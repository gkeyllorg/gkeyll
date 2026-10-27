#!/usr/bin/env python3
"""Coordinate persistent Jenkins CI trees without deleting a live run's files."""

import argparse
import fcntl
import getpass
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request


PLATFORMS = {"personal", "team-workstation", "stellar-cpu", "perlmutter-gpu"}
SLURM_PLATFORMS = {"stellar-cpu", "perlmutter-gpu"}
SHA = re.compile(r"[0-9a-f]{40}\Z")
TAG = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
WORKSPACE_SUFFIX = re.compile(r"(?:tmp|script|\d+(?:@tmp)?)\Z")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def remove_tree(path):
    require(not path.is_symlink(), f"refusing symlink: {path}")
    if path.is_dir():
        print(f"cache-manager: removing {path}", file=sys.stderr)
        shutil.rmtree(path)


def write_json(path, value):
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, sort_keys=True) + "\n")
    os.replace(temporary, path)


class CacheManager:
    def __init__(self, args):
        self.root = Path(args.root).absolute()
        require(self.root != Path("/") and ".." not in self.root.parts, "unsafe CI root")
        require(args.platform in PLATFORMS, "invalid platform")
        require(TAG.fullmatch(args.tag) is not None and ".." not in args.tag, "invalid build tag")
        require(args.number.isdecimal() and int(args.number) > 0, "invalid build number")
        require(args.job and all(part not in ("", ".", "..") for part in args.job.split("/")), "invalid job name")
        self.platform = args.platform
        self.tag = args.tag
        self.job = args.job
        self.number = args.number
        self.url = args.url.rstrip("/")
        parsed = urllib.parse.urlsplit(self.url)
        require(parsed.scheme in ("http", "https") and parsed.netloc and "/job/" in parsed.path, "invalid BUILD_URL")
        build_controller = urllib.parse.urlunsplit((parsed.scheme, parsed.netloc, parsed.path.split("/job/", 1)[0], "", ""))
        controller = urllib.parse.urlsplit(args.controller_url or build_controller)
        require(controller.scheme in ("http", "https") and controller.netloc
                and not controller.username and not controller.password
                and not controller.query and not controller.fragment and "/job/" not in controller.path,
                "invalid controller URL")
        self.controller = urllib.parse.urlunsplit((controller.scheme, controller.netloc,
                                                   controller.path.rstrip("/"), "", ""))
        self.auth_file = Path(args.auth_file) if args.auth_file else None
        self.candidates = self.root / "candidate-cache" / self.platform
        self.baselines = self.root / "baseline-cache" / self.platform
        self.runs = self.root / "cache-runs" / self.platform
        self.workspaces = self.root / "workspaces" / self.job
        self.record_path = self.runs / f"{self.tag}.json"
        self.latest_path = self.runs / "latest-run"
        self.live_cache = {}
        self.slurm_dirs = None

    def locked(self, operation):
        lock_dir = self.root / "cache-locks"
        lock_dir.mkdir(parents=True, exist_ok=True)
        lock_path = lock_dir / f"{self.platform}.lock"
        with lock_path.open("a+") as lock_file:
            deadline = time.monotonic() + 60
            while True:
                try:
                    fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except BlockingIOError:
                    require(time.monotonic() < deadline, f"cache lock timed out: {lock_path}")
                    time.sleep(0.25)
            try:
                return operation()
            finally:
                fcntl.flock(lock_file, fcntl.LOCK_UN)

    def read_record(self, tag):
        path = self.runs / f"{tag}.json"
        if path.is_file():
            data = json.loads(path.read_text())
            require(data.get("tag") == tag and isinstance(data.get("job"), str)
                    and str(data.get("number", "")).isdecimal(), f"invalid run record: {path}")
            return data
        # Pre-migration cluster workspaces have no run record.
        prefix = f"jenkins-{self.job.replace('/', '-')}-"
        if tag.startswith(prefix) and tag[len(prefix):].isdecimal():
            return {"tag": tag, "job": self.job, "number": tag[len(prefix):]}
        return None

    def active_slurm_directories(self):
        if self.platform not in SLURM_PLATFORMS:
            return set()
        if self.slurm_dirs is None:
            result = subprocess.run(
                ["squeue", "--noheader", "--user", getpass.getuser(), "--format", "%Z"],
                text=True, capture_output=True, check=False,
            )
            require(result.returncode == 0, f"cannot query Slurm: {result.stderr.strip()}")
            directories = {Path(line.strip()) for line in result.stdout.splitlines() if line.strip()}
            require(all(path.is_absolute() and str(path) != "/" for path in directories),
                    "unknown Slurm working directory; preserving caches")
            self.slurm_dirs = directories
        return self.slurm_dirs

    def slurm_uses_tag(self, tag):
        for path in self.active_slurm_directories():
            if self.root in path.parents and tag in path.relative_to(self.root).parts:
                return True
        return False

    def is_live(self, tag):
        if tag == self.tag:
            return True
        if tag in self.live_cache:
            return self.live_cache[tag]
        record = self.read_record(tag)
        require(record is not None, f"cannot identify owner of {tag}; preserve its files")
        job_path = "".join("/job/" + urllib.parse.quote(part, safe="") for part in record["job"].split("/"))
        url = f"{self.controller}{job_path}/{record['number']}/api/json?tree=building,result"
        request = urllib.request.Request(url)
        if self.auth_file:
            import base64
            credentials = self.auth_file.read_text().strip()
            require(":" in credentials, "invalid Jenkins read credential file")
            request.add_header("Authorization", "Basic " + base64.b64encode(credentials.encode()).decode())
        try:
            with urllib.request.urlopen(request, timeout=15) as response:
                building = json.load(response).get("building")
            require(isinstance(building, bool), f"Jenkins returned no build state: {url}")
        except urllib.error.HTTPError as error:
            if error.code == 404:
                building = False
            else:
                raise RuntimeError(f"Jenkins build-state query failed ({error.code}): {url}") from error
        except (urllib.error.URLError, TimeoutError) as error:
            raise RuntimeError(f"Jenkins build-state query unavailable: {url}") from error
        live = building or self.slurm_uses_tag(tag)
        self.live_cache[tag] = live
        return live

    def candidate_tags(self):
        if not self.candidates.is_dir():
            return []
        require(not self.candidates.is_symlink(), "candidate root is a symlink")
        return sorted(path for path in self.candidates.iterdir() if path.is_dir() and not path.is_symlink())

    def prune_candidates(self, keep_latest=False):
        paths = self.candidate_tags()
        latest = ""
        if keep_latest:
            latest = json.loads(self.latest_path.read_text()).get("tag", "")
            require(TAG.fullmatch(latest) is not None, "invalid latest run marker")
        for path in paths:
            tag = path.name
            if SHA.fullmatch(tag):
                # Legacy cache has one candidate per platform. Do not remove it
                # while a pre-migration Slurm allocation may still use it.
                require(not self.active_slurm_directories(), f"legacy candidate may be active: {path}")
                manifest = path / "candidate-manifest.txt"
                require(manifest.is_file() and "result=building" not in manifest.read_text(),
                        f"legacy candidate may be active: {path}")
                remove_tree(path)
            elif TAG.fullmatch(tag):
                record = self.read_record(tag)
                require(record is not None, f"cannot identify owner of {path}; preserve its files")
                if record.get("state") == "finished":
                    if tag == latest or self.slurm_uses_tag(tag):
                        print(f"cache-manager: preserving selected candidate {path}", file=sys.stderr)
                    else:
                        remove_tree(path)
                elif self.is_live(tag):
                    print(f"cache-manager: preserving active candidate {path}", file=sys.stderr)
                else:
                    remove_tree(path)
            elif not TAG.fullmatch(tag):
                raise RuntimeError(f"unexpected candidate entry: {path}")

    def baseline_valid(self, path, sha):
        manifest = path / "cache-manifest.txt"
        if not manifest.is_file():
            return False
        values = dict(line.split("=", 1) for line in manifest.read_text().splitlines() if "=" in line)
        if values.get("format_version") != "1" or values.get("platform") != self.platform or values.get("baseline_commit") != sha:
            return False
        if not (path / "gkeyll" / ".git").is_dir() or not os.access(path / "gkylsoft" / "gkeyll" / "bin" / "gkeyll", os.X_OK):
            return False
        for layer in ("moments", "vlasov", "gyrokinetic", "pkpm"):
            for subdir in ("", "parallel-c-4"):
                if not (path / "gkylsoft" / "gkeyll-results" / subdir / layer / "creg-accepted").is_dir():
                    return False
        result = subprocess.run(["git", "-C", str(path / "gkeyll"), "rev-parse", "HEAD"], capture_output=True, text=True)
        return result.returncode == 0 and result.stdout.strip() == sha

    def baseline_users(self):
        users = set()
        if self.runs.is_dir():
            for path in self.runs.glob("*.json"):
                if path.is_symlink():
                    raise RuntimeError(f"run record is a symlink: {path}")
                record = self.read_record(path.stem)
                if record.get("baseline") and (self.slurm_uses_tag(path.stem)
                    if record.get("state") == "finished" else self.is_live(path.stem)):
                    users.add(record["baseline"])
        return users

    def prune_baselines(self, selected="", keep_previous=True):
        if not self.baselines.is_dir():
            return
        require(not self.baselines.is_symlink(), "baseline root is a symlink")
        users = self.baseline_users()
        entries = sorted((path for path in self.baselines.iterdir() if path.is_dir() and not path.is_symlink()),
                         key=lambda path: path.stat().st_mtime, reverse=True)
        keep = selected if selected and self.baseline_valid(self.baselines / selected, selected) else ""
        if not keep and keep_previous:
            keep = next((path.name for path in entries if SHA.fullmatch(path.name) and self.baseline_valid(path, path.name)), "")
        for path in entries:
            require(SHA.fullmatch(path.name), f"unexpected baseline entry: {path}")
            if path.name not in users and path.name != keep:
                remove_tree(path)

    def prune_workspaces(self):
        if self.platform not in SLURM_PLATFORMS or not self.workspaces.is_dir():
            return
        require(not self.workspaces.is_symlink(), "workspace root is a symlink")
        for path in self.workspaces.iterdir():
            if not path.is_dir() or path.is_symlink():
                continue
            tag, separator, suffix = path.name.partition("@")
            require(TAG.fullmatch(tag) and (not separator or WORKSPACE_SUFFIX.fullmatch(suffix)),
                    f"unexpected workspace entry: {path}")
            if self.is_live(tag):
                print(f"cache-manager: preserving active workspace {path}", file=sys.stderr)
            else:
                remove_tree(path)

    def prune_controls(self):
        parent = self.root / "cache-control"
        if not parent.is_dir():
            return
        require(not parent.is_symlink(), "control root is a symlink")
        for path in parent.iterdir():
            if path.is_dir() and not path.is_symlink() and path.name != self.tag:
                require(TAG.fullmatch(path.name), f"unexpected control entry: {path}")
                if self.read_record(path.name) is not None and not self.is_live(path.name):
                    remove_tree(path)

    def begin(self):
        def operation():
            self.runs.mkdir(parents=True, exist_ok=True)
            require(not self.record_path.exists(), f"run already registered: {self.tag}")
            write_json(self.record_path, {"tag": self.tag, "job": self.job, "number": self.number,
                                          "state": "building", "started_ns": time.time_ns()})
            write_json(self.latest_path, {"tag": self.tag})
            self.active_slurm_directories()
            self.prune_candidates()
            self.prune_workspaces()
            self.prune_controls()
            self.prune_baselines()
        self.locked(operation)

    def claim_baseline(self, sha):
        require(SHA.fullmatch(sha), "invalid baseline SHA")
        def operation():
            record = self.read_record(self.tag)
            require(record is not None, "run is not registered")
            record["baseline"] = sha
            write_json(self.record_path, record)
            self.baselines.mkdir(parents=True, exist_ok=True)
            self.active_slurm_directories()
            path = self.baselines / sha
            if self.baseline_valid(path, sha):
                os.utime(path, None)
                self.prune_baselines(sha)
                return "hit"
            if path.exists():
                owner = path / "cache-owner.txt"
                if owner.is_file() and self.is_live(owner.read_text().strip()):
                    return "wait"
                remove_tree(path)
            self.prune_baselines(keep_previous=False)
            (path / "gkeyll").mkdir(parents=True)
            (path / "cache-owner.txt").write_text(self.tag + "\n")
            return "build"
        return self.locked(operation)

    def publish_baseline(self, sha):
        require(SHA.fullmatch(sha), "invalid baseline SHA")
        def operation():
            path = self.baselines / sha
            require((path / "cache-owner.txt").read_text().strip() == self.tag, "baseline has another owner")
            require(self.baseline_valid(path, sha), "baseline build/results are incomplete")
            (path / "cache-owner.txt").unlink()
            self.prune_baselines(sha)
        self.locked(operation)

    def abort_baseline(self, sha):
        require(SHA.fullmatch(sha), "invalid baseline SHA")
        def operation():
            path = self.baselines / sha
            owner = path / "cache-owner.txt"
            if owner.is_file() and owner.read_text().strip() == self.tag and not self.baseline_valid(path, sha):
                require(not self.slurm_uses_tag(self.tag), "baseline still used by Slurm")
                remove_tree(path)
        self.locked(operation)

    def finish(self):
        def operation():
            record = self.read_record(self.tag)
            if record:
                record["state"] = "finished"
                write_json(self.record_path, record)
            self.active_slurm_directories()
            self.prune_candidates(keep_latest=True)
            self.prune_baselines()
            for path in self.runs.glob("*.json"):
                if path.stem != self.tag and not self.is_live(path.stem) and not (self.candidates / path.stem).exists():
                    path.unlink()
        self.locked(operation)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("begin", "claim-baseline", "publish-baseline", "abort-baseline", "finish"))
    parser.add_argument("--root", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--job", required=True)
    parser.add_argument("--number", required=True)
    parser.add_argument("--url", required=True)
    parser.add_argument("--controller-url", default="")
    parser.add_argument("--auth-file", default="")
    parser.add_argument("--sha", default="")
    args = parser.parse_args()
    manager = CacheManager(args)
    try:
        result = getattr(manager, args.command.replace("-", "_"))(args.sha) if args.command in ("claim-baseline", "publish-baseline", "abort-baseline") else getattr(manager, args.command)()
        if result is not None:
            print(result)
    except (OSError, ValueError, RuntimeError) as error:
        print(f"cache-manager: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
