#!/usr/bin/env python3
"""Prepare formatted three-way merge candidates without changing the checkout."""

import argparse
import difflib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ours", default="HEAD")
    parser.add_argument("--theirs", default="MERGE_HEAD")
    parser.add_argument("--base", help="Explicit base for histories with multiple merge bases")
    parser.add_argument("--style-ref", help="Revision providing formatting tools; defaults to theirs")
    parser.add_argument("--clang-format", default="clang-format")
    parser.add_argument("--output", type=Path, help="New output directory; defaults to a temporary one")
    args = parser.parse_args()
    root = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip())

    def git(*command):
        return subprocess.check_output(["git", *command], cwd=root)

    def revision(ref):
        return git("rev-parse", "--verify", "--end-of-options", ref + "^{commit}").decode().strip()

    refs = {"ours": revision(args.ours), "theirs": revision(args.theirs)}
    bases = [revision(args.base)] if args.base else git(
        "merge-base", "--all", refs["ours"], refs["theirs"]
    ).decode().splitlines()
    if len(bases) != 1:
        parser.error("Choose --base explicitly: the histories do not have one unique merge base")
    refs["base"] = bases[0]
    style_ref = revision(args.style_ref) if args.style_ref else refs["theirs"]
    config = git("show", style_ref + ":.clang-format")
    strip_source = git("show", style_ref + ":ci/strip-trailing-commas.py")
    format_script = git("show", style_ref + ":ci/format-all.sh").decode()
    pinned = re.search(r'PINNED_VERSION="([0-9.]+)"', format_script)
    if not pinned:
        parser.error("Cannot read PINNED_VERSION from the selected revision's ci/format-all.sh")
    version = subprocess.check_output([args.clang_format, "--version"], text=True)
    actual = re.search(r"\b\d+\.\d+\.\d+\b", version)
    if not actual or actual.group() != pinned.group(1):
        parser.error(f"Use clang-format {pinned.group(1)}; found {version.strip()}")

    output = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix="gkeyll-merge-"))
    if args.output:
        output.mkdir(parents=True, exist_ok=False)
    (output / ".clang-format").write_bytes(config)
    strip_path = output / "strip_trailing_commas.py"
    strip_path.write_bytes(strip_source)
    spec = importlib.util.spec_from_file_location("merge_strip_commas", strip_path)
    strip = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(strip)

    def normalize(source, path):
        for _ in range(5):
            result = subprocess.run(
                [args.clang_format, "--style=file:" + str(output / ".clang-format"),
                 "--assume-filename=" + str(root / path)],
                input=strip.strip_trailing_commas(source), text=True, capture_output=True, check=True
            ).stdout
            if result == source:
                return result
            source = result
        raise RuntimeError("Formatting did not stabilize: " + path)

    def write(path, text):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    trees = {}
    for label, ref in refs.items():
        trees[label] = {}
        for entry in git("ls-tree", "-r", "-z", ref).decode().split("\0"):
            if entry:
                metadata, path = entry.split("\t", 1)
                trees[label][path] = metadata.split()[0]
    paths = git("diff", "--name-only", "--no-renames", "-z", refs["base"], refs["ours"])
    paths = [path for path in paths.decode().split("\0") if path]
    rows, excluded = [], []
    for path in paths:
        if (Path(path).suffix not in {".c", ".h", ".cpp", ".hpp", ".cu", ".cuh"}
                or "ker" in Path(path).parts[:-1] or path.startswith("core/minus/")
                or any(tree.get(path, "100644") not in {"100644", "100755"}
                       for tree in trees.values())):
            excluded.append(path)
            continue
        versions, normalized = {}, {}
        for label, ref in refs.items():
            versions[label] = None
            if path in trees[label]:
                source = git("show", ref + ":" + path).decode()
                if source not in normalized:
                    normalized[source] = normalize(source, path)
                versions[label] = normalized[source]
            write(output / label / path, versions[label] or "")
        base, ours, theirs = (versions[label] for label in ("base", "ours", "theirs"))
        row = {"path": path, "incoming_format_only": base == theirs, "conflicts": 0}
        candidate = None
        if ours is None:
            row["status"] = "delete-clean" if theirs is None or theirs == base else "delete-conflict"
        elif theirs is None:
            if base is None:
                row["status"], candidate = "ours-added", ours
            else:
                row["status"] = "delete-clean" if ours == base else "delete-conflict"
        else:
            merge = subprocess.run(
                ["git", "merge-file", "-p", "--diff3", "-L", "ours", "-L", "base", "-L", "theirs",
                 str(output / "ours" / path), str(output / "base" / path), str(output / "theirs" / path)],
                capture_output=True, text=True
            )
            if merge.returncode < 0 or merge.returncode > 127:
                raise RuntimeError(merge.stderr)
            row["status"] = "conflict" if merge.returncode else "clean"
            row["conflicts"] = sum(line == "<<<<<<< ours" for line in merge.stdout.splitlines())
            candidate = merge.stdout
            if not merge.returncode:
                candidate = normalize(candidate, path)
        if candidate is not None:
            write(output / "merged" / path, candidate)
        for label in ("ours", "theirs"):
            diff = "".join(difflib.unified_diff(
                (base or "").splitlines(True), (versions[label] or "").splitlines(True),
                fromfile="base/" + path, tofile=label + "/" + path
            ))
            folder = "ours-changes" if label == "ours" else "incoming-changes"
            write(output / folder / (path + ".diff"), diff)
        rows.append(row)
        if row["status"] in {"conflict", "delete-conflict"}:
            print(f"{row['status']}: {path}")

    manifest = {"refs": refs, "style_ref": style_ref, "clang_format": pinned.group(1),
                "files": rows, "excluded_paths": excluded}
    write(output / "manifest.json", json.dumps(manifest, indent=2) + "\n")
    manual = sum(row["status"] in {"conflict", "delete-conflict"} for row in rows)
    print(f"Prepared {len(rows)} paths; {manual} need conflict review. Output: {output}")


if __name__ == "__main__":
    main()
