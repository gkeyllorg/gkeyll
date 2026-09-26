#!/usr/bin/env python3
"""Build a Markdown CI report from Jenkins workspace artifacts and publish it
to GitHub as a pull-request comment (or a commit comment when the run has no
pull request). Standard library only.

    github_report.py build   --platform NAME --context CTX --result STATE [--pr N] [--output ci-report.md]
    github_report.py publish --report ci-report.md --commit SHA --context CTX [--pr N] [--repo OWNER/NAME]

`build` reads the summary files the Jenkinsfiles already write into the
workspace root (ci-selection.txt, ci-candidate-commit.txt, ci-failure-summary.txt,
ci-regression-summary.txt, ...). `publish` needs GITHUB_TOKEN in the environment
and updates the previous report comment for the same context in place, so each
CI machine owns exactly one comment per pull request or commit.

The Jenkinsfiles run `publish` only from a trusted copy of this file (main, or
the team-workstation trusted checkout), never from the candidate checkout, so
that the token is never bound while candidate code runs.
"""

import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request

API = "https://api.github.com"
DEFAULT_REPO = "gkeyllorg/gkeyll"
MARKER_FORMAT = '<!-- gkeyll-ci-report context="{}" -->'
COMMENT_LIMIT = 60000  # GitHub allows 65536 characters; keep headroom.
DETAIL_LIMIT = 6000


# ---- helpers ---------------------------------------------------------------

def read_text(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            return f.read()
    except OSError:
        return ""


def read_kv(path):
    """Parse key=value lines. Repeated keys accumulate into a list."""
    values = {}
    for line in read_text(path).splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        key, value = key.strip(), value.strip()
        values.setdefault(key, []).append(value)
    return values


def first(values, key, default=""):
    return values.get(key, [default])[0] if values.get(key) else default


def short(sha):
    return sha[:7] if re.fullmatch(r"[0-9a-fA-F]{40}", sha or "") else (sha or "?")


def code(text):
    return "`" + text.replace("`", "'") + "`"


def truncate(text, limit):
    if len(text) <= limit:
        return text
    return text[:limit] + "\n... [truncated {} characters]".format(len(text) - limit)


def fenced(text, lang="text"):
    return "```{}\n{}\n```".format(lang, text.rstrip("\n").replace("```", "'''"))


# ---- build -----------------------------------------------------------------

STATE_LABEL = {
    "success": ":white_check_mark: passed",
    "failure": ":x: failed",
    "error": ":warning: errored or aborted",
    "pending": ":hourglass: running",
}


def regression_section(title, summary_file):
    values = read_kv(summary_file)
    if not values:
        return ""
    passed = first(values, "c_regression_passed", "?")
    acked = first(values, "c_regression_acknowledged", "0")
    unacked = first(values, "c_regression_unacknowledged", "0")
    lines = ["**{}:** {} passed, {} acknowledged, {} unacknowledged".format(
        title, passed, acked, unacked)]
    failures = values.get("c_regression_failure", [])
    if failures:
        lines += ["", "| Failing test | Status |", "| --- | --- |"]
        for entry in failures:
            name, _, status = entry.rpartition(":")
            lines.append("| {} | {} |".format(code(name or entry), status or "fail"))
    acked_tests = values.get("c_regression_acknowledged_test", [])
    if acked_tests:
        lines += ["", "Acknowledged diffs (listed in expected_regression_diffs.txt): "
                  + ", ".join(code(t) for t in acked_tests)]
    return "\n".join(lines)


def failure_detail():
    detail = read_text("ci-failure-detail.txt").strip()
    if detail:
        return detail
    # HPC pipelines: fall back to the most recent Slurm output, if any.
    outs = [p for p in os.listdir(".") if p.startswith("slurm-") and p.endswith(".out")]
    if not outs:
        return ""
    newest = max(outs, key=lambda p: os.path.getmtime(p))
    tail = read_text(newest).splitlines()[-40:]
    return "{} (last 40 lines)\n{}".format(newest, "\n".join(tail))


def timing_section():
    values = read_kv("ci-timing-summary.txt")
    rows = [(k, first(values, k)) for k in values if k.endswith("_seconds")]
    rows = [(k, v) for k, v in rows if v and v != "not-recorded"]
    if not rows:
        return ""
    body = ["| Step | Seconds |", "| --- | ---: |"]
    body += ["| {} | {} |".format(k[:-len("_seconds")].replace("_", " "), v) for k, v in rows]
    return "<details><summary>Timings</summary>\n\n" + "\n".join(body) + "\n\n</details>"


def build_report(args):
    selection = read_kv("ci-selection.txt")
    candidate = read_text("ci-candidate-commit.txt").strip()
    baseline = read_text("ci-baseline-commit.txt").strip()
    failure = read_kv("ci-failure-summary.txt")
    build_number = os.environ.get("BUILD_NUMBER", "?")
    build_url = os.environ.get("BUILD_URL", "")
    node = os.environ.get("NODE_NAME") or os.environ.get("HOSTNAME") or ""

    parts = [MARKER_FORMAT.format(args.context)]
    parts.append("### Gkeyll CI on {}: {}".format(
        code(args.platform), STATE_LABEL.get(args.result, args.result)))

    candidate_sel = first(selection, "candidate_selector", "PR #{}".format(args.pr) if args.pr else "?")
    baseline_sel = first(selection, "baseline_selector", "?")
    meta = ["**Candidate:** {} @ {}".format(candidate_sel, code(short(candidate))),
            "**Baseline:** {} @ {}".format(baseline_sel, code(short(baseline)))]
    where = "Jenkins build #{}".format(build_number)
    if node:
        where += " on {}".format(code(node))
    if build_url:
        where += " ({}, reachable only from that controller)".format(build_url)
    meta.append(where)
    parts.append("  \n".join(meta))

    if args.result != "success":
        stage = first(failure, "stage", "unknown")
        message = first(failure, "message", "")
        parts.append("**Failed at stage:** {}{}".format(stage, " — " + message if message else ""))
        detail = failure_detail()
        if detail:
            parts.append("<details><summary>Failure detail</summary>\n\n"
                         + fenced(truncate(detail, DETAIL_LIMIT)) + "\n\n</details>")

    for title, path in (("Serial C regressions", "ci-regression-summary.txt"),
                        ("Parallel C regressions", "ci-parallel-regression-summary.txt")):
        section = regression_section(title, path)
        if section:
            parts.append(section)

    timings = timing_section()
    if timings:
        parts.append(timings)

    parts.append("_Artifacts: `./ci/jenkins/gkeyll-ci.sh {} artifact --build {} --fetch` "
                 "for anyone with access to this controller._".format(args.platform, build_number))

    report = "\n\n".join(parts) + "\n"
    report = truncate(report, COMMENT_LIMIT)
    with open(args.output, "w", encoding="utf-8") as f:
        f.write(report)
    print("Wrote {} ({} characters)".format(args.output, len(report)))


# ---- publish ---------------------------------------------------------------

def api(method, url, token, payload=None):
    data = json.dumps(payload).encode() if payload is not None else None
    request = urllib.request.Request(url, data=data, method=method, headers={
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
        "Authorization": "Bearer " + token,
        "Content-Type": "application/json",
        "User-Agent": "gkeyll-ci-report",
    })
    try:
        with urllib.request.urlopen(request, timeout=45) as response:
            body = response.read().decode("utf-8", "replace")
            link = response.headers.get("Link", "")
            return (json.loads(body) if body else None), link
    except urllib.error.HTTPError as err:
        body = err.read().decode("utf-8", "replace")[:500]
        sys.exit("GitHub API {} {} failed: HTTP {} {}".format(method, url, err.code, body))


def next_link(link_header):
    for part in link_header.split(","):
        match = re.search(r'<([^>]+)>;\s*rel="next"', part)
        if match:
            return match.group(1)
    return None


def find_existing(url, token, marker):
    """Return the id of the first comment whose body contains the marker."""
    url += ("&" if "?" in url else "?") + "per_page=100"
    for _ in range(20):
        comments, link = api("GET", url, token)
        for comment in comments or []:
            if marker in (comment.get("body") or ""):
                return comment["id"]
        url = next_link(link)
        if not url:
            break
    return None


def publish_report(args):
    token = os.environ.get("GITHUB_TOKEN")
    if not token:
        sys.exit("GITHUB_TOKEN is not set")
    if not re.fullmatch(r"[0-9a-fA-F]{40}", args.commit):
        sys.exit("--commit must be a full 40-character SHA")
    body = read_text(args.report)
    marker = MARKER_FORMAT.format(args.context)
    if marker not in body:
        sys.exit("report {} does not carry the marker for context {}".format(args.report, args.context))
    repo = args.repo
    payload = {"body": body}

    if args.pr:
        list_url = "{}/repos/{}/issues/{}/comments".format(API, repo, args.pr)
        update_url = "{}/repos/{}/issues/comments/{{}}".format(API, repo)
        target = "PR #{}".format(args.pr)
    else:
        list_url = "{}/repos/{}/commits/{}/comments".format(API, repo, args.commit)
        update_url = "{}/repos/{}/comments/{{}}".format(API, repo)
        target = "commit {}".format(short(args.commit))

    existing = find_existing(list_url, token, marker)
    if existing:
        result, _ = api("PATCH", update_url.format(existing), token, payload)
        action = "Updated"
    else:
        result, _ = api("POST", list_url, token, payload)
        action = "Created"
    print("{} CI report comment on {}: {}".format(action, target, (result or {}).get("html_url", "")))


# ---- main ------------------------------------------------------------------

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    build = sub.add_parser("build", help="write the Markdown report from workspace artifacts")
    build.add_argument("--platform", required=True, help="CI platform name, e.g. personal, team, stellar_cpu")
    build.add_argument("--context", required=True, help="GitHub status context this report belongs to")
    build.add_argument("--result", required=True, choices=sorted(STATE_LABEL))
    build.add_argument("--pr", default="", help="pull-request number, if the run tests a PR")
    build.add_argument("--output", default="ci-report.md")
    build.set_defaults(func=build_report)

    publish = sub.add_parser("publish", help="create or update the GitHub comment carrying the report")
    publish.add_argument("--report", default="ci-report.md")
    publish.add_argument("--commit", required=True, help="candidate commit SHA")
    publish.add_argument("--context", required=True)
    publish.add_argument("--pr", default="", help="pull-request number; omitted for commit comments")
    publish.add_argument("--repo", default=DEFAULT_REPO)
    publish.set_defaults(func=publish_report)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
