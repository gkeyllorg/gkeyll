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
import datetime
import json
import os
import re
import socket
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


def describe(selector, commit):
    """'main @ `abc1234`', or just the short SHA when the selector is a SHA,
    or 'main (not reached)' when the run ended before that tree was built."""
    if not commit:
        return "{} (not reached)".format(selector)
    if re.fullmatch(r"[0-9a-fA-F]{40}", selector or ""):
        return code(short(commit))
    return "{} @ {}".format(selector, code(short(commit)))


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
    layers = [entry.split(":") for entry in values.get("c_regression_layer", [])]
    layers = [row for row in layers if len(row) == 4]
    if layers:
        lines += ["", "| Layer | Passed | Acknowledged | Failed |", "| --- | ---: | ---: | ---: |"]
        for layer, p, a, f in layers:
            mark = " :x:" if f != "0" else ""
            lines.append("| {} | {} | {} | {}{} |".format(layer, p, a, f, mark))
    failures = values.get("c_regression_failure", [])
    if failures:
        details = {}
        for entry in values.get("c_regression_failure_detail", []):
            name, _, line = entry.partition("|")
            details.setdefault(name, []).append(line)
        lines += ["", "| Failing test | Status |", "| --- | --- |"]
        for entry in failures:
            name, _, status = entry.rpartition(":")
            lines.append("| {} | {} |".format(code(name or entry), status or "fail"))
        for entry in failures:
            name = entry.rpartition(":")[0] or entry
            if details.get(name):
                lines += ["", "<details><summary>{}: files that differ from the baseline</summary>".format(code(name)),
                          "", fenced("\n".join(details[name])), "", "</details>"]
    acked_tests = values.get("c_regression_acknowledged_test", [])
    if acked_tests:
        lines += ["", "Acknowledged diffs (listed in expected_regression_diffs.txt): "
                  + ", ".join(code(t) for t in acked_tests)]
    return "\n".join(lines)


ERROR_LINE = re.compile(
    r"\[ FAILED \]|\.\.\. failed$|^FAILED:|^Failed tests:|^FAIL |error:|undefined reference"
    r"|^make(\[\d+\])?: \*\*\*|[Ss]egmentation fault|Abort trap|Bus error|[Tt]imed? ?out")


def extract_errors(log_path, max_lines=40, tail_lines=15):
    """Pull the informative lines out of a build/test log: matches of ERROR_LINE
    plus the indented lines that follow a 'Failed tests:' summary."""
    lines = read_text(log_path).splitlines()
    picked, in_failed_block = [], False
    for number, line in enumerate(lines, 1):
        if in_failed_block and line.startswith("  "):
            picked.append("{}: {}".format(number, line))
            continue
        in_failed_block = line.startswith("Failed tests:")
        if ERROR_LINE.search(line):
            picked.append("{}: {}".format(number, line))
        if len(picked) >= max_lines:
            picked.append("... [more matches omitted]")
            break
    out = ["--- error lines ---"] + (picked or ["(no recognised error lines)"])
    out += ["--- last {} lines ---".format(tail_lines)] + lines[-tail_lines:]
    return "\n".join(out)


LAYER_ORDER = ["core", "moments", "vlasov", "gyrokinetic", "pkpm"]


def unit_results(results_file):
    """Parse 'PASS <layer>: <test>' / 'FAIL <layer>: <test>' lines into
    {layer: {"passed": n, "failed": [test, ...]}}, or None if absent."""
    text = read_text(results_file)
    if not text:
        return None
    layers = {}
    for line in text.splitlines():
        match = re.match(r"^(PASS|FAIL) (\S+): (\S+)", line)
        if not match:
            continue
        entry = layers.setdefault(match.group(2), {"passed": 0, "failed": []})
        if match.group(1) == "PASS":
            entry["passed"] += 1
        else:
            entry["failed"].append(match.group(3))
    return layers


def assertion_lines(test, log_file, limit=8):
    """Lines the unit-test framework printed for a failing test: the
    '[ FAILED ]' check names and 'file.c:NN: Check ... failed' lines."""
    picked, want = [], test + ".c:"
    lines = read_text(log_file).splitlines()
    for index, line in enumerate(lines):
        if want in line:
            if index and "[ FAILED ]" in lines[index - 1] and (not picked or picked[-1] != lines[index - 1].strip()):
                picked.append(lines[index - 1].strip())
            picked.append(line.strip())
        if len(picked) >= limit:
            picked.append("...")
            break
    return picked


def unit_section(label, results_file, log_file):
    layers = unit_results(results_file)
    if layers is None:
        return "", 0, []
    order = LAYER_ORDER + sorted(set(layers) - set(LAYER_ORDER))
    rows, nfail, failing = [], 0, []
    for layer in order:
        if layer not in layers:
            continue
        entry = layers[layer]
        failed = entry["failed"]
        nfail += len(failed)
        failing += ["{}: {}".format(layer, t) for t in failed]
        mark = " :x:" if failed else ""
        rows.append("| {} | {} | {}{} | {} |".format(
            layer, entry["passed"], len(failed), mark, ", ".join(code(t) for t in failed)))
    total = sum(e["passed"] for e in layers.values())
    head = "{} unit tests: {} passed, {} failed".format(label, total, nfail)
    table = ["| Layer | Passed | Failed | Failing tests |", "| --- | ---: | ---: | --- |"] + rows
    if nfail:
        body = ["**{}**".format(head), ""] + table
    else:  # all green: keep the per-layer table one click away; no markup in <summary>
        body = ["<details><summary>{}</summary>".format(head), ""] + table + ["", "</details>"]
    for layer in order:
        for test in layers.get(layer, {}).get("failed", []):
            detail = assertion_lines(test, log_file)
            if detail:
                body += ["", "<details><summary>{}: failed checks</summary>".format(code(test)),
                         "", fenced("\n".join(detail)), "", "</details>"]
    return "\n".join(body), nfail, failing


def status_description(result, stage, unit_fail, unit_failing, regression_files):
    """One line for the GitHub status: name what failed, or what passed."""
    if result == "success":
        return "Passed: unit tests and C regressions."
    if result == "error":
        return "Aborted or errored at stage: {}.".format(stage)
    if unit_fail:
        by_layer = {}
        for item in unit_failing:
            by_layer[item.split(":")[0]] = by_layer.get(item.split(":")[0], 0) + 1
        return "Failed: unit tests ({}).".format(
            ", ".join("{} {}".format(l, n) for l, n in by_layer.items()))
    failing = []
    for path in regression_files:
        failing += [e.rpartition(":")[0] for e in read_kv(path).get("c_regression_failure", [])]
    if failing:
        by_layer = {}
        for name in failing:
            by_layer[name.split("/")[0]] = by_layer.get(name.split("/")[0], 0) + 1
        return "Failed: C regressions ({}).".format(
            ", ".join("{} {}".format(l, n) for l, n in by_layer.items()))
    return "Failed at stage: {}.".format(stage)


def failure_detail():
    detail = read_text("ci-failure-detail.txt").strip()
    if detail:
        # The pipeline records which log the failed command wrote; re-scan it
        # here so the extraction heuristics live in one testable place.
        match = re.search(r"^Full log: (\S+)", detail, re.M)
        command = re.search(r"^Command: (.*)$", detail, re.M)
        if match and os.path.isfile(match.group(1)):
            head = ["Command: " + command.group(1)] if command else []
            head.append("Full log: {} (archived with the build)".format(match.group(1)))
            return "\n".join(head) + "\n" + extract_errors(match.group(1))
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
    node = os.environ.get("NODE_NAME", "")
    if node in ("", "built-in", "master"):
        node = os.environ.get("HOSTNAME") or socket.gethostname().split(".")[0]

    parts = [MARKER_FORMAT.format(args.context)]
    parts.append("### Gkeyll CI on {}: {}".format(
        code(args.platform), STATE_LABEL.get(args.result, args.result)))

    candidate_sel = first(selection, "candidate_selector", "PR #{}".format(args.pr) if args.pr else "?")
    baseline_sel = first(selection, "baseline_selector", "?")
    meta = ["**Candidate:** " + describe(candidate_sel, candidate),
            "**Baseline:** " + describe(baseline_sel, baseline)]
    # No controller URL: every controller is loopback-only, so a link would be
    # dead for everyone but the machine owner. The build number is what that
    # owner needs to fetch artifacts (see the footer).
    finished = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    meta.append("**Run:** {} build #{} on {}, finished {}".format(
        code(args.platform), build_number, code(node), finished))
    parts.append("  \n".join(meta))

    stage = first(failure, "stage", "unknown")
    if args.result != "success":
        message = re.sub(r"^(?:[A-Za-z_$][\w$]*\.)+[A-Z]\w*(?:Exception|Error): ", "", first(failure, "message", ""))
        parts.append("**Failed at stage:** {}{}".format(stage, " — " + message if message else ""))

    unit_fail, unit_failing = 0, []
    for label, results_file, log_file in (("Candidate", "candidate-unit-results.txt", "candidate-unit-test.log"),
                                          ("Baseline", "baseline-unit-results.txt", "baseline-unit-test.log")):
        section, nfail, failing = unit_section(label, results_file, log_file)
        if section:
            parts.append(section)
            unit_fail += nfail
            unit_failing += failing

    if args.result != "success":
        detail = failure_detail()
        # With a structured unit-test table above, the raw log excerpt is only
        # needed when the failure was something else (compile error, crash...).
        if detail and not unit_fail:
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

    parts.append("_Logs and regression databases stay on that machine; its owner can fetch them with "
                 "`./ci/jenkins/gkeyll-ci.sh {} artifact --build {} --fetch`._".format(args.platform, build_number))

    report = "\n\n".join(parts) + "\n"
    report = truncate(report, COMMENT_LIMIT)
    with open(args.output, "w", encoding="utf-8") as f:
        f.write(report)
    print("Wrote {} ({} characters)".format(args.output, len(report)))
    line = status_description(args.result, stage, unit_fail, unit_failing,
                              ["ci-regression-summary.txt", "ci-parallel-regression-summary.txt"])
    with open("ci-status-description.txt", "w", encoding="utf-8") as f:
        f.write(line + "\n")
    print("Status description: " + line)


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
