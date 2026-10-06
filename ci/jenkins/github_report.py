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
from contextlib import closing
import datetime
import glob
import html
import json
import os
import re
import socket
import sqlite3
import sys
import urllib.error
import urllib.parse
import urllib.request

API = "https://api.github.com"
DEFAULT_REPO = "gkeyllorg/gkeyll"
MARKER_FORMAT = '<!-- gkeyll-ci-report context="{}" -->'
COMMENT_LIMIT = 60000  # GitHub allows 65536 characters; keep headroom.


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
    r"|^make(\[\d+\])?: \*\*\*|[Ss]egmentation fault|Abort trap|Bus error|[Tt]imed? ?out"
    r"|\berror\b(?:\s*#?\d+)?\s*:|\bfatal\b|Traceback \(most recent call last\)|command not found|No such file or directory", re.I)


# Compiler and runtime diagnostics, including GCC/Clang, NVCC and Intel forms.
WARNING_LINE = re.compile(r"\bwarning\b(?:\s*#?\d+)?\s*[:\[]|\bWARNING\b", re.I)
ANSI = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")


def log_paths():
    patterns = ('candidate-*.log', 'baseline-*.log', 'ci-command-logs/*.log',
                'slurm-*.out', 'build/**/*.log', 'cuda-build/**/*.log',
                '_baseline/build/**/*.log', '_baseline/cuda-build/**/*.log')
    return sorted({p for pattern in patterns for p in glob.glob(pattern, recursive=True)
                   if os.path.isfile(p) and not os.path.islink(p)})


def captured_logs():
    for path in log_paths():
        yield path, os.path.basename(path), read_text(path), False
    # runregression stores individual compiler/runtime output in SQLite rather
    # than printing it to the command log. Open these artifacts read-only.
    for root, label in (('gkylsoft', 'candidate'), ('_baseline/gkylsoft', 'baseline')):
        for path in sorted(glob.glob(root + '/gkeyll-results/**/regressiondb', recursive=True)):
            try:
                uri = 'file:' + urllib.parse.quote(os.path.abspath(path)) + '?mode=ro'
                with closing(sqlite3.connect(uri, uri=True)) as db:
                    rows = db.execute('SELECT name, runlog, status FROM RegressionData').fetchall()
                for name, output, status in rows:
                    step = '-regression/' + os.path.relpath(path, root) + '/' + str(name)
                    yield path + ':' + str(name), label + step, output or '', status in (-2, 1)
            except sqlite3.Error as err:
                yield path, label + '-regression-output', 'ERROR: Could not read regression diagnostics: ' + str(err), False


def warning_key(line):
    """Ignore checkout roots and source line shifts, but retain file and message."""
    line = ANSI.sub('', line)
    root = os.environ.get('WORKSPACE', os.getcwd()).rstrip('/')
    line = line.replace(root + '/_baseline/', '').replace(root + '/', '')
    line = re.sub(r'(?<!\w)(?:\.\./|\./|_baseline/)+', '', line)
    line = re.sub(r':\d+(?::\d+)?(?=:)', '', line)
    line = re.sub(r'\(\d+(?:,\d+)?\)(?=\s*:)', '', line)
    return line.strip()


def detail_sections(title, text):
    """Split inside the log, keeping every comment's HTML and fences balanced."""
    chunks = []
    while text:
        end = min(len(text), 48000)
        if end < len(text):
            end = text.rfind('\n', 0, end) + 1 or end
        chunks.append(text[:end])
        text = text[end:]
    return ['<details><summary>{}{}</summary>\n\n{}\n\n</details>'.format(
        html.escape(title), ' (part {})'.format(i + 1) if len(chunks) > 1 else '',
        fenced(chunk)) for i, chunk in enumerate(chunks or ['None.'])]


def diagnostic_sections():
    warnings, errors = [], []
    candidate, baseline = {}, set()
    comparable = set()
    paths = []
    for path, name, output, complete in captured_logs():
        paths.append(path)
        is_baseline = name.startswith('baseline-') or path.startswith('_baseline/')
        # Only compare equivalent completed commands; failed/missing steps are unknown.
        if name.startswith('baseline-'):
            complete = complete or (read_text(path + '.exit').strip() == '0' or
                        (name.endswith('.log') and os.path.isfile(name[:-4] + '-seconds.txt')))
            if complete:
                comparable.add(name[len('baseline-'):])
        lines = ANSI.sub('', output).splitlines()
        for i, line in enumerate(lines):
            entry = '{}:{}: {}'.format(path, i + 1, line)
            # Include source/caret/continuation lines with the diagnostic.
            context = []
            for following in lines[i + 1:i + 7]:
                if not following.strip() or not following[:1].isspace() or WARNING_LINE.search(following) or ERROR_LINE.search(following):
                    break
                context.append(following)
            if context:
                entry += '\n' + '\n'.join(context)
            if WARNING_LINE.search(line):
                warnings.append(entry)
                key = warning_key(line)
                if is_baseline:
                    baseline.add(key)
                elif name.startswith('candidate-'):
                    candidate[(name[len('candidate-'):], key)] = entry
            if ERROR_LINE.search(line) and not WARNING_LINE.search(line):
                errors.append(entry)

    selection = read_kv('ci-selection.txt')
    ref = first(selection, 'baseline_selector', 'main')
    new = [entry for (step, key), entry in candidate.items() if step in comparable and key not in baseline]
    unknown = [entry for (step, key), entry in candidate.items() if step not in comparable]
    note = 'Compared equivalent completed build steps against {}. Source line/column shifts are ignored.'.format(ref)
    if ref != 'main':
        note += ' This run selected a baseline other than main.'
    if not comparable:
        note = 'Comparison unavailable: no matching completed baseline build logs. Warnings cannot be classified as new.'
    elif unknown:
        note += '\n{} distinct candidate warnings could not be compared because their baseline step did not complete.'.format(len(unknown))
    sections = detail_sections('New warnings vs {} ({})'.format(ref, len(new) if comparable else 'unknown'),
                               note + '\n\n' + ('\n\n'.join(new) if new else 'No new warnings in comparable steps.' if comparable else ''))
    sections += detail_sections('All warnings ({})'.format(len(warnings)), '\n\n'.join(warnings) or 'No warnings found in captured logs.')
    sections += detail_sections('All errors ({})'.format(len(errors)), '\n\n'.join(errors) or 'No recognised error lines in captured logs; see failure details for command exits and infrastructure failures.')
    sections += detail_sections('Captured logs ({})'.format(len(paths)), '\n'.join(paths) or 'No command logs were produced before this run ended.')
    return sections


def report_pages(summary, sections, context):
    pages, current = [], summary
    for section in sections:
        if len(current) + len(section) + 2 > COMMENT_LIMIT:
            pages.append(current)
            current = MARKER_FORMAT.format(context + '/part-{}'.format(len(pages) + 1)) + '\n\n### CI diagnostics (continued)'
        current += '\n\n' + section
    pages.append(current)
    return pages


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


def failure_sections():
    """Show unfiltered tails first, keeping extracted warnings/errors separate."""
    detail = read_text('ci-failure-detail.txt').strip()
    match = re.search(r'^Full log: (.+?)(?: \(archived with the build\))?$', detail, re.M)
    command = re.search(r'^Command: (.*)$', detail, re.M)
    logs = {}
    if match and os.path.isfile(match.group(1)):
        logs[match.group(1)] = ['Command: ' + command.group(1)] if command else []
    for path in log_paths():
        exit_code = read_text(path + '.exit').strip()
        if exit_code and exit_code != '0':
            logs.setdefault(path, []).append('Exit code: ' + exit_code)
    # A Slurm launcher log describes scheduling, while the .out contains the
    # compiler/test output. Include both when a Slurm command failed.
    if not logs or any('shared-slurm-' in path for path in logs):
        outs = glob.glob('slurm-*.out')
        if outs:
            logs.setdefault(max(outs, key=os.path.getmtime), [])
    sections = []
    for path, metadata in logs.items():
        tail = '\n'.join(ANSI.sub('', read_text(path)).splitlines()[-100:])
        heading = metadata + ['Full log: ' + path, '', tail or '(log is empty)']
        sections += detail_sections('Failed build log — last 100 lines: ' + path,
                                    '\n'.join(heading))
    if not sections:
        sections += detail_sections('Failure details — build log unavailable',
                                    detail or 'No build log was captured before this failure; see the failed stage above.')
    return sections


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
        parts.extend(failure_sections())

    unit_fail, unit_failing = 0, []
    for label, results_file, log_file in (("Candidate", "candidate-unit-results.txt", "candidate-unit-test.log"),
                                          ("Baseline", "baseline-unit-results.txt", "baseline-unit-test.log")):
        section, nfail, failing = unit_section(label, results_file, log_file)
        if section:
            parts.append(section)
            unit_fail += nfail
            unit_failing += failing

    for title, path in (("Serial C regressions", "ci-regression-summary.txt"),
                        ("Parallel C regressions", "ci-parallel-regression-summary.txt")):
        section = regression_section(title, path)
        if section:
            parts.append(section)

    timings = timing_section()
    if timings:
        parts.append(timings)

    parts.append("_Full raw logs and regression databases stay on that machine; its owner can fetch them with "
                 "`./ci/jenkins/gkeyll-ci.sh {} artifact --build {} --fetch`._".format(args.platform, build_number))

    summary_sections = []
    for part in parts[1:]:
        summary_sections.extend([part] if len(part) <= 48000 else
                                detail_sections('Extended CI summary (Markdown)', part))
    pages = report_pages(parts[0], summary_sections + diagnostic_sections(), args.context)
    with open(args.output + '.json', 'w', encoding='utf-8') as f:
        json.dump(pages, f)
    report = '\n\n'.join(pages)
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


def publish_report(args):
    token = os.environ.get("GITHUB_TOKEN")
    if not token:
        sys.exit("GITHUB_TOKEN is not set")
    if not args.pr and not args.commit and args.ref:
        commit, _ = api('GET', '{}/repos/{}/commits/{}'.format(
            API, args.repo, urllib.parse.quote(args.ref, safe='')), token)
        args.commit = commit.get('sha', '')
    if not args.pr and not re.fullmatch(r"[0-9a-fA-F]{40}", args.commit):
        sys.exit("--commit must be a full 40-character SHA")
    body = read_text(args.report)
    marker = MARKER_FORMAT.format(args.context)
    if marker not in body:
        sys.exit("report {} does not carry the marker for context {}".format(args.report, args.context))
    repo = args.repo

    if args.pr:
        list_url = "{}/repos/{}/issues/{}/comments".format(API, repo, args.pr)
        update_url = "{}/repos/{}/issues/comments/{{}}".format(API, repo)
        target = "PR #{}".format(args.pr)
    else:
        list_url = "{}/repos/{}/commits/{}/comments".format(API, repo, args.commit)
        update_url = "{}/repos/{}/comments/{{}}".format(API, repo)
        target = "commit {}".format(short(args.commit))

    pages_path = args.report + '.json'
    pages = json.loads(read_text(pages_path)) if os.path.isfile(pages_path) else [body]
    if not isinstance(pages, list) or not pages:
        sys.exit('Report has no pages')
    markers = [MARKER_FORMAT.format(args.context if i == 0 else args.context + '/part-{}'.format(i + 1))
               for i in range(len(pages))]
    for page, page_marker in zip(pages, markers):
        if not isinstance(page, str) or page_marker not in page or len(page) > COMMENT_LIMIT:
            sys.exit('Invalid or oversized report page')
    # List once, including older overflow pages so shorter later runs cannot
    # leave stale errors and warnings on the PR.
    comments = []
    url = list_url + '?per_page=100'
    while url:
        batch, link = api('GET', url, token)
        comments.extend(batch or [])
        url = next_link(link)
    for page, page_marker in zip(pages, markers):
        existing = next((c['id'] for c in comments if page_marker in (c.get('body') or '')), None)
        result, _ = api('PATCH' if existing else 'POST',
                        update_url.format(existing) if existing else list_url, token, {'body': page})
        print('Published CI report on {}: {}'.format(target, (result or {}).get('html_url', '')))
    prefix = MARKER_FORMAT.format(args.context + '/part-').split('" -->')[0]
    for comment in comments:
        old = comment.get('body') or ''
        if prefix in old and not any(m in old for m in markers):
            old_marker = re.search(re.escape(prefix) + r'\d+" -->', old)
            if old_marker:
                api('PATCH', update_url.format(comment['id']), token,
                    {'body': old_marker.group() + '\n\nDiagnostics for this part were cleared by the latest CI run; see the main report.'})



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
    publish.add_argument("--commit", default="", help="candidate commit SHA")
    publish.add_argument("--ref", default="", help="resolve this ref if checkout failed before recording a SHA")
    publish.add_argument("--context", required=True)
    publish.add_argument("--pr", default="", help="pull-request number; omitted for commit comments")
    publish.add_argument("--repo", default=DEFAULT_REPO)
    publish.set_defaults(func=publish_report)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
