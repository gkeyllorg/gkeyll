#!/usr/bin/env python3
"""Build Jenkins reports and publish GitHub statuses and comments on tested commits.

Standard library only; uses the existing platform username/token credential.

    github_report.py build --platform NAME --context CTX --result STATE [--pr N]
    github_report.py update --platform NAME --context CTX --result STATE --commit SHA [--pr N]
    github_report.py publish --report ci-report.md --context CTX --commit SHA [--pr N]

`build` reads workspace artifacts without credentials. `update` publishes stage
statuses or the completed report, recording delivery in ci-report-delivery.json.
`publish` supports older jobs that publish their commit status separately.
The Pipelines load this file from reviewed, trusted CI code before checkout.
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
import time
import urllib.error
import urllib.parse
import urllib.request

API = "https://api.github.com"
DEFAULT_REPO = "gkeyllorg/gkeyll"
MARKER_FORMAT = '<!-- gkeyll-ci-report context="{}" -->'
COMMENT_LIMIT = 60000  # UTF-8 bytes; leave room for run metadata and navigation.
RUN_MARKER = '<!-- gkeyll-ci-run {} -->'


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


def dropdown(title, body):
    return '<details><summary>{}</summary>\n\n{}\n\n</details>'.format(
        html.escape(title), body)


def wrapped_log(text):
    # GitHub's inline code preserves spaces and wraps to the available width.
    # Put breaks outside code (GitHub hides br inside code), and keep this HTML
    # block free of blank source lines so log contents cannot become Markdown.
    return '<div>\n' + '<br>\n'.join(
        '<code>{}</code>'.format(html.escape(line, quote=False))
        for line in text.split('\n')) + '\n</div>'


# ---- build -----------------------------------------------------------------

STATE_LABEL = {
    "success": ":white_check_mark: passed",
    "failure": ":x: failed",
    "error": ":warning: errored or aborted",
    "pending": ":hourglass: running",
    "cancelled": ":stop_sign: cancelled",
    "timed_out": ":alarm_clock: timed out",
}


def regression_section(title, summary_file):
    values = read_kv(summary_file)
    if not values:
        return ""
    passed = first(values, "c_regression_passed", "?")
    candidate_only = first(values, "c_regression_candidate_only", "0")
    acked = first(values, "c_regression_acknowledged", "0")
    unacked = first(values, "c_regression_unacknowledged", "0")
    head = "{}: {} passed, {} candidate-only, {} acknowledged, {} unacknowledged".format(
        title, passed, candidate_only, acked, unacked)
    lines = []
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
                lines += ["", dropdown(name + ': files that differ from the baseline',
                                       wrapped_log("\n".join(details[name])))]
    acked_tests = values.get("c_regression_acknowledged_test", [])
    if acked_tests:
        lines += ["", "Acknowledged diffs (new or updated versus the baseline): "
                  + ", ".join(code(t) for t in acked_tests)]
    candidate_only_tests = values.get("c_regression_candidate_only_test", [])
    if candidate_only_tests:
        lines += ["", "Candidate-only tests (executed but not compared to a baseline): "
                  + ", ".join(code(t) for t in candidate_only_tests)]
    return dropdown(head, "\n".join(lines).strip() or 'No per-test details recorded.')


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
    # Baseline logs survive across runs and retain their original absolute root.
    line = re.sub(r'(?<!\w)/[^\s:]*?/(?:baseline-cache/[^/]+/[0-9a-f]{40}'
                  r'|runs/[^/]+/[^/]+/(?:candidate|baseline)/[0-9a-f]{40})/gkeyll/', '', line)
    line = re.sub(r'(?<!\w)(?:\.\./|\./|_baseline/)+', '', line)
    line = re.sub(r':\d+(?::\d+)?(?=:)', '', line)
    line = re.sub(r'\(\d+(?:,\d+)?\)(?=\s*:)', '', line)
    return line.strip()


def byte_prefix(text, limit):
    return text.encode('utf-8')[:limit].decode('utf-8', errors='ignore')


def detail_sections(title, text):
    """Split inside the log, keeping rendered HTML within the comment limit."""
    chunks = []
    while text:
        end = len(byte_prefix(text, 48000))
        while True:
            if end < len(text):
                end = text.rfind('\n', 0, end) + 1 or end
            chunk = wrapped_log(text[:end])
            if len(chunk.encode('utf-8')) <= 48000:
                break
            # Escaping and per-line tags can expand the log substantially.
            end //= 2
        chunks.append(chunk)
        text = text[end:]
    return ['<details><summary>{}{}</summary>\n\n{}\n\n</details>'.format(
        html.escape(title), ' (part {})'.format(i + 1) if len(chunks) > 1 else '',
        chunk) for i, chunk in enumerate(chunks or [wrapped_log('None.')])]


def diagnostic_sections(summary_output=None):
    warnings, errors = [], []
    candidate, baseline = {}, {}
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
                    step = name[len('baseline-'):] if name.startswith('baseline-') else name
                    baseline.setdefault(step, set()).add(key)
                elif name.startswith('candidate-'):
                    candidate[(name[len('candidate-'):], key)] = entry
            if ERROR_LINE.search(line) and not WARNING_LINE.search(line):
                errors.append(entry)

    selection = read_kv('ci-selection.txt')
    ref = first(selection, 'baseline_selector', 'main')
    new = [entry for (step, key), entry in candidate.items()
           if step in comparable and key not in baseline.get(step, set())]
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
    if summary_output:
        with open(summary_output, 'w', encoding='utf-8') as stream:
            json.dump({'warnings': len(warnings), 'new_warnings': len(new) if comparable else None,
                       'unclassified_warnings': len(unknown), 'errors': len(errors),
                       'captured_logs': len(paths), 'captured_log_paths': paths}, stream)
    return sections


def report_pages(summary, sections, context):
    pages, current = [], summary
    for section in sections:
        if len((current + '\n\n' + section).encode('utf-8')) > COMMENT_LIMIT:
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
    body = table
    for layer in order:
        for test in layers.get(layer, {}).get("failed", []):
            detail = assertion_lines(test, log_file)
            if detail:
                body += ["", dropdown(test + ': failed checks', wrapped_log("\n".join(detail)))]
    return dropdown(head, "\n".join(body)), nfail, failing


def status_description(result, stage, unit_fail, unit_failing, regression_files):
    """One line for the GitHub status: name what failed, or what passed."""
    if result == "success":
        return "Passed: unit tests and C regressions."
    if result == "pending":
        return "Running: {}.".format(stage)
    if result in ("error", "cancelled", "timed_out"):
        return "{} at stage: {}.".format(
            {'error': 'Errored', 'cancelled': 'Cancelled', 'timed_out': 'Timed out'}[result], stage)
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

    def log_metadata(path):
        # Failure records can use absolute paths while discovery uses relative
        # paths. Resolve aliases before merging metadata and rendering one tail.
        return logs.setdefault(os.path.relpath(os.path.realpath(path)), [])

    if match and os.path.isfile(match.group(1)):
        metadata = log_metadata(match.group(1))
        if command:
            metadata.append('Command: ' + command.group(1))
    for path in log_paths():
        exit_code = read_text(path + '.exit').strip()
        if exit_code and exit_code != '0':
            log_metadata(path).append('Exit code: ' + exit_code)
    # A Slurm launcher log describes scheduling, while the .out contains the
    # compiler/test output. Include both when a Slurm command failed.
    if not logs or any('shared-slurm-' in path for path in logs):
        outs = glob.glob('slurm-*.out')
        if outs:
            log_metadata(max(outs, key=os.path.getmtime))
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


def timing_section(elapsed=None):
    values = read_kv("ci-timing-summary.txt")
    rows = [(k, first(values, k)) for k in values if k.endswith("_seconds")]
    rows = [(k, v) for k, v in rows if v and v != "not-recorded"]
    if not rows and elapsed is None:
        return ""
    body = ["| Step | Seconds |", "| --- | ---: |"]
    body += ["| {} | {} |".format(k[:-len("_seconds")].replace("_", " "), v) for k, v in rows]
    if elapsed is not None:
        body.append("| **Total elapsed** | **{}** |".format(elapsed))
    return dropdown('Timings', "\n".join(body))


def run_timing():
    """Use controller timestamps so agent clock differences do not affect elapsed time."""
    end_ms = int(os.environ.get('CI_REPORT_END_MS') or time.time() * 1000)
    finished = datetime.datetime.fromtimestamp(end_ms / 1000, datetime.timezone.utc)
    try:
        start_ms = int(os.environ.get('CI_REPORT_START_MS', ''))
        if not 0 <= start_ms <= end_ms:
            return None, finished, None
        started = datetime.datetime.fromtimestamp(start_ms / 1000, datetime.timezone.utc)
    except (ValueError, OverflowError, OSError):
        return None, finished, None
    return started, finished, (end_ms - start_ms) // 1000


def build_report(args):
    selection = read_kv("ci-selection.txt")
    preflight = read_kv("ci-baseline-preflight.txt")
    cache = read_kv("ci-baseline-cache.txt")
    candidate = (read_text("ci-candidate-commit.txt").strip() or first(preflight, "candidate_commit")
                 or os.environ.get("CI_REPORT_COMMIT", ""))
    baseline = read_text("ci-baseline-commit.txt").strip() or first(preflight, "baseline_commit")
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
    cache_label = {'hit': 'loaded from cache', 'saved': 'saved to cache',
                   'miss': 'cache miss; not saved'}.get(first(cache, 'status'), '')
    if first(cache, 'baseline_commit') != baseline:
        cache_label = ''
    meta = ["**Candidate:** " + describe(candidate_sel, candidate),
            "**Baseline:** " + describe(baseline_sel, baseline)
            + (" ({})".format(cache_label) if cache_label else '')]
    pipeline_commit = (read_text('ci-trusted-ci-commit.txt').strip()
                       or os.environ.get('CI_TRUSTED_CI_COMMIT', ''))
    pipeline_ref = os.environ.get('CI_TRUSTED_CI_REF', '')
    reporting_commit = read_text('ci-reporting-commit.txt').strip()
    checker_commit = read_text('ci-trusted-checker-commit.txt').strip()
    for label, commit in [('Jenkinsfile', pipeline_commit), ('Reporting tools', reporting_commit),
                          ('Regression checker', checker_commit)]:
        if re.fullmatch(r'[0-9a-fA-F]{40}', commit):
            source = '[{}](https://github.com/{}/commit/{})'.format(code(commit), DEFAULT_REPO, commit)
            if label == 'Jenkinsfile' and pipeline_ref:
                source = code(pipeline_ref) + ' @ ' + source
            meta.append('**{}:** {}'.format(label, source))
    # No controller URL: every controller is loopback-only, so a link would be
    # dead for everyone but the machine owner. The build number is what that
    # owner needs to fetch artifacts (see the footer).
    started, finished, elapsed = run_timing()
    timestamp_format = "%Y-%m-%d %H:%M:%S UTC"
    meta.append("**Run:** {} build #{} on {}".format(code(args.platform), build_number, code(node)))
    meta.append("**Start:** {} · **{}:** {} · **Elapsed:** {}".format(
        started.strftime(timestamp_format) if started else 'not recorded',
        'Updated' if args.result == 'pending' else 'End', finished.strftime(timestamp_format),
        '{} ({} s)'.format(datetime.timedelta(seconds=elapsed), elapsed) if elapsed is not None else 'not recorded'))
    queued_ms = os.environ.get('CI_QUEUE_ENQUEUED_MS', '')
    if queued_ms.isdigit() and started:
        queue_seconds = max(0, (int(started.timestamp() * 1000) - int(queued_ms)) // 1000)
        meta.append('**Queue wait:** {} s (excluded from execution elapsed)'.format(queue_seconds))
    retry_ms = os.environ.get('CI_BOOTSTRAP_RETRY_WAIT_MS', '')
    if retry_ms.isdigit():
        meta.append('**Setup retry waiting:** {} s (included in execution elapsed)'.format(int(retry_ms) // 1000))
    if os.environ.get('CI_QUEUE_ID'):
        meta.append('**Jenkins queue:** #' + os.environ['CI_QUEUE_ID'])
    meta.insert(0, '**Status context:** ' + code(args.context))
    if preflight:
        behind = first(preflight, "behind_by", "?")
        override = first(preflight, "override", "false")
        state = "override enabled" if override == "true" else "verified"
        meta.append("**Baseline preflight:** {} ({} commit(s) behind)".format(state, behind))
    parts.append(dropdown('Run details', "  \n".join(meta)))

    stage = first(failure, "stage", os.environ.get('CI_FAILURE_STAGE', 'unknown'))
    if args.result == 'pending':
        parts.append(dropdown('Current stage: ' + stage, '**Current stage:** ' + stage))
    elif args.result != "success":
        message = re.sub(r"^(?:[A-Za-z_$][\w$]*\.)+[A-Z]\w*(?:Exception|Error): ", "", first(failure, "message", ""))
        parts.append(dropdown('Failed at stage: ' + stage,
                              "**Failed at stage:** {}{}".format(stage, " — " + message if message else "")))
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

    stage_timings = read_text('ci-stage-timings.json')
    if stage_timings:
        rows = json.loads(stage_timings).get('stages', [])
        parts.append(dropdown('Stage durations', '| Stage | Duration | Result |\n| --- | ---: | --- |\n' +
                     '\n'.join('| {} | {:.3f} s | {} |'.format(
                         code(row['stage']), row['elapsed_ms'] / 1000, row['result']) for row in rows)))
    timings = timing_section(elapsed)
    if timings:
        parts.append(timings)

    parts.append(dropdown('Fetch full logs and regression databases',
                 "_Full raw logs and regression databases stay on that machine; its owner can fetch them with "
                 "`./ci/jenkins/gkeyll-ci.sh {} artifact --build {} --fetch`._".format(args.platform, build_number)))

    progress = read_text('ci-stage-history.txt').strip()
    if progress:
        parts += detail_sections('Stage history (UTC)', progress)
    for path in sorted(glob.glob('slurm-*-job-status.txt')):
        parts += detail_sections('Slurm allocation: ' + path, read_text(path))

    summary_sections = []
    for part in parts[1:]:
        summary_sections.extend([part] if len(part.encode('utf-8')) <= 48000 else
                                detail_sections('Extended CI summary (Markdown)', part))
    pages = report_pages(parts[0], summary_sections + diagnostic_sections('ci-diagnostic-summary.json'), args.context)
    with open(args.output + '.json', 'w', encoding='utf-8') as f:
        json.dump(pages, f)
    report = '\n\n'.join(pages)
    with open(args.output, "w", encoding="utf-8") as f:
        f.write(report)
    print("Wrote {} ({} characters)".format(args.output, len(report)))
    line = status_description(args.result, stage, unit_fail, unit_failing,
                              ["ci-regression-summary.txt", "ci-parallel-regression-summary.txt"])
    if elapsed is not None and args.result != 'pending':
        outcome, separator, detail = line.partition(':')
        outcome = re.sub(r' at stage$', '', outcome.rstrip('.'))
        line = '{} after {} s{}{}'.format(outcome, elapsed, separator, detail)
    diagnostics = json.loads(read_text('ci-diagnostic-summary.json'))
    if diagnostics['warnings']:
        if diagnostics['new_warnings'] is None:
            line += ' Warning comparison unavailable.'
        elif diagnostics['new_warnings']:
            line += ' {} new warning(s).'.format(diagnostics['new_warnings'])
    with open("ci-status-description.txt", "w", encoding="utf-8") as f:
        f.write(line + "\n")
    print("Status description: " + line)


# ---- publish ---------------------------------------------------------------

class GitHubError(RuntimeError):
    def __init__(self, method, status, retryable=False):
        super().__init__('GitHub {} failed: HTTP {}'.format(method, status))
        self.retryable = retryable


def api(method, url, token, payload=None):
    # Pagination links must never redirect the credential to another host.
    parsed = urllib.parse.urlsplit(url)
    if parsed.scheme != 'https' or parsed.netloc != 'api.github.com':
        raise ValueError('Unexpected GitHub API URL')
    data = json.dumps(payload).encode() if payload is not None else None
    request = urllib.request.Request(url, data=data, method=method, headers={
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
        "Authorization": "Bearer " + token,
        "Content-Type": "application/json",
        "User-Agent": "gkeyll-ci-report",
    })
    try:
        with urllib.request.urlopen(request, timeout=25) as response:
            body = response.read().decode("utf-8", "replace")
            return (json.loads(body) if body else None), response.headers.get("Link", "")
    except urllib.error.HTTPError as err:
        # Response bodies can contain private data; archive only method/status.
        raise GitHubError(method, err.code, err.code in (408, 429, 500, 502, 503, 504)) from None
    except (OSError, TimeoutError):
        raise OSError('GitHub {} request failed (network or timeout)'.format(method)) from None


def retry(operation):
    # Retry the entire operation, including comment discovery. If a POST was
    # accepted before a connection failed, the next attempt updates that comment.
    for attempt in range(3):
        try:
            return operation()
        except (GitHubError, OSError) as err:
            if attempt == 2 or (isinstance(err, GitHubError) and not err.retryable):
                raise
            time.sleep(2 ** attempt)


def next_link(link_header):
    for part in link_header.split(","):
        match = re.search(r'<([^>]+)>;\s*rel="next"', part)
        if match:
            return match.group(1)
    return None


def run_metadata(args):
    return dict(commit=args.commit, queue_id=getattr(args, 'queue_id', '') or os.environ.get('CI_QUEUE_ID', ''),
                build_id=os.environ.get('BUILD_TAG', ''), job=os.environ.get('JOB_NAME', ''),
                build_number=os.environ.get('BUILD_NUMBER', ''), result=getattr(args, 'result', ''))


def comment_metadata(body):
    match = re.search(r'^<!-- gkeyll-ci-run (.+) -->$', body, re.M)
    try:
        value = json.loads(match.group(1)) if match else {}
        return value if isinstance(value, dict) else {}
    except ValueError:
        return {}


def superseded(previous, current):
    """Queue IDs increase on a controller; contexts must be unique per machine."""
    before, after = str(previous.get('queue_id', '')), str(current.get('queue_id', ''))
    if before.isdigit() and after.isdigit():
        return int(before) > int(after)
    before, after = str(previous.get('build_number', '')), str(current.get('build_number', ''))
    return (bool(current.get('job')) and previous.get('job') == current['job'] and
            before.isdigit() and after.isdigit() and int(before) > int(after))


def publish_report(args, token=None):
    token = token or os.environ.get("GITHUB_TOKEN")
    if not token:
        raise ValueError("GITHUB_TOKEN is not set")
    # The PR head or branch can move during a run. Only the recorded candidate
    # SHA identifies the commit whose results this report describes.
    if not re.fullmatch(r"[0-9a-fA-F]{40}", args.commit):
        raise ValueError("--commit must be a full 40-character SHA")
    body = read_text(args.report)
    marker = MARKER_FORMAT.format(args.context)
    if not body.startswith(marker + '\n'):
        raise ValueError('Report does not carry the expected context marker')
    base = '{}/repos/{}'.format(API, args.repo)
    list_url = base + '/commits/{}/comments'.format(args.commit)
    update_url = base + '/comments/{}'

    pages_path = args.report + '.json'
    pages = json.loads(read_text(pages_path)) if os.path.isfile(pages_path) else [body]
    if not isinstance(pages, list) or not pages:
        raise ValueError('Report has no pages')
    markers = [MARKER_FORMAT.format(args.context if i == 0 else args.context + '/part-{}'.format(i + 1))
               for i in range(len(pages))]
    for page, page_marker in zip(pages, markers):
        if not isinstance(page, str) or not page.startswith(page_marker + '\n') or len(page.encode('utf-8')) > COMMENT_LIMIT:
            raise ValueError('Invalid or oversized report page')
    # Never edit a contributor's comment just because it quotes a report marker.
    user, _ = api('GET', API + '/user', token)
    comments = []
    url = list_url + '?per_page=100'
    while url:
        batch, link = api('GET', url, token)
        comments.extend(c for c in batch or [] if c.get('user', {}).get('id') == user['id'])
        url = next_link(link)
    metadata = run_metadata(args)
    main = next((c for c in reversed(comments) if (c.get('body') or '').startswith(marker + '\n')), None)
    if main and superseded(comment_metadata(main['body']), metadata):
        print('A newer run owns the current report; this report remains in Jenkins artifacts.')
        return {'skipped': True, 'html_url': ''}

    published = {}
    # Publish continuations first so a completed main report never points to
    # missing diagnostics. Reuse those comments and clear obsolete pages.
    for index in reversed(range(len(pages))):
        page, page_marker = pages[index], markers[index]
        existing = next((c for c in reversed(comments) if (c.get('body') or '').startswith(page_marker + '\n')), None)
        content = page.replace(page_marker, page_marker + '\n' + RUN_MARKER.format(json.dumps(metadata)), 1)
        if index > 0 and index + 1 in published:
            content += '\n\n[Next diagnostics page]({})'.format(published[index + 1]['html_url'])
        if index == 0 and published:
            navigation = '\n\n**More diagnostics:** ' + ' · '.join(
                '[Part {}]({})'.format(i + 1, published[i]['html_url']) for i in sorted(published))
            if len((content + navigation).encode('utf-8')) > 65000:
                navigation = '\n\n[Continue diagnostics]({}) ({} more pages).'.format(
                    published[1]['html_url'], len(published))
            content += navigation
        if len(content.encode('utf-8')) > 65000:
            raise ValueError('Report page exceeds the publication limit after adding metadata')
        result, _ = api('PATCH' if existing else 'POST',
                        update_url.format(existing['id']) if existing else list_url, token, {'body': content})
        published[index] = result
    prefix = MARKER_FORMAT.format(args.context + '/part-').split('" -->')[0]
    for comment in comments:
        old = comment.get('body') or ''
        if old.startswith(prefix) and not any(old.startswith(m + '\n') for m in markers):
            old_marker = re.match(re.escape(prefix) + r'\d+" -->', old)
            if old_marker:
                api('PATCH', update_url.format(comment['id']), token, {'body': old_marker.group() +
                    '\n\nThis continuation is no longer needed; see the [current report]({}).'.format(published[0]['html_url'])})
    print('Published CI report: ' + published[0]['html_url'])
    return published[0]


def latest_status(args, token):
    url = '{}/repos/{}/commits/{}/status?per_page=100'.format(API, args.repo, args.commit)
    while url:
        response, link = api('GET', url, token)
        current = next((s for s in response['statuses'] if s['context'] == args.context), None)
        if current:
            return current
        url = next_link(link)
    return None


def publish_status(args, token, target_url='', report_failed=False):
    state = {'cancelled': 'error', 'timed_out': 'error'}.get(args.result, args.result)
    description = args.description or status_description(args.result, args.stage, 0, [], [])
    current = latest_status(args, token)
    previous_description = (current or {}).get('description') or ''
    # Retain compatibility with installed controller hooks and older Pipelines.
    owner = re.search(r' \(Jenkins (?:queue #|(?:build #\d+; )?ID )(\d+)\)\.$', previous_description)
    if owner and args.queue_id.isdigit():
        if int(owner.group(1)) > int(args.queue_id):
            return {'skipped': True}
        if owner.group(1) == args.queue_id and current['state'] != 'pending' and state == 'pending':
            return {'skipped': True}
    build_number = os.environ.get('BUILD_NUMBER', '')
    build = 'build #{}; '.format(build_number) if build_number.isdigit() else ''
    suffix = ' (Jenkins {}ID {}).'.format(build, args.queue_id) if args.queue_id else ''
    if args.result == 'pending':
        activity = (description[9:] if description.startswith('Running: ') else description).rstrip('.')
        # Keep the most recent Jenkins timing snapshot through stage changes.
        # The controller refreshes both elapsed time and ETA every minute.
        timing = None
        if owner and owner.group(1) == args.queue_id and current['state'] == 'pending':
            timing = re.search(r'(?:^Gkeyll CI running: |; )(elapsed .*; ETA .*?) \(Jenkins ', previous_description)
        timing = '; ' + timing.group(1) if timing else ''
        prefix = 'Gkeyll CI running: '
        description = prefix + activity[:max(0, 140 - len(prefix + timing + suffix))] + timing
    if report_failed:
        suffix = '; report unavailable' + suffix
    description = description[:140 - len(suffix)] + suffix
    payload = {'state': state, 'context': args.context, 'description': description}
    # Controllers are loopback-only. Always link to a GitHub page people can read.
    payload['target_url'] = target_url or (current or {}).get('target_url') or \
        'https://github.com/{}/commit/{}'.format(args.repo, args.commit)
    if current and all(current.get(k) == v for k, v in payload.items()):
        return current
    return api('POST', '{}/repos/{}/statuses/{}'.format(API, args.repo, args.commit), token, payload)[0]


def update_report(args):
    delivery = dict(commit=args.commit, context=args.context, queue_id=args.queue_id,
                    comment=False, status=False, errors=[])
    try:
        token = os.environ.get('GITHUB_TOKEN')
        if not token:
            raise ValueError('GITHUB_TOKEN is not set')
        if not re.fullmatch(r'[0-9a-f]{40}', args.commit):
            raise ValueError('Reporting requires the exact candidate commit SHA')
        target = ''
        # Stage progress uses statuses. One comment is updated at completion,
        # avoiding a stream of comment notifications for every command or heartbeat.
        if args.result != 'pending':
            try:
                comment = retry(lambda: publish_report(args, token))
                delivery['comment'] = not comment.get('skipped', False)
                delivery['comment_skipped'] = comment.get('skipped', False)
                target = comment.get('html_url', '')
                delivery['report_url'] = target
            except (RuntimeError, OSError, ValueError) as err:
                delivery['errors'].append('Comment: ' + str(err))
        try:
            status = retry(lambda: publish_status(args, token, target, bool(delivery['errors'])))
            delivery['status'] = not status.get('skipped', False)
            delivery['status_skipped'] = status.get('skipped', False)
        except (RuntimeError, OSError, ValueError) as err:
            delivery['errors'].append('Commit status: ' + str(err))
    except (RuntimeError, OSError, ValueError) as err:
        delivery['errors'].append(str(err))
    finally:
        with open('ci-report-delivery.json', 'w', encoding='utf-8') as stream:
            json.dump(delivery, stream, indent=2)
    if delivery['errors']:
        raise RuntimeError('CI reporting incomplete: ' + '; '.join(delivery['errors']) +
                           '; see ci-report-delivery.json')


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

    publish = sub.add_parser("publish", help="create or update the GitHub commit comment carrying the report")
    publish.add_argument("--report", default="ci-report.md")
    publish.add_argument("--commit", required=True, help="exact tested commit SHA")
    publish.add_argument("--ref", default="", help="candidate ref (informational only; --commit selects the target)")
    publish.add_argument("--context", required=True)
    publish.add_argument("--pr", default="", help="pull-request number (informational only; --commit selects the target)")
    publish.add_argument("--repo", default=DEFAULT_REPO)
    publish.set_defaults(func=publish_report)

    update = sub.add_parser('update', help='publish commit status and the completed report on the tested commit')
    update.add_argument('--report', default='ci-report.md')
    update.add_argument('--commit', required=True)
    update.add_argument('--context', required=True)
    update.add_argument('--platform', required=True)
    update.add_argument('--result', required=True, choices=sorted(STATE_LABEL))
    update.add_argument('--description', default='')
    update.add_argument('--stage', default=os.environ.get('CI_FAILURE_STAGE', 'starting'))
    update.add_argument('--pr', default='', help='pull-request number (informational only; --commit selects the target)')
    update.add_argument('--ref', default='', help='candidate ref (informational only; --commit selects the target)')
    update.add_argument('--repo', default=DEFAULT_REPO)
    update.add_argument('--queue-id', default=os.environ.get('CI_QUEUE_ID', ''))
    update.set_defaults(func=update_report)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except (RuntimeError, OSError, ValueError) as err:
        sys.exit('ERROR: ' + str(err))
