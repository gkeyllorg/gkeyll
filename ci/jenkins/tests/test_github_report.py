"""Offline checks for report timing and log presentation."""
import argparse
import contextlib
from html.parser import HTMLParser
import io
import json
import os
from pathlib import Path
import re
import sqlite3
import tempfile
import unittest
from unittest.mock import patch

from ci.jenkins import github_report as report


def log_text(markup):
    """Read displayed log text, ignoring report headings and HTML layout."""
    class LogParser(HTMLParser):
        in_code = False

        def handle_starttag(self, tag, attrs):
            if tag == 'code':
                self.in_code = True
            elif tag == 'br':
                parts.append('\n')

        def handle_endtag(self, tag):
            if tag == 'code':
                self.in_code = False

        def handle_data(self, data):
            if self.in_code:
                parts.append(data)

    parts = []
    LogParser().feed(markup)
    return ''.join(parts)


class FailureReportTests(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory()
        self.previous = os.getcwd()
        os.chdir(self.workspace.name)

    def tearDown(self):
        os.chdir(self.previous)
        self.workspace.cleanup()

    def write(self, path, text):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(text)

    def test_regression_frames_count_unique_nonzero_frames_without_field(self):
        root = 'gkylsoft/gkeyll-results/'
        self.assertEqual(report.regression_frame_sections(), [])
        for filename in ('field_0.gkyl', 'field_1.gkyl', 'elc_0.gkyl', 'elc_1.gkyl',
                         'ion_1.gkyl', 'integrated_mom.gkyl'):
            self.write(root + 'vlasov/creg-runs/cpu_serial/rt_single/' + filename, '')
        self.assertEqual(report.regression_frame_sections(), [])
        for filename in ('elc_0.gkyl', 'elc_1.gkyl', 'elc_003.gkyl', 'ion_3.gkyl',
                         'ion_8.gkyl', 'integrated_mom.gkyl', 'notes_9.txt'):
            self.write(root + 'vlasov/creg-runs/cpu_serial/rt_many/' + filename, '')
        # Accepted/baseline outputs, subdirectories, and symlinks are not run frames.
        for prefix in ('vlasov/creg-accepted/rt_single/',
                       'vlasov/creg-runs/cpu_serial/rt_single/_cpu_output/'):
            for frame in (2, 3):
                self.write(root + prefix + 'field_{}.gkyl'.format(frame), '')
        for frame in (2, 3):
            self.write('_baseline/' + root + 'vlasov/creg-runs/cpu_serial/rt_base/field_{}.gkyl'.format(frame), '')
        Path(root + 'vlasov/creg-runs/cpu_serial/rt_single/field_4.gkyl').symlink_to('field_1.gkyl')
        body = '\n'.join(report.regression_frame_sections())
        self.assertIn('| `vlasov/creg/rt_many` | `vlasov` | `cpu_serial` | 3 |', body)
        self.assertNotIn('rt_single', body)
        self.assertNotIn('rt_base', body)

    def test_regression_frames_from_staged_inventory_in_completed_report(self):
        self.write('gkylsoft/gkeyll-results/ci-output-files.txt', '\0'.join([
            './parallel-c-4/moments/creg-runs/gpu_parallel/rt_many/field_{}.gkyl'.format(frame)
            for frame in (0, 1, 2)
        ] + [
            './vlasov/luareg-runs/cpu_serial/no_field/elc-{}.gkyl'.format(frame)
            for frame in (0, 1, 2, 3)
        ]) + '\0')
        args = argparse.Namespace(platform='personal', context='ci/test', result='success', pr='', output='report.md')
        with contextlib.redirect_stdout(io.StringIO()):
            report.build_report(args)
        body = Path('report.md').read_text()
        self.assertIn('<summary>Regression tests outputting more than 1 frame</summary>', body)
        self.assertIn('| `moments/creg/rt_many` | `parallel-c-4/moments` | `gpu_parallel` | 2 |', body)
        self.assertIn('| `vlasov/luareg/no_field` | `vlasov` | `cpu_serial` | 3 |', body)
        self.write('gkylsoft/gkeyll-results/ci-output-files.txt', '')
        self.assertEqual(report.regression_frame_sections(), [])

    def test_expected_diff_entries_for_passing_tests(self):
        self.write('summary.txt', 'c_regression_passed=2\n'
                   'c_regression_listed_passing=moments/creg/rt_pass|inherited|rt_pass # old | reason\n'
                   'c_regression_listed_passing=moments/creg/rt_pass|new or updated|rt_pass cpu_serial out.gkyl\n')
        text = '\n'.join(report.regression_sections('Serial C regressions', 'summary.txt')[0])
        self.assertIn('2 expected-diff entries with passing tests', text)
        self.assertIn('Listed in expected_regression_diffs.txt but passing', text)
        self.assertIn('inherited', text)
        self.assertIn('new or updated', text)
        self.assertIn(r'old \| reason', text)
        self.assertIn('This notice does not fail CI.', text)
        self.write('summary.txt', 'c_regression_passed=2\n')
        self.assertNotIn('but passing', '\n'.join(report.regression_sections('Serial', 'summary.txt')[0]))

    def regression_fixture(self, files=15):
        lines = ['c_regression_passed=3', 'c_regression_acknowledged=2',
                 'c_regression_unacknowledged=2', 'c_regression_layer=moments:3:2:2',
                 'c_regression_acknowledged_test=moments/creg/rt_ack',
                 'c_regression_acknowledged_status=moments/creg/rt_ack:fail',
                 'c_regression_acknowledged_test=moments/creg/rt_timeout',
                 'c_regression_acknowledged_status=moments/creg/rt_timeout:timeout',
                 'c_regression_failure=moments/creg/rt_diff:fail',
                 'c_regression_failure=moments/creg/rt_crash:crash']
        for kind, name in [('acknowledged', 'ack'), ('failure', 'diff')]:
            for index in range(files):
                lines.append('c_regression_{}_detail=moments/creg/rt_{}|'
                             '{}_{}.gkyl  [DIFF]  max_abs=1 max_rel=2 <&>'.format(
                                 kind, name, name, index))
            lines.append('c_regression_{}_detail=moments/creg/rt_{}|'
                         '{}_missing.gkyl  [MISSING]'.format(kind, name, name))
        self.write('ci-regression-summary.txt', '\n'.join(lines))

    def test_regressions_have_two_groups_with_all_files_and_instructions(self):
        self.regression_fixture()
        sections, full_lists = report.regression_sections('Serial C regressions', 'ci-regression-summary.txt')
        text = '\n'.join(sections)
        self.assertEqual(text.count('<details>'), 3)
        self.assertIn('| moments | 3 | 2 | 2 :x: |', text)
        self.assertIn('Acknowledged tests (2)', sections[1])
        self.assertIn('Unacknowledged tests (2)', sections[2])
        for section, name, failure, status in [(sections[1], 'ack', 'timeout', 'timeout'),
                                              (sections[2], 'diff', 'crash', 'crash')]:
            displayed = log_text(section)
            self.assertIn('rt_' + failure, displayed)
            self.assertIn(status, displayed)
            self.assertIn('Comparison reasons:', section)
            self.assertIn('15 differing files, 1 missing file.', section)
            self.assertIn('max_abs=1 ({}_0.gkyl)'.format(name), displayed)
            self.assertIn('missing output', displayed)
            self.assertLess(section.index('rt_' + failure), section.index('rt_' + name))
        for name in ('ack', 'diff'):
            for index in range(15):
                self.assertIn('{}_{}.gkyl  [DIFF]  max_abs=1 max_rel=2 <&>'.format(name, index),
                              log_text('\n'.join(full_lists)))
        self.assertNotIn('diff_0.gkyl', sections[1])
        self.assertNotIn('ack_0.gkyl', sections[2])
        self.assertIn('ci/jenkins/expected_regression_diffs.txt', sections[2])
        self.assertIn('moments/creg/rt_diff   # <why this changed>', log_text(sections[2]))
        self.assertIn('moments/creg/rt_crash   # <why this changed>', log_text(sections[2]))
        self.assertNotIn('rt_ack   #', log_text(text))
        self.assertIn('Only entries new or changed relative to the baseline', text)

    def test_regression_groups_with_empty_or_legacy_summaries(self):
        self.assertEqual(report.regression_sections('Serial', 'absent'), ([], []))
        self.write('summary.txt', 'c_regression_passed=1\n')
        sections, full_lists = report.regression_sections('Serial', 'summary.txt')
        self.assertEqual(full_lists, [])
        self.assertEqual(len(sections), 3)
        for section in sections[1:]:
            self.assertIn('None.', section)
        self.write('summary.txt', 'c_regression_acknowledged=1\n'
                   'c_regression_acknowledged_test=moments/creg/rt_legacy\n')
        text = '\n'.join(report.regression_sections('Serial', 'summary.txt')[0])
        self.assertIn('rt_legacy', text)
        self.assertIn('status unavailable', text)
        self.assertIn('No per-file comparison details recorded.', text)

    def test_large_regression_report_keeps_all_files_in_artifact(self):
        self.regression_fixture(files=1000)
        args = argparse.Namespace(platform='personal', context='ci/test', result='failure', pr='', output='ci-report.md')
        with contextlib.redirect_stdout(io.StringIO()):
            report.build_report(args)
        pages = json.loads(Path(args.output + '.json').read_text())
        self.assertEqual(len(pages), 1)
        for page in pages:
            self.assertLessEqual(len(page.encode('utf-8')), report.COMMENT_LIMIT)
            for tag in ('details', 'div', 'code'):
                self.assertEqual(page.count('<' + tag + '>'), page.count('</' + tag + '>'))
        text = '\n'.join(pages)
        self.assertEqual(pages[0].count('<details><summary>Serial C regressions:'), 3)
        self.assertIn('1000 differing files', pages[0])
        self.assertNotIn('ack_500.gkyl', pages[0])
        self.assertNotIn('full file lists</summary>', text)
        self.assertNotIn('Extended CI summary (Markdown)', text)
        self.assertNotIn('Fetch full logs', text)
        self.assertNotIn('artifact --build', text)
        displayed = log_text(Path('ci-report-regressions.md').read_text())
        for name in ('ack', 'diff'):
            for index in range(1000):
                self.assertEqual(displayed.count('{}_{}.gkyl  [DIFF]'.format(name, index)), 1)
        self.assertIn('moments/creg/rt_diff   # <why this changed>', log_text(text))

    def test_many_accepted_diffs_fit_one_comment(self):
        lines = ['c_regression_passed=281', 'c_regression_acknowledged=37',
                 'c_regression_unacknowledged=2']
        for index in range(39):
            name = 'gyrokinetic/creg/rt_gk_helical_adapt_src_sheath_surr_3x2v_p1_{}'.format(index)
            if index < 37:
                lines += ['c_regression_acknowledged_test=' + name,
                          'c_regression_acknowledged_status=' + name + ':fail']
                kind = 'acknowledged'
            else:
                lines += ['c_regression_failure=' + name + ':fail']
                kind = 'failure'
            for field in range(70):
                lines += ['c_regression_{}_detail={}|{}-moment_{}_1.gkyl [DIFF] '
                          'max_abs={} max_rel=1'.format(kind, name, name.rsplit('/', 1)[-1], field, field)]
        self.write('ci-regression-summary.txt', '\n'.join(lines))
        args = argparse.Namespace(platform='personal', context='ci/test', result='failure', pr='', output='ci-report.md')
        with contextlib.redirect_stdout(io.StringIO()):
            report.build_report(args)
        pages = json.loads(Path(args.output + '.json').read_text())
        self.assertEqual(len(pages), 1)
        self.assertLess(len(pages[0].encode('utf-8')), 40000)
        self.assertEqual(Path('ci-report-regressions.md').read_text().count('[DIFF]'), 39 * 70)
        for index in range(39):
            self.assertIn('3x2v_p1_{}'.format(index), pages[0])

    def test_compact_reasons_and_distinct_maxima(self):
        lines = ['test-elc_1.gkyl [DIFF] max_abs=10 max_rel=0.1',
                 'test-ion_1.gkyl [DIFF] max_abs=2 max_rel=1',
                 'test-dt.gkyl [DIFF] dynvec not compatible',
                 'test-moms.gkyl [DIFF] dynvec not compatible',
                 'test-field.gkyl [MISSING]']
        text = report.regression_preview('test', 'fail', lines, compact=True)
        self.assertIn('4 differing files, 1 missing file', text)
        self.assertIn('2 × <code>dynvec not compatible', text)
        self.assertIn('2 × <code>numerical difference', text)
        self.assertIn('max_abs=10 (test-elc_1.gkyl)', text)
        self.assertIn('max_rel=1 (test-ion_1.gkyl)', text)
        self.assertNotIn('| Files |', text)

    def test_comparison_log_is_summarized_only_with_checker_results(self):
        path = 'ci-command-logs/candidate-regression-check.log'
        self.write(path, 'UNACKNOWLEDGED FAILURES:\nrepeated file listing')
        self.write(path + '.exit', '1')
        self.assertIn('repeated file listing', '\n'.join(report.failure_sections()))
        self.regression_fixture()
        text = '\n'.join(report.failure_sections())
        self.assertIn(path, text)
        self.assertIn('Exit code: 1', text)
        self.assertNotIn('repeated file listing', text)

    def test_timing_omissions_are_preserved_in_artifact(self):
        issues = ['database: invalid execution time for test{}'.format(i) for i in range(200)]
        with patch.object(report, 'regression_timings', side_effect=[({}, issues), ({}, [])]):
            text = '\n'.join(report.regression_timing_sections())
        self.assertIn('200 candidate and 0 baseline', text)
        self.assertNotIn('test199', text)
        artifact = Path('ci-report-timing-details.md').read_text()
        for issue in issues:
            self.assertIn(issue, artifact)

    def test_frame_ranges_are_exact_and_preserve_padding(self):
        for frames, expected in [(['0', '20', '10', '40'], '0, 10, …, 20, 40'),
                                 (['001', '002', '003', '007'], '001, 002, …, 003, 007'),
                                 (['1', '3'], '1, 3'), (['9'], '9')]:
            self.assertEqual(report.frame_ranges(frames), expected)

    def test_grouping_preserves_failure_kinds_reasons_and_distinct_extrema(self):
        lines = ['elc_0.gkyl [DIFF] max_abs=9 max_rel=0.1',
                 'elc_10.gkyl [DIFF] max_abs=1 max_rel=0.9',
                 'elc_20.gkyl [MISSING]', 'elc_30.gkyl [DIFF] shape mismatch',
                 'elc_001.gkyl [DIFF] max_abs=2 max_rel=0.2',
                 'elc_40.gkyl [DIFF] max_abs=nan max_rel=inf',
                 'unrecognized <diagnostic> | reason']
        groups = report.regression_file_groups(lines)
        self.assertEqual(len(groups), 6)
        self.assertEqual(groups[0]['kind'], 'MISSING')
        combined = next(group for group in groups if len(group['files']) == 2)
        self.assertEqual(combined['files'], ['elc_0.gkyl', 'elc_10.gkyl'])
        self.assertEqual(combined['metrics']['max_abs'][2], 'elc_0.gkyl')
        self.assertEqual(combined['metrics']['max_rel'][2], 'elc_10.gkyl')
        preview = report.regression_preview('rt_test', 'fail', lines)
        self.assertIn('max_abs=9 (elc_0.gkyl)', preview)
        self.assertIn('max_rel=0.9 (elc_10.gkyl)', preview)
        self.assertIn('max_rel=inf (elc_40.gkyl)', preview)
        self.assertIn('max_abs=nan', preview)
        self.assertNotIn('<diagnostic>', preview)
        self.assertIn('&#124;', preview)

    def test_preview_limits_groups_and_includes_omitted_worst_file(self):
        lines = ['species_{}_frame_0.gkyl [DIFF] max_abs={} max_rel={}'.format(i, i, i / 100)
                 for i in range(20)] + ['geometry.gkyl [MISSING]']
        preview = report.regression_preview('rt_test', 'fail', lines)
        rows = [line for line in preview.splitlines() if line.startswith('| <code>')]
        self.assertEqual(len(rows), 10)
        self.assertIn('geometry.gkyl', rows[0])
        self.assertIn('11 additional files in 11 groups omitted', preview)
        self.assertIn('20 differing files, 1 missing file.', preview)
        self.assertIn('max_abs=19 (species_19_frame_0.gkyl)', preview)
        self.assertIn('max_rel=0.19 (species_19_frame_0.gkyl)', preview)

    def test_irregular_frame_groups_are_bounded_but_full_lists_remain_complete(self):
        lines = ['elc_{}.gkyl [DIFF] max_abs=1'.format(i * i) for i in range(500)]
        preview = report.regression_preview('rt_test', 'fail', lines)
        self.assertLess(len(preview), 1500)
        self.assertIn('(500 files)', preview)
        self.assertIn('…', preview)

    def test_many_test_previews_paginate_without_splitting_tables(self):
        lines = ['species_{}_0.gkyl [DIFF] {}'.format(i, 'long reason ' * 30) for i in range(10)]
        blocks = [report.regression_preview('rt_test_{}'.format(i), 'fail', lines) for i in range(40)]
        sections = report.regression_dropdowns('Unacknowledged tests', blocks)
        pages = report.report_pages('Summary', sections, 'ci/test')
        self.assertGreater(len(pages), 1)
        for block in blocks:
            self.assertEqual(sum(block in page for page in pages), 1)
        for page in pages:
            self.assertLessEqual(len(page.encode('utf-8')), report.COMMENT_LIMIT)
            self.assertEqual(page.count('<details>'), page.count('</details>'))

    def test_statistics_preserves_cached_timings_and_failed_command_status(self):
        self.write('ci-timing-summary.txt', 'baseline_install_seconds=0\n'
                   'cached_baseline_install_seconds=3661\ncandidate_install_seconds=61\n')
        self.write('baseline-install-seconds.txt', '3661\n')
        self.write('candidate-install-seconds.txt', '61\n')
        self.write('candidate-install.log', 'failure details')
        self.write('candidate-install.log.exit', '2\n')
        self.write('ci-command-logs/configure.log.exit', '0\n')
        report.main(['statistics'])
        text = Path('ci-statistics.md').read_text()
        self.assertIn('| baseline install seconds | 0 | 0m 0s |', text)
        self.assertIn('| cached baseline install seconds | 3661 | 1h 1m 1s |', text)
        self.assertIn('| candidate install | 61 | 1m 1s |', text)
        self.assertIn('[candidate-install.log](candidate-install.log) | 2 | Failed', text)
        self.assertIn('ci-command-logs/configure.log | 0 | Success', text)
        self.assertEqual(Path('baseline-install-seconds.txt').read_text(), '3661\n')
        self.assertEqual(Path('candidate-install.log.exit').read_text(), '2\n')

    def test_statistics_preserves_run_metadata(self):
        metadata = {
            'baseline-slurm-parallel-regression-job-id.txt': '2963900\n',
            'slurm-regression-job-id.txt': '2963976\n',
            'slurm-parallel-regression-job-status.txt': 'State=FAILED\nExitCode=1:0\n',
            'ci-regression-summary.txt': 'passed=3\nfailed=1\n',
            'ci-parallel-regression-summary.txt': 'passed=2\n',
            'ci-failure-summary.txt': 'failure <detail> | context\n',
            'ci-report-delivery.json': '{"comment": false, "errors": ["offline"]}',
            'ci-candidate-commit.txt': 'a' * 40 + '\n',
            'ci-run-path.txt': '/runs/11\n',
            'ci-baseline-cache.txt': 'hit\n',
        }
        for path, value in metadata.items():
            self.write(path, value)
        report.main(['statistics'])
        text = Path('ci-statistics.md').read_text()
        for path, value in metadata.items():
            escaped = report.html.escape(value.rstrip('\n')).replace('|', '&#124;').replace('\n', '<br>')
            self.assertIn('| {} | {} |'.format(path, escaped), text)
            self.assertEqual(Path(path).read_text(), value)
        self.assertIn('State=FAILED<br>ExitCode=1:0', text)

    def test_statistics_handles_early_failure_and_malformed_measurements(self):
        report.main(['statistics'])
        self.assertIn('# CI statistics', Path('ci-statistics.md').read_text())
        self.write('candidate-install-seconds.txt', '<unknown>|\nmore')
        self.write('candidate-install.log.exit', '')
        report.main(['statistics'])
        text = Path('ci-statistics.md').read_text()
        self.assertIn('&lt;unknown&gt;&#124;<br>more | —', text)
        self.assertIn('| candidate-install.log |  | Unknown |', text)

    def test_baseline_cache_status_is_on_baseline_line(self):
        sha = 'a' * 40
        self.write('ci-baseline-commit.txt', sha)
        for status, label in [('hit', 'loaded from cache'), ('saved', 'saved to cache'),
                              ('miss', 'cache miss; not saved')]:
            with self.subTest(status=status):
                self.write('ci-baseline-cache.txt', f'status={status}\nbaseline_commit={sha}\n')
                args = argparse.Namespace(platform='personal', context='ci/test', result='success', pr='', output='ci-report.md')
                with contextlib.redirect_stdout(io.StringIO()):
                    report.build_report(args)
                line = next(line for line in Path(args.output).read_text().splitlines() if '**Baseline:**' in line)
                self.assertIn(f'({label})', line)
                self.assertIn(sha[:7], line)

    def test_cached_compiler_warnings_compare_across_persistent_roots(self):
        baseline = '/ci/baseline-cache/personal/' + 'a' * 40 + '/gkeyll/'
        candidate = '/ci/runs/personal/jenkins-run-2/candidate/' + 'b' * 40 + '/gkeyll/'
        self.write('baseline-unit-build.log', baseline + 'core/example.c:10:2: warning: existing warning\n')
        self.write('baseline-unit-build.log.exit', '0\n')
        self.write('candidate-unit-build.log', candidate + 'core/example.c:20:3: warning: existing warning\n'
                   + candidate + 'core/example.c:30:4: warning: introduced warning\n')
        report.diagnostic_sections('diagnostics.json')
        summary = json.loads(Path('diagnostics.json').read_text())
        self.assertEqual(summary['new_warnings'], 1)
        self.assertEqual(summary['unclassified_warnings'], 0)

    def test_large_captured_log_inventory_is_archived_without_comment_pages(self):
        paths = ['gkylsoft/gkeyll-results/' + 'long-directory/' * 10 +
                 'regressiondb:test_{}'.format(i) for i in range(674)]
        logs = [(path, 'candidate-regression/test_{}'.format(i), '', True)
                for i, path in enumerate(paths)]
        logs.append(('candidate-unit-build.log', 'candidate-unit-build.log',
                     'example.c:1: warning: example warning\nexample.c:2: error: example error\n', False))
        args = argparse.Namespace(platform='personal', context='ci/test', result='success', pr='', output='ci-report.md')
        with patch.object(report, 'captured_logs', return_value=iter(logs)), \
                contextlib.redirect_stdout(io.StringIO()):
            report.build_report(args)
        summary = json.loads(Path('ci-diagnostic-summary.json').read_text())
        self.assertEqual(summary['captured_log_paths'], paths + ['candidate-unit-build.log'])
        self.assertEqual(summary['captured_logs'], len(logs))
        pages = json.loads(Path(args.output + '.json').read_text())
        self.assertEqual(len(pages), 1)
        self.assertNotIn('Captured logs', pages[0])
        self.assertNotIn('regressiondb:test_', pages[0])
        self.assertIn('All compiler and other warnings (1)', pages[0])
        self.assertIn('example warning', pages[0])
        self.assertIn('example warning', Path('ci-report-warnings.md').read_text())
        self.assertIn('All errors (1)', pages[0])
        self.assertIn('example error', pages[0])

    def test_large_timestep_warning_report_preserves_artifact_without_comment_flood(self):
        warnings = ['WARNING: Time-step dt = {}e-16 is below 0.001*dt_init ... num_failures = 1'.format(i)
                    for i in range(1, 17001)]
        logs = [('baseline:test', 'baseline-regression/test', warnings[0], True),
                ('candidate:test', 'candidate-regression/test', '\n'.join(warnings), True),
                ('baseline-build.log', 'baseline-build.log', 'a.c:1: warning: existing\n', True),
                ('candidate-build.log', 'candidate-build.log',
                 'a.c:2: warning: existing\na.c:3: warning: introduced\n', True)]
        args = argparse.Namespace(platform='personal', context='ci/test', result='success', pr='', output='ci-report.md')
        with patch.object(report, 'captured_logs', return_value=iter(logs)), \
                contextlib.redirect_stdout(io.StringIO()):
            report.build_report(args)
        summary = json.loads(Path('ci-diagnostic-summary.json').read_text())
        self.assertEqual(summary['warnings'], 17004)
        self.assertEqual(summary['new_warnings'], 17000)
        pages = json.loads(Path(args.output + '.json').read_text())
        self.assertEqual(len(pages), 1)
        self.assertIn('16999 distinct new warnings grouped', pages[0])
        self.assertIn('warning: introduced', pages[0])
        self.assertIn('warning: existing', pages[0])
        self.assertIn('All compiler and other warnings (3)', pages[0])
        self.assertIn('All time-step warnings (17001)', pages[0])
        self.assertIn('New time-step warnings vs main (16999)', pages[0])
        archived = log_text(Path('ci-report-warnings.md').read_text())
        for warning in warnings:
            self.assertIn(warning, archived)
        self.assertIn('a.c:1: warning: existing', archived)

    def test_warning_categories_are_separate_in_comments_and_archive(self):
        compiler = 'a.c:7: warning: compiler diagnostic'
        timestep = 'WARNING: Time-step dt = nan is below 0.001*dt_init ... num_failures = 20'
        error = 'ERROR: Time-step was below 0.001*dt_init 20 consecutive times. Aborting simulation.'
        self.write('candidate-build.log', '\n'.join([compiler, timestep, error]))
        self.write('baseline-build.log', '')
        self.write('baseline-build.log.exit', '0\n')
        sections = report.diagnostic_sections('diagnostics.json', 'ci-report-warnings.md')
        for body in ('\n\n'.join(sections), Path('ci-report-warnings.md').read_text()):
            for section in body.split('</details>'):
                if '<summary>New compiler' in section or '<summary>All compiler' in section:
                    self.assertIn(compiler, log_text(section))
                    self.assertNotIn(timestep, log_text(section))
                elif '<summary>New time-step' in section or '<summary>All time-step' in section:
                    self.assertIn(timestep, log_text(section))
                    self.assertNotIn(compiler, log_text(section))
            self.assertIn('New time-step warnings vs main (1)', body)
            self.assertIn('New compiler and other warnings vs main (1)', body)
        self.assertIn(error, log_text(sections[-1]))
        summary = json.loads(Path('diagnostics.json').read_text())
        self.assertEqual((summary['warnings'], summary['new_warnings'], summary['errors']), (2, 2, 1))

    def test_kernel_warnings_are_excluded_but_errors_and_other_warnings_remain(self):
        kernel = [
            'ker/example.c:12:3: warning: generated kernel warning',
            './core/ker/example.h:12: warning: generated kernel warning',
            '../moments/ker/example.c:12:3: warning: generated kernel warning',
            '/ci/run/vlasov/ker/example.c:12:3: warning: generated kernel warning',
            'gyrokinetic/ker/example.cu(12): warning #177-D: generated kernel warning',
            '"pkpm/ker/example.c"(12, 3): warning #167: generated kernel warning',
            '\x1b[33mcustom/nested/ker/example.cpp:12:3: warning: generated kernel warning\x1b[0m',
        ]
        retained = [
            'core/zero/example.c:12:3: warning: refers to ker/example.h',
            'worker/example.c:12:3: warning: worker diagnostic',
            'core/kernels/example.c:12:3: warning: different directory',
            'core/zero/example.c:14:3: warning: outside kernel',
            'warning: unknown source location',
        ]
        output = '\n'.join(kernel + retained + [
            'In file included from core/ker/example.h:1:',
            'core/zero/example.c:15:3: warning: included nonkernel header',
            '  #include "ker/example.h"',
            'core/ker/example.c:16:3: error: kernel compilation failed',
        ])
        for tree in ('candidate', 'baseline'):
            self.write(tree + '-build.log', output)
        self.write('baseline-build.log.exit', '0\n')
        body = '\n'.join(report.diagnostic_sections('diagnostics.json', 'warnings.md'))
        for rendered in (body, Path('warnings.md').read_text()):
            self.assertNotIn('generated kernel warning', rendered)
            for warning in retained:
                self.assertIn(warning, log_text(rendered))
            self.assertIn('included nonkernel header', rendered)
        self.assertIn('kernel compilation failed', body)
        summary = json.loads(Path('diagnostics.json').read_text())
        self.assertEqual((summary['warnings'], summary['new_warnings'], summary['errors']), (12, 0, 2))
        self.assertEqual(Path('candidate-build.log').read_text(), output)

    def test_candidate_only_kernel_warning_is_not_new(self):
        self.write('baseline-build.log', '')
        self.write('baseline-build.log.exit', '0\n')
        self.write('candidate-build.log', 'ker/example.c:1: warning: generated warning\n')
        report.diagnostic_sections('diagnostics.json')
        summary = json.loads(Path('diagnostics.json').read_text())
        self.assertEqual((summary['warnings'], summary['new_warnings']), (0, 0))

    def test_fast_math_warnings_are_archived_but_not_posted_or_counted(self):
        suppressed = 'zero/math.c:51:13: warning: use of infinity is undefined behavior [-Wnan-infinity-disabled]'
        retained = 'zero/util.c:269:12: warning: discards qualifiers [-Wincompatible-pointer-types-discards-qualifiers]'
        self.write('baseline-build.log', suppressed)
        self.write('baseline-build.log.exit', '0\n')
        self.write('candidate-build.log', suppressed + '\n' + retained)
        self.write('candidate-unmatched.log', suppressed)
        self.write('candidate-test.log', 'WARNING: runtime NaN detected\nERROR: infinity detected\n')
        body = '\n'.join(report.diagnostic_sections('diagnostics.json', 'warnings.md'))
        self.assertNotIn('[-Wnan-infinity-disabled]', body)
        self.assertIn(retained, log_text(body))
        self.assertIn('WARNING: runtime NaN detected', body)
        self.assertIn('ERROR: infinity detected', body)
        self.assertIn('All compiler and other warnings (2)', body)
        self.assertIn('New compiler and other warnings vs main (1)', body)
        self.assertIn(suppressed, log_text(Path('warnings.md').read_text()))
        summary = json.loads(Path('diagnostics.json').read_text())
        self.assertEqual((summary['warnings'], summary['new_warnings'],
                          summary['unclassified_warnings'], summary['errors']), (2, 1, 1, 1))

    def test_timestep_groups_keep_tests_thresholds_and_other_warnings_separate(self):
        def warning(path, line, dt, threshold='0.001', failures=1):
            return '{}:{}: WARNING: Time-step dt = {} is below {}*dt_init ... num_failures = {}'.format(
                path, line, dt, threshold, failures)
        entries = [warning('db:test-a', 3, '2e-12'), 'a.c:7: warning: keep this compiler message',
                   warning('db:test-b', 4, '3e-12'), warning('db:test-a', 8, '1e-12', failures=4),
                   warning('db:test-a', 9, '4e-12', threshold='0.01'), warning('db:test-a', 10, 'nan')]
        grouped = report.group_timestep_warnings(entries)
        self.assertEqual(len(grouped), 5)
        self.assertIn('db:test-a:3-8:', grouped[0])
        self.assertIn('2 distinct new warnings grouped', grouped[0])
        self.assertIn('dt range 1e-12 to 2e-12', grouped[0])
        self.assertIn('maximum num_failures = 4', grouped[0])
        self.assertEqual(grouped[1:], [entries[i] for i in (1, 2, 4, 5)])

    def test_all_compiler_warnings_are_posted_and_archived(self):
        warnings = ['a.c:{}: warning: distinct message {}'.format(i + 1, i) for i in range(2000)]
        warnings.append(warnings[0])  # All warnings must retain repeated occurrences too.
        self.write('candidate-build.log', '\n'.join(warnings) + '\nERROR: useful failure\n')
        self.write('baseline-build.log', '')
        self.write('baseline-build.log.exit', '0\n')
        args = argparse.Namespace(platform='personal', context='ci/test', result='failure', pr='', output='ci-report.md')
        with contextlib.redirect_stdout(io.StringIO()):
            report.build_report(args)
        pages = json.loads(Path(args.output + '.json').read_text())
        posted = log_text('\n'.join(pages))
        all_sections = re.findall(r'<details><summary>All compiler and other warnings.*?</details>',
                                  '\n'.join(pages), re.S)
        all_posted = log_text('\n'.join(all_sections))
        self.assertIn('All compiler and other warnings (2001)', '\n'.join(pages))
        self.assertEqual(all_posted.count(warnings[0]), 2)
        self.assertIn('New compiler and other warnings vs main (2000)', pages[0])
        self.assertIn('ERROR: useful failure', posted)
        for warning in warnings:
            self.assertIn(warning, posted)
            self.assertIn(warning, all_posted)
        for page in pages:
            self.assertLessEqual(len(page.encode('utf-8')), report.COMMENT_LIMIT)
        archived = log_text(Path('ci-report-warnings.md').read_text())
        for warning in warnings:
            self.assertIn(warning, archived)
        self.assertEqual(json.loads(Path('ci-diagnostic-summary.json').read_text())['new_warnings'], 2000)

    def test_timestep_grouping_preserves_compiler_warnings(self):
        path = 'gkylsoft/gkeyll-results/gyrokinetic/regressiondb:gyrokinetic/creg/rt_gk_cbc_2x2v_p1'
        entries = ['{}:{}: WARNING: Time-step dt = {}e-14 is below 0.001*dt_init ... num_failures = 1'.format(
            path, i + 607, i + 1) for i in range(100)]
        entries.append('a.c:1: warning: useful compiler warning')
        preview = '\n\n'.join(report.group_timestep_warnings(entries))
        self.assertIn(path + ':607-706:', preview)
        self.assertIn('100 distinct new warnings grouped', preview)
        self.assertIn('dt range 1e-14 to 1e-12', preview)
        self.assertIn(entries[-1], preview)
        self.assertNotIn('omitted', preview)

    def test_unavailable_baseline_shows_warnings_in_report_and_artifact(self):
        self.write('candidate-build.log', 'a.c:1: warning: cannot classify\n')
        body = '\n'.join(report.diagnostic_sections('diagnostics.json', 'ci-report-warnings.md'))
        self.assertIn('New compiler and other warnings vs main (unknown)', body)
        self.assertIn('Comparison unavailable', body)
        self.assertIn('All compiler and other warnings (1)', body)
        self.assertIn('warning: cannot classify', body)
        self.assertIn('warning: cannot classify', Path('ci-report-warnings.md').read_text())
        self.assertIsNone(json.loads(Path('diagnostics.json').read_text())['new_warnings'])

    def test_repeated_identical_timesteps_without_baseline_are_grouped(self):
        warning = 'WARNING: Time-step dt = 1.5263e-14 is below 0.001*dt_init ... num_failures = 1'
        self.write('candidate-build.log', '\n'.join([warning] * 100) + '\na.c:1: warning: compiler detail\n')
        body = '\n'.join(report.diagnostic_sections('diagnostics.json', 'ci-report-warnings.md'))
        self.assertIn('All compiler and other warnings (1)', body)
        self.assertIn('All time-step warnings (100)', body)
        self.assertIn('New time-step warnings vs main (unknown)', body)
        self.assertIn('100 warning occurrences grouped', body)
        self.assertEqual(body.count('Time-step dt'), 1)
        self.assertIn('a.c:1: warning: compiler detail', log_text(body))
        self.assertEqual(log_text(Path('ci-report-warnings.md').read_text()).count(warning), 100)

    def test_run_header_and_total_use_wall_time_instead_of_overlapping_steps(self):
        self.write('ci-timing-summary.txt', 'candidate_unit_build_seconds=3000\n'
                   'unit_slurm_elapsed_seconds=3500\nmissing_seconds=not-recorded\n')
        args = argparse.Namespace(platform='personal', context='ci/test', result='success', pr='', output='ci-report.md')
        with patch.dict(os.environ, CI_REPORT_START_MS='1700000000000', CI_REPORT_END_MS='1700003661000'), \
                contextlib.redirect_stdout(io.StringIO()):
            report.build_report(args)
        body = Path(args.output).read_text()
        self.assertIn('**Start:** 2023-11-14 22:13:20 UTC', body)
        self.assertIn('**End:** 2023-11-14 23:14:21 UTC', body)
        self.assertIn('**Elapsed:** 1:01:01 (3661 s)', body)
        self.assertIn('| unit slurm elapsed | 3500 |\n| **Total elapsed** | **3661** |\n\n</details>', body)
        self.assertNotIn('| missing |', body)

    def test_early_failure_and_pending_reports_have_timing_without_step_artifacts(self):
        for result in ('failure', 'cancelled', 'timed_out', 'pending'):
            with self.subTest(result=result):
                args = argparse.Namespace(platform='personal', context='ci/test', result=result, pr='', output='ci-report.md')
                with patch.dict(os.environ, CI_REPORT_START_MS='1700000000000', CI_REPORT_END_MS='1700000000000'), \
                        contextlib.redirect_stdout(io.StringIO()):
                    report.build_report(args)
                body = Path(args.output).read_text()
                self.assertIn('**Elapsed:** 0:00:00 (0 s)', body)
                self.assertIn('| **Total elapsed** | **0** |', body)
                self.assertIn('**Updated:**' if result == 'pending' else '**End:**', body)
                if result == 'pending':
                    self.assertNotIn('**End:**', body)

    def test_missing_or_invalid_start_does_not_invent_elapsed_time(self):
        for start in ('', 'invalid', '-1', '1700000001000'):
            with patch.dict(os.environ, CI_REPORT_START_MS=start, CI_REPORT_END_MS='1700000000000'):
                started, finished, elapsed = report.run_timing()
            self.assertIsNone(started)
            self.assertIsNone(elapsed)
        self.assertEqual(report.timing_section(), '')

    def test_log_lines_allow_wrapping_without_losing_text_or_indentation(self):
        lines = ['warning: ' + 'long diagnostic ' * 30, '/workspace/' + 'long-path/' * 40,
                 '   56 |     bad', '      |     ^~~', '',
                 '<script>bad & worse</script> ``` **literal**', '',
                 '</code></div></details>', '\tindented', '']
        text = '\n'.join(lines)
        rendered = report.wrapped_log(text)
        self.assertEqual(log_text(rendered), text)
        self.assertNotIn('<pre', rendered)
        self.assertNotIn('<script>', rendered)
        self.assertNotIn('</details>', rendered)
        self.assertNotIn('\n\n', rendered)
        self.assertIn('<code>   56 |     bad</code>', rendered)

    def test_wrapped_large_diagnostics_fit_comment_limits(self):
        for text in ('x' * 120000, '<&>λ漢字' * 20000, '\n' * 10000):
            with self.subTest(sample=text[:10]):
                sections = report.detail_sections('Long diagnostic', text)
                pages = report.report_pages('Summary', sections, 'ci/test')
                self.assertGreater(len(pages), 1)
                self.assertTrue(all(len(page.encode('utf-8')) <= report.COMMENT_LIMIT for page in pages))
                self.assertEqual(log_text('\n'.join(pages)), text)
                for page in pages:
                    for tag in ('details', 'div', 'code'):
                        self.assertEqual(page.count('<' + tag + '>'), page.count('</' + tag + '>'))

    def test_pipeline_provenance_is_distinct_from_candidate_and_baseline(self):
        self.write('ci-candidate-commit.txt', 'a' * 40)
        self.write('ci-baseline-commit.txt', 'b' * 40)
        self.write('ci-trusted-ci-commit.txt', 'c' * 40)
        self.write('ci-reporting-commit.txt', 'd' * 40)
        self.write('ci-trusted-checker-commit.txt', 'e' * 40)
        args = argparse.Namespace(platform='personal', context='ci/test', result='success', pr='', output='ci-report.md')
        with patch.dict(os.environ, CI_TRUSTED_CI_REF='ci-feature'), contextlib.redirect_stdout(io.StringIO()):
            report.build_report(args)
        body = Path(args.output).read_text()
        self.assertIn('**Jenkinsfile:** `ci-feature` @ [`' + 'c' * 40 + '`]', body)
        self.assertIn('**Reporting tools:** [`' + 'd' * 40 + '`]', body)
        self.assertIn('**Regression checker:** [`' + 'e' * 40 + '`]', body)
        self.assertIn('https://github.com/gkeyllorg/gkeyll/commit/' + 'c' * 40, body)

    def test_pipeline_provenance_survives_early_failure_and_never_guesses_candidate(self):
        self.write('ci-candidate-commit.txt', 'a' * 40)
        args = argparse.Namespace(platform='personal', context='ci/test', result='failure', pr='', output='ci-report.md')
        for commit in ['', 'f' * 40]:
            with patch.dict(os.environ, CI_TRUSTED_CI_COMMIT=commit), contextlib.redirect_stdout(io.StringIO()):
                report.build_report(args)
            body = Path(args.output).read_text()
            self.assertEqual('**Jenkinsfile:**' in body, bool(commit))
            self.assertNotIn('**Reporting tools:**', body)

    def test_failed_build_tail_is_collapsed_and_before_test_results(self):
        lines = ['build line {}'.format(i) for i in range(110)] + [
            'cc -c -O3 zero/array_ops.c',
            "zero/array_ops.c:56:5: error: undeclared identifier 'bad'",
            '   56 |     bad',
            '      |     ^~~',
            'zero/array_ops.c:724:19: warning: NaN disabled',
            '  724 |   if (isnan(a)) {',
            '      |       ^~~~~~~~',
            'make: *** [core] Error 2',
        ]
        self.write('candidate-unit-build.log', '\n'.join(lines) + '\n')
        self.write('ci-failure-detail.txt', 'Command: make -j8 unit\nFull log: candidate-unit-build.log (archived with the build)\n')
        self.write('candidate-unit-results.txt', 'FAIL core: ctest_example\n')
        args = argparse.Namespace(platform='personal', context='ci/test', result='failure', pr='1', output='ci-report.md')
        with contextlib.redirect_stdout(io.StringIO()):
            report.build_report(args)
        body = Path(args.output).read_text()
        tail = body.split('<summary>Failed build log — last 100 lines:')[1].split('</details>')[0]
        self.assertIn('\n'.join(lines[-100:]), log_text(tail))
        self.assertNotIn('build line 17\n', tail)
        self.assertNotIn('<details open', body)
        self.assertLess(body.index('Failed build log'), body.index('Candidate unit tests'))
        warnings = body.split('<summary>All compiler and other warnings (1)')[1].split('</details>')[0]
        errors = body.split('<summary>All errors (2)')[1].split('</details>')[0]
        self.assertNotIn('undeclared identifier', warnings)
        self.assertNotIn('NaN disabled', errors)
        self.assertIn('^~~', errors)

    def test_failed_logged_command_uses_same_tail_length(self):
        path = 'ci-command-logs/candidate-setup.log'
        lines = ['line {}'.format(i) for i in range(125)]
        self.write(path, '\n'.join(lines))
        self.write(path + '.exit', '7\n')
        body = '\n'.join(report.failure_sections())
        self.assertIn('Exit code: 7', body)
        self.assertIn('\n'.join(lines[-100:]), log_text(body))
        self.assertNotIn('line 24\n', body)

    def test_baseline_cache_rejection_reports_crashed_tests(self):
        path = 'ci-command-logs/baseline-cache-publish.log'
        self.write(path, 'baseline-cache: incomplete or failed baseline creation\n'
                   '  gyrokinetic/creg/rt_gk_cbc_2x2v_p1: crashed (-6)\n'
                   '  gyrokinetic/creg/rt_gk_rosenbluth_hinton_2x2v_p1: crashed (-6)\n')
        self.write(path + '.exit', '1\n')
        body = '\n'.join(report.failure_sections())
        self.assertIn('rt_gk_cbc_2x2v_p1: crashed (-6)', body)
        self.assertIn('rt_gk_rosenbluth_hinton_2x2v_p1: crashed (-6)', body)
        self.assertIn('Exit code: 1', body)
        self.assertNotIn('build log unavailable', body)

    def execution_database(self, rows, baseline=True, suite='gyrokinetic', guid='create', finalized=True):
        root = '_baseline/' if baseline else ''
        path = Path(root + 'gkylsoft/gkeyll-results/' + suite + '/regressiondb')
        path.parent.mkdir(parents=True, exist_ok=True)
        mode = 'cpu_parallel' if suite.startswith('parallel-') else 'cpu_serial'
        with contextlib.closing(sqlite3.connect(path)) as db, db:
            db.executescript('''
                CREATE TABLE IF NOT EXISTS RegressionMeta (guid TEXT, run_mode TEXT, ntotal INTEGER);
                CREATE TABLE IF NOT EXISTS RegressionData
                    (guid TEXT, name TEXT, test_type TEXT, status INTEGER, runtime REAL, runlog TEXT);
            ''')
            if finalized:
                db.execute('INSERT INTO RegressionMeta VALUES (?, ?, ?)', (guid, mode, len(rows)))
            db.executemany('INSERT INTO RegressionData VALUES (?, ?, ?, ?, ?, ?)',
                           [(guid, name, 'c', status, 1., output) for name, status, output in rows])
        return path

    def test_baseline_numerical_aborts_are_prominent_without_a_failed_command_log(self):
        names = ['gyrokinetic/creg/rt_gk_cbc_2x2v_p1',
                 'gyrokinetic/creg/rt_gk_rosenbluth_hinton_2x2v_p1']
        fatal = 'ERROR: Time-step was below 0.001*dt_init 20 consecutive times. Aborting simulation ....'
        output = 'WARNING: Time-step dt = 1e-14 is below 0.001*dt_init\n' * 200
        output += fatal + '\n' + 'Timing statistic: 1.0\n' * 150
        self.execution_database([(name, -6, output) for name in names])
        self.write('ci-failure-summary.txt', 'stage=Build baseline cache\nmessage=script returned exit code 1\n')
        self.write('baseline-c-regression-create.log.exit', '0\n')
        self.write('candidate-unit-results.txt', 'PASS core: ctest_array\n')
        args = argparse.Namespace(platform='personal', context='ci/test', result='failure', pr='', output='ci-report.md')
        with contextlib.redirect_stdout(io.StringIO()):
            report.build_report(args)
        body = Path(args.output).read_text()
        section = body.split('<summary>Baseline regression execution failures (2)</summary>')[1].split('</details>')[0]
        for name in names:
            self.assertIn(name + ' [cpu_serial]: Simulation aborted: time-step collapse', log_text(section))
        self.assertIn(fatal, log_text(section))
        self.assertIn('cannot be published as a baseline cache', section)
        self.assertIn('RegressionData.runlog', section)
        self.assertNotIn('Timing statistic', section)
        self.assertNotIn('WARNING:', section)
        self.assertNotIn('build log unavailable', body)
        self.assertLess(body.index('Baseline regression execution failures'), body.index('Candidate unit tests'))
        self.assertIn('baseline regressions (2 execution failures)', Path('ci-status-description.txt').read_text())

    def test_execution_failure_modes_and_fallback_output(self):
        self.execution_database([('timeout', -3, 'last progress <1%>'), ('empty', -5, '')],
                                suite='parallel-c-4/gyrokinetic')
        self.execution_database([('compile', -4, 'a.c:1: error: broken'),
                                 ('crash', -6, 'Segmentation fault'),
                                 ('difference', 0, 'different arrays'),
                                 ('passed', 1, ''), ('skipped', -1, '')], baseline=False)
        sections, counts = report.regression_execution_failure_sections()
        body = '\n'.join(sections)
        self.assertEqual(counts, {'Baseline': 2, 'Candidate': 2})
        self.assertIn('timeout [cpu_parallel]: Timed out', body)
        self.assertIn('last progress &lt;1%&gt;', body)
        self.assertIn('No per-test output was captured', body)
        self.assertIn('compile [cpu_serial]: Compilation failed', body)
        self.assertIn('Segmentation fault', body)
        self.assertNotIn('difference [', body)
        self.assertNotIn('passed [', body)
        self.assertNotIn('skipped [', body)

    def test_execution_failures_do_not_reuse_historical_results(self):
        self.execution_database([('old_crash', -6, 'ERROR: old')], guid='old')
        self.execution_database([('new_pass', -2, '')], guid='new')
        self.assertEqual(report.regression_execution_failure_sections(), ([], {}))
        self.execution_database([('in_progress', -6, 'ERROR: partial')], guid='unfinished', finalized=False)
        sections, counts = report.regression_execution_failure_sections()
        self.assertEqual(counts, {})
        self.assertIn('latest invocation is unfinished', '\n'.join(sections))
        self.assertNotIn('old_crash', '\n'.join(sections))

    def test_unreadable_execution_database_does_not_break_reporting(self):
        self.write('_baseline/gkylsoft/gkeyll-results/gyrokinetic/regressiondb', 'not sqlite')
        sections, counts = report.regression_execution_failure_sections()
        self.assertEqual(counts, {})
        self.assertIn('details unavailable', '\n'.join(sections))

    def test_failed_build_path_aliases_share_one_tail_and_metadata(self):
        path = 'candidate-unit-build.log'
        lines = ['build line {}'.format(i) for i in range(125)] + ['make: *** [unit] Error 2']
        self.write(path, '\n'.join(lines) + '\n')
        self.write(path + '.exit', '2\n')
        Path('build-log-alias').symlink_to(path)
        for recorded_path in (str(Path(path).resolve()), './' + path, 'build-log-alias'):
            with self.subTest(recorded_path=recorded_path):
                self.write('ci-failure-detail.txt', 'Command: make -j8 unit\n'
                           'Full log: {} (archived with the build)\n'.format(recorded_path))
                body = '\n'.join(report.failure_sections())
                self.assertEqual(body.count('<summary>Failed build log'), 1)
                self.assertIn('last 100 lines: ' + path + '</summary>', body)
                displayed = log_text(body)
                self.assertIn('Command: make -j8 unit\nExit code: 2', displayed)
                tail = displayed.split('Full log: ' + path + '\n\n', 1)[1]
                self.assertEqual(tail.splitlines(), lines[-100:])

    def test_distinct_failed_logs_with_same_basename_keep_separate_tails(self):
        for directory in ('build', '_baseline/build'):
            path = directory + '/compile.log'
            self.write(path, directory + ' failure\n')
            self.write(path + '.exit', '2\n')
        body = '\n'.join(report.failure_sections())
        self.assertEqual(body.count('<summary>Failed build log'), 2)
        for directory in ('build', '_baseline/build'):
            self.assertIn('last 100 lines: ' + directory + '/compile.log</summary>', body)

    def test_slurm_fallback_and_missing_logs(self):
        self.assertIn('build log unavailable', report.failure_sections()[0])
        self.write('slurm-unit-123.out', 'compiler command\nerror: failed\n')
        body = '\n'.join(report.failure_sections())
        self.assertIn('last 100 lines: slurm-unit-123.out', body)
        self.assertIn('compiler command\nerror: failed', log_text(body))

    def test_warning_mentioning_error_is_not_an_error(self):
        self.write('candidate-unit-build.log', 'a.c:4: warning: error: in diagnostic message\n')
        body = '\n'.join(report.diagnostic_sections())
        self.assertIn('All compiler and other warnings (1)', body)
        self.assertIn('All errors (0)', body)


if __name__ == '__main__':
    unittest.main()
