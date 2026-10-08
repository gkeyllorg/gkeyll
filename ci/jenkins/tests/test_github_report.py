"""Offline checks for report timing and log presentation."""
import argparse
import contextlib
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from ci.jenkins import github_report as report


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

    def test_long_log_lines_wrap_without_losing_text_or_short_line_indentation(self):
        lines = ['warning: ' + 'long diagnostic ' * 30, '/workspace/' + 'long-path/' * 40,
                 '   56 |     bad', '      |     ^~~', '']
        for line in lines:
            rendered = report.fenced(line).split('\n')[1:-1]
            self.assertTrue(all(len(part) <= 100 for part in rendered))
            self.assertEqual(rendered[0] + ''.join(part[2:] for part in rendered[1:]), line)

    def test_wrapped_large_diagnostics_fit_comment_limits(self):
        sections = report.detail_sections('Long diagnostic', 'x' * 120000)
        pages = report.report_pages('Summary', sections, 'ci/test')
        self.assertGreater(len(pages), 1)
        self.assertTrue(all(len(page.encode('utf-8')) <= report.COMMENT_LIMIT for page in pages))
        self.assertTrue(all(page.count('```') % 2 == 0 for page in pages))

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
        self.assertIn('\n'.join(lines[-100:]), tail)
        self.assertNotIn('build line 17\n', tail)
        self.assertNotIn('<details open', body)
        self.assertLess(body.index('Failed build log'), body.index('Candidate unit tests'))
        warnings = body.split('<summary>All warnings (1)')[1].split('</details>')[0]
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
        self.assertIn('\n'.join(lines[-100:]), body)
        self.assertNotIn('line 24\n', body)

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
                self.assertIn('Command: make -j8 unit\nExit code: 2', body)
                tail = body.split('Full log: ' + path + '\n\n', 1)[1].split('\n```', 1)[0]
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
        self.assertIn('compiler command\nerror: failed', body)

    def test_warning_mentioning_error_is_not_an_error(self):
        self.write('candidate-unit-build.log', 'a.c:4: warning: error: in diagnostic message\n')
        body = '\n'.join(report.diagnostic_sections())
        self.assertIn('All warnings (1)', body)
        self.assertIn('All errors (0)', body)


if __name__ == '__main__':
    unittest.main()
