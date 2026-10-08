"""Offline checks for failure log presentation."""
import argparse
import contextlib
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import github_report as report


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
