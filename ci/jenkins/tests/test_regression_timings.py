"""Offline regression timing reports using runregression's SQLite schema."""
import argparse
import contextlib
import io
import os
from pathlib import Path
import sqlite3
import tempfile
import unittest

from ci.jenkins import github_report as report


class RegressionTimingTests(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory()
        self.previous = os.getcwd()
        os.chdir(self.workspace.name)

    def tearDown(self):
        os.chdir(self.previous)
        self.workspace.cleanup()

    def database(self, rows, baseline=False, suite='moments', mode='cpu_serial', guid='run', finalized=True):
        root = '_baseline/' if baseline else ''
        path = Path(root + 'gkylsoft/gkeyll-results/' + suite + '/regressiondb')
        path.parent.mkdir(parents=True, exist_ok=True)
        with contextlib.closing(sqlite3.connect(path)) as db, db:
            db.executescript('''
                CREATE TABLE IF NOT EXISTS RegressionMeta (guid TEXT, run_mode TEXT, ntotal INTEGER);
                CREATE TABLE IF NOT EXISTS RegressionData
                    (guid TEXT, name TEXT, test_type TEXT, status INTEGER, runtime REAL, runlog TEXT);
            ''')
            if finalized:
                db.execute('INSERT INTO RegressionMeta VALUES (?, ?, ?)', (guid, mode, len(rows)))
            db.executemany('INSERT INTO RegressionData VALUES (?, ?, ?, ?, ?, ?)',
                           [(guid, name, kind, status, seconds, '') for name, kind, status, seconds in rows])
        return path

    def sections(self):
        return report.regression_timing_sections()

    def test_top_twenty_and_build_integration(self):
        self.database([('test{:02}'.format(i), 'c', 1, i) for i in range(1, 26)])
        # Execution failures and skipped tests must not crowd out slow tests.
        self.database([('timeout', 'c', -3, 1000), ('skip', 'c', -1, 1000),
                       ('different', 'c', 0, 26)], suite='vlasov')
        args = argparse.Namespace(platform='personal', context='ci/test', result='success', pr='', output='report.md')
        with contextlib.redirect_stdout(io.StringIO()):
            report.build_report(args)
        body = Path('report.md').read_text()
        table = body.split('<summary>20 slowest candidate regression tests</summary>')[1].split('</details>')[0]
        self.assertEqual(table.count('\n| `'), 20)
        self.assertIn('`different`', table)
        self.assertIn('| 26.000 | numerical difference |', table)
        self.assertIn('`test07`', table)
        self.assertNotIn('`test06`', table)
        self.assertNotIn('`timeout`', table)
        self.assertNotIn('`skip`', table)
        self.assertLess(table.index('`test25`'), table.index('`test24`'))
        self.assertIn('not statistically significant', body)
        self.assertEqual(Path('ci-status-description.txt').read_text().strip(),
                         'Passed: unit tests and C regressions.')

    def test_thresholds_in_both_directions_and_all_changes_reported(self):
        pairs = [('slower', 5, 10), ('faster', 10, 5), ('fast_base', 4.999, 20),
                 ('fast_candidate', 20, 4.999), ('below_ratio', 5, 9.999), ('equal', 8, 8)]
        pairs += [('change{:02}'.format(i), 5, 20 + i) for i in range(25)]
        self.database([(name, 'c', 1, cand) for name, base, cand in pairs])
        self.database([(name, 'c', -2, base) for name, base, cand in pairs], baseline=True)
        body = self.sections()[1]
        self.assertIn('26 slower, 1 faster', body)
        self.assertIn('| 5.000 | 10.000 | +5.000 | 2.000x | slower |', body)
        self.assertIn('| 10.000 | 5.000 | -5.000 | 0.500x | faster |', body)
        for name in ('fast_base', 'fast_candidate', 'below_ratio', 'equal'):
            self.assertNotIn('`' + name + '`', body)
        self.assertEqual(body.count('\n| `'), 27)
        self.assertLess(body.index('`change24`'), body.index('`change00`'))

    def test_latest_invocation_identity_and_unsuccessful_results(self):
        self.database([('old', 'c', 1, 100)], guid='old')
        self.database([('current', 'c', 1, 10), ('diff', 'c', 0, 100),
                       ('crash', 'c', -6, 100), ('kind', 'lua', 1, 50)], guid='new')
        self.database([('old', 'c', -2, 10), ('current', 'c', -2, 5),
                       ('diff', 'c', -2, 5), ('crash', 'c', -2, 5),
                       ('kind', 'c', -2, 5)], baseline=True)
        self.database([('current', 'c', 1, 100)], suite='parallel-c-4/moments', mode='cpu_parallel')
        self.database([('current', 'c', -2, 5)], baseline=True,
                      suite='parallel-c-4/moments', mode='gpu_parallel')
        body = self.sections()[1]
        self.assertEqual(body.count('\n| `'), 1)
        self.assertIn('`current` | `moments` | `c` | `cpu_serial`', body)
        self.assertNotIn('`old`', '\n'.join(self.sections()))

    def test_invalid_duplicate_and_missing_timings(self):
        self.assertEqual(self.sections(), [])
        self.database([('zero', 'c', 1, 0), ('negative', 'c', 1, -3),
                       ('missing', 'c', 1, None), ('text', 'c', 1, 'bad'),
                       ('infinite', 'c', 1, float('inf')), ('nan', 'c', 1, float('nan')),
                       ('duplicate', 'c', 1, 10), ('duplicate', 'c', 1, 20),
                       ('good', 'c', 1, 10)])
        values, issues = report.regression_timings('gkylsoft/gkeyll-results')
        self.assertEqual(len(values), 1)
        self.assertEqual(len(issues), 7)
        self.assertIn('comparison unavailable', self.sections()[1])

    def test_unfinished_invocation_does_not_reuse_old_timings(self):
        self.database([('old', 'c', 1, 10)], guid='old')
        self.database([('new', 'c', 1, 50)], guid='new', finalized=False)
        values, issues = report.regression_timings('gkylsoft/gkeyll-results')
        self.assertEqual(values, {})
        self.assertIn('unfinished', issues[0])

    def test_incomplete_and_legacy_databases_do_not_break_report(self):
        path = self.database([('test', 'c', 1, 10)])
        with contextlib.closing(sqlite3.connect(path)) as db, db:
            db.execute('UPDATE RegressionMeta SET ntotal=2')
        self.assertIn('incomplete results', '\n'.join(self.sections()))
        path.write_bytes(b'not a database')
        self.assertIn('timing data unavailable', '\n'.join(self.sections()))
        path.unlink()
        with contextlib.closing(sqlite3.connect(path)) as db, db:
            db.execute('CREATE TABLE RegressionMeta (guid TEXT)')
        self.assertIn('timing data unavailable', '\n'.join(self.sections()))

    def test_cache_caveat_and_separate_parallel_comparison(self):
        for suite, mode in [('moments', 'cpu_serial'), ('parallel-c-4/moments', 'cpu_parallel')]:
            self.database([('test', 'c', 1, 20)], suite=suite, mode=mode)
            self.database([('test', 'c', -2, 5)], baseline=True, suite=suite, mode=mode)
        Path('ci-baseline-cache.txt').write_text('status=hit\n')
        body = self.sections()[1]
        self.assertIn('2 slower, 0 faster', body)
        self.assertIn('measured in an earlier run', body)
        self.assertIn('`cpu_serial`', body)
        self.assertIn('`cpu_parallel`', body)


if __name__ == '__main__':
    unittest.main()
