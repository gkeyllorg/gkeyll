"""Narrow acknowledgments for corrupt historical arrays must fail closed."""
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
GKEYLL = os.environ.get('GKEYLL') or shutil.which('gkeyll')
TEST = 'gyrokinetic/creg/rt_fixture'
MODE = 'cpu_parallel'
FILES = ['rt_fixture-elc.gkyl', 'rt_fixture-ion.gkyl']
DETAIL = 'baseline array read failed (candidate readable)'


@unittest.skipUnless(GKEYLL, 'Set GKEYLL to run the checker fixtures')
class BaselineArrayAcknowledgmentTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix='gkeyll-baseline-array-')
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.results = self.root / 'results'
        (self.results / 'gyrokinetic').mkdir(parents=True)
        self.db = sqlite3.connect(self.results / 'gyrokinetic/regressiondb')
        self.addCleanup(self.db.close)
        self.db.executescript('''
            CREATE TABLE RegressionMeta (guid TEXT, run_mode TEXT);
            INSERT INTO RegressionMeta VALUES ('fixture', 'cpu_parallel');
            CREATE TABLE RegressionData (guid TEXT, name TEXT, status INTEGER, runlog TEXT);
        ''')
        self.ack = self.root / 'ack.txt'
        self.baseline_ack = self.root / 'baseline-ack.txt'
        self.baseline_ack.write_text('')
        self.ack.write_text('\n'.join(
            f'{TEST} {MODE} {file} baseline-array-read-failed # repaired writer'
            for file in FILES))

    def check(self, lines=None, status=0, passed=False):
        if lines is None:
            lines = [f'{file}  [DIFF]  {DETAIL}' for file in FILES]
        log = '--- Comparison failures ---\n' + '\n'.join(lines)
        self.db.execute('DELETE FROM RegressionData')
        self.db.execute('INSERT INTO RegressionData VALUES (?, ?, ?, ?)',
                        ('fixture', TEST, status, log))
        self.db.commit()
        summary = self.root / 'summary.txt'
        result = subprocess.run([
            GKEYLL, '-S', str(ROOT / 'ci/jenkins/check_regression_results.lua'),
            str(self.results), str(self.ack), str(summary), str(self.baseline_ack)],
            cwd=ROOT, capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode, 0 if passed else 1, result.stdout + result.stderr)
        self.assertTrue(summary.exists(), result.stdout + result.stderr)
        self.assertIn(f'c_regression_unacknowledged={0 if passed else 1}', summary.read_text())

    def test_only_named_baseline_arrays(self):
        self.check(passed=True)

    def test_inherited_acknowledgment_is_inert(self):
        self.baseline_ack.write_text(self.ack.read_text())
        self.check()

    def test_wrong_mode(self):
        self.db.execute("UPDATE RegressionMeta SET run_mode='cpu_serial'")
        self.check()

    def test_wrong_test(self):
        self.ack.write_text(self.ack.read_text().replace(TEST, TEST + '_other'))
        self.check()

    def test_missing_acknowledgment(self):
        self.ack.write_text(self.ack.read_text().splitlines()[0])
        self.check()

    def test_other_failure_in_same_test(self):
        for line in ('other.gkyl  [DIFF]  ' + DETAIL,
                     FILES[0] + '  [DIFF]  candidate array read failed',
                     FILES[0] + '  [DIFF]  max_abs=1 max_rel=1',
                     FILES[0] + '  [MISSING]',
                     FILES[0] + '  [DIFF]  crash: error'):
            with self.subTest(line=line):
                self.check([f'{FILES[1]}  [DIFF]  {DETAIL}', line])

    def test_empty_comparison_log(self):
        self.check([])
        self.check(['', ''])

    def test_execution_failures(self):
        for status in (-3, -4, -5, -6):
            with self.subTest(status=status):
                self.check(status=status)


if __name__ == '__main__':
    unittest.main()
