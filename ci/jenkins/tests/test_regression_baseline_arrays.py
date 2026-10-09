"""Regression acknowledgments must respect test, run-mode, and file scope."""
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
class RegressionAcknowledgmentTests(unittest.TestCase):
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
            f'{TEST} {MODE} {file} # repaired writer'
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
        summary.unlink(missing_ok=True)
        result = subprocess.run([
            GKEYLL, '-S', str(ROOT / 'ci/jenkins/check_regression_results.lua'),
            str(self.results), str(self.ack), str(summary), str(self.baseline_ack)],
            cwd=ROOT, capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode, 0 if passed else 1, result.stdout + result.stderr)
        self.assertTrue(summary.exists(), result.stdout + result.stderr)
        self.assertIn(f'c_regression_unacknowledged={0 if passed else 1}', summary.read_text())

    def test_named_baseline_arrays(self):
        self.check(passed=True)

    def test_all_forms_names_and_modes(self):
        modes = ('cpu_serial', 'cpu_parallel', 'gpu_serial', 'gpu_parallel')
        for name in (TEST, 'gyrokinetic/rt_fixture', 'rt_fixture'):
            for mode in modes:
                for scope in ('test', 'mode', 'file', 'mode_file'):
                    with self.subTest(name=name, mode=mode, scope=scope):
                        prefix = f'{name} {mode}' if scope in ('mode', 'mode_file') else name
                        entries = [f'{prefix} {file}' for file in FILES] if scope in ('file', 'mode_file') else [prefix]
                        self.ack.write_text('\n'.join(entries))
                        for run_mode in modes:
                            self.db.execute('UPDATE RegressionMeta SET run_mode=?', (run_mode,))
                            self.check(passed=scope in ('test', 'file') or run_mode == mode)

    def test_whitespace_and_comments(self):
        self.ack.write_text('# ignored\n\n' + '\n'.join(
            f'  {TEST}\t {MODE}   {file}  # expected change' for file in FILES))
        self.check(passed=True)

    def test_inherited_acknowledgment_is_inert(self):
        for entry in (TEST, f'{TEST} {MODE}', f'{TEST} {FILES[0]}', f'{TEST} {MODE} {FILES[0]}'):
            with self.subTest(entry=entry):
                self.ack.write_text(entry + ' # old reason')
                self.baseline_ack.write_text(self.ack.read_text())
                self.check([f'{FILES[0]}  [DIFF]'])
                self.ack.write_text(entry + ' # new reason')
                self.check([f'{FILES[0]}  [DIFF]'], passed=True)

    def test_wrong_mode(self):
        self.db.execute("UPDATE RegressionMeta SET run_mode='cpu_serial'")
        self.check()

    def test_wrong_test(self):
        for entry in (TEST, f'{TEST} {MODE}', f'{TEST} {FILES[0]}', f'{TEST} {MODE} {FILES[0]}'):
            with self.subTest(entry=entry):
                self.ack.write_text(entry.replace(TEST, TEST + '_other'))
                self.check([f'{FILES[0]}  [DIFF]'])

    def test_missing_acknowledgment(self):
        self.ack.write_text(self.ack.read_text().splitlines()[0])
        self.check()

    def test_other_failure_in_same_test(self):
        for line in ('other.gkyl  [DIFF]  ' + DETAIL, 'other.gkyl  [MISSING]',
                     FILES[0] + '.extra  [DIFF]', 'unrecognized failure',
                     FILES[0] + '  [UNKNOWN]', FILES[0] + '  [DIFF]malformed'):
            with self.subTest(line=line):
                self.check([f'{FILES[1]}  [DIFF]  {DETAIL}', line])

    def test_named_file_failure_types(self):
        for line in (FILES[0] + '  [DIFF]',
                     FILES[0] + '  [DIFF]  ' + DETAIL,
                     FILES[0] + '  [DIFF]  candidate array read failed',
                     FILES[0] + '  [DIFF]  max_abs=1 max_rel=1',
                     FILES[0] + '  [MISSING]',
                     FILES[0] + '  [DIFF]  crash: error'):
            with self.subTest(line=line):
                self.check([f'{FILES[1]}  [DIFF]  {DETAIL}', line], passed=True)

    def test_mixed_file_scopes(self):
        self.ack.write_text(f'{TEST} {FILES[0]}\n{TEST} {MODE} {FILES[1]}')
        self.check(passed=True)
        self.db.execute("UPDATE RegressionMeta SET run_mode='gpu_parallel'")
        self.check()

    def test_extra_fields_do_not_acknowledge(self):
        self.ack.write_text('\n'.join(
            f'{TEST} {MODE} {file} unexpected' for file in FILES))
        self.check()

    def test_empty_comparison_log(self):
        self.check([])
        self.check(['', ''])

    def test_execution_failures(self):
        for mode in ('cpu_parallel', 'gpu_parallel'):
            self.db.execute('UPDATE RegressionMeta SET run_mode=?', (mode,))
            for entry, passed in ((TEST, True), (f'{TEST} {MODE}', mode == MODE),
                                  (f'{TEST} {FILES[0]}', False),
                                  (f'{TEST} {MODE} {FILES[0]}', False)):
                self.ack.write_text(entry)
                for status in (-3, -4, -5, -6):
                    with self.subTest(mode=mode, entry=entry, status=status):
                        self.check(status=status, passed=passed)


if __name__ == '__main__':
    unittest.main()
