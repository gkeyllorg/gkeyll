"""Exercise installed regression tools without running any simulations.

Set GKEYLL to a built executable. Its shared libraries must be discoverable
(set LD_LIBRARY_PATH when using a build-tree executable).
"""
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
GKEYLL = os.environ.get('GKEYLL') or shutil.which('gkeyll')
LAYERS = ('moments', 'vlasov', 'gyrokinetic', 'pkpm')


@unittest.skipUnless(GKEYLL, 'Set GKEYLL to run the installed-tool fixtures')
class RegressionConfigTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix='gkeyll-regression-config-')
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.source = self.root / 'source'
        for layer in LAYERS:
            creg = self.source / layer / 'creg'
            creg.mkdir(parents=True)
            (creg / 'rt_config_fixture.c').touch()
        # A home config cannot be written here. Installed tools must keep
        # candidate and baseline configuration under their own install prefixes.
        self.test_home = self.root / 'home'
        self.test_home.mkdir()
        (self.test_home / 'runregression.config.lua').mkdir()
        self.env = dict(os.environ, HOME=str(self.test_home))

    def install(self, name, config=None):
        prefix = self.root / name / 'gkylsoft'
        bindir = prefix / 'gkeyll/bin'
        bindir.mkdir(parents=True)
        shutil.copy2(GKEYLL, bindir / 'gkeyll')
        shutil.copytree(ROOT / 'gkeyll/lua', bindir / 'lua')
        if config is not None:
            share = prefix / 'gkeyll/share'
            share.mkdir()
            (share / 'config.mak').write_text(config.format(prefix=prefix))
        return prefix, bindir / 'gkeyll'

    def run_tool(self, executable, *args, success=True):
        result = subprocess.run(
            [str(executable), '-S', *map(str, args)], cwd=self.root,
            env=self.env, capture_output=True, text=True, timeout=20)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertNotIn('*** ERROR:', result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def configure(self, executable, *args):
        return self.run_tool(executable, 'runregression', 'configure',
                             '--source-dir', self.source, *args)

    def check_loaded(self, executable):
        listed = self.run_tool(executable, 'runregression', 'list', '--c-only')
        for layer in LAYERS:
            self.assertIn(layer + '/creg/rt_config_fixture', listed.stdout)
        self.run_tool(executable, 'queryrdb', '--layer', 'moments', 'summary')

    def test_configure_and_reload_installed_prefix(self):
        configs = (
            'PREFIX={prefix}\n',
            'PREFIX := {prefix}\n',
            'GKEYLL_CI_DEPENDENCY_PREFIX=/unavailable/dependencies\n'
            'override PREFIX := {prefix}\noverride INSTALL_PREFIX := {prefix}\n',
            'PREFIX=/unavailable/dependencies\nINSTALL_PREFIX={prefix}\n',
            None,
        )
        for index, config in enumerate(configs):
            with self.subTest(config=config):
                prefix, executable = self.install(str(index), config)
                self.configure(executable)
                self.assertTrue((prefix / 'gkeyll-results/runregression.config.lua').is_file())
                for layer in LAYERS:
                    self.assertTrue((prefix / 'gkeyll-results' / layer / 'regressiondb').is_file())
                self.check_loaded(executable)

    def test_explicit_results_prefix_and_quoted_mpi_args(self):
        prefix, executable = self.install('custom', 'PREFIX={prefix}\n')
        results_prefix = self.root / 'custom results'
        results_prefix.mkdir()
        self.configure(executable, '--prefix', results_prefix,
                       '--mpi-args=--label="two words"')
        self.assertTrue((prefix / 'gkeyll-results/runregression.config.lua').is_file())
        self.assertTrue((results_prefix / 'gkeyll-results/moments/regressiondb').is_file())
        self.check_loaded(executable)

    def test_candidate_baseline_and_parallel_isolation(self):
        candidate, candidate_exec = self.install('candidate')
        baseline, baseline_exec = self.install('baseline')
        self.configure(candidate_exec, '--prefix', candidate)
        self.configure(baseline_exec, '--prefix', baseline,
                       '--results-subdir', 'parallel-c-4')
        baseline_db = baseline / 'gkeyll-results/parallel-c-4/moments/regressiondb'
        candidate_db = candidate / 'gkeyll-results/moments/regressiondb'
        self.assertTrue(baseline_db.is_file())
        self.assertTrue(candidate_db.is_file())
        # Removing one DB must not affect queries against the other install.
        baseline_db.unlink()
        self.run_tool(candidate_exec, 'queryrdb', '--layer', 'moments', 'summary')
        result = self.run_tool(baseline_exec, 'queryrdb', '--layer', 'moments',
                               'summary', success=False)
        self.assertIn(str(baseline_db), result.stdout)

    def test_configuration_write_failure(self):
        prefix, executable = self.install('blocked', 'PREFIX={prefix}\n')
        config_file = prefix / 'gkeyll-results/runregression.config.lua'
        config_file.mkdir(parents=True)
        result = self.run_tool(executable, 'runregression', 'configure',
                               '--source-dir', self.source, success=False)
        self.assertIn(str(config_file), result.stderr)
        self.assertNotIn("attempt to index local 'fn'", result.stderr)

    def test_parallel_execution_uses_shared_worker_pool(self):
        prefix, executable = self.install('parallel')
        launcher = self.root / 'mpiexec'
        launcher.write_text('#!/bin/sh\nset -eu\n'
                            'test "$1" = -n && test "$2" = 4\nshift 2\nexec "$@"\n')
        launcher.chmod(0o755)
        self.configure(executable, '--mpiexec', launcher)
        names = [f'rt_mpi_fixture_{idx}' for idx in range(4)]
        creg = self.source / 'moments/creg'
        (creg / 'c_test_manifest.lua').write_text(
            'return {parallel={' + ','.join(
                '{name="' + name + '", cuts={2,2}}' for name in names) + '}}\n')
        run_root = prefix / 'gkeyll-results/moments/creg-runs/cpu_parallel'
        for idx, name in enumerate(names):
            (creg / (name + '.c')).touch()
            run_dir = run_root / name
            run_dir.mkdir(parents=True)
            binary = run_dir / name
            binary.write_text('#!/bin/sh\nset -eu\n'
                              'test "$*" = "-M -c 2 -d 2"\n'
                              f'sleep {1.5 if idx == 0 else 0.1}\n'
                              'printf fixture > result.gkyl\n')
            binary.chmod(0o755)
        result = self.run_tool(executable, 'runregression', 'run', 'moments',
                               '--parallel', '--execute-only', '--jobs', '8', 'create')
        intervals = []
        for name in names:
            raw = (run_root / name / '_parallel_out.txt').read_text()
            markers = dict(line.split(':', 1) for line in raw.splitlines()
                           if line.startswith(('__START__:', '__END__:')))
            intervals.append((float(markers['__START__']), float(markers['__END__'])))
        self.assertLess(intervals[1][0], intervals[0][1], result.stdout)
        self.assertLess(intervals[3][0], intervals[0][1], result.stdout)
        self.assertGreaterEqual(intervals[2][0], intervals[1][1])
        self.assertGreaterEqual(intervals[3][0], intervals[2][1])
        with sqlite3.connect(prefix / 'gkeyll-results/moments/regressiondb') as db:
            rows = db.execute('SELECT name, test_type, status FROM RegressionData').fetchall()
        self.assertCountEqual(rows, [('moments/creg/' + name, 'c', -2) for name in names])

    def test_legacy_home_configuration(self):
        prefix, executable = self.install('legacy')
        self.configure(executable)
        legacy_config = self.test_home / 'runregression.config.lua'
        legacy_config.rmdir()
        (prefix / 'gkeyll-results/runregression.config.lua').rename(legacy_config)
        self.check_loaded(executable)

    def test_lua_errors_reach_the_shell(self):
        _, executable = self.install('errors')
        script = self.root / 'failure.lua'
        script.write_text('error("fixture input failure")\n')
        result = self.run_tool(executable, script, success=False)
        self.assertIn('fixture input failure', result.stderr)
        script.write_text('print("input should not run")\n')
        result = self.run_tool(executable, '-e', 'error("fixture chunk failure")',
                               script, success=False)
        self.assertIn('fixture chunk failure', result.stderr)
        self.assertNotIn('input should not run', result.stdout)
        script.write_text('this is not valid Lua\n')
        self.run_tool(executable, script, success=False)
        self.run_tool(executable, 'nonexistent-input.lua', success=False)
        self.run_tool(executable, '-e', 'assert(1 + 1 == 2)')


if __name__ == '__main__':
    unittest.main()
