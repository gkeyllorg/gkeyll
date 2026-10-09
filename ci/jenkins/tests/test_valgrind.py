"""CPU configuration, memory-check failure propagation, and complete reports."""
import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from ci.jenkins import github_report, valgrind
from ci.jenkins.tests.test_github_report import log_text


ROOT = Path(__file__).resolve().parents[3]
RUNNER = ROOT / 'core/minus/valcheck.sh'


class ValgrindTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.previous = Path.cwd()
        os.chdir(self.temp.name)
        Path('config.mak').write_text('CC=cc\nARCH_FLAGS=-march=native\nBUILD_APP=core\n')

    def tearDown(self):
        os.chdir(self.previous)
        self.temp.cleanup()

    def fixture(self, names):
        Path('core/unit').mkdir(parents=True)
        Path('build/core/unit').mkdir(parents=True)
        for name in names:
            Path(f'core/unit/ctest_{name}.c').touch()
            test = Path(f'build/core/unit/ctest_{name}')
            test.write_text('#!/bin/sh\nexit 0\n')
            test.chmod(0o755)
        Path('Makefile').write_text('.PHONY: core-unit\ncore-unit:\n\t@true\n')

    def fake_valgrind(self):
        Path('bin').mkdir()
        executable = Path('bin/valgrind')
        executable.write_text('''#!/bin/sh
case "$*" in
  *ctest_crash*) echo 'illegal instruction' >&2; exit 132 ;;
  *ctest_error*) echo 'Invalid write of size 4' >&2; echo 'ERROR SUMMARY: 1 errors' >&2; exit 99 ;;
  *ctest_assertion*) echo 'Assertion failed'; echo 'ERROR SUMMARY: 0 errors' >&2; exit 1 ;;
  *) echo 'ERROR SUMMARY: 0 errors' >&2; exit 0 ;;
esac
''')
        executable.chmod(0o755)
        return str(Path('bin').resolve()) + os.pathsep + os.environ['PATH']

    def test_flag_validation_and_disabled_noop(self):
        for value in ('', 'true', '2', ' 1'):
            with patch.dict(os.environ, GKEYLL_USE_VALGRIND=value):
                with self.assertRaises(ValueError):
                    valgrind.enabled()
        for value in ('0', '1'):
            with patch.dict(os.environ, GKEYLL_USE_VALGRIND=value):
                self.assertEqual(valgrind.enabled(), value == '1')
        with patch.dict(os.environ, GKEYLL_USE_VALGRIND='0'):
            result = subprocess.run(['python3', str(ROOT / 'ci/jenkins/valgrind.py'), 'run'])
        self.assertEqual(result.returncode, 0)
        self.assertFalse(Path('ci-valgrind').exists())

    def test_config_preserves_defaults_and_disables_instruction_sets(self):
        # Configure an older source whose Makefile lacks the new flag.
        Path('Makefile').write_text('CFLAGS ?= -O3 -g -fPIC\ninclude config.mak\n'
                                   'probe:\n\t@echo $(CFLAGS) $(ARCH_FLAGS)\n')
        valgrind.configure()
        flags = subprocess.check_output(['make', '-s', 'probe'], text=True)
        self.assertIn('-O3 -g -fPIC', flags)
        self.assertNotIn('-march=native', flags)
        self.assertIn(valgrind.FLAGS, flags)
        if shutil.which('cc'):
            macros = subprocess.check_output(['cc', *flags.split(), '-dM', '-E', '-'], input='', text=True)
            self.assertNotIn('__AVX', macros)

    def test_reject_gpu_configuration(self):
        for config in ('CC=nvcc', 'CC=ccache /opt/cuda/bin/nvcc', 'CC=cc\nUSE_CUDSS=1',
                       'CC=cc\nUSE_NCCL=1'):
            Path('config.mak').write_text(config + '\n')
            with self.assertRaisesRegex(ValueError, 'CPU'):
                valgrind.configure()

    def test_configured_layers_and_duplicate_names(self):
        self.fixture(['shared'])
        with Path('config.mak').open('a') as stream:
            stream.write('BUILD_APP=moments\n')
        for app in ('moments', 'vlasov'):
            Path(app, 'unit').mkdir(parents=True)
            Path(app, 'unit/ctest_shared.c').touch()
            Path(app, 'unit/mctest_mpi.c').touch()
            Path(app, 'unit/lctest_lua.c').touch()
        Path('Makefile').write_text('moments-unit:\n\t@true\n')
        with patch.dict(os.environ, PATH=self.fake_valgrind()):
            self.assertEqual(valgrind.run(Path('ci-valgrind/candidate')), 0)
        summary = json.loads(Path('ci-valgrind/candidate/summary.json').read_text())
        self.assertEqual(summary['layers'], ['core', 'moments'])
        self.assertEqual([t['name'] for t in summary['tests']],
                         ['core/ctest_shared', 'moments/ctest_shared'])
        for test in summary['tests']:
            self.assertTrue((Path('ci-valgrind/candidate') / test['log']).is_file())
        # Completed passing logs are omitted, but an interrupted nested log is kept.
        Path('ci-valgrind/candidate/moments/ctest_interrupted.log').write_text('partial trace')
        sections = ''.join(github_report.valgrind_sections())
        self.assertIn('partial trace', sections)
        self.assertNotIn('ERROR SUMMARY', sections)

    def test_make_profiles_and_default_build(self):
        # Read the actual top-level Makefile in isolation, without local config.
        probe = Path('probe.mak')
        probe.write_text(f'include {ROOT}/Makefile\n'
                         'probe:\n\t@echo $(CFLAGS)\n\t@echo $(SQL_CFLAGS)\n'
                         '\t@echo $(ARCH_FLAGS)\n\t@echo $(VALGRIND_APPS)\n')
        def flags(*args):
            return subprocess.check_output(['make', '-s', '-f', str(probe), 'probe',
                                            'CC=cc', 'GIT_TIP=fixture', *args], text=True)
        normal = flags('GKEYLL_USE_VALGRIND=0')
        self.assertIn('-march=native', normal)
        checked = flags('GKEYLL_USE_VALGRIND=1', 'BUILD_APP=vlasov',
                        'USE_MPI=1', 'ARCH_FLAGS=-march=native -mavx512f')
        self.assertNotIn('-march=native', checked)
        self.assertNotIn('-mavx512f', checked)
        self.assertIn('-DGKYL_HAVE_MPI', checked)
        self.assertIn('-DGKYL_HAVE_VLASOV', checked)
        self.assertIn('-g -gdwarf-4', checked)
        self.assertEqual(checked.splitlines()[-1], 'core moments vlasov')
        for row in checked.splitlines()[:3]:
            self.assertIn(valgrind.FLAGS, row)

    def test_recursive_make_cannot_restore_native_flags(self):
        output = subprocess.check_output([
            'make', '-n', 'GKEYLL_USE_VALGRIND=1', 'CC=cc',
            'ARCH_FLAGS=-march=native', 'CFLAGS=-O2 -mavx512f',
            'BUILD_DIR=valgrind-flags-fixture', 'core-unit'], cwd=ROOT, text=True)
        compiles = [line for line in output.splitlines() if line.startswith('cc ') and ' -shared ' not in line]
        self.assertTrue(compiles)
        for line in compiles:
            self.assertNotIn('-march=native', line)
            self.assertNotIn('-mavx512f', line)
            self.assertIn(valgrind.FLAGS, line)

    def test_arm_profile_uses_compiler_target(self):
        compiler = Path('arm-cc').resolve()
        compiler.write_text('#!/bin/sh\necho aarch64-linux-gnu\n')
        compiler.chmod(0o755)
        Path('config.mak').write_text(f'CC={compiler}\n')
        valgrind.configure()
        config = Path('config.mak').read_text()
        self.assertIn('-march=armv8-a', config)
        self.assertNotIn('-march=x86-64', config)

    def test_script_valgrind_options_and_precedence(self):
        for script, output in ((ROOT / 'configure', 'config.mak'),
                               (ROOT / 'install-deps/mkdeps.sh', 'build-opts.sh')):
            cases = (
                ({}, [], '0'),
                ({}, ['--use-valgrind=yes'], '1'),
                ({}, ['--use-valgrind=no'], '0'),
                ({'GKEYLL_USE_VALGRIND': '1'}, [], '1'),
                ({'GKEYLL_USE_VALGRIND': '0'}, ['--use-valgrind=yes'], '1'),
                ({'GKEYLL_USE_VALGRIND': '1'}, ['--use-valgrind=no'], '0'),
                ({}, ['--use-valgrind=yes', '--use-valgrind=no'], '0'),
                ({}, ['GKEYLL_USE_VALGRIND=1'], '1'),
                ({}, ['GKEYLL_USE_VALGRIND=1', '--use-valgrind=no'], '0'),
            )
            env = dict(os.environ)
            env.pop('GKEYLL_USE_VALGRIND', None)
            for extra_env, args, expected in cases:
                with self.subTest(script=script.name, env=extra_env, args=args):
                    subprocess.run(['sh', str(script), *args], env=env | extra_env,
                                   check=True, capture_output=True)
                    self.assertIn(f'GKEYLL_USE_VALGRIND={expected}\n', Path(output).read_text())
            saved = Path(output).read_text()
            for arg in ('--use-valgrind', '--use-valgrind=', '--use-valgrind=1',
                        '--use-valgrind=true', '--use-valgrind=maybe'):
                with self.subTest(script=script.name, arg=arg):
                    result = subprocess.run(['sh', str(script), arg], env=env,
                                            text=True, capture_output=True)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn('--use-valgrind requires yes or no', result.stdout)
                    self.assertEqual(Path(output).read_text(), saved)

    def test_configure_saves_make_profile(self):
        Path('probe.mak').write_text(f'include {ROOT}/Makefile\n'
                                    'probe:\n\t@echo $(ARCH_FLAGS) $(CFLAGS)\n')
        for value, environment in (('yes', '0'), ('no', '1')):
            with self.subTest(value=value), patch.dict(os.environ, GKEYLL_USE_VALGRIND=environment):
                subprocess.run(['sh', str(ROOT / 'configure'), f'--use-valgrind={value}'],
                               check=True, capture_output=True)
                flags = subprocess.check_output(['make', '-s', '-f', 'probe.mak',
                                                 'GIT_TIP=fixture', 'probe'], text=True)
                if value == 'yes':
                    self.assertIn(valgrind.FLAGS, flags)
                    self.assertIn('-gdwarf-4', flags)
                    self.assertNotIn('-march=native', flags)
                else:
                    self.assertIn('-march=native', flags)
                    self.assertNotIn('-gdwarf-4', flags)

    def test_openblas_default_and_valgrind_build_options(self):
        # Exercise the real dependency scripts without an expensive build.
        for script in ('mkdeps.sh', 'build-openblas.sh'):
            shutil.copy(ROOT / 'install-deps' / script, script)
        Path('bin').mkdir()
        commands = {
            'make': 'echo "$*" >> "$BUILD_CALLS"',
            'tar': 'mkdir -p OpenBLAS-0.3.30',
            'cc': 'echo x86_64-linux-gnu',
        }
        for name, body in commands.items():
            tool = Path('bin', name)
            tool.write_text('#!/bin/sh\n' + body + '\n')
            tool.chmod(0o755)
        calls = Path('calls').resolve()
        prefix = Path('deps').resolve()
        with patch.dict(os.environ, PATH=str(Path('bin').resolve()) + os.pathsep + os.environ['PATH'],
                        BUILD_CALLS=str(calls), NPROC='1', GKEYLL_USE_VALGRIND='0'):
            base = ['sh', './mkdeps.sh', f'--prefix={prefix}', '--download=no',
                    '--build-openblas=yes', 'CC=cc']
            subprocess.run(base, check=True, capture_output=True)
            self.assertNotIn('NO_AVX', calls.read_text())
            self.assertNotIn('TARGET=', calls.read_text())
            calls.unlink()
            subprocess.run(base + ['--use-valgrind=yes'], check=True, capture_output=True)
            for line in calls.read_text().splitlines():
                self.assertIn('TARGET=GENERIC', line)
                self.assertIn('NO_AVX512=1', line)
                self.assertIn('DYNAMIC_ARCH=0', line)
            calls.unlink()
            # Subsequent builds retain the selected profile despite the environment.
            subprocess.run(['bash', './build-openblas.sh'], check=True, capture_output=True)
            self.assertIn('TARGET=GENERIC', calls.read_text())
            calls.unlink()
            with patch.dict(os.environ, GKEYLL_USE_VALGRIND='1'):
                subprocess.run(base + ['--use-valgrind=no'], check=True, capture_output=True)
            self.assertNotIn('NO_AVX', calls.read_text())
            self.assertNotIn('TARGET=', calls.read_text())

    def test_success_and_stale_log_replacement(self):
        self.fixture(['clean'])
        Path('build/core/unit/ctest_clean_val_err').write_text('old errors')
        with patch.dict(os.environ, PATH=self.fake_valgrind()):
            self.assertEqual(valgrind.run(Path('ci-valgrind/candidate')), 0)
        summary = json.loads(Path('ci-valgrind/candidate/summary.json').read_text())
        self.assertEqual(summary['status'], 'passed')
        self.assertNotIn('old errors', Path('ci-valgrind/candidate/core/ctest_clean.log').read_text())

    def test_all_failures_run_and_retain_complete_logs(self):
        self.fixture(['error', 'crash', 'assertion', 'clean'])
        with patch.dict(os.environ, PATH=self.fake_valgrind()):
            self.assertEqual(valgrind.run(Path('ci-valgrind/candidate')), 1)
        summary = json.loads(Path('ci-valgrind/candidate/summary.json').read_text())
        tests = {row['name']: row for row in summary['tests']}
        self.assertEqual(len(tests), 4)
        self.assertEqual(tests['core/ctest_clean']['status'], 'passed')
        self.assertEqual(tests['core/ctest_error']['errors'], 1)
        self.assertIsNone(tests['core/ctest_crash']['errors'])
        self.assertEqual(tests['core/ctest_assertion']['status'], 'failed')
        self.assertIn('Assertion failed', Path('ci-valgrind/candidate/core/ctest_assertion.log').read_text())

    def test_make_target_propagates_memory_error(self):
        self.fixture(['error', 'clean'])
        Path('core/minus').mkdir()
        Path('core/corelinkobjs.mak').touch()
        shutil.copyfile(RUNNER, 'core/minus/valcheck.sh')
        shutil.copyfile(ROOT / 'Makefile.valgrind', 'Makefile.valgrind')
        # Use the real make recipe while skipping compilation of fixture tests.
        with patch.dict(os.environ, PATH=self.fake_valgrind()):
            result = subprocess.run(['make', '-f', str(ROOT / 'core/Makefile-core'),
                                     '-o', 'unit', 'valcheck',
                                     'BUILD_DIR=build', 'KERNELS_DIR=ker',
                                     'UNITS=../build/core/unit/ctest_error ../build/core/unit/ctest_clean'],
                                    cwd='core', capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('FAIL valgrind:', result.stdout, result.stderr)
        self.assertIn('PASS valgrind:', result.stdout)
        self.assertEqual(Path('build/core/unit/ctest_error_val_err.exit').read_text(), '99\n')

    def test_missing_tool_and_build_failure_are_not_success(self):
        self.fixture(['clean'])
        with patch.object(shutil, 'which', return_value=None):
            self.assertEqual(valgrind.run(Path('missing-tool')), 1)
        self.assertIn('not installed', Path('missing-tool/summary.json').read_text())
        Path('Makefile').write_text('core-unit:\n\t@exit 2\n')
        with patch.dict(os.environ, PATH=self.fake_valgrind()):
            self.assertEqual(valgrind.run(Path('failed-build')), 1)
        summary = json.loads(Path('failed-build/summary.json').read_text())
        self.assertNotEqual(summary['make_exit'], 0)
        self.assertIsNone(summary['tests'][0]['log'])

    def test_large_error_report_round_trips_every_character(self):
        root = Path('ci-valgrind/candidate')
        root.mkdir(parents=True)
        output = ''.join(f'==123== frame {i}: <&> λ\n' for i in range(12000))
        (root / 'ctest_leak.log').write_text(output)
        (root / 'summary.json').write_text(json.dumps({'status': 'failed', 'tests': [
            dict(name='ctest_leak', status='failed', errors=1, exit='99', log='ctest_leak.log')]}))
        sections = github_report.valgrind_sections()
        displayed = ''.join(log_text(section) for section in sections[1:])
        self.assertEqual(displayed, output)
        pages = github_report.report_pages('report', sections, 'ci/test')
        self.assertGreater(len(pages), 1)
        self.assertTrue(all(len(page.encode()) < github_report.COMMENT_LIMIT for page in pages))
        args = argparse.Namespace(platform='personal', context='ci/test', result='failure', pr='', output='ci-report.md')
        with contextlib.redirect_stdout(io.StringIO()):
            github_report.build_report(args)
        self.assertIn('Valgrind errors', Path('ci-report.md').read_text())
        self.assertIn('frame 11999', Path('ci-report.md').read_text())

    def test_incomplete_report_keeps_partial_log(self):
        root = Path('ci-valgrind/candidate')
        root.mkdir(parents=True)
        (root / 'ctest_crash.log').write_text('Valgrind crashed before its summary\n')
        sections = github_report.valgrind_sections()
        self.assertIn('incomplete', sections[0])
        self.assertIn('Valgrind crashed', ''.join(sections))

    def test_cache_separates_modes(self):
        def key(value):
            with patch.dict(os.environ, GKEYLL_USE_VALGRIND=value):
                return subprocess.check_output(['bash', str(ROOT / 'ci/jenkins/baseline_cache.sh'), 'context'])
        self.assertNotEqual(key('0'), key('1'))
        self.assertEqual(key('1'), key('1'))

    @unittest.skipUnless(shutil.which('valgrind') and shutil.which('cc'), 'requires Valgrind and a C compiler')
    def test_real_memcheck_detects_leak_invalid_access_and_assertion(self):
        cases = {
            'clean': '#include <stdlib.h>\nint main(void) { void *ptr = malloc(16); free(ptr); return 0; }\n',
            'leak': '#include <stdlib.h>\nint main(void) { void *ptr = malloc(16); return ptr == NULL; }\n',
            'invalid': '#include <stdlib.h>\nint main(void) { int *ptr = malloc(4); free(ptr); *ptr = 1; return 0; }\n',
            'assertion': '#include <assert.h>\nint main(void) { assert(0); }\n',
        }
        self.fixture(cases)
        for name, source in cases.items():
            filename = Path(f'core/unit/ctest_{name}.c')
            filename.write_text(source)
            subprocess.run(['cc', '-g', '-O0', str(filename), '-o', f'build/core/unit/ctest_{name}'], check=True)
        self.assertEqual(valgrind.run(Path('ci-valgrind/candidate')), 1)
        summary = json.loads(Path('ci-valgrind/candidate/summary.json').read_text())
        states = {row['name']: row['status'] for row in summary['tests']}
        self.assertEqual(states, {'core/ctest_' + name: 'passed' if name == 'clean' else 'failed' for name in cases})


if __name__ == '__main__':
    unittest.main()
