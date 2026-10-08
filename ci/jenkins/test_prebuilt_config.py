"""Check private dependency copies and Make configuration without simulations."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from ci.jenkins.prebuilt_config import evaluate_config, write_config


class PrebuiltConfigTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory(prefix='gkeyll-prebuilt-')
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.original = self.root / 'original'
        self.library = self.original / 'OpenBLAS/lib'
        self.library.mkdir(parents=True)
        (self.original / 'OpenBLAS/include').mkdir()
        (self.library / 'libfixture.so.1').write_text('original library')
        (self.library / 'libfixture.so').symlink_to(self.library / 'libfixture.so.1')
        (self.original / 'gkeyll-results').mkdir()
        (self.original / 'gkeyll-results/old-result.gkyl').write_text('old result')
        self.config = self.root / 'source-config.mak'
        self.config.write_text(
            f'PREFIX={self.original}\nCC=cc\n'
            'LAPACK_INC_DIR=$(PREFIX)/OpenBLAS/include\n'
            'LAPACK_LIB_DIR=${PREFIX}/OpenBLAS/lib\n')
        self.run = self.root / 'run-1'
        self.run.mkdir()
        self.output = self.run / 'config.mak'
        self.dependencies = self.run / 'dependencies'
        self.prefix = self.run / 'candidate/gkylsoft'

    def prepare(self):
        return write_config(self.config, self.prefix, self.dependencies, self.output)

    def test_copy_is_private_and_reused_only_within_run(self):
        self.prepare()
        values = evaluate_config(self.output)
        self.assertEqual(values['PREFIX'], str(self.prefix))
        copied = Path(values['LAPACK_LIB_DIR']) / 'libfixture.so'
        self.assertTrue(copied.is_relative_to(self.dependencies))
        self.assertFalse(copied.is_symlink())
        self.assertEqual(copied.read_text(), 'original library')
        copied.write_text('private modification')
        self.assertEqual((self.library / 'libfixture.so').read_text(), 'original library')
        self.assertFalse((self.dependencies / 'gkylsoft/gkeyll-results').exists())
        write_config(self.config, self.run / 'baseline/gkylsoft', self.dependencies, self.output)
        self.assertEqual(copied.read_text(), 'private modification')
        second = self.root / 'run-2/dependencies'
        write_config(self.config, self.root / 'run-2/gkylsoft', second, self.output)
        self.assertEqual((second / 'gkylsoft/OpenBLAS/lib/libfixture.so').read_text(),
                         'original library')

    def test_external_dependency_and_source_config_preserved(self):
        external = self.root / 'external/superlu'
        (external / 'include').mkdir(parents=True)
        (external / 'lib').mkdir()
        (external / 'lib/libsuperlu.so').write_text('external library')
        with self.config.open('a') as config:
            config.write(f'SUPERLU_INC_DIR={external}/include\nSUPERLU_LIB_DIR={external}/lib\n')
        original = self.config.read_bytes()
        self.prepare()
        copied_dir = Path(evaluate_config(self.output)['SUPERLU_LIB_DIR'])
        self.assertTrue(copied_dir.is_relative_to(self.dependencies))
        self.assertEqual((copied_dir / 'libsuperlu.so').read_text(), 'external library')
        self.assertEqual(self.config.read_bytes(), original)

    def test_config_change_cannot_mix_candidate_and_baseline_dependencies(self):
        self.prepare()
        self.config.write_text(self.config.read_text() + '\nCC=gcc\n')
        with self.assertRaisesRegex(ValueError, 'changed during this run'):
            self.prepare()

    def test_missing_dependency_fails_before_building(self):
        shutil.rmtree(self.library)
        with self.assertRaisesRegex(ValueError, 'LAPACK_LIB_DIR'):
            self.prepare()
        self.assertFalse(self.output.exists())

    def test_system_libraries_remain_machine_provided(self):
        self.config.write_text('PREFIX=/usr\nLAPACK_INC_DIR=/usr/include\n'
                               'LAPACK_LIB_DIR=/usr/lib\n')
        self.prepare()
        values = evaluate_config(self.output)
        self.assertEqual(values['PREFIX'], str(self.prefix))
        self.assertEqual(values['LAPACK_INC_DIR'], '/usr/include')
        self.assertEqual(values['LAPACK_LIB_DIR'], '/usr/lib')

    def run_cli(self, *options, library_path=''):
        return subprocess.run(
            [sys.executable, '-I', str(Path(__file__).with_name('prebuilt_config.py')),
             '--config', str(self.config), '--prefix', str(self.prefix),
             '--dependencies', str(self.dependencies), '--output', str(self.output), *options],
            env=dict(os.environ, LD_LIBRARY_PATH=library_path),
            capture_output=True, text=True)

    def test_cli_env_matches_json_with_empty_and_literal_values(self):
        library_path = '/existing/lib=with spaces:/literal/"quotes"/and\\backslash'
        default = self.run_cli(library_path=library_path)
        self.assertEqual(default.returncode, 0, default.stderr)
        expected = json.loads(default.stdout)
        self.assertEqual(expected['MPI_HOME'], '')
        self.assertEqual(expected['OPAL_PREFIX'], '')
        self.assertEqual(expected['MPIEXEC'], '')
        self.assertTrue(expected['LD_LIBRARY_PATH'].endswith(':' + library_path))
        result = self.run_cli('--format', 'env', library_path=library_path)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(dict(line.split('=', 1) for line in result.stdout.splitlines()), expected)

    def test_cli_env_rejects_line_breaks_without_partial_output(self):
        for line_break in ('\n', '\r'):
            with self.subTest(line_break=line_break):
                result = self.run_cli('--format', 'env',
                                      library_path='/existing/lib' + line_break + 'MPI_HOME=/wrong')
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('must not contain line breaks', result.stderr)
                self.assertEqual(result.stdout, '')

    @unittest.skipUnless(shutil.which('cc'), 'A C compiler is needed for the runtime fixture')
    def test_runtime_uses_copied_transitive_libraries(self):
        # The first shared library has an absolute RUNPATH back to the original
        # install. The saved environment must still load its dependency locally.
        second = self.root / 'second.c'
        first = self.root / 'first.c'
        main = self.root / 'main.c'
        second.write_text('int second(void) { return 42; }\n')
        first.write_text('int second(void); int first(void) { return second(); }\n')
        main.write_text('int first(void); int main(void) { return first() != 42; }\n')
        subprocess.run(['cc', '-shared', '-fPIC', str(second), '-o',
                        str(self.library / 'libsecond.so')], check=True)
        subprocess.run(['cc', '-shared', '-fPIC', str(first), '-o',
                        str(self.library / 'libfirst.so'), '-L' + str(self.library),
                        '-lsecond', '-Wl,-rpath,' + str(self.library)], check=True)
        runtime = self.prepare()
        copied_dir = evaluate_config(self.output)['LAPACK_LIB_DIR']
        executable = self.run / 'fixture'
        subprocess.run(['cc', str(main), '-o', str(executable), '-L' + copied_dir,
                        '-lfirst', '-Wl,-rpath,' + copied_dir], check=True)
        self.original.rename(self.root / 'original-unavailable')
        subprocess.run([str(executable)], env=dict(os.environ, **runtime), check=True)
        subprocess.run(['bash', '-c', '. "$1"; "$2"', 'fixture',
                        str(self.dependencies / 'env.sh'), str(executable)], check=True)


if __name__ == '__main__':
    unittest.main()
