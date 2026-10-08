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
        self.adas = self.original / 'gkeyll/share/adas'
        self.adas.mkdir(parents=True)
        for name in ('ioniz_h.npy', 'recomb_h.npy', 'logT_h.npy', 'logN_h.npy'):
            (self.adas / name).write_bytes(b'\x00ADAS fixture\xff')
        (self.adas / 'radiation_fit_parameters.txt').write_text('original radiation fits')
        (self.original / 'gkeyll-results').mkdir()
        (self.original / 'gkeyll-results/old-result.gkyl').write_text('old result')
        self.config = self.root / 'source-config.mak'
        self.config.write_text(
            f'PREFIX={self.original}\nCC=cc\nBUILD_APP=pkpm\n'
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

    def test_adas_snapshot_and_build_trees_are_isolated(self):
        # Both an installation alias and file symlinks must be dereferenced.
        external = self.root / 'external-adas'
        self.adas.rename(external)
        self.adas.symlink_to(external, target_is_directory=True)
        table = external / 'ioniz_h.npy'
        table.rename(external / 'table.npy')
        table.symlink_to(external / 'table.npy')
        self.prepare()
        manifest = json.loads((self.dependencies / 'manifest.json').read_text())
        snapshot = Path(manifest['adas_dir'])
        self.assertTrue(snapshot.is_relative_to(self.dependencies))
        candidate = self.prefix / 'gkeyll/share/adas'
        for directory in (snapshot, candidate):
            self.assertFalse(directory.is_symlink())
            self.assertFalse((directory / 'ioniz_h.npy').is_symlink())
            self.assertEqual((directory / 'ioniz_h.npy').read_bytes(), table.read_bytes())
        (candidate / 'ioniz_h.npy').write_bytes(b'candidate edit')
        (candidate / 'radiation_fit_parameters.txt').write_text('candidate radiation fits')
        # Reusing the run must not need the original installation or data.
        self.original.rename(self.root / 'original-unavailable')
        external.rename(self.root / 'adas-unavailable')
        baseline_prefix = self.run / 'baseline/gkylsoft'
        write_config(self.config, baseline_prefix, self.dependencies, self.output)
        baseline = baseline_prefix / 'gkeyll/share/adas'
        self.assertEqual((baseline / 'ioniz_h.npy').read_bytes(), b'\x00ADAS fixture\xff')
        (baseline / 'radiation_fit_parameters.txt').write_text('baseline radiation fits')
        self.assertEqual((candidate / 'radiation_fit_parameters.txt').read_text(),
                         'candidate radiation fits')
        self.assertEqual((snapshot / 'radiation_fit_parameters.txt').read_text(),
                         'original radiation fits')
        self.assertEqual((snapshot / 'ioniz_h.npy').read_bytes(), b'\x00ADAS fixture\xff')
        (self.root / 'original-unavailable').rename(self.original)
        (self.root / 'adas-unavailable').rename(external)
        self.assertEqual(table.read_bytes(), b'\x00ADAS fixture\xff')
        table.write_bytes(b'updated source table')
        second_prefix = self.root / 'run-2/candidate/gkylsoft'
        write_config(self.config, second_prefix, self.root / 'run-2/dependencies', self.output)
        self.assertEqual((second_prefix / 'gkeyll/share/adas/ioniz_h.npy').read_bytes(),
                         b'updated source table')
        self.assertEqual((snapshot / 'ioniz_h.npy').read_bytes(), b'\x00ADAS fixture\xff')

    def test_adas_copied_when_library_root_is_entire_prefix(self):
        (self.original / 'include').mkdir()
        (self.original / 'lib').mkdir()
        (self.original / 'gkeyll/bin').mkdir()
        (self.original / 'gkeyll/bin/gkeyll').write_text('old executable')
        with self.config.open('a') as config:
            config.write('LAPACK_INC_DIR=$(PREFIX)/include\nLAPACK_LIB_DIR=$(PREFIX)/lib\n')
        self.prepare()
        self.assertFalse((self.dependencies / 'gkylsoft/gkeyll').exists())
        self.assertFalse((self.dependencies / 'gkylsoft/gkeyll-results').exists())
        self.assertFalse((self.prefix / 'gkeyll/bin/gkeyll').exists())
        self.assertEqual((self.prefix / 'gkeyll/share/adas/recomb_h.npy').read_bytes(),
                         (self.adas / 'recomb_h.npy').read_bytes())

    def test_missing_adas_fails_before_building_for_gyrokinetic_and_pkpm(self):
        for table in self.adas.glob('*.npy'):
            table.unlink()
        for app in ('gyrokinetic', 'pkpm'):
            with self.subTest(app=app):
                with self.config.open('a') as config:
                    config.write(f'BUILD_APP={app}\n')
                with self.assertRaisesRegex(ValueError, 'Missing ADAS .npy data'):
                    self.prepare()
                self.assertFalse(self.output.exists())
                self.assertFalse(self.dependencies.exists())

    def test_non_gyrokinetic_build_does_not_require_adas(self):
        shutil.rmtree(self.adas)
        with self.config.open('a') as config:
            config.write('BUILD_APP=vlasov\n')
        self.prepare()
        self.assertFalse((self.prefix / 'gkeyll/share/adas').exists())

    def test_old_snapshot_requires_new_run(self):
        self.prepare()
        path = self.dependencies / 'manifest.json'
        manifest = json.loads(path.read_text())
        del manifest['adas_dir']
        path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, 'start a new CI run'):
            self.prepare()

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
            [sys.executable, '-I', str(Path(__file__).resolve().parents[1] / 'prebuilt_config.py'),
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
