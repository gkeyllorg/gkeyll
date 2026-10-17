"""Round-trip CI artifacts without compiling or running simulations."""
import fnmatch
import io
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tarfile
import tempfile
import unittest


CI = Path(__file__).resolve().parents[1]
SHA = 'a' * 40
RESULTS = Path('gkylsoft/gkeyll-results')


class ArtifactStagingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='gkeyll-artifacts-')
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.candidate = self.root / 'runs/personal/fixture/candidate' / SHA
        self.baseline = self.root / 'baseline-cache/personal' / SHA
        self.workspace = self.root / 'workspace with spaces'
        for path in (self.candidate, self.baseline, self.workspace):
            path.mkdir(parents=True)
        # Exercise gzip even on developer machines with zstd installed.
        self.bin = self.root / 'bin'
        self.bin.mkdir()
        for name in ('bash', 'mkdir', 'cp', 'rm', 'dirname', 'find', 'tar', 'gzip', 'mv', 'wc'):
            executable = shutil.which(name)
            self.assertIsNotNone(executable, name)
            (self.bin / name).symlink_to(executable)
        self.env = dict(os.environ, PATH=str(self.bin), BUILD_TAG='fixture')

    def write(self, root, relative, data):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        return path

    def populate(self):
        outputs = {
            RESULTS / 'moments/creg-runs/space dir/-field.gkyl': b'\0candidate\xff',
            RESULTS / 'parallel-c-4/vlasov/creg-runs/line\nbreak.gkyl': b'parallel',
            RESULTS / 'gyrokinetic/creg-runs/empty.gkyl': b'',
        }
        for relative, data in outputs.items():
            self.write(self.candidate, relative, data)
            self.write(self.baseline, relative, b'baseline:' + data)
        for tree in (self.candidate, self.baseline):
            self.write(tree, RESULTS / 'moments/regressiondb', b'database')
            self.write(tree, RESULTS / 'runregression.config.lua', b'return {}')
            self.write(tree, RESULTS / 'moments/creg-runs/_rr_failures.txt', b'failure details')
            self.write(tree, RESULTS / 'moments/creg-runs/line\nbreak.txt', b'diagnostic')
            self.write(tree, Path('gkeyll/build/nested/compile.log'), b'CPU log')
            self.write(tree, Path('gkeyll/cuda-build/compile.log'), b'CUDA log')
            self.write(tree, Path('gkeyll/build/not-an-artifact.o'), b'object')
            outside = self.write(self.root, Path('outside/secret.gkyl'), b'outside')
            (tree / RESULTS / 'linked.gkyl').symlink_to(outside)
            (tree / RESULTS / 'linked.log').symlink_to(outside)
            (tree / RESULTS / 'linked-directory').symlink_to(outside.parent)
            os.mkfifo(tree / RESULTS / 'pipe.gkyl')
        return outputs

    def stage(self, baseline=True, success=True):
        command = [str(CI / 'baseline_cache.sh'), 'stage-candidate-artifacts',
                   str(self.root), 'personal', SHA, str(self.workspace)]
        if baseline:
            command.append(SHA)
        result = subprocess.run(command, env=self.env, capture_output=True, text=True, timeout=30)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
        return result

    def read_archive(self, role, extension):
        path = self.workspace / 'ci-numerical' / (role + '.tar.' + extension)
        if extension == 'zst':
            data = subprocess.check_output([shutil.which('zstd'), '-q', '-d', '-c', str(path)])
            archive = tarfile.open(fileobj=io.BytesIO(data))
        else:
            archive = tarfile.open(path)
        with archive:
            members = archive.getmembers()
            self.assertTrue(all(member.isfile() for member in members))
            return {Path(member.name): archive.extractfile(member).read() for member in members}

    def test_gzip_round_trip_diagnostics_and_independent_snapshot(self):
        outputs = self.populate()
        # Candidate-controlled destination symlinks must be replaced, not used.
        for name in ('gkylsoft', 'build', 'cuda-build', '_baseline', 'ci-numerical'):
            (self.workspace / name).symlink_to(self.root / 'outside')
        (self.candidate / '_baseline').symlink_to(self.root / 'outside')
        self.stage()
        self.assertEqual(self.read_archive('candidate', 'gz'), outputs)
        self.assertEqual(self.read_archive('baseline', 'gz'),
                         {path: b'baseline:' + data for path, data in outputs.items()})
        self.assertEqual(sorted(p.name for p in (self.root / 'outside').iterdir()), ['secret.gkyl'])
        self.assertFalse(list(self.workspace.rglob('*.gkyl')))
        for root in (self.workspace, self.workspace / '_baseline'):
            self.assertEqual((root / RESULTS / 'moments/regressiondb').read_bytes(), b'database')
            self.assertEqual((root / RESULTS / 'runregression.config.lua').read_bytes(), b'return {}')
            self.assertEqual((root / RESULTS / 'moments/creg-runs/line\nbreak.txt').read_bytes(), b'diagnostic')
            self.assertEqual((root / 'build/nested/compile.log').read_bytes(), b'CPU log')
            self.assertEqual((root / 'cuda-build/compile.log').read_bytes(), b'CUDA log')
            self.assertFalse((root / 'build/not-an-artifact.o').exists())
            self.assertFalse((root / RESULTS / 'linked.log').exists())
            self.assertFalse((root / RESULTS / 'linked-directory').exists())
        shutil.rmtree(self.baseline)
        for path, data in outputs.items():
            self.assertEqual((self.candidate / path).read_bytes(), data)
            self.assertEqual((self.candidate / '_baseline' / path).read_bytes(), b'baseline:' + data)

    @unittest.skipUnless(shutil.which('zstd'), 'zstd is not installed')
    def test_zstd_round_trip_and_restage_removes_old_archives(self):
        outputs = self.populate()
        self.stage()
        (self.bin / 'zstd').symlink_to(shutil.which('zstd'))
        self.stage()
        self.assertEqual(self.read_archive('candidate', 'zst'), outputs)
        self.assertEqual(self.read_archive('baseline', 'zst'),
                         {path: b'baseline:' + data for path, data in outputs.items()})
        self.assertFalse(list((self.workspace / 'ci-numerical').glob('*.gz')))

    def test_missing_and_empty_results_after_early_failure(self):
        self.stage(baseline=False)
        self.assertFalse(list((self.workspace / 'ci-numerical').iterdir()))
        (self.candidate / RESULTS).mkdir(parents=True)
        self.stage()
        self.assertEqual(self.read_archive('candidate', 'gz'), {})
        self.assertFalse((self.workspace / 'ci-numerical/baseline.tar.gz').exists())
        # A symlink used as the results root is also excluded.
        (self.candidate / RESULTS).rmdir()
        (self.candidate / RESULTS).symlink_to(self.baseline)
        self.stage()
        self.assertFalse(list((self.workspace / 'ci-numerical').iterdir()))

    def test_compressor_failure_does_not_publish_partial_archive(self):
        self.populate()
        self.stage()
        (self.bin / 'gzip').unlink()
        (self.bin / 'gzip').write_text('#!/bin/bash\nprintf partial\nexit 1\n')
        (self.bin / 'gzip').chmod(0o755)
        result = self.stage(success=False)
        self.assertIn('could not package numerical outputs', result.stderr)
        self.assertFalse(list((self.workspace / 'ci-numerical').iterdir()))
        self.assertTrue((self.workspace / RESULTS / 'moments/regressiondb').exists())
        self.assertTrue((self.candidate / '_baseline' / RESULTS / 'moments/regressiondb').exists())

    def test_incomplete_file_scan_does_not_publish_partial_archive(self):
        self.populate()
        (self.bin / 'find').unlink()
        (self.bin / 'find').write_text(
            '#!/bin/bash\n' + shlex.quote(shutil.which('find')) + ' "$@"\n'
            'if [[ "$1" == gkylsoft/gkeyll-results ]]; then exit 1; fi\n')
        (self.bin / 'find').chmod(0o755)
        result = self.stage(success=False)
        self.assertIn('could not package numerical outputs', result.stderr)
        self.assertFalse(list((self.workspace / 'ci-numerical').iterdir()))

    def test_pipeline_masks_include_both_formats_but_not_temporary_files(self):
        for platform in ('personal', 'stellar_cpu', 'perlmutter_gpu'):
            masks = re.findall(r"artifacts: '([^']+)'", (CI / ('jenkinsfile.' + platform)).read_text())
            patterns = ','.join(masks).split(',')
            for role in ('candidate', 'baseline'):
                for extension in ('zst', 'gz'):
                    path = 'ci-numerical/{}.tar.{}'.format(role, extension)
                    self.assertTrue(any(fnmatch.fnmatchcase(path, mask) for mask in patterns), platform)
                    self.assertFalse(any(fnmatch.fnmatchcase(path + '.tmp', mask) for mask in patterns), platform)


if __name__ == '__main__':
    unittest.main()
