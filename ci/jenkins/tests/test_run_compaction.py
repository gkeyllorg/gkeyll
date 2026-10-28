"""Verify run compaction preserves payloads and refuses unfinished runs."""
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
import unittest

CI = Path(__file__).resolve().parents[1]


class RunCompactionTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='gkeyll-compaction-')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.run = self.root / 'runs/personal/fixture'
        self.candidate = self.run / 'candidate' / ('a' * 40)
        self.candidate.mkdir(parents=True)
        self.manifest = self.candidate / 'candidate-manifest.txt'
        self.manifest.write_text('result=success\nfinished_utc=2026-10-15T00:00:00Z\n')
        (self.candidate / 'payload').write_bytes(b'\x00simulation data\xff')
        self.bin = self.root / 'bin'
        self.bin.mkdir()
        for name in ('bash', 'mkdir', 'rmdir', 'cat', 'grep', 'tar', 'gzip', 'rm', 'mv'):
            (self.bin / name).symlink_to(shutil.which(name))
        self.env = dict(os.environ, PATH=str(self.bin), BUILD_TAG='fixture')

    def compact(self):
        return subprocess.run([str(CI / 'baseline_cache.sh'), 'compact-run',
                               str(self.root), 'personal'], env=self.env,
                              capture_output=True, text=True)

    def test_round_trip_and_repeat(self):
        dependencies = self.run / 'dependencies'
        dependencies.mkdir()
        (dependencies / 'library').write_bytes(b'library')
        (self.candidate / 'link').symlink_to('payload')
        result = self.compact()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(self.candidate.exists())
        self.assertFalse(dependencies.exists())
        with tarfile.open(self.run / 'run.tar.gz') as archive:
            prefix = 'candidate/' + 'a' * 40
            self.assertEqual(archive.extractfile(prefix + '/payload').read(), b'\x00simulation data\xff')
            self.assertTrue(archive.getmember(prefix + '/link').issym())
            self.assertEqual(archive.extractfile('dependencies/library').read(), b'library')
        self.assertEqual(self.compact().returncode, 0)

    def test_unfinished_run_is_preserved(self):
        self.manifest.write_text('result=building\n')
        self.assertNotEqual(self.compact().returncode, 0)
        self.assertTrue((self.candidate / 'payload').exists())
        self.assertFalse((self.run / 'run.tar.gz').exists())

    def test_other_cache_owner_is_rejected(self):
        lock = self.root / 'baseline-cache/personal/.lock'
        lock.mkdir(parents=True)
        (lock / 'owner').write_text('another-build')
        self.assertNotEqual(self.compact().returncode, 0)
        self.assertTrue(self.candidate.exists())

    def test_compressor_failure_preserves_payload(self):
        (self.bin / 'gzip').unlink()
        (self.bin / 'gzip').write_text('#!/bin/bash\nexit 1\n')
        (self.bin / 'gzip').chmod(0o755)
        self.assertNotEqual(self.compact().returncode, 0)
        self.assertTrue((self.candidate / 'payload').exists())
        self.assertFalse((self.run / 'run.tar.gz').exists())

    def test_verification_failure_preserves_payload(self):
        real_tar = shutil.which('tar')
        (self.bin / 'tar').unlink()
        (self.bin / 'tar').write_text(
            '#!/bin/bash\nif [[ "$1" == -df ]]; then exit 1; fi\n'
            f'exec "{real_tar}" "$@"\n')
        (self.bin / 'tar').chmod(0o755)
        self.assertNotEqual(self.compact().returncode, 0)
        self.assertTrue((self.candidate / 'payload').exists())
        self.assertFalse((self.run / 'run.tar.gz').exists())

    def test_pipelines_require_explicit_opt_in(self):
        for platform in ('personal', 'stellar_cpu', 'perlmutter_gpu'):
            source = (CI / ('jenkinsfile.' + platform)).read_text()
            self.assertIn("env.GKEYLL_CI_COMPACT_RUNS == '1'", source)
            self.assertNotIn("env.GKEYLL_CI_COMPACT_RUNS != '0'", source)


if __name__ == '__main__':
    unittest.main()
