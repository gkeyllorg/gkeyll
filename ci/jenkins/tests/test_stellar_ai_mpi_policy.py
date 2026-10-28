"""Exercise the Stellar AI MPI environment and its baseline-cache fingerprint."""
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
PAYLOAD = ROOT / 'ci/jenkins/slurm-parallel-regression-tests.stellar-ai-cpu.sh'
PIPELINE = (ROOT / 'ci/jenkins/jenkinsfile.stellar-ai-cpu').read_text()
CONTEXT_FUNCTION = PIPELINE.split('def baselineCacheContext() {', 1)[1].split('\ndef ', 1)[0]
CONTEXT_SCRIPT = re.search(r"script: '''(.*?)'''", CONTEXT_FUNCTION, re.S).group(1)


class StellarAiMpiPolicyTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix='stellar-mpi-policy-')
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        (self.root / 'machines').mkdir()
        (self.root / 'machines/module_load.stellarai-cpu.sh').write_text(':\n')
        self.env = dict(os.environ, CI_WORKSPACE=str(self.root),
                        CI_BASELINE_DIR=str(self.root / 'baseline/gkeyll'),
                        CI_BASELINE_PREFIX=str(self.root / 'baseline/prefix'),
                        CI_CANDIDATE_DIR=str(self.root / 'candidate/gkeyll'),
                        CI_CANDIDATE_PREFIX=str(self.root / 'candidate/prefix'),
                        CI_TRUSTED_CHECKER=str(self.root / 'checker.lua'),
                        SLURM_MPI_TYPE='none',
                        OMPI_MCA_coll='han,hcoll',
                        OMPI_MCA_coll_hcoll_enable='1',
                        OMPI_MCA_coll_tuned_use_dynamic_rules='0',
                        OMPI_MCA_coll_tuned_allreduce_algorithm='0')
        (self.root / 'checker.lua').write_text('-- fixture\n')
        for kind in ('BASELINE', 'CANDIDATE'):
            Path(self.env[f'CI_{kind}_DIR']).mkdir(parents=True)
            binary = Path(self.env[f'CI_{kind}_PREFIX']) / 'gkeyll/bin/gkeyll'
            binary.parent.mkdir(parents=True)
            binary.write_text('''#!/bin/bash
set -eu
printf '%s %s %s %s %s\\n' "$OMPI_MCA_coll_hcoll_enable" \\
  "$OMPI_MCA_coll_tuned_use_dynamic_rules" "$OMPI_MCA_coll_tuned_allreduce_algorithm" "$OMPI_MCA_coll" "$SLURM_MPI_TYPE" \\
  >> "$CI_WORKSPACE/observed-policy.txt"
''')
            binary.chmod(0o755)

    def test_baseline_and_candidate_use_identical_fixed_collectives(self):
        for mode in ('baseline-create', 'candidate-check'):
            with self.subTest(mode=mode):
                observed = self.root / 'observed-policy.txt'
                observed.unlink(missing_ok=True)
                result = subprocess.run(['bash', str(PAYLOAD)], cwd=self.root,
                                        env=dict(self.env, CI_REGRESSION_MODE=mode),
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual(observed.read_text().splitlines(),
                                 ['0 1 3 ^hcoll,han pmix'] * (1 if mode == 'baseline-create' else 2))

    def test_runtime_payload_and_toolchain_both_affect_cache_key(self):
        payload = self.root / 'ci/jenkins' / PAYLOAD.name
        payload.parent.mkdir(parents=True)
        payload.write_text(PAYLOAD.read_text())
        for name in ('module_load', 'configure', 'mkdeps'):
            (self.root / f'machines/{name}.stellarai-cpu.sh').write_text(':\n')
        tool = self.root / 'cache-tool'
        tool.write_text('#!/bin/bash\necho "${FIXTURE_CONTEXT:-toolchain-v1}"\n')
        tool.chmod(0o755)
        env = dict(self.env, CI_BASELINE_CACHE_TOOL=str(tool),
                   CI_ARCH_FLAGS='-march=x86-64-v4 -mtune=znver5')

        def context():
            return subprocess.check_output(['bash', '-c', CONTEXT_SCRIPT],
                                           cwd=self.root, env=env, text=True).strip()

        original = context()
        self.assertEqual(original, context())
        payload.write_text(payload.read_text().replace('allreduce_algorithm=3',
                                                      'allreduce_algorithm=1'))
        self.assertNotEqual(original, context())
        payload.write_text(PAYLOAD.read_text())
        self.assertEqual(original, context())
        module = self.root / 'machines/module_load.stellarai-cpu.sh'
        module.write_text('# upgraded modules\n')
        self.assertNotEqual(original, context())
        module.write_text(':\n')
        self.assertEqual(original, context())
        env['FIXTURE_CONTEXT'] = 'toolchain-v2'
        self.assertNotEqual(original, context())
        del env['FIXTURE_CONTEXT']
        self.assertEqual(original, context())
        env['CI_ARCH_FLAGS'] = '-march=x86-64-v3 -mtune=generic'
        self.assertNotEqual(original, context())


if __name__ == '__main__':
    unittest.main()
