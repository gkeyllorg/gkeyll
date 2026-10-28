"""Exercise the actual submission shell with a contaminated host environment."""
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
PIPELINE = (ROOT / 'ci/jenkins/jenkinsfile.stellar-ai-cpu').read_text()
SUBMIT = PIPELINE.split('def submitSlurmJob(', 1)[1].split('\ndef ', 1)[0]
SCRIPT = re.search(r"loggedSh\(.*?, '''(.*?)'''", SUBMIT, re.S).group(1)


class StellarAiBatchEnvironmentTests(unittest.TestCase):
    def test_submission_removes_host_state_and_retains_runtime_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            observed = root / 'environment.json'
            sbatch = root / 'sbatch'
            sbatch.write_text('''#!/usr/bin/python3
import json, os, sys
from pathlib import Path
Path(%r).write_text(json.dumps(dict(environment=dict(os.environ), args=sys.argv[1:])))
print('4481')
''' % str(observed))
            sbatch.chmod(0o755)
            for name, output in [('squeue', ''), ('sacct', 'COMPLETED|0:0|1')]:
                executable = root / name
                executable.write_text('#!/bin/sh\nprintf "%%s\\n" "%s"\n' % output)
                executable.chmod(0o755)
            env = dict(os.environ, PATH=str(root) + ':/usr/bin:/bin',
                       HOME=str(root), USER='fixture', LOGNAME='fixture',
                       CI_WORKSPACE=str(root), CI_SLURM_JOB_KIND='unit',
                       CI_SLURM_JOB_ID_FILE=str(root / 'id'),
                       CI_SLURM_JOB_STATUS_FILE=str(root / 'status'),
                       CI_STELLAR_AI_CPU_QOS='test', CI_STELLAR_AI_CPU_TIME='00:02:00',
                       CI_STELLAR_AI_CPU_CPUS='4', CI_SLURM_TASKS='1',
                       CI_SLURM_OUTPUT='slurm.out', CI_STELLAR_AI_CPU_ACCOUNT='csd',
                       CI_SLURM_PAYLOAD='ci/jenkins/slurm-unit-tests.stellar-ai-cpu.sh',
                       CI_SLURM_EXPORTS='CI_WORKSPACE=%s,CI_CANDIDATE_DIR=%s/candidate' % (root, root),
                       GKEYLL_CI_DEPENDENCY_ENV=str(root / 'dependencies/env.sh'),
                       BUILD_NUMBER='7', LOADEDMODULES='nvhpc/26.5:cudatoolkit/13.3',
                       _LMFILES_='/missing/gpu/modulefile', MODULEPATH='/host/modules',
                       LD_LIBRARY_PATH='/host/gpu/lib', PYTHONPATH='/host/anaconda',
                       GITHUB_TOKEN='must-not-reach-batch-job')
            # Exported Bash functions and BASH_ENV must not reach the payload.
            script = SCRIPT.replace('            job_id=',
                                    '            export BASH_ENV=/missing/host/bashrc\n'
                                    '            module() { return 1; }; export -f module\n'
                                    '            job_id=', 1)
            result = subprocess.run(['bash', '--noprofile', '--norc', '-c', script],
                                    env=env, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            data = json.loads(observed.read_text())
            batch = data['environment']
            for key in ('LOADEDMODULES', '_LMFILES_', 'MODULEPATH', 'LD_LIBRARY_PATH',
                        'PYTHONPATH', 'GITHUB_TOKEN', 'BASH_ENV', 'BASH_FUNC_module%%'):
                self.assertNotIn(key, batch)
            for key in ('HOME', 'USER', 'LOGNAME', 'GKEYLL_CI_DEPENDENCY_ENV'):
                self.assertEqual(batch[key], env[key])
            self.assertEqual(batch['PATH'], '/usr/bin:/bin')
            self.assertEqual(batch['SLURM_EXPORT_ENV'], 'ALL')
            self.assertIn('--export=ALL,' + env['CI_SLURM_EXPORTS'], data['args'])
            self.assertEqual(data['args'][-1], env['CI_SLURM_PAYLOAD'])
            self.assertEqual((root / 'id').read_text(), '4481\n')

    def test_all_payloads_skip_login_profiles(self):
        for path in (ROOT / 'ci/jenkins').glob('slurm-*.stellar-ai-cpu.sh'):
            with self.subTest(payload=path.name):
                self.assertEqual(path.read_text().splitlines()[0], '#!/bin/bash')


if __name__ == '__main__':
    unittest.main()
