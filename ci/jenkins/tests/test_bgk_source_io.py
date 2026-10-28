"""Compare actual BGK source output from serial and four-rank CBC runs.

Set GKEYLL, GKEYLL_CBC, and MPIEXEC, and expose the build's shared libraries.
"""
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

GKEYLL = os.environ.get('GKEYLL')
CBC = os.environ.get('GKEYLL_CBC')
MPIEXEC = os.environ.get('MPIEXEC')


@unittest.skipUnless(GKEYLL and CBC and MPIEXEC,
                     'Set GKEYLL, GKEYLL_CBC, and MPIEXEC for the BGK MPI regression')
class BgkSourceIOTests(unittest.TestCase):
    def test_parallel_source_matches_serial(self):
        with tempfile.TemporaryDirectory(prefix='gkeyll-bgk-io-') as directory:
            root = Path(directory)
            for mode, command in (
                    ('serial', [CBC, '-s1']),
                    ('parallel', shlex.split(MPIEXEC) + ['-n', '4', CBC, '-M', '-e', '4', '-s1'])):
                work = root / mode
                work.mkdir()
                result = subprocess.run(command, cwd=work, capture_output=True,
                                        text=True, timeout=120)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            script = root / 'compare.lua'
            script.write_text('''
local root = %s
for _, species in ipairs({"ion", "elc"}) do
   for _, suffix in ipairs({"source_bgk_rate", "source_bgk_feq"}) do
      local name = "rt_gk_cbc_3x2v_p1-" .. species .. "_" .. suffix .. "_0.gkyl"
      local serialGrid, serialArray = G0.Zero.arrayNewFromFile(root .. "/serial/" .. name)
      local mpiGrid, mpiArray = G0.Zero.arrayNewFromFile(root .. "/parallel/" .. name)
      assert(serialGrid and mpiGrid, "array read failed: " .. name)
      assert(G0.Zero.rectGridCmp(serialGrid, mpiGrid), "grid mismatch: " .. name)
      local range = G0.Zero.createGridRanges(serialGrid, {0,0,0,0,0,0})
      local diff = G0.Zero.arrayDiff(serialArray, mpiArray, range)
      assert(diff.is_compatible, "incompatible arrays: " .. name)
      local maxAbs = math.max(math.abs(diff.min_abs_diff), math.abs(diff.max_abs_diff))
      assert(maxAbs <= 1e-12 or diff.max_rel_diff <= 1e-12, "data mismatch: " .. name)
   end
end
''' % json.dumps(str(root)))
            result = subprocess.run([GKEYLL, '-S', str(script)], cwd=root,
                                    capture_output=True, text=True, timeout=20)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertNotIn('*** ERROR:', result.stderr)


if __name__ == '__main__':
    unittest.main()
