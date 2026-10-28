"""Check actual DG kernel code generation with the Stellar AI CI flags."""
import ctypes
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]


class StellarAiCodegenTests(unittest.TestCase):
    def test_conf_phase_product_matches_basis_products(self):
        compiler = shutil.which('cc')
        cpuinfo = Path('/proc/cpuinfo')
        if not compiler or not cpuinfo.exists():
            self.skipTest('Requires a C compiler and an x86 Linux host')
        pipeline = (ROOT / 'ci/jenkins/jenkinsfile.stellar-ai-cpu').read_text()
        flags = shlex.split(re.search(
            r"def ciArchitectureFlags\(\)\s*\{\s*return '([^']+)'", pipeline).group(1))
        features = set(cpuinfo.read_text().split())
        required = {'avx2', 'fma', 'bmi1', 'bmi2', 'f16c', 'movbe', 'xsave'}
        if '-march=x86-64-v4' in flags:
            required.update(('avx512f', 'avx512bw', 'avx512cd', 'avx512dq', 'avx512vl'))
        if not required.issubset(features):
            self.skipTest('Host cannot execute the configured CI instruction set')

        with tempfile.TemporaryDirectory(prefix='stellar-ai-codegen-') as directory:
            library = Path(directory) / 'cross_mul.so'
            result = subprocess.run([
                compiler, '-shared', '-fPIC', '-O3', '-ffast-math', *flags,
                '-Icore/zero', '-Icore/minus', '-Icore/ker/bin_op',
                'core/ker/bin_op/binop_cross_mul_1d_3d_ser_p1.c', '-o', str(library),
            ], cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            module = ctypes.CDLL(str(library))
            multiply = module.binop_cross_mul_1d_3d_ser_p1
            multiply.argtypes = [ctypes.POINTER(ctypes.c_double)] * 3
            multiply.restype = None

            # Multiplication by the constant mode preserves each phase mode.
            # Multiplication by the x mode exchanges 1<->x, y<->xy,
            # z<->xz, and yz<->xyz in the orthonormal p1 basis.
            x_partner = [1, 0, 4, 5, 2, 3, 7, 6]
            for conf_mode in range(2):
                for phase_mode in range(8):
                    with self.subTest(conf_mode=conf_mode, phase_mode=phase_mode):
                        conf = (ctypes.c_double * 2)()
                        phase = (ctypes.c_double * 8)()
                        actual = (ctypes.c_double * 8)()
                        conf[conf_mode] = phase[phase_mode] = 1.0
                        multiply(conf, phase, actual)
                        target = phase_mode if conf_mode == 0 else x_partner[phase_mode]
                        for mode in range(8):
                            expected = 2 ** -0.5 if mode == target else 0.0
                            self.assertAlmostEqual(actual[mode], expected, delta=1e-14)


if __name__ == '__main__':
    unittest.main()
