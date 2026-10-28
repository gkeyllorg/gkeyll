"""Check field-line derivatives in the actual CBC and Rosenbluth-Hinton inputs.

Set GKEYLL to a matching CPU installation and expose its shared libraries.
These small compiled fixtures evaluate geometry only; they do not run a solver.
"""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
GKEYLL = os.environ.get('GKEYLL')


@unittest.skipUnless(GKEYLL, 'Set GKEYLL to a matching CPU installation')
class PeriodicGeometryInputsTests(unittest.TestCase):
    def test_field_line_derivative_across_poloidal_boundary(self):
        for name in ('rt_gk_cbc_2x2v_p1', 'rt_gk_rosenbluth_hinton_2x2v_p1'):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                work = Path(directory)
                shutil.copy(Path(GKEYLL).resolve().parents[1] / 'share/Makefile', work)
                shutil.copy(ROOT / 'gyrokinetic/creg/rt_arg_parse.h', work)
                source = ROOT / 'gyrokinetic/creg' / (name + '.c')
                # The inputs tabulate their radial integral; two nearby entries suffice
                # to test the poloidal derivative at a fixed radius.
                setup = '''
  double radii[] = {ctx.r0 - 0.01, ctx.r0 + 0.01};
  double integrals[] = {intdPsidr(radii[0], &ctx), intdPsidr(radii[1], &ctx)};
  ctx.psi_lut_size = 2;
  ctx.r_lut = radii;
  ctx.dPsidr_int_lut = integrals;
'''
                (work / 'boundary.c').write_text(
                    '#define main input_main\n#include "' + str(source) + '"\n#undef main\n'
                    'int main(void) {\n  struct gk_app_ctx ctx = create_ctx();\n' + setup + r'''
  double h = 1e-4;
  double angles[] = {0., -M_PI-h, -M_PI, -M_PI+h, M_PI-h, M_PI, M_PI+h};
  for (int i = 0; i < sizeof(angles)/sizeof(angles[0]); ++i) {
    double theta = angles[i], r = ctx.r0;
    double radius = R_rtheta(r, theta, &ctx);
    // Fundamental theorem of calculus applied to the field-line integral.
    double expected = -Bphi(radius, &ctx) * Jr(r, theta, &ctx)
      / (radius * dPsidr(r, theta, &ctx));
    double actual = (alpha(r, theta+h, 0., &ctx) - alpha(r, theta-h, 0., &ctx))/(2*h);
    if (!(fabs(actual-expected) < 1e-3*fabs(expected))) {
      fprintf(stderr, "theta=%g: derivative %g, expected %g\n", theta, actual, expected);
      return 1;
    }
  }
  return 0;
}
''')
                build = subprocess.run(['make', '-j1', 'boundary'], cwd=work,
                                       capture_output=True, text=True, timeout=60)
                self.assertEqual(build.returncode, 0, build.stdout + build.stderr)
                result = subprocess.run([str(work / 'boundary')], cwd=work,
                                        capture_output=True, text=True, timeout=20)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == '__main__':
    unittest.main()
