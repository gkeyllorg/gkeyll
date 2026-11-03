"""Run with python -m unittest discover -s gyrokinetic/data/eqdsk/miller -v."""
from dataclasses import replace
import contextlib
import io
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import RectBivariateSpline

from gkeyll_gen_miller_eqdsk import (
    MillerGeometry, MillerEquilibrium, load_input, main, sample_eqdsk,
    summary, write_eqdsk,
)

HERE = Path(__file__).resolve().parent


def read_geqdsk(path):
  """Independent fixed-width parser: checks records and coordinate pairing."""
  lines = Path(path).read_text().splitlines()
  nr, nz = int(lines[0][52:56]), int(lines[0][56:60])
  cursor = 1

  def array(count):
    """Read exactly count values from complete fixed-width GEQDSK records."""
    nonlocal cursor
    values = []
    while len(values) < count:
      line = lines[cursor]
      cursor += 1
      if len(line) % 16 or len(line) > 80:
        raise AssertionError("malformed fixed-width record")
      values.extend(float(line[i:i+16]) for i in range(0, len(line), 16))
    if len(values) != count:
      raise AssertionError("array did not end on a record boundary")
    return np.array(values)

  header = array(20)
  profiles = [array(nr) for _ in range(4)]
  psi = array(nr*nz).reshape(nz, nr)
  qpsi = array(nr)
  nboundary, nlimiter = map(int, lines[cursor].split())
  cursor += 1
  boundary = array(2*nboundary).reshape(nboundary, 2)
  if nlimiter != 0 or cursor != len(lines):
    raise AssertionError("unexpected trailing records")
  return header, profiles, psi, qpsi, boundary


def reconstruction_errors(eq, sampled):
  """Field and contour errors from an independent tensor cubic interpolant."""
  spline = RectBivariateSpline(sampled["Z"], sampled["R"], sampled["psi"])
  theta = np.linspace(-np.pi, np.pi, 1025)
  lower = float(eq.geometry.r_at_R_omp(eq.bounds[0])) if eq.bounds else 0.2*eq.a
  upper = float(eq.geometry.r_at_R_omp(eq.bounds[1])) if eq.bounds else eq.a
  errors = dict(psi=0.0, B=0.0, Bp=0.0, q=0.0)
  for radius in np.linspace(max(lower, 0.1*eq.a), upper, 7):
    R, Z = eq.geometry.RZ_at_rtheta(radius, theta)
    dpsi_R = spline.ev(Z, R, dx=0, dy=1)
    dpsi_Z = spline.ev(Z, R, dx=1, dy=0)
    Bp = np.hypot(dpsi_R, dpsi_Z)/R
    Br, Bz, Bt = eq.Bfield_at_rtheta(radius, theta)
    Bp_ref = np.hypot(Br, Bz)
    B_ref = np.hypot(Bp_ref, Bt)
    _, Rt, _, Zt = eq.geometry.RZ_derivatives_at_rtheta(radius, theta)
    q_reconstructed = eq.geometry.F/(2*np.pi)*simpson(np.hypot(Rt, Zt)/(R*R*Bp), x=theta)
    errors["psi"] = max(errors["psi"], np.max(np.abs(spline.ev(Z, R)-eq.psi_at_r(radius)))/eq.psi_boundary)
    errors["B"] = max(errors["B"], np.max(np.abs(np.hypot(Bp, Bt)/B_ref-1)))
    errors["Bp"] = max(errors["Bp"], np.max(np.abs(Bp/Bp_ref-1)))
    errors["q"] = max(errors["q"], abs(q_reconstructed/(np.sign(eq.geometry.F)*eq.q_at_r(radius))-1))
  return {key: float(value) for key, value in errors.items()}


class TestMiller(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    """Construct shared circular and LTX, DIII-D, and TCV equilibria for the tests."""
    cls.circle = MillerEquilibrium(MillerGeometry(2.0, 1.5, 2.5), lambda R: 2.0,
                                   guard_fraction=0.25)
    cls.examples = {}
    for name in ("ltx", "d3d", "tcv"):
      geometry, q, bounds, _ = load_input(HERE/f"miller_inp_{name}.py")
      cls.examples[name] = MillerEquilibrium(geometry, q, bounds)

  def test_circular_flux_and_inverse(self):
    """Check circular constant-q flux and its inverse against the analytic solution."""
    eq = self.circle
    for radius in np.linspace(0, eq.rminor_guard, 101):
      # Rationalized expression avoids cancellation near the axis.
      reference = eq.geometry.F/2*radius**2/(2+math.sqrt(4-radius**2))
      self.assertAlmostEqual(float(eq.psi_at_r(radius))/eq.psi_boundary,
                             reference/eq.psi_boundary, delta=1e-8)
      self.assertAlmostEqual(eq.r_at_psi(eq.psi_at_r(radius))/eq.rminor_guard,
                             radius/eq.rminor_guard, delta=1e-8)
    self.assertEqual(float(eq.dpsidr_at_r(0)), 0)
    with self.assertRaises(ValueError):
      eq.r_at_psi(eq.psi_guard*1.01)

  def test_signed_toroidal_field(self):
    """Changing toroidal orientation preserves flux and the poloidal field."""
    positive = self.circle
    negative = MillerEquilibrium(replace(positive.geometry, bmag_axis=-1.5),
                                 positive.q, guard_fraction=0.25)
    radius = np.linspace(0.01, positive.a, 17)
    theta = np.linspace(-np.pi, np.pi, 17)
    np.testing.assert_array_equal(positive.psi_at_r(radius), negative.psi_at_r(radius))
    bp = positive.Bfield_at_rtheta(radius, theta)
    bn = negative.Bfield_at_rtheta(radius, theta)
    np.testing.assert_array_equal(bp[:2], bn[:2])
    np.testing.assert_array_equal(bp[2], -bn[2])
    for eq in (positive, negative):
      sampled = sample_eqdsk(eq, nr=129, nz=129)
      with tempfile.TemporaryDirectory() as directory:
        path = Path(directory)/"signed.geqdsk"
        write_eqdsk(path, eq, sampled)
        header, profiles, psi, qpsi, _ = read_geqdsk(path)
      self.assertEqual(header[9], eq.geometry.bmag_axis)
      np.testing.assert_allclose(profiles[0], eq.geometry.F, rtol=1e-9)
      np.testing.assert_allclose(qpsi, np.sign(eq.geometry.F)*2, rtol=1e-9)
      self.assertGreater(header[8], header[7])

  def test_derivatives_and_coordinate_round_trips(self):
    """Check surface derivatives, coordinate inversions, and outward flux monotonicity."""
    for eq in self.examples.values():
      g = eq.geometry
      theta = np.linspace(-np.pi, np.pi, 137)
      radius = eq.a*0.81
      step = 1e-6*eq.a
      Rr, Rt, Zr, Zt = g.RZ_derivatives_at_rtheta(radius, theta)
      Rp, Zp = g.RZ_at_rtheta(radius+step, theta)
      Rm, Zm = g.RZ_at_rtheta(radius-step, theta)
      np.testing.assert_allclose(Rr, (Rp-Rm)/(2*step), atol=2e-9)
      np.testing.assert_allclose(Zr, (Zp-Zm)/(2*step), atol=2e-9)
      Rp, Zp = g.RZ_at_rtheta(radius, theta+1e-6)
      Rm, Zm = g.RZ_at_rtheta(radius, theta-1e-6)
      np.testing.assert_allclose(Rt, (Rp-Rm)/2e-6, atol=2e-9)
      np.testing.assert_allclose(Zt, (Zp-Zm)/2e-6, atol=2e-9)
      R, Z = g.RZ_at_rtheta(radius, theta)
      recovered, angle = eq.rtheta_at_RZ(R, Z)
      np.testing.assert_allclose(recovered, radius, rtol=1e-12)
      np.testing.assert_allclose(np.cos(angle), np.cos(theta), atol=1e-12)
      radii = np.linspace(0, eq.rminor_guard, 101)
      np.testing.assert_allclose(g.r_at_R_omp(g.R_omp_at_r(radii)), radii, atol=1e-14)
      psi = eq.psi_at_r(radii)
      self.assertTrue(np.all(np.diff(psi) > 0))

  def test_ltx_replacement_q_and_shear(self):
    """Verify that the replacement LTX profile preserves LCFS q and magnetic shear."""
    eq = self.examples["ltx"]
    h = eq.a*1e-5
    shear = eq.a/eq.q_at_r(eq.a)*(eq.q_at_r(eq.a+h)-eq.q_at_r(eq.a-h))/(2*h)
    self.assertAlmostEqual(eq.q_at_r(eq.a), 4.3153848, places=10)
    self.assertAlmostEqual(shear, 2.6899871, places=8)

  def test_geqdsk_format_and_determinism(self):
    """Check repeatable GEQDSK output, array ordering, profiles, and boundary metadata."""
    eq = self.examples["tcv"]
    sampled = sample_eqdsk(eq, nr=257, nz=259)
    with tempfile.TemporaryDirectory() as directory:
      path = Path(directory)/"test.geqdsk"
      write_eqdsk(path, eq, sampled)
      original = path.read_bytes()
      write_eqdsk(path, eq, sampled)
      self.assertEqual(path.read_bytes(), original)
      header, profiles, psi, qpsi, boundary = read_geqdsk(path)
    np.testing.assert_allclose(psi, sampled["psi"], rtol=1e-9)
    np.testing.assert_allclose(qpsi, sampled["qpsi"], rtol=1e-9)
    np.testing.assert_allclose(boundary, sampled["boundary"], atol=5e-10)
    np.testing.assert_allclose(profiles[0], eq.geometry.F, rtol=1e-9)
    self.assertTrue(all(np.all(profile == 0) for profile in profiles[1:]))
    self.assertAlmostEqual(header[6], eq.geometry.Zvert_axis, places=9)
    self.assertAlmostEqual(header[8], eq.psi_boundary, places=9)

  def test_device_field_reconstruction(self):
    """Check reconstructed flux, magnetic fields, and contour q for each device."""
    for name, eq in self.examples.items():
      with self.subTest(device=name):
        sampled = sample_eqdsk(eq)
        errors = reconstruction_errors(eq, sampled)
        for quantity in ("B", "Bp", "q"):
          self.assertLess(errors[quantity], 1e-3, (name, errors))
        self.assertLess(errors["psi"], 1e-5, (name, errors))

  def test_refinement(self):
    """Verify decreasing flux and poloidal-field errors as the R-Z mesh is refined."""
    # A larger physical guard isolates mesh convergence even on the 65 grid.
    g, q, bounds, _ = load_input(HERE/"miller_inp_tcv.py")
    eq = MillerEquilibrium(g, q, bounds, guard_fraction=0.25, rtol=1e-8)
    errors = [reconstruction_errors(eq, sample_eqdsk(eq, n, n)) for n in (65, 129, 257)]
    for coarse, fine in zip(errors, errors[1:]):
      self.assertLess(fine["psi"], coarse["psi"])
      self.assertLess(fine["Bp"], coarse["Bp"])

  def test_three_domain_arrangements_and_omitted_bounds(self):
    """Check coverage metadata for core, SOL, crossing, and unspecified domains."""
    for bounds in (None, (2.15, 2.3), (2.55, 2.6), (2.4, 2.6)):
      eq = MillerEquilibrium(self.circle.geometry, lambda R: 2, bounds)
      sampled = sample_eqdsk(eq)
      metadata = summary(eq, sampled)
      self.assertEqual(metadata["requested_bounds"] is None, bounds is None)
      self.assertGreater(eq.rminor_guard, eq.a)
      if bounds:
        converted = metadata["requested_bounds"]
        self.assertLess(converted[0]["psi"], converted[1]["psi"])

  def test_explicit_bounds_and_exterior(self):
    """Check explicit R-Z bounds and radial derivative continuity at the guard."""
    eq = self.circle
    data = sample_eqdsk(eq, RmajorZvert_bounds=(1.3, 2.7, -0.7, 0.7))
    self.assertTrue(np.any(data["exterior"]))
    eps = eq.rminor_guard*1e-7
    derivative_lo = float(eq.dpsidr_at_r(eq.rminor_guard-eps))
    derivative_hi = float(eq.dpsidr_at_r(eq.rminor_guard+eps))
    self.assertAlmostEqual(derivative_lo/derivative_hi, 1, delta=1e-6)

  def test_invalid_inputs(self):
    """Verify rejection of invalid geometry, q profiles, coverage, and grid settings."""
    for kwargs in (dict(delta=1), dict(kappa=0), dict(bmag_axis=0), dict(bmag_axis=math.inf), dict(bmag_axis=math.nan), dict(a_shift=10)):
      fields = dict(Rmajor_axis=2, bmag_axis=1, Rmajor_LCFS=2.5)
      fields.update(kwargs)
      with self.assertRaises(ValueError):
        MillerGeometry(**fields)
    for q in (lambda R: 0, lambda R: -1, lambda R: math.nan):
      with self.assertRaisesRegex(ValueError, "q must"):
        MillerEquilibrium(self.circle.geometry, q)
    with self.assertRaises(ValueError):
      MillerEquilibrium(self.circle.geometry, lambda R: 2, (2.5, 2.3))
    with self.assertRaisesRegex(ValueError, "clearance"):
      sample_eqdsk(self.examples["tcv"], 17, 17)
    with self.assertRaises(ValueError):
      sample_eqdsk(self.circle, 64, 65)
    with self.assertRaises(ValueError):
      sample_eqdsk(self.circle, RmajorZvert_bounds=(1.9, 2.1, -0.1, 0.1))
    eq = self.examples["ltx"]
    with self.assertRaisesRegex(ValueError, "q must"):
      MillerEquilibrium(eq.geometry,
                        lambda R: 4.3153848/(1-2.6899871*((R-eq.geometry.Rmajor_axis)/eq.a-1)),
                        eq.bounds)

  def test_cli_and_input_loading(self):
    """Check CLI generation, input provenance, and reporting of malformed inputs."""
    with tempfile.TemporaryDirectory() as directory:
      input_path = Path(directory)/"circle.py"
      input_path.write_text("geometry=dict(Rmajor_axis=2.,bmag_axis=1.,Rmajor_LCFS=2.5)\n"
                            "def q(R_omp): return 2.\n")
      output = Path(directory)/"output.geqdsk"
      with contextlib.redirect_stdout(io.StringIO()):
        self.assertEqual(main([str(input_path), "--output", str(output)]), 0)
      meta = json.loads(Path(str(output)+".json").read_text())
      self.assertEqual(meta["source"]["path"], str(input_path.resolve()))
      self.assertIsNone(meta["requested_bounds"])
      input_path.write_text("geometry={}\nq=4\n")
      with contextlib.redirect_stderr(io.StringIO()) as error:
        self.assertEqual(main([str(input_path), "--output", str(output)]), 1)
      self.assertIn(str(input_path.resolve()), error.getvalue())


if __name__ == "__main__":
  unittest.main()
