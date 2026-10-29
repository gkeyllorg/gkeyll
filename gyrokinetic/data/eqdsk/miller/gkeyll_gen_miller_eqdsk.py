#!/usr/bin/env python3
"""Generate a synthetic Miller GEQDSK from a Python equilibrium input file.

See README.md and doc/notes/miller_eq/gkeyll_gen_miller_eqdsk.tex. All lengths are
meters, bmag_axis is the signed toroidal field in tesla, and psi is poloidal flux / (2*pi).
"""

import argparse
from dataclasses import asdict, dataclass
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import runpy
import sys
import tempfile

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class MillerGeometry:
  """Describe a family of axisymmetric Miller flux surfaces and its toroidal field.

  Rmajor_axis and Zvert_axis locate the magnetic axis [m]; bmag_axis is the toroidal field
  there [T], signed in the cylindrical toroidal direction. Rmajor_LCFS specifies the LCFS by its outboard midplane major radius [m].
  The dimensionless elongation kappa and triangularity delta are constant across
  surfaces. The dimensionless a_shift sets the radial displacement of each
  surface center: R_center(r) = Rmajor_axis - a_shift*r**2/(2*Rmajor_axis).

  Surfaces are labeled by minor radius r [m] and poloidal angle theta
  [rad], with theta=0 at the outboard midplane. Methods evaluate R(r,theta),
  Z(r,theta), their derivatives, the Jacobian factor, and conversions between
  r and outboard major radius. The toroidal field is B_phi=Rmajor_axis*bmag_axis/R.
  MillerEquilibrium combines this immutable geometry with a q profile to
  construct the poloidal flux and magnetic field.
  """

  Rmajor_axis: float
  bmag_axis: float
  Rmajor_LCFS: float
  Zvert_axis: float = 0.0
  kappa: float = 1.0
  delta: float = 0.0
  a_shift: float = 0.0

  def __post_init__(self):
    """Validate geometry parameters and the increasing outboard LCFS branch."""
    if not all(np.isfinite(v) for v in asdict(self).values()):
      raise ValueError("all geometry parameters must be finite")
    if self.Rmajor_axis <= 0 or self.bmag_axis == 0 or self.kappa <= 0:
      raise ValueError("Rmajor_axis and kappa must be positive; bmag_axis must be nonzero")
    if self.Rmajor_LCFS <= self.Rmajor_axis or abs(self.delta) >= 1:
      raise ValueError("require Rmajor_LCFS > Rmajor_axis and abs(delta) < 1")
    self.r_at_R_omp(self.Rmajor_LCFS)

  @property
  def F(self):
    """Return the constant toroidal field function F = R*B_phi [T m]."""
    return self.Rmajor_axis * self.bmag_axis

  def R_omp_at_r(self, r):
    """Map minor radius r to outboard midplane major radius [m]."""
    r = np.asarray(r)
    return self.Rmajor_axis - self.a_shift * r*r / (2*self.Rmajor_axis) + r

  def r_at_R_omp(self, R_omp):
    """Invert the increasing outboard branch, stably also at zero shift."""
    offset = np.asarray(R_omp) - self.Rmajor_axis
    disc = 1 - 2*self.a_shift*offset/self.Rmajor_axis
    if np.any(~np.isfinite(offset)) or np.any(offset < 0) or np.any(disc <= 0):
      raise ValueError("R_omp must lie on the increasing branch above Rmajor_axis")
    return 2*offset / (1 + np.sqrt(disc))

  def RZ_at_rtheta(self, r, theta):
    """Return cylindrical R,Z [m] at minor radius r and poloidal angle theta."""
    r, theta = np.asarray(r), np.asarray(theta)
    angle = theta + math.asin(self.delta)*np.sin(theta)
    return (self.Rmajor_axis - self.a_shift*r*r/(2*self.Rmajor_axis) + r*np.cos(angle),
            self.Zvert_axis + self.kappa*r*np.sin(theta))

  def RZ_derivatives_at_rtheta(self, r, theta):
    """Return R_r, R_theta, Z_r, Z_theta; arrays broadcast normally."""
    r, theta = np.asarray(r), np.asarray(theta)
    angle = theta + math.asin(self.delta)*np.sin(theta)
    return (-self.a_shift*r/self.Rmajor_axis + np.cos(angle),
            -r*(1 + math.asin(self.delta)*np.cos(theta))*np.sin(angle),
            self.kappa*np.sin(theta), self.kappa*r*np.cos(theta))

  def RZdet_over_r(self, r, theta):
    """(R_r Z_theta - R_theta Z_r)/r, including its regular axis limit."""
    angle = theta + math.asin(self.delta)*np.sin(theta)
    return self.kappa*((np.cos(angle) - self.a_shift*r/self.Rmajor_axis)*np.cos(theta)
                       + (1 + math.asin(self.delta)*np.cos(theta))
                       * np.sin(angle)*np.sin(theta))


def load_input(filename):
  """Execute a Python input and return geometry, scalar q, bounds, provenance."""
  path = Path(filename).resolve()
  try:
    namespace = runpy.run_path(str(path))
    geometry = MillerGeometry(**namespace["geometry"])
    q = namespace["q"]
    if not callable(q):
      raise ValueError("q must be callable as q(R_omp)")
    bounds = namespace.get("Rmajor_omp_bounds")
    if bounds is not None:
      bounds = tuple(float(v) for v in bounds)
    source = {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    return geometry, q, bounds, source
  except Exception as exc:
    raise ValueError(f"{path}: {exc}") from exc


class MillerEquilibrium:
  """Physical flux table plus explicit exterior continuation.

  q is a scalar callback q(R_omp). Rmajor_omp_bounds labels optional simulation
  surfaces; the physical construction always includes the axis and LCFS.
  """

  def __init__(self, geometry, q, Rmajor_omp_bounds=None, *, rtol=1e-10, guard_fraction=0.05):
    """Validate inputs, choose the physical guard, and construct the radial flux table."""
    self.geometry, self.q = geometry, q
    if not callable(q):
      raise ValueError("q must be callable")
    if not np.isfinite(rtol) or not 1e-13 <= rtol <= 1e-3:
      raise ValueError("rtol must be between 1e-13 and 1e-3")
    if not np.isfinite(guard_fraction) or guard_fraction <= 0:
      raise ValueError("guard_fraction must be finite and positive")
    self.rtol, self.guard_fraction = rtol, guard_fraction
    self.a = float(geometry.r_at_R_omp(geometry.Rmajor_LCFS))
    self.bounds = None
    self.rminor_required = self.a
    if Rmajor_omp_bounds is not None:
      bounds = np.asarray(Rmajor_omp_bounds, dtype=float)
      if bounds.shape != (2,) or not np.all(np.isfinite(bounds)) or bounds[0] >= bounds[1]:
        raise ValueError("Rmajor_omp_bounds must contain two finite increasing radii")
      r_bounds = geometry.r_at_R_omp(bounds)
      self.bounds = tuple(bounds)
      self.rminor_required = max(self.a, float(r_bounds[1]))
    self.rminor_guard = self.rminor_required*(1 + guard_fraction)
    self.s_guard = self.rminor_guard**2
    self._validate_surfaces()
    self._build_flux_table()
    self.psi_boundary = float(self.psi_at_r(self.a))

  def q_at_r(self, r):
    """Evaluate q at minor radius r and require a finite positive value."""
    R_omp = float(self.geometry.R_omp_at_r(r))
    try:
      value = float(self.q(R_omp))
    except Exception as exc:
      raise ValueError(f"q failed at R_omp={R_omp:.16g}: {exc}") from exc
    if not math.isfinite(value) or value <= 0:
      raise ValueError(f"q must be finite and positive at R_omp={R_omp:.16g}; got {value}")
    return value

  def _validate_surfaces(self):
    """Check sampled surfaces for folds, positive R, and valid q through the guard."""
    g = self.geometry
    theta = np.linspace(0, 2*np.pi, 1025)
    for radius in np.linspace(0, self.rminor_guard, 257):
      R, _ = g.RZ_at_rtheta(radius, theta)
      det = g.RZdet_over_r(radius, theta)
      if np.min(R) <= 0 or np.min(det) <= 0:
        raise ValueError(f"surfaces cross R=0 or fold at minor radius r={radius:.16g}")
      if 1-g.a_shift*radius/g.Rmajor_axis <= 0:
        raise ValueError("outboard midplane mapping turns before the guard surface")
      self.q_at_r(radius)

  def _flux_derivative(self, s):
    """Direct quadrature for dpsi/ds; independent of the radial interpolant."""
    g, radius = self.geometry, math.sqrt(s)
    triangularity = math.asin(g.delta)

    def integrand(theta):
      """Return the dimensionless poloidal integrand for the flux derivative."""
      sint, cost = math.sin(theta), math.cos(theta)
      angle = theta + triangularity*sint
      cosa, sina = math.cos(angle), math.sin(angle)
      R = g.Rmajor_axis - g.a_shift*s/(2*g.Rmajor_axis) + radius*cosa
      det = g.kappa*((cosa-g.a_shift*radius/g.Rmajor_axis)*cost
                    + (1+triangularity*cost)*sina*sint)
      return det*g.Rmajor_axis/R  # Dimensionless; finite at the axis.

    integral, error = quad(integrand, 0, 2*math.pi, epsabs=self.rtol*0.05,
                           epsrel=self.rtol*0.05, limit=200)
    if error > self.rtol*max(1, abs(integral)):
      raise ValueError(f"poloidal quadrature did not converge at r={radius:.16g}")
    value = abs(g.F)*integral/(4*math.pi*g.Rmajor_axis*self.q_at_r(radius))
    if not math.isfinite(value) or value <= 0:
      raise ValueError("flux derivative must be finite and positive")
    return value

  def _build_flux_table(self):
    """Adaptively tabulate dpsi/ds for s=r^2 and integrate its PCHIP interpolant."""
    # Triangular surfaces can produce sqrt(s) terms at the axis. Cluster
    # initial samples there and cache unchanged midpoint quadratures.
    evaluate = lru_cache(maxsize=None)(self._flux_derivative)
    nodes = self.s_guard*np.linspace(0, 1, 65)**4
    values = np.array([evaluate(s) for s in nodes])
    for _ in range(120):
      derivative = PchipInterpolator(nodes, values, extrapolate=False)
      mids = (nodes[:-1]+nodes[1:])/2
      exact = np.array([evaluate(s) for s in mids])
      scale = np.maximum(np.abs(exact), np.max(values)*1e-8)
      errors = np.abs(derivative(mids)-exact)/scale
      bad = errors > self.rtol
      if not np.any(bad):
        self.s_nodes, self.g_nodes = nodes, values
        self._derivative = derivative
        self._integral = derivative.antiderivative()
        self.radial_error = float(np.max(errors))
        self.psi_guard = float(self._integral(self.s_guard))
        return
      if len(nodes)+np.count_nonzero(bad) > 262145:
        break
      order = np.argsort(np.concatenate((nodes, mids[bad])))
      nodes = np.concatenate((nodes, mids[bad]))[order]
      values = np.concatenate((values, exact[bad]))[order]
    raise ValueError(f"radial flux interpolation did not converge ({len(nodes)} nodes, "
                     f"relative error {np.max(errors):.3g}); check q smoothness or rtol")

  def psi_at_r(self, r):
    """Flux at minor radius, including the documented linear-in-r² exterior."""
    radius = np.asarray(r)
    if np.any(~np.isfinite(radius)) or np.any(radius < 0):
      raise ValueError("Minor radius must be finite and nonnegative")
    s = radius**2
    return self._integral(np.minimum(s, self.s_guard)) + self.g_nodes[-1]*np.maximum(s-self.s_guard, 0)

  def dpsidr_at_r(self, r):
    """Return the radial flux derivative, including the exterior continuation."""
    radius = np.asarray(r)
    return 2*radius*self._derivative(np.minimum(radius**2, self.s_guard))

  def r_at_psi(self, psi):
    """Invert physical flux; do not silently return a background surface."""
    value = float(psi)
    if not math.isfinite(value) or not 0 <= value <= self.psi_guard:
      raise ValueError("psi is outside the physical axis-to-guard interval")
    s = brentq(lambda s: float(self._integral(s))-value, 0, self.s_guard,
               xtol=self.s_guard*1e-14, rtol=1e-14)
    return math.sqrt(s)

  def fs_coords_at_R_omp(self, R_omp):
    """Return r, psi, normalized flux radius, and q for a physical R_omp surface."""
    radius = float(self.geometry.r_at_R_omp(R_omp))
    if radius > self.rminor_guard:
      raise ValueError("R_omp exceeds the physical guard surface")
    flux = float(self.psi_at_r(radius))
    return {"R_omp": float(R_omp), "r": radius, "psi": flux,
            "rho_psi": math.sqrt(flux/self.psi_boundary), "q": self.q_at_r(radius)}

  def _R_uv(self, u, v):
    """Evaluate major radius using u=r*cos(theta) and v=r*sin(theta)."""
    radius = np.hypot(u, v)
    ratio = np.divide(v, radius, out=np.zeros_like(radius), where=radius > 0)
    angle = math.asin(self.geometry.delta)*ratio
    return (self.geometry.Rmajor_axis-self.geometry.a_shift*radius**2/(2*self.geometry.Rmajor_axis)
            + u*np.cos(angle)-v*np.sin(angle))

  def _dR_du(self, u, v):
    """Evaluate the analytic derivative of R(u,v) with respect to u at fixed v."""
    radius = np.hypot(u, v)
    ratio = np.divide(v, radius, out=np.zeros_like(radius), where=radius > 0)
    angle = math.asin(self.geometry.delta)*ratio
    angle_deriv = np.divide(-math.asin(self.geometry.delta)*v*u, radius**3,
                           out=np.zeros_like(radius), where=radius > 0)
    return (np.cos(angle)-self.geometry.a_shift*u/self.geometry.Rmajor_axis
            -(u*np.sin(angle)+v*np.cos(angle))*angle_deriv)

  def rtheta_at_RZ(self, R, Z):
    """Return r, theta by vectorized bracketed bisection of the analytic map."""
    R, Z = np.broadcast_arrays(np.asarray(R, dtype=float), np.asarray(Z, dtype=float))
    if np.any(~np.isfinite(R)) or np.any(~np.isfinite(Z)) or np.any(R <= 0):
      raise ValueError("R,Z must be finite and R positive")
    v = (Z-self.geometry.Zvert_axis)/self.geometry.kappa
    bound = np.full(R.shape, self.rminor_guard)
    for _ in range(12):
      needs = (self._R_uv(-bound, v) > R) | (self._R_uv(bound, v) < R)
      if not np.any(needs):
        break
      bound = np.where(needs, 1.5*bound, bound)
    else:
      raise ValueError("cannot bracket R–Z inverse; reduce the rectangular bounds")
    for fraction in np.linspace(-1, 1, 65):
      if np.any(self._dR_du(fraction*bound, v) <= 0):
        raise ValueError("R–Z rectangle reaches a folded inverse branch; reduce its bounds")
    lo, hi = -bound, bound
    for _ in range(50):
      mid = (lo+hi)/2
      left = self._R_uv(mid, v) < R
      lo, hi = np.where(left, mid, lo), np.where(left, hi, mid)
    u = (lo+hi)/2
    radius, theta = np.hypot(u, v), np.arctan2(v, u)
    forward_R, forward_Z = self.geometry.RZ_at_rtheta(radius, theta)
    residual = np.hypot(forward_R-R, forward_Z-Z)
    if np.any(residual > self.geometry.Rmajor_axis*1e-11):
      index = np.unravel_index(np.argmax(residual), R.shape)
      raise ValueError(f"R–Z inverse failed at ({R[index]}, {Z[index]})")
    return radius, theta

  def Bfield_at_rtheta(self, r, theta):
    """Return cylindrical (B_R, B_Z, B_phi), with B_p tangent to +theta."""
    R, _ = self.geometry.RZ_at_rtheta(r, theta)
    _, Rt, _, Zt = self.geometry.RZ_derivatives_at_rtheta(r, theta)
    det = np.asarray(r)*self.geometry.RZdet_over_r(r, theta)
    factor = np.divide(self.dpsidr_at_r(r), R*det, out=np.zeros_like(R), where=det != 0)
    return factor*Rt, factor*Zt, self.geometry.F/R


def sample_eqdsk(equilibrium, nr=257, nz=257, RmajorZvert_bounds=None):
  """Sample psi on a rectangular mesh and validate the certified surfaces."""
  if any(not isinstance(n, (int, np.integer)) or n < 5 or n % 2 == 0 or n > 999 for n in (nr, nz)):
    raise ValueError("nr and nz must be odd integers between 5 and 999")
  eq, g = equilibrium, equilibrium.geometry
  theta = np.linspace(-np.pi, np.pi, 4097)
  Rg, Zg = g.RZ_at_rtheta(eq.rminor_guard, theta)
  if RmajorZvert_bounds is None:
    rpad, zpad = 0.02*np.ptp(Rg), 0.02*np.ptp(Zg)
    RmajorZvert_bounds = (float(Rg.min()-rpad), float(Rg.max()+rpad),
                 float(Zg.min()-zpad), float(Zg.max()+zpad))
  bounds = np.asarray(RmajorZvert_bounds, dtype=float)
  if bounds.shape != (4,) or not np.all(np.isfinite(bounds)):
    raise ValueError("RmajorZvert_bounds must contain four finite values: RMIN RMAX ZMIN ZMAX")
  rmin, rmax, zmin, zmax = bounds
  if not (0 < rmin < g.Rmajor_axis < rmax and zmin < g.Zvert_axis < zmax):
    raise ValueError("R–Z bounds must have RMIN>0 and contain the magnetic axis")
  Rgrid, Zgrid = np.linspace(rmin, rmax, nr), np.linspace(zmin, zmax, nz)
  spacing = np.array([(rmax-rmin)/(nr-1), (zmax-zmin)/(nz-1)])
  guard_tree = cKDTree(np.column_stack((Rg, Zg))/spacing)
  clearance = math.inf
  for radius in (eq.a, eq.rminor_required):
    Rc, Zc = g.RZ_at_rtheta(radius, theta)
    distances, _ = guard_tree.query(np.column_stack((Rc, Zc))/spacing)
    edge = np.minimum.reduce(((Rc-rmin)/spacing[0], (rmax-Rc)/spacing[0],
                              (Zc-zmin)/spacing[1], (zmax-Zc)/spacing[1]))
    clearance = min(clearance, float(distances.min()), float(edge.min()))
  if clearance < 3:
    raise ValueError(f"only {clearance:.3g} cells of contour clearance; need 3: "
                     "increase nr/nz or guard_fraction, or enlarge R–Z bounds")
  RR, ZZ = np.meshgrid(Rgrid, Zgrid)
  radius, _ = eq.rtheta_at_RZ(RR, ZZ)
  psi = eq.psi_at_r(radius)
  if not np.all(np.isfinite(psi)):
    raise ValueError("nonfinite flux on the R–Z mesh")
  flux_nodes = np.linspace(0, eq.psi_boundary, nr)
  qpsi = np.sign(g.F)*np.array([eq.q_at_r(eq.r_at_psi(value)) for value in flux_nodes])
  boundary = np.column_stack(g.RZ_at_rtheta(eq.a, np.linspace(0, 2*np.pi, 257)))
  return {"R": Rgrid, "Z": Zgrid, "psi": psi, "qpsi": qpsi,
          "boundary": boundary, "exterior": radius > eq.rminor_guard,
          "clearance_cells": clearance}


def _write_values(stream, values):
  """Write flattened values as GEQDSK records of up to five E16.9 fields."""
  array = np.asarray(values).ravel()
  for first in range(0, len(array), 5):
    fields = [f"{float(value):16.9E}" for value in array[first:first+5]]
    if any(len(field) != 16 for field in fields):
      raise ValueError("value exceeds GEQDSK E16.9 width")
    stream.write("".join(fields)+"\n")


def _atomic_text(path, writer):
  """Write in the destination directory and replace only after successful close."""
  path = Path(path)
  temporary = None
  try:
    with tempfile.NamedTemporaryFile(mode="w", encoding="ascii", newline="\n",
                                     dir=path.parent, prefix=f".{path.name}.", delete=False) as stream:
      temporary = Path(stream.name)
      writer(stream)
    temporary.replace(path)
  finally:
    if temporary is not None:
      temporary.unlink(missing_ok=True)


def write_eqdsk(filename, equilibrium, sampled):
  """Write a deterministic GEQDSK with explicit zero pressure/current placeholders."""
  eq, g = equilibrium, equilibrium.geometry
  R, Z = sampled["R"], sampled["Z"]
  nr, nz = len(R), len(Z)

  def write(stream):
    """Serialize the GEQDSK header, profiles, flux mesh, and LCFS into the stream."""
    stream.write(f"{'GKEYLL SYNTHETIC MILLER':48s}{0:4d}{nr:4d}{nz:4d}\n")
    _write_values(stream, [R[-1]-R[0], Z[-1]-Z[0], g.Rmajor_axis, R[0], (Z[0]+Z[-1])/2])
    _write_values(stream, [g.Rmajor_axis, g.Zvert_axis, 0, eq.psi_boundary, g.bmag_axis])
    _write_values(stream, [0, 0, 0, g.Rmajor_axis, 0])
    _write_values(stream, [g.Zvert_axis, 0, eq.psi_boundary, 0, 0])
    for profile in (np.full(nr, g.F), np.zeros(nr), np.zeros(nr), np.zeros(nr)):
      _write_values(stream, profile)
    _write_values(stream, sampled["psi"])
    _write_values(stream, sampled["qpsi"])
    stream.write(f"{len(sampled['boundary']):5d}{0:5d}\n")
    _write_values(stream, sampled["boundary"])

  _atomic_text(filename, write)


def summary(equilibrium, sampled, source=None):
  """Build JSON metadata describing inputs, coverage, numerics, and radial conversions."""
  eq = equilibrium
  table = [eq.fs_coords_at_R_omp(eq.geometry.R_omp_at_r(radius)) for radius in np.linspace(0, eq.rminor_guard, 129)]
  return {"generator": "gkeyll_gen_miller_eqdsk", "geometry": asdict(eq.geometry),
          "source": source, "psi_axis": 0.0, "psi_boundary": eq.psi_boundary,
          "guard": {"r": eq.rminor_guard, "R_omp": float(eq.geometry.R_omp_at_r(eq.rminor_guard)),
                    "psi": eq.psi_guard},
          "certified_outer_surface": eq.fs_coords_at_R_omp(eq.geometry.R_omp_at_r(eq.rminor_required)),
          "requested_bounds": None if eq.bounds is None else [eq.fs_coords_at_R_omp(R) for R in eq.bounds],
          "numerics": {"nr": len(sampled["R"]), "nz": len(sampled["Z"]), "rtol": eq.rtol,
                       "guard_fraction": eq.guard_fraction, "num_nodes_radial": len(eq.s_nodes),
                       "radial_relative_error": eq.radial_error,
                       "clearance_cells": sampled["clearance_cells"]},
          "RmajorZvert_bounds": [float(sampled[k][i]) for k, i in (("R", 0), ("R", -1), ("Z", 0), ("Z", -1))],
          "exterior_nodes": int(np.count_nonzero(sampled["exterior"])),
          "conventions": {"psi": "poloidal flux / (2*pi), increasing outwards; axis=0",
                          "q_argument": "outboard midplane major radius [m]",
                          "q": "positive input magnitude; signed GEQDSK q follows bmag_axis",
                          "bmag_axis": "signed toroidal field at the magnetic axis [T]",
                          "exterior": "linear in r^2 beyond guard; not physical coverage",
                          "pressure_and_current": "zero placeholders; not a Grad-Shafranov solution",
                          "limiter": "not specified"},
          "radial_table": table}


def plot_eqdsk(filename, equilibrium, sampled):
  """Save diagnostic panels showing flux contours, q(R_omp), and psi(r)."""
  import matplotlib
  matplotlib.use("Agg")
  import matplotlib.pyplot as plt

  eq = equilibrium
  fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
  contour = axes[0].contour(sampled["R"], sampled["Z"], sampled["psi"], levels=20)
  axes[0].contourf(sampled["R"], sampled["Z"], sampled["exterior"].astype(float),
                   levels=[0.5, 1.5], colors=["lightgray"], alpha=0.5)
  theta = np.linspace(0, 2*np.pi, 1025)
  for radius, label, style in ((eq.a, "LCFS", "-"), (eq.rminor_guard, "guard", "--")):
    axes[0].plot(*eq.geometry.RZ_at_rtheta(radius, theta), style, label=label)
  if eq.bounds:
    for index, R in enumerate(eq.bounds):
      axes[0].plot(*eq.geometry.RZ_at_rtheta(eq.geometry.r_at_R_omp(R), theta), ":",
                   label=f"requested {'lower' if index == 0 else 'upper'}")
  axes[0].set(xlabel="R [m]", ylabel="Z [m]", aspect="equal")
  axes[0].legend(fontsize="small")
  fig.colorbar(contour, ax=axes[0], label=r"$\psi$ [Wb/rad]")
  radius = np.linspace(0, eq.rminor_guard, 513)
  axes[1].plot(eq.geometry.R_omp_at_r(radius), [eq.q_at_r(r) for r in radius])
  axes[1].set(xlabel="R_omp [m]", ylabel="q")
  axes[2].plot(radius, eq.psi_at_r(radius))
  axes[2].set(xlabel="Minor radius, r [m]", ylabel=r"$\psi$ [Wb/rad]")
  fig.savefig(filename, dpi=160)
  plt.close(fig)


def main(argv=None):
  """Parse CLI options, generate output files, and return an exit status."""
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("input", type=Path, help="Python file defining geometry and q(R_omp)")
  parser.add_argument("--output", type=Path,
    help="GEQDSK output path (default: <cwd>/<input-stem>.geqdsk)")
  parser.add_argument("--nr", type=int, default=257,
    help="R grid nodes; odd, 5..999 (default: %(default)s)")
  parser.add_argument("--nz", type=int, default=257,
    help="Z grid nodes; odd, 5..999 (default: %(default)s)")
  parser.add_argument("--rtol", type=float, default=1e-10,
    help="Relative tolerance for flux quadrature/interpolation (default: %(default)s)")
  parser.add_argument("--rz-bounds", dest="RmajorZvert_bounds", type=float, nargs=4, metavar=("RMIN", "RMAX", "ZMIN", "ZMAX"),
    help="Rectangular R-Z bounds in meters (default: padded guard contour)")
  parser.add_argument("--guard-fraction", type=float, default=0.05,
    help="Fractional minor-radius margin beyond required surfaces (default: %(default)s)")
  parser.add_argument("--plot", type=Path,
    help="Save diagnostic plot to this path (default: no plot)")
  args = parser.parse_args(argv)
  output = args.output or Path(args.input.stem+".geqdsk")
  try:
    geometry, q, bounds, source = load_input(args.input)
    eq = MillerEquilibrium(geometry, q, bounds, rtol=args.rtol, guard_fraction=args.guard_fraction)
    sampled = sample_eqdsk(eq, args.nr, args.nz, args.RmajorZvert_bounds)
    metadata = summary(eq, sampled, source)
    write_eqdsk(output, eq, sampled)
    _atomic_text(str(output)+".json", lambda stream: json.dump(metadata, stream, indent=2, allow_nan=False))
    if args.plot:
      plot_eqdsk(args.plot, eq, sampled)
    print(f"Wrote {output}; psi_axis=0, psi_boundary={eq.psi_boundary:.16g}")
    if eq.bounds:
      print("Requested psi bounds:", *(f"{eq.fs_coords_at_R_omp(R)['psi']:.16g}" for R in eq.bounds))
    print(f"Physical guard: R_omp={geometry.R_omp_at_r(eq.rminor_guard):.16g}; "
          f"radial interpolation error={eq.radial_error:.3g}")
    return 0
  except (ValueError, OSError, TypeError, ImportError) as exc:
    print(f"{parser.prog}: {exc}", file=sys.stderr)
    return 1


if __name__ == "__main__":
  sys.exit(main())
