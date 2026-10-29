#!/usr/bin/env python3
"""Generate a TCV GEQDSK and JSON summary using the Miller Python API."""
import json
from pathlib import Path

from gkeyll_gen_miller_eqdsk import (
  MillerGeometry, MillerEquilibrium, sample_eqdsk, summary, write_eqdsk,
)

# TCV geometry: lengths in meters and toroidal magnetic field in tesla.
geometry = MillerGeometry(
  Rmajor_axis=0.8867856264, Zvert_axis=0.1414361745, bmag_axis=-1.4,
  Rmajor_LCFS=1.0870056099999, kappa=1.4, delta=-0.38, a_shift=0.5,
)

# Convert the requested minor-radius endpoints to outboard major radii.
a = geometry.r_at_R_omp(geometry.Rmajor_LCFS)
Rmajor_omp_bounds = (geometry.R_omp_at_r(a-0.04), geometry.R_omp_at_r(a+0.08))

def q(R_omp):
  """Return the TCV safety factor at outboard midplane major radius R_omp [m]."""
  return ((497.3420166252413*R_omp-1408.736172826569)*R_omp
          + 1331.4134861681464)*R_omp-419.00692601227627

# Build the equilibrium and sample its flux on the rectangular R-Z grid.
equilibrium = MillerEquilibrium(geometry, q, Rmajor_omp_bounds)
sampled = sample_eqdsk(equilibrium, nr=257, nz=257)

output = Path("tcv_miller.geqdsk")
write_eqdsk(output, equilibrium, sampled)
metadata = summary(equilibrium, sampled)
Path(str(output)+".json").write_text(json.dumps(metadata, indent=2, allow_nan=False)+"\n")

# Optional diagnostic plot:
from gkeyll_gen_miller_eqdsk import plot_eqdsk
plot_eqdsk("tcv_miller.png", equilibrium, sampled)

print(f"Wrote {output} and {output}.json")
for name, radius in zip(("lower", "upper"), Rmajor_omp_bounds):
  coordinates = equilibrium.fs_coords_at_R_omp(radius)
  print(f"{name}: R_omp={radius:.12g} m, psi={coordinates['psi']:.12g} Wb/rad")
