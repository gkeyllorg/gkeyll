"""LTX shape with a positive q profile preserving the old LCFS q and shear."""
import math

geometry = dict(Rmajor_axis=0.406051632, bmag_axis=-0.240606108, Rmajor_LCFS=0.61183,
                kappa=1.3, delta=0.4)
Rmajor_omp_bounds = (0.58, 0.7)


def q(R_omp):
  # Safety factor at outboard midplane major radius R_omp [m].
  # Positive exponential preserving LCFS q=4.3153848 and shear r*q'/q=2.6899871.
  radius = R_omp-geometry["Rmajor_axis"]  # Zero shift in this example.
  a = geometry["Rmajor_LCFS"]-geometry["Rmajor_axis"]
  return 4.3153848*math.exp(0.5*2.6899871*((radius/a)**2-1))
