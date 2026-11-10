#pragma once

#include <math.h>
#include <stdbool.h>
#include <gkyl_util.h>

// Excision is a lapse decision, not a repair for invalid live geometry.
// Retain the all-zero representation of an already-cleared cell at threshold 0.
GKYL_CU_DH static inline bool
gkyl_gr_einstein_is_excised(double threshold, const double *q, int ncomp)
{
  if (isfinite(q[9]) && q[9] < threshold) return true;
  if (q[9] != 0.0) return false;
  for (int m = 0; m < ncomp; ++m)
    if (q[m] != 0.0) return false;
  return true;
}

GKYL_CU_DH static inline bool
gkyl_gr_conformal_scalars_valid(const double *q)
{
  return isfinite(q[9]) && q[9] > 0.0 && isfinite(q[64]) && q[64] > 0.0;
}
