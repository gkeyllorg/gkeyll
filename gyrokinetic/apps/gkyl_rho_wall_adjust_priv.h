#pragma once

#include <gkyl_tok_geo.h>

#include <stdbool.h>
#include <stdint.h>
#include <string.h>

// Inward rho increments for the wall adjustment, searched in two phases.
//
// The search walks inward from the requested boundary and takes the FIRST
// admissible rho, so the increment sets the PRECISION of the answer, not just
// the speed: with a 1e-3 step the search can overshoot the true limit by up to
// 1e-3 of rho and discard SOL width that would have fit.
//
// But a single fine step is slow when the request starts far outside. A
// widened ASDEX declaration beginning 0.02 of rho beyond the admissible
// boundary needs ~220 fine trials, each of which REBUILDS the geometry, and
// was killed by a 1800 s per-case timeout mid-convergence.
//
// So: approach on the COARSE lattice, then polish on the FINE one. When a
// coarse step first lands admissible, the answer is bracketed within one
// coarse increment; rewinding one coarse step and re-approaching finely
// recovers the full 1e-4 precision for ~10 extra trials instead of ~220. The
// coarse phase is a search accelerator and never the reported answer.
#define GKYL_RHO_WALL_STEP 0.0001
#define GKYL_RHO_WALL_COARSE_STEP 0.001
#define GKYL_RHO_WALL_MAX_STEPS 10000

// Material-boundary adjustment families. Core, IWL, and the unimplemented
// USN type deliberately have no automatic adjustment policy.
static inline int
gkyl_rho_wall_family(enum gkyl_tok_geo_type type)
{
  switch (type) {
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP:
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL:
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO:
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID:
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP:
      return 1;
    case GKYL_GEOMETRY_TOKAMAK_PF_UP_L:
    case GKYL_GEOMETRY_TOKAMAK_PF_UP_R:
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
      return 2;
    default:
      return 0;
  }
}

// isfinite() can be optimized away by this project's -ffast-math builds.
static inline bool
gkyl_rho_wall_finite(double value)
{
  uint64_t bits;
  memcpy(&bits, &value, sizeof bits);
  return (bits & UINT64_C(0x7ff0000000000000)) !=
    UINT64_C(0x7ff0000000000000);
}

// Return trial number step measured from the ORIGINAL requested boundary.
// SOL shrinks downward toward max(1, other_rho); PF shrinks upward toward
// min(1, other_rho). The fixed other boundary and separatrix are never crossed
// or reached. psi_sep must be the same native separatrix used by construction.
// This computes a candidate only; it does not certify wall containment.
// Both output objects must be distinct, non-NULL, and remain untouched on error.
static inline bool
gkyl_rho_wall_next_phased(double requested_rho, double other_rho,
  double psi_axis, double psi_sep, int family, int coarse_steps, int fine_steps,
  double *rho, double *psi)
{
  const int step = coarse_steps*10 + fine_steps;   // in units of the FINE step
  if (!rho || !psi || rho == psi || coarse_steps < 0 || fine_steps < 0 ||
      step < 1 ||
      step > GKYL_RHO_WALL_MAX_STEPS || (family != 1 && family != 2) ||
      !gkyl_rho_wall_finite(requested_rho) ||
      !gkyl_rho_wall_finite(other_rho) || other_rho < 0.0 ||
      !gkyl_rho_wall_finite(psi_axis) ||
      !gkyl_rho_wall_finite(psi_sep) || psi_axis == psi_sep)
    return false;

  // COARSE_STEP is exactly 10*STEP, so the two phases share one lattice and a
  // polished answer can never land between representable fine points.
  const double delta = GKYL_RHO_WALL_COARSE_STEP*coarse_steps
                     + GKYL_RHO_WALL_STEP*fine_steps;
  double next;
  if (family == 1) {
    if (!(requested_rho > 1.0)) return false;
    const double stop = other_rho > 1.0 ? other_rho : 1.0;
    next = requested_rho-delta;
    if (!(next > stop && next < requested_rho)) return false;
  }
  else {
    if (!(requested_rho >= 0.0 && requested_rho < 1.0)) return false;
    const double stop = other_rho < 1.0 ? other_rho : 1.0;
    next = requested_rho+delta;
    if (!(next < stop && next > requested_rho)) return false;
  }
  if (!gkyl_rho_wall_finite(next)) return false;

  const double flux_span = psi_sep-psi_axis;
  if (!gkyl_rho_wall_finite(flux_span) || flux_span == 0.0) return false;
  const double next_psi = psi_axis+next*next*flux_span;
  if (!gkyl_rho_wall_finite(next_psi)) return false;

  *rho = next;
  *psi = next_psi;
  return true;
}

// Single-phase form, kept so existing callers and tests read unchanged: a
// pure fine-step walk is the coarse count set to zero.
static inline bool
gkyl_rho_wall_next(double requested_rho, double other_rho,
  double psi_axis, double psi_sep, int family, int step,
  double *rho, double *psi)
{
  return gkyl_rho_wall_next_phased(requested_rho, other_rho, psi_axis, psi_sep,
    family, 0, step, rho, psi);
}
