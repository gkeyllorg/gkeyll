// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_bc_pkpm_species_reflect.h>
#include <gkyl_pkpm_layout.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct gkyl_bc_pkpm_species_reflect {
  int dir; // Configuration-space direction of the wall normal.
  enum gkyl_edge_loc edge; // Lower or upper wall.
  int cdim; // Configuration-space dimensions.
  int num_basis; // Phase-space basis size (one block).
  int conf_num_basis; // Configuration-space basis size (layout of bvar).
  int num_blocks; // Number of blocks of the distribution array.
  const struct gkyl_basis *basis; // Phase-space basis (device pointer on device).
  struct gkyl_range skin_r, ghost_r; // Phase-space skin and ghost ranges of the wall.
  struct gkyl_range phase_range_ext, conf_range_ext; // Ranges indexing f and bvar.
  bool use_gpu; // Whether the updater runs on the GPU.
  uint32_t flags;
  struct gkyl_bc_pkpm_species_reflect *on_dev; // pointer to itself or device data
};

// Mirror image of the ghost phase-space cell gidx: fills sidx with the skin cell it reflects
// (mirrored in the wall direction, and in the normalized parallel velocity when flip_w) and
// returns flip_w, decided by the parity of the field across the wall from the cell averages of
// the ghost and skin field directions b', b: flip when the normal field is even (b'_n b_n > 0:
// symmetry wall, copied field) or the image is the reversed mirror (b' . R b < 0, which decides
// when b_n = 0); keep w when the normal field is odd (PEC: b' = R b).
GKYL_CU_DH static inline bool
bc_pkpm_species_reflect_skin(
  const struct gkyl_bc_pkpm_species_reflect *up, const int *gidx, const struct gkyl_array *bvar,
  int *sidx
)
{
  int dir = up->dir, cdim = up->cdim, nc = up->conf_num_basis;
  for (int d = 0; d < cdim + 1; ++d) {
    sidx[d] = gidx[d];
  }
  sidx[dir] = up->edge == GKYL_LOWER_EDGE ? 2 * up->skin_r.lower[dir] - 1 - gidx[dir] :
                                            2 * up->skin_r.upper[dir] + 1 - gidx[dir];
  const double *bg =
    (const double *)gkyl_array_cfetch(bvar, gkyl_range_idx(&up->conf_range_ext, gidx));
  const double *bs =
    (const double *)gkyl_array_cfetch(bvar, gkyl_range_idx(&up->conf_range_ext, sidx));
  double sdot = 0.0;
  for (int i = 0; i < 3; ++i) {
    sdot += bg[i * nc] * (i == dir ? -bs[i * nc] : bs[i * nc]);
  }
  bool flip_w = sdot < 0.0 || bg[dir * nc] * bs[dir * nc] > 0.0;
  if (flip_w) {
    sidx[cdim] = up->ghost_r.lower[cdim] + up->ghost_r.upper[cdim] - gidx[cdim];
  }
  return flip_w;
}

// Ghost cell from its skin cell: every block mirrored in the wall direction (and in the
// normalized parallel velocity when flip_w); the wall-normal Cartesian components of the
// first-harmonic vectors (blocks 2-4 and 5-7) negated.
GKYL_CU_DH static inline void
bc_pkpm_species_reflect_cell(
  const struct gkyl_bc_pkpm_species_reflect *up, bool flip_w, const double *inp, double *out
)
{
  int np = up->num_basis;
  for (int k = 0; k < up->num_blocks; ++k) {
    int loc = k * np;
    up->basis->flip_odd_sign(up->dir, &inp[loc], &out[loc]);
    if (flip_w) {
      up->basis->flip_odd_sign(up->cdim, &out[loc], &out[loc]);
    }
    if (k >= PKPM_DIST_M0 && (k - PKPM_DIST_M0) % 3 == up->dir) {
      for (int m = 0; m < np; ++m) {
        out[loc + m] = -out[loc + m];
      }
    }
  }
}
