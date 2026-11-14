// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_dg_pkpm_harm_vol.h>
#include <gkyl_pkpm_harm_ops.h>
#include <gkyl_pkpm_harm_slots.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_pkpm_layout.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>

// Tensor Legendre modes of degree <= 2 per configuration direction (3^cdim, k0 + 3 k1 + 9 k2).
#define PKPM_HM_N2 27
// (Output, input) block pairs of the source-rate estimate.
#define PKPM_HM_MAXSRC (PKPM_DIST_NUM_BLOCKS * PKPM_DIST_NUM_BLOCKS)
// The nc x nc matrices of the stored operator of a configuration cell (pkpm_harm_mat_*).
#define PKPM_HM_MAT (PKPM_HO_NOP_MAX * PKPM_MAX_CONF_BASIS * PKPM_MAX_CONF_BASIS)

// Generated kernels: the stored operator of a configuration cell from the p=1 fields (with the raw
// w-flux values at the w-face nodes), its matrices, and the volume terms of a phase cell from them.
typedef void (*pkpm_harm_build_t)(
  double qbym, const double *bvar, const double *grad_b, const double *u, const double *grad_u,
  const double *xi, const double *grad_xi, const double *grad_ln_n, const double *theta,
  const double *dtb, const double *em, double *GKYL_RESTRICT ho, double *GKYL_RESTRICT wn
);
typedef void (*pkpm_harm_mat_t)(
  const double *dxv, const double *jacob_pos, const double *ho, double *GKYL_RESTRICT mat
);
typedef void (*pkpm_harm_vol_t)(
  const double *w, const double *dxv, const double *mat, const double *f, double *GKYL_RESTRICT out
);

struct gkyl_dg_pkpm_harm_vol {
  struct gkyl_rect_grid phase_grid; // Phase-space grid.
  int cdim; // Configuration-space dimensions.
  int n2; // Number of degree <= 2 tensor modes (3^cdim).
  int nop; // Number of entries of the stored operator (gkyl_pkpm_harm_ops.h).
  double qbym; // Charge to mass ratio (gyration).
  double bmax[PKPM_HM_N2]; // Largest |value| over the cell of each degree <= 2 mode.
  double wmax; // Larger |w| edge of the velocity domain.
  // the (output, input) source pairs of the rate estimate: the source entry of each power of w
  // (-1 if absent); the gyration is rated apart (pkpm_hm_src_bound)
  int nsrc, src_o[PKPM_HM_MAXSRC], src_i[PKPM_HM_MAXSRC], src_ent[PKPM_HM_MAXSRC][3];
  pkpm_harm_build_t build_kern; // stored operator of a configuration cell
  pkpm_harm_mat_t mat_kern; // matrices of the stored operator of a configuration cell
  pkpm_harm_vol_t vol_kern; // volume terms of a phase cell
  const struct gkyl_vlasov_position_map
    *pos_map; // Position map (acquired; host pointer, 0 on the device copy).
  const struct gkyl_array
    *jacob_pos; // Per-direction Jacobians (borrowed from pos_map; device pointer on device).
  bool use_gpu; // Whether the updater runs on the GPU.
  uint32_t flags;
  struct gkyl_dg_pkpm_harm_vol *on_dev; // pointer to itself or device data
};

// Fill the updater (sizes, the largest values of the degree <= 2 modes, the source pairs of the
// rate) from the input and the generated tables; host function shared by the host and device
// constructors (the kernel pointers are set by each).
void gkyl_dg_pkpm_harm_vol_setup(
  struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_dg_pkpm_harm_vol_inp *inp
);

// Generated kernels (p=1 serendipity; indexed by cdim-1).
GKYL_CU_D static const pkpm_harm_build_t ser_pkpm_harm_build_kernels[] = {
  pkpm_harm_build_1x_ser_p1, pkpm_harm_build_2x_ser_p1, pkpm_harm_build_3x_ser_p1
};
GKYL_CU_D static const pkpm_harm_mat_t ser_pkpm_harm_mat_kernels[] = {
  pkpm_harm_mat_1x_ser_p1, pkpm_harm_mat_2x_ser_p1, pkpm_harm_mat_3x_ser_p1
};
GKYL_CU_D static const pkpm_harm_vol_t ser_pkpm_harm_vol_kernels[] = {
  pkpm_harm_vol_1x1v_ser_p1, pkpm_harm_vol_2x1v_ser_p1, pkpm_harm_vol_3x1v_ser_p1
};

// Bound of the source rate of the column of a configuration cell (time-step rate, plan 17.27, 17.44):
// the largest row sum (Gershgorin) of the source matrix with every entry bounded over the column by
// its absolute coefficients, sum_j W^j sum_k |e_jk| max|mode_k|, W the larger |w| edge (no w of the
// column and no point of the cell exceeds it), plus the bound of |Omega| = |(q/m) B| over the cell,
// the spectral radius of the gyration (explicit in SSP-RK3; its row sums would depend on the
// direction of B relative to the axes).
GKYL_CU_DH static inline double
pkpm_hm_src_bound(const struct gkyl_dg_pkpm_harm_vol *up, const double *ho)
{
  int n2 = up->n2;
  double row[PKPM_DIST_NUM_BLOCKS] = {0.0}, wj[3] = {1.0, up->wmax, up->wmax * up->wmax};
  for (int p = 0; p < up->nsrc; ++p) {
    for (int j = 0; j < 3; ++j) {
      if (up->src_ent[p][j] >= 0) {
        const double *c = &ho[up->src_ent[p][j] * n2];
        double v = 0.0;
        for (int k = 0; k < n2; ++k) {
          v += up->bmax[k] * fabs(c[k]);
        }
        row[up->src_o[p]] += wj[j] * v;
      }
    }
  }
  double rate = 0.0;
  for (int o = 0; o < PKPM_DIST_NUM_BLOCKS; ++o) {
    rate = fmax(rate, row[o]);
  }
  double om2 = 0.0;
  for (int e = 0; e < 3; ++e) {
    const double *c = &ho[PKPM_HO_GYR(e) * n2];
    double v = 0.0;
    for (int k = 0; k < n2; ++k) {
      v += up->bmax[k] * fabs(c[k]);
    }
    om2 += v * v;
  }
  return rate + sqrt(om2);
}
