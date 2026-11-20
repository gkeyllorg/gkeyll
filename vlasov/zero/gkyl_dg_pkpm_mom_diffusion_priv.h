// Private header: not for direct use
#pragma once

#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_dg_pkpm_mom_diffusion.h>
#include <gkyl_pkpm_vars_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>

// Generated kernel of one direction: increments the momentum blocks of out, returns the rate.
typedef double (*pkpm_md_kern_t)(
  const double *dxv, double D, const double *jac_l, const double *jac_c, const double *jac_r,
  const double *fl, const double *fc, const double *fr, double *GKYL_RESTRICT out
);

struct gkyl_dg_pkpm_mom_diffusion {
  struct gkyl_rect_grid conf_grid; // Configuration-space grid (computational).
  int cdim; // Configuration-space dimensions.
  int order; // Order of the diffusion (2 or 4).
  int num_dir; // Number of directions with a nonzero coefficient.
  int dirs[GKYL_MAX_CDIM]; // Those directions.
  double D[GKYL_MAX_CDIM]; // Diffusion coefficient of each of those directions.
  pkpm_md_kern_t kern[GKYL_MAX_CDIM]; // Kernel of each of those directions.
  const struct gkyl_vlasov_position_map *pos_map; // Position map (acquired; host pointer, 0 on the device copy).
  const struct gkyl_array *jacob_pos; // Per-direction Jacobians (borrowed from pos_map; device pointer on device).
  bool use_gpu; // Whether the updater runs on the GPU.
  uint32_t flags;
  struct gkyl_dg_pkpm_mom_diffusion *on_dev; // pointer to itself or device data
};

// Fill the updater (sizes, directions, coefficients) from the input; host function shared by the
// host and device constructors (the kernel pointers are set by each).
void gkyl_dg_pkpm_mom_diffusion_setup(
  struct gkyl_dg_pkpm_mom_diffusion *up, const struct gkyl_dg_pkpm_mom_diffusion_inp *inp
);

// Generated kernels (p=1 serendipity): [order 2, order 4][cdim - 1][direction].
GKYL_CU_D static const pkpm_md_kern_t ser_pkpm_md_kernels[2][3][3] = {
  {
    {pkpm_mom_diffusion_order2_x_1x_ser_p1, 0, 0},
    {pkpm_mom_diffusion_order2_x_2x_ser_p1, pkpm_mom_diffusion_order2_y_2x_ser_p1, 0},
    {pkpm_mom_diffusion_order2_x_3x_ser_p1, pkpm_mom_diffusion_order2_y_3x_ser_p1, pkpm_mom_diffusion_order2_z_3x_ser_p1},
  },
  {
    {pkpm_mom_diffusion_order4_x_1x_ser_p1, 0, 0},
    {pkpm_mom_diffusion_order4_x_2x_ser_p1, pkpm_mom_diffusion_order4_y_2x_ser_p1, 0},
    {pkpm_mom_diffusion_order4_x_3x_ser_p1, pkpm_mom_diffusion_order4_y_3x_ser_p1, pkpm_mom_diffusion_order4_z_3x_ser_p1},
  },
};

// Kernel of an order, configuration dimension and direction.
GKYL_CU_DH static inline pkpm_md_kern_t
pkpm_md_choose_kern(int order, int cdim, int dir)
{
  return ser_pkpm_md_kernels[order == 2 ? 0 : 1][cdim - 1][dir];
}

// The diffusion of one cell: every direction with a nonzero coefficient, with the neighbours of
// the cell along it; returns the summed rate.
GKYL_CU_DH static inline double
pkpm_md_cell(
  const struct gkyl_dg_pkpm_mom_diffusion *up, const struct gkyl_range *conf_range, const int *idx,
  const struct gkyl_array *fluid, double *GKYL_RESTRICT out
)
{
  int cdim = up->cdim, idx_n[GKYL_MAX_DIM];
  long loc = gkyl_range_idx(conf_range, idx);
  const double *fc = (const double *)gkyl_array_cfetch(fluid, loc);
  const double *jc = (const double *)gkyl_array_cfetch(up->jacob_pos, loc);
  double rate = 0.0;
  for (int k = 0; k < up->num_dir; ++k) {
    int d = up->dirs[k];
    gkyl_copy_int_arr(cdim, idx, idx_n);
    idx_n[d] = idx[d] - 1;
    long loc_l = gkyl_range_idx(conf_range, idx_n);
    idx_n[d] = idx[d] + 1;
    long loc_r = gkyl_range_idx(conf_range, idx_n);
    rate += up->kern[k](up->conf_grid.dx, up->D[k], (const double *)gkyl_array_cfetch(up->jacob_pos, loc_l), jc,
      (const double *)gkyl_array_cfetch(up->jacob_pos, loc_r), (const double *)gkyl_array_cfetch(fluid, loc_l), fc,
      (const double *)gkyl_array_cfetch(fluid, loc_r), out);
  }
  return rate;
}
