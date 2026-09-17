#pragma once

// Private header for positivity_shift_gyrokinetic updater, not for direct use in user code.

#include <gkyl_positivity_shift_gyrokinetic.h>
#include <gkyl_positivity_shift_gyrokinetic_kernels.h>
#include <gkyl_mom_gyrokinetic_kernels.h>
#include <gkyl_mom_gyrokinetic_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_util.h>
#include <assert.h>

enum gkyl_positivity_shift_type {
  GKYL_POSITIVITY_SHIFT_TYPE_SHIFT_ONLY = 0,
  GKYL_POSITIVITY_SHIFT_TYPE_MRS_LIMITER // Moe-Rossmanith-Seal limiter.
};

// Function pointer type for sheath reflection kernels.
typedef bool (*m0_pos_check_t)(const double *m0);
typedef bool (*shift_t)(double ffloor, double *distf);
typedef void (*m0_t)(
  const double *dxv, const double *vmap, double m_, const double *bmag, const double *f,
  double *GKYL_RESTRICT out
);

typedef struct {
  m0_pos_check_t kernels[3];
} pos_shift_gk_kern_list_m0_pos_check; // For use in kernel tables.
typedef struct {
  shift_t kernels[3];
} pos_shift_gk_kern_list_shift; // For use in kernel tables.
typedef struct {
  m0_t kernels[3];
} pos_shift_gk_kern_list_m0; // For use in kernel tables.

// Serendipity  kernels.
GKYL_CU_D static const pos_shift_gk_kern_list_m0_pos_check
  pos_shift_gk_kern_list_m0_pos_check_ser[] = {
    {positivity_shift_gyrokinetic_conf_pos_check_1x_ser_p1, NULL, NULL},
    {positivity_shift_gyrokinetic_conf_pos_check_2x_ser_p1, NULL, NULL},
    {positivity_shift_gyrokinetic_conf_pos_check_3x_ser_p1, NULL, NULL}
};

GKYL_CU_D static const pos_shift_gk_kern_list_shift pos_shift_gk_kern_list_shift_ser[] = {
  {positivity_shift_gyrokinetic_shift_only_1x1v_ser_p1, NULL, NULL},
  {positivity_shift_gyrokinetic_shift_only_1x2v_ser_p1, NULL, NULL},
  {positivity_shift_gyrokinetic_shift_only_2x2v_ser_p1, NULL, NULL},
  {positivity_shift_gyrokinetic_shift_only_3x2v_ser_p1, NULL, NULL}
};

GKYL_CU_D static const pos_shift_gk_kern_list_shift pos_shift_gk_kern_list_MRSlimiter_ser[] = {
  {positivity_shift_gyrokinetic_MRS_limiter_1x1v_ser_p1, NULL, NULL},
  {positivity_shift_gyrokinetic_MRS_limiter_1x2v_ser_p1, NULL, NULL},
  {positivity_shift_gyrokinetic_MRS_limiter_2x2v_ser_p1, NULL, NULL},
  {positivity_shift_gyrokinetic_MRS_limiter_3x2v_ser_p1, NULL, NULL}
};

GKYL_CU_D static const pos_shift_gk_kern_list_m0 pos_shift_gk_kern_list_m0_ser[] = {
  {gyrokinetic_M0_1x1v_ser_p1, gyrokinetic_M0_1x1v_ser_p2, NULL},
  {gyrokinetic_M0_1x2v_ser_p1, gyrokinetic_M0_1x2v_ser_p2, NULL},
  {gyrokinetic_M0_2x2v_ser_p1, gyrokinetic_M0_2x2v_ser_p2, NULL},
  {gyrokinetic_M0_3x2v_ser_p1, NULL, NULL}
};

struct gkyl_positivity_shift_gyrokinetic_kernels {
  m0_pos_check_t is_m0_positive; // Kernels that checks if m0 is positive.
  shift_t shift; // Kernel that shifts f to enforce positivity if needed.
  m0_t m0; // Kernel that computes the number density.
  inv_op_t conf_inv_op; // Conf-space weak inversion (1/A) kernel (p=1 only).
  mul_op_t conf_mul_op; // Conf-space weak multiplication kernel.
  mul_op_t conf_phase_mul_op; // Conf-phase weak multiplication kernel.
};

// Configuration-space regions to restrict the shift to, expressed as inclusive
// cell-index intervals along the field-line direction (z). This is the
// cell-index counterpart of gkyl_positivity_shift_gyrokinetic_regions, computed
// once from the user's z-coordinates in _new().
struct gkyl_positivity_shift_gyrokinetic_shift_region_idx {
  int z_dir; // Field-line (z) direction, the last conf-space direction.
  int num_regions; // Number of regions; 0 applies the shift everywhere.
  int idx_lo[GKYL_MAX_POSITIVITY_SHIFT_REGIONS]; // Lower z cell index of each region (inclusive).
  int idx_hi[GKYL_MAX_POSITIVITY_SHIFT_REGIONS]; // Upper z cell index of each region (inclusive).
};

// Return true if the conf-space cell with index conf_idx should be shifted,
// i.e. its z cell index falls within one of the restriction regions. With no
// regions specified the shift is applied everywhere.
GKYL_CU_DH static inline bool pos_shift_gk_apply_in_cell(
  const struct gkyl_positivity_shift_gyrokinetic_shift_region_idx *regions, const int *conf_idx
)
{
  if (regions->num_regions == 0) {
    return true;
  }

  int z_idx = conf_idx[regions->z_dir];
  for (int r = 0; r < regions->num_regions; r++) {
    if (regions->idx_lo[r] <= z_idx && z_idx <= regions->idx_hi[r]) {
      return true;
    }
  }
  return false;
}

// Primary struct in this updater.
struct gkyl_positivity_shift_gyrokinetic {
  int num_cbasis; // Number of conf-space basis monomials.
  struct gkyl_rect_grid grid; // Phase-space grid.
  double mass; // Species mass.
  double *ffloor; // Minimum f to shift distribution to when it's <0.
  double ffloor_fac; // ffloor = max(f)*ffloor_fac.
  double cellav_fac; // Factor multiplying 0th DG coefficient to give cellav.
  const struct gk_geometry *gk_geom; // Pointer to geometry object.
  const struct gkyl_velocity_map *vel_map; // Pointer to velocity mapping object.
  struct gkyl_positivity_shift_gyrokinetic_shift_region_idx shift_regions; // Conf-cells to shift.
  bool use_gpu;
  struct gkyl_positivity_shift_gyrokinetic_kernels *kernels;
  struct gkyl_array *shiftedf; // Marks if a shift occured at a given conf-cell.
};

#ifdef GKYL_HAVE_CUDA
// Declaration of cuda device functions.

void pos_shift_gk_choose_shift_kernel_cu(
  struct gkyl_positivity_shift_gyrokinetic_kernels *kernels, struct gkyl_basis cbasis,
  struct gkyl_basis pbasis, enum gkyl_positivity_shift_type stype
);

void gkyl_positivity_shift_gyrokinetic_advance_cu(
  gkyl_positivity_shift_gyrokinetic *up, const struct gkyl_range *conf_rng,
  const struct gkyl_range *phase_rng, struct gkyl_array *GKYL_RESTRICT distf,
  struct gkyl_array *GKYL_RESTRICT m0, struct gkyl_array *GKYL_RESTRICT delta_m0
);

void gkyl_positivity_shift_gyrokinetic_quasineutrality_scale_cu(
  gkyl_positivity_shift_gyrokinetic *up, const struct gkyl_range *conf_rng,
  const struct gkyl_range *phase_rng, const struct gkyl_array *GKYL_RESTRICT delta_m0s,
  const struct gkyl_array *GKYL_RESTRICT delta_m0s_tot,
  const struct gkyl_array *GKYL_RESTRICT delta_m0r, const struct gkyl_array *GKYL_RESTRICT m0s,
  struct gkyl_array *GKYL_RESTRICT fs
);
#endif

GKYL_CU_D static void pos_shift_gk_choose_shift_kernel(
  struct gkyl_positivity_shift_gyrokinetic_kernels *kernels, struct gkyl_basis cbasis,
  struct gkyl_basis pbasis, enum gkyl_positivity_shift_type stype, bool use_gpu
)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    pos_shift_gk_choose_shift_kernel_cu(kernels, cbasis, pbasis, stype);
    return;
  }
#endif

  enum gkyl_basis_type cbasis_type = cbasis.b_type, pbasis_type = pbasis.b_type;
  int cdim = cbasis.ndim, pdim = pbasis.ndim;
  int poly_order = pbasis.poly_order;

  switch (pbasis_type) {
  case GKYL_BASIS_MODAL_GKHYBRID:
  case GKYL_BASIS_MODAL_SERENDIPITY:
    kernels->is_m0_positive =
      pos_shift_gk_kern_list_m0_pos_check_ser[cdim - 1].kernels[poly_order - 1];
    kernels->shift = stype == GKYL_POSITIVITY_SHIFT_TYPE_SHIFT_ONLY ?
                       pos_shift_gk_kern_list_shift_ser[pdim - 2].kernels[poly_order - 1] :
                       pos_shift_gk_kern_list_MRSlimiter_ser[pdim - 2].kernels[poly_order - 1];
    kernels->m0 = pos_shift_gk_kern_list_m0_ser[pdim - 2].kernels[poly_order - 1];
    kernels->conf_phase_mul_op =
      choose_mul_conf_phase_kern(pbasis_type, cdim, pdim - cdim, poly_order);
    break;
  default:
    assert(false);
    break;
  }

  switch (cbasis_type) {
  case GKYL_BASIS_MODAL_SERENDIPITY:
    kernels->conf_inv_op = choose_ser_inv_kern(cdim, poly_order);
    kernels->conf_mul_op = choose_ser_mul_kern(cdim, poly_order);
    break;
  default:
    assert(false);
    break;
  }
}
