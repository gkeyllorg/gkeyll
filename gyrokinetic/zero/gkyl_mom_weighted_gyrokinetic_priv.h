#pragma once

// Private header for mom_weighted_gyrokinetic updater, not for direct use in user code.

#include <gkyl_mom_weighted_gyrokinetic.h>
#include <gkyl_mom_weighted_gyrokinetic_kernels.h>
#include <gkyl_util.h>
#include <assert.h>

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
GKYL_CU_D static struct {
  int vdim[3];
} cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense
  {-1, 0, 1}, // 1x kernel indices
  {-1, -1, 2}, // 2x kernel indices
  {-1, -1, 3} // 3x kernel indices
};

// Kernel signature.
typedef void (*wmomf_t)(
  const double *dxv, const double *vmap, double m_, double q_, const double *bmag,
  const double *phi, const double *f, const double *wgt, double *GKYL_RESTRICT out
);

// for use in kernel tables
typedef struct {
  wmomf_t kernels[3]; // One for each poly_order.
} gkyl_gyrokinetic_mom_kern_list;

// Generated selection tables and CUDA call targets.
#include <gkyl_mom_weighted_gyrokinetic_registry_priv.h>

struct gkyl_mom_weighted_gyrokinetic_kernels {
  wmomf_t wmom; // Kernels that computes weighted moment.
};

// Object definition.
struct gkyl_mom_weighted_gyrokinetic {
  struct gkyl_rect_grid phase_grid; // Phase-space grid.
  double mass; // Species mass.
  double charge; // Species charge.
  enum gkyl_mom_weight_type wgt_type; // Type of weight.
  const struct gk_geometry *gk_geom; // Pointer to geometry object.
  const struct gkyl_velocity_map *vel_map; // Pointer to velocity mapping object.
  bool use_gpu; // Whether to run on the GPU.
  struct gkyl_mom_weighted_gyrokinetic_kernels *kernels;
  int num_mom; // Number of moments to compute.
};

#ifdef GKYL_HAVE_CUDA
// Declaration of cuda device functions.

void mom_weighted_gk_choose_kernel_cu(
  struct gkyl_mom_weighted_gyrokinetic_kernels *kernels, int *num_mom,
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  enum gkyl_distribution_moments mom_type, enum gkyl_mom_weight_type wgt_type, bool is_integrated
);

void gkyl_mom_weighted_gyrokinetic_advance_cu(
  struct gkyl_mom_weighted_gyrokinetic *up, const struct gkyl_range *phase_rng,
  const struct gkyl_range *conf_rng, const struct gkyl_range *wgt_rng,
  const struct gkyl_array *GKYL_RESTRICT wgt, const struct gkyl_array *GKYL_RESTRICT phi,
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT mout
);
#endif

GKYL_CU_D static void
mom_weighted_gk_choose_kernel(
  struct gkyl_mom_weighted_gyrokinetic_kernels *kernels, int *num_mom,
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  enum gkyl_distribution_moments mom_type, enum gkyl_mom_weight_type wgt_type, bool is_integrated
)
{
  int index = cv_index[cbasis->ndim].vdim[pbasis->ndim - cbasis->ndim];
  int order = pbasis->poly_order;
  const gkyl_gyrokinetic_mom_kern_list *table = NULL;
  if (mom_type == GKYL_F_MOMENT_M0) {
    *num_mom = 1;
    if (wgt_type == GKYL_F_MOMENT_WEIGHT_CONF) {
      table = is_integrated ? ser_int_m0_wgtx_kernels : ser_m0_wgtx_kernels;
    } else {
      table = is_integrated ? ser_int_m0_wgtv_kernels : ser_m0_wgtv_kernels;
    }
  }
  if (mom_type == GKYL_F_MOMENT_M1) {
    *num_mom = 1;
    if (wgt_type == GKYL_F_MOMENT_WEIGHT_CONF) {
      table = is_integrated ? ser_int_m1_wgtx_kernels : ser_m1_wgtx_kernels;
    } else {
      table = is_integrated ? ser_int_m1_wgtv_kernels : ser_m1_wgtv_kernels;
    }
  }
  if (mom_type == GKYL_F_MOMENT_M2) {
    *num_mom = 1;
    if (wgt_type == GKYL_F_MOMENT_WEIGHT_CONF) {
      table = is_integrated ? ser_int_m2_wgtx_kernels : ser_m2_wgtx_kernels;
    } else {
      table = is_integrated ? ser_int_m2_wgtv_kernels : ser_m2_wgtv_kernels;
    }
  }
  if (mom_type == GKYL_F_MOMENT_M2PAR) {
    *num_mom = 1;
    if (wgt_type == GKYL_F_MOMENT_WEIGHT_CONF) {
      table = is_integrated ? ser_int_m2par_wgtx_kernels : ser_m2par_wgtx_kernels;
    } else {
      table = is_integrated ? ser_int_m2par_wgtv_kernels : ser_m2par_wgtv_kernels;
    }
  }
  if (mom_type == GKYL_F_MOMENT_M2PERP) {
    *num_mom = 1;
    if (wgt_type == GKYL_F_MOMENT_WEIGHT_CONF) {
      table = is_integrated ? ser_int_m2perp_wgtx_kernels : ser_m2perp_wgtx_kernels;
    } else {
      table = is_integrated ? ser_int_m2perp_wgtv_kernels : ser_m2perp_wgtv_kernels;
    }
  }
  if (mom_type == GKYL_F_MOMENT_M0M1M2) {
    *num_mom = 3;
    if (wgt_type == GKYL_F_MOMENT_WEIGHT_CONF) {
      table = is_integrated ? ser_int_m0m1m2_wgtx_kernels : ser_m0m1m2_wgtx_kernels;
    } else {
      table = is_integrated ? ser_int_m0m1m2_wgtv_kernels : ser_m0m1m2_wgtv_kernels;
    }
  }
  if (mom_type == GKYL_F_MOMENT_M0M1M2PARM2PERP) {
    *num_mom = 4;
    if (wgt_type == GKYL_F_MOMENT_WEIGHT_CONF) {
      table = is_integrated ? ser_int_m0m1m2parm2perp_wgtx_kernels :
                              ser_m0m1m2parm2perp_wgtx_kernels;
    } else {
      table = is_integrated ? ser_int_m0m1m2parm2perp_wgtv_kernels :
                              ser_m0m1m2parm2perp_wgtv_kernels;
    }
  }
  assert(table && table[index].kernels[order]);
  kernels->wmom = table[index].kernels[order];
}

/**
 * Reject selectors outside density, momentum and second moments.
 *
 * @param moment Requested weighted moment.
 */
static inline void
mom_weighted_gk_validate_weighted(enum gkyl_distribution_moments moment)
{
  switch (moment) {
    case GKYL_F_MOMENT_M0:
    case GKYL_F_MOMENT_M1:
    case GKYL_F_MOMENT_M2:
    case GKYL_F_MOMENT_M2PAR:
    case GKYL_F_MOMENT_M2PERP:
    case GKYL_F_MOMENT_M0M1M2:
    case GKYL_F_MOMENT_M0M1M2PARM2PERP:
      return;
    default:
      gkyl_exit("Unsupported weighted gyrokinetic moment selector.");
  }
}
