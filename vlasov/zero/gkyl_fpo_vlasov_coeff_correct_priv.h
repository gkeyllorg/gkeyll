#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_mom_fpo_vlasov_kernels.h>
#include <gkyl_mat.h>
#include <gkyl_rect_grid.h>
#include <gkyl_range.h>

static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];
  for (int i=cdim; i<pdim; ++i)
    out[i] = vidx[i-cdim];
}

// Kernel function pointers
typedef void (*fpo_correct_mat_set_t)(struct gkyl_mat *lhs, struct gkyl_mat *rhs, 
  const double *fpo_moms, const double *boundary_corrections, const double *moms);

typedef void (*fpo_correct_accum_t)(const double *drag_diff_coeff_corrs, double *drag_coeff,
  double *drag_coeff_surf, double *diff_coeff, double *diff_coeff_surf);

// Object type
struct gkyl_fpo_coeff_correct {
  int num_conf_basis;

  const struct gkyl_rect_grid *grid;
  const struct gkyl_basis *conf_basis;

  bool is_first; // flag to indicate first call to update
  struct gkyl_nmat *As, *xs; // matrices used for LHS and RHS
  gkyl_nmat_mem *mem; // memory for used in batched linear solve

  fpo_correct_mat_set_t mat_set_kernel;
  fpo_correct_accum_t accum_kernel;

  uint32_t flags;

  struct gkyl_fpo_coeff_correct *on_dev; // Pointer to device copy of struct
};

// For use in kernel tables
typedef struct { fpo_correct_mat_set_t kernels[3]; } gkyl_fpo_coeff_correct_mat_set_kern_list;
typedef struct { fpo_correct_accum_t kernels[3]; } gkyl_fpo_coeff_correct_accum_kern_list;

GKYL_CU_D
static const gkyl_fpo_coeff_correct_mat_set_kern_list ser_fpo_coeff_correct_mat_set_kernels[] = {
  {NULL, NULL, NULL },
  { NULL, mom_fpo_vlasov_coeff_correct_mat_1x3v_ser_p1, mom_fpo_vlasov_coeff_correct_mat_1x3v_ser_p2 },
  { NULL, mom_fpo_vlasov_coeff_correct_mat_2x3v_ser_p1, NULL },
  { NULL, NULL, NULL }
};

GKYL_CU_D
static const gkyl_fpo_coeff_correct_accum_kern_list ser_fpo_coeff_correct_accum_kernels[] = 
{
  { NULL, NULL, NULL },
  { NULL, mom_fpo_vlasov_coeff_correct_accum_1x3v_ser_p1, mom_fpo_vlasov_coeff_correct_accum_1x3v_ser_p2 },
  { NULL, mom_fpo_vlasov_coeff_correct_accum_2x3v_ser_p1, NULL },
  { NULL, NULL, NULL }
};
