#pragma once

// Private header for positivity_shift_vlasov updater, not for direct use in user code.

#include <gkyl_positivity_shift_vlasov.h>
#include <gkyl_positivity_shift_vlasov_kernels.h>
#include <gkyl_mom_vlasov_kernels.h>
#include <gkyl_mom_vlasov_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_util.h>
#include <assert.h>

enum gkyl_positivity_shift_type {
  GKYL_POSITIVITY_SHIFT_TYPE_SHIFT_ONLY = 0,
  GKYL_POSITIVITY_SHIFT_TYPE_MRS_LIMITER, // Moe-Rossmanith-Seal limiter.
};

// Function pointer type for sheath reflection kernels.
typedef bool (*m0_pos_check_t)(const double *m0);
typedef bool (*shift_t)(double ffloor, double *distf);
typedef void (*m0_t)(const double *xc, const double *dx,
  const int *idx, const double *fIn, double* GKYL_RESTRICT out);

typedef struct { m0_pos_check_t kernels[3]; } pos_shift_vlasov_kern_list_m0_pos_check; // For use in kernel tables.
typedef struct { shift_t kernels[3]; } pos_shift_vlasov_kern_list_shift; // For use in kernel tables.
typedef struct { m0_t kernels[3]; } pos_shift_vlasov_kern_list_m0; // For use in kernel tables.

// Serendipity  kernels.
GKYL_CU_D
static const pos_shift_vlasov_kern_list_m0_pos_check pos_shift_vlasov_kern_list_m0_pos_check_tensor[] = {
  { positivity_shift_vlasov_conf_pos_check_1x_tensor_p1, NULL, NULL },
  { positivity_shift_vlasov_conf_pos_check_2x_tensor_p1, NULL, NULL },
  { positivity_shift_vlasov_conf_pos_check_3x_tensor_p1, NULL, NULL },
};

GKYL_CU_D
static const pos_shift_vlasov_kern_list_shift pos_shift_vlasov_kern_list_shift_tensor[] = {
  { positivity_shift_vlasov_shift_only_1x1v_tensor_p1, NULL, NULL },
  { positivity_shift_vlasov_shift_only_1x2v_tensor_p1, NULL, NULL },
  { positivity_shift_vlasov_shift_only_1x3v_tensor_p1, NULL, NULL },
  { positivity_shift_vlasov_shift_only_2x2v_tensor_p1, NULL, NULL },
  { positivity_shift_vlasov_shift_only_2x3v_tensor_p1, NULL, NULL },
  { positivity_shift_vlasov_shift_only_3x3v_tensor_p1, NULL, NULL },
};

GKYL_CU_D
static const pos_shift_vlasov_kern_list_shift pos_shift_vlasov_kern_list_MRSlimiter_tensor[] = {
  { positivity_shift_vlasov_MRS_limiter_1x1v_tensor_p1, NULL, NULL },
  { positivity_shift_vlasov_MRS_limiter_1x2v_tensor_p1, NULL, NULL },
  { positivity_shift_vlasov_MRS_limiter_1x3v_tensor_p1, NULL, NULL },
  { positivity_shift_vlasov_MRS_limiter_2x2v_tensor_p1, NULL, NULL },
  { positivity_shift_vlasov_MRS_limiter_2x3v_tensor_p1, NULL, NULL },
  { positivity_shift_vlasov_MRS_limiter_3x3v_tensor_p1, NULL, NULL },
};

GKYL_CU_D
static const pos_shift_vlasov_kern_list_m0 pos_shift_vlasov_kern_list_m0_tensor[] = {
  { vlasov_M0_1x1v_tensor_p1, NULL, NULL },
  { vlasov_M0_1x2v_tensor_p1, NULL, NULL },
  { vlasov_M0_1x3v_tensor_p1, NULL, NULL },
  { vlasov_M0_2x2v_tensor_p1, NULL, NULL },
  { vlasov_M0_2x3v_tensor_p1, NULL, NULL },
  { vlasov_M0_3x3v_tensor_p1, NULL, NULL },
};

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
GKYL_CU_D
static struct { int vdim[4]; } pos_shift_vlasov_cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices
};

struct gkyl_positivity_shift_vlasov_kernels {
  m0_pos_check_t is_m0_positive; // Kernels that checks if m0 is positive.
  shift_t shift;  // Kernel that shifts f to enforce positivity if needed.
  m0_t m0;  // Kernel that computes the number density.
  inv_op_t conf_inv_op; // Conf-space weak inversion (1/A) kernel (p=1 only).
  mul_op_t conf_mul_op; // Conf-space weak multiplication kernel.
  mul_op_t conf_phase_mul_op; // Conf-phase weak multiplication kernel.
};

// Primary struct in this updater.
struct gkyl_positivity_shift_vlasov {
  int num_cbasis;  // Number of conf-space basis monomials.
  struct gkyl_rect_grid grid;  // Phase-space grid.
  double *ffloor;  // Minimum f to shift distribution to when it's <0.
  double ffloor_fac;  // ffloor = max(f)*ffloor_fac.
  double cellav_fac; // Factor multiplying 0th DG coefficient to give cellav.
  bool use_gpu;
  struct gkyl_positivity_shift_vlasov_kernels *kernels;
  struct gkyl_array *shiftedf; // Marks if a shift occured at a given conf-cell.
};

#ifdef GKYL_HAVE_CUDA
// Declaration of cuda device functions.

void
pos_shift_vlasov_choose_shift_kernel_cu(struct gkyl_positivity_shift_vlasov_kernels *kernels,
  struct gkyl_basis cbasis, struct gkyl_basis pbasis, enum gkyl_positivity_shift_type stype);

void
gkyl_positivity_shift_vlasov_advance_cu(gkyl_positivity_shift_vlasov* up,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  struct gkyl_array *GKYL_RESTRICT distf, struct gkyl_array *GKYL_RESTRICT m0,
  struct gkyl_array *GKYL_RESTRICT delta_m0);
#endif

GKYL_CU_D
static void pos_shift_vlasov_choose_shift_kernel(struct gkyl_positivity_shift_vlasov_kernels *kernels,
  struct gkyl_basis cbasis, struct gkyl_basis pbasis, enum gkyl_positivity_shift_type stype,
  bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    pos_shift_vlasov_choose_shift_kernel_cu(kernels, cbasis, pbasis, stype);
    return;
  }
#endif

  enum gkyl_basis_type cbasis_type = cbasis.b_type, pbasis_type = pbasis.b_type;
  int cdim = cbasis.ndim, pdim = pbasis.ndim;
  int vdim = pdim-cdim;
  int poly_order = pbasis.poly_order;

  int plin = pos_shift_vlasov_cv_index[cdim].vdim[vdim];

  switch (pbasis_type) {
    case GKYL_BASIS_MODAL_TENSOR:
      kernels->is_m0_positive = pos_shift_vlasov_kern_list_m0_pos_check_tensor[cdim-1].kernels[poly_order-1];
      kernels->shift = stype == GKYL_POSITIVITY_SHIFT_TYPE_SHIFT_ONLY?
        pos_shift_vlasov_kern_list_shift_tensor[plin].kernels[poly_order-1] :
        pos_shift_vlasov_kern_list_MRSlimiter_tensor[plin].kernels[poly_order-1];
      kernels->m0 = pos_shift_vlasov_kern_list_m0_tensor[plin].kernels[poly_order-1];
      kernels->conf_phase_mul_op = choose_mul_conf_phase_kern(pbasis_type, cdim, vdim, poly_order);
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
