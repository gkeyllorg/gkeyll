#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_bcorr_lbo_vlasov_kernels.h>
#include <gkyl_vlasov_hyb_build.h>

typedef void (*lbo_vlasov_momf_t)(
  const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv,
  const double *fIn, double *GKYL_RESTRICT out
);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct {
  int vdim[4];
} cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1, 0, 1, 2}, // 1x kernel indices
  {-1, 3, 4, 5}, // 2x kernel indices
  {-1, -1, -1, 6} // 3x kernel indices
};

// for use in kernel tables
typedef struct {
  lbo_vlasov_momf_t kernels[4];
} gkyl_mom_bcorr_lbo_vlasov_kern_list;

//
// Serendipity basis kernels
//

// boundary integral moment correction kernel lists (both momentum and energy)
GKYL_CU_D static const gkyl_mom_bcorr_lbo_vlasov_kern_list ser_mom_bcorr_lbo_vlasov_kernels[] = {
  // 1x kernels
  {NULL, NULL, mom_bcorr_lbo_vlasov_1x1v_ser_p2, mom_bcorr_lbo_vlasov_1x1v_ser_p3}, // 0
  {NULL, NULL, mom_bcorr_lbo_vlasov_1x2v_ser_p2, NULL}, // 1
  {NULL, NULL, mom_bcorr_lbo_vlasov_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, NULL, mom_bcorr_lbo_vlasov_2x1v_ser_p2, mom_bcorr_lbo_vlasov_2x1v_ser_p3}, // 3
  {NULL, NULL, mom_bcorr_lbo_vlasov_2x2v_ser_p2, NULL}, // 4
  {NULL, NULL, mom_bcorr_lbo_vlasov_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, NULL, NULL, NULL} // 6
};

//
// Tensor basis kernels
//

// boundary integral moment correction kernel lists (both momentum and energy)
GKYL_CU_D static const gkyl_mom_bcorr_lbo_vlasov_kern_list ten_mom_bcorr_lbo_vlasov_kernels[] = {
  // 1x kernels
  {NULL, mom_bcorr_lbo_vlasov_1x1v_tensor_p1, mom_bcorr_lbo_vlasov_1x1v_tensor_p2,
   mom_bcorr_lbo_vlasov_1x1v_tensor_p3}, // 0
  {NULL, mom_bcorr_lbo_vlasov_1x2v_tensor_p1, mom_bcorr_lbo_vlasov_1x2v_tensor_p2, NULL}, // 1
  {NULL, mom_bcorr_lbo_vlasov_1x3v_tensor_p1, mom_bcorr_lbo_vlasov_1x3v_tensor_p2, NULL}, // 2
  // 2x kernels
  {NULL, mom_bcorr_lbo_vlasov_2x1v_tensor_p1, mom_bcorr_lbo_vlasov_2x1v_tensor_p2,
   mom_bcorr_lbo_vlasov_2x1v_tensor_p3}, // 3
  {NULL, mom_bcorr_lbo_vlasov_2x2v_tensor_p1, mom_bcorr_lbo_vlasov_2x2v_tensor_p2, NULL}, // 4
  {NULL, GKYL_HYB_2X3V(mom_bcorr_lbo_vlasov_2x3v_tensor_p1), NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(mom_bcorr_lbo_vlasov_3x3v_tensor_p1), NULL, NULL} // 6
};

struct mom_type_bcorr_lbo_vlasov {
  struct gkyl_mom_type momt;
  lbo_vlasov_momf_t kernel; // moment calculation kernel
  double vBoundary[2 * GKYL_MAX_DIM];
};

void gkyl_mom_bcorr_vlasov(const struct gkyl_ref_count *ref);

GKYL_CU_D static void
kernel(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, void *param
)
{
  struct mom_type_bcorr_lbo_vlasov *mom_bcorr =
    container_of(momt, struct mom_type_bcorr_lbo_vlasov, momt);
  enum gkyl_vel_edge edge = *(enum gkyl_vel_edge *)param;

  return mom_bcorr->kernel(idx, edge, mom_bcorr->vBoundary, dx, f, out);
}
