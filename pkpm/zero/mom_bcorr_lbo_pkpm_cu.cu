/* -*- c++ -*- */

#include <assert.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_util.h>
#include <gkyl_mom_bcorr_lbo_pkpm.h>
#include <gkyl_mom_bcorr_lbo_pkpm_priv.h>
}

#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

__global__
static void
gkyl_mom_bcorr_lbo_pkpm_set_cu_dev_ptrs(struct mom_type_bcorr_lbo_pkpm* mom_bcorr, 
  enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  // choose kernel tables based on basis-function type
  const gkyl_mom_bcorr_lbo_pkpm_kern_list *mom_bcorr_lbo_pkpm_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      mom_bcorr_lbo_pkpm_kernels = ser_mom_bcorr_lbo_pkpm_kernels;

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      mom_bcorr_lbo_pkpm_kernels = ten_mom_bcorr_lbo_pkpm_kernels;
      
      break;

    default:
      assert(false);
      break;
  }
  mom_bcorr->momt.kernel = CK(mom_bcorr_lbo_pkpm_kernels, cdim, poly_order);
}

struct gkyl_mom_type*
gkyl_mom_bcorr_lbo_pkpm_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const double *vBoundary, double mass)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_bcorr_lbo_pkpm *mom_bcorr = (struct mom_type_bcorr_lbo_pkpm*) gkyl_malloc(sizeof(struct mom_type_bcorr_lbo_pkpm));

  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int poly_order = cbasis->poly_order;

  mom_bcorr->momt.cdim = cdim;
  mom_bcorr->momt.pdim = pdim;
  mom_bcorr->momt.poly_order = poly_order;
  mom_bcorr->momt.num_config = cbasis->num_basis;
  mom_bcorr->momt.num_phase = pbasis->num_basis;

  mom_bcorr->vBoundary[0] = vBoundary[0];
  mom_bcorr->vBoundary[1] = vBoundary[1];
  mom_bcorr->mass = mass;
  mom_bcorr->momt.num_mom = 2; // number of moments

  mom_bcorr->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_bcorr->momt.flags);
  mom_bcorr->momt.ref_count = gkyl_ref_count_init(gkyl_mom_bcorr_lbo_pkpm_free);

  // copy struct to device
  struct mom_type_bcorr_lbo_pkpm *mom_bcorr_cu = (struct mom_type_bcorr_lbo_pkpm*)
    gkyl_cu_malloc(sizeof(struct mom_type_bcorr_lbo_pkpm));
  gkyl_cu_memcpy(mom_bcorr_cu, mom_bcorr, sizeof(struct mom_type_bcorr_lbo_pkpm), GKYL_CU_MEMCPY_H2D);


  gkyl_mom_bcorr_lbo_pkpm_set_cu_dev_ptrs<<<1,1>>>(mom_bcorr_cu, cbasis->b_type, cdim, poly_order);

  mom_bcorr->momt.on_dev = &mom_bcorr_cu->momt;

  return &mom_bcorr->momt;
}
