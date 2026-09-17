#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_bcorr_fpo_vlasov.h>
#include <gkyl_mom_bcorr_fpo_vlasov_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

void
gkyl_mom_bcorr_fpo_vlasov_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_mom_type *momt = container_of(ref, struct gkyl_mom_type, ref_count);
  if (GKYL_IS_CU_ALLOC(momt->flags))
    gkyl_cu_free(momt->on_dev);
  gkyl_free(momt);
}

void
gkyl_mom_bcorr_fpo_vlasov_set_auxfields(const struct gkyl_mom_type *momt, struct gkyl_mom_bcorr_fpo_vlasov_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_mom_type_is_cu_dev(momt)) {
    return gkyl_mom_bcorr_fpo_vlasov_set_auxfields_cu(momt, auxin);
  }
#endif

  struct mom_type_bcorr_fpo_vlasov *mom_bcorr = container_of(momt, struct mom_type_bcorr_fpo_vlasov, momt);
  mom_bcorr->auxfields.D = auxin.D;
}

struct gkyl_mom_type*
gkyl_mom_bcorr_fpo_vlasov_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range, const double* vBoundary, bool use_gpu)
{
  assert(cbasis->poly_order == pbasis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_mom_bcorr_fpo_vlasov_cu_dev_new(cbasis, pbasis, phase_range, vBoundary);
  } 
#endif  
  struct mom_type_bcorr_fpo_vlasov *mom_bcorr = gkyl_malloc(sizeof(struct mom_type_bcorr_fpo_vlasov));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_bcorr->momt.cdim = cdim;
  mom_bcorr->momt.pdim = pdim;
  mom_bcorr->momt.poly_order = poly_order;
  mom_bcorr->momt.num_config = cbasis->num_basis;
  mom_bcorr->momt.num_phase = pbasis->num_basis;
  mom_bcorr->use_gpu = false;

  // FPO is 3V by default
  mom_bcorr->vBoundary[0] = vBoundary[0];
  mom_bcorr->vBoundary[1] = vBoundary[1];
  mom_bcorr->vBoundary[2] = vBoundary[2];
  mom_bcorr->vBoundary[3] = vBoundary[3];
  mom_bcorr->vBoundary[4] = vBoundary[4];
  mom_bcorr->vBoundary[5] = vBoundary[5];

  // choose kernel tables based on basis-function type
  const gkyl_mom_bcorr_fpo_vlasov_kern_list *mom_bcorr_fpo_vlasov_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      mom_bcorr_fpo_vlasov_kernels = ser_mom_bcorr_fpo_vlasov_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  mom_bcorr->momt.kernel = CK(mom_bcorr_fpo_vlasov_kernels, cdim, poly_order);
  mom_bcorr->momt.num_mom = 8; // 8 component field (2 sets of 3 momentum corrections, 2 sets of energy corrections)
                               // One set of corrections is identical to the LBO, the other involves integrations of
                               // the diffusion tensor at the edge of velocity space.

  mom_bcorr->phase_range = *phase_range;

  mom_bcorr->auxfields.D = 0;

  mom_bcorr->momt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mom_bcorr->momt.flags);
  mom_bcorr->momt.ref_count = gkyl_ref_count_init(gkyl_mom_bcorr_fpo_vlasov_free);

  mom_bcorr->momt.on_dev = &mom_bcorr->momt;
 
  return &mom_bcorr->momt;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_mom_type*
gkyl_mom_bcorr_fpo_vlasov_cu_dev_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range, 
  const double* vBoundary)
{
  assert(false);
}

#endif
