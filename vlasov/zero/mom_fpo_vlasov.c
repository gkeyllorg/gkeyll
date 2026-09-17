#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_fpo_vlasov.h>
#include <gkyl_mom_fpo_vlasov_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

void
gkyl_mom_fpo_vlasov_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_mom_type *momt = container_of(ref, struct gkyl_mom_type, ref_count);
  struct mom_type_fpo_vlasov *fpo_momt = container_of(momt, struct mom_type_fpo_vlasov, momt);
  if (GKYL_IS_CU_ALLOC(momt->flags))
    gkyl_cu_free(momt->on_dev);
  gkyl_free(momt);
  gkyl_free(fpo_momt);
}

void
gkyl_mom_fpo_vlasov_set_auxfields(const struct gkyl_mom_type *momt, struct gkyl_mom_fpo_vlasov_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_mom_type_is_cu_dev(momt)) {
    gkyl_mom_fpo_vlasov_set_auxfields_cu(momt->on_dev, auxin);
    return;
  }
#endif

  struct mom_type_fpo_vlasov *mom_fpo_vlasov = container_of(momt, struct mom_type_fpo_vlasov, momt);
  mom_fpo_vlasov->auxfields.a = auxin.a;
  mom_fpo_vlasov->auxfields.D = auxin.D;
}

struct gkyl_mom_type*
gkyl_mom_fpo_vlasov_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range, bool use_gpu)
{
  assert(cbasis->poly_order == pbasis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_mom_fpo_vlasov_cu_dev_new(cbasis, pbasis, phase_range);
  } 
#endif    
  struct mom_type_fpo_vlasov *mom_fpo_vlasov = gkyl_malloc(sizeof(struct mom_type_fpo_vlasov));
  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int poly_order = cbasis->poly_order;

  mom_fpo_vlasov->momt.cdim = cdim;
  mom_fpo_vlasov->momt.pdim = pdim;
  mom_fpo_vlasov->momt.poly_order = poly_order;
  mom_fpo_vlasov->momt.num_config = cbasis->num_basis;
  mom_fpo_vlasov->momt.num_phase = pbasis->num_basis;

  // choose kernel tables based on basis-function type
  const gkyl_mom_fpo_vlasov_kern_list *mom_fpo_vlasov_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      mom_fpo_vlasov_kernels = ser_mom_fpo_vlasov_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  mom_fpo_vlasov->momt.kernel = CK(mom_fpo_vlasov_kernels, cdim, poly_order);
  mom_fpo_vlasov->momt.num_mom = 4; // 4 component field ((a + div(D)), v . a + div(D . v))

  mom_fpo_vlasov->phase_range = *phase_range;

  mom_fpo_vlasov->auxfields.a = 0;
  mom_fpo_vlasov->auxfields.D = 0;
    
  mom_fpo_vlasov->momt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mom_fpo_vlasov->momt.flags);
  mom_fpo_vlasov->momt.ref_count = gkyl_ref_count_init(gkyl_mom_fpo_vlasov_free);
  
  mom_fpo_vlasov->momt.on_dev = &mom_fpo_vlasov->momt; // on host, self-reference
    
  return &mom_fpo_vlasov->momt;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_mom_type*
gkyl_mom_fpo_vlasov_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range)
{
  assert(false);
  return 0;
}

#endif
