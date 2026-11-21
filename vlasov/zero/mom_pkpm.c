#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_pkpm.h>
#include <gkyl_mom_pkpm_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst, cdim, poly_order) lst[cdim - 1].kernels[poly_order]

void
gkyl_mom_pkpm_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_mom_type *momt = container_of(ref, struct gkyl_mom_type, ref_count);
  if (GKYL_IS_CU_ALLOC(momt->flags)) {
    gkyl_cu_free(momt->on_dev);
  }
  gkyl_free(momt);
}

struct gkyl_mom_type *
gkyl_mom_pkpm_inew(const struct gkyl_mom_pkpm_inp *inp)
{
  assert(inp->conf_basis->poly_order == inp->phase_basis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_mom_pkpm_cu_dev_inew(inp);
  }
#endif
  struct mom_type_pkpm *mom_pkpm = gkyl_malloc(sizeof(struct mom_type_pkpm));
  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim;
  int poly_order = inp->conf_basis->poly_order;

  mom_pkpm->momt.cdim = cdim;
  mom_pkpm->momt.pdim = pdim;
  mom_pkpm->momt.poly_order = poly_order;
  mom_pkpm->momt.num_config = inp->conf_basis->num_basis;
  mom_pkpm->momt.num_phase = inp->phase_basis->num_basis;

  // choose kernel tables based on basis-function type
  const gkyl_mom_pkpm_kern_list *mom_pkpm_kernels, *mom_pkpm_diag_kernels;

  switch (inp->conf_basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
    case GKYL_BASIS_MODAL_TENSOR: // the same basis at polynomial order 1
      mom_pkpm_kernels = ser_mom_pkpm_kernels;
      mom_pkpm_diag_kernels = ser_mom_pkpm_diag_kernels;

      break;

    default:
      assert(false);
      break;
  }

  if (inp->is_diag) {
    mom_pkpm->momt.kernel = CK(mom_pkpm_diag_kernels, cdim, poly_order);
    // mass times [M0, M1, M2, M0_G, M3, M1_G, M4, M2_G] of the input distributions
    mom_pkpm->momt.num_mom = 8;
  } else {
    mom_pkpm->momt.kernel = CK(mom_pkpm_kernels, cdim, poly_order);
    // M0, M2, M0_G, M1, M3, M1_G, int (1, w, w^2) M0_perp dw (3 each), int N_perp dw (3)
    mom_pkpm->momt.num_mom = 18;
  }

  mom_pkpm->mass = inp->mass;

  mom_pkpm->momt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mom_pkpm->momt.flags);
  mom_pkpm->momt.ref_count = gkyl_ref_count_init(gkyl_mom_pkpm_free);

  mom_pkpm->momt.on_dev = &mom_pkpm->momt; // on host, self-reference

  return &mom_pkpm->momt;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_mom_type *
gkyl_mom_pkpm_cu_dev_inew(const struct gkyl_mom_pkpm_inp *inp)
{
  assert(false);
  return 0;
}

#endif
