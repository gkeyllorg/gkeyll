/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_pkpm.h>
#include <gkyl_mom_pkpm_priv.h>
#include <gkyl_util.h>
}

#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

__global__
static void
set_cu_ptrs(struct mom_type_pkpm *mom_pkpm,
  enum gkyl_basis_type b_type, int cdim, int poly_order, bool diag)
{
  // choose kernel tables based on basis-function type
  const gkyl_mom_pkpm_kern_list *mom_pkpm_kernels, *mom_pkpm_diag_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      mom_pkpm_kernels = ser_mom_pkpm_kernels;
      mom_pkpm_diag_kernels = ser_mom_pkpm_diag_kernels;

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      mom_pkpm_kernels = ten_mom_pkpm_kernels;
      mom_pkpm_diag_kernels = ten_mom_pkpm_diag_kernels;
      
      break;

    default:
      assert(false);
      break;    
  }

  if (diag) {
    mom_pkpm->momt.kernel = CK(mom_pkpm_diag_kernels, cdim, poly_order);
    mom_pkpm->momt.num_mom = 8; // rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp
  }
  else {
    mom_pkpm->momt.kernel = CK(mom_pkpm_kernels, cdim, poly_order);
    mom_pkpm->momt.num_mom = 3; // rho, p_par, p_perp
  }
}

struct gkyl_mom_type*
gkyl_mom_pkpm_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  double mass, bool diag)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_pkpm *mom_pkpm = (struct mom_type_pkpm*)
    gkyl_malloc(sizeof(struct mom_type_pkpm));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int poly_order = cbasis->poly_order;

  mom_pkpm->momt.cdim = cdim;
  mom_pkpm->momt.pdim = pdim;
  mom_pkpm->momt.poly_order = poly_order;
  mom_pkpm->momt.num_config = cbasis->num_basis;
  mom_pkpm->momt.num_phase = pbasis->num_basis;

  if (diag)
    mom_pkpm->momt.num_mom = 8; // rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp
  else
    mom_pkpm->momt.num_mom = 4; // rho, p_par, p_perp, M1
  mom_pkpm->mass = mass;

  mom_pkpm->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_pkpm->momt.flags);
  mom_pkpm->momt.ref_count = gkyl_ref_count_init(gkyl_mom_pkpm_free);
  
  // copy struct to device
  struct mom_type_pkpm *mom_pkpm_cu = (struct mom_type_pkpm*)
    gkyl_cu_malloc(sizeof(struct mom_type_pkpm));
  gkyl_cu_memcpy(mom_pkpm_cu, mom_pkpm, sizeof(struct mom_type_pkpm), GKYL_CU_MEMCPY_H2D);

  set_cu_ptrs<<<1,1>>>(mom_pkpm_cu, cbasis->b_type, cdim, poly_order, diag);

  mom_pkpm->momt.on_dev = &mom_pkpm_cu->momt;
  
  return &mom_pkpm->momt;
}
