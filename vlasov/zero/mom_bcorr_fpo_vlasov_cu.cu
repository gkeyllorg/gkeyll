extern "C" {
  #include <gkyl_alloc.h>
  #include <gkyl_alloc_flags_priv.h>
  #include <gkyl_mom_bcorr_fpo_vlasov.h>
  #include <gkyl_mom_bcorr_fpo_vlasov_priv.h>
}

#include <assert.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

// CUDA kernel to set pointer to diffusion tensor.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ 
static void
gkyl_mom_bcorr_fpo_vlasov_set_auxfields_cu_kernel(const struct gkyl_mom_type *momt, const struct gkyl_array *diff_coeff)
{
  struct mom_type_bcorr_fpo_vlasov *mom_fpo_vlasov = container_of(momt, struct mom_type_bcorr_fpo_vlasov, momt);
  mom_fpo_vlasov->auxfields.D = diff_coeff;
}

//// Host-side wrapper for device kernels setting drag coefficient and diffusion tensor.
void 
gkyl_mom_bcorr_fpo_vlasov_set_auxfields_cu(const struct gkyl_mom_type *momt, struct gkyl_mom_bcorr_fpo_vlasov_auxfields auxin)
{
  gkyl_mom_bcorr_fpo_vlasov_set_auxfields_cu_kernel<<<1,1>>>(momt->on_dev, auxin.D->on_dev);
}

// CUDA kernel to set device pointers to range object and vlasov fpo kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__
static void
set_cu_ptrs(struct mom_type_bcorr_fpo_vlasov *mom_bcorr, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  // choose kernel tables based on basis-function type
  const gkyl_mom_bcorr_fpo_vlasov_kern_list *mom_bcorr_fpo_vlasov_kernels;

  switch (b_type) {
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
}


struct gkyl_mom_type*
gkyl_mom_bcorr_fpo_vlasov_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range, const double* vBoundary)
{
  struct mom_type_bcorr_fpo_vlasov *mom_bcorr = (struct mom_type_bcorr_fpo_vlasov *)
    gkyl_malloc(sizeof(struct mom_type_bcorr_fpo_vlasov));
  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int poly_order = cbasis->poly_order;

  mom_bcorr->momt.cdim = cdim;
  mom_bcorr->momt.pdim = pdim;
  mom_bcorr->momt.poly_order = poly_order;
  mom_bcorr->momt.num_config = cbasis->num_basis;
  mom_bcorr->momt.num_phase = pbasis->num_basis;
  mom_bcorr->use_gpu = true;

  mom_bcorr->phase_range = *phase_range;

  // FPO is 3V by default
  mom_bcorr->vBoundary[0] = vBoundary[0];
  mom_bcorr->vBoundary[1] = vBoundary[1];
  mom_bcorr->vBoundary[2] = vBoundary[2];
  mom_bcorr->vBoundary[3] = vBoundary[3];
  mom_bcorr->vBoundary[4] = vBoundary[4];
  mom_bcorr->vBoundary[5] = vBoundary[5];

  mom_bcorr->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_bcorr->momt.flags);
  mom_bcorr->momt.ref_count = gkyl_ref_count_init(gkyl_mom_bcorr_fpo_vlasov_free);

  // copy struct to device
  struct mom_type_bcorr_fpo_vlasov *mom_bcorr_cu = (struct mom_type_bcorr_fpo_vlasov *)
    gkyl_cu_malloc(sizeof(struct mom_type_bcorr_fpo_vlasov));
  gkyl_cu_memcpy(mom_bcorr_cu, mom_bcorr, sizeof(struct mom_type_bcorr_fpo_vlasov), GKYL_CU_MEMCPY_H2D);

  set_cu_ptrs<<<1,1>>>(mom_bcorr_cu, cbasis->b_type, cdim, poly_order);

  mom_bcorr->momt.on_dev = &mom_bcorr_cu->momt;
    
  return &mom_bcorr->momt;
}



