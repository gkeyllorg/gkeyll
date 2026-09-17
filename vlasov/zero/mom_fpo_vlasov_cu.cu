extern "C" {
  #include <gkyl_alloc.h>
  #include <gkyl_alloc_flags_priv.h>
  #include <gkyl_mom_fpo_vlasov.h>
  #include <gkyl_mom_fpo_vlasov_priv.h>
}

#include <assert.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

// CUDA kernel to set pointer to drag coefficient and diffusion tensor.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_mom_fpo_vlasov_set_auxfields_cu_kernel(const struct gkyl_mom_type *momt, const struct gkyl_array *drag_coeff, const struct gkyl_array *diff_coeff)
{
  struct mom_type_fpo_vlasov *mom_fpo_vlasov = container_of(momt, struct mom_type_fpo_vlasov, momt);
  mom_fpo_vlasov->auxfields.a = drag_coeff;
  mom_fpo_vlasov->auxfields.D = diff_coeff;
}

//// Host-side wrapper for device kernels setting drag coefficient and diffusion tensor.
void
gkyl_mom_fpo_vlasov_set_auxfields_cu(const struct gkyl_mom_type *momt, struct gkyl_mom_fpo_vlasov_auxfields auxin)
{
  gkyl_mom_fpo_vlasov_set_auxfields_cu_kernel<<<1,1>>>(momt, auxin.a->on_dev, auxin.D->on_dev);
}

// CUDA kernel to set device pointers to range object and vlasov fpo kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__
static void
set_cu_ptrs(struct mom_type_fpo_vlasov *mom_fpo_vlasov, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  // choose kernel tables based on basis-function type
  const gkyl_mom_fpo_vlasov_kern_list *mom_fpo_vlasov_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      mom_fpo_vlasov_kernels = ser_mom_fpo_vlasov_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  mom_fpo_vlasov->momt.kernel = CK(mom_fpo_vlasov_kernels, cdim, poly_order);
  mom_fpo_vlasov->momt.num_mom = 4; // 4 component field ((a + div(D)), v . a + div(D . v))
}


struct gkyl_mom_type*
gkyl_mom_fpo_vlasov_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range)
{
  struct mom_type_fpo_vlasov *mom_fpo_vlasov = (struct mom_type_fpo_vlasov*) gkyl_malloc(sizeof(struct mom_type_fpo_vlasov));
  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int poly_order = cbasis->poly_order;

  mom_fpo_vlasov->momt.cdim = cdim;
  mom_fpo_vlasov->momt.pdim = pdim;
  mom_fpo_vlasov->momt.poly_order = poly_order;
  mom_fpo_vlasov->momt.num_config = cbasis->num_basis;
  mom_fpo_vlasov->momt.num_phase = pbasis->num_basis;

  mom_fpo_vlasov->momt.num_mom = 4; // 4 component field ((a + div(D)), v . a + div(D . v))

  mom_fpo_vlasov->phase_range = *phase_range;

  mom_fpo_vlasov->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_fpo_vlasov->momt.flags);
  mom_fpo_vlasov->momt.ref_count = gkyl_ref_count_init(gkyl_mom_fpo_vlasov_free);
 
  // copy struct to device
  struct mom_type_fpo_vlasov *mom_fpo_vlasov_cu = (struct mom_type_fpo_vlasov*)
    gkyl_cu_malloc(sizeof(struct mom_type_fpo_vlasov));
  gkyl_cu_memcpy(mom_fpo_vlasov_cu, mom_fpo_vlasov, sizeof(struct mom_type_fpo_vlasov), GKYL_CU_MEMCPY_H2D);

  // set pointers
  set_cu_ptrs<<<1,1>>>(mom_fpo_vlasov_cu, cbasis->b_type, cdim, poly_order);

  mom_fpo_vlasov->momt.on_dev = &mom_fpo_vlasov_cu->momt;

  return &mom_fpo_vlasov->momt;
}



