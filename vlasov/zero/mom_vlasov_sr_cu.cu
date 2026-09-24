/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_vlasov_sr.h>
#include <gkyl_mom_vlasov_sr_priv.h>
#include <gkyl_util.h>
}

static int
v_num_mom(int vdim, enum gkyl_distribution_moments mom_type)
{
  int num_mom = 0;
  
  switch (mom_type) {
    case GKYL_F_MOMENT_M0:
    case GKYL_F_MOMENT_M2:
      num_mom = 1;
      break;

    case GKYL_F_MOMENT_M1:
    case GKYL_F_MOMENT_M3:
      num_mom = vdim;
      break;   

    case GKYL_F_MOMENT_NI:
      num_mom = vdim+1;
      break;   
    
    case GKYL_F_MOMENT_TIJ:
      num_mom = 1+vdim+(vdim*(vdim+1))/2;
      break;   

    case GKYL_F_MOMENT_M0ENERGYM3:
      num_mom = vdim+2;
      break;  

    default: // Can't happen.
      fprintf(stderr,"Moment option %d not available.\n",mom_type);
      assert(false);
      break;
  }

  return num_mom;
}

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_mom_vlasov_sr_set_auxfields_cu_kernel(const struct gkyl_mom_type *momt, 
  const struct gkyl_array *gamma)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);
  mom_vm_sr->auxfields.gamma = gamma;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_mom_vlasov_sr_set_auxfields_cu(const struct gkyl_mom_type *momt, struct gkyl_mom_vlasov_sr_auxfields auxin)
{
  gkyl_mom_vlasov_sr_set_auxfields_cu_kernel<<<1,1>>>(momt, auxin.gamma->on_dev);
}


__global__
static void
set_cu_ptrs(struct mom_type_vlasov_sr* mom_vm_sr, enum gkyl_distribution_moments mom_type,
  enum gkyl_basis_type b_type, int vdim, int poly_order, int tblidx)
{
  mom_vm_sr->auxfields.gamma = 0;
  
  // choose kernel tables based on basis-function type
  const gkyl_vlasov_sr_mom_kern_list *m0_kernels, *m1i_kernels, 
    *m2_kernels, *m3i_kernels, *Ni_kernels, *Tij_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      m0_kernels = ser_m0_kernels;
      m1i_kernels = ser_m1i_kernels;
      m2_kernels = ser_m2_kernels;
      m3i_kernels = ser_m3i_kernels;
      Ni_kernels = ser_Ni_kernels;
      Tij_kernels = ser_Tij_kernels;
      break;

    default:
      assert(false);
      break;    
  }
  
  switch (mom_type) {
    case GKYL_F_MOMENT_M0:
      mom_vm_sr->momt.kernel = m0_kernels[tblidx].kernels[poly_order];
      mom_vm_sr->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M1:
      mom_vm_sr->momt.kernel = m1i_kernels[tblidx].kernels[poly_order];
      mom_vm_sr->momt.num_mom = vdim;
      break;

    case GKYL_F_MOMENT_M2:
      mom_vm_sr->momt.kernel = m2_kernels[tblidx].kernels[poly_order];
      mom_vm_sr->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M3:
      mom_vm_sr->momt.kernel = m3i_kernels[tblidx].kernels[poly_order];
      mom_vm_sr->momt.num_mom = vdim;
      break;

    case GKYL_F_MOMENT_NI:
      mom_vm_sr->momt.kernel = Ni_kernels[tblidx].kernels[poly_order];
      mom_vm_sr->momt.num_mom = 1+vdim;
      break;

    case GKYL_F_MOMENT_TIJ:
      mom_vm_sr->momt.kernel = Tij_kernels[tblidx].kernels[poly_order];
      mom_vm_sr->momt.num_mom = 1+vdim+(vdim*(vdim+1))/2;
      break;

    default: // can't happen
      assert(false);
      break;
  }
}

struct gkyl_mom_type*
gkyl_mom_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, 
  enum gkyl_distribution_moments mom_type)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_vlasov_sr *mom_vm_sr = (struct mom_type_vlasov_sr*)
    gkyl_malloc(sizeof(struct mom_type_vlasov_sr));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_vm_sr->momt.cdim = cdim;
  mom_vm_sr->momt.pdim = pdim;
  mom_vm_sr->momt.poly_order = poly_order;
  mom_vm_sr->momt.num_config = cbasis->num_basis;
  mom_vm_sr->momt.num_phase = pbasis->num_basis;

  mom_vm_sr->conf_range = *conf_range;
  mom_vm_sr->vel_range = *vel_range;

  mom_vm_sr->momt.num_mom = v_num_mom(vdim, mom_type); // number of moments

  mom_vm_sr->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_vm_sr->momt.flags);
  mom_vm_sr->momt.ref_count = gkyl_ref_count_init(gkyl_mom_vm_sr_free);
  
  // copy struct to device
  struct mom_type_vlasov_sr *momt_cu = (struct mom_type_vlasov_sr*)
    gkyl_cu_malloc(sizeof(struct mom_type_vlasov_sr));
  gkyl_cu_memcpy(momt_cu, mom_vm_sr, sizeof(struct mom_type_vlasov_sr), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);

  set_cu_ptrs<<<1,1>>>(momt_cu, mom_type, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  mom_vm_sr->momt.on_dev = &momt_cu->momt;
  
  return &mom_vm_sr->momt;
}

__global__
static void
set_int_cu_ptrs(struct mom_type_vlasov_sr* mom_vm_sr, enum gkyl_distribution_moments mom_type,
  enum gkyl_basis_type b_type, int vdim, int poly_order, int tblidx)
{
  mom_vm_sr->auxfields.gamma = 0;

  // choose kernel tables based on basis-function type
  const gkyl_vlasov_sr_mom_kern_list *int_five_moments_kernels;  
  
  // set kernel pointer
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      int_five_moments_kernels = ser_int_five_moments_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  switch (mom_type) {
    case GKYL_F_MOMENT_M0ENERGYM3:
      mom_vm_sr->momt.kernel = int_five_moments_kernels[tblidx].kernels[poly_order];
      mom_vm_sr->momt.num_mom = 2+vdim;
      break;

    default:
      assert(false);
      break;
  }
}

struct gkyl_mom_type *
gkyl_int_mom_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, enum gkyl_distribution_moments mom_type)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_vlasov_sr *mom_vm_sr = (struct mom_type_vlasov_sr*)
    gkyl_malloc(sizeof(struct mom_type_vlasov_sr));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_vm_sr->momt.cdim = cdim;
  mom_vm_sr->momt.pdim = pdim;
  mom_vm_sr->momt.poly_order = poly_order;
  mom_vm_sr->momt.num_config = cbasis->num_basis;
  mom_vm_sr->momt.num_phase = pbasis->num_basis;

  mom_vm_sr->momt.num_mom = v_num_mom(vdim, mom_type); // number of moments

  mom_vm_sr->conf_range = *conf_range;
  mom_vm_sr->vel_range = *vel_range;

  mom_vm_sr->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_vm_sr->momt.flags);
  mom_vm_sr->momt.ref_count = gkyl_ref_count_init(gkyl_mom_vm_sr_free);
  
  // copy struct to device
  struct mom_type_vlasov_sr *momt_cu = (struct mom_type_vlasov_sr*)
    gkyl_cu_malloc(sizeof(struct mom_type_vlasov_sr));
  gkyl_cu_memcpy(momt_cu, mom_vm_sr, sizeof(struct mom_type_vlasov_sr), GKYL_CU_MEMCPY_H2D);

  set_int_cu_ptrs<<<1,1>>>(momt_cu, mom_type, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  mom_vm_sr->momt.on_dev = &momt_cu->momt;
  
  return &mom_vm_sr->momt;
}
