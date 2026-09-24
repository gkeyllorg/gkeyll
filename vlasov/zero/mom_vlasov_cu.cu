/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_mom_vlasov_priv.h>
#include <gkyl_util.h>
}

static int
v_num_mom(int vdim, enum gkyl_distribution_moments mom_type)
{
  int m3ijk_count[] = { 1, 4, 10 };
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

    case GKYL_F_MOMENT_M2IJ:
      num_mom = vdim*(vdim+1)/2;
      break;

    case GKYL_F_MOMENT_M3IJK:
      num_mom = m3ijk_count[vdim-1];
      break;

    case GKYL_F_MOMENT_M0M1M2:
      num_mom = vdim+2;
      break;      
      
    default: // Can't happen.
      fprintf(stderr,"Moment option %d not available.\n",mom_type);
      assert(false);
      break;
  }

  return num_mom;
}

__global__
static void
set_cu_ptrs(struct mom_type_vlasov* momt, enum gkyl_distribution_moments mom_type,
  enum gkyl_basis_type b_type, int vdim, int poly_order, int tblidx)
{
  int m3ijk_count[] = { 1, 4, 10 };

  momt->momt.kernel = kernel;
  
  // choose kernel tables based on basis-function type
  const gkyl_mom_kern_list *m0_kernels, *m1i_kernels,
    *m2_kernels, *m2ij_kernels, *m3i_kernels, *m3ijk_kernels, *five_moments_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      m0_kernels = ser_m0_kernels;
      m1i_kernels = ser_m1i_kernels;
      m2_kernels = ser_m2_kernels;
      m2ij_kernels = ser_m2ij_kernels;
      m3i_kernels = ser_m3i_kernels;
      m3ijk_kernels = ser_m3ijk_kernels;
      five_moments_kernels = ser_five_moments_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      m0_kernels = tensor_m0_kernels;
      m1i_kernels = tensor_m1i_kernels;
      m2_kernels = tensor_m2_kernels;
      m2ij_kernels = tensor_m2ij_kernels;
      m3i_kernels = tensor_m3i_kernels;
      m3ijk_kernels = tensor_m3ijk_kernels;
      five_moments_kernels = tensor_five_moments_kernels;
      break;

    default:
      assert(false);
      break;    
  }  
  
  switch (mom_type) {
    case GKYL_F_MOMENT_M0:
      momt->kernel = m0_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M1:
      momt->kernel = m1i_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = vdim;
      break;

    case GKYL_F_MOMENT_M2:
      momt->kernel = m2_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M2IJ:
      momt->kernel = m2ij_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = vdim*(vdim+1)/2;
      break;

    case GKYL_F_MOMENT_M3:
      momt->kernel = m3i_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = vdim;
      break;

    case GKYL_F_MOMENT_M3IJK:
      momt->kernel = m3ijk_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = m3ijk_count[vdim-1];
      break;

    case GKYL_F_MOMENT_M0M1M2:
      momt->kernel = five_moments_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = vdim+2;
      break;
      
    default: // can't happen
      break;
  }
}

struct gkyl_mom_type*
gkyl_mom_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, enum gkyl_distribution_moments mom_type)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_vlasov *momt = (struct mom_type_vlasov*)
    gkyl_malloc(sizeof(struct mom_type_vlasov));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  momt->momt.cdim = cdim;
  momt->momt.pdim = pdim;
  momt->momt.poly_order = poly_order;
  momt->momt.num_config = cbasis->num_basis;
  momt->momt.num_phase = pbasis->num_basis;

  momt->momt.num_mom = v_num_mom(vdim, mom_type); // Number of moments.

  momt->momt.flags = 0;
  GKYL_SET_CU_ALLOC(momt->momt.flags);
  momt->momt.ref_count = gkyl_ref_count_init(gkyl_mom_free);
  
  // copy struct to device
  struct mom_type_vlasov *momt_cu = (struct mom_type_vlasov*)
    gkyl_cu_malloc(sizeof(struct mom_type_vlasov));
  gkyl_cu_memcpy(momt_cu, momt, sizeof(struct mom_type_vlasov), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);

  set_cu_ptrs<<<1,1>>>(momt_cu, mom_type, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  momt->momt.on_dev = &momt_cu->momt;
  
  return &momt->momt;
}

__global__
static void
set_int_cu_ptrs(struct mom_type_vlasov* momt, enum gkyl_distribution_moments mom_type,
  enum gkyl_basis_type b_type, int vdim, int poly_order, int tblidx)
{
  momt->momt.kernel = kernel;

  // Choose kernel tables based on basis-function type.
  const gkyl_mom_kern_list *int_five_moments_kernels;

  // Set kernel pointer.
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      int_five_moments_kernels = ser_int_five_moments_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:      
      int_five_moments_kernels = tensor_int_five_moments_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  switch (mom_type) {
    case GKYL_F_MOMENT_M0M1M2:
      momt->kernel = int_five_moments_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = 2+vdim;
      break;

    default:
      assert(false);
      break;
  }
}

struct gkyl_mom_type *
gkyl_int_mom_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, enum gkyl_distribution_moments mom_type)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_vlasov *momt = (struct mom_type_vlasov*)
    gkyl_malloc(sizeof(struct mom_type_vlasov));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  momt->momt.cdim = cdim;
  momt->momt.pdim = pdim;
  momt->momt.poly_order = poly_order;
  momt->momt.num_config = cbasis->num_basis;
  momt->momt.num_phase = pbasis->num_basis;

  momt->momt.num_mom = v_num_mom(vdim, mom_type); // Number of moments.

  momt->momt.flags = 0;
  GKYL_SET_CU_ALLOC(momt->momt.flags);
  momt->momt.ref_count = gkyl_ref_count_init(gkyl_mom_free);
  
  // copy struct to device
  struct mom_type_vlasov *momt_cu = (struct mom_type_vlasov*)
    gkyl_cu_malloc(sizeof(struct mom_type_vlasov));
  gkyl_cu_memcpy(momt_cu, momt, sizeof(struct mom_type_vlasov), GKYL_CU_MEMCPY_H2D);

  set_int_cu_ptrs<<<1,1>>>(momt_cu, mom_type, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  momt->momt.on_dev = &momt_cu->momt;
  
  return &momt->momt;
}
