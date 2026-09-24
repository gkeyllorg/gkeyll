/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_canonical_pb.h>
#include <gkyl_mom_canonical_pb_priv.h>
#include <gkyl_util.h>
}

#include <cassert>

static int
v_num_mom(int vdim, enum gkyl_distribution_moments mom_type)
{
  int num_mom = 0;
  
  switch (mom_type) {
    case GKYL_F_MOMENT_ENERGY:
      num_mom = 1;
      break;    

    case GKYL_F_MOMENT_M1_FROM_H:
      num_mom = vdim;
      break;   

    case GKYL_F_MOMENT_M0M1M2:
      num_mom = 2+vdim;
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
gkyl_mom_canonical_pb_set_auxfields_cu_kernel(const struct gkyl_mom_type *momt, 
  const struct gkyl_array *hamil)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  mom_can_pb->auxfields.hamil = hamil;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_mom_canonical_pb_set_auxfields_cu(const struct gkyl_mom_type *momt, struct gkyl_mom_canonical_pb_auxfields auxin)
{
  gkyl_mom_canonical_pb_set_auxfields_cu_kernel<<<1,1>>>(momt, auxin.hamil->on_dev);
}


__global__
static void
set_cu_ptrs(struct mom_type_canonical_pb* mom_can_pb, enum gkyl_distribution_moments mom_type,
  enum gkyl_basis_type b_type, int vdim, int poly_order, int tblidx)
{
  mom_can_pb->auxfields.hamil = 0;
  
  // choose kernel tables based on basis-function type
  const gkyl_canonical_pb_mom_kern_list *menergy_kernels;
  const gkyl_canonical_pb_mom_kern_list *m1i_from_h_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      // Verify that the poly-order is 2 for ser case
      assert(poly_order == 2);
      menergy_kernels = ser_menergy_kernels;
      m1i_from_h_kernels = ser_m1i_from_h_kernels;
      break;

    case GKYL_BASIS_MODAL_HYBRID:
      // Verify that the poly-order is 1 for hybrid case
      assert(poly_order == 1);
      menergy_kernels = ser_menergy_kernels;
      m1i_from_h_kernels = ser_m1i_from_h_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      menergy_kernels = tensor_menergy_kernels;
      m1i_from_h_kernels = tensor_m1i_from_h_kernels;
      break;

    default:
      assert(false);
      break;    
  }
  
  switch (mom_type) {
    case GKYL_F_MOMENT_ENERGY:
      mom_can_pb->momt.kernel = menergy_kernels[tblidx].kernels[poly_order];
      mom_can_pb->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M1_FROM_H:
      mom_can_pb->momt.kernel = m1i_from_h_kernels[tblidx].kernels[poly_order];
      mom_can_pb->momt.num_mom = vdim;
      break;

    default: // can't happen
      break;
  }
}

struct gkyl_mom_type*
gkyl_mom_canonical_pb_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range, enum gkyl_distribution_moments mom_type)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_canonical_pb *mom_can_pb = (struct mom_type_canonical_pb*)
    gkyl_malloc(sizeof(struct mom_type_canonical_pb));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_can_pb->momt.cdim = cdim;
  mom_can_pb->momt.pdim = pdim;
  mom_can_pb->momt.poly_order = poly_order;
  mom_can_pb->momt.num_config = cbasis->num_basis;
  mom_can_pb->momt.num_phase = pbasis->num_basis;

  mom_can_pb->phase_range = *phase_range;

  mom_can_pb->momt.num_mom = v_num_mom(vdim, mom_type); // Number of moments.

  mom_can_pb->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_can_pb->momt.flags);
  mom_can_pb->momt.ref_count = gkyl_ref_count_init(gkyl_mom_can_pb_free);
  
  // copy struct to device
  struct mom_type_canonical_pb *momt_cu = (struct mom_type_canonical_pb*)
    gkyl_cu_malloc(sizeof(struct mom_type_canonical_pb));
  gkyl_cu_memcpy(momt_cu, mom_can_pb, sizeof(struct mom_type_canonical_pb), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);

  set_cu_ptrs<<<1,1>>>(momt_cu, mom_type, pbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  mom_can_pb->momt.on_dev = &momt_cu->momt;
  
  return &mom_can_pb->momt;
}

__global__
static void
set_int_cu_ptrs(struct mom_type_canonical_pb* mom_can_pb, enum gkyl_distribution_moments mom_type,
  enum gkyl_basis_type b_type, int vdim, int poly_order, int tblidx)
{
  mom_can_pb->auxfields.hamil = 0;

  // Coose kernel tables based on basis-function type.
  const gkyl_canonical_pb_mom_kern_list *int_five_moments_kernels;
  
  // Set kernel pointer.
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      // Verify that the poly-order is 2 for ser case
      assert(poly_order == 2);
      int_five_moments_kernels = ser_int_five_moments_kernels;
      break;

    case GKYL_BASIS_MODAL_HYBRID:
      // Verify that the poly-order is 1 for hybrid case
      assert(poly_order == 1);
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
      mom_can_pb->momt.kernel = int_five_moments_kernels[tblidx].kernels[poly_order];
      mom_can_pb->momt.num_mom = 2+vdim;
      break;

    default:
      assert(false);
      break;
  }
}

struct gkyl_mom_type *
gkyl_int_mom_canonical_pb_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range, enum gkyl_distribution_moments mom_type)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_canonical_pb *mom_can_pb = (struct mom_type_canonical_pb*)
    gkyl_malloc(sizeof(struct mom_type_canonical_pb));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_can_pb->momt.cdim = cdim;
  mom_can_pb->momt.pdim = pdim;
  mom_can_pb->momt.poly_order = poly_order;
  mom_can_pb->momt.num_config = cbasis->num_basis;
  mom_can_pb->momt.num_phase = pbasis->num_basis;

  mom_can_pb->momt.num_mom = v_num_mom(vdim, mom_type); // Number of moments.

  mom_can_pb->phase_range = *phase_range;

  mom_can_pb->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_can_pb->momt.flags);
  mom_can_pb->momt.ref_count = gkyl_ref_count_init(gkyl_mom_can_pb_free);
  
  // copy struct to device
  struct mom_type_canonical_pb *momt_cu = (struct mom_type_canonical_pb*)
    gkyl_cu_malloc(sizeof(struct mom_type_canonical_pb));
  gkyl_cu_memcpy(momt_cu, mom_can_pb, sizeof(struct mom_type_canonical_pb), GKYL_CU_MEMCPY_H2D);

  set_int_cu_ptrs<<<1,1>>>(momt_cu, mom_type, pbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  mom_can_pb->momt.on_dev = &momt_cu->momt;
  
  return &mom_can_pb->momt;
}

