// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_gk_neut_hamil_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*hamil_t)(const double *w, const double *dv, 
  const double* gij, double* GKYL_RESTRICT hamil);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
GKYL_CU_D
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { hamil_t kernels[3]; } gkyl_dg_calc_gk_neut_hamil_kern_list;

struct gkyl_dg_calc_gk_neut_hamil {
  struct gkyl_rect_grid phase_grid; // Phase-space grid for cell spacing and cell center 

  hamil_t calc_hamil; // kernel for computing hamiltonian.
 
  uint32_t flags;
  struct gkyl_dg_calc_gk_neut_hamil *on_dev; // pointer to itself or device data.
};

// Calculate hamiltonian
GKYL_CU_D
static const gkyl_dg_calc_gk_neut_hamil_kern_list tensor_gk_neut_hamil_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_neut_hamil_1x3v_tensor_p1, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, gk_neut_hamil_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, gk_neut_hamil_3x3v_tensor_p1, NULL }, // 5
};

GKYL_CU_D
static hamil_t
choose_kern(enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_TENSOR:  
      return tensor_gk_neut_hamil_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
  return 0;
}

#ifdef GKYL_HAVE_CUDA
/**
 * Create new updater to compute relativistic variables on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_gk_neut_hamil* 
gkyl_dg_calc_gk_neut_hamil_cu_dev_new(const struct gkyl_rect_grid *phase_grid,
  const struct gkyl_basis *basis, int cdim);

/**
 * Host-side wrappers for sr vars operations on device
 */

 void gkyl_dg_calc_gk_neut_hamil_calc_cu(struct gkyl_dg_calc_gk_neut_hamil *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_array* gij, struct gkyl_array* hamil);
#endif
