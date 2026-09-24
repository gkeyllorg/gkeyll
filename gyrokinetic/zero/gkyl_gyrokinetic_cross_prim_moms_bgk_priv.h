// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gyrokinetic_cross_prim_moms_bgk_kernels.h>
#include <gkyl_mat.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

// Function pointer type for cross moments calculation
typedef void (*gyrokinetic_cross_prim_moms_bgk_t)(const double delta_sr, const double betap1, 
  const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, 
  double *prim_moms_cross);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below.
GKYL_CU_D
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense
  {-1,  0,  1}, // 1x kernel indices
  {-1, -1,  2}, // 2x kernel indices
  {-1, -1,  3}, // 3x kernel indices
};

// For use in kernel tables
typedef struct { gyrokinetic_cross_prim_moms_bgk_t kernels[3]; } gkyl_gyrokinetic_cross_prim_moms_bgk_kern_list;

// Cross moments kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_cross_prim_moms_bgk_kern_list ser_gyrokinetic_cross_prim_moms_bgk_list[] = {
  // 1x kernels 
  { NULL, gyrokinetic_cross_prim_moms_bgk_1x1v_ser_p1, NULL }, //0
  { NULL, gyrokinetic_cross_prim_moms_bgk_1x2v_ser_p1, NULL }, //1
  // 2x kernels
  { NULL, gyrokinetic_cross_prim_moms_bgk_2x2v_ser_p1, NULL }, // no gyrokinetic_cross_prim_moms_bgk_2x2v_ser_p2 due to the lack of gkyl_basis_ser_2x_p2_inv.h //2
  // 3x kernels
  { NULL, gyrokinetic_cross_prim_moms_bgk_3x2v_ser_p1, NULL }, // no gyrokinetic_cross_prim_moms_bgk_2x2v_ser_p2 due to the lack of gkyl_basis_ser_3x_p2_inv.h //3
};

GKYL_CU_D
static gyrokinetic_cross_prim_moms_bgk_t
choose_gyrokinetic_cross_prim_moms_bgk_kern(int cdim, int vdim, int poly_order){
  return ser_gyrokinetic_cross_prim_moms_bgk_list[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}

// Primary struct in this updater
struct gkyl_gyrokinetic_cross_prim_moms_bgk {
  bool use_gpu;
  gyrokinetic_cross_prim_moms_bgk_t cross_prim_moms_calc; // a pointer to the cross primitive moments kernel

  struct gkyl_gyrokinetic_cross_prim_moms_bgk *on_dev; 
};

#ifdef GKYL_HAVE_CUDA
/**
 * Create new updater to compute cross BGK moments on
 * NV-GPU. See new() method for documentation.
 */
gkyl_gyrokinetic_cross_prim_moms_bgk* gkyl_gyrokinetic_cross_prim_moms_bgk_cu_dev_new(
  const struct gkyl_basis *phase_basis, const struct gkyl_basis *conf_basis);

/**
 * Host-side wrappers for cross BGK moments operations on device
 */
void gkyl_gyrokinetic_cross_prim_moms_bgk_advance_cu(gkyl_gyrokinetic_cross_prim_moms_bgk *up,
  const struct gkyl_range *conf_rng, double delta_sr, double betap1,
  double m_self, const struct gkyl_array *prim_moms_self, double m_other, const struct gkyl_array *prim_moms_other,
  struct gkyl_array *prim_moms_cross);
#endif
