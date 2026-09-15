// Private header: not for direct use
#pragma once

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>

#include <gkyl_basis.h>
#include <gkyl_dg_gr_maxwell_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_hyb_build.h>

typedef void (*dg_gr_maxwell_lorentz_conf_t)(const gkyl_dg_gr_maxwell_inp *meq, const double *dx,
  const double *jacob_pos,
  const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal,
  const double *h_ij_inv_nodal, const double *J_c, const double *fields,
  const double *vierb_cov_nodal, const double *vierb_con_nodal,
  double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force);

// The cv_index[cd].vdim[vd] is used to index the kernel tables below.
GKYL_CU_TABLE
static const struct { int vdim[4]; } lorentz_cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense.
  {-1,  0,  1,  2}, // 1x kernel indices.
  {-1, -1,  3,  4}, // 2x kernel indices.
  {-1, -1, -1,  5}, // 3x kernel indices.
};

typedef struct { dg_gr_maxwell_lorentz_conf_t kernels[4]; } gkyl_dg_gr_maxwell_lorentz_conf_kern_list;

struct gkyl_dg_gr_maxwell_lorentz_conf {
  struct gkyl_rect_grid conf_grid; // Configuration-space grid.
  gkyl_dg_gr_maxwell_inp gr_maxwell_data; // Parameters needed by GR-Maxwell kernels.
  int cdim; // Configuration-space dimensions.
  int vdim; // Velocity-space dimensions.
  int num_basis; // Number of configuration-space basis functions.
  double qbym; // Charge over mass.
  const struct gkyl_vlasov_position_map *pos_map; // acquired; jacob_pos borrowed
  const struct gkyl_array *jacob_pos;
  dg_gr_maxwell_lorentz_conf_t lorentz_conf; // Kernel computing local Lorentz-force fields.

  uint32_t flags;
  bool use_gpu;
  struct gkyl_dg_gr_maxwell_lorentz_conf *on_dev; // Pointer to itself or device data.
};

GKYL_CU_D
static const gkyl_dg_gr_maxwell_lorentz_conf_kern_list ser_lorentz_conf_kernels[] = {
  { NULL, dg_gr_maxwell_lorentz_conf_1x1v_ser_p1, dg_gr_maxwell_lorentz_conf_1x1v_ser_p2, NULL }, // 0
  { NULL, dg_gr_maxwell_lorentz_conf_1x2v_ser_p1, dg_gr_maxwell_lorentz_conf_1x2v_ser_p2, NULL }, // 1
  { NULL, dg_gr_maxwell_lorentz_conf_1x3v_ser_p1, dg_gr_maxwell_lorentz_conf_1x3v_ser_p2, NULL }, // 2
  { NULL, dg_gr_maxwell_lorentz_conf_2x2v_ser_p1, dg_gr_maxwell_lorentz_conf_2x2v_ser_p2, NULL }, // 3
  { NULL, dg_gr_maxwell_lorentz_conf_2x3v_ser_p1, dg_gr_maxwell_lorentz_conf_2x3v_ser_p2, NULL }, // 4
  { NULL, dg_gr_maxwell_lorentz_conf_3x3v_ser_p1, dg_gr_maxwell_lorentz_conf_3x3v_ser_p2, NULL }, // 5
};

GKYL_CU_D
static const gkyl_dg_gr_maxwell_lorentz_conf_kern_list tensor_lorentz_conf_kernels[] = {
  { NULL, dg_gr_maxwell_lorentz_conf_1x1v_tensor_p1, dg_gr_maxwell_lorentz_conf_1x1v_tensor_p2, NULL }, // 0
  { NULL, dg_gr_maxwell_lorentz_conf_1x2v_tensor_p1, dg_gr_maxwell_lorentz_conf_1x2v_tensor_p2, NULL }, // 1
  { NULL, dg_gr_maxwell_lorentz_conf_1x3v_tensor_p1, dg_gr_maxwell_lorentz_conf_1x3v_tensor_p2, NULL }, // 2
  { NULL, dg_gr_maxwell_lorentz_conf_2x2v_tensor_p1, dg_gr_maxwell_lorentz_conf_2x2v_tensor_p2, NULL }, // 3
  { NULL, GKYL_HYB_2X3V(dg_gr_maxwell_lorentz_conf_2x3v_tensor_p1), dg_gr_maxwell_lorentz_conf_2x3v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_3X3V(dg_gr_maxwell_lorentz_conf_3x3v_tensor_p1), dg_gr_maxwell_lorentz_conf_3x3v_tensor_p2, NULL }, // 5
};
