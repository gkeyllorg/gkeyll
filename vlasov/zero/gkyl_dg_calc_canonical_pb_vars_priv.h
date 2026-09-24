// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_canonical_pb.h>
#include <gkyl_canonical_pb_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef int (*canonical_pb_alpha_surf_t)(const double *w, const double *dxv, const double *hamil,
  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
typedef void (*canonical_pb_m1i_contra_to_cov_t)(const double *h_ij, const double *v_i, const double *nv_i,
 double* GKYL_RESTRICT v_i_cov, double* GKYL_RESTRICT nv_i_cov); 
typedef void (*canonical_pb_pressure_t)(const double *h_ij_inv, const double *MEnergy, const double *v_i,
  const double *nv_i, double* GKYL_RESTRICT d_Jv_P); 

// for use in kernel tables
typedef struct { canonical_pb_alpha_surf_t kernels[3]; } gkyl_dg_canonical_pb_alpha_surf_kern_list;
typedef struct { canonical_pb_m1i_contra_to_cov_t kernels[3]; } gkyl_dg_canonical_pb_m1i_contra_to_cov_kern_list;
typedef struct { canonical_pb_pressure_t kernels[3]; } gkyl_dg_canonical_pb_pressure_kern_list;

struct gkyl_dg_calc_canonical_pb_vars {
  struct gkyl_rect_grid phase_grid; // Phase space grid for cell spacing and cell center
  int cdim; // Configuration space dimensionality
  int pdim; // Phase space dimensionality
  canonical_pb_alpha_surf_t alpha_surf[6]; // kernel for computing surface expansion of phase space flux alpha
  canonical_pb_alpha_surf_t alpha_edge_surf[3]; // kernel for computing surface expansion of phase space flux alpha
                                               // at upper configuration space edge
  canonical_pb_m1i_contra_to_cov_t canonical_pb_covariant_u_i; // Canonical pb covariant u_i components
  canonical_pb_pressure_t canonical_pb_pressure; // Canonical pb pressure
  uint32_t flags;
  struct gkyl_dg_calc_canonical_pb_vars *on_dev; // pointer to itself or device data
};

// The cv_index[cd].vdim[cd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

//
// Serendipity surface kernels general geometry
//
// canonical_pb general geometry phase space flux alpha surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list ser_canonical_pb_alpha_surfx_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_alpha_surfx_1x1v_ser_p1, canonical_pb_alpha_surfx_1x1v_ser_p2 }, // 0
  { NULL, canonical_pb_alpha_surfx_1x2v_ser_p1, canonical_pb_alpha_surfx_1x2v_ser_p2 }, // 1
  { NULL, canonical_pb_alpha_surfx_1x3v_ser_p1, canonical_pb_alpha_surfx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_alpha_surfx_2x2v_ser_p1, canonical_pb_alpha_surfx_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_alpha_surfx_2x3v_ser_p1, canonical_pb_alpha_surfx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// canonical_pb general geometry phase space flux alpha edge surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list ser_canonical_pb_alpha_edge_surfx_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_alpha_edge_surfx_1x1v_ser_p1, canonical_pb_alpha_edge_surfx_1x1v_ser_p2 }, // 0
  { NULL, canonical_pb_alpha_edge_surfx_1x2v_ser_p1, canonical_pb_alpha_edge_surfx_1x2v_ser_p2 }, // 1
  { NULL, canonical_pb_alpha_edge_surfx_1x3v_ser_p1, canonical_pb_alpha_edge_surfx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_alpha_edge_surfx_2x2v_ser_p1, canonical_pb_alpha_edge_surfx_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_alpha_edge_surfx_2x3v_ser_p1, canonical_pb_alpha_edge_surfx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// canonical_pb general geometry phase space flux alpha surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list ser_canonical_pb_alpha_surfy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, canonical_pb_alpha_surfy_2x2v_ser_p1, canonical_pb_alpha_surfy_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_alpha_surfy_2x3v_ser_p1, canonical_pb_alpha_surfy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// canonical_pb general geometry phase space flux alpha edge surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list ser_canonical_pb_alpha_edge_surfy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, canonical_pb_alpha_edge_surfy_2x2v_ser_p1, canonical_pb_alpha_edge_surfy_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_alpha_edge_surfy_2x3v_ser_p1, canonical_pb_alpha_edge_surfy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// canonical_pb general geometry phase space flux alpha surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list ser_canonical_pb_alpha_surfz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// canonical_pb general geometry phase space flux alpha edge surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list ser_canonical_pb_alpha_edge_surfz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

//
// Serendipity surface kernels general geometry
//
// canonical_pb general geometry phase space flux alpha surface expansions in vx (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list ser_canonical_pb_alpha_surfvx_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_alpha_surfvx_1x1v_ser_p1, canonical_pb_alpha_surfvx_1x1v_ser_p2 }, // 0
  { NULL, canonical_pb_alpha_surfvx_1x2v_ser_p1, canonical_pb_alpha_surfvx_1x2v_ser_p2 }, // 1
  { NULL, canonical_pb_alpha_surfvx_1x3v_ser_p1, canonical_pb_alpha_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_alpha_surfvx_2x2v_ser_p1, canonical_pb_alpha_surfvx_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_alpha_surfvx_2x3v_ser_p1, canonical_pb_alpha_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};


// canonical_pb general geometry phase space flux alpha surface expansions in vy (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list ser_canonical_pb_alpha_surfvy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, canonical_pb_alpha_surfvy_1x2v_ser_p1, canonical_pb_alpha_surfvy_1x2v_ser_p2 }, // 1
  { NULL, canonical_pb_alpha_surfvy_1x3v_ser_p1, canonical_pb_alpha_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_alpha_surfvy_2x2v_ser_p1, canonical_pb_alpha_surfvy_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_alpha_surfvy_2x3v_ser_p1, canonical_pb_alpha_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// canonical_pb general geometry phase space flux alpha surface expansions in vz (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list ser_canonical_pb_alpha_surfvz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, canonical_pb_alpha_surfvz_1x3v_ser_p1, canonical_pb_alpha_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, canonical_pb_alpha_surfvz_2x3v_ser_p1, canonical_pb_alpha_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// canonical_pb contravaraint to covariant conversion (Serendipity kernels)
// (Jnu_i = h_{ij}Jnu^i and u_i = h_{ij}u^j) 
GKYL_CU_D
static const gkyl_dg_canonical_pb_m1i_contra_to_cov_kern_list ser_canonical_pb_m1i_contra_to_cov_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_vars_m1i_contra_to_cov_1x1v_ser_p1, canonical_pb_vars_m1i_contra_to_cov_1x1v_ser_p2 }, // 0
  { NULL, canonical_pb_vars_m1i_contra_to_cov_1x2v_ser_p1, canonical_pb_vars_m1i_contra_to_cov_1x2v_ser_p2 }, // 1
  { NULL, canonical_pb_vars_m1i_contra_to_cov_1x3v_ser_p1, canonical_pb_vars_m1i_contra_to_cov_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_vars_m1i_contra_to_cov_2x2v_ser_p1, canonical_pb_vars_m1i_contra_to_cov_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_vars_m1i_contra_to_cov_2x3v_ser_p1, canonical_pb_vars_m1i_contra_to_cov_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};


// canonical_pb Pressure (d*P*Jv = 2*E - n*h^{ij}*u_i*u_j) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_pressure_kern_list ser_canonical_pb_pressure_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_vars_pressure_1x1v_ser_p1, canonical_pb_vars_pressure_1x1v_ser_p2 }, // 0
  { NULL, canonical_pb_vars_pressure_1x2v_ser_p1, canonical_pb_vars_pressure_1x2v_ser_p2 }, // 1
  { NULL, canonical_pb_vars_pressure_1x3v_ser_p1, canonical_pb_vars_pressure_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_vars_pressure_2x2v_ser_p1, canonical_pb_vars_pressure_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_vars_pressure_2x3v_ser_p1, canonical_pb_vars_pressure_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

//
// Tensor surface kernels general geometry
//
// canonical_pb general geometry phase space flux alpha surface expansions in x (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list tensor_canonical_pb_alpha_surfx_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_alpha_surfx_1x1v_tensor_p1, canonical_pb_alpha_surfx_1x1v_tensor_p2 }, // 0
  { NULL, canonical_pb_alpha_surfx_1x2v_tensor_p1, canonical_pb_alpha_surfx_1x2v_tensor_p2 }, // 1
  { NULL, canonical_pb_alpha_surfx_1x3v_tensor_p1, canonical_pb_alpha_surfx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_alpha_surfx_2x2v_tensor_p1, canonical_pb_alpha_surfx_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_alpha_surfx_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_alpha_surfx_3x3v_tensor_p1, NULL }, // 5
};

// canonical_pb general geometry phase space flux alpha edge surface expansions in x (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list tensor_canonical_pb_alpha_edge_surfx_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_alpha_edge_surfx_1x1v_tensor_p1, canonical_pb_alpha_edge_surfx_1x1v_tensor_p2 }, // 0
  { NULL, canonical_pb_alpha_edge_surfx_1x2v_tensor_p1, canonical_pb_alpha_edge_surfx_1x2v_tensor_p2 }, // 1
  { NULL, canonical_pb_alpha_edge_surfx_1x3v_tensor_p1, canonical_pb_alpha_edge_surfx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_alpha_edge_surfx_2x2v_tensor_p1, canonical_pb_alpha_edge_surfx_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_alpha_edge_surfx_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_alpha_edge_surfx_3x3v_tensor_p1, NULL }, // 5
};

// canonical_pb general geometry phase space flux alpha surface expansions in y (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list tensor_canonical_pb_alpha_surfy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, canonical_pb_alpha_surfy_2x2v_tensor_p1, canonical_pb_alpha_surfy_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_alpha_surfy_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_alpha_surfy_3x3v_tensor_p1, NULL }, // 5
};

// canonical_pb general geometry phase space flux alpha edge surface expansions in y (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list tensor_canonical_pb_alpha_edge_surfy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, canonical_pb_alpha_edge_surfy_2x2v_tensor_p1, canonical_pb_alpha_edge_surfy_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_alpha_edge_surfy_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_alpha_edge_surfy_3x3v_tensor_p1, NULL }, // 5
};

// canonical_pb general geometry phase space flux alpha surface expansions in z (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list tensor_canonical_pb_alpha_surfz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_alpha_surfz_3x3v_tensor_p1, NULL }, // 5
};

// canonical_pb general geometry phase space flux alpha edge surface expansions in z (tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list tensor_canonical_pb_alpha_edge_surfz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_alpha_edge_surfz_3x3v_tensor_p1, NULL }, // 5
};

//
// Tensor surface kernels general geometry
//
// canonical_pb general geometry phase space flux alpha surface expansions in vx (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list tensor_canonical_pb_alpha_surfvx_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_alpha_surfvx_1x1v_tensor_p1, canonical_pb_alpha_surfvx_1x1v_tensor_p2 }, // 0
  { NULL, canonical_pb_alpha_surfvx_1x2v_tensor_p1, canonical_pb_alpha_surfvx_1x2v_tensor_p2 }, // 1
  { NULL, canonical_pb_alpha_surfvx_1x3v_tensor_p1, canonical_pb_alpha_surfvx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_alpha_surfvx_2x2v_tensor_p1, canonical_pb_alpha_surfvx_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_alpha_surfvx_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_alpha_surfvx_3x3v_tensor_p1, NULL }, // 5
};


// canonical_pb general geometry phase space flux alpha surface expansions in vy (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list tensor_canonical_pb_alpha_surfvy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, canonical_pb_alpha_surfvy_1x2v_tensor_p1, canonical_pb_alpha_surfvy_1x2v_tensor_p2 }, // 1
  { NULL, canonical_pb_alpha_surfvy_1x3v_tensor_p1, canonical_pb_alpha_surfvy_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_alpha_surfvy_2x2v_tensor_p1, canonical_pb_alpha_surfvy_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_alpha_surfvy_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_alpha_surfvy_3x3v_tensor_p1, NULL }, // 5
};

// canonical_pb general geometry phase space flux alpha surface expansions in vz (tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_alpha_surf_kern_list tensor_canonical_pb_alpha_surfvz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, canonical_pb_alpha_surfvz_1x3v_tensor_p1, canonical_pb_alpha_surfvz_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, canonical_pb_alpha_surfvz_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_alpha_surfvz_3x3v_tensor_p1, NULL }, // 5
};

// canonical_pb contravaraint to covariant conversion (Tensor kernels)
// (Jnu_i = h_{ij}Jnu^i and u_i = h_{ij}u^j) 
GKYL_CU_D
static const gkyl_dg_canonical_pb_m1i_contra_to_cov_kern_list tensor_canonical_pb_m1i_contra_to_cov_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_vars_m1i_contra_to_cov_1x1v_tensor_p1, canonical_pb_vars_m1i_contra_to_cov_1x1v_tensor_p2 }, // 0
  { NULL, canonical_pb_vars_m1i_contra_to_cov_1x2v_tensor_p1, canonical_pb_vars_m1i_contra_to_cov_1x2v_tensor_p2 }, // 1
  { NULL, canonical_pb_vars_m1i_contra_to_cov_1x3v_tensor_p1, canonical_pb_vars_m1i_contra_to_cov_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_vars_m1i_contra_to_cov_2x2v_tensor_p1, canonical_pb_vars_m1i_contra_to_cov_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_vars_m1i_contra_to_cov_2x3v_tensor_p1, NULL }, //4
  // 3x kernels
  { NULL, canonical_pb_vars_m1i_contra_to_cov_3x3v_tensor_p1, NULL }, // 5
};


// canonical_pb Pressure (d*P*Jv = h^{ij}*M2_{ij} - n*h^{ij}*u_i*u_j) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_pressure_kern_list tensor_canonical_pb_pressure_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_vars_pressure_1x1v_tensor_p1, canonical_pb_vars_pressure_1x1v_tensor_p2 }, // 0
  { NULL, canonical_pb_vars_pressure_1x2v_tensor_p1, canonical_pb_vars_pressure_1x2v_tensor_p2 }, // 1
  { NULL, canonical_pb_vars_pressure_1x3v_tensor_p1, canonical_pb_vars_pressure_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_vars_pressure_2x2v_tensor_p1, canonical_pb_vars_pressure_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_vars_pressure_2x3v_tensor_p1, NULL }, //4
  // 3x kernels
  { NULL, canonical_pb_vars_pressure_3x3v_tensor_p1, NULL }, // 5
};


GKYL_CU_D
static canonical_pb_alpha_surf_t
choose_canonical_pb_alpha_surf_kern(enum gkyl_basis_type b_type, int dir, int cv_index, int cdim, int vdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      // Verify that the poly-order is 2 for ser case
      assert(poly_order == 2);
      if (dir == 0)
        return ser_canonical_pb_alpha_surfx_kernels[cv_index].kernels[poly_order];
      else if (dir == 1)
        return ser_canonical_pb_alpha_surfy_kernels[cv_index].kernels[poly_order];
      else if (dir == 2)
        return ser_canonical_pb_alpha_surfz_kernels[cv_index].kernels[poly_order];
      else
        return NULL;
      break;
    case GKYL_BASIS_MODAL_HYBRID:
      // Verify that the poly-order is 1 for hybrid case
      assert(poly_order == 1);
      if (dir == 0)
        return ser_canonical_pb_alpha_surfx_kernels[cv_index].kernels[poly_order];
      else if (dir == 1)
        return ser_canonical_pb_alpha_surfy_kernels[cv_index].kernels[poly_order];
      else if (dir == 2)
        return ser_canonical_pb_alpha_surfz_kernels[cv_index].kernels[poly_order];
      else
        return NULL;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      if (dir == 0)
        return tensor_canonical_pb_alpha_surfx_kernels[cv_index].kernels[poly_order];
      else if (dir == 1)
        return tensor_canonical_pb_alpha_surfy_kernels[cv_index].kernels[poly_order];
      else if (dir == 2)
        return tensor_canonical_pb_alpha_surfz_kernels[cv_index].kernels[poly_order];
      else
        return NULL;
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static canonical_pb_alpha_surf_t
choose_canonical_pb_alpha_edge_surf_kern(enum gkyl_basis_type b_type, int dir, int cv_index, int cdim, int vdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      // Verify that the poly-order is 2 for ser case
      assert(poly_order == 2);
      if (dir == 0)
        return ser_canonical_pb_alpha_edge_surfx_kernels[cv_index].kernels[poly_order];
      else if (dir == 1)
        return ser_canonical_pb_alpha_edge_surfy_kernels[cv_index].kernels[poly_order];
      else if (dir == 2)
        return ser_canonical_pb_alpha_edge_surfz_kernels[cv_index].kernels[poly_order];
      else
        return NULL;
      break;
    case GKYL_BASIS_MODAL_HYBRID:
      // Verify that the poly-order is 1 for hybrid case
      assert(poly_order == 1);
      if (dir == 0)
        return ser_canonical_pb_alpha_edge_surfx_kernels[cv_index].kernels[poly_order];
      else if (dir == 1)
        return ser_canonical_pb_alpha_edge_surfy_kernels[cv_index].kernels[poly_order];
      else if (dir == 2)
        return ser_canonical_pb_alpha_edge_surfz_kernels[cv_index].kernels[poly_order];
      else
        return NULL;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      if (dir == 0)
        return tensor_canonical_pb_alpha_edge_surfx_kernels[cv_index].kernels[poly_order];
      else if (dir == 1)
        return tensor_canonical_pb_alpha_edge_surfy_kernels[cv_index].kernels[poly_order];
      else if (dir == 2)
        return tensor_canonical_pb_alpha_edge_surfz_kernels[cv_index].kernels[poly_order];
      else
        return NULL;
      break;
    default:
      assert(false);
      break;  
  }
}


GKYL_CU_D
static canonical_pb_alpha_surf_t
choose_canonical_pb_alpha_surf_v_kern(enum gkyl_basis_type b_type, int dir, int cv_index, int cdim, int vdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      // Verify that the poly-order is 2 for ser case
      assert(poly_order == 2);
      if (dir == 0)
        return ser_canonical_pb_alpha_surfvx_kernels[cv_index].kernels[poly_order];
      else if (dir == 1)
        return ser_canonical_pb_alpha_surfvy_kernels[cv_index].kernels[poly_order];
      else if (dir == 2)
        return ser_canonical_pb_alpha_surfvz_kernels[cv_index].kernels[poly_order];
      else
        return NULL;
      break;
    case GKYL_BASIS_MODAL_HYBRID:
      // Verify that the poly-order is 1 for hybrid case
      assert(poly_order == 1);
      if (dir == 0)
        return ser_canonical_pb_alpha_surfvx_kernels[cv_index].kernels[poly_order];
      else if (dir == 1)
        return ser_canonical_pb_alpha_surfvy_kernels[cv_index].kernels[poly_order];
      else if (dir == 2)
        return ser_canonical_pb_alpha_surfvz_kernels[cv_index].kernels[poly_order];
      else
        return NULL;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      if (dir == 0)
        return tensor_canonical_pb_alpha_surfvx_kernels[cv_index].kernels[poly_order];
      else if (dir == 1)
        return tensor_canonical_pb_alpha_surfvy_kernels[cv_index].kernels[poly_order];
      else if (dir == 2)
        return tensor_canonical_pb_alpha_surfvz_kernels[cv_index].kernels[poly_order];
      else
        return NULL;
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static canonical_pb_m1i_contra_to_cov_t
choose_canonical_pb_m1i_contra_to_cov_kern(enum gkyl_basis_type b_type, int cv_index, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      // Verify that the poly-order is 2 for ser case
      assert(poly_order == 2);
      return ser_canonical_pb_m1i_contra_to_cov_kernels[cv_index].kernels[poly_order];
      break; 
    case GKYL_BASIS_MODAL_HYBRID:
      // Verify that the poly-order is 1 for hybrid case
      assert(poly_order == 1);
      return ser_canonical_pb_m1i_contra_to_cov_kernels[cv_index].kernels[poly_order];
      break; 
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_canonical_pb_m1i_contra_to_cov_kernels[cv_index].kernels[poly_order];
      break; 
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static canonical_pb_pressure_t
choose_canonical_pb_pressure_kern(enum gkyl_basis_type b_type, int cv_index, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      // Verify that the poly-order is 2 for ser case
      assert(poly_order == 2);
      return ser_canonical_pb_pressure_kernels[cv_index].kernels[poly_order];
      break; 
    case GKYL_BASIS_MODAL_HYBRID:
      // Verify that the poly-order is 1 for hybrid case
      assert(poly_order == 1);
      return ser_canonical_pb_pressure_kernels[cv_index].kernels[poly_order];
      break; 
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_canonical_pb_pressure_kernels[cv_index].kernels[poly_order];
      break; 
    default:
      assert(false);
      break;  
  }
}