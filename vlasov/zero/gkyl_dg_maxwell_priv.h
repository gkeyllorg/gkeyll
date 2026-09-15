#pragma once

#include <gkyl_maxwell_kernels.h>
#include <gkyl_dg_gr_maxwell_kernels.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

#include <assert.h>

// private header for use in Maxwell DG equation object creation
// functions

// Types for various kernels
typedef double (*maxwell_surf_t)(const gkyl_maxwell_inp *meq, const double *w, const double *dx,
  const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r,
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);

typedef double (*maxwell_surf_from_flux_t)(const double *w, const double *dx,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out);

typedef void (*maxwell_gr_maxwell_vol_t)(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
  const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *det_h_nodal, const double *fields_no_J, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[4]; } gkyl_dg_maxwell_vol_kern_list;
typedef struct { maxwell_surf_t kernels[4]; } gkyl_dg_maxwell_surf_kern_list;
typedef struct { maxwell_gr_maxwell_vol_t kernels[4]; } gkyl_dg_maxwell_gr_maxwell_vol_kern_list;
typedef struct { maxwell_surf_from_flux_t kernels[4]; } gkyl_dg_maxwell_surf_from_flux_kern_list;

struct dg_maxwell {
  struct gkyl_dg_eqn eqn; // Base object
  gkyl_maxwell_inp maxwell_data; // Parameters needed by kernels
  gkyl_dg_gr_maxwell_inp gr_maxwell_data; // Parameters needed by GR-Maxwell kernels
  struct gkyl_range crange; // Configuration-space range for use in indexing conf_flux/jacob_pos
  const struct gkyl_array *jacob_pos; // Configuration-space position-map Jacobian (borrowed; never NULL:
                                      // callers without a nonuniform map pass an identity position map).

  maxwell_surf_t surf[3]; // pointers to surface kernels
  bool use_conf_flux_surf; // If using configuration-space fluxes for streaming terms (GR Maxwell).
  const struct gkyl_array *conf_flux_surf; // Nodal expansion of fluxes at configuration space surfaces. 
  maxwell_surf_from_flux_t surf_from_flux[3]; // Surface terms for streaming using fluxes.

  maxwell_gr_maxwell_vol_t vol; // pointers to volume kernels (for nodal geometry)
  const struct gkyl_surf_and_vol_node_arrays *lapse; // nodal expansion of lapse
  const struct gkyl_surf_and_vol_node_arrays *shift; // nodal expansion of shift
  const struct gkyl_surf_and_vol_node_arrays *h_ij; // nodal expansion of spatial metric
  const struct gkyl_surf_and_vol_node_arrays *h_ij_inv; // nodal expansion of inverse spatial metric
  const struct gkyl_surf_and_vol_node_arrays *det_h; // nodal expansion of spatial metric determinant
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_maxwell_vol_1x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  long cidx = gkyl_range_idx(&maxwell->crange, idx);
  const double *jacob_pos_d = (const double*) gkyl_array_cfetch(maxwell->jacob_pos, cidx);
  return maxwell_vol_1x_ser_p1(&maxwell->maxwell_data, xc, dx, jacob_pos_d, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_1x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  long cidx = gkyl_range_idx(&maxwell->crange, idx);
  const double *jacob_pos_d = (const double*) gkyl_array_cfetch(maxwell->jacob_pos, cidx);
  return maxwell_vol_1x_ser_p2(&maxwell->maxwell_data, xc, dx, jacob_pos_d, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_1x_ser_p3(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  long cidx = gkyl_range_idx(&maxwell->crange, idx);
  const double *jacob_pos_d = (const double*) gkyl_array_cfetch(maxwell->jacob_pos, cidx);
  return maxwell_vol_1x_ser_p3(&maxwell->maxwell_data, xc, dx, jacob_pos_d, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_2x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  long cidx = gkyl_range_idx(&maxwell->crange, idx);
  const double *jacob_pos_d = (const double*) gkyl_array_cfetch(maxwell->jacob_pos, cidx);
  return maxwell_vol_2x_ser_p1(&maxwell->maxwell_data, xc, dx, jacob_pos_d, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_2x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  long cidx = gkyl_range_idx(&maxwell->crange, idx);
  const double *jacob_pos_d = (const double*) gkyl_array_cfetch(maxwell->jacob_pos, cidx);
  return maxwell_vol_2x_ser_p2(&maxwell->maxwell_data, xc, dx, jacob_pos_d, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_2x_ser_p3(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  long cidx = gkyl_range_idx(&maxwell->crange, idx);
  const double *jacob_pos_d = (const double*) gkyl_array_cfetch(maxwell->jacob_pos, cidx);
  return maxwell_vol_2x_ser_p3(&maxwell->maxwell_data, xc, dx, jacob_pos_d, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_2x_tensor_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  long cidx = gkyl_range_idx(&maxwell->crange, idx);
  const double *jacob_pos_d = (const double*) gkyl_array_cfetch(maxwell->jacob_pos, cidx);
  return maxwell_vol_2x_tensor_p2(&maxwell->maxwell_data, xc, dx, jacob_pos_d, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_2x_tensor_p3(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  long cidx = gkyl_range_idx(&maxwell->crange, idx);
  const double *jacob_pos_d = (const double*) gkyl_array_cfetch(maxwell->jacob_pos, cidx);
  return maxwell_vol_2x_tensor_p3(&maxwell->maxwell_data, xc, dx, jacob_pos_d, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_3x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  long cidx = gkyl_range_idx(&maxwell->crange, idx);
  const double *jacob_pos_d = (const double*) gkyl_array_cfetch(maxwell->jacob_pos, cidx);
  return maxwell_vol_3x_ser_p1(&maxwell->maxwell_data, xc, dx, jacob_pos_d, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_3x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  long cidx = gkyl_range_idx(&maxwell->crange, idx);
  const double *jacob_pos_d = (const double*) gkyl_array_cfetch(maxwell->jacob_pos, cidx);
  return maxwell_vol_3x_ser_p2(&maxwell->maxwell_data, xc, dx, jacob_pos_d, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_maxwell_vol_3x_tensor_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  long cidx = gkyl_range_idx(&maxwell->crange, idx);
  const double *jacob_pos_d = (const double*) gkyl_array_cfetch(maxwell->jacob_pos, cidx);
  return maxwell_vol_3x_tensor_p2(&maxwell->maxwell_data, xc, dx, jacob_pos_d, qIn, qRhsOut);
}

// Volume kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_maxwell_vol_kern_list ser_vol_kernels[] = {
  { NULL, kernel_maxwell_vol_1x_ser_p1, kernel_maxwell_vol_1x_ser_p2, kernel_maxwell_vol_1x_ser_p3 }, // 0
  { NULL, kernel_maxwell_vol_2x_ser_p1, kernel_maxwell_vol_2x_ser_p2, kernel_maxwell_vol_2x_ser_p3 }, // 1
  { NULL, kernel_maxwell_vol_3x_ser_p1, kernel_maxwell_vol_3x_ser_p2, NULL },              // 2
};

// Volume kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_dg_maxwell_vol_kern_list ten_vol_kernels[] = {
  { NULL, kernel_maxwell_vol_1x_ser_p1, kernel_maxwell_vol_1x_ser_p2, kernel_maxwell_vol_1x_ser_p3 }, // 0
  { NULL, kernel_maxwell_vol_2x_ser_p1, kernel_maxwell_vol_2x_tensor_p2, kernel_maxwell_vol_2x_tensor_p3 }, // 1
  { NULL, kernel_maxwell_vol_3x_ser_p1, kernel_maxwell_vol_3x_tensor_p2, NULL },              // 2
};

// Volume kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_maxwell_gr_maxwell_vol_kern_list ser_gr_maxwell_vol_kernels[] = {
  { NULL, gr_maxwell_vol_1x_ser_p1, gr_maxwell_vol_1x_ser_p2, NULL }, // 0
  { NULL, gr_maxwell_vol_2x_ser_p1, gr_maxwell_vol_2x_ser_p2, NULL }, // 1
  { NULL, gr_maxwell_vol_3x_ser_p1, gr_maxwell_vol_3x_ser_p2, NULL },              // 2
};

// Volume kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_dg_maxwell_gr_maxwell_vol_kern_list ten_gr_maxwell_vol_kernels[] = {
  { NULL, gr_maxwell_vol_1x_tensor_p1, gr_maxwell_vol_1x_tensor_p2, NULL }, // 0
  { NULL, gr_maxwell_vol_2x_tensor_p1, gr_maxwell_vol_2x_tensor_p2, NULL }, // 1
  { NULL, gr_maxwell_vol_3x_tensor_p1, gr_maxwell_vol_3x_tensor_p2, NULL },              // 2
};

// Surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, maxwell_surfx_1x_ser_p1, maxwell_surfx_1x_ser_p2, maxwell_surfx_1x_ser_p3 }, // 0
  { NULL, maxwell_surfx_2x_ser_p1, maxwell_surfx_2x_ser_p2, maxwell_surfx_2x_ser_p3 }, // 1
  { NULL, maxwell_surfx_3x_ser_p1, maxwell_surfx_3x_ser_p2, NULL },                 // 2
};

// Surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ten_surf_x_kernels[] = {
  { NULL, maxwell_surfx_1x_ser_p1, maxwell_surfx_1x_ser_p2, maxwell_surfx_1x_ser_p3 }, // 0
  { NULL, maxwell_surfx_2x_ser_p1, maxwell_surfx_2x_tensor_p2, maxwell_surfx_2x_tensor_p3 }, // 1
  { NULL, maxwell_surfx_3x_ser_p1, maxwell_surfx_3x_tensor_p2, NULL },                 // 2
};

// Surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, maxwell_surfy_2x_ser_p1, maxwell_surfy_2x_ser_p2, maxwell_surfy_2x_ser_p3 }, // 1
  { NULL, maxwell_surfy_3x_ser_p1, maxwell_surfy_3x_ser_p2, NULL },                 // 2
};

// Surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ten_surf_y_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, maxwell_surfy_2x_ser_p1, maxwell_surfy_2x_tensor_p2, maxwell_surfy_2x_tensor_p3 }, // 1
  { NULL, maxwell_surfy_3x_ser_p1, maxwell_surfy_3x_tensor_p2, NULL },                 // 2
};

// Surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL, NULL },                 // 0
  { NULL, NULL, NULL, NULL },                 // 1
  { NULL, maxwell_surfz_3x_ser_p1, maxwell_surfz_3x_ser_p2, NULL }, // 2
};

// Surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ten_surf_z_kernels[] = {
  { NULL, NULL, NULL, NULL },                 // 0
  { NULL, NULL, NULL, NULL },                 // 1
  { NULL, maxwell_surfz_3x_ser_p1, maxwell_surfz_3x_tensor_p2, NULL }, // 2
};

// DG GR Maxwell surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_from_flux_kern_list ser_gr_maxwell_surf_from_flux_x_kernels[] = {
  // 1x kernels
  { NULL, gr_maxwell_surfx_1x_ser_p1, gr_maxwell_surfx_1x_ser_p2, NULL }, // 0
  // 2x kernels
  { NULL, gr_maxwell_surfx_2x_ser_p1, gr_maxwell_surfx_2x_ser_p2, NULL }, // 1
  // 3x kernels
  { NULL, gr_maxwell_surfx_3x_ser_p1, gr_maxwell_surfx_3x_ser_p2, NULL }, // 2
};

// DG GR Maxwell surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_from_flux_kern_list ser_gr_maxwell_surf_from_flux_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, gr_maxwell_surfy_2x_ser_p1, gr_maxwell_surfy_2x_ser_p2, NULL }, // 1
  // 3x kernels
  { NULL, gr_maxwell_surfy_3x_ser_p1, gr_maxwell_surfy_3x_ser_p2, NULL }, // 2
};

// DG GR Maxwell surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_from_flux_kern_list ser_gr_maxwell_surf_from_flux_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
// 2x kernels
  { NULL, NULL, NULL, NULL }, // 1
  // 3x kernels
  { NULL, gr_maxwell_surfz_3x_ser_p1, gr_maxwell_surfz_3x_ser_p2, NULL }, // 2
};


// DG GR Maxwell surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_from_flux_kern_list tensor_gr_maxwell_surf_from_flux_x_kernels[] = {
  // 1x kernels
  { NULL, gr_maxwell_surfx_1x_tensor_p1, gr_maxwell_surfx_1x_tensor_p2, NULL }, // 0
  // 2x kernels
  { NULL, gr_maxwell_surfx_2x_tensor_p1, gr_maxwell_surfx_2x_tensor_p2, NULL }, // 1
  // 3x kernels
  { NULL, gr_maxwell_surfx_3x_tensor_p1, gr_maxwell_surfx_3x_tensor_p2, NULL }, // 2
};

// DG GR Maxwell surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_from_flux_kern_list tensor_gr_maxwell_surf_from_flux_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, gr_maxwell_surfy_2x_tensor_p1, gr_maxwell_surfy_2x_tensor_p2, NULL }, // 1
  // 3x kernels
  { NULL, gr_maxwell_surfy_3x_tensor_p1, gr_maxwell_surfy_3x_tensor_p2, NULL }, // 2
};

// DG GR Maxwell surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_maxwell_surf_from_flux_kern_list tensor_gr_maxwell_surf_from_flux_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
// 2x kernels
  { NULL, NULL, NULL, NULL }, // 1
  // 3x kernels
  { NULL, gr_maxwell_surfz_3x_tensor_p1, gr_maxwell_surfz_3x_tensor_p2, NULL }, // 2
};



/**
 * Free Maxwell equation object
 *
 * @param ref Reference counter for Maxwell equation
 */
void gkyl_maxwell_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn *eqn, const double* xc, const double*  dx,
  const int*  idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);

  // For GR problems, use conf_flux_surf
  if (maxwell->use_conf_flux_surf) {
    // Each cell owns the *lower* edge surface flux
    long cidx = gkyl_range_idx(&maxwell->crange, idx);
    const double* lapse = (const double*) gkyl_array_cfetch(maxwell->lapse->nodal_arr_vol, cidx);
    const double* shift = (const double*) gkyl_array_cfetch(maxwell->shift->nodal_arr_vol, cidx);
    const double* h_ij = (const double*) gkyl_array_cfetch(maxwell->h_ij->nodal_arr_vol, cidx);
    const double* h_ij_inv = (const double*) gkyl_array_cfetch(maxwell->h_ij_inv->nodal_arr_vol, cidx);
    const double* det_h = (const double*) gkyl_array_cfetch(maxwell->det_h->nodal_arr_vol, cidx);
    const double* jacob_pos_d = (const double*) gkyl_array_cfetch(maxwell->jacob_pos, cidx);

    // For GR Maxwell the volume term does not contribute to CFL, only the conf-flux.
    maxwell->vol(&maxwell->gr_maxwell_data, xc, dx, jacob_pos_d, lapse, shift, h_ij, h_ij_inv, det_h, qIn, qRhsOut);

    return 0.0;
  }
  else {
    // For non-GR problems, use regular volume kernel (which includes streaming terms)
    assert(false);
  }
}

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);

  // For GR problems, use conf_flux_surf
  if (maxwell->use_conf_flux_surf) {
    // Each cell owns the *lower* edge surface flux
    long cidxC = gkyl_range_idx(&maxwell->crange, idxC);
    long cidxR = gkyl_range_idx(&maxwell->crange, idxR);
    const double* conf_flux_surf_l = (const double*) gkyl_array_cfetch(maxwell->conf_flux_surf, cidxC);
    const double* conf_flux_surf_r = (const double*) gkyl_array_cfetch(maxwell->conf_flux_surf, cidxR);

    return maxwell->surf_from_flux[dir](xcC, dxC,
      conf_flux_surf_l, conf_flux_surf_r, qRhsOut);
  }
  else {
    // Position-map Jacobian of the left, center, and right cells; the surface
    // kernel divides each side's stored J*field by its own normal-direction
    // Jacobian to recover the physical field for the interface flux. With no
    // position map this is the owned identity array (J = 1), bit-identical to
    // the uniform-grid kernels.
    long cidxL = gkyl_range_idx(&maxwell->crange, idxL);
    long cidxC = gkyl_range_idx(&maxwell->crange, idxC);
    long cidxR = gkyl_range_idx(&maxwell->crange, idxR);
    const double *jpl = (const double*) gkyl_array_cfetch(maxwell->jacob_pos, cidxL);
    const double *jpc = (const double*) gkyl_array_cfetch(maxwell->jacob_pos, cidxC);
    const double *jpr = (const double*) gkyl_array_cfetch(maxwell->jacob_pos, cidxR);
    return maxwell->surf[dir](&maxwell->maxwell_data, xcC, dxC,
      jpl, jpc, jpr, qInL, qInC, qInR, qRhsOut);
  }
}

GKYL_CU_D
static double
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  return 0.;
}
