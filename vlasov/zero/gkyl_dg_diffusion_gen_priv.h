#pragma once

#include <gkyl_dg_diffusion_gen_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in diffusion_gen DG equation object creation
// functions

// Types for various kernels
typedef double (*diffusion_gen_surf_t)(const double *w, const double *dx,
  const double* D, const double *q[27],
  double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_diffusion_gen_vol_kern_list;
typedef struct { diffusion_gen_surf_t kernels[3]; } gkyl_dg_diffusion_gen_surf_kern_list;

struct dg_diffusion_gen {
  struct gkyl_dg_eqn eqn;
  diffusion_gen_surf_t surf[3][3];
  struct gkyl_range conf_range;
  struct gkyl_dg_diffusion_gen_auxfields auxfields;
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_dg_diffusion_gen_vol_2x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gen* diffusion_gen = container_of(eqn, struct dg_diffusion_gen, eqn);
  
  long cidx = gkyl_range_idx(&diffusion_gen->conf_range, idx);
  
  return dg_diffusion_gen_vol_2x_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion_gen->auxfields.Dij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_gen_vol_2x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gen* diffusion_gen = container_of(eqn, struct dg_diffusion_gen, eqn);
  
  long cidx = gkyl_range_idx(&diffusion_gen->conf_range, idx);
  
  return dg_diffusion_gen_vol_2x_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion_gen->auxfields.Dij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_gen_vol_3x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gen* diffusion_gen = container_of(eqn, struct dg_diffusion_gen, eqn);
  
  long cidx = gkyl_range_idx(&diffusion_gen->conf_range, idx);
  
  return dg_diffusion_gen_vol_3x_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion_gen->auxfields.Dij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_gen_vol_3x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gen* diffusion_gen = container_of(eqn, struct dg_diffusion_gen, eqn);
  
  long cidx = gkyl_range_idx(&diffusion_gen->conf_range, idx);
  
  return dg_diffusion_gen_vol_3x_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion_gen->auxfields.Dij, cidx),
    qIn, qRhsOut);
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_diffusion_gen_vol_kern_list ser_vol_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, kernel_dg_diffusion_gen_vol_2x_ser_p1, kernel_dg_diffusion_gen_vol_2x_ser_p2 },
  { NULL, kernel_dg_diffusion_gen_vol_3x_ser_p1, kernel_dg_diffusion_gen_vol_3x_ser_p2 },
};

// Surface kernel list: xx-direction
GKYL_CU_D
static const gkyl_dg_diffusion_gen_surf_kern_list ser_surf_xx_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, dg_diffusion_gen_surfxx_2x_ser_p1, dg_diffusion_gen_surfxx_2x_ser_p2 },
  { NULL, dg_diffusion_gen_surfxx_3x_ser_p1, dg_diffusion_gen_surfxx_3x_ser_p2 },
};

// Surface kernel list: xy-direction
GKYL_CU_D
static const gkyl_dg_diffusion_gen_surf_kern_list ser_surf_xy_kernels[] = {
  { NULL, NULL, NULL },// general diffusion is not implemented for 1D, use diffusion instead
  { NULL, dg_diffusion_gen_surfxy_2x_ser_p1, dg_diffusion_gen_surfxy_2x_ser_p2 },
  { NULL, dg_diffusion_gen_surfxy_3x_ser_p1, dg_diffusion_gen_surfxy_3x_ser_p2 },
};

// Surface kernel list: xz-direction
GKYL_CU_D
static const gkyl_dg_diffusion_gen_surf_kern_list ser_surf_xz_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion_gen_surfxz_3x_ser_p1, dg_diffusion_gen_surfxz_3x_ser_p2 },
};

// Surface kernel list: yx-direction
GKYL_CU_D
static const gkyl_dg_diffusion_gen_surf_kern_list ser_surf_yx_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, dg_diffusion_gen_surfyx_2x_ser_p1, dg_diffusion_gen_surfyx_2x_ser_p2 },
  { NULL, dg_diffusion_gen_surfyx_3x_ser_p1, dg_diffusion_gen_surfyx_3x_ser_p2 },
};

// Surface kernel list: yy-direction
GKYL_CU_D
static const gkyl_dg_diffusion_gen_surf_kern_list ser_surf_yy_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, dg_diffusion_gen_surfyy_2x_ser_p1, dg_diffusion_gen_surfyy_2x_ser_p2 },
  { NULL, dg_diffusion_gen_surfyy_3x_ser_p1, dg_diffusion_gen_surfyy_3x_ser_p2 },
};

// Surface kernel list: yz-direction
GKYL_CU_D
static const gkyl_dg_diffusion_gen_surf_kern_list ser_surf_yz_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion_gen_surfyz_3x_ser_p1, dg_diffusion_gen_surfyz_3x_ser_p2 },
};

// Surface kernel list: zx-direction
GKYL_CU_D
static const gkyl_dg_diffusion_gen_surf_kern_list ser_surf_zx_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion_gen_surfzx_3x_ser_p1, dg_diffusion_gen_surfzx_3x_ser_p2 },
};

// Surface kernel list: zy-direction
GKYL_CU_D
static const gkyl_dg_diffusion_gen_surf_kern_list ser_surf_zy_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion_gen_surfzy_3x_ser_p1, dg_diffusion_gen_surfzy_3x_ser_p2 },
};

// Surface kernel list: zz-direction
GKYL_CU_D
static const gkyl_dg_diffusion_gen_surf_kern_list ser_surf_zz_kernels[] = {
  { NULL, NULL, NULL }, // general diffusion is not implemented for 1D, use diffusion instead
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion_gen_surfzz_3x_ser_p1, dg_diffusion_gen_surfzz_3x_ser_p2 },
};

/**
 * Free diffusion_gen equation object
 *
 * @param ref Reference counter for constant diffusion_gen equation
 */
void gkyl_diffusion_gen_free(const struct gkyl_ref_count* ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn* eqn, int dir1, int dir2,
  const double* xc, const double* dxc, const int* idxc,
  long sz_dim, const int idx[27][GKYL_MAX_DIM], const double* qIn[27],
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gen* diffusion_gen = container_of(eqn, struct dg_diffusion_gen, eqn);
  long cidx = gkyl_range_idx(&diffusion_gen->conf_range, idxc);
  
  return diffusion_gen->surf[dir1][dir2](xc, dxc,
    (const double*) gkyl_array_cfetch(diffusion_gen->auxfields.Dij, cidx), 
    qIn, qRhsOut);
}
