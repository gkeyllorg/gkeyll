// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_collisionless_flux_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef double (*gk_collisionless_flux_surf_t)(const double *w, const double *dxv, 
  const double *vmap, const double *vmapSq, const double q_, const double m_, 
  const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
  const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR, const double *phi,
  const double *JfL, const double *JfR, double* GKYL_RESTRICT flux_surf); 

typedef double (*gk_collisionless_flux_surfvpar_t)( 
  const double *w, const double *dxv, 
  const double *vmap_prime_l, const double *vmap_prime_r,
  const double *vmap, const double *vmapSq, const double q_, const double m_, 
  const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
  const double *bmag, const double *phi, const double *JfL, const double *JfR, double* GKYL_RESTRICT flux_surf); 

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below.
GKYL_CU_D
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense.
  {-1,  0,  1}, // 1x kernel indices.
  {-1, -1,  2}, // 2x kernel indices.
  {-1, -1,  3}, // 3x kernel indices.
};

// for use in kernel tables
typedef struct { gk_collisionless_flux_surf_t kernels[3]; } gkyl_gk_collisionless_flux_surf_kern_list;
typedef struct { gk_collisionless_flux_surfvpar_t kernels[3]; } gkyl_gk_collisionless_flux_surfvpar_kern_list;

struct gkyl_gk_collisionless_flux {
  struct gkyl_rect_grid phase_grid; // Phase space grid for cell spacing and cell center
  int cdim; // Configuration space dimensionality
  int pdim; // Phase space dimensionality
  gk_collisionless_flux_surf_t flux_surf[GKYL_MAX_CDIM]; // kernel for computing surface expansion of phase space flux.
  gk_collisionless_flux_surf_t flux_surf_edge_lo[GKYL_MAX_CDIM]; // kernel for computing surface expansion of phase space flux.
  gk_collisionless_flux_surf_t flux_surf_edge_up[GKYL_MAX_CDIM]; // kernel for computing surface expansion of phase space flux.
                                                        // at upper configuration space edge
  gk_collisionless_flux_surfvpar_t flux_surfvpar[1]; // kernel for computing surface expansion of phase space flux alpha
  double charge, mass;
  const struct gk_geometry *gk_geom; // Pointer to geometry struct.
  const struct gkyl_dg_geom *dg_geom; // Pointer to vol dg geometry struct.
  const struct gkyl_gk_dg_geom *gk_dg_geom; // Pointer to vol gk dg geometry struct.
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.

  uint32_t flags;
  struct gkyl_gk_collisionless_flux *on_dev; // pointer to itself or device data
};

//
// Serendipity surface kernels.
//
// Gyrokinetic phase space flux alpha surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_surfx_kernels[] = {
  { NULL, gk_collisionless_flux_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_surfx_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_surfx_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_edge_surfx_kernels[] = {
  { NULL, gk_collisionless_flux_edge_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_edge_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_edge_surfx_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_edge_surfx_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_surfy_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux edge surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_edge_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_edge_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_edge_surfy_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_surfz_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux edge surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_edge_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_edge_surfz_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in vpar (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surfvpar_kern_list ser_gk_collisionless_flux_surfvpar_kernels[] = {
  { NULL, gk_collisionless_flux_surfvpar_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_surfvpar_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_surfvpar_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_surfvpar_3x2v_ser_p1, NULL }, // 3
};

//
// Kernels used at multiblock boundaries.
//

// Gyrokinetic phase space flux alpha surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_multib_boundary_surfx_kernels[] = {
  { NULL, gk_collisionless_flux_multib_boundary_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_multib_boundary_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_multib_boundary_surfx_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_multib_boundary_surfx_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_multib_boundary_edge_surfx_kernels[] = {
  { NULL, gk_collisionless_flux_multib_boundary_edge_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_multib_boundary_edge_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_multib_boundary_edge_surfx_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_multib_boundary_edge_surfx_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_multib_boundary_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_multib_boundary_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_multib_boundary_surfy_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux edge surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_multib_boundary_edge_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_multib_boundary_edge_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_multib_boundary_edge_surfy_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_multib_boundary_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_multib_boundary_surfz_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux edge surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_multib_boundary_edge_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_multib_boundary_edge_surfz_3x2v_ser_p1, NULL }, // 3
};


//
// Serendipity surface kernels without toroidal field (by=0)
//
// Gyrokinetic phase space flux alpha surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_surfx_kernels[] = {
  { NULL, gk_collisionless_flux_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_no_by_surfx_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_no_by_surfx_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_edge_surfx_kernels[] = {
  { NULL, gk_collisionless_flux_edge_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_edge_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_no_by_edge_surfx_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_no_by_edge_surfx_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_no_by_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_no_by_surfy_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_edge_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_no_by_edge_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_no_by_edge_surfy_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_no_by_surfz_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_edge_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_no_by_edge_surfz_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in vpar (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surfvpar_kern_list ser_gk_collisionless_flux_no_by_surfvpar_kernels[] = {
  { NULL, gk_collisionless_flux_surfvpar_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_surfvpar_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_no_by_surfvpar_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_no_by_surfvpar_3x2v_ser_p1, NULL }, // 3
};

//
// Kernels used at multiblock boundaries.
//

// Gyrokinetic phase space flux alpha surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_multib_boundary_surfx_kernels[] = {
  { NULL, gk_collisionless_flux_multib_boundary_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_multib_boundary_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_no_by_multib_boundary_surfx_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_no_by_multib_boundary_surfx_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_multib_boundary_edge_surfx_kernels[] = {
  { NULL, gk_collisionless_flux_multib_boundary_edge_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_multib_boundary_edge_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_no_by_multib_boundary_edge_surfx_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_no_by_multib_boundary_edge_surfx_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_multib_boundary_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_no_by_multib_boundary_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_no_by_multib_boundary_surfy_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_multib_boundary_edge_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_no_by_multib_boundary_edge_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_no_by_multib_boundary_edge_surfy_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_multib_boundary_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_no_by_multib_boundary_surfz_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_multib_boundary_edge_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_no_by_multib_boundary_edge_surfz_3x2v_ser_p1, NULL }, // 3
};


GKYL_CU_D
static gk_collisionless_flux_surf_t
choose_gk_collisionless_flux_surf_conf_kern(int dir, int cdim, int vdim, int poly_order, enum gkyl_gyrokinetic_bc_type bc)
{
  if (bc == GKYL_BC_GK_SKIP) {
    if (dir == 0)
      return ser_gk_collisionless_flux_multib_boundary_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return ser_gk_collisionless_flux_multib_boundary_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return ser_gk_collisionless_flux_multib_boundary_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else
      return NULL;
  }
  else {
    if (dir == 0)
      return ser_gk_collisionless_flux_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return ser_gk_collisionless_flux_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return ser_gk_collisionless_flux_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else
      return NULL;
  }
}

GKYL_CU_D
static gk_collisionless_flux_surf_t
choose_gk_collisionless_flux_edge_surf_conf_kern(int dir, int cdim, int vdim, int poly_order, enum gkyl_gyrokinetic_bc_type bc)
{
  if (bc == GKYL_BC_GK_SKIP) {
    if (dir == 0)
      return ser_gk_collisionless_flux_multib_boundary_edge_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return ser_gk_collisionless_flux_multib_boundary_edge_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return ser_gk_collisionless_flux_multib_boundary_edge_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else
      return NULL;
  }
  else {
    if (dir == 0)
      return ser_gk_collisionless_flux_edge_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return ser_gk_collisionless_flux_edge_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return ser_gk_collisionless_flux_edge_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else
      return NULL;
  }
}

GKYL_CU_D
static gk_collisionless_flux_surfvpar_t
choose_gk_collisionless_flux_surf_vpar_kern(int cdim, int vdim, int poly_order)
{
  return ser_gk_collisionless_flux_surfvpar_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}

GKYL_CU_D
static gk_collisionless_flux_surf_t
choose_gk_collisionless_flux_no_by_surf_conf_kern(int dir, int cdim, int vdim, int poly_order, enum gkyl_gyrokinetic_bc_type bc)
{
  if (bc == GKYL_BC_GK_SKIP) {
    if (dir == 0)
      return ser_gk_collisionless_flux_no_by_multib_boundary_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)                                       
      return ser_gk_collisionless_flux_no_by_multib_boundary_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)                                       
      return ser_gk_collisionless_flux_no_by_multib_boundary_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else
      return NULL;
  }
  else {
    if (dir == 0)
      return ser_gk_collisionless_flux_no_by_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return ser_gk_collisionless_flux_no_by_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return ser_gk_collisionless_flux_no_by_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else
      return NULL;
  }
}

GKYL_CU_D
static gk_collisionless_flux_surf_t
choose_gk_collisionless_flux_no_by_edge_surf_conf_kern(int dir, int cdim, int vdim, int poly_order, enum gkyl_gyrokinetic_bc_type bc)
{
  if (bc == GKYL_BC_GK_SKIP) {
    if (dir == 0)
      return ser_gk_collisionless_flux_no_by_multib_boundary_edge_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return ser_gk_collisionless_flux_no_by_multib_boundary_edge_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return ser_gk_collisionless_flux_no_by_multib_boundary_edge_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else
      return NULL;
  }
  else {
    if (dir == 0)
      return ser_gk_collisionless_flux_no_by_edge_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return ser_gk_collisionless_flux_no_by_edge_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return ser_gk_collisionless_flux_no_by_edge_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else
      return NULL;
  }
}

GKYL_CU_D
static gk_collisionless_flux_surfvpar_t
choose_gk_collisionless_flux_no_by_surf_vpar_kern(int cdim, int vdim, int poly_order)
{
  return ser_gk_collisionless_flux_no_by_surfvpar_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}

#ifdef GKYL_HAVE_CUDA
/**
 * Create new updater to compute gyrokinetic variables on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_gk_collisionless_flux* 
gkyl_gk_collisionless_flux_cu_dev_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  double charge, double mass,
  enum gkyl_gk_collisionless_type collless_type,
  const struct gk_geometry *gk_geom, const struct gkyl_dg_geom *dg_geom, 
  const struct gkyl_gk_dg_geom *gk_dg_geom, const struct gkyl_velocity_map *vel_map,
  const enum gkyl_gyrokinetic_bc_type *bctype_conf);

/**
 * Host-side wrappers for gyrokinetic vars operations on device
 */
void gkyl_gk_collisionless_flux_surf_cu(struct gkyl_gk_collisionless_flux *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_ext_range, const struct gkyl_range *phase_ext_range, const struct gkyl_array *phi, 
  const struct gkyl_array* fin, struct gkyl_array* flux_surf, struct gkyl_array* cflrate);
#endif
