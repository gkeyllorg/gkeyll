// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_collisionless_passive_flux_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef double (*gk_collisionless_passive_flux_surf_t)(const double *w, const double *dxv,
  const double *vmap, const double *vmapSq, const double q_, const double m_,
  const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs,
  const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
  const double *speeds_L, const double *speeds_R,
  const double *JfL, const double *JfR, double *GKYL_RESTRICT flux_surf);

// The cv_index[cd].vdim[vd] is used to index the various lists of kernels below.
GKYL_CU_D
static struct { int vdim[3]; } cv_index_passive[] = {
  {-1, -1, -1}, // 0x makes no sense.
  {-1,  0,  1}, // 1x kernel indices.
  {-1, -1,  2}, // 2x kernel indices.
  {-1, -1,  3}, // 3x kernel indices.
};

// For use in kernel tables.
typedef struct { gk_collisionless_passive_flux_surf_t kernels[3]; } gkyl_gk_collisionless_passive_flux_surf_kern_list;

struct gkyl_gk_collisionless_passive_flux {
  struct gkyl_rect_grid phase_grid; // Phase space grid for cell spacing and cell center.
  int cdim; // Configuration space dimensionality.
  int pdim; // Phase space dimensionality.
  double charge, mass;
  // Pointers to projected conf-space passive-speed arrays (NOT owned by this struct).
  const struct gkyl_array *passive_speeds;
  const struct gk_geometry *gk_geom; // Pointer to geometry struct.
  const struct gkyl_dg_geom *dg_geom; // Pointer to vol DG geometry struct.
  const struct gkyl_gk_dg_geom *gk_dg_geom; // Pointer to GK DG geometry struct.
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.
  // Kernels for interior, lower-edge, and upper-edge surfaces in each conf direction.
  gk_collisionless_passive_flux_surf_t flux_surf[GKYL_MAX_CDIM];
  gk_collisionless_passive_flux_surf_t flux_surf_edge_lo[GKYL_MAX_CDIM];
  gk_collisionless_passive_flux_surf_t flux_surf_edge_up[GKYL_MAX_CDIM];
  // No vpar kernel — passive advection is conf-space only.
  uint32_t flags;
  struct gkyl_gk_collisionless_passive_flux *on_dev; // pointer to itself or device data
};

//
// Serendipity surface kernels — x-direction (interior)
//
GKYL_CU_D
static const gkyl_gk_collisionless_passive_flux_surf_kern_list ser_gk_collisionless_passive_flux_surfx_kernels[] = {
  { NULL, gk_collisionless_passive_flux_surfx_1x1v_ser_p1, NULL }, // 0: 1x1v
  { NULL, gk_collisionless_passive_flux_surfx_1x2v_ser_p1, NULL }, // 1: 1x2v
  { NULL, gk_collisionless_passive_flux_surfx_2x2v_ser_p1, NULL }, // 2: 2x2v
  { NULL, gk_collisionless_passive_flux_surfx_3x2v_ser_p1, NULL }, // 3: 3x2v
};

// x-direction edge (lower boundary)
GKYL_CU_D
static const gkyl_gk_collisionless_passive_flux_surf_kern_list ser_gk_collisionless_passive_flux_edge_surfx_kernels[] = {
  { NULL, gk_collisionless_passive_flux_edge_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_passive_flux_edge_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_passive_flux_edge_surfx_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_passive_flux_edge_surfx_3x2v_ser_p1, NULL }, // 3
};

//
// y-direction (interior)
//
GKYL_CU_D
static const gkyl_gk_collisionless_passive_flux_surf_kern_list ser_gk_collisionless_passive_flux_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0: 1x — no y
  { NULL, NULL, NULL }, // 1: 1x — no y
  { NULL, gk_collisionless_passive_flux_surfy_2x2v_ser_p1, NULL }, // 2: 2x2v
  { NULL, gk_collisionless_passive_flux_surfy_3x2v_ser_p1, NULL }, // 3: 3x2v
};

// y-direction edge
GKYL_CU_D
static const gkyl_gk_collisionless_passive_flux_surf_kern_list ser_gk_collisionless_passive_flux_edge_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_passive_flux_edge_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_passive_flux_edge_surfy_3x2v_ser_p1, NULL }, // 3
};

//
// z-direction (interior)
//
GKYL_CU_D
static const gkyl_gk_collisionless_passive_flux_surf_kern_list ser_gk_collisionless_passive_flux_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_passive_flux_surfz_3x2v_ser_p1, NULL }, // 3: 3x2v
};

// z-direction edge
GKYL_CU_D
static const gkyl_gk_collisionless_passive_flux_surf_kern_list ser_gk_collisionless_passive_flux_edge_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_passive_flux_edge_surfz_3x2v_ser_p1, NULL }, // 3
};

//
// Kernels used at multiblock boundaries.
//

// Gyrokinetic phase space flux alpha surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_passive_flux_surf_kern_list ser_gk_collisionless_passive_flux_multib_boundary_surfx_kernels[] = {
  { NULL, gk_collisionless_passive_flux_multib_boundary_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_passive_flux_multib_boundary_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_passive_flux_multib_boundary_surfx_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_passive_flux_multib_boundary_surfx_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_passive_flux_surf_kern_list ser_gk_collisionless_passive_flux_multib_boundary_edge_surfx_kernels[] = {
  { NULL, gk_collisionless_passive_flux_multib_boundary_edge_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_passive_flux_multib_boundary_edge_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_passive_flux_multib_boundary_edge_surfx_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_passive_flux_multib_boundary_edge_surfx_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_passive_flux_surf_kern_list ser_gk_collisionless_passive_flux_multib_boundary_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_passive_flux_multib_boundary_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_passive_flux_multib_boundary_surfy_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux edge surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_passive_flux_surf_kern_list ser_gk_collisionless_passive_flux_multib_boundary_edge_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_passive_flux_multib_boundary_edge_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_passive_flux_multib_boundary_edge_surfy_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_passive_flux_surf_kern_list ser_gk_collisionless_passive_flux_multib_boundary_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_passive_flux_multib_boundary_surfz_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux edge surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_passive_flux_surf_kern_list ser_gk_collisionless_passive_flux_multib_boundary_edge_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_passive_flux_multib_boundary_edge_surfz_3x2v_ser_p1, NULL }, // 3
};

// Helper: choose interior conf-surface kernel for a given direction.
GKYL_CU_D
static gk_collisionless_passive_flux_surf_t
choose_gk_collisionless_passive_flux_surf_conf_kern(int dir, int cdim, int vdim, int poly_order)
{
  if (dir == 0)
    return ser_gk_collisionless_passive_flux_surfx_kernels[cv_index_passive[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 1)
    return ser_gk_collisionless_passive_flux_surfy_kernels[cv_index_passive[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 2)
    return ser_gk_collisionless_passive_flux_surfz_kernels[cv_index_passive[cdim].vdim[vdim]].kernels[poly_order];

  return NULL;
}

// Helper: choose lower-edge conf-surface kernel for a given direction.
GKYL_CU_D
static gk_collisionless_passive_flux_surf_t
choose_gk_collisionless_passive_flux_edge_lo_surf_conf_kern(int dir, int cdim, int vdim, int poly_order, enum gkyl_gyrokinetic_bc_type bc)
{
  if (bc == GKYL_BC_GK_SKIP) {
    if (dir == 0)
      return ser_gk_collisionless_passive_flux_multib_boundary_surfx_kernels[cv_index_passive[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return ser_gk_collisionless_passive_flux_multib_boundary_surfy_kernels[cv_index_passive[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return ser_gk_collisionless_passive_flux_multib_boundary_surfz_kernels[cv_index_passive[cdim].vdim[vdim]].kernels[poly_order];
    else
      return NULL;
  }
  else {
    if (dir == 0)
      return ser_gk_collisionless_passive_flux_surfx_kernels[cv_index_passive[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return ser_gk_collisionless_passive_flux_surfy_kernels[cv_index_passive[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return ser_gk_collisionless_passive_flux_surfz_kernels[cv_index_passive[cdim].vdim[vdim]].kernels[poly_order];
    else
      return NULL;
  }

  return NULL;
}

// Helper: choose upper-edge conf-surface kernel for a given direction.
GKYL_CU_D
static gk_collisionless_passive_flux_surf_t
choose_gk_collisionless_passive_flux_edge_up_surf_conf_kern(int dir, int cdim, int vdim, int poly_order, enum gkyl_gyrokinetic_bc_type bc)
{
  if (bc == GKYL_BC_GK_SKIP) {
    if (dir == 0)
      return ser_gk_collisionless_passive_flux_multib_boundary_edge_surfx_kernels[cv_index_passive[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return ser_gk_collisionless_passive_flux_multib_boundary_edge_surfy_kernels[cv_index_passive[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return ser_gk_collisionless_passive_flux_multib_boundary_edge_surfz_kernels[cv_index_passive[cdim].vdim[vdim]].kernels[poly_order];
    else
      return NULL;
  }
  else {
    if (dir == 0)
      return ser_gk_collisionless_passive_flux_edge_surfx_kernels[cv_index_passive[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return ser_gk_collisionless_passive_flux_edge_surfy_kernels[cv_index_passive[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return ser_gk_collisionless_passive_flux_edge_surfz_kernels[cv_index_passive[cdim].vdim[vdim]].kernels[poly_order];
    return NULL;
  }

  return NULL;
}

// GPU-only functions (defined in gk_collisionless_passive_flux_cu.cu).
gkyl_gk_collisionless_passive_flux*
gkyl_gk_collisionless_passive_flux_cu_dev_new(const struct gkyl_rect_grid *phase_grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_array *passive_speeds,
  const double charge, const double mass,
  const struct gk_geometry *gk_geom, const struct gkyl_dg_geom *dg_geom,
  const struct gkyl_gk_dg_geom *gk_dg_geom, const struct gkyl_velocity_map *vel_map,
  const enum gkyl_gyrokinetic_bc_type *bctype_conf);

void gkyl_gk_collisionless_passive_flux_surf_cu(gkyl_gk_collisionless_passive_flux *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_ext_range, const struct gkyl_range *phase_ext_range,
  const struct gkyl_array *fin,
  struct gkyl_array *flux_surf, struct gkyl_array *cflrate);
