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
  const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
  const double *phiL, const double *phiR, const double *JfL, const double *JfR,
  double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, double* GKYL_RESTRICT flux_surf); 

typedef double (*gk_collisionless_flux_surfvpar_t)( 
  const double *w, const double *dxv, 
  const double *vmap_prime_l, const double *vmap_prime_r,
  const double *vmap, const double *vmapSq, const double q_, const double m_, 
  const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
  const double *bmag, const double *yfield, const double *JfL, const double *JfR, double* GKYL_RESTRICT flux_surf); 

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

// Gyrokinetic phase space flux flux surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_surfy_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_surfz_3x2v_ser_p1, NULL }, // 3
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
// Boundary kernels: lower/upper domain or block edge, split by periodic,
// nonperiodic (physical BC), and multib (multiblock/SKIP) since the
// treatment of the non-owning side's data differs in each case.
//

// x direction, lower/upper edge.
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfx_lower_periodic_kernels[] = {
  { NULL, gk_collisionless_flux_boundary_surfx_lower_periodic_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_boundary_surfx_lower_periodic_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfx_lower_periodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfx_lower_periodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfx_upper_periodic_kernels[] = {
  { NULL, gk_collisionless_flux_boundary_surfx_upper_periodic_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_boundary_surfx_upper_periodic_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfx_upper_periodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfx_upper_periodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfx_lower_nonperiodic_kernels[] = {
  { NULL, gk_collisionless_flux_boundary_surfx_lower_nonperiodic_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_boundary_surfx_lower_nonperiodic_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfx_lower_nonperiodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfx_lower_nonperiodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfx_upper_nonperiodic_kernels[] = {
  { NULL, gk_collisionless_flux_boundary_surfx_upper_nonperiodic_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_boundary_surfx_upper_nonperiodic_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfx_upper_nonperiodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfx_upper_nonperiodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfx_lower_multib_kernels[] = {
  { NULL, gk_collisionless_flux_boundary_surfx_lower_multib_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_boundary_surfx_lower_multib_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfx_lower_multib_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfx_lower_multib_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfx_upper_multib_kernels[] = {
  { NULL, gk_collisionless_flux_boundary_surfx_upper_multib_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_boundary_surfx_upper_multib_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfx_upper_multib_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfx_upper_multib_3x2v_ser_p1, NULL }, // 3
};

// y direction, lower/upper edge.
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfy_lower_periodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfy_lower_periodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfy_lower_periodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfy_upper_periodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfy_upper_periodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfy_upper_periodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfy_lower_nonperiodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfy_lower_nonperiodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfy_lower_nonperiodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfy_upper_nonperiodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfy_upper_nonperiodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfy_upper_nonperiodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfy_lower_multib_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfy_lower_multib_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfy_lower_multib_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfy_upper_multib_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfy_upper_multib_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfy_upper_multib_3x2v_ser_p1, NULL }, // 3
};

// z direction, lower/upper edge.
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfz_lower_periodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfz_lower_periodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfz_upper_periodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfz_upper_periodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfz_lower_nonperiodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfz_lower_nonperiodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfz_upper_nonperiodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfz_upper_nonperiodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfz_lower_multib_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfz_lower_multib_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_boundary_surfz_upper_multib_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfz_upper_multib_3x2v_ser_p1, NULL }, // 3
};

//
// Serendipity surface kernels without toroidal field (by=0)
//
// Gyrokinetic phase space flux alpha surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_surfx_kernels[] = {
  { NULL, gk_collisionless_flux_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_surfx_no_by_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_surfx_no_by_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_surfy_no_by_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_surfy_no_by_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_surfz_no_by_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in vpar (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surfvpar_kern_list ser_gk_collisionless_flux_no_by_surfvpar_kernels[] = {
  { NULL, gk_collisionless_flux_surfvpar_1x1v_ser_p1, NULL }, // 0
  { NULL, gk_collisionless_flux_surfvpar_1x2v_ser_p1, NULL }, // 1
  { NULL, gk_collisionless_flux_surfvpar_no_by_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_surfvpar_no_by_3x2v_ser_p1, NULL }, // 3
};

//
// no_by boundary kernels: lower/upper domain or block edge, split by
// periodic, nonperiodic, and multib. Not yet generated (2x2v/3x2v to
// follow); all entries NULL for now.
//
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfx_lower_periodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfx_no_by_lower_periodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfx_no_by_lower_periodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfx_upper_periodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfx_no_by_upper_periodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfx_no_by_upper_periodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfx_lower_nonperiodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfx_no_by_lower_nonperiodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfx_no_by_lower_nonperiodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfx_upper_nonperiodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfx_no_by_upper_nonperiodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfx_no_by_upper_nonperiodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfx_lower_multib_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfx_no_by_lower_multib_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfx_no_by_lower_multib_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfx_upper_multib_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfx_no_by_upper_multib_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfx_no_by_upper_multib_3x2v_ser_p1, NULL }, // 3
};

GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfy_lower_periodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfy_no_by_lower_periodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfy_no_by_lower_periodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfy_upper_periodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfy_no_by_upper_periodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfy_no_by_upper_periodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfy_lower_nonperiodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfy_no_by_lower_nonperiodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfy_no_by_lower_nonperiodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfy_upper_nonperiodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfy_no_by_upper_nonperiodic_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfy_no_by_upper_nonperiodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfy_lower_multib_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfy_no_by_lower_multib_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfy_no_by_lower_multib_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfy_upper_multib_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gk_collisionless_flux_boundary_surfy_no_by_upper_multib_2x2v_ser_p1, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfy_no_by_upper_multib_3x2v_ser_p1, NULL }, // 3
};

GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfz_lower_periodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfz_no_by_lower_periodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfz_upper_periodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfz_no_by_upper_periodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfz_lower_nonperiodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfz_no_by_lower_nonperiodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfz_upper_nonperiodic_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfz_no_by_upper_nonperiodic_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfz_lower_multib_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfz_no_by_lower_multib_3x2v_ser_p1, NULL }, // 3
};
GKYL_CU_D
static const gkyl_gk_collisionless_flux_surf_kern_list ser_gk_collisionless_flux_no_by_boundary_surfz_upper_multib_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gk_collisionless_flux_boundary_surfz_no_by_upper_multib_3x2v_ser_p1, NULL }, // 3
};

GKYL_CU_D
static gk_collisionless_flux_surf_t
choose_gk_collisionless_flux_surf_conf_kern(int dir, int cdim, int vdim, int poly_order, enum gkyl_gyrokinetic_bc_type bc)
{
  if (dir == 0)
    return ser_gk_collisionless_flux_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 1)
    return ser_gk_collisionless_flux_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 2)
    return ser_gk_collisionless_flux_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else
    return NULL;
}

// Chooses the boundary (lower- or upper-edge) kernel for direction 'dir',
// discriminating between periodic, nonperiodic (physical BC), and
// multiblock (SKIP) treatment of the non-owning side's data.
GKYL_CU_D
static gk_collisionless_flux_surf_t
choose_gk_collisionless_flux_boundary_conf_kern(int dir, int cdim, int vdim, int poly_order,
  enum gkyl_gyrokinetic_bc_type bc, enum gkyl_edge_loc edge)
{
  bool is_lower = edge == GKYL_LOWER_EDGE;
  if (bc == GKYL_BC_GK_SKIP) {
    if (dir == 0)
      return (is_lower? ser_gk_collisionless_flux_boundary_surfx_lower_multib_kernels
                      : ser_gk_collisionless_flux_boundary_surfx_upper_multib_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return (is_lower? ser_gk_collisionless_flux_boundary_surfy_lower_multib_kernels
                      : ser_gk_collisionless_flux_boundary_surfy_upper_multib_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return (is_lower? ser_gk_collisionless_flux_boundary_surfz_lower_multib_kernels
                      : ser_gk_collisionless_flux_boundary_surfz_upper_multib_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else
      return NULL;
  }
  else if (bc == GKYL_BC_GK_SPECIES_PERIODIC) {
    if (dir == 0)
      return (is_lower? ser_gk_collisionless_flux_boundary_surfx_lower_periodic_kernels
                      : ser_gk_collisionless_flux_boundary_surfx_upper_periodic_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return (is_lower? ser_gk_collisionless_flux_boundary_surfy_lower_periodic_kernels
                      : ser_gk_collisionless_flux_boundary_surfy_upper_periodic_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return (is_lower? ser_gk_collisionless_flux_boundary_surfz_lower_periodic_kernels
                      : ser_gk_collisionless_flux_boundary_surfz_upper_periodic_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else
      return NULL;
  }
  else {
    if (dir == 0)
      return (is_lower? ser_gk_collisionless_flux_boundary_surfx_lower_nonperiodic_kernels
                      : ser_gk_collisionless_flux_boundary_surfx_upper_nonperiodic_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return (is_lower? ser_gk_collisionless_flux_boundary_surfy_lower_nonperiodic_kernels
                      : ser_gk_collisionless_flux_boundary_surfy_upper_nonperiodic_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return (is_lower? ser_gk_collisionless_flux_boundary_surfz_lower_nonperiodic_kernels
                      : ser_gk_collisionless_flux_boundary_surfz_upper_nonperiodic_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
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
  if (dir == 0)
    return ser_gk_collisionless_flux_no_by_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 1)
    return ser_gk_collisionless_flux_no_by_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 2)
    return ser_gk_collisionless_flux_no_by_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else
    return NULL;
}

// no_by counterpart of choose_gk_collisionless_flux_boundary_conf_kern().
GKYL_CU_D
static gk_collisionless_flux_surf_t
choose_gk_collisionless_flux_boundary_conf_kern_no_by(int dir, int cdim, int vdim, int poly_order,
  enum gkyl_gyrokinetic_bc_type bc, enum gkyl_edge_loc edge)
{
  bool is_lower = edge == GKYL_LOWER_EDGE;
  if (bc == GKYL_BC_GK_SKIP) {
    if (dir == 0)
      return (is_lower? ser_gk_collisionless_flux_no_by_boundary_surfx_lower_multib_kernels
                       : ser_gk_collisionless_flux_no_by_boundary_surfx_upper_multib_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return (is_lower? ser_gk_collisionless_flux_no_by_boundary_surfy_lower_multib_kernels
                       : ser_gk_collisionless_flux_no_by_boundary_surfy_upper_multib_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return (is_lower? ser_gk_collisionless_flux_no_by_boundary_surfz_lower_multib_kernels
                       : ser_gk_collisionless_flux_no_by_boundary_surfz_upper_multib_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else
      return NULL;
  }
  else if (bc == GKYL_BC_GK_SPECIES_PERIODIC) {
    if (dir == 0)
      return (is_lower? ser_gk_collisionless_flux_no_by_boundary_surfx_lower_periodic_kernels
                       : ser_gk_collisionless_flux_no_by_boundary_surfx_upper_periodic_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return (is_lower? ser_gk_collisionless_flux_no_by_boundary_surfy_lower_periodic_kernels
                       : ser_gk_collisionless_flux_no_by_boundary_surfy_upper_periodic_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return (is_lower? ser_gk_collisionless_flux_no_by_boundary_surfz_lower_periodic_kernels
                       : ser_gk_collisionless_flux_no_by_boundary_surfz_upper_periodic_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else
      return NULL;
  }
  else {
    if (dir == 0)
      return (is_lower? ser_gk_collisionless_flux_no_by_boundary_surfx_lower_nonperiodic_kernels
                       : ser_gk_collisionless_flux_no_by_boundary_surfx_upper_nonperiodic_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 1)
      return (is_lower? ser_gk_collisionless_flux_no_by_boundary_surfy_lower_nonperiodic_kernels
                       : ser_gk_collisionless_flux_no_by_boundary_surfy_upper_nonperiodic_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    else if (dir == 2)
      return (is_lower? ser_gk_collisionless_flux_no_by_boundary_surfz_lower_nonperiodic_kernels
                       : ser_gk_collisionless_flux_no_by_boundary_surfz_upper_nonperiodic_kernels)[cv_index[cdim].vdim[vdim]].kernels[poly_order];
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
