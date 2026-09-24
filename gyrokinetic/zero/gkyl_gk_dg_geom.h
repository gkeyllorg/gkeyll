#pragma once

#include <gkyl_gk_geometry.h>
#include <gkyl_dg_geom.h>
#include <gkyl_array.h>
#include <gkyl_evalf_def.h>
#include <gkyl_math.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>


// Geometry information for lower surfaces (left, bottom, front) for a cell
struct gkyl_gk_dg_surf_geom {
  double B3; // n^3 \cdot \vec{B}
  double normcurlbhat; // normal to face perp to direction 'd' dotted with curl(bhat)
  double bmag; // |B| 
  struct gkyl_vec3 bhat; // Covariant components of bhat (b_i)
  double Jc; // J_c
};

// Geometry information for a single cell
struct gkyl_gk_dg_vol_geom {
  double bmag; // |B|
  double B3; // n^3 \cdot \vec{B}
  struct gkyl_vec3 dualcurlbhat; // duals dotted with curl(bhat)
};


// geometry information over a range of cells: 
struct gkyl_gk_dg_geom {
  struct gkyl_range range; // range over which geometry is defined
  struct gkyl_range surf_quad_range; // range for indexing surface nodes
  struct gkyl_range vol_quad_range;  // range for indexing volume nodes

  
  struct gkyl_array *surf_geom[GKYL_MAX_CDIM]; // surface geometry in dir 'd' in each cell
  struct gkyl_array *vol_geom; // cell geometry
  
  uint32_t flags;
  struct gkyl_ref_count ref_count;  
  struct gkyl_gk_dg_geom *on_dev; // pointer to itself or device object
};

// Input for constructing geometry with mapc2p
struct gkyl_gk_dg_geom_inp {
  const struct gkyl_rect_grid *grid; // grid over which geometry needs to be defined
  struct gkyl_range *range; // range of grid: this needs to be extended grid
  bool skip_geometry_creation; // true if only allocation is needed (geom must be computed eleswhere)
  int nquad; // number of quadrature points in each direction
  evalf_t mapc2p; // mapping function: set mapc2p = 0 to use identity map
  void *ctx; // context for use in mapc2p
  bool use_gpu; // does this live on the GPU?
};

/**
 * Create a new DG geometry object from mapc2p. 
 *
 * @param inp Inputs for use in constructing geometry
 */
struct gkyl_gk_dg_geom* gkyl_gk_dg_geom_new(const struct gkyl_gk_dg_geom_inp *inp);

/**
 * Create a new DG geometry object from host object. 
 *
 * @param up_host geometry object on host
 * @param inp Inputs for use in constructing geometry
 * @param use_gpu whether to use gpu
 */
struct gkyl_gk_dg_geom* gkyl_gk_dg_geom_new_from_host(const struct gkyl_gk_dg_geom_inp *inp, struct gkyl_gk_dg_geom *up_host, bool use_gpu);

struct gkyl_gk_dg_geom* gkyl_gk_dg_geom_cu_dev_new_from_host(const struct gkyl_gk_dg_geom_inp *inp, struct gkyl_gk_dg_geom *up_host);

/**
 * Acquire pointer to geometry object. The pointer must be released
 * using gkyl_gk_dg_geom_release method.
 *
 * @param dgg Geometry to which a pointer is needed
 * @return Pointer to acquired geometry
 */
struct gkyl_gk_dg_geom* gkyl_gk_dg_geom_acquire(const struct gkyl_gk_dg_geom* dgg);

/**
 * Write out geometry data to file. The "fprefix" is the prefix of the
 * set of file names.
 *
 * @param dgg Geometry to which a pointer is needed
 * @param fname Name of output file to write
 */
void gkyl_gk_dg_geom_write(const struct gkyl_gk_dg_geom* dgg, const char *fprefix);

/**
 * Get pointer to geometry on the surface normal to 'd' given by idx
 * into the range over which the geometry was constructed. The
 * returned array needs to be further indexed to fetch the specific
 * struct at a quadrature point
 *
 * @param dgg DG geometry object
 * @param d Direction to which surface is perpendicular
 * @param idx Index into grid
 * @return Pointer to surface geometry at all quadrature nodes in cell @a idx
 */
GKYL_CU_DH
static inline const struct gkyl_gk_dg_surf_geom*
gkyl_gk_dg_geom_get_surf(const struct gkyl_gk_dg_geom *dgg, int d, const int *idx)
{
  return (const struct gkyl_gk_dg_surf_geom*) gkyl_array_cfetch(dgg->surf_geom[d], gkyl_range_idx(&dgg->range, idx));
}

/**
 * Get linear index into surface quadrature node corresponding to
 * multi-dimensional index @a idx.
 *
 * @param dgg DG geometry object
 * @param idx Index (ndim-1) of surface quadrature node
 * @return Linear index for indexing surface quadrature array
 */
GKYL_CU_DH
static inline long
gkyl_gk_dg_geom_surf_quad_idx(const struct gkyl_gk_dg_geom *dgg, const int *idx)
{
  return gkyl_range_idx(&dgg->surf_quad_range, idx);
}

/**
 * Get pointer to geometry in cell given by idx into the range over
 * which the geometry was constructed. The returned array needs to be
 * further indexed to fetch the specific struct at a quadrature point
 *
 * @param dgg DG geometry object
 * @param idx Index into grid
 * @return Pointer to cell geometry at all quadrature nodes in cell @a idx
 */
GKYL_CU_DH
static inline const struct gkyl_gk_dg_vol_geom*
gkyl_gk_dg_geom_get_vol(const struct gkyl_gk_dg_geom *dgg, const int *idx)
{
  return (const struct gkyl_gk_dg_vol_geom*) gkyl_array_cfetch(dgg->vol_geom, gkyl_range_idx(&dgg->range, idx));
}

/**
 * Get linear index into volume quadrature node corresponding to
 * multi-dimensional index @a vidx.
 *
 * @param dgg DG geometry object
 * @param vidx Index of volume quadrature node
 * @return Linear index for indexing volume quadrature array
 */
GKYL_CU_DH
static inline long
gkyl_gk_dg_geom_vol_quad_idx(const struct gkyl_gk_dg_geom *dgg, const int *vidx)
{
  return gkyl_range_idx(&dgg->vol_quad_range, vidx);
}

/**
 * Release geometry object.
 *
 * @param dgg Dg geometry object to release.
 */
void gkyl_gk_dg_geom_release(const struct gkyl_gk_dg_geom *dgg);

/**
 * Release DG geometry object from mapc2p. 
 *
 * @param inp Inputs for use in constructing geometry
 */
void gk_dg_geom_free(const struct gkyl_ref_count *ref);

/**
 * Populate volume quad point DG geometry objects 
 * using geometry from gk_geom . 
 *
 * @param dg_geom dg_geom object (stores geo common to all simulation types)
 * @param gk_dg_geom gk_dg_geom object (stores geo specific to gyrokinetics)
 * @param gk_geom gk_geom object from which to populate dg geom objects
 */
void gkyl_gk_dg_geom_populate_vol(struct gkyl_dg_geom *dg_geom, struct gkyl_gk_dg_geom *gk_dg_geom, struct gk_geometry* gk_geom);

/**
 * Populate surface quad point DG geometry objects 
 * using geometry from gk_geom . 
 *
 * @param dg_geom dg_geom object (stores geo common to all simulation types)
 * @param gk_dg_geom gk_dg_geom object (stores geo specific to gyrokinetics)
 * @param gk_geom gk_geom object from which to populate dg geom objects
 */
void gkyl_gk_dg_geom_populate_surf(struct gkyl_dg_geom *dg_geom, struct gkyl_gk_dg_geom *gk_dg_geom, struct gk_geometry* gk_geom);
