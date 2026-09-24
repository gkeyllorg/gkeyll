#pragma once

#include <gkyl_array.h>
#include <gkyl_evalf_def.h>
#include <gkyl_math.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Geometry information for lower surfaces (left, bottom, front) for a cell
struct gkyl_dg_surf_geom {
  double area_elem; // area element: Jc*norm(e^d), where e^d is dual normal
  struct gkyl_vec3 norm; // norm is the normal to face perp to direction 'd'
  // tau1 X tau2 = norm are tangents to face perp to direction 'd'
  struct gkyl_vec3 tau1, tau2;
};

// Geometry information for a single cell
struct gkyl_dg_vol_geom {
  struct gkyl_vec3 tang[3]; // tangent vectors, e_i
  struct gkyl_vec3 dual[3]; // dual vectors, e^i
  double Jc; // Jacobian = e_1*(e_2 X e_3)  = 1/e^1*(e^2 X e^3)
};

// geometry information over a range of cells: 
struct gkyl_dg_geom {
  struct gkyl_range range; // range over which geometry is defined
  struct gkyl_range surf_quad_range; // range for indexing surface nodes
  struct gkyl_range vol_quad_range;  // range for indexing volume nodes

  // index weights and ords arrays below using the appropriate
  // methods below
  
  // weights and ordinates for surface quadrature
  double *surf_weights, *surf_ords;
  // weights and ordinates for surface quadrature
  double *vol_weights, *vol_ords;
  
  struct gkyl_array *surf_geom[GKYL_MAX_CDIM]; // surface geometry in dir 'd' in each cell
  struct gkyl_array *vol_geom; // cell geometry
  
  uint32_t flags;
  struct gkyl_ref_count ref_count;  
  struct gkyl_dg_geom *on_dev; // pointer to itself or device object
};

// Input for constructing geometry with mapc2p
struct gkyl_dg_geom_inp {
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
struct gkyl_dg_geom* gkyl_dg_geom_new(const struct gkyl_dg_geom_inp *inp);

/**
 * Create a new DG geometry object from host object. 
 *
 * @param up_host geometry object on host
 * @param inp Inputs for use in constructing geometry
 * @param use_gpu whether to use gpu
 */
struct gkyl_dg_geom* gkyl_dg_geom_new_from_host(const struct gkyl_dg_geom_inp *inp, struct gkyl_dg_geom *up_host, bool use_gpu);

struct gkyl_dg_geom* gkyl_dg_geom_cu_dev_new_from_host(const struct gkyl_dg_geom_inp *inp, struct gkyl_dg_geom *up_host);


/**
 * Acquire pointer to geometry object. The pointer must be released
 * using gkyl_dg_geom_release method.
 *
 * @param dgg Geometry to which a pointer is needed
 * @return Pointer to acquired geometry
 */
struct gkyl_dg_geom* gkyl_dg_geom_acquire(const struct gkyl_dg_geom* dgg);

/**
 * Write out geometry data to file. The "fprefix" is the prefix of the
 * set of file names.
 *
 * @param dgg Geometry to which a pointer is needed
 * @param fname Name of output file to write
 */
void gkyl_dg_geom_write(const struct gkyl_dg_geom* dgg, const char *fprefix);

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
static inline const struct gkyl_dg_surf_geom*
gkyl_dg_geom_get_surf(const struct gkyl_dg_geom *dgg, int d, const int *idx)
{
  return (const struct gkyl_dg_surf_geom*) gkyl_array_cfetch(dgg->surf_geom[d], gkyl_range_idx(&dgg->range, idx));
}

/**
 * Get linear index into surface quadrature node corresponding to
 * multi-dimensional index @a sidx.
 *
 * @param dgg DG geometry object
 * @param sidx Index (ndim-1) of surface quadrature node
 * @return Weight at surface quadrature node
 */
GKYL_CU_DH
static inline long
gkyl_dg_geom_surf_quad_idx(const struct gkyl_dg_geom *dgg, const int *sidx)
{
  return gkyl_range_idx(&dgg->surf_quad_range, sidx);
}

/**
 * Get weight for surface quadrature node corresponding to
 * multi-dimensional surface index @a sidx.
 *
 * @param dgg DG geometry object
 * @param sidx Index (ndim-1) of surface quadrature node
 * @return Weight at node
 */
GKYL_CU_DH
static inline double
gkyl_dg_geom_surf_quad_weight(const struct gkyl_dg_geom *dgg, const int *sidx)
{
  return dgg->surf_weights[gkyl_dg_geom_surf_quad_idx(dgg, sidx)];
}

/**
 * Get ordinates for surface quadrature node corresponding to
 * multi-dimensional surface index @a sidx.
 *
 * @param dgg DG geometry object
 * @param sidx Index of surface quadrature node
 * @return Ordinates at nodes: ndim-1 size array
 */
GKYL_CU_DH
static inline const double*
gkyl_dg_geom_surf_quad_ords(const struct gkyl_dg_geom *dgg, const int *sidx)
{
  int sdim = dgg->surf_quad_range.ndim;
  return &dgg->surf_ords[sdim*gkyl_dg_geom_surf_quad_idx(dgg, sidx)];
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
static inline const struct gkyl_dg_vol_geom*
gkyl_dg_geom_get_vol(const struct gkyl_dg_geom *dgg, const int *idx)
{
  return (const struct gkyl_dg_vol_geom*) gkyl_array_cfetch(dgg->vol_geom, gkyl_range_idx(&dgg->range, idx));
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
gkyl_dg_geom_vol_quad_idx(const struct gkyl_dg_geom *dgg, const int *vidx)
{
  return gkyl_range_idx(&dgg->vol_quad_range, vidx);
}

/**
 * Get weight for volume quadrature node corresponding to
 * multi-dimensional vol index @a vidx.
 *
 * @param dgg DG geometry object
 * @param vidx Index of volume quadrature node
 * @return Weight at node
 */
GKYL_CU_DH
static inline double
gkyl_dg_geom_vol_quad_weight(const struct gkyl_dg_geom *dgg, const int *vidx)
{
  return dgg->vol_weights[gkyl_dg_geom_vol_quad_idx(dgg, vidx)];
}

/**
 * Get ordinates for volume quadrature node corresponding to
 * multi-dimensional volume index @a sidx.
 *
 * @param dgg DG geometry object
 * @param vidx Index of volume quadrature node
 * @return Ordinates at nodes: ndim size array
 */
GKYL_CU_DH
static inline const double*
gkyl_dg_geom_vol_quad_ords(const struct gkyl_dg_geom *dgg, const int *vidx)
{
  int ndim = dgg->vol_quad_range.ndim;
  return &dgg->vol_ords[ndim*gkyl_dg_geom_vol_quad_idx(dgg, vidx)];
}

/**
 * Release geometry object.
 *
 * @param dgg Dg geometry object to release.
 */
void gkyl_dg_geom_release(const struct gkyl_dg_geom *dgg);

/**
 * Release DG geometry object from mapc2p. 
 *
 * @param inp Inputs for use in constructing geometry
 */
void dg_geom_free(const struct gkyl_ref_count *ref);

