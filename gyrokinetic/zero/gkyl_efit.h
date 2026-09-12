#pragma once

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_efit gkyl_efit;

struct gkyl_efit_inp {
  // Inputs to get psiRZ and related inputs from efit
  char filepath[1024]; // filepath path to eqdsk file
  int rz_poly_order; // poly order for DG rep of psi, psi/R, and psi/R^2
  int flux_poly_order; // poly order to use for DG rep of F(psi)
  bool reflect; // whether to reflect across R axis to preserve symmetry
  bool use_gpu; // whether to use the GPU

  // Optional lower-half vessel outline (Z <= 0) used to restrict the X-point search.
  // psi has real saddle points outside the machine (coil-induced); the search picks
  // the critical point closest in flux to sibry, which without a spatial constraint
  // can select one of those instead of the plasma X-point.  When xpt_bound_n > 0,
  // candidates outside this polygon are rejected.  Leave unset (0/NULL) to search
  // the whole domain, which is the previous behaviour.
  int xpt_bound_n;             // number of polygon vertices (0 disables the filter)
  const double *xpt_bound_R;   // vertex R coordinates
  const double *xpt_bound_Z;   // vertex Z coordinates (all <= 0)
};

struct gkyl_efit{
  const char* filepath;
  char name[1024]; // name (filepath stripped of full path and extension)
  int nr, nz;
  double rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry, bcentr, current, xdum;
  double rmin, rmax, zmin, zmax;

  // Actual EQDSK limiter polyline, in file order (not xpt_bound_*).
  // The outline closes last-to-first; reflection does not alter these vertices.
  // Owned by this object.
  // 1: usable -- >=3 finite vertices, enough to bound a region.
  // 0: absent -- the file supplies no limiter record, or declares zero vertices.
  // 2: degenerate -- a readable record with 1 or 2 vertices, which cannot bound
  //    a region. Benign input, not corrupt data; separated from -1 so that an
  //    explicit no-vessel-outline declaration may cover it while corrupt data
  //    stays un-declarable.
  // -1: malformed -- counts or coordinates could not be read, were out of range,
  //    or were non-finite. Says nothing about polygon SHAPE: convexity,
  //    self-intersection and winding are not examined here or anywhere else.
  int limiter_status;
  int limiter_n;
  double *limiter_R, *limiter_Z;

  // Outline SHAPE, examined separately from the vertex-count classification
  // above. Reported, never enforced: a self-intersecting outline is still
  // limiter_status = 1, because the one real instance we have measured
  // (tcv_upper_SN.geqdsk) is a degenerate pinch that encloses no area and
  // changes no containment answer. Refusing on that evidence would reject
  // usable data; saying nothing would let a genuinely overlapping outline be
  // handed an arbitrary interior by the crossing-number test.
  int limiter_self_intersections; // 0 for a well-formed outline

  double psisep; // Separatrix psi for our DG representation
  double psisep_cubic; // Separatrix psi for our cubic DG representation
                 // Can differ from sibry, but we need to keep sibry
                 // because fpol, q, etc. are defined based on it

  struct gkyl_basis rzbasis;
  struct gkyl_basis rzbasis_cubic;
  struct gkyl_rect_grid rzgrid;
  struct gkyl_range rzlocal;
  struct gkyl_range rzlocal_ext;
  struct gkyl_rect_grid rzgrid_cubic;
  struct gkyl_range rzlocal_cubic;
  struct gkyl_range rzlocal_cubic_ext;
  struct gkyl_array *psizr;
  struct gkyl_array *psizr_cubic;
  struct gkyl_array *bmagzr;
  struct gkyl_basis_ops_evalf *evf ;

  struct gkyl_basis fluxbasis;
  struct gkyl_rect_grid fluxgrid;
  struct gkyl_range fluxlocal;
  struct gkyl_range fluxlocal_ext;
  struct gkyl_array* fpolflux;
  struct gkyl_array* fpolprimeflux;
  struct gkyl_array* qflux;

  int num_xpts; // Number of X-points
  double *Rxpt; // R coordinates of X points
  double *Zxpt; // Z coordinates of X-points

  int num_xpts_cubic; // Number of X-points of cubic rep
  double *Rxpt_cubic; // R coordinates of X points of cubic rep
  double *Zxpt_cubic; // Z coordinates of X-points of cubic rep

  bool xpt_diag_quadratic_found;
  bool xpt_diag_fallback_to_cubic;
  double xpt_diag_quad_R;
  double xpt_diag_quad_Z;
  double xpt_diag_quad_psi;
  double xpt_diag_quad_dist_cell;

  bool reflect;
  bool use_gpu;

  // Copy of the optional X-point search bounds (see struct gkyl_efit_inp).
  int xpt_bound_n;
  double *xpt_bound_R;
  double *xpt_bound_Z;
};

/**
 * Create new updater which reads in magnetic equilibrium 
 * parameters from a geqdsk file,
 * projects the poloidal flux psi(R,Z), psi/R, psi/R^2 on 
 * the RZ grid,
 * and projects F(psi)= R*B_phi on a poloidal flux grid.
 * Also finds X-points of the EQDSK file.
 *
 * @param inp efit_inp input to create efit updayer
 * @return New updater pointer.
 */

gkyl_efit* gkyl_efit_new(const struct gkyl_efit_inp *inp);

/**
 * Number of properly crossing non-adjacent edge pairs in the vessel outline.
 *
 * Zero for a well-formed outline. A nonzero count means the outline overlaps
 * itself, so the crossing-number interior test is answering a question with no
 * single right answer. Reported, never enforced -- see limiter_self_intersections.
 *
 * @param e EFIT object with a usable outline (limiter_status == 1).
 * @return Number of crossing edge pairs; 0 if there is no usable outline.
 */
int gkyl_efit_limiter_self_intersections(const struct gkyl_efit *e);

/**
 * Fetch magnetic-axis and separatrix psi from an EFIT object.
 *
 * @param up EFIT object
 * @param simag On output, magnetic-axis psi
 * @param psisep On output, separatrix psi for the DG representation
 */
void gkyl_efit_get_psi_bounds(const gkyl_efit *up, double *simag, double *psisep);

void gkyl_efit_release(gkyl_efit* up);
