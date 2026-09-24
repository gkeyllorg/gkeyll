#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_tok_geo.h>

// Object type
typedef struct gkyl_calc_bmag gkyl_calc_bmag;

// Context object for computing bmag
struct gkyl_bmag_ctx{
   const struct gkyl_rect_grid* grid; // Physical RZ grid
   const struct gkyl_rect_grid* cgrid; // Computational grid
   const struct gkyl_range* range; // Physical RZ range
   const struct gkyl_range* crange; // Computational range
   const struct gkyl_range* crange_global; // Global computational range
   const struct gkyl_basis* basis; // Physical RZ basis
   const struct gkyl_basis* cbasis; // Computational basis 
   const struct gkyl_array* bmagdg; // DG representation of bmag in physical RZ coordinates
   const struct gkyl_array* bmag; // DG representation of bmag in computational coordinates
   const struct gkyl_array* mapc2p; // DG representation of mapc2p
};

/**
 * Create new updater to compute the bmag on the compuational grid 
 *
 * @param cbasis Basis object ( computational configuration space).
 * @ param pbasis (physical RZ basis)
 * @param cgrid computational grid
 * @param pgrid physical RZ grid
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_calc_bmag* 
gkyl_calc_bmag_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_rect_grid *cgrid, const struct gkyl_rect_grid *pgrid, bool use_gpu);

/**
 * Computes the magnitude of the magnetic field using a global field aligned representation.
 * 
 * @param t Time at which to compute the magnetic field.
 * @param xn Coordinates at which to compute the magnetic field.
 * @param fout Output array.
 * @param ctx Context object. Of type bmag_ctx
*/
void gkyl_calc_bmag_global(double t, const double *xn, double *fout, void *ctx);


/**
 * Advance calc_bmag (Convert B to computational grid from RZ grid)
 *
 * @param up calc_bmag updater object.
 * @param crange, _ext computational local Config-space range and extended range.
 * @param crange_global, _ext computational global Config-space range and extended range.
 * @param prange, _ext physical RZ range and extended range.
 * @param bmagrz DG B(R,Z) on the RZ grid
 * @param bmag_compdg output field where DG bmag on the computational basis/grid will be plaed
 * @param mapc2p DG rep of mapc2p on the computational grid
* @param bool use_quad flag to indicate whether to use quadrature points or not
 */

void gkyl_calc_bmag_advance(const gkyl_calc_bmag *up, const struct gkyl_range *crange,
    const struct gkyl_range *crange_ext,  const struct gkyl_range *crange_global,
    const struct gkyl_range *prange, const struct gkyl_range *prange_ext, 
    const struct gkyl_array *bmagrz, struct gkyl_array* bmag_compdg, struct gkyl_array* mapc2p, bool use_quad);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_calc_bmag_release(gkyl_calc_bmag* up);
