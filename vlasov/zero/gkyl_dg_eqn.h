#pragma once

#include <gkyl_ref_count.h>
#include <gkyl_util.h>

#include <stdbool.h>

// Forward declare for use in function pointers
struct gkyl_dg_eqn;

// Function pointer type for volume kernel
typedef double (*vol_termf_t)(const struct gkyl_dg_eqn *eqn, 
  const double*  xc, const double*  dx, const int*  idx,
  const double* qIn, double* GKYL_RESTRICT qRhsOut);

// Function pointer type for surface kernel
typedef double (*surf_termf_t)(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut);

// Function pointer type for boundary surface kernel
typedef double (*boundary_surf_termf_t)(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut);

// Function pointer type for generic stencil kernel
// Similar to surface kernel, but size of input arrays unspecified
// Could be 9 (for a generic 2D stencil) or 27 (for a generic 3D stencil)
// Hard code to 27 for 3D for GPU compatibility
// NOTE: ASSUMES UNIFORM GRIDS FOR NOW
// NOTE: Takes the index of the cell being updated (idxc) and array of indices
//       (idx) so we can fetch auxiliary variables easily for neighbors or just
//       the cell being updated. Need size of integer array (sz_dim)
typedef double (*gen_termf_t)(const struct gkyl_dg_eqn *eqn, 
  int dir1, int dir2,
  const double* xc, const double* dxc, const int* idxc,
  long sz_dim, const int idx[27][GKYL_MAX_DIM], const double* qIn[27], 
  double* GKYL_RESTRICT qRhsOut);

struct gkyl_dg_eqn {
  int num_equations; // Number of equations in system.
  vol_termf_t vol_term; // Volume term kernel.
  surf_termf_t surf_term; // Surface term kernel.
  boundary_surf_termf_t boundary_surf_term; // Boundary surface term kernel.
  boundary_surf_termf_t boundary_diag_term; // Boundary diagnotic kernel.
  gen_termf_t gen_surf_term; // Generic stencil kernel with input variable size unspecified.
  gen_termf_t gen_boundary_surf_term; // Generic stencil kernel with input variable size unspecified
                                      // for boundary surface updates.
  uint32_t flags;
  struct gkyl_ref_count ref_count; // Reference count.
  struct  gkyl_dg_eqn *on_dev; // Pointer to itself or device data.
};

// context for use in BCs
struct dg_bc_ctx {
  int dir; // direction for BCs
  int cdim; // config-space dimensions
  const struct gkyl_basis *basis; // basis function
};

/**
 * Check if equation is on device.
 *
 * @param eqn Equation to check
 * @return true if eqn on device, false otherwise
 */
bool gkyl_dg_eqn_is_cu_dev(const struct gkyl_dg_eqn *eqn);

/**
 * Acquire pointer to equation object. Delete using the release()
 * method
 *
 * @param eqn Equation object.
 */
struct gkyl_dg_eqn* gkyl_dg_eqn_acquire(const struct gkyl_dg_eqn* eqn);

/**
 * Compute the volume term in a DG discretization.
 *
 * @param eqn Equation object (contains the volume kernel function pointer vol_termf_t).
 * @param xc  Cell center coordinates
 * @param dx  Grid spacing
 * @param idx Cell index
 * @param qIn Input state vector
 * @param qRhsOut Output RHS for use in an explicit time-stepping scheme
 * @return cfl frequency *if volume term is computing stable time step* (otherwise returns 0.0)
 */
GKYL_CU_DH
static inline double
gkyl_dg_eqn_vol_update(const struct gkyl_dg_eqn *eqn, 
  const double* xc, const double* dx, const int* idx,
  const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return eqn->vol_term(eqn, xc, dx, idx, qIn, qRhsOut);
}

/**
 * Compute the surface term in a DG discretization in direction dir.
 *
 * @param eqn Equation object (contains the surface kernel function pointers surf_termf_t).
 * @param dir Direction of surface update to be computed
 * @param xcL Cell center coordinates in left (or lower) cell
 * @param xcC Cell center coordinates in center cell (the cell being updated)
 * @param xcR Cell center coordinates in right (or upper) cell
 * @param dxL Grid spacing in left (or lower) cell
 * @param dxC Grid spacing in center cell (the cell being updated)
 * @param dxR Grid spacing in right (or upper) cell
 * @param idxL Cell index in left (or lower) cell
 * @param idxC Cell index in center cell (the cell being updated)
 * @param idxR Cell index in right (or upper) cell
 * @param qInL Input state vector in left (or lower) cell
 * @param qInC Input state vector in center cell (the cell being updated)
 * @param qInR Input state vector in right (or upper) cell
 * @param qRhsOut Output RHS for use in an explicit time-stepping scheme in center cell (the cell being updated)
 * @return cfl frequency *if surface term is computing stable time step* (otherwise returns 0.0)
 */
GKYL_CU_DH
static inline double
gkyl_dg_eqn_surf_update(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  return eqn->surf_term(eqn, dir, 
    xcL, xcC, xcR, 
    dxL, dxC, dxR, 
    idxL, idxC, idxR, 
    qInL, qInC, qInR, qRhsOut);
}

/**
 * Compute the boundary surface term in a DG discretization in direction dir.
 *
 * @param eqn Equation object (contains the surface kernel function pointers surf_termf_t).
 * @param dir Direction of surface update to be computed
 * @param xcEdge Cell center coordinates in edge cell (next cell in interior relative to skin cell)
 * @param xcSkin Cell center coordinates in skin cell (the cell being updated)
 * @param dxEdge Grid spacing in edge cell (next cell in interior relative to skin cell)
 * @param dxSkin Grid spacing in skin cell (the cell being updated)
 * @param idxEdge Cell index in edge cell (next cell in interior relative to skin cell)
 * @param idxSkin Cell index in skin cell (the cell being updated)
 * @param edge Parameter for whether we are at lower or upper edge for determining relative location of edge cell
 * @param qInEdge Input state vector in edge cell (next cell in interior relative to skin cell)
 * @param qInSkin Input state vector in skin cell (the cell being updated)
 * @param qRhsOut Output RHS for use in an explicit time-stepping scheme in skin cell (the cell being updated)
 * @return cfl frequency *if boundary surface term is computing stable time step* (otherwise returns 0.0)
 */
GKYL_CU_DH
static inline double
gkyl_dg_eqn_boundary_surf_update(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  return eqn->boundary_surf_term(eqn, dir, 
    xcEdge, xcSkin, 
    dxEdge, dxSkin, 
    idxEdge, idxSkin, edge, 
    qInEdge, qInSkin, qRhsOut);
}

/**
 * Compute the surface term of a diffusion operator at the boundary
 * in a DG discretization in direction dir.
 *
 * @param eqn Equation object (contains the surface kernel function pointers surf_termf_t).
 * @param dir Direction of surface update to be computed
 * @param xcEdge Cell center coordinates in edge cell (next cell in interior relative to skin cell)
 * @param xcSkin Cell center coordinates in skin cell (the cell being updated)
 * @param dxEdge Grid spacing in edge cell (next cell in interior relative to skin cell)
 * @param dxSkin Grid spacing in skin cell (the cell being updated)
 * @param idxEdge Cell index in edge cell (next cell in interior relative to skin cell)
 * @param idxSkin Cell index in skin cell (the cell being updated)
 * @param edge Parameter for whether we are at lower or upper edge for determining relative location of edge cell
 * @param qInEdge Input state vector in edge cell (next cell in interior relative to skin cell)
 * @param qInSkin Input state vector in skin cell (the cell being updated)
 * @param qRhsOut Output RHS for use in an explicit time-stepping scheme in skin cell (the cell being updated)
 * @return cfl frequency *if boundary surface term is computing stable time step* (otherwise returns 0.0)
 */
GKYL_CU_DH
static inline double
gkyl_dg_eqn_boundary_diag_update(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  return eqn->boundary_diag_term(eqn, dir, 
    xcEdge, xcSkin, 
    dxEdge, dxSkin, 
    idxEdge, idxSkin, edge, 
    qInEdge, qInSkin, qRhsOut);
}

/**
 * Delete equation object
 *
 * @param eqn Equation object to delete.
 */
void gkyl_dg_eqn_release(const struct gkyl_dg_eqn *eqn);
