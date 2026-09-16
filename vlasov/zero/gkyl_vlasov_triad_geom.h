#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Triad geometry context if needed for preset geometries.
struct gkyl_triad_geom_ctx {
  double mass_bh; // Mass of the black hole (M)
  double spin_bh; // Spin of the black hole (a = J/M)
};

enum gkyl_triad_preset_geom_type {
  GKYL_TRIAD_NONE = 0,
  GKYL_TRIAD_FLAT,
  GKYL_TRIAD_ANNULUS,
  GKYL_TRIAD_CYLINDRICAL_RZ,
  GKYL_TRIAD_SPHERICAL_RTHETA,
  GKYL_TRIAD_GR_KERR_SCHILD_RPHI,
  GKYL_TRIAD_GR_KERR_SCHILD_R,
  GKYL_TRIAD_GR_KERR_SCHILD_RTHETA,
  GKYL_TRIAD_GR_KERR_SCHILD_3V,
  GKYL_TRIAD_CART_GR_KERR_SCHILD_3V,
};

// Function transforming a set of cdim computational coordinates to physical
// ones (e.g. the position map on non-uniform conf meshes). NULL => identity.
typedef void (*vlasov_triad_geom_c2p_t)(const double *xcomp, double *xphys, void *ctx);

struct gkyl_vlasov_triad_geom_inp {
  evalf_t eval_cov_tangent_basis; // The covariant tangent basis to be evaluated within each configuration space cell.
  evalf_t eval_triad_basis; // The triad basis to be evaluated within each configuration space cell.
  evalf_t eval_triad_basis_gradient; // The triad basis gradient to be evaluated within each configuration space cell.
  void *eval_cov_tangent_basis_ctx; // Context for function evaluation. Can be NULL.
  void *eval_triad_basis_ctx; // Context for function evaluation. Can be NULL.
  void *eval_triad_basis_gradient_ctx; // Context for function evaluation. Can be NULL.

  evalf_t eval_vierbein; // The vierbein to be evaluated within each configuration space cell.
  evalf_t eval_vierbein_gradient; // The vierbein gradient to be evaluated within each configuration space cell.
  void *eval_vierbein_ctx; // Context for function evaluation. Can be NULL.
  void *eval_vierbein_gradient_ctx; // Context for function evaluation. Can be NULL.

  bool use_preset_geom; // bool for determining if we have specified a preset geometry.
  bool use_vierbein; // bool for determining which geometry convention we are constructing PT from
  enum gkyl_triad_preset_geom_type triad_preset_geom_type;  // geom type for preset geometries for triads

  vlasov_triad_geom_c2p_t c2p_func; // comp->phys conf-coordinate map for sampling the
                                    // geometry on non-uniform meshes. NULL => identity.
  void *c2p_func_ctx; // Context for c2p_func.
};

/**
 * Construct the triad geometry objects
 *
 * @param cgrid Configuration-space grid object
 * @param crange Configuration-space range
 * @param cbasis Configuration-space basis
 * @param pgrid Phase-space grid object
 * @param prange Phase-space range
 * @param pbasis Phase-space basis
 * @param inp_triad_geom Basis mapping input (function and context) cov_tangent_basis, traid_basis
 * @param conf_poisson_tensor The configuration component of the Poisson tensor
 */
void gkyl_vlasov_triad_geom_new(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, const struct gkyl_basis cbasis, 
  const struct gkyl_rect_grid *pgrid, const struct gkyl_range *prange, const struct gkyl_basis pbasis, 
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom, struct gkyl_array *conf_poisson_tensor);

/**
 * Construct the configuration-space Poisson tensor directly from nodal samples of
 * the geometry, bypassing the evalf callbacks and the comp->phys coordinate map
 * (for callers, such as gyrokinetics, that already hold the geometry as nodal
 * arrays). Each input array is indexed on crange with a per-cell, node-major
 * layout in cbasis.node_list order: for node i the components are contiguous,
 * so component c of node i lives at [i*ncomp_per_node + c].
 *
 * @param cgrid Configuration-space grid object
 * @param crange Configuration-space range
 * @param cbasis Configuration-space basis
 * @param pgrid Phase-space grid object
 * @param prange Phase-space range
 * @param pbasis Phase-space basis
 * @param cov_tangent_basis_nodal Cartesian components of the covariant tangent basis
 *   at conf-basis nodes, vdim*vdim per node (row a = tangent e_a)
 * @param triad_basis_nodal Cartesian components of the orthonormal triad at
 *   conf-basis nodes, vdim*vdim per node (row a = frame vector e_a)
 * @param triad_basis_gradient_nodal d(triad)/dz^i at conf-basis nodes,
 *   vdim*vdim*vdim per node (gradient direction slowest)
 * @param conf_poisson_tensor The configuration component of the Poisson tensor (output)
 */
void gkyl_vlasov_triad_geom_from_nodal(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, const struct gkyl_basis cbasis,
  const struct gkyl_rect_grid *pgrid, const struct gkyl_range *prange, const struct gkyl_basis pbasis,
  const struct gkyl_array *cov_tangent_basis_nodal, const struct gkyl_array *triad_basis_nodal,
  const struct gkyl_array *triad_basis_gradient_nodal, struct gkyl_array *conf_poisson_tensor);

/**
 * Same as gkyl_vlasov_triad_geom_from_nodal, but for nodal samples living at the
 * per-cell Gauss-Legendre interior nodes of a global nodal range (the layout the
 * gyrokinetic geometry stores, cf. gkyl_nodal_ops_n2m_interior): 2 nodes per
 * direction per cell, none shared, node i of a cell decomposed lexicographically
 * with the last dimension fastest. The Poisson tensor is assembled pointwise at
 * the nodes and converted to modal with quad_nodal_to_modal. p=1 only.
 *
 * @param cgrid Configuration-space grid object
 * @param crange Configuration-space range (cells to update)
 * @param cbasis Configuration-space basis (p=1)
 * @param pgrid Phase-space grid object
 * @param prange Phase-space range
 * @param pbasis Phase-space basis
 * @param nrange Nodal range holding the interior nodes (2*cells per direction)
 * @param cov_tangent_basis_nodal Cartesian covariant tangents on nrange, vdim*vdim per node
 * @param triad_basis_nodal Orthonormal triad on nrange, vdim*vdim per node
 * @param triad_basis_gradient_nodal d(triad)/dz^i on nrange, vdim*vdim*vdim per node
 * @param conf_poisson_tensor The configuration component of the Poisson tensor (output)
 */
void gkyl_vlasov_triad_geom_from_nodal_interior(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, const struct gkyl_basis cbasis,
  const struct gkyl_rect_grid *pgrid, const struct gkyl_range *prange, const struct gkyl_basis pbasis,
  const struct gkyl_range *nrange,
  const struct gkyl_array *cov_tangent_basis_nodal, const struct gkyl_array *triad_basis_nodal,
  const struct gkyl_array *triad_basis_gradient_nodal, struct gkyl_array *conf_poisson_tensor);

/**
 * Build the b-aligned orthonormal triad and its gradients from nodal Cartesian
 * covariant tangents in FIELD-ALIGNED coordinates (the magnetic field along the
 * third tangent, as in gyrokinetic geometry), then assemble the conf-space
 * Poisson tensor. Per node: bhat = e_3/|e_3|, e1 = normalize(t_1 - (t_1.bhat)bhat),
 * e2 = bhat x e1 (right-handed). Triad gradients are second-order finite
 * differences across the interior node grid. Construction-time checks: the
 * Gram-Schmidt degeneracy (t_1 parallel to bhat) and, when bhat_check_nodal is
 * given, agreement of bhat with the independently computed field direction
 * (e.g. gyrokinetic bcart_nodal); failures print and, with exit_at_checks,
 * assert. 3x3v and p=1 only; deflated axisymmetric conf grids need analytic
 * d/dalpha terms and are not supported yet.
 *
 * @param cgrid Configuration-space grid object
 * @param crange Configuration-space range (cells to update)
 * @param cbasis Configuration-space basis (p=1)
 * @param pgrid Phase-space grid object
 * @param prange Phase-space range
 * @param pbasis Phase-space basis
 * @param nrange Nodal range holding the interior nodes (2*cells per direction)
 * @param cov_tangent_basis_nodal Cartesian covariant tangents on nrange, 9 per node
 *   (e.g. gyrokinetic dxdz_nodal)
 * @param bhat_check_nodal Cartesian field direction on nrange, 3 per node, or NULL
 *   to skip the cross-check
 * @param exit_at_checks Assert when a construction-time check fails
 * @param conf_poisson_tensor The configuration component of the Poisson tensor (output)
 */
void gkyl_vlasov_triad_geom_from_tangents_interior(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, const struct gkyl_basis cbasis,
  const struct gkyl_rect_grid *pgrid, const struct gkyl_range *prange, const struct gkyl_basis pbasis,
  const struct gkyl_range *nrange,
  const struct gkyl_array *cov_tangent_basis_nodal, const struct gkyl_array *bhat_check_nodal,
  bool exit_at_checks, struct gkyl_array *conf_poisson_tensor);

/**
 * Preset function for the triad hamil
 *
 * @param cdim Number of configuration space dimenions
 * @param vdim Number of velocity space dimenions
 * @param preset_geom_type Preset geometry type
 */
evalf_t gkyl_vlasov_triad_preset_hamil(const int cdim, const int vdim, enum gkyl_triad_preset_geom_type preset_geom_type);

/**
 * Preset function for the covariant vierbein.
 *
 * @param vdim Number of velocity-space dimensions
 * @param preset_geom_type Preset geometry type
 */
evalf_t gkyl_vlasov_triad_preset_vierbein(const int vdim, enum gkyl_triad_preset_geom_type preset_geom_type);

/**
 * Preset function for the contravariant vierbein.
 *
 * @param vdim Number of velocity-space dimensions
 * @param preset_geom_type Preset geometry type
 */
evalf_t gkyl_vlasov_triad_preset_vierbein_inv(const int vdim, enum gkyl_triad_preset_geom_type preset_geom_type);
