#pragma once

#include <gkyl_comm.h>

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
  GKYL_TRIAD_CART_GR_KERR_SCHILD_3V
};

// Function transforming a set of cdim computational coordinates to physical
// ones (e.g. the position map on non-uniform conf meshes). NULL => identity.
typedef void (*vlasov_triad_geom_c2p_t)(const double *xcomp, double *xphys, void *ctx);

struct gkyl_vlasov_triad_geom_inp {
  evalf_t
    eval_cov_tangent_basis; // The covariant tangent basis to be evaluated within each configuration space cell.
  evalf_t eval_triad_basis; // The triad basis to be evaluated within each configuration space cell.
  evalf_t
    eval_triad_basis_gradient; // The triad basis gradient to be evaluated within each configuration space cell.
  void *eval_cov_tangent_basis_ctx; // Context for function evaluation. Can be NULL.
  void *eval_triad_basis_ctx; // Context for function evaluation. Can be NULL.
  void *eval_triad_basis_gradient_ctx; // Context for function evaluation. Can be NULL.

  evalf_t eval_vierbein; // The vierbein to be evaluated within each configuration space cell.
  evalf_t
    eval_vierbein_gradient; // The vierbein gradient to be evaluated within each configuration space cell.
  void *eval_vierbein_ctx; // Context for function evaluation. Can be NULL.
  void *eval_vierbein_gradient_ctx; // Context for function evaluation. Can be NULL.

  bool use_preset_geom; // bool for determining if we have specified a preset geometry.
  bool use_vierbein; // bool for determining which geometry convention we are constructing PT from
  enum gkyl_triad_preset_geom_type
    triad_preset_geom_type; // geom type for preset geometries for triads

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
void gkyl_vlasov_triad_geom_new(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange,
  const struct gkyl_basis cbasis, const struct gkyl_rect_grid *pgrid,
  const struct gkyl_range *prange, const struct gkyl_basis pbasis,
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom, struct gkyl_array *conf_poisson_tensor
);

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
void gkyl_vlasov_triad_geom_from_nodal(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange,
  const struct gkyl_basis cbasis, const struct gkyl_rect_grid *pgrid,
  const struct gkyl_range *prange, const struct gkyl_basis pbasis,
  const struct gkyl_array *cov_tangent_basis_nodal, const struct gkyl_array *triad_basis_nodal,
  const struct gkyl_array *triad_basis_gradient_nodal, struct gkyl_array *conf_poisson_tensor
);

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
void gkyl_vlasov_triad_geom_from_nodal_interior(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange,
  const struct gkyl_basis cbasis, const struct gkyl_rect_grid *pgrid,
  const struct gkyl_range *prange, const struct gkyl_basis pbasis, const struct gkyl_range *nrange,
  const struct gkyl_array *cov_tangent_basis_nodal, const struct gkyl_array *triad_basis_nodal,
  const struct gkyl_array *triad_basis_gradient_nodal, struct gkyl_array *conf_poisson_tensor
);

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
void gkyl_vlasov_triad_geom_from_tangents_interior(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange,
  const struct gkyl_basis cbasis, const struct gkyl_rect_grid *pgrid,
  const struct gkyl_range *prange, const struct gkyl_basis pbasis, const struct gkyl_range *nrange,
  const struct gkyl_array *cov_tangent_basis_nodal, const struct gkyl_array *bhat_check_nodal,
  bool exit_at_checks, struct gkyl_array *conf_poisson_tensor
);

/**
 * Same b-aligned triad as gkyl_vlasov_triad_geom_from_tangents_interior, for
 * configuration grids that hold only cdim <= 3 of the three field-aligned
 * coordinates (e.g. a deflated 1x or 2x gyrokinetic geometry), where the triad
 * gradients along the coordinates absent from the grid cannot be taken by finite
 * differences. The caller supplies, at the interior nodes, the Cartesian
 * tangents of all three coordinates and their gradients along all three
 * coordinates (e.g. the gyrokinetic geometry's dxdz_nodal and dxdz_grad_nodal),
 * and the triad gradients follow pointwise from the chain rule of the
 * Gram-Schmidt construction; nothing is differentiated here. The tangent rows
 * and gradient directions are in the 3D coordinate order, in which the third
 * tangent is field-aligned; dir_map gives, for each coordinate index of the
 * Poisson tensor (grid directions first, then the absent ones), its index in
 * that 3D order. Construction-time checks as in the finite-difference variant.
 * p=1 only.
 *
 * @param cgrid Configuration-space grid object (cdim <= 3)
 * @param crange Configuration-space range (cells to update)
 * @param cbasis Configuration-space basis (p=1)
 * @param pgrid Phase-space grid object
 * @param prange Phase-space range
 * @param pbasis Phase-space basis
 * @param nrange Nodal range holding the interior nodes (2*cells per direction)
 * @param cov_tangent_basis_nodal Cartesian covariant tangents on nrange, 9 per node
 *   (row a = tangent of 3D coordinate a)
 * @param cov_tangent_basis_grad_nodal d(tangents)/dz^j on nrange, 27 per node,
 *   gradient direction j (3D order) slowest: [9 j + 3 a + c]
 * @param bhat_check_nodal Cartesian field direction on nrange, 3 per node, or NULL
 *   to skip the cross-check
 * @param dir_map 3D coordinate index of each Poisson-tensor coordinate index, or
 *   NULL for the identity (cdim = 3)
 * @param exit_at_checks Assert when a construction-time check fails
 * @param conf_poisson_tensor The configuration component of the Poisson tensor (output)
 */
void gkyl_vlasov_triad_geom_from_tangents_grad_interior(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange,
  const struct gkyl_basis cbasis, const struct gkyl_rect_grid *pgrid,
  const struct gkyl_range *prange, const struct gkyl_basis pbasis, const struct gkyl_range *nrange,
  const struct gkyl_array *cov_tangent_basis_nodal,
  const struct gkyl_array *cov_tangent_basis_grad_nodal, const struct gkyl_array *bhat_check_nodal,
  const int *dir_map, bool exit_at_checks, struct gkyl_array *conf_poisson_tensor
);

/**
 * Poisson tensor of a triad given by evaluation functions (vierbein, preset vierbein,
 * or tangent basis + triad + triad gradient, as in gkyl_vlasov_triad_geom_new), stored
 * at the configuration-space Gauss-Legendre nodes for the nodal (tensor hybrid p=1,
 * sparse velocity Hamiltonian) triad kernels, vdim = 2 or 3. pt_vol and pt_surf have
 * the layout of gkyl_vlasov_triad_geom_nodal_from_tangents; every value is evaluated
 * pointwise (at the physical coordinates of the node if a c2p map is given), so the
 * face values are single valued.
 *
 * @param cgrid Configuration-space grid
 * @param crange Configuration-space range to fill
 * @param vdim Number of velocity dimensions
 * @param inp_triad_geom Geometry evaluation functions
 * @param pt_vol Volume storage (output, host)
 * @param pt_surf Face storage (output, host)
 */
void gkyl_vlasov_triad_geom_nodal_new(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, int vdim,
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom, struct gkyl_array *pt_vol,
  struct gkyl_array *pt_surf
);

/**
 * Poisson tensor of the b-aligned triad (e1, e2, bhat) of a field-aligned geometry,
 * vdim = 3, stored at the configuration-space Gauss-Legendre nodes for the nodal
 * (tensor hybrid p=1, sparse velocity Hamiltonian) triad kernels:
 *   pt_vol:  per cell, [node][component] at the 2^cdim interior nodes (18 components),
 *            nodes in lexicographic order with the first direction slowest and the
 *            -1/sqrt(3) node first;
 *   pt_surf: per cell, [dir][face][face node][b] with face 0 the lower and 1 the upper
 *            face of the cell, the 2^(cdim-1) face nodes of the transverse directions
 *            in the same order, and b the conf-velocity row of dir (Pi^{dir,b}),
 *            evaluated from the geometry at the face (single valued across the face).
 * Triad gradients along the grid directions are finite differences across the
 * interior nodes; one cell of nodes is exchanged across every internal decomposition
 * boundary, so the result does not depend on the decomposition. Along the coordinates
 * removed from a deflated geometry they come from the tangent gradients.
 *
 * @param cgrid Configuration-space grid
 * @param local Local configuration-space range
 * @param local_ext Local extended configuration-space range
 * @param global Global configuration-space range
 * @param comm Communicator used to exchange the halo
 * @param use_gpu Whether comm exchanges device arrays (the construction itself is on host)
 * @param nrange_int Nodal range of the interior nodes (2*cells of local per direction)
 * @param tan_nodal Cartesian covariant tangents on nrange_int, 9 per node (3D order)
 * @param tan_grad_nodal d(tangents)/dz^j on nrange_int, 27 per node, or NULL when
 *   cdim = 3 (only read along removed coordinates)
 * @param bhat_nodal Cartesian field direction on nrange_int, or NULL to skip the check
 * @param nrange_surf Nodal ranges of the face nodes of each grid direction
 * @param tan_surf_nodal Cartesian covariant tangents at the face nodes of each direction
 * @param dir_map 3D coordinate index of each Poisson-tensor coordinate index (grid
 *   directions first), or NULL for the identity (cdim = 3)
 * @param exit_at_checks Assert when a construction-time check fails
 * @param pt_vol Volume storage (output, host, on local_ext)
 * @param pt_surf Face storage (output, host, on local_ext)
 */
void gkyl_vlasov_triad_geom_nodal_from_tangents(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *local,
  const struct gkyl_range *local_ext, const struct gkyl_range *global, struct gkyl_comm *comm,
  bool use_gpu, const struct gkyl_range *nrange_int, const struct gkyl_array *tan_nodal,
  const struct gkyl_array *tan_grad_nodal, const struct gkyl_array *bhat_nodal,
  const struct gkyl_range *nrange_surf, struct gkyl_array *const *tan_surf_nodal,
  const int *dir_map, bool exit_at_checks, struct gkyl_array *pt_vol, struct gkyl_array *pt_surf
);

/**
 * Preset function for the triad hamil
 *
 * @param cdim Number of configuration space dimenions
 * @param vdim Number of velocity space dimenions
 * @param preset_geom_type Preset geometry type
 */
evalf_t gkyl_vlasov_triad_preset_hamil(
  const int cdim, const int vdim, enum gkyl_triad_preset_geom_type preset_geom_type
);

/**
 * Preset function for the covariant vierbein.
 *
 * @param vdim Number of velocity-space dimensions
 * @param preset_geom_type Preset geometry type
 */
evalf_t gkyl_vlasov_triad_preset_vierbein(
  const int vdim, enum gkyl_triad_preset_geom_type preset_geom_type
);

/**
 * Preset function for the contravariant vierbein.
 *
 * @param vdim Number of velocity-space dimensions
 * @param preset_geom_type Preset geometry type
 */
evalf_t gkyl_vlasov_triad_preset_vierbein_inv(
  const int vdim, enum gkyl_triad_preset_geom_type preset_geom_type
);
