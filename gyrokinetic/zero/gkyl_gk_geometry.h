#pragma once

#include <gkyl_array.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_eqn_type.h>
#include <gkyl_tok_geo.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_position_map.h>


typedef struct gk_geometry gk_geometry;

struct gk_geom_surf {

  struct gkyl_array *jacobgeo; // 1 component. Configuration space jacobian J.
  struct gkyl_array *jacobgeo_ratio; // 1 component. Ratio of jacobgeo from neigboring blocks (valid only at the boundary, stored in ghost cell).
  struct gkyl_array *bmag; // 1 component. B Magnitude of magnetic field.
  struct gkyl_array *b_i; // 3 components. Contravariant components of magnetic field vector b_1, b_2, b_3.
  struct gkyl_array *cmag; // 1 component. C = JB/sqrt(g_33).
  struct gkyl_array *jacobtot_inv; // 1 component. 1/(JB).
  struct gkyl_array *B3; // 1 component n^3 \dot \vec{B} = 1/g_33.
  struct gkyl_array *normcurlbhat; // 1 component, n^m \dot curl(bhat).
  struct gkyl_array *normals; // 9 components Cartesian components of normal. vectors in order n^1, n^2, n^3.
  struct gkyl_array *lenr; // 1 components Jc|n^i|.
  struct gkyl_array *bimpactangle; // 1 component arcsin(1/sqrt(g_33 * g^33))
  struct gkyl_array *deltats; // 1 component. Function of psi-alpha.

  // Arrays below are just for computation of arrays above
  struct gkyl_array *mc2p_nodal_fd; // 3 components. Cartesian X,Y, and Z at surf quad nodes and nodes epsilon away.
  struct gkyl_array *mc2p_nodal; // 3 components. Cartesian X,Y, and Z at surf quad nodes.
  struct gkyl_array *bmag_nodal; // 1 component. B Magnitude of magnetic field.
  struct gkyl_array *curlbhat_nodal; // Cartesian components of curl(bhat).
  struct gkyl_array *normcurlbhat_nodal; // 1 component, n^m \dot curl(bhat).
  struct gkyl_array *jacobgeo_nodal; // 1 component. Configuration space jacobian J.
  struct gkyl_array *b_i_nodal; // 3 components. Contravariant components of magnetic field vector b_1, b_2, b_3.
  struct gkyl_array *b_i_nodal_fd; // 3 components. b_i at surf quad nodes and nodes epsilon away.
  struct gkyl_array *cmag_nodal; // 1 component. C = JB/sqrt(g_33).
  struct gkyl_array *jacobtot_inv_nodal; // 1 component. 1/(JB).
  struct gkyl_array *ddtheta_nodal; // dphi/dtheta, dR/dtheta, dz/dtheta at surf quad nodes.
  struct gkyl_array *ddpsi_nodal; // dPsi/dpsi at surf quad nodes.
  struct gkyl_array *g_ij_nodal; // g_{ij}.
  struct gkyl_array *dxdz_nodal; // 9 components. Cartesian components of tangent Vectors stored in order e_1, e_2, e_3.
  struct gkyl_array *dzdx_nodal; // 9 components. Cartesian components of dual vectors stroed in order e^1, e^2, e^3.
  struct gkyl_array *dualmag_nodal; // 3 components norms of the dual vectors : sqrt(e^i.e^i).
  struct gkyl_array *normals_nodal; // 9 components Cartesian components of normal vectors in order n^1,, n^2, n^3.
  struct gkyl_array *bcart_nodal; // 3 components. Cartesian components of magnetic field unit vector b_X, b_Y, b_Z.
                     
  struct gkyl_array *B3_nodal; // 1 component n^3 \dot \vec{B} = 1/g_33.
  struct gkyl_array *lenr_nodal; // 1 components Jc|n^i|.
  struct gkyl_array *bimpactangle_nodal; // 1 component arcsin(1/sqrt(g_33 * g^33))
  struct gkyl_array *deltats_nodal; // 1 component. Function of psi-alpha.
};

struct gk_geom_corn {
  struct gkyl_array *mc2p; // 3 components. Cartesian X,Y, and Z.
  struct gkyl_array *mc2p_deflated; // cdim components. Component removed (Z in 1x, R,Z in 2x, R,Z,phi in 3x).
  struct gkyl_array *mc2nu_pos; // 3 components. Uniform computational space to non-uniform computational space mapping.
  struct gkyl_array *mc2nu_pos_deflated; // cdim components. Uniform computational space to non-uniform computational space mapping.
  struct gkyl_array *bmag; // 1 component. B Magnitude of magnetic field.
  struct gkyl_array *bmag_inv; // 1 component. 1/bmag.

  // Arrays below are just for computation of arrays above.
  struct gkyl_array *mc2p_nodal; // 3 components. Cartesian X,Y, and Z.
  struct gkyl_array *mc2nu_pos_nodal; // 3 components. Uniform computational space to non-uniform computational space mapping
  struct gkyl_array *bmag_nodal; // 1 components. Magnitude of Magnetic Field.

};

struct gk_geom_int {
  struct gkyl_array *mc2p; // 3 components. Cartesian X, Y and Z,
  struct gkyl_array *bmag; // 1 component. B Magnitude of magnetic field.
  struct gkyl_array *g_ij; // 6 components. Metric coefficients g_{ij} Stored in order g_11, g12, g_13, g_22, g_23, g_33.
  struct gkyl_array *g_ij_neut; // 6 components.  Metric coefficients g_{ij} Stored in order g_11, g12, g_13, g_22, g_23, g_33.
                                // Calculated with coord definition alpha = phi for tokamak geometry.
  struct gkyl_array *dxdz; // 9 components. Cartesian components of tangent Vectors stored in order e_1, e_2, e_3.
  struct gkyl_array *dzdx; // 9 components. Cartesian components of dual vectors stored in order e^1, e^2, e^3.
  struct gkyl_array *dualmag; // 3 components. Norms of the dual vectors: sqrt(e^i.e^i).
  struct gkyl_array *normals; // 9 components. Cartesian components of normal vectors in order n^1,, n^2, n^3.
  struct gkyl_array *jacobgeo; // 1 component. Configuration space jacobian J.
  struct gkyl_array *jacobgeo_ghost; // 1 component. Configuration space jacobian J.
  struct gkyl_array *jacobgeo_inv; // 1 component. 1/J
  struct gkyl_array *gij; // Metric coefficients g^{ij}. See g_ij for order.
  struct gkyl_array *gij_neut; // Metric coefficients g^{ij}. See g_ij for order.
                               // Calculated with coord definition alpha = phi for tokamak geometry.
  struct gkyl_array *b_i; // 3 components. Covariant components of magnetic field unit vector b_1, b_2, b_3.
  struct gkyl_array *bcart; // 3 components. Cartesian components of magnetic field unit vector b_X, b_Y, b_Z.
  struct gkyl_array *cmag; // 1 component. C = JB/sqrt(g_33).
  struct gkyl_array *jacobtot; // 1 component. Phase space Jacobian = JB.
  struct gkyl_array *jacobtot_inv; // 1 component. 1/(JB).
  struct gkyl_array *gxxj; // 1 component. g^{xx} * J. For poisson solve.
  struct gkyl_array *gxyj; // 1 component. g^{xy} * J. For poisson solve.
  struct gkyl_array *gyyj; // 1 component. g^{yy} * J. For poisson solve.
  struct gkyl_array *gxzj; // 1 component. g^{xz} * J. For poisson solve if z derivatives are kept.
  struct gkyl_array *eps2; // 1 component. eps2 = Jg^33 - J/g_33. For poisson if z derivatives are kept.
  struct gkyl_array *dualcurlbhat; // 3 components, e^m \dot curl(bhat).
  struct gkyl_array *dualcurlbhatoverB; // 3 components, e^m \dot curl(bhat)/|B|.
  struct gkyl_array *rtg33inv; // 1 component 1/sqrt(g_33).
  struct gkyl_array *bioverJB; // 1 component b_i/J/|B|.
  struct gkyl_array *B3; // 1 component e^3 \dot \vec{B} = 1/g_33.
  struct gkyl_array *qprofile; // 1 component. Flux surface averaged q profle q(psi).

  // Arrays below are just for computation of arrays above
  struct gkyl_array *bmag_nodal;
  struct gkyl_array *ddtheta_nodal;
  struct gkyl_array *ddpsi_nodal;// dPsi/dpsi at interior quad nodes
  struct gkyl_array *mc2p_nodal; // 3 components. Cartesian X,Y, and Z
  struct gkyl_array *mc2p_nodal_fd; // 39 components. Cartesian X,Y, and Z at nodes and FD nodes.
  /* Array containing cartesian coordinates at nodes and nearby nodes (epsilon and 2 epsilon away) used for FD
  *    At each array location 39 values are stored.
  *    The arrangement is as follows: X_c, Y_c, Z_c, 
  *    X_L1, Y_L1, Z_L1, X_R1, Y_R1, Z_R1,
  *    X_L2, Y_L2, Z_L2, X_R2, Y_R2, Z_R2,
  *    X_L3, Y_L3, Z_L3, X_R3, Y_R3, Z_R3,
  *    X_LL1, Y_LL1, Z_LL1, X_RR1, Y_RR1, Z_RR1,
  *    X_LL2, Y_LL2, Z_LL2, X_RR2, Y_RR2, Z_RR2,
  *    X_LL3, Y_LL3, Z_LL3, X_RR3, Y_RR3, Z_RR3
  *    where L#/R# indicates a node shifted to the left/right by epsilon in coordinate #
  *    and LL#/RR# indicates a node shifted to the left/right by 2 epsilon in coordinate #
  */
  struct gkyl_array *curlbhat_nodal; // Cartesian components of curl(bhat)
  struct gkyl_array *dualcurlbhat_nodal; // 3 components, e^m \dot curl(bhat)
  struct gkyl_array *jacobgeo_nodal; // jacobian 
  struct gkyl_array *g_ij_nodal; // 6 components. 
                           // Metric coefficients g_{ij} Stored in order g_11, g12, g_13, g_22, g_23, g_33
  struct gkyl_array *g_ij_neut_nodal; // 6 components. 
                           // Metric coefficients g_{ij} Stored in order g_11, g12, g_13, g_22, g_23, g_33
                           // Calculated with coord definition alpha = phi for tokamak geometry
  struct gkyl_array *dxdz_nodal; // 9 components.
                           // Cartesian components of tangent Vectors stored in order e_1, e_2, e_3
  struct gkyl_array *dzdx_nodal; // 9 components.
                           // Cartesian components of dual vectors stroed in order e^1, e^2, e^3
  struct gkyl_array *dualmag_nodal; // 3 components
                              // norms of the dual vectors : sqrt(e^i.e^i)
  struct gkyl_array *normals_nodal; // 9 components
                              // Cartesian components of normal vectors in order n^1,, n^2, n^3
  struct gkyl_array *gij_neut_nodal; // Metric coefficients g^{ij}. See g_ij for order. 
                               // Calculated with coord definition alpha = phi for tokamak geometry
  struct gkyl_array *b_i_nodal; // 3 components. Covariant components of magnetic field unit vector b_1, b_2, b_3.
  struct gkyl_array *b_i_nodal_fd; // 3 components. b_i at interior quad nodes and nodes epsilon away
  struct gkyl_array *bcart_nodal; // 3 components. Cartesian components of magnetic field unit vector b_X, b_Y, b_Z.
  struct gkyl_array *B3_nodal; // 1 component e^3 \dot \vec{B} = 1/g_33 
  struct gkyl_array *dualcurlbhatoverB_nodal; // 3 components, e^m \dot curl(bhat)/|B|
  struct gkyl_array *rtg33inv_nodal; // 1 component 1/sqrt(g_33)
  struct gkyl_array *bioverJB_nodal; // 3 components b_i/J/|B|
};

struct gk_geometry {
  enum gkyl_geometry_id geometry_id;

  // Objects for mapc2p and finite differences array.
  struct gkyl_range local;
  struct gkyl_range local_ext;
  struct gkyl_range global;
  struct gkyl_range global_ext;
  struct gkyl_basis basis;
  struct gkyl_basis surf_basis;
  int num_surf_basis;
  struct gkyl_rect_grid grid;
  double dzc[3];

  // Nodal Ranges.
  struct gkyl_range nrange_corn;
  struct gkyl_range nrange_int;
  struct gkyl_range nrange_surf[3];

  // The fields in these structs contain the geometric quantities needed to solve the
  // GK Equation and Poisson Equation and to apply certain BC's.
  struct gk_geom_corn geo_corn; // Volume geometry from corner Nodes.
  struct gk_geom_int geo_int; // Volume geometry from interior nodes.
  struct gk_geom_surf geo_surf[3]; // Surface geometry.

  int geqdsk_sign_convention; // 0 if psi increases away from magnetic axis.
                              // 1 if psi increases toward magnetic axis.

  int half_domain; // For double null geometry. 0 for full domain 1 for lower half domain.

  bool has_LCFS; // Whether the geometry has an LCFS.
  double x_LCFS; // For mapc2p IWL geometry, the user has to provide the
                 // location of the LCFS. For numerical IWL, it may be stored
                 // in the eqdsk.
  int idx_LCFS_lo; // Index of the cell that abuts the LCFS from below.

  // Functions defining the twistshift for parallel BCs.
  void (*parallel_lower_bc_shift_func)(double t, const double *xn, double *fout, void *ctx);
  void (*parallel_upper_bc_shift_func)(double t, const double *xn, double *fout, void *ctx);
  void *parallel_lower_bc_shift_ctx; // Context for parallel_lower_bc_shift_func.
  void *parallel_upper_bc_shift_ctx; // Context for parallel_upper_bc_shift_func.

  struct gkyl_msgpack_map_elem* io_meta_basic; // Metadata for I/O.
  int io_meta_basic_len; // Number of elements in io_meta_basic.
  
  uint32_t flags;
  struct gkyl_ref_count ref_count;  
  struct gk_geometry *on_dev; // Pointer to itself or device object.
};


// Inputs to create geometry for a specific computational grid.
struct gkyl_mirror_geo_grid_inp {
  char filename_psi[256]; // File with psi(R,Z) data.
  double rclose; // Closest R to discrimate.
  double rright; // Closest R to discrimate.
  double zmin, zmax; // Extents of Z for integration.

  bool include_axis; // Add nodes on r=0 axis (the axis is assumed be psi=0).
  enum gkyl_mirror_grid_gen_field_line_coord fl_coord; // Field-line coordinate to use.
};

// Input struct for geometry creation
struct gkyl_gk_geometry_inp {
  enum gkyl_geometry_id geometry_id;
  char geometry_path[128]; // Path to geometry files

  void *c2p_ctx; // Context for mapc2p function.
  // Pointer to mapc2p function: xc are the computational space
  // coordinates and on output xp are the corresponding physical space
  // coordinates.
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);

  void *bfield_ctx; // Context for bfield function.
  // Pointer to bfield function.
  void (*bfield_func)(double t, const double *xc, double *xp, void *ctx);

  // Functions defining the twistshift for parallel BCs.
  void (*parallel_lower_bc_shift_func)(double t, const double *xn, double *fout, void *ctx);
  void (*parallel_upper_bc_shift_func)(double t, const double *xn, double *fout, void *ctx);
  void *parallel_lower_bc_shift_ctx; // Context for parallel_lower_bc_shift_func.
  void *parallel_upper_bc_shift_ctx; // Context for parallel_upper_bc_shift_func.

  struct gkyl_efit_inp efit_info; // Context with RZ data such as efit file for a tokamak or mirror.
  struct gkyl_tok_geo_grid_inp tok_grid_info; // Context for tokamak geometry with computational domain info.
  struct gkyl_mirror_geo_grid_inp mirror_grid_info; // Context for mirror geometry with computational domain info.
  struct gkyl_position_map *position_map; // Position map object.
  struct gkyl_comm *comm; // Communicator object.

  double world[3]; // Extra computational coordinates for cases with reduced dimensionality.

  bool has_LCFS; // Whether the geometry has a last closed flux surface (LCFS).
  double x_LCFS; // x location of the LCFS.

  // 3D grid ranges and basis.
  struct gkyl_rect_grid geo_grid;
  struct gkyl_range geo_local;
  struct gkyl_range geo_local_ext;
  struct gkyl_range geo_global;
  struct gkyl_range geo_global_ext;
  struct gkyl_basis geo_basis;

  // Grid ranges and basis with cdim of simulation.
  struct gkyl_rect_grid grid;
  struct gkyl_range local;
  struct gkyl_range local_ext;
  struct gkyl_range global;
  struct gkyl_range global_ext;
  struct gkyl_basis basis;

};

/**
 * Create a new geometry object
 * If use-gpu is false, the geometric quantities will be empty, so quantities should be read 
 * from file later on
 * If use-gpu is true, the geometric quantities will be copied from the pre-populated host object
 *
 * @param geo_host gk_geometry object on the host
 *   Unused if calling with use_gpu=false
 *   If use_gpu=true, geo_host must already be initialized.
 * @param geometry_inp geometry input struct containing grid, range, and other geo info
 * @param use_gpu whether or not to use gpu
 */
struct gk_geometry*
gkyl_gk_geometry_new(struct gk_geometry* geo_host, struct gkyl_gk_geometry_inp *geometry_inp, bool use_gpu);

/**
 * Create a new gk geometry object that lives on NV-GPU from a host geometry object: see new() method
 * above for documentation.
 */
struct gk_geometry* 
gkyl_gk_geometry_cu_dev_new(struct gk_geometry* geo_host, struct gkyl_gk_geometry_inp *geometry_inp);

/**
 * Augment a grid with dim < 3 to 3d by adding 1 cell in the other directions
 * If dim=1, the input grid is assumed to be in z
 * If dim =2, the input grid is assumed to be in xz
 *
 * @param grid Input grid with dim <3.
 * @param geometry Geometry input struct with context for augmenting grid.
 */
struct gkyl_rect_grid 
gkyl_gk_geometry_augment_grid(struct gkyl_rect_grid grid, struct gkyl_gk_geometry_inp geometry);


/**
 * Augment a range with dim < 3 to 3d by adding 1 cell in the other directions
 * If dim=1, the input range is assumed to be in z
 * If dim =2, the input range is assumed to be in xz
 *
 * @param inrange Input range with dim <3.
 * @param nghost Number of ghost cells.
 * @param ext_range Output, augmented extended range.
 * @param range Output, augmented range.
 */
void 
gkyl_gk_geometry_augment_local(const struct gkyl_range *inrange, const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range);

/**
 * Reduce bmag to get min or max value, by evaluating bmag at basis nodes.
 * Only to be used during initialization because it allocates memory.
 *
 *  @param up gk_geometry object.
 *  @param op Pperation to perform (GKYL_MAX or GKYL_MIN).
 */
double gkyl_gk_geometry_reduce_bmag(struct gk_geometry* up, enum gkyl_array_op op);

/**
 * Reduce bmag to get min or max value, by evaluating bmag at basis nodes,
 * and the location of the extrema.
 * Only to be used during initialization because it allocates memory.
 *
 * @param up gk_geometry object.
 * @param op Operation to perform (GKYL_MAX or GKYL_MIN).
 * @param op Coordinate where extrema occurs.
 */
double gkyl_gk_geometry_reduce_arg_bmag(struct gk_geometry* up, enum gkyl_array_op op, double *coord);

/**
 * Init nodal range from modal range
 *
 * @param nrange nodal Range to be initialized.
 * @param range modal Range.
 * @param poly_order Polynomial order.
 */
void
gkyl_gk_geometry_init_nodal_range( struct gkyl_range *nrange, struct gkyl_range *range, int poly_order);

/**
 * Init nodal grid from modal grid
 *
 * @param ngrid Nodal grid to be initialized.
 * @param grid Modal grid.
 * @param nrange Nodal range.
 */
void
gkyl_gk_geometry_init_nodal_grid(struct gkyl_rect_grid *ngrid, struct gkyl_rect_grid *grid, struct gkyl_range *nrange);

/**
 * Deflate geometry to lower dimensionality.
 *
 * @param up_3d 3d geometry object to deflate
 * @param grid deflated grid
 * @param local deflated local range
 * @param local_ext deflated local extended range
 * @param basis deflated basis
 * @param use_gpu whether or not to use gpu
 */
struct gk_geometry* gkyl_gk_geometry_deflate(const struct gk_geometry* up_3d, struct gkyl_gk_geometry_inp *geometry_inp);

/**
 * Populate nodal arrays from modal geometry
 * @param gk_geom geometry object to deflate
 */
void gkyl_gk_geometry_populate_nodal(struct gk_geometry *gk_geom);

/**
 * Write psi(R,Z)
* @param geometry_inp Geometry input struct.
* @param io_meta_basic_basic basic metadata
* @param io_meta_basic_basic_len length of basic metadata
 */
void gkyl_gk_geometry_write_efit(struct gkyl_gk_geometry_inp *geometry_inp, struct gkyl_msgpack_map_elem* io_meta_basic_basic, int io_meta_basic_basic_len);

/**
 * Reset the metadata values with corresponding values in GK geometry object
 * (they may have been updated).
 *
 * @param up GK geometry object.
 */
void gkyl_gk_geometry_reset_io_meta(struct gk_geometry *up);

/**
 * Acquire pointer to gk geometry object. The pointer must be released
 * using gkyl_gk_geometry_release method.
 *
 * @param up Geometry to which a pointer is needed
 * @return Pointer to acquired geometry
 */
struct gk_geometry* gkyl_gk_geometry_acquire(const struct gk_geometry* up);

void gkyl_gk_geometry_free(const struct gkyl_ref_count *ref);

/**
 * Release gk geometry object.
 *
 * @param up gk geometry object to release.
 */
void gkyl_gk_geometry_release(const struct gk_geometry *up);
