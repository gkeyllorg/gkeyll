#pragma once

/* Basis function identifiers */
enum gkyl_basis_type {
  GKYL_BASIS_MODAL_SERENDIPITY,
  GKYL_BASIS_MODAL_TENSOR,
  GKYL_BASIS_MODAL_HYBRID,
  GKYL_BASIS_MODAL_GKHYBRID,
  GKYL_BASIS_MODAL_GKHYBRID_VEL,
};

typedef void (*nodal_to_modal_quad_surf_t)(const double *fnodal, double *fmodal);
typedef void (*node_quad_surf_list_t)(double *node_coords);

/**
 * Basis function object
 */
struct gkyl_basis {
  unsigned ndim, poly_order, num_basis, num_quad;
  char id[64]; // "serendipity", "tensor", "hybrid, "gkhybrid", "gkhybrid_vel"
  enum gkyl_basis_type b_type; // identifier for basis function
    
/**
 * Evaluate basis in unit cell (i.e. a hypercube with each side
 * [-1,1])
 *
 * @param z Location to evaluate basis. z \in [-1,1]^n
 * @param b On output, value of basis at 'z'
 */
  void (*eval)(const double *z, double *b);

/**
 * Evaluate expansion at point in the logical cell (hypercube)
 *
 * @param z Location to evaluate exansion. z \in [-1,1]^n
 * @param f Expansion coefficients
 * @return Expansion evaluated at z
 */
  double (*eval_expand)(const double *z, const double *f);

/**
 * Evaluate gradient, given expansion at point in the logical cell
 * (hypercube)
 *
 * @param dir Direction to compute gradient
 * @param z Location to evaluate exansion. z \in [-1,1]^n
 * @param f Expansion coefficients
 * @return Expansion evaluated at z
 */
  double (*eval_grad_expand)(int dir, const double *z, const double *f);

/**
 * Flip-sign function: changes signs of input expansion cofficients by
 * changing sign of odd monomial powers in specified direction. So if
 * dir=0, all odd powers of x appearing in the expansion will have
 * sign flipped.
 *
 * @param dir Direction to flip sign
 * @param f Input expansion
 * @param fout On output, flipped version of @a f
 */
  void (*flip_odd_sign)(int dir, const double *f, double *fout);

  /**
 * Flip-sign function: changes signs of input expansion cofficients by
 * changing sign of even monomial powers in specified direction. So if dir=0, all
 * even powers of x appearing in the expansion will have sign flipped.
 *
 * @param dir Direction to flip sign
 * @param f Input expansion
 * @param fout On output, flipped version of @a f
 */
  void (*flip_even_sign)(int dir, const double *f, double *fout);

/**
 * Construct list of nodes corresponding to this basis set. The nodes
 * coordinates are in the unit cell [-1,1]^ndim and stored such that
 * the coodinates of a node are contiguous, starting at index ndim*n,
 * n = 0, ... num_basis-1. The node_coords array must be pre-allocated
 * by the caller.
 */
  void (*node_list)(double *node_coords);

/**
 * Given expansion coefficients on nodal basis (nodes specified by the
 * node_list method), compute modal expansion coefficients.
 *
 * @param fnodal Coefficients of nodal expansion
 * @param fmodal On output, coefficients of modal expansion
 */
  void (*nodal_to_modal)(const double *fnodal, double *fmodal);

/**
 * Given expansion coefficients on nodes at the surface in one direction,
 * and Gauss-Legendre nodes in the other, compute modal expansion coefficients.
 *
 * @param fnodal Coefficients of nodal expansion
 * @param fmodal On output, coefficients of modal expansion
 */
  nodal_to_modal_quad_surf_t nodal_to_modal_quad_surf[3];

/**
 * Construct list of nodes that are on the surface in one direction
 * and on Gauss-Legendre coordinates in the other. The nodes
 * coordinates are in the unit cell [-1,1]^ndim and stored such that
 * the coodinates of a node are contiguous, starting at index ndim*n,
 * n = 0, ... num_basis-1. The node_coords array must be pre-allocated
 * by the caller.
 */
  node_quad_surf_list_t node_quad_surf_list[3];

/**
 * Given expansion coefficients on nodal basis defined by Gauss-Legendre
 * quadrature points, compute modal expansion coefficients.
 *
 * @param fquad Coefficients of nodal expansion in quadrature node basis
 * @param fmodal On output, coefficients of modal expansion
 * @param linc2 Modal component being updated (allows for parallelization over basis functions)
 */
  void (*quad_nodal_to_modal)(const double *fquad, double *fmodal, long linc2);  

/**
 * Given expansion coefficients of DG modal basis, evaluate basis at Gauss-Legendre
 * quadrature points of order p+1 (the Gauss-Legendre nodal basis).
 * 
 * @param fmodal Coefficients of modal expansion
 * @param fquad On output, evaluation of modal expansion at Gauss-Legendre quadrature nodal basis
 * @param linc2 Quadrature node being updated (allows for parallelization over quadrature points)
 */
  void (*modal_to_quad_nodal)(const double *fmodal, double *fquad, long linc2);   
}; 

/**
 * Assign object members in modal serendipity basis object.
 *
 * @param basis Basis object to initialize
 * @param ndim Dimension of reference element.
 * @param poly_order Polynomial order.
 */
void gkyl_cart_modal_serendip(struct gkyl_basis *basis, int ndim,
  int poly_order);
void gkyl_cart_modal_serendip_cu_dev(struct gkyl_basis *basis, int ndim,
  int poly_order);


/**
 * Create new modal serendipity basis function object.
 * This basis needs to be deallocated with free/release methods.
 *
 * @param ndim Dimension of reference element.
 * @param poly_order Polynomial order.
 * @return new basis struct.
 */
struct gkyl_basis * gkyl_cart_modal_serendip_new(int ndim, int poly_order);
struct gkyl_basis * gkyl_cart_modal_serendip_cu_dev_new(int ndim, int poly_order);

/**
 * Assign object members in modal tensor-product basis object.
 *
 * @param basis Basis object to initialize
 * @param ndim Dimension of reference element.
 * @param poly_order Polynomial order.
 */
void gkyl_cart_modal_tensor(struct gkyl_basis *basis, int ndim,
  int poly_order);
void gkyl_cart_modal_tensor_cu_dev(struct gkyl_basis *basis, int ndim,
  int poly_order);

/**
 * Create new modal tensor-product basis function object.
 * This basis needs to be deallocated with free/release methods.
 *
 * @param ndim Dimension of reference element.
 * @param poly_order Polynomial order.
 * @return new basis struct.
 */
struct gkyl_basis * gkyl_cart_modal_tensor_new(int ndim, int poly_order);
struct gkyl_basis * gkyl_cart_modal_tensor_cu_dev_new(int ndim, int poly_order);

/**
 * Assign object members in hybrid basis. These are p=1 in configuration space
 * and p=2 in velocity space.
 *
 * @param basis Basis object to initialize
 * @param cdim dimension of configuration space.
 * @param vdim dimension of velocity space.
 */
void gkyl_cart_modal_hybrid(struct gkyl_basis *basis, int cdim, int vdim);
void gkyl_cart_modal_hybrid_cu_dev(struct gkyl_basis *basis, int cdim, int vdim);

/**
 * Create new hybrid basis. These are p=1 in configuration space
 * and p=2 in velocity space.
 * This basis needs to be deallocated with free/release methods.
 *
 * @param cdim dimension of configuration space.
 * @param vdim dimension of velocity space.
 * @return new basis struct.
 */
struct gkyl_basis * gkyl_cart_modal_hybrid_new(int cdim, int vdim);
struct gkyl_basis * gkyl_cart_modal_hybrid_cu_dev_new(int cdim, int vdim);

/**
 * Assign object members in hybrid basis for use in gyrokinetics p=1
 * simulations. These basis have the v_par^2 monomial and its tensor
 * product with other monomials included. 
 *
 * @param basis Basis object to initialize
 * @param cdim dimension of configuration space.
 * @param vdim dimension of velocity space.
 */
void gkyl_cart_modal_gkhybrid(struct gkyl_basis *basis, int cdim, int vdim);
void gkyl_cart_modal_gkhybrid_cu_dev(struct gkyl_basis *basis, int cdim, int vdim);

/**
 * Create new hybrid basis for use in gyrokinetics p=1
 * simulations. These basis have the v_par^2 monomial and its tensor
 * product with other monomials included. 
 * This basis needs to be deallocated with free/release methods.
 *
 * @param cdim dimension of configuration space.
 * @param vdim dimension of velocity space.
 * @return new basis struct.
 */
struct gkyl_basis * gkyl_cart_modal_gkhybrid_new(int cdim, int vdim);
struct gkyl_basis * gkyl_cart_modal_gkhybrid_cu_dev_new(int cdim, int vdim);

/**
 * Assign object members in hybrid basis for use in gyrokinetics p=1
 * velocity space. These basis have the v_par^2 monomial and its tensor
 * product with mu monomials (for vdim=2). 
 *
 * @param basis Basis object to initialize
 * @param cdim dimension of configuration space.
 * @param vdim dimension of velocity space.
 */
void gkyl_cart_modal_gkhybrid_vel(struct gkyl_basis *basis, int vdim);
void gkyl_cart_modal_gkhybrid_vel_cu_dev(struct gkyl_basis *basis, int vdim);

/**
 * Create new hybrid basis for use in gyrokinetics p=1 velocity space.
 * These basis have the v_par^2 monomial and its tensor
 * product with mu monomials (for vdim=2). 
 * This basis needs to be deallocated with free/release methods.
 *
 * @param cdim dimension of configuration space.
 * @param vdim dimension of velocity space.
 * @return new basis struct.
 */
struct gkyl_basis * gkyl_cart_modal_gkhybrid_vel_new(int vdim);
struct gkyl_basis * gkyl_cart_modal_gkhybrid_vel_cu_dev_new(int vdim);

/**
 * Get the dimensionality of a basis (for opaque pointers).
 *
 * @param basis Basis object. 
 * @return Dimensionality of the basis.
 */
unsigned gkyl_cart_modal_basis_get_ndim(const struct gkyl_basis *basis);

/**
 * Get the polynomial order of a basis (for opaque pointers).
 *
 * @param basis Basis object. 
 * @return Polynomial order of the basis.
 */
unsigned gkyl_cart_modal_basis_get_poly_order(const struct gkyl_basis *basis);

/**
 * Get the num_basis of a basis (for opaque pointers).
 *
 * @param basis Basis object. 
 * @return Number of basis monomials.
 */
unsigned gkyl_cart_modal_basis_get_num_basis(const struct gkyl_basis *basis);

/**
 * Get the type of a basis (for opaque pointers).
 *
 * @param basis Basis object. 
 * @return Type of the basis.
 */
enum gkyl_basis_type gkyl_cart_modal_basis_get_type(const struct gkyl_basis *basis);

/**
 * Free the memory associated with basis objects
 * created with _new methods above.
 *
 * @param basis Basis object to free.
 */
void gkyl_cart_modal_basis_release(struct gkyl_basis *basis);
void gkyl_cart_modal_basis_release_cu(struct gkyl_basis *basis);
