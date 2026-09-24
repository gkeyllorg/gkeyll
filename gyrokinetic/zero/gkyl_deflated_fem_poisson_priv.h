#pragma once
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_fem_poisson_bctype.h>

// Struct containing all fields and solver
// needed to solve the poisson equation on a
// surface which is constant in the last coordinate
// of the configuration space grid
struct deflated_fem_data {
  struct gkyl_array *deflated_rhs;
  struct gkyl_array *deflated_phibc;
  struct gkyl_array *deflated_phi;
  struct gkyl_array *deflated_epsilon;
  struct gkyl_array *deflated_kSq;
  struct gkyl_array *deflated_nodal_fld;
  struct gkyl_fem_poisson *fem_poisson;
};

// Updater type
struct gkyl_deflated_fem_poisson {
  bool isdirichletvar; // True if using a spatially varying Dirichlet BC.
  bool ishelmholtz; // if solving Helmholtz equation (kSq is not zero/NULL).

  struct gkyl_rect_grid grid; // Conf space grid
  struct gkyl_rect_grid deflated_grid; // Conf space grid with last
                                       // dimension removed
  struct gkyl_basis basis; // Basis object
  struct gkyl_basis *basis_on_dev; // Basis object on GPU
  struct gkyl_basis deflated_basis; // Basis object with one
                                    // less dimension than basis
  struct gkyl_basis *deflated_basis_on_dev; // deflated basis on GPU
  struct gkyl_range local; // Local range where poisson problem
                           // should be solved
  struct gkyl_range deflated_local; // Local range with last dim removed
  struct gkyl_range global_sub_range; // Local range as a
                                      // sub range of the global range
  struct gkyl_range deflated_local_ext; // extended deflated local range
  struct gkyl_range nrange; // nodal range corresponding to
                            // local range
  struct gkyl_range deflated_nrange; // nodal range corresponding to
                                     // deflated local range
  struct gkyl_poisson_bc poisson_bc; // Boundary conditions
  struct gkyl_poisson_bias_plane_list *bias_plane_list; // Biased planes
  struct deflated_fem_data *d_fem_data; // Array of deflated_dem_data
                                        // to be used for individual surface solves
  struct gkyl_array *nodal_fld; // Nodal field which holds solution
  int num_solves_z; // Number of surfaces to solve on
  int cdim; // Dimension of configureation space
  struct gkyl_nodal_ops *n2m; // Nodal to modal operator to be
                              // used to construct the final DG solution
  struct gkyl_nodal_ops *n2m_deflated; // Nodal to modal operator to
                                       // be used at each surface
  struct gkyl_deflate_zsurf *deflator_lo; // Deflation operator used to
                                          // evaluate fields at a lower surface
  struct gkyl_deflate_zsurf *deflator_up; // Deflation operator used to
                                          // evaluate fields at an upper surface
  bool use_gpu; // Whether to use the GPU
};
