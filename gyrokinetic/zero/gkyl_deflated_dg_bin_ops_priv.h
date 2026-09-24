#pragma once
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_basis.h>

// Struct containing all fields 
// needed to do the div or mul on a
// surface which is constant in the last coordinate
// of the configuration space grid
struct deflated_dg_bin_ops_data {
  struct gkyl_array *deflated_lop;
  struct gkyl_array *deflated_rop;
  struct gkyl_array *deflated_out;
};

enum deflated_dg_bin_ops_type { 
  GKYL_DEFLATED_MUL = 0, 
  GKYL_DEFLATED_DIV = 1, 
};


// Updater type
struct gkyl_deflated_dg_bin_ops {
  struct gkyl_rect_grid grid; // Conf space grid
  struct gkyl_rect_grid deflated_grid; // Conf space grid with last
                                       // dimension removed
  struct gkyl_basis basis; // Basis object
  struct gkyl_basis *basis_on_dev; // Basis object on GPU
  struct gkyl_basis deflated_basis; // Basis object with one
                                    // less dimension than basis
  struct gkyl_basis *deflated_basis_on_dev; // deflated basis on GPU
  struct gkyl_range local; // Local range where operation should hapen
  struct gkyl_range deflated_local; // Local range with last dim removed
  struct gkyl_range deflated_local_ext; // extended deflated local range
  struct gkyl_range nrange; // nodal range corresponding to
                            // local range
  struct gkyl_range deflated_nrange; // nodal range corresponding to
                                     // deflated local range

  struct gkyl_dg_bin_op_mem *mem; // memory fir division of deflated fields

  struct deflated_dg_bin_ops_data *d_bop_data; // Array of deflated_dem_data
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
