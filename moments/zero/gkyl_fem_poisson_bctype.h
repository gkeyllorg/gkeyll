#pragma once

#include <gkyl_util.h>
#include <math.h>

// Boundary condition types.
enum gkyl_poisson_bc_type {
  GKYL_POISSON_PERIODIC = 0,
  GKYL_POISSON_DIRICHLET, // Sets the value.
  GKYL_POISSON_NEUMANN, // Sets the slope normal to the boundary.
  GKYL_POISSON_ROBIN, // A combination of dirichlet and neumann.  
  GKYL_POISSON_DIRICHLET_VARYING, // Sets the value, spatially varying.
};

// Boundary condition values. Dirichlet and Neumann use only one value,
// Robin uses 3, and periodic ignores the value.
struct gkyl_poisson_bc_value { double v[3]; };

struct gkyl_poisson_bias_plane {
  int dir; // Direction perpendicular to the plane.
  double loc; // Location of the plane in the 'dir' dimension.
  double val; // Biasing value.
};

struct gkyl_poisson_bias_plane_list {
  int num_bias_plane; // Number of bias planes.
  struct gkyl_poisson_bias_plane *bp;
};

struct gkyl_poisson_bias_line {
  // Directions perpendicular to the line and coordinates in those directions
  // (in 3D space; for 2x sims value in 3rd direction is ignored).
  int perp_dirs[GKYL_MAX_CDIM-1];
  double perp_coords[GKYL_MAX_CDIM-1];
  double val; // Biasing value.
};

struct gkyl_poisson_bias_line_list {
  int num_bias_line; // Number of bias planes.
  struct gkyl_poisson_bias_line *bl;
};

struct gkyl_poisson_bc {
  enum gkyl_poisson_bc_type lo_type[GKYL_MAX_CDIM], up_type[GKYL_MAX_CDIM];
  struct gkyl_poisson_bc_value lo_value[GKYL_MAX_CDIM], up_value[GKYL_MAX_CDIM];
  // Function specifying a spatially varying BC.
  void (*bc_value_func)(double t, const double *xn, double *phi_wall_up_out, void *ctx);
  void *bc_value_func_ctx;

  // Additional attributes to apply a bias plane at the extremal z values only.
  bool contains_lower_z_edge, contains_upper_z_edge;
};
