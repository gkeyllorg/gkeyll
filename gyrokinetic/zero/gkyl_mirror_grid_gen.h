#pragma once

#include <gkyl_array.h>
#include <gkyl_math.h>
#include <gkyl_rect_grid.h>
#include <gkyl_range.h>
#include <gkyl_position_map.h>

// Forward declare internal object
struct gkyl_mirror_grid_gen_x;

// Geometric quantities: various vectors are stored as contravariant
// components
struct gkyl_mirror_grid_gen_geom {
  struct gkyl_vec3 tang[3]; // tangent vectors, e_i
  struct gkyl_vec3 dual[3]; // dual vectors, e^i
  struct gkyl_vec3 B; // Magnetic field
  double Jc; // Jacobian = e_1*(e_2 X e_3)  = 1/e^1*(e^2 X e^3)
  struct gkyl_vec3 curlbhat; // \nabla X \hat{b}
};

struct gkyl_mirror_grid_gen {
  struct gkyl_array *nodes_rza; // r,z,phi coordinates of corner nodes of cells
  struct gkyl_array *nodes_psi; // psi values at nodes
  struct gkyl_array *nodes_geom; // geometric quantities at nodes:
  // this is an array of gkyl_mirror_grid_gen_geom objects

  struct gkyl_mirror_grid_gen_x *gg_x; // Internal object that contains flags
  // describing the type of geometry being generated
};

// flag to indicate what field-line coordinate to use
enum gkyl_mirror_grid_gen_field_line_coord {
  GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z, // use psi and Cartesian Z coordinate
  GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z // use sqrt(psi) and Cartesian Z coordinate
};

// input struct to construct the mirror geometry
struct gkyl_mirror_grid_gen_inp {
  const struct gkyl_rect_grid *comp_grid; // Computational space grid (psi, phi, z)
  const struct gkyl_range nrange; // Nodal range
  const struct gkyl_range local; // local range
  const struct gkyl_range global; //global range
  int dir; // Surface direction
  const struct gkyl_position_map *position_map; // position map
  enum gkyl_mirror_grid_gen_field_line_coord fl_coord; // field-line coordinate to use
  bool include_axis; // add nodes on r=0 axis (the axis is assumed be psi=0)

  double R[2], Z[2]; // extents of R,Z grid on which psi(R,Z) is given
  int nrcells, nzcells; // number of cells in R and Z
  const struct gkyl_array *psiRZ; // nodal values of psi(R,Z)

  bool write_psi_cubic; // set to true to write the cubic fit to file
  const char *psi_cubic_fname; // name for cubic fit file
};

/**
 * Create new grid generator for mirror geometries. Return a NULL
 * pointer if the grid generator failed.
 *
 * @param inp Input parameters to the grid generator
 * @return newly create mirror geometry object
 */
struct gkyl_mirror_grid_gen *gkyl_mirror_grid_gen_inew(const struct gkyl_mirror_grid_gen_inp *inp);

struct gkyl_mirror_grid_gen *
gkyl_mirror_grid_gen_int_inew(const struct gkyl_mirror_grid_gen_inp *inp);

struct gkyl_mirror_grid_gen *
gkyl_mirror_grid_gen_surf_inew(const struct gkyl_mirror_grid_gen_inp *inp);

/**
 * Does the grid include the axis?
 *
 * @param geom Geometry object
 * @return true if axis is included, false otherwise.
 */
bool gkyl_mirror_grid_gen_is_include_axis(const struct gkyl_mirror_grid_gen *geom);

/**
 * Get field-line coordinate used in grid
 *
 * @param geom Geometry object
 * @return field-line coordinate used in grid
 */
enum gkyl_mirror_grid_gen_field_line_coord
gkyl_mirror_grid_gen_fl_coord(const struct gkyl_mirror_grid_gen *geom);

/**
 * Release the mirror grid object.
 *
 * @param geom Mirror grid object to release
 */
void gkyl_mirror_grid_gen_release(struct gkyl_mirror_grid_gen *geom);
