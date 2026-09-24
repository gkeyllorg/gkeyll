#pragma once

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_gr_spacetime.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct gr_neutronstar {
  struct gkyl_gr_spacetime spacetime; // Base spacetime object.

  double mass; // Mass of the neutron star.
  double spin; // Spin of the neutron star.

  double mass_quadrupole; // Mass quadrupole of the neutron star.
  double spin_octupole; // Spin octupole of the neutron star.
  double mass_hexadecapole; // Mass hexadecapole of the neutron star.

  double pos_x; // Position of the neutron star (x-direction).
  double pos_y; // Position of the neutron star (y-direction).
  double pos_z; // Position of the neutron star (z-direction).
};

// Input context, packaged as a struct.
struct gkyl_gr_neutronstar_inp {
  bool use_gpu; // Whether the spacetime object is on the host (false) or the device (true).

  double mass; // Mass of the neutron star.
  double spin; // Spin of the neutron star.

  double mass_quadrupole; // Mass quadrupole of the neutron star.
  double spin_octupole; // Spin octupole of the neutron star.
  double mass_hexadecapole; // Mass hexadecapole of the neutron star.

  double pos_x; // Position of the neutron star (x-direction).
  double pos_y; // Position of the neutron star (y-direction).
  double pos_z; // Position of the neutron star (z-direction).
};

/**
* Compute the scalar quantity A appearing within the Weyl-Lewis-Papapetrou metric functions, at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @return The Weyl-Lewis-Papapetrou scalar A.
*/
double
neutronstar_A_scalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z);

/**
* Compute the scalar quantity B appearing within the Weyl-Lewis-Papapetrou metric functions, at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @return The Weyl-Lewis-Papapetrou scalar B.
*/
double
neutronstar_B_scalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z);

/**
* Compute the scalar quantity H appearing within the Weyl-Lewis-Papapetrou metric functions, at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @return The Weyl-Lewis-Papapetrou scalar H.
*/
double
neutronstar_H_scalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z);

/**
* Compute the scalar quantity G appearing within the Weyl-Lewis-Papapetrou metric functions, at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @return The Weyl-Lewis-Papapetrou scalar G.
*/
double
neutronstar_G_scalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z);

/**
* Compute the scalar quantity F appearing within the Weyl-Lewis-Papapetrou metric functions, at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @return The Weyl-Lewis-Papapetrou scalar F.
*/
double
neutronstar_F_scalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z);

/**
* Compute the metric function f appearing within the Weyl-Lewis-Papapetrou metric, at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @return The Weyl-Lewis-Papapetrou function f.
*/
double
neutronstar_f_function(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z);

/**
* Compute the metric function omega appearing within the Weyl-Lewis-Papapetrou metric, at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @return The Weyl-Lewis-Papapetrou function omega.
*/
double
neutronstar_omega_function(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z);

/**
* Compute the metric function gamma appearing within the Weyl-Lewis-Papapetrou metric, at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @return The Weyl-Lewis-Papapetrou function gamma.
*/
double
neutronstar_gamma_function(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z);

/**
* Compute the spatial transformation tensor from Weyl-Lewis-Papapetrou cylindrical coordinates to Kerr-Schild Cartesian coordinates, at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @return The rank-2 (spatial) transformation tensor from Weyl-Lewis-Papapetrou to Kerr-Schild coordinates.
*/
double**
neutronstar_spatial_transformation_tensor(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z);

/**
* Compute the spacetime transformation tensor from Weyl-Lewis-Papapetrou cylindrical coordinates to Kerr-Schild Cartesian coordinates, at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @return The rank-2 (spacetime) transformation tensor from Weyl-Lewis-Papapetrou to Kerr-Schild coordinates.
*/
double**
neutronstar_spacetime_transformation_tensor(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z);

/**
* Compute the (inverse) spatial transformation tensor from Kerr-Schild Cartesian coordinates to Weyl-Lewis-Papapetrou cylindrical coordinates, at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @return The rank-2 (spatial, inverse) transformation tensor from Kerr-Schild to Weyl-Lewis-Papapetrou coordinates.
*/
double**
neutronstar_spatial_inv_transformation_tensor(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z);

/**
* Compute the (inverse) spacetime transformation tensor from Kerr-Schild Cartesian coordinates to Weyl-Lewis-Papapetrou cylindrical coordinates, at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @return The rank-2 (spacetime, inverse) transformation tensor from Kerr-Schild to Weyl-Lewis-Papapetrou coordinates.
*/
double**
neutronstar_spacetime_inv_transformation_tensor(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z);

/**
* Compute the rank-2 spatial metric tensor at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param spatial_metric_tensor Rank-2 spatial metric tensor (output).
*/
GKYL_CU_D
static void
neutronstar_spatial_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_metric_tensor);

/**
* Compute the rank-2 spacetime metric tensor at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordiante (y-direction).
* @param z Spatial coordinate (z-direction).
* @param spacetime_metric_tensor Rank-2 spacetime metric tensor (output).
*/
GKYL_CU_D
static void
neutronstar_spacetime_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_metric_tensor);

/**
* Compute the rank-2 inverse spatial metric tensor at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param spatial_inv_metric_tensor Rank-2 inverse spatial metric tensor (output).
*/
GKYL_CU_D
static void
neutronstar_spatial_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_inv_metric_tensor);

/**
* Compute the rank-2 inverse spacetime metric tensor at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param spacetime_inv_metric_tensor Rank-2 inverse spacetime metric tensor (output).
*/
GKYL_CU_D
static void
neutronstar_spacetime_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_inv_metric_tensor);

/**
* Compute the (scalar) spatial metric determinant at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param spatial_metric_det Spatial metric determinant (output).
*/
GKYL_CU_D
static void
neutronstar_spatial_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spatial_metric_det);

/**
* Compute the (scalar) spacetime metric determinant at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param spacetime_metric_det Spacetime metric determinant (output).
*/
GKYL_CU_D
static void
neutronstar_spacetime_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spacetime_metric_det);

/**
* Compute the rank-3 (spatial) partial derivative of the spatial metric tensor at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param spatial_metric_tensor_der Rank-3 partial derivative of the spatial metric tensor (output).
*/
GKYL_CU_D
static void
neutronstar_spatial_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_metric_tensor_der);

/**
* Compute the rank-3 (spacetime) partial derivative of the spacetime metric tensor at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dt Time coordinate spacing.
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param spacetime_metric_tensor_der Rank-3 partial derivative of the spacetime metric tensor (output).
*/
GKYL_CU_D
static void
neutronstar_spacetime_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_metric_tensor_der);

/**
* Compute the (scalar) lapse function at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param lapse_function Lapse function (output).
*/
GKYL_CU_D
static void
neutronstar_lapse_function(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* lapse_function);

/**
* Compute the rank-1 shift vector at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param shift_vector Rank-1 shift vector (output).
*/
GKYL_CU_D
static void
neutronstar_shift_vector(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double** shift_vector);

/**
* Compute the rank-1 (spatial) partial derivative of the lapse function at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param lapse_function_der Rank-1 partial derivative of the lapse function (output).
*/
GKYL_CU_D
static void
neutronstar_lapse_function_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** lapse_function_der);

/**
* Compute the rank-2 (spatial) partial derivative of the shift vector at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param shift_vector_der Rank-2 partial derivative of the shift vector (output).
*/
GKYL_CU_D
static void
neutronstar_shift_vector_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** shift_vector_der);

/**
* Compute the rank-3 (spatial) Christoffel symbols at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param spatial_christoffel Rank-3 spatial Christoffel symbols (output).
*/
GKYL_CU_D
static void
neutronstar_spatial_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_christoffel);

/**
* Compute the rank-3 (spacetime) Christoffel symbols at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dt Time coordinate spacing.
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param spacetime_christoffel Rank-3 spacetime Christoffel symbols (output).
*/
GKYL_CU_D
static void
neutronstar_spacetime_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_christoffel);

/**
* Compute the rank-4 (spatial) Riemann curvature tensor at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param spatial_riemann_tensor Rank-4 spatial Riemann curvature tensor (output).
*/
GKYL_CU_D
static void
neutronstar_spatial_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_riemann_tensor);

/**
* Compute the rank-4 (spacetime) Riemann curvature tensor at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dt Time coordinate spacing.
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param spacetime_riemann_tensor Rank-4 spacetime Riemann curvature tensor (output).
*/
GKYL_CU_D
static void
neutronstar_spacetime_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_riemann_tensor);

/**
* Compute the rank-2 (spatial) Ricci curvature tensor at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param spatial_ricci_tensor Rank-2 spatial Ricci curvature tensor (output).
*/
GKYL_CU_D
static void
neutronstar_spatial_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** spatial_ricci_tensor);

/**
* Compute the rank-2 (spacetime) Ricci curvature tensor at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dt Time coordinate spacing.
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param spacetime_ricci_tensor Rank-2 spacetime Ricci curvature tensor (output).
*/
GKYL_CU_D
static void
neutronstar_spacetime_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double*** spacetime_ricci_tensor);

/**
* Compute the (spatial) Ricci scalar curvature at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param spatial_ricci_scalar Spatial Ricci scalar curvature (output).
*/
GKYL_CU_D
static void
neutronstar_spatial_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double* spatial_ricci_scalar);

/**
* Compute the (spacetime) Ricci scalar curvature at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dt Time coordinate spacing.
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param spacetime_ricci_tensor Spacetime Ricci scalar curvature (output).
*/
GKYL_CU_D
static void
neutronstar_spacetime_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double* spacetime_ricci_scalar);

/**
* Compute the rank-4 (spatial) Weyl curvature tensor at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param spatial_weyl_tensor Rank-4 spatial Weyl curvature tensor (output).
*/
GKYL_CU_D
static void
neutronstar_spatial_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_weyl_tensor);

/**
* Compute the rank-4 (spacetime) Weyl curvature tensor at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dt Time coordinate spacing.
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param spacetime_weyl_tensor Rank-4 spacetime Weyl curvature tensor (output).
*/
GKYL_CU_D
static void
neutronstar_spacetime_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_weyl_tensor);

/**
* Compute the rank-2 extrinsic curvature tensor at a given point in a neutron star spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param extrinsic_curvature_tensor Rank-2 extrinsic curvature tensor (output).
*/
GKYL_CU_D
static void
neutronstar_extrinsic_curvature_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** extrinsic_curvature_tensor);

/**
* Determine whether a given point in a neutron star spacetime lies within an excision region.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param in_excision_region Whether the point lies in an excision region (output).
*/
GKYL_CU_D
static void
neutronstar_excision_region(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* in_excision_region);

/**
* Free neutron star spacetime object.
*
* @param ref Reference counter for neutron star spacetime.
*/
void
gkyl_gr_neutronstar_free(const struct gkyl_ref_count* ref);

/**
* Create a new neutron star spacetime object.
*
* @param use_gpu Whether the spacetime object is on the host (false) or the device (true).
* @return Pointer to the neutron star spacetime object.
*/
struct gkyl_gr_spacetime*
gkyl_gr_neutronstar_new(bool use_gpu, double mass, double spin, double mass_quadrupole, double spin_octupole, double mass_hexadecapole, double pos_x, double pos_y, double pos_z);

/**
* Create a new neutron star spacetime object, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the neutron star spacetime object.
*/
struct gkyl_gr_spacetime*
gkyl_gr_neutronstar_inew(const struct gkyl_gr_neutronstar_inp* inp);