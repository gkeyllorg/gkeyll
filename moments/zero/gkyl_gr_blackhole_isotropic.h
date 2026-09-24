#pragma once

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_gr_spacetime.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct gr_blackhole_isotropic {
  struct gkyl_gr_spacetime spacetime; // Base spacetime object.

  double mass; // Mass of the black hole in isotropic gauge.
  double spin; // Spin of the black hole in isotropic gauge.

  double pos_x; // Position of the black hole in isotropic gauge (x-direction).
  double pos_y; // Position of the black hole in isotropic gauge (y-direction).
  double pos_z; // Position of the black hole in isotropic gauge (z-direction).
};

// Input context, packaged as a struct.
struct gkyl_gr_blackhole_isotropic_inp {
  bool use_gpu; // Whether the spacetime object is on the host (false) or the device (true).

  double mass; // Mass of the black hole in isotropic gauge.
  double spin; // Spin of the black hole in isotropic gauge.

  double pos_x; // Position of the black hole in isotropic gauge (x-direction).
  double pos_y; // Position of the black hole in isotropic gauge (y-direction).
  double pos_z; // Position of the black hole in isotropic gauge (z-direction).
};

/**
* Compute the rank-2 spatial metric tensor at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spatial_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_metric_tensor);

/**
* Compute the rank-2 spacetime metric tensor at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spacetime_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_metric_tensor);

/**
* Compute the rank-2 inverse spatial metric tensor at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spatial_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_inv_metric_tensor);

/**
* Compute the rank-2 inverse spacetime metric tensor at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spacetime_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_inv_metric_tensor);

/**
* Compute the (scalar) spatial metric determinant at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spatial_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spatial_metric_det);

/**
* Compute the (scalar) spacetime metric determinant at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spacetime_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spacetime_metric_det);

/**
* Compute the rank-3 (spatial) partial derivative of the spatial metric tensor at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spatial_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_metric_tensor_der);

/**
* Compute the rank-3 (spacetime) partial derivative of the spacetime metric tensor at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spacetime_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_metric_tensor_der);

/**
* Compute the (scalar) lapse function at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_lapse_function(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* lapse_function);

/**
* Compute the rank-1 shift vector at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_shift_vector(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double** shift_vector);

/**
* Compute the rank-1 (spatial) partial derivative of the lapse function at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_lapse_function_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** lapse_function_der);

/**
* Compute the rank-2 (spatial) partial derivative of the shift vector at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_shift_vector_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** shift_vector_der);

/**
* Compute the rank-3 (spatial) Christoffel symbols at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spatial_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_christoffel);

/**
* Compute the rank-3 (spacetime) Christoffel symbols at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spacetime_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_christoffel);

/**
* Compute the rank-4 (spatial) Riemann curvature tensor at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spatial_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_riemann_tensor);

/**
* Compute the rank-4 (spacetime) Riemann curvature tensor at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spacetime_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_riemann_tensor);

/**
* Compute the rank-2 (spatial) Ricci curvature tensor at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spatial_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** spatial_ricci_tensor);

/**
* Compute the rank-2 (spacetime) Ricci curvature tensor at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spacetime_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double*** spacetime_ricci_tensor);

/**
* Compute the (spatial) Ricci scalar curvature at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spatial_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double* spatial_ricci_scalar);

/**
* Compute the (spacetime) Ricci scalar curvature at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spacetime_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double* spacetime_ricci_scalar);

/**
* Compute the rank-4 (spatial) Weyl curvature tensor at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spatial_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_weyl_tensor);

/**
* Compute the rank-4 (spacetime) Weyl curvature tensor at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_spacetime_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_weyl_tensor);

/**
* Compute the rank-2 extrinsic curvature tensor at a given point in a black hole spacetime in isotropic gauge.
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
blackhole_isotropic_extrinsic_curvature_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** extrinsic_curvature_tensor);

/**
* Determine whether a given point in a black hole spacetime in isotropic gauge lies within an excision region.
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
blackhole_isotropic_excision_region(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* in_excision_region);

/**
* Free black hole spacetime in isotropic gauge object.
*
* @param ref Reference counter for black hole spacetime in isotropic gauge.
*/
void
gkyl_gr_blackhole_isotropic_free(const struct gkyl_ref_count* ref);

/**
* Create a new black hole spacetime in isotropic gauge object.
*
* @param use_gpu Whether the spacetime object is on the host (false) or the device (true).
* @return Pointer to the black hole spacetime in isotropic gauge object.
*/
struct gkyl_gr_spacetime*
gkyl_gr_blackhole_isotropic_new(bool use_gpu, double mass, double spin, double pos_x, double pos_y, double pos_z);

/**
* Create a new black hole spacetime in isotropic gauge object, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the black hole spacetime in isotropic gauge object.
*/
struct gkyl_gr_spacetime*
gkyl_gr_blackhole_isotropic_inew(const struct gkyl_gr_blackhole_isotropic_inp* inp);