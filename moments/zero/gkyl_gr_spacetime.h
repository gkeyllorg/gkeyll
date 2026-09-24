#pragma once

#include <gkyl_ref_count.h>
#include <gkyl_util.h>
#include <gkyl_evalf_def.h>

// Spacetime gauge choices.
enum gkyl_spacetime_gauge {
  GKYL_STATIC_GAUGE = 0, // Default.
  GKYL_BLACKHOLE_COLLAPSE_GAUGE
};

// Dynamic spacetime slicing condition.
enum gkyl_spacetime_slicing {
  GKYL_GEODESIC_SLICING = 0, // Default.
  GKYL_HARMONIC_SLICING,
  GKYL_1PLUSLOG_SLICING
};

// Dynamic spacetime evolution system.
enum gkyl_spacetime_evolution {
  GKYL_RICCI_EVOLUTION = 0, // Default.
  GKYL_EINSTEIN_EVOLUTION
};

// Forward declare the spacetime struct, for use in future function pointers.
struct gkyl_gr_spacetime;

// Function pointer to compute the rank-2 spatial metric tensor at a given point in spacetime.
typedef void (*gr_spatial_metric_tensor_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, double*** spatial_metric_tensor);

// Function pointer to compute the rank-2 spacetime metric tensor at a given point in spacetime.
typedef void (*gr_spacetime_metric_tensor_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, double*** spacetime_metric_tensor);

// Function pointer to compute the rank-2 inverse spatial metric tensor at a given point in spacetime.
typedef void (*gr_spatial_inv_metric_tensor_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, double*** spatial_inv_metric_tensor);

// Function pointer to compute the rank-2 inverse spacetime metric tensor at a given point in spacetime.
typedef void (*gr_spacetime_inv_metric_tensor_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, double*** spacetime_inv_metric_tensor);

// Function pointer to compute the (scalar) spatial metric determinant at a given point in spacetime.
typedef void (*gr_spatial_metric_det_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, double* spatial_metric_det);

// Function pointer to compute the (scalar) spacetime metric determinant at a given point in spacetime.
typedef void (*gr_spacetime_metric_det_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, double* spacetime_metric_det);

// Function pointer to compute the rank-3 (spatial) partial derivative of the spatial metric tensor at a given point in spacetime.
typedef void (*gr_spatial_metric_tensor_der_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dx, const double dy, const double dz, double**** spatial_metric_tensor_der);

// Function pointer to compute the rank-3 (spacetime) partial derivative of the spacetime metric tensor at a given point in spacetime.
typedef void (*gr_spacetime_metric_tensor_der_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dt, const double dx, const double dy, const double dz, double**** spacetime_metric_tensor_der);

// Function pointer to compute the (scalar) lapse function at a given point in spacetime.
typedef void (*gr_lapse_function_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, double* lapse_function);

// Function pointer to compute the rank-1 shift vector at a given point in spacetime.
typedef void (*gr_shift_vector_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, double** shift_vector);

// Function pointer to compute the rank-1 (spatial) partial derivative of the lapse function at a given point in spacetime.
typedef void (*gr_lapse_function_der_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dx, const double dy, const double dz, double** lapse_function_der);

// Function pointer to compute the rank-2 (spatial) partial derivative of the shift vector at a given point in spacetime.
typedef void (*gr_shift_vector_der_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dx, const double dy, const double dz, double*** shift_vector_der);

// Function pointer to compute the rank-3 (spatial) Christoffel symbols at a given point in spacetime.
typedef void (*gr_spatial_christoffel_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dx, const double dy, const double dz, double**** spatial_christoffel);

// Function pointer to compute the rank-3 (spacetime) Christoffel symbols at a given point in spacetime.
typedef void (*gr_spacetime_christoffel_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dt, const double dx, const double dy, const double dz, double**** spacetime_christoffel);

// Function pointer to compute the rank-4 (spatial) Riemann curvature tensor at a given point in spacetime.
typedef void (*gr_spatial_riemann_tensor_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dx, const double dy, const double dz, double***** spatial_riemann_tensor);

// Function pointer to compute the rank-4 (spacetime) Riemann curvature tensor at a given point in spacetime.
typedef void (*gr_spacetime_riemann_tensor_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dt, const double dx, const double dy, const double dz, double***** spacetime_riemann_tensor);

// Function pointer to compute the rank-2 (spatial) Ricci curvature tensor at a given point in spacetime.
typedef void (*gr_spatial_ricci_tensor_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dx, const double dy, const double dz, double*** spatial_ricci_tensor);

// Function pointer to compute the rank-2 (spacetime) Ricci curvature tensor at a given point in spacetime.
typedef void (*gr_spacetime_ricci_tensor_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dt, const double dx, const double dy, const double dz, double*** spacetime_ricci_tensor);

// Function pointer to compute the (spatial) Ricci scalar curvature at a given point in spacetime.
typedef void (*gr_spatial_ricci_scalar_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dx, const double dy, const double dz, double* spatial_ricci_scalar);

// Function pointer to compute the (spacetime) Ricci scalar curvature at a given point in spacetime.
typedef void (*gr_spacetime_ricci_scalar_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dt, const double dx, const double dy, const double dz, double* spacetime_ricci_scalar);

// Function pointer to compute the rank-4 (spatial) Weyl curvature tensor at a given point in spacetime.
typedef void (*gr_spatial_weyl_tensor_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dx, const double dy, const double dz, double***** spatial_weyl_tensor);

// Function pointer to compute the rank-4 (spacetime) Weyl curvature tensor at a given point in spacetime.
typedef void (*gr_spacetime_weyl_tensor_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dt, const double dx, const double dy, const double dz, double***** spacetime_weyl_tensor);

// Function pointer to compute the rank-2 extrinsic curvature tensor at a given point in spacetime.
typedef void (*gr_extrinsic_curvature_tensor_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dx, const double dy, const double dz, double*** extrinsic_curvature_tensor);

// Function pointer to compute the conformal factor at a given point in spacetime.
typedef void (*gr_conformal_factor_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, double* conformal_factor);

// Function pointer to compute the BSSN conformal factor at a given point in spacetime.
typedef void (*gr_bssn_conformal_factor_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, double* bssn_conformal_factor);

// Function pointer to compute the rank-1 (spatial) partial derivative of the conformal factor function at a given point in spacetime.
typedef void (*gr_conformal_factor_der_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dx, const double dy, const double dz, double** conformal_factor_der);

// Function pointer to compute the rank-1 (spatial) partial derivative of the BSSN conformal factor function at a given point in spacetime.
typedef void (*gr_bssn_conformal_factor_der_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dx, const double dy, const double dz, double** bssn_conformal_factor_der);

// Function pointer to compute the rank-2 (spatial) second partial derivative of the conformal factor function at a given point in spacetime.
typedef void (*gr_conformal_factor_der2_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dx, const double dy, const double dz, double*** conformal_factor_der2);

// Function pointer to compute the rank-2 (spatial) second partial derivative of the BSSN conformal factor function at a given point in spacetime.
typedef void (*gr_bssn_conformal_factor_der2_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, const double dx, const double dy, const double dz, double*** bssn_conformal_factor_der2);

// Function pointer to determine whether a given point in spacetime lies inside an excision region.
typedef void (*gr_excision_region_t)(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y,
  const double z, bool* in_excision_region);

struct gkyl_gr_spacetime {
  gr_spatial_metric_tensor_t spatial_metric_tensor_func; // Function to compute spatial metric tensor.
  gr_spacetime_metric_tensor_t spacetime_metric_tensor_func; // Function to compute spacetime metric tensor.

  gr_spatial_inv_metric_tensor_t spatial_inv_metric_tensor_func; // Function to compute inverse spatial metric tensor.
  gr_spacetime_inv_metric_tensor_t spacetime_inv_metric_tensor_func; // Function to compute inverse spacetime metric tensor.

  gr_spatial_metric_det_t spatial_metric_det_func; // Function to compute spatial metric determinant.
  gr_spacetime_metric_det_t spacetime_metric_det_func; // Function to compute spacetime metric determinant.

  gr_spatial_metric_tensor_der_t spatial_metric_tensor_der_func; // Function to compute partial derivative of spatial metric tensor.
  gr_spacetime_metric_tensor_der_t spacetime_metric_tensor_der_func; // Function to compute partial derivative of spacetime metric tensor.

  gr_lapse_function_t lapse_function_func; // Function to compute lapse function.
  gr_shift_vector_t shift_vector_func; // Function to compute shift vector.

  gr_lapse_function_der_t lapse_function_der_func; // Function to compute partial derivative of lapse function.
  gr_shift_vector_der_t shift_vector_der_func; // Function to compute partial derivative of shift vector.

  gr_spatial_christoffel_t spatial_christoffel_func; // Function to compute spatial Christoffel symbols.
  gr_spacetime_christoffel_t spacetime_christoffel_func; // Function to compute spacetime Christoffel symbols.

  gr_spatial_riemann_tensor_t spatial_riemann_tensor_func; // Function to compute spatial Riemann curvature tensor.
  gr_spacetime_riemann_tensor_t spacetime_riemann_tensor_func; // Function to compute spacetime Riemann curvature tensor.

  gr_spatial_ricci_tensor_t spatial_ricci_tensor_func; // Function to compute spatial Ricci curvature tensor.
  gr_spacetime_ricci_tensor_t spacetime_ricci_tensor_func; // Function to compute spacetime Ricci curvature tensor.

  gr_spatial_ricci_scalar_t spatial_ricci_scalar_func; // Function to compute spatial Ricci scalar curvature.
  gr_spacetime_ricci_scalar_t spacetime_ricci_scalar_func; // Function to compute spacetime Ricci scalar curvature.

  gr_spatial_weyl_tensor_t spatial_weyl_tensor_func; // Function to compute spatial Weyl curvature tensor.
  gr_spacetime_weyl_tensor_t spacetime_weyl_tensor_func; // Function to compute spacetime Weyl curvature tensor.

  gr_extrinsic_curvature_tensor_t extrinsic_curvature_tensor_func; // Function to compute extrinsic curvature tensor.

  gr_conformal_factor_t conformal_factor_func; // Function to compute conformal factor.
  gr_bssn_conformal_factor_t bssn_conformal_factor_func; // Function to compute BSSN conformal factor.

  gr_conformal_factor_der_t conformal_factor_der_func; // Function to compute partial derivative of conformal factor.
  gr_bssn_conformal_factor_der_t bssn_conformal_factor_der_func; // Function to compute partial derivative of BSSN conformal factor.

  gr_conformal_factor_der2_t conformal_factor_der2_func; // Function to compute second partial derivative of conformal factor.
  gr_bssn_conformal_factor_der2_t bssn_conformal_factor_der2_func; // Function to compute second partial derivative of BSSN conformal factor.

  gr_excision_region_t excision_region_func; // Function to determine whether point lies within excision region.

  uint32_t flags;
  struct gkyl_ref_count ref_count; // Reference count.
  struct gkyl_gr_spacetime *on_dev; // Pointer to itself, or device
};

/**
* Compute the rank-2 spatial metric tensor at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param spatial_metric_tensor Rank-2 spatial metric tensor (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_spatial_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_metric_tensor);

/**
* Compute the rank-2 spacetime metric tensor at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param spacetime_metric_tensor Rank-2 spacetime metric tensor (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_spacetime_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_metric_tensor);

/**
* Compute the rank-2 inverse spatial metric tensor at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param spatial_inv_metric_tensor Rank-2 inverse spatial metric tensor (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_spatial_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_inv_metric_tensor);

/**
* Compute the rank-2 inverse spacetime metric tensor at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param spacetime_inv_metric_tensor Rank-2 inverse spacetime metric tensor (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_spacetime_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_inv_metric_tensor);

/**
* Compute the (scalar) spatial metric determinant at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param spatial_metric_det Spatial metric determinant (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_spatial_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spatial_metric_det);

/**
* Compute the (scalar) spacetime metric determinant at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param spacetime_metric_det Spacetime metric determinant (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_spacetime_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spacetime_metric_det);

/**
* Compute the rank-3 (spatial) partial derivative of the spatial metric tensor at a given point in spacetime.
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
GKYL_CU_DH
static inline void
gkyl_gr_spatial_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_metric_tensor_der);

/**
* Compute the rank-3 (spacetime) partial derivative of the spacetime metric tensor at a given point in spacetime.
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
GKYL_CU_DH
static inline void
gkyl_gr_spacetime_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_metric_tensor_der);

/**
* Compute the (scalar) lapse function at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param lapse_function Lapse function (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_lapse_function(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* lapse_function);

/**
* Compute the rank-1 shift vector at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param shift_vector Rank-1 shift vector (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_shift_vector(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double** shift_vector);

/**
* Compute the rank-1 (spatial) partial derivative of the lapse function at a given point in spacetime.
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
GKYL_CU_DH
static inline void
gkyl_gr_lapse_function_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** lapse_function_der);

/**
* Compute the rank-2 (spatial) partial derivative of the shift vector at a given point in spacetime.
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
GKYL_CU_DH
static inline void
gkyl_gr_shift_vector_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** shift_vector_der);

/**
* Compute the rank-3 (spatial) Christoffel symbols at a given point in spacetime.
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
GKYL_CU_DH
static inline void
gkyl_gr_spatial_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_christoffel);

/**
* Compute the rank-3 (spacetime) Christoffel symbols at a given point in spacetime.
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
GKYL_CU_DH
static inline void
gkyl_gr_spacetime_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_christoffel);

/**
* Compute the rank-4 (spatial) Riemann curvature tensor at a given point in spacetime.
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
GKYL_CU_DH
static inline void
gkyl_gr_spatial_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_riemann_tensor);

/**
* Compute the rank-4 (spacetime) Riemann curvature tensor at a given point in spacetime.
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
GKYL_CU_DH
static inline void
gkyl_gr_spacetime_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_riemann_tensor);

/**
* Compute the rank-2 (spatial) Ricci curvature tensor at a given point in spacetime.
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
GKYL_CU_DH
static inline void
gkyl_gr_spatial_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** spatial_ricci_tensor);

/**
* Compute the rank-2 (spacetime) Ricci curvature tensor at a given point in spacetime.
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
* @param spacetime_ricci_tensor Rank-2 spatial Ricci curvature tensor (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_spacetime_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double*** spacetime_ricci_tensor);

/**
* Compute the (spatial) Ricci scalar curvature at a given point in spacetime.
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
GKYL_CU_DH
static inline void
gkyl_gr_spatial_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double* spatial_ricci_scalar);

/**
* Compute the (spacetime) Ricci scalar curvature at a given point in spacetime.
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
* @param spacetime_ricci_scalar Spatial Ricci scalar curvature (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_spacetime_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double* spacetime_ricci_scalar);

/**
* Compute the rank-4 (spatial) Weyl curvature tensor at a given point in spacetime.
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
GKYL_CU_DH
static inline void
gkyl_gr_spatial_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_weyl_tensor);

/**
* Compute the rank-4 (spacetime) Weyl curvature tensor at a given point in spacetime.
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
GKYL_CU_DH
static inline void
gkyl_gr_spacetime_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_weyl_tensor);

/**
* Compute the rank-2 extrinsic curvature tensor at a given point in spacetime.
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
GKYL_CU_DH
static inline void
gkyl_gr_extrinsic_curvature_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** extrinsic_curvature_tensor);

/**
* Compute the conformal factor at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param conformal_factor Conformal factor (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_conformal_factor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* conformal_factor);

/**
* Compute the BSSN conformal factor at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param conformal_factor BSSN conformal factor (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_bssn_conformal_factor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* bssn_conformal_factor);

/**
* Compute the rank-1 (spatial) partial derivative of the conformal factor at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param conformal_factor_der Rank-1 partial derivative of the conformal factor (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_conformal_factor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** conformal_factor_der);

/**
* Compute the rank-1 (spatial) partial derivative of the BSSN conformal factor at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param bssn_conformal_factor_der Rank-1 partial derivative of the BSSN conformal factor (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_bssn_conformal_factor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** bssn_conformal_factor_der);

/**
* Compute the rank-2 (spatial) second partial derivative of the conformal factor at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param conformal_factor_der2 Rank-2 second partial derivative of the conformal factor (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_conformal_factor_der2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** conformal_factor_der2);

/**
* Compute the rank-2 (spatial) second partial derivative of the BSSN conformal factor at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param bssn_conformal_factor_der2 Rank-2 second partial derivative of the BSSN conformal factor (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_bssn_conformal_factor_der2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** bssn_conformal_factor_der2);

/**
* Determine whether a given point in spacetime lies within an excision region.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param in_excision_region Whether the spacetime point lies in an excision region (output).
*/
GKYL_CU_DH
static inline void
gkyl_gr_excision_region(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* in_excision_region);

/**
* Check whether the spacetime is on device.
*
* @param spacetime Spacetime to check.
* @return Whether the spacetime is on device.
*/
bool gkyl_gr_spacetime_is_cu_dev(const struct gkyl_gr_spacetime* spacetime);

/**
* Acquire pointer to the spacetime object. Delete using the release() method.
*
* @param spacetime Spacetime object.
* @return Acquired spacetime object pointer.
*/
struct gkyl_gr_spacetime* gkyl_gr_spacetime_acquire(const struct gkyl_gr_spacetime* spacetime);

/**
* Delete spacetime object.
*
* @param spacetime Spacetime object to delete.
*/
void gkyl_gr_spacetime_release(const struct gkyl_gr_spacetime* spacetime);