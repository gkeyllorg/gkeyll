#include <gkyl_gr_spacetime.h>
#include <gkyl_alloc_flags_priv.h>

// These ensure that inline functions are only defined once.

static inline void
gkyl_gr_spatial_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_metric_tensor)
{
  return spacetime->spatial_metric_tensor_func(spacetime, t, x, y, z, spatial_metric_tensor);
}

static inline void
gkyl_gr_spacetime_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_metric_tensor)
{
  return spacetime->spacetime_metric_tensor_func(spacetime, t, x, y, z, spacetime_metric_tensor);
}

static inline void
gkyl_gr_spatial_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_inv_metric_tensor)
{
  return spacetime->spatial_inv_metric_tensor_func(spacetime, t, x, y, z, spatial_inv_metric_tensor);
}

static inline void
gkyl_gr_spacetime_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_inv_metric_tensor)
{
  return spacetime->spacetime_inv_metric_tensor_func(spacetime, t, x, y, z, spacetime_inv_metric_tensor);
}

static inline void
gkyl_gr_spatial_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spatial_metric_det)
{
  return spacetime->spatial_metric_det_func(spacetime, t, x, y, z, spatial_metric_det);
}

static inline void
gkyl_gr_spacetime_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spacetime_metric_det)
{
  return spacetime->spacetime_metric_det_func(spacetime, t, x, y, z, spacetime_metric_det);
}

static inline void
gkyl_gr_spatial_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_metric_tensor_der)
{
  return spacetime->spatial_metric_tensor_der_func(spacetime, t, x, y, z, dx, dy, dz, spatial_metric_tensor_der);
}

static inline void
gkyl_gr_spacetime_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_metric_tensor_der)
{
  return spacetime->spacetime_metric_tensor_der_func(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_metric_tensor_der);
}

static inline void
gkyl_gr_lapse_function(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* lapse_function)
{
  return spacetime->lapse_function_func(spacetime, t, x, y, z, lapse_function);
}

static inline void
gkyl_gr_shift_vector(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double** shift_vector)
{
  return spacetime->shift_vector_func(spacetime, t, x, y, z, shift_vector);
}

static inline void
gkyl_gr_lapse_function_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** lapse_function_der)
{
  return spacetime->lapse_function_der_func(spacetime, t, x, y, z, dx, dy, dz, lapse_function_der);
}

static inline void
gkyl_gr_shift_vector_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** shift_vector_der)
{
  return spacetime->shift_vector_der_func(spacetime, t, x, y, z, dx, dy, dz, shift_vector_der);
}

static inline void
gkyl_gr_spatial_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_christoffel)
{
  return spacetime->spatial_christoffel_func(spacetime, t, x, y, z, dx, dy, dz, spatial_christoffel);
}

static inline void
gkyl_gr_spacetime_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_christoffel)
{
  return spacetime->spacetime_christoffel_func(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_christoffel);
}

static inline void
gkyl_gr_spatial_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_riemann_tensor)
{
  return spacetime->spatial_riemann_tensor_func(spacetime, t, x, y, z, dx, dy, dz, spatial_riemann_tensor);
}

static inline void
gkyl_gr_spacetime_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_riemann_tensor)
{
  return spacetime->spacetime_riemann_tensor_func(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_riemann_tensor);
}

static inline void
gkyl_gr_spatial_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** spatial_ricci_tensor)
{
  return spacetime->spatial_ricci_tensor_func(spacetime, t, x, y, z, dx, dy, dz, spatial_ricci_tensor);
}

static inline void
gkyl_gr_spacetime_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double*** spacetime_ricci_tensor)
{
  return spacetime->spacetime_ricci_tensor_func(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_ricci_tensor);
}

static inline void
gkyl_gr_spatial_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double* spatial_ricci_scalar)
{
  return spacetime->spatial_ricci_scalar_func(spacetime, t, x, y, z, dx, dy, dz, spatial_ricci_scalar);
}

static inline void
gkyl_gr_spacetime_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double* spacetime_ricci_scalar)
{
  return spacetime->spacetime_ricci_scalar_func(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_ricci_scalar);
}

static inline void
gkyl_gr_spatial_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_weyl_tensor)
{
  return spacetime->spatial_weyl_tensor_func(spacetime, t, x, y, z, dx, dy, dz, spatial_weyl_tensor);
}

static inline void
gkyl_gr_spacetime_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_weyl_tensor)
{
  return spacetime->spacetime_weyl_tensor_func(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_weyl_tensor);
}

static inline void
gkyl_gr_extrinsic_curvature_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** extrinsic_curvature_tensor)
{
  return spacetime->extrinsic_curvature_tensor_func(spacetime, t, x, y, z, dx, dy, dz, extrinsic_curvature_tensor);
}

static inline void
gkyl_gr_conformal_factor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* conformal_factor)
{
  return spacetime->conformal_factor_func(spacetime, t, x, y, z, conformal_factor);
}

static inline void
gkyl_gr_bssn_conformal_factor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* bssn_conformal_factor)
{
  return spacetime->bssn_conformal_factor_func(spacetime, t, x, y, z, bssn_conformal_factor);
}

static inline void
gkyl_gr_conformal_factor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** conformal_factor_der)
{
  return spacetime->conformal_factor_der_func(spacetime, t, x, y, z, dx, dy, dz, conformal_factor_der);
}

static inline void
gkyl_gr_bssn_conformal_factor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** bssn_conformal_factor_der)
{
  return spacetime->bssn_conformal_factor_der_func(spacetime, t, x, y, z, dx, dy, dz, bssn_conformal_factor_der);
}

static inline void
gkyl_gr_conformal_factor_der2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** conformal_factor_der2)
{
  return spacetime->conformal_factor_der2_func(spacetime, t, x, y, z, dx, dy, dz, conformal_factor_der2);
}

static inline void
gkyl_gr_bssn_conformal_factor_der2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** bssn_conformal_factor_der2)
{
  return spacetime->bssn_conformal_factor_der2_func(spacetime, t, x, y, z, dx, dy, dz, bssn_conformal_factor_der2);
}

static inline void
gkyl_gr_excision_region(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* in_excision_region)
{
  return spacetime->excision_region_func(spacetime, t, x, y, z, in_excision_region);
}

bool
gkyl_gr_spacetime_is_cu_dev(const struct gkyl_gr_spacetime* spacetime)
{
  return GKYL_IS_CU_ALLOC(spacetime->flags);
}

struct gkyl_gr_spacetime*
gkyl_gr_spacetime_acquire(const struct gkyl_gr_spacetime* spacetime)
{
  gkyl_ref_count_inc(&spacetime->ref_count);

  return (struct gkyl_gr_spacetime*) spacetime;
}

void
gkyl_gr_spacetime_release(const struct gkyl_gr_spacetime* spacetime)
{
  gkyl_ref_count_dec(&spacetime->ref_count);
}