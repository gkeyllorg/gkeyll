#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gr_minkowski.h>

static void
minkowski_spatial_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_metric_tensor)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        (*spatial_metric_tensor)[i][j] = 1.0;
      }
      else {
        (*spatial_metric_tensor)[i][j] = 0.0;
      }
    }
  }
}

static void
minkowski_spacetime_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_metric_tensor)
{
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      if (i == j) {
        if (i == 0) {
          (*spacetime_metric_tensor)[i][j] = -1.0;
        }
        else {
          (*spacetime_metric_tensor)[i][j] = 1.0;
        }
      }
      else {
        (*spacetime_metric_tensor)[i][j] = 0.0;
      }
    }
  }
}

static void
minkowski_spatial_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_inv_metric_tensor)
{
  minkowski_spatial_metric_tensor(spacetime, t, x, y, z, spatial_inv_metric_tensor);
}

static void
minkowski_spacetime_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_inv_metric_tensor)
{
  minkowski_spacetime_metric_tensor(spacetime, t, x, y, z, spacetime_inv_metric_tensor);
}

static void
minkowski_spatial_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spatial_metric_det)
{
  *spatial_metric_det = 1.0;
}

static void
minkowski_spacetime_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spacetime_metric_det)
{
  *spacetime_metric_det = -1.0;
}

static void
minkowski_spatial_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
   const double dx, const double dy, const double dz, double**** spatial_metric_tensor_der)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        (*spatial_metric_tensor_der)[i][j][k] = 0.0;
      }
    }
  }
}

static void
minkowski_spacetime_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_metric_tensor_der)
{
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        (*spacetime_metric_tensor_der)[i][j][k] = 0.0;
      }
    }
  }
}

static void
minkowski_lapse_function(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* lapse_function)
{
  *lapse_function = 1.0;
}

static void
minkowski_shift_vector(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double** shift_vector)
{
  for (int i = 0; i < 3; i++) {
    (*shift_vector)[i] = 0.0;
  }
}

static void
minkowski_lapse_function_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** lapse_function_der)
{
  for (int i = 0; i < 3; i++) {
    (*lapse_function_der)[i] = 0.0;
  }
}

static void
minkowski_shift_vector_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** shift_vector_der)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*shift_vector_der)[i][j] = 0.0;
    }
  }
}

static void
minkowski_spatial_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_christoffel)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        (*spatial_christoffel)[i][j][k] = 0.0;
      }
    }
  }
}

static void
minkowski_spacetime_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_christoffel)
{
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        (*spacetime_christoffel)[i][j][k] = 0.0;
      }
    }
  }
}

static void
minkowski_spatial_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_riemann_tensor)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          (*spatial_riemann_tensor)[i][j][k][l] = 0.0;
        }
      }
    }
  }
}

static void
minkowski_spacetime_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_riemann_tensor)
{
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        for (int l = 0; l < 4; l++) {
          (*spacetime_riemann_tensor)[i][j][k][l] = 0.0;
        }
      }
    }
  }
}

static void
minkowski_spatial_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** spatial_ricci_tensor)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*spatial_ricci_tensor)[i][j] = 0.0;
    }
  }
}

static void
minkowski_spacetime_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double*** spacetime_ricci_tensor)
{
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      (*spacetime_ricci_tensor)[i][j] = 0.0;
    }
  }
}

static void
minkowski_spatial_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double* spatial_ricci_scalar)
{
  *spatial_ricci_scalar = 0.0;
}

static void
minkowski_spacetime_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double* spacetime_ricci_scalar)
{
  *spacetime_ricci_scalar = 0.0;
}

static void
minkowski_spatial_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_weyl_tensor)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          (*spatial_weyl_tensor)[i][j][k][l] = 0.0;
        }
      }
    }
  }
}

static void
minkowski_spacetime_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_weyl_tensor)
{
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        for (int l = 0; l < 4; l++) {
          (*spacetime_weyl_tensor)[i][j][k][l] = 0.0;
        }
      }
    }
  }
}

static void
minkowski_extrinsic_curvature_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** extrinsic_curvature_tensor)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*extrinsic_curvature_tensor)[i][j] = 0.0;
    }
  }
}

static void
minkowski_conformal_factor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* conformal_factor)
{
  *conformal_factor = 1.0;
}

static void
minkowski_bssn_conformal_factor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* bssn_conformal_factor)
{
  *bssn_conformal_factor = 1.0;
}

static void
minkowski_conformal_factor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** conformal_factor_der)
{
  for (int i = 0; i < 3; i++) {
    (*conformal_factor_der)[i] = 0.0;
  }
}

static void
minkowski_bssn_conformal_factor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** bssn_conformal_factor_der)
{
  for (int i = 0; i < 3; i++) {
    (*bssn_conformal_factor_der)[i] = 0.0;
  }
}

static void
minkowski_conformal_factor_der2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** conformal_factor_der2)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*conformal_factor_der2)[i][j] = 0.0;
    }
  }
}

static void
minkowski_bssn_conformal_factor_der2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** bssn_conformal_factor_der2)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*bssn_conformal_factor_der2)[i][j] = 0.0;
    }
  }
}

static void
minkowski_excision_region(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* in_excision_region)
{
  *in_excision_region = false;
}

void
gkyl_gr_minkowski_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_gr_spacetime* base = container_of(ref, struct gkyl_gr_spacetime, ref_count);

  if (gkyl_gr_spacetime_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct gr_minkowski *gr_minkowski = container_of(base->on_dev, struct gr_minkowski, spacetime);
    gkyl_cu_free(gr_minkowski);
  }

  struct gr_minkowski *gr_minkowski = container_of(base, struct gr_minkowski, spacetime);
  gkyl_free(gr_minkowski);
}

struct gkyl_gr_spacetime*
gkyl_gr_minkowski_new(bool use_gpu)
{
  return gkyl_gr_minkowski_inew(&(struct gkyl_gr_minkowski_inp) {
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_gr_spacetime*
gkyl_gr_minkowski_inew(const struct gkyl_gr_minkowski_inp* inp)
{
  struct gr_minkowski *gr_minkowski = gkyl_malloc(sizeof(struct gr_minkowski));

  gr_minkowski->spacetime.spatial_metric_tensor_func = minkowski_spatial_metric_tensor;
  gr_minkowski->spacetime.spacetime_metric_tensor_func = minkowski_spacetime_metric_tensor;

  gr_minkowski->spacetime.spatial_inv_metric_tensor_func = minkowski_spatial_inv_metric_tensor;
  gr_minkowski->spacetime.spacetime_inv_metric_tensor_func = minkowski_spacetime_inv_metric_tensor;

  gr_minkowski->spacetime.spatial_metric_det_func = minkowski_spatial_metric_det;
  gr_minkowski->spacetime.spacetime_metric_det_func = minkowski_spacetime_metric_det;

  gr_minkowski->spacetime.spatial_metric_tensor_der_func = minkowski_spatial_metric_tensor_der;
  gr_minkowski->spacetime.spacetime_metric_tensor_der_func = minkowski_spacetime_metric_tensor_der;

  gr_minkowski->spacetime.lapse_function_func = minkowski_lapse_function;
  gr_minkowski->spacetime.shift_vector_func = minkowski_shift_vector;

  gr_minkowski->spacetime.lapse_function_der_func = minkowski_lapse_function_der;
  gr_minkowski->spacetime.shift_vector_der_func = minkowski_shift_vector_der;

  gr_minkowski->spacetime.spatial_christoffel_func = minkowski_spatial_christoffel;
  gr_minkowski->spacetime.spacetime_christoffel_func = minkowski_spacetime_christoffel;

  gr_minkowski->spacetime.spatial_riemann_tensor_func = minkowski_spatial_riemann_tensor;
  gr_minkowski->spacetime.spacetime_riemann_tensor_func = minkowski_spacetime_riemann_tensor;

  gr_minkowski->spacetime.spatial_ricci_tensor_func = minkowski_spatial_ricci_tensor;
  gr_minkowski->spacetime.spacetime_ricci_tensor_func = minkowski_spacetime_ricci_tensor;

  gr_minkowski->spacetime.spatial_ricci_scalar_func = minkowski_spatial_ricci_scalar;
  gr_minkowski->spacetime.spacetime_ricci_scalar_func = minkowski_spacetime_ricci_scalar;

  gr_minkowski->spacetime.spatial_weyl_tensor_func = minkowski_spatial_weyl_tensor;
  gr_minkowski->spacetime.spacetime_weyl_tensor_func = minkowski_spacetime_weyl_tensor;

  gr_minkowski->spacetime.extrinsic_curvature_tensor_func = minkowski_extrinsic_curvature_tensor;

  gr_minkowski->spacetime.conformal_factor_func = minkowski_conformal_factor;
  gr_minkowski->spacetime.bssn_conformal_factor_func = minkowski_bssn_conformal_factor;

  gr_minkowski->spacetime.conformal_factor_der_func = minkowski_conformal_factor_der;
  gr_minkowski->spacetime.bssn_conformal_factor_der_func = minkowski_bssn_conformal_factor_der;

  gr_minkowski->spacetime.conformal_factor_der2_func = minkowski_conformal_factor_der2;
  gr_minkowski->spacetime.bssn_conformal_factor_der2_func = minkowski_bssn_conformal_factor_der2;

  gr_minkowski->spacetime.excision_region_func = minkowski_excision_region;

  gr_minkowski->spacetime.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_minkowski->spacetime.flags);
  gr_minkowski->spacetime.ref_count = gkyl_ref_count_init(gkyl_gr_minkowski_free);
  gr_minkowski->spacetime.on_dev = &gr_minkowski->spacetime; // On the CPU, the spacetime object points to itself.

  return &gr_minkowski->spacetime;
}