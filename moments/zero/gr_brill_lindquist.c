#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gr_brill_lindquist.h>
#include <gkyl_gr_spacetime_diff.h>

double
brill_lindquist_phi(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_brill_lindquist *brill_lindquist = container_of(spacetime, struct gr_brill_lindquist, spacetime);

  double mass1 = brill_lindquist->mass1;
  double mass2 = brill_lindquist->mass2;

  double pos_x1 = brill_lindquist->pos_x1;
  double pos_y1 = brill_lindquist->pos_y1;
  double pos_z1 = brill_lindquist->pos_z1;

  double pos_x2 = brill_lindquist->pos_x2;
  double pos_y2 = brill_lindquist->pos_y2;
  double pos_z2 = brill_lindquist->pos_z2;

  double radial_12 = sqrt(((pos_x1 - pos_x2) * (pos_x1 - pos_x2)) + ((pos_y1 - pos_y2) * (pos_y1 - pos_y2)) + ((pos_z1 - pos_z2) * (pos_z1 - pos_z2)));
  double radial_1 = sqrt(((x - pos_x1) * (x - pos_x1)) + ((y - pos_y1) * (y - pos_y1)) + ((z - pos_z1) * (z - pos_z1)));
  double radial_2 = sqrt(((x - pos_x2) * (x - pos_x2)) + ((y - pos_y2) * (y - pos_y2)) + ((z - pos_z2) * (z - pos_z2)));

  double alpha1 = -(0.25 * ((2.0 * radial_12) + mass2 - mass1)) +
    ((0.25 * radial_12) * sqrt(4.0 + ((4.0 / radial_12) * (mass1 + mass2)) + (((mass1 - mass2) / radial_12) * (mass1 - mass2) / radial_12)));
  double alpha2 = -(0.25 * ((2.0 * radial_12) + mass1 - mass2)) +
    ((0.25 * radial_12) * sqrt(4.0 + ((4.0 / radial_12) * (mass1 + mass2)) + (((mass1 - mass2) / radial_12) * (mass1 - mass2) / radial_12)));
  
  return 8.0 * ((alpha1 / radial_1) + (alpha2 / radial_2));
}

double
brill_lindquist_psi(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_brill_lindquist *brill_lindquist = container_of(spacetime, struct gr_brill_lindquist, spacetime);

  double mass1 = brill_lindquist->mass1;
  double mass2 = brill_lindquist->mass2;

  double pos_x1 = brill_lindquist->pos_x1;
  double pos_y1 = brill_lindquist->pos_y1;
  double pos_z1 = brill_lindquist->pos_z1;

  double pos_x2 = brill_lindquist->pos_x2;
  double pos_y2 = brill_lindquist->pos_y2;
  double pos_z2 = brill_lindquist->pos_z2;

  double radial_12 = sqrt(((pos_x1 - pos_x2) * (pos_x1 - pos_x2)) + ((pos_y1 - pos_y2) * (pos_y1 - pos_y2)) + ((pos_z1 - pos_z2) * (pos_z1 - pos_z2)));
  double radial_1 = sqrt(((x - pos_x1) * (x - pos_x1)) + ((y - pos_y1) * (y - pos_y1)) + ((z - pos_z1) * (z - pos_z1)));
  double radial_2 = sqrt(((x - pos_x2) * (x - pos_x2)) + ((y - pos_y2) * (y - pos_y2)) + ((z - pos_z2) * (z - pos_z2)));

  double alpha1 = -(0.25 * ((2.0 * radial_12) + mass2 - mass1)) +
    ((0.25 * radial_12) * sqrt(4.0 + ((4.0 / radial_12) * (mass1 + mass2)) + (((mass1 - mass2) / radial_12) * (mass1 - mass2) / radial_12)));
  double alpha2 = -(0.25 * ((2.0 * radial_12) + mass1 - mass2)) +
    ((0.25 * radial_12) * sqrt(4.0 + ((4.0 / radial_12) * (mass1 + mass2)) + (((mass1 - mass2) / radial_12) * (mass1 - mass2) / radial_12)));

  double beta1 = alpha1 * ((radial_12 + alpha1 - alpha2) / (radial_12 + alpha1 + alpha2));
  double beta2 = alpha2 * ((radial_12 + alpha2 - alpha1) / (radial_12 + alpha1 + alpha2));
  
  return 8.0 * ((beta1 / radial_1) + (beta2 / radial_2));
}

static void
brill_lindquist_spatial_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_metric_tensor)
{
  double phi = brill_lindquist_phi(spacetime, x, y, z);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        (*spatial_metric_tensor)[i][j] = (1.0 + (0.125 * phi)) * (1.0 + (0.125 * phi)) * (1.0 + (0.125 * phi)) * (1.0 + (0.125 * phi));
      }
      else {
        (*spatial_metric_tensor)[i][j] = 0.0;
      }
    }
  }
}

static void
brill_lindquist_spacetime_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_metric_tensor)
{
  double phi = brill_lindquist_phi(spacetime, x, y, z);
  double psi = brill_lindquist_psi(spacetime, x, y, z);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      if (i == j) {
        if (i == 0) {
          (*spacetime_metric_tensor)[i][j] = -((1.0 - (0.125 * psi)) / (1.0 + (0.5 * phi))) * ((1.0 - (0.125 * psi)) / (1.0 + (0.5 * phi)));
        }
        else {
          (*spacetime_metric_tensor)[i][j] = (1.0 + (0.125 * phi)) * (1.0 + (0.125 * phi)) * (1.0 + (0.125 * phi)) * (1.0 + (0.125 * phi));
        }
      }
      else {
        (*spacetime_metric_tensor)[i][j] = 0.0;
      }
    }
  }
}

static void
brill_lindquist_spatial_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_inv_metric_tensor)
{
  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  brill_lindquist_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);
  double spatial_metric_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += spatial_metric[i][i];
  }

  double **spatial_metric_sq = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_sq[i] = gkyl_malloc(sizeof(double[3]));
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      spatial_metric_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_sq[i][j] += spatial_metric[i][k] * spatial_metric[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += spatial_metric_sq[i][i];
  }

  double **euclidean_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    euclidean_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        euclidean_metric[i][j] = 1.0;
      }
      else {
        euclidean_metric[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*spatial_inv_metric_tensor)[i][j] = (1.0 / spatial_metric_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * spatial_metric[i][j]) + spatial_metric_sq[i][j]);
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(spatial_metric_sq[i]);
    gkyl_free(euclidean_metric[i]);
  }
  gkyl_free(spatial_metric);
  gkyl_free(spatial_metric_sq);
  gkyl_free(euclidean_metric);
}

static void
brill_lindquist_spacetime_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_inv_metric_tensor)
{
  double** spatial_metric = gkyl_malloc(sizeof(double*[3]));
  double** inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  brill_lindquist_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);
  brill_lindquist_spatial_inv_metric_tensor(spacetime, t, x, y, z, &inv_spatial_metric);

  double lapse_function;
  double* shift_vector = gkyl_malloc(sizeof(double[3]));
  brill_lindquist_lapse_function(spacetime, t, x, y, z, &lapse_function);
  brill_lindquist_shift_vector(spacetime, t, x, y, z, &shift_vector);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      (*spacetime_inv_metric_tensor)[i][j] = 0.0;
    }
  }

  (*spacetime_inv_metric_tensor)[0][0] = -1.0 / (lapse_function * lapse_function);
  for (int i = 0; i < 3; i++) {
    (*spacetime_inv_metric_tensor)[0][i + 1] = shift_vector[i] / (lapse_function * lapse_function);
    (*spacetime_inv_metric_tensor)[i + 1][0] = shift_vector[i] / (lapse_function * lapse_function);

    for (int j = 0; j < 3; j++) {
      (*spacetime_inv_metric_tensor)[i + 1][j + 1] = inv_spatial_metric[i][j] - (shift_vector[i] * shift_vector[j]) / (lapse_function * lapse_function);
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(spatial_metric);
  gkyl_free(inv_spatial_metric);
  gkyl_free(shift_vector);
}

static void
brill_lindquist_spatial_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spatial_metric_det)
{
  double** spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  brill_lindquist_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);

  *spatial_metric_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
  }
  gkyl_free(spatial_metric);
}

static void
brill_lindquist_spacetime_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spacetime_metric_det)
{
  double spatial_metric_det;
  double lapse_function;
  brill_lindquist_spatial_metric_det(spacetime, t, x, y, z, &spatial_metric_det);
  brill_lindquist_lapse_function(spacetime, t, x, y, z, &lapse_function);

  *spacetime_metric_det = - (lapse_function * lapse_function) * spatial_metric_det;
}

static void
brill_lindquist_spatial_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
   const double dx, const double dy, const double dz, double**** spatial_metric_tensor_der)
{
  gkyl_gr_spatial_metric_tensor_diff(spacetime, t, x, y, z, dx, dy, dz, spatial_metric_tensor_der);
}

static void
brill_lindquist_spacetime_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_metric_tensor_der)
{
  gkyl_gr_spacetime_metric_tensor_diff(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_metric_tensor_der);
}

static void
brill_lindquist_lapse_function(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* lapse_function)
{
  double phi = brill_lindquist_phi(spacetime, x, y, z);
  double psi = brill_lindquist_psi(spacetime, x, y, z);

  *lapse_function = (1.0 - (0.125 * psi)) / (1.0 + (0.5 * phi));
}

static void
brill_lindquist_shift_vector(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double** shift_vector)
{
  for (int i = 0; i < 3; i++) {
    (*shift_vector)[i] = 0.0;
  }
}

static void
brill_lindquist_lapse_function_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** lapse_function_der)
{
  gkyl_gr_lapse_function_diff(spacetime, t, x, y, z, dx, dy, dz, lapse_function_der);
}

static void
brill_lindquist_shift_vector_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** shift_vector_der)
{
  gkyl_gr_shift_vector_diff(spacetime, t, x, y, z, dx, dy, dz, shift_vector_der);
}

static void
brill_lindquist_spatial_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_christoffel)
{
  gkyl_gr_spatial_christoffel_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_christoffel);
}

static void
brill_lindquist_spacetime_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_christoffel)
{
  gkyl_gr_spacetime_christoffel_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_christoffel);
}

static void
brill_lindquist_spatial_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_riemann_tensor)
{
  gkyl_gr_spatial_riemann_tensor_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_riemann_tensor);
}

static void
brill_lindquist_spacetime_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_riemann_tensor)
{
  gkyl_gr_spacetime_riemann_tensor_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_riemann_tensor);
}

static void
brill_lindquist_spatial_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** spatial_ricci_tensor)
{
  gkyl_gr_spatial_ricci_tensor_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_ricci_tensor);
}

static void
brill_lindquist_spacetime_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double*** spacetime_ricci_tensor)
{
  gkyl_gr_spacetime_ricci_tensor_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_ricci_tensor);
}

static void
brill_lindquist_spatial_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double* spatial_ricci_scalar)
{
  gkyl_gr_spatial_ricci_scalar_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_ricci_scalar);
}

static void
brill_lindquist_spacetime_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double* spacetime_ricci_scalar)
{
  gkyl_gr_spacetime_ricci_scalar_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_ricci_scalar);
}

static void
brill_lindquist_spatial_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_weyl_tensor)
{
  gkyl_gr_spatial_weyl_tensor_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_weyl_tensor);
}

static void
brill_lindquist_spacetime_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_weyl_tensor)
{
  gkyl_gr_spacetime_weyl_tensor_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_weyl_tensor);
}

static void
brill_lindquist_extrinsic_curvature_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** extrinsic_curvature_tensor)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*extrinsic_curvature_tensor)[i][j] = 0.0;
    }
  }
}

static void
brill_lindquist_conformal_factor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* conformal_factor)
{
  double spatial_metric_det;
  brill_lindquist_spatial_metric_det(spacetime, t, x, y, z, &spatial_metric_det);

  *conformal_factor = pow(spatial_metric_det, 1.0 / 12.0);
}

static void
brill_lindquist_bssn_conformal_factor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* bssn_conformal_factor)
{
  double spatial_metric_det;
  brill_lindquist_spatial_metric_det(spacetime, t, x, y, z, &spatial_metric_det);

  *bssn_conformal_factor = 1.0 / pow(spatial_metric_det, 1.0 / 6.0);
}

static void
brill_lindquist_conformal_factor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** conformal_factor_der)
{
  gkyl_gr_conformal_factor_diff(spacetime, t, x, y, z, dx, dy, dz, conformal_factor_der);
}

static void
brill_lindquist_bssn_conformal_factor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** bssn_conformal_factor_der)
{
  gkyl_gr_bssn_conformal_factor_diff(spacetime, t, x, y, z, dx, dy, dz, bssn_conformal_factor_der);
}

static void
brill_lindquist_conformal_factor_der2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** conformal_factor_der2)
{
  gkyl_gr_conformal_factor_diff2(spacetime, t, x, y, z, dx, dy, dz, conformal_factor_der2);
}

static void
brill_lindquist_bssn_conformal_factor_der2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** bssn_conformal_factor_der2)
{
  gkyl_gr_bssn_conformal_factor_diff2(spacetime, t, x, y, z, dx, dy, dz, bssn_conformal_factor_der2);
}

static void
brill_lindquist_excision_region(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* in_excision_region)
{
  const struct gr_brill_lindquist *brill_lindquist = container_of(spacetime, struct gr_brill_lindquist, spacetime);

  double mass1 = brill_lindquist->mass1;
  double mass2 = brill_lindquist->mass2;

  double pos_x1 = brill_lindquist->pos_x1;
  double pos_y1 = brill_lindquist->pos_y1;
  double pos_z1 = brill_lindquist->pos_z1;

  double pos_x2 = brill_lindquist->pos_x2;
  double pos_y2 = brill_lindquist->pos_y2;
  double pos_z2 = brill_lindquist->pos_z2;

  double r1 = sqrt(((x - pos_x1) * (x - pos_x1)) + ((y - pos_y1) * (y - pos_y1)) + ((z - pos_z1) * (z - pos_z1)));
  double r2 = sqrt(((x - pos_x2) * (x - pos_x2)) + ((y - pos_y2) * (y - pos_y2)) + ((z - pos_z2) * (z - pos_z2)));

  if (r1 <= 2.0 * mass1 || r2 <= 2.0 * mass2) {
    *in_excision_region = true;
  }
  else {
    *in_excision_region = false;
  }
}

void
gkyl_gr_brill_lindquist_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_gr_spacetime* base = container_of(ref, struct gkyl_gr_spacetime, ref_count);

  if (gkyl_gr_spacetime_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct gr_brill_lindquist *gr_brill_lindquist = container_of(base->on_dev, struct gr_brill_lindquist, spacetime);
    gkyl_cu_free(gr_brill_lindquist);
  }

  struct gr_brill_lindquist *gr_brill_lindquist = container_of(base, struct gr_brill_lindquist, spacetime);
  gkyl_free(gr_brill_lindquist);
}

struct gkyl_gr_spacetime*
gkyl_gr_brill_lindquist_new(bool use_gpu, double mass1, double mass2, double pos_x1, double pos_y1, double pos_z1, double pos_x2, double pos_y2, double pos_z2)
{
  return gkyl_gr_brill_lindquist_inew(&(struct gkyl_gr_brill_lindquist_inp) {
      .use_gpu = use_gpu,
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
    }
  );
}

struct gkyl_gr_spacetime*
gkyl_gr_brill_lindquist_inew(const struct gkyl_gr_brill_lindquist_inp* inp)
{
  struct gr_brill_lindquist *gr_brill_lindquist = gkyl_malloc(sizeof(struct gr_brill_lindquist));

  gr_brill_lindquist->mass1 = inp->mass1;
  gr_brill_lindquist->mass2 = inp->mass2;

  gr_brill_lindquist->pos_x1 = inp->pos_x1;
  gr_brill_lindquist->pos_y1 = inp->pos_y1;
  gr_brill_lindquist->pos_z1 = inp->pos_z1;

  gr_brill_lindquist->pos_x2 = inp->pos_x2;
  gr_brill_lindquist->pos_y2 = inp->pos_y2;
  gr_brill_lindquist->pos_z2 = inp->pos_z2;

  gr_brill_lindquist->spacetime.spatial_metric_tensor_func = brill_lindquist_spatial_metric_tensor;
  gr_brill_lindquist->spacetime.spacetime_metric_tensor_func = brill_lindquist_spacetime_metric_tensor;

  gr_brill_lindquist->spacetime.spatial_inv_metric_tensor_func = brill_lindquist_spatial_inv_metric_tensor;
  gr_brill_lindquist->spacetime.spacetime_inv_metric_tensor_func = brill_lindquist_spacetime_inv_metric_tensor;

  gr_brill_lindquist->spacetime.spatial_metric_det_func = brill_lindquist_spatial_metric_det;
  gr_brill_lindquist->spacetime.spacetime_metric_det_func = brill_lindquist_spacetime_metric_det;

  gr_brill_lindquist->spacetime.spatial_metric_tensor_der_func = brill_lindquist_spatial_metric_tensor_der;
  gr_brill_lindquist->spacetime.spacetime_metric_tensor_der_func = brill_lindquist_spacetime_metric_tensor_der;

  gr_brill_lindquist->spacetime.lapse_function_func = brill_lindquist_lapse_function;
  gr_brill_lindquist->spacetime.shift_vector_func = brill_lindquist_shift_vector;

  gr_brill_lindquist->spacetime.lapse_function_der_func = brill_lindquist_lapse_function_der;
  gr_brill_lindquist->spacetime.shift_vector_der_func = brill_lindquist_shift_vector_der;

  gr_brill_lindquist->spacetime.spatial_christoffel_func = brill_lindquist_spatial_christoffel;
  gr_brill_lindquist->spacetime.spacetime_christoffel_func = brill_lindquist_spacetime_christoffel;

  gr_brill_lindquist->spacetime.spatial_riemann_tensor_func = brill_lindquist_spatial_riemann_tensor;
  gr_brill_lindquist->spacetime.spacetime_riemann_tensor_func = brill_lindquist_spacetime_riemann_tensor;

  gr_brill_lindquist->spacetime.spatial_ricci_tensor_func = brill_lindquist_spatial_ricci_tensor;
  gr_brill_lindquist->spacetime.spacetime_ricci_tensor_func = brill_lindquist_spacetime_ricci_tensor;

  gr_brill_lindquist->spacetime.spatial_ricci_scalar_func = brill_lindquist_spatial_ricci_scalar;
  gr_brill_lindquist->spacetime.spacetime_ricci_scalar_func = brill_lindquist_spacetime_ricci_scalar;

  gr_brill_lindquist->spacetime.spatial_weyl_tensor_func = brill_lindquist_spatial_weyl_tensor;
  gr_brill_lindquist->spacetime.spacetime_weyl_tensor_func = brill_lindquist_spacetime_weyl_tensor;

  gr_brill_lindquist->spacetime.extrinsic_curvature_tensor_func = brill_lindquist_extrinsic_curvature_tensor;

  gr_brill_lindquist->spacetime.conformal_factor_func = brill_lindquist_conformal_factor;
  gr_brill_lindquist->spacetime.bssn_conformal_factor_func = brill_lindquist_bssn_conformal_factor;

  gr_brill_lindquist->spacetime.conformal_factor_der_func = brill_lindquist_conformal_factor_der;
  gr_brill_lindquist->spacetime.bssn_conformal_factor_der_func = brill_lindquist_bssn_conformal_factor_der;

  gr_brill_lindquist->spacetime.conformal_factor_der2_func = brill_lindquist_conformal_factor_der2;
  gr_brill_lindquist->spacetime.bssn_conformal_factor_der2_func = brill_lindquist_bssn_conformal_factor_der2;

  gr_brill_lindquist->spacetime.excision_region_func = brill_lindquist_excision_region;

  gr_brill_lindquist->spacetime.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_brill_lindquist->spacetime.flags);
  gr_brill_lindquist->spacetime.ref_count = gkyl_ref_count_init(gkyl_gr_brill_lindquist_free);
  gr_brill_lindquist->spacetime.on_dev = &gr_brill_lindquist->spacetime; // On the CPU, the spacetime object points to itself.

  return &gr_brill_lindquist->spacetime;
}