#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gr_blackhole_isotropic.h>
#include <gkyl_gr_spacetime_diff.h>

static void
blackhole_isotropic_spatial_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_metric_tensor)
{
  const struct gr_blackhole_isotropic *blackhole_isotropic = container_of(spacetime, struct gr_blackhole_isotropic, spacetime);

  double mass = blackhole_isotropic->mass;
  double spin = blackhole_isotropic->spin;

  double pos_x = blackhole_isotropic->pos_x;
  double pos_y = blackhole_isotropic->pos_y;
  double pos_z = blackhole_isotropic->pos_z;

  double Mspin = mass * spin;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
  double sin_theta = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double cos_theta = (z - pos_z) / rho;

  double r_BL = rho * (1.0 + ((mass + Mspin) / (2.0 * rho))) * (1.0 + ((mass - Mspin) / (2.0 * rho)));
  double sigma = (r_BL * r_BL) + ((Mspin * Mspin) * (cos_theta * cos_theta));
  double delta = (r_BL * r_BL) - (2.0 * mass * r_BL) + (Mspin * Mspin);

  double A_BL = (((r_BL * r_BL) + (Mspin * Mspin)) * ((r_BL * r_BL) + (Mspin + Mspin))) - ((Mspin * Mspin) * delta * (sin_theta * sin_theta));
  double A_func = sqrt(sigma) / rho;
  double B_func = sqrt(A_BL) / (rho * sqrt(sigma));

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*spatial_metric_tensor)[i][j] = 0.0;
    }
  }

  (*spatial_metric_tensor)[0][0] = ((A_func * A_func) + (((B_func * B_func) - (A_func * A_func)) * ((y - pos_y) * (y - pos_y)))) /
    (((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  (*spatial_metric_tensor)[1][1] = ((A_func * A_func) + (((B_func * B_func) - (A_func * A_func)) * ((y - pos_y) * (y - pos_y)))) /
    (((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  (*spatial_metric_tensor)[2][2] = (A_func * A_func);

  (*spatial_metric_tensor)[0][1] = -(((B_func * B_func) - (A_func * A_func)) * ((x - pos_x) * (y - pos_y))) /
    (((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  (*spatial_metric_tensor)[1][0] = -(((B_func * B_func) - (A_func * A_func)) * ((x - pos_x) * (y - pos_y))) /
    (((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
}

static void
blackhole_isotropic_spacetime_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_metric_tensor)
{
  double** spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  blackhole_isotropic_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);

  double lapse_function;
  double* shift_vector = gkyl_malloc(sizeof(double[3]));
  blackhole_isotropic_lapse_function(spacetime, t, x, y, z, &lapse_function);
  blackhole_isotropic_shift_vector(spacetime, t, x, y, z, &shift_vector);

  double cov_shift_vector[3];
  for (int i = 0; i < 3; i++) {
    cov_shift_vector[0] = 0.0;


    for (int j = 0; j < 3; j++) {
      cov_shift_vector[i] += spatial_metric[i][j] * shift_vector[j];
    }
  }

  (*spacetime_metric_tensor)[0][0] = -(lapse_function * lapse_function);
  for (int i = 0; i < 3; i++) {
    (*spacetime_metric_tensor)[0][0] += cov_shift_vector[i] * shift_vector[i];

    (*spacetime_metric_tensor)[0][i + 1] = cov_shift_vector[i];
    (*spacetime_metric_tensor)[i + 1][0] = cov_shift_vector[i];

    for (int j = 0; j < 3; j++) {
      (*spacetime_metric_tensor)[i + 1][j + 1] = spatial_metric[i][j];
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
  }
  gkyl_free(spatial_metric);
  gkyl_free(shift_vector);
}

static void
blackhole_isotropic_spatial_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_inv_metric_tensor)
{
  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  blackhole_isotropic_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);
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
blackhole_isotropic_spacetime_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_inv_metric_tensor)
{
  double** inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  blackhole_isotropic_spatial_inv_metric_tensor(spacetime, t, x, y, z, &inv_spatial_metric);

  double lapse_function;
  double* shift_vector = gkyl_malloc(sizeof(double[3]));
  blackhole_isotropic_lapse_function(spacetime, t, x, y, z, &lapse_function);
  blackhole_isotropic_shift_vector(spacetime, t, x, y, z, &shift_vector);

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
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(inv_spatial_metric);
  gkyl_free(shift_vector);
}

static void
blackhole_isotropic_spatial_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spatial_metric_det)
{
  double** spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  blackhole_isotropic_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);

  *spatial_metric_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
  }
  gkyl_free(spatial_metric);
}

static void
blackhole_isotropic_spacetime_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spacetime_metric_det)
{
  double spatial_metric_det;
  double lapse_function;
  blackhole_isotropic_spatial_metric_det(spacetime, t, x, y, z, &spatial_metric_det);
  blackhole_isotropic_lapse_function(spacetime, t, x, y, z, &lapse_function);

  *spacetime_metric_det = - (lapse_function * lapse_function) * spatial_metric_det;
}

static void
blackhole_isotropic_spatial_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
   const double dx, const double dy, const double dz, double**** spatial_metric_tensor_der)
{
  gkyl_gr_spatial_metric_tensor_diff(spacetime, t, x, y, z, dx, dy, dz, spatial_metric_tensor_der);
}

static void
blackhole_isotropic_spacetime_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_metric_tensor_der)
{
  gkyl_gr_spacetime_metric_tensor_diff(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_metric_tensor_der);
}

static void
blackhole_isotropic_lapse_function(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* lapse_function)
{
  const struct gr_blackhole_isotropic *blackhole_isotropic = container_of(spacetime, struct gr_blackhole_isotropic, spacetime);

  double mass = blackhole_isotropic->mass;
  double spin = blackhole_isotropic->spin;

  double pos_x = blackhole_isotropic->pos_x;
  double pos_y = blackhole_isotropic->pos_y;
  double pos_z = blackhole_isotropic->pos_z;

  double Mspin = mass * spin;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
  double sin_theta = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double cos_theta = (z - pos_z) / rho;

  double r_BL = rho * (1.0 + ((mass + Mspin) / (2.0 * rho))) * (1.0 + ((mass - Mspin) / (2.0 * rho)));
  double sigma = (r_BL * r_BL) + ((Mspin * Mspin) * (cos_theta * cos_theta));
  double delta = (r_BL * r_BL) - (2.0 * mass * r_BL) + (Mspin * Mspin);

  double A_BL = (((r_BL * r_BL) + (Mspin * Mspin)) * ((r_BL * r_BL) + (Mspin + Mspin))) - ((Mspin * Mspin) * delta * (sin_theta * sin_theta));

  *lapse_function = sqrt((sigma * delta) / A_BL);
}

static void
blackhole_isotropic_shift_vector(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double** shift_vector)
{
  const struct gr_blackhole_isotropic *blackhole_isotropic = container_of(spacetime, struct gr_blackhole_isotropic, spacetime);

  double mass = blackhole_isotropic->mass;
  double spin = blackhole_isotropic->spin;

  double pos_x = blackhole_isotropic->pos_x;
  double pos_y = blackhole_isotropic->pos_y;
  double pos_z = blackhole_isotropic->pos_z;

  double Mspin = mass * spin;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
  double sin_theta = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double cos_theta = (z - pos_z) / rho;

  double r_BL = rho * (1.0 + ((mass + Mspin) / (2.0 * rho))) * (1.0 + ((mass - Mspin) / (2.0 * rho)));
  double delta = (r_BL * r_BL) - (2.0 * mass * r_BL) + (Mspin * Mspin);

  double A_BL = (((r_BL * r_BL) + (Mspin * Mspin)) * ((r_BL * r_BL) + (Mspin + Mspin))) - ((Mspin * Mspin) * delta * (sin_theta * sin_theta));
  double omega = (2.0 * mass * Mspin * r_BL) / (A_BL);

  (*shift_vector)[0] = omega * (y - pos_y);
  (*shift_vector)[1] = -omega * (x - pos_x);
  (*shift_vector)[2] = 0.0;
}

static void
blackhole_isotropic_lapse_function_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** lapse_function_der)
{
  gkyl_gr_lapse_function_diff(spacetime, t, x, y, z, dx, dy, dz, lapse_function_der);
}

static void
blackhole_isotropic_shift_vector_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** shift_vector_der)
{
  gkyl_gr_shift_vector_diff(spacetime, t, x, y, z, dx, dy, dz, shift_vector_der);
}

static void
blackhole_isotropic_spatial_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_christoffel)
{
  gkyl_gr_spatial_christoffel_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_christoffel);
}

static void
blackhole_isotropic_spacetime_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_christoffel)
{
  gkyl_gr_spacetime_christoffel_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_christoffel);
}

static void
blackhole_isotropic_spatial_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_riemann_tensor)
{
  gkyl_gr_spatial_riemann_tensor_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_riemann_tensor);
}

static void
blackhole_isotropic_spacetime_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_riemann_tensor)
{
  gkyl_gr_spacetime_riemann_tensor_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_riemann_tensor);
}

static void
blackhole_isotropic_spatial_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** spatial_ricci_tensor)
{
  gkyl_gr_spatial_ricci_tensor_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_ricci_tensor);
}

static void
blackhole_isotropic_spacetime_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double*** spacetime_ricci_tensor)
{
  gkyl_gr_spacetime_ricci_tensor_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_ricci_tensor);
}

static void
blackhole_isotropic_spatial_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double* spatial_ricci_scalar)
{
  gkyl_gr_spatial_ricci_scalar_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_ricci_scalar);
}

static void
blackhole_isotropic_spacetime_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double* spacetime_ricci_scalar)
{
  gkyl_gr_spacetime_ricci_scalar_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_ricci_scalar);
}

static void
blackhole_isotropic_spatial_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_weyl_tensor)
{
  gkyl_gr_spatial_weyl_tensor_fd(spacetime, t, x, y, z, dx, dy, dx, spatial_weyl_tensor);
}

static void
blackhole_isotropic_spacetime_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_weyl_tensor)
{
  gkyl_gr_spacetime_weyl_tensor_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_weyl_tensor);
}

static void
blackhole_isotropic_extrinsic_curvature_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** extrinsic_curvature_tensor)
{
  double lapse_function;
  double *shift_vector = gkyl_malloc(sizeof(double[3]));

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  double **shift_vector_der = gkyl_malloc(sizeof(double*[3]));
  double **shift_vector_cov_der = gkyl_malloc(sizeof(double*[3]));
  double **shift_covector_cov_der = gkyl_malloc(sizeof(double*[3]));

  double ***spatial_christoffel = gkyl_malloc(sizeof(double**[3]));

  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
    shift_vector_der[i] = gkyl_malloc(sizeof(double[3]));
    shift_vector_cov_der[i] = gkyl_malloc(sizeof(double[3]));
    shift_covector_cov_der[i] = gkyl_malloc(sizeof(double[3]));

    spatial_christoffel[i] = gkyl_malloc(sizeof(double*[3]));
    for (int j = 0; j < 3; j++) {
      spatial_christoffel[i][j] = gkyl_malloc(sizeof(double[3]));

      shift_covector_cov_der[i][j] = 0.0;
    }
  }

  blackhole_isotropic_lapse_function(spacetime, t, x, y, z, &lapse_function);
  blackhole_isotropic_shift_vector(spacetime, t, x, y, z, &shift_vector);
  blackhole_isotropic_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);
  blackhole_isotropic_shift_vector_der(spacetime, t, x, y, z, dx, dy, dz, &shift_vector_der);
  blackhole_isotropic_spatial_christoffel(spacetime, t, x, y, z, dx, dy, dz, &spatial_christoffel);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      shift_vector_cov_der[i][j] = shift_vector_der[i][j];

      for (int k = 0; k < 3; k++) {
        shift_vector_cov_der[i][j] += spatial_christoffel[j][i][k] * shift_vector[k];
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        shift_covector_cov_der[i][j] += spatial_metric[j][k] * shift_vector_cov_der[i][k];
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*extrinsic_curvature_tensor)[i][j] = -(1.0 / (2.0 * lapse_function)) * (shift_covector_cov_der[j][i] + shift_covector_cov_der[i][j]);
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(shift_vector_der[i]);
    gkyl_free(shift_vector_cov_der[i]);
    gkyl_free(shift_covector_cov_der[i]);

    for (int j = 0; j < 3; j++) {
      gkyl_free(spatial_christoffel[i][j]);
    }
    gkyl_free(spatial_christoffel[i]);
  }
  gkyl_free(shift_vector);
  gkyl_free(spatial_metric);
  gkyl_free(shift_vector_der);
  gkyl_free(shift_vector_cov_der);
  gkyl_free(shift_covector_cov_der);
  gkyl_free(spatial_christoffel);
}

static void
blackhole_isotropic_excision_region(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* in_excision_region)
{
  const struct gr_blackhole_isotropic *blackhole_isotropic = container_of(spacetime, struct gr_blackhole_isotropic, spacetime);

  double mass = blackhole_isotropic->mass;

  double pos_x = blackhole_isotropic->pos_x;
  double pos_y = blackhole_isotropic->pos_y;
  double pos_z = blackhole_isotropic->pos_z;

  double eta = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));

  if (eta <= mass * 2.0) {
    *in_excision_region = true;
  }
  else {
    *in_excision_region = false;
  }
}

void
gkyl_gr_blackhole_isotropic_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_gr_spacetime* base = container_of(ref, struct gkyl_gr_spacetime, ref_count);

  if (gkyl_gr_spacetime_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct gr_blackhole_isotropic *gr_blackhole_isotropic = container_of(base->on_dev, struct gr_blackhole_isotropic, spacetime);
    gkyl_cu_free(gr_blackhole_isotropic);
  }

  struct gr_blackhole_isotropic *gr_blackhole_isotropic = container_of(base, struct gr_blackhole_isotropic, spacetime);
  gkyl_free(gr_blackhole_isotropic);
}

struct gkyl_gr_spacetime*
gkyl_gr_blackhole_isotropic_new(bool use_gpu, double mass, double spin, double pos_x, double pos_y, double pos_z)
{
  return gkyl_gr_blackhole_isotropic_inew(&(struct gkyl_gr_blackhole_isotropic_inp) {
      .use_gpu = use_gpu,
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
    }
  );
}

struct gkyl_gr_spacetime*
gkyl_gr_blackhole_isotropic_inew(const struct gkyl_gr_blackhole_isotropic_inp* inp)
{
  struct gr_blackhole_isotropic *gr_blackhole_isotropic = gkyl_malloc(sizeof(struct gr_blackhole_isotropic));

  gr_blackhole_isotropic->mass = inp->mass;
  gr_blackhole_isotropic->spin = inp->spin;

  gr_blackhole_isotropic->pos_x = inp->pos_x;
  gr_blackhole_isotropic->pos_y = inp->pos_y;
  gr_blackhole_isotropic->pos_z = inp->pos_z;

  gr_blackhole_isotropic->spacetime.spatial_metric_tensor_func = blackhole_isotropic_spatial_metric_tensor;
  gr_blackhole_isotropic->spacetime.spacetime_metric_tensor_func = blackhole_isotropic_spacetime_metric_tensor;

  gr_blackhole_isotropic->spacetime.spatial_inv_metric_tensor_func = blackhole_isotropic_spatial_inv_metric_tensor;
  gr_blackhole_isotropic->spacetime.spacetime_inv_metric_tensor_func = blackhole_isotropic_spacetime_inv_metric_tensor;

  gr_blackhole_isotropic->spacetime.spatial_metric_det_func = blackhole_isotropic_spatial_metric_det;
  gr_blackhole_isotropic->spacetime.spacetime_metric_det_func = blackhole_isotropic_spacetime_metric_det;

  gr_blackhole_isotropic->spacetime.spatial_metric_tensor_der_func = blackhole_isotropic_spatial_metric_tensor_der;
  gr_blackhole_isotropic->spacetime.spacetime_metric_tensor_der_func = blackhole_isotropic_spacetime_metric_tensor_der;

  gr_blackhole_isotropic->spacetime.lapse_function_func = blackhole_isotropic_lapse_function;
  gr_blackhole_isotropic->spacetime.shift_vector_func = blackhole_isotropic_shift_vector;

  gr_blackhole_isotropic->spacetime.lapse_function_der_func = blackhole_isotropic_lapse_function_der;
  gr_blackhole_isotropic->spacetime.shift_vector_der_func = blackhole_isotropic_shift_vector_der;

  gr_blackhole_isotropic->spacetime.spatial_christoffel_func = blackhole_isotropic_spatial_christoffel;
  gr_blackhole_isotropic->spacetime.spacetime_christoffel_func = blackhole_isotropic_spacetime_christoffel;

  gr_blackhole_isotropic->spacetime.spatial_riemann_tensor_func = blackhole_isotropic_spatial_riemann_tensor;
  gr_blackhole_isotropic->spacetime.spacetime_riemann_tensor_func = blackhole_isotropic_spacetime_riemann_tensor;

  gr_blackhole_isotropic->spacetime.spatial_ricci_tensor_func = blackhole_isotropic_spatial_ricci_tensor;
  gr_blackhole_isotropic->spacetime.spacetime_ricci_tensor_func = blackhole_isotropic_spacetime_ricci_tensor;

  gr_blackhole_isotropic->spacetime.spatial_ricci_scalar_func = blackhole_isotropic_spatial_ricci_scalar;
  gr_blackhole_isotropic->spacetime.spacetime_ricci_scalar_func = blackhole_isotropic_spacetime_ricci_scalar;

  gr_blackhole_isotropic->spacetime.spatial_weyl_tensor_func = blackhole_isotropic_spatial_weyl_tensor;
  gr_blackhole_isotropic->spacetime.spacetime_weyl_tensor_func = blackhole_isotropic_spacetime_weyl_tensor;

  gr_blackhole_isotropic->spacetime.extrinsic_curvature_tensor_func = blackhole_isotropic_extrinsic_curvature_tensor;

  gr_blackhole_isotropic->spacetime.excision_region_func = blackhole_isotropic_excision_region;

  gr_blackhole_isotropic->spacetime.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_blackhole_isotropic->spacetime.flags);
  gr_blackhole_isotropic->spacetime.ref_count = gkyl_ref_count_init(gkyl_gr_blackhole_isotropic_free);
  gr_blackhole_isotropic->spacetime.on_dev = &gr_blackhole_isotropic->spacetime; // On the CPU, the spacetime object points to itself.

  return &gr_blackhole_isotropic->spacetime;
}