#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gr_spacetime_diff.h>
#include <gkyl_gr_tov.h>

// at Cartesian (x,y,z), return the areal radius r, the radial metric factor A(r) = 1/(1 - 2 m(r)/r), the lapse alpha = exp(Phi(r)), and the unit radial vector n^i = x^i / r (zeroed at the origin). 
// Phi(r), m(r) come from the frozen TOV table by areal-radius lookup (via the gkyl_tov_eval_bl), so interior and exterior are both covered.
static void
tov_geometry(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z, double* r_out, double* A_out, double* alpha_out, double n[3])
{
  const struct gr_tov *gr_tov = container_of(spacetime, struct gr_tov, spacetime);

  double dx = x - gr_tov->pos_x, dy = y - gr_tov->pos_y, dz = z - gr_tov->pos_z;
  double r = sqrt((dx * dx) + (dy * dy) + (dz * dz));

  struct tov_eval_bl bl = { 0 };
  gkyl_tov_eval_bl(gr_tov->tov, r, &bl);

  double A = 1.0;
  if (r > 1.0e-12) {
    double comp = 2.0 * bl.m / r;
    if (comp > 1.0 - 1.0e-12) {
      comp = 1.0 - 1.0e-12; // for a regular star 2m/r < 1 everywhere; collapsed though?
    }
    A = 1.0 / (1.0 - comp);
    n[0] = dx / r;
    n[1] = dy / r;
    n[2] = dz / r;
  }
  else {
    n[0] = 0.0;
    n[1] = 0.0;
    n[2] = 0.0;
  }

  *r_out = r;
  *A_out = A;
  *alpha_out = exp(bl.Phi);
}

static void
tov_spatial_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_metric_tensor)
{
  const struct gr_tov *gr_tov = container_of(spacetime, struct gr_tov, spacetime);

  double r, A, alpha, n[3];
  tov_geometry(spacetime, x, y, z, &r, &A, &alpha, n);

  if (gr_tov->use_kerr_schild) {
    // Kerr-Schild: gamma_ij = delta_ij + (2m/r) n_i n_j, comp = 2m/r = 1 - 1/A
    double comp = 1.0 - (1.0 / A);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (i == j) {
          (*spatial_metric_tensor)[i][j] = 1.0 + (comp * n[i] * n[j]);
        }
        else {
          (*spatial_metric_tensor)[i][j] = (comp * n[i] * n[j]);
        }
      }
    }
  }
  else {
    // Areal (static) gauge: gamma_ij = delta_ij + (A - 1) n_i n_j.
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (i == j) {
          (*spatial_metric_tensor)[i][j] = 1.0 + ((A - 1.0) * n[i] * n[j]);
        }
        else {
          (*spatial_metric_tensor)[i][j] = (A - 1.0) * n[i] * n[j];
        }
      }
    }
  }
}

static void
tov_lapse_function(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* lapse_function)
{
  const struct gr_tov *gr_tov = container_of(spacetime, struct gr_tov, spacetime);

  double r, A, alpha, n[3];
  tov_geometry(spacetime, x, y, z, &r, &A, &alpha, n);

  if (gr_tov->use_kerr_schild) {
    // Kerr-Schild alpha_KS = e^Phi / sqrt((1 - 2m/r)(1 + 2m/r)) = e^Phi / sqrt(1 - comp^2).
    // Schwarzschild-KS 1/sqrt(1+2M/r) in the exterior, where e^Phi = sqrt(1-2M/r)
    double comp = 1.0 - (1.0 / A);
    *lapse_function = alpha / sqrt((1.0 - comp) * (1.0 + comp));
  }
  else {
    *lapse_function = alpha; // areal gauge: alpha = e^Phi.
  }
}

static void
tov_shift_vector(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double** shift_vector)
{
  const struct gr_tov *gr_tov = container_of(spacetime, struct gr_tov, spacetime);

  if (gr_tov->use_kerr_schild) {
    // Ingoing Kerr-Schild shift (outward, radial): beta^i = e^Phi (2m/r) / [(1+2m/r) sqrt(1-2m/r)] n^i.
    // Schwarzschild-KS (2M/r)/(1+2M/r) n^i in the exterior
    double r, A, alpha, n[3];
    tov_geometry(spacetime, x, y, z, &r, &A, &alpha, n);
    double comp = 1.0 - (1.0 / A);
    double beta_mag = alpha * comp / ((1.0 + comp) * sqrt(1.0 - comp));
    for (int i = 0; i < 3; i++) {
      (*shift_vector)[i] = beta_mag * n[i];
    }
  }
  else {
    // Areal (static) gauge: g_ti = 0 => beta^i = 0.
    for (int i = 0; i < 3; i++) {
      (*shift_vector)[i] = 0.0;
    }
  }
}

static void
tov_spacetime_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_metric_tensor)
{
  double** spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }
  tov_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);

  double lapse_function;
  tov_lapse_function(spacetime, t, x, y, z, &lapse_function);

  double* shift_vector = gkyl_malloc(sizeof(double[3]));
  tov_shift_vector(spacetime, t, x, y, z, &shift_vector);

  double shift_cov[3] = { 0.0, 0.0, 0.0 };
  double shift_sq = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      shift_cov[i] += spatial_metric[i][j] * shift_vector[j];
    }
    shift_sq += shift_cov[i] * shift_vector[i];
  }

  // ADM line element: g_00 = -alpha^2 + beta_i beta^i, g_0i = beta_i, g_ij = gamma_ij.
  (*spacetime_metric_tensor)[0][0] = -(lapse_function * lapse_function) + shift_sq;
  for (int i = 0; i < 3; i++) {
    (*spacetime_metric_tensor)[0][i + 1] = shift_cov[i];
    (*spacetime_metric_tensor)[i + 1][0] = shift_cov[i];

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
tov_spatial_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_inv_metric_tensor)
{
  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  tov_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);
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
tov_spacetime_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_inv_metric_tensor)
{
  double** inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  tov_spatial_inv_metric_tensor(spacetime, t, x, y, z, &inv_spatial_metric);

  double lapse_function;
  double* shift_vector = gkyl_malloc(sizeof(double[3]));
  tov_lapse_function(spacetime, t, x, y, z, &lapse_function);
  tov_shift_vector(spacetime, t, x, y, z, &shift_vector);

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
tov_spatial_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spatial_metric_det)
{
  double** spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  tov_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);

  *spatial_metric_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
  }
  gkyl_free(spatial_metric);
}

static void
tov_spacetime_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spacetime_metric_det)
{
  double spatial_metric_det;
  double lapse_function;
  tov_spatial_metric_det(spacetime, t, x, y, z, &spatial_metric_det);
  tov_lapse_function(spacetime, t, x, y, z, &lapse_function);

  *spacetime_metric_det = -(lapse_function * lapse_function) * spatial_metric_det;
}

static void
tov_spatial_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_metric_tensor_der)
{
  const struct gr_tov *gr_tov = container_of(spacetime, struct gr_tov, spacetime);

  double r, A, alpha, n[3];
  tov_geometry(spacetime, x, y, z, &r, &A, &alpha, n);
  if (r <= 1.0e-12) {
    for (int k = 0; k < 3; k++)
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          (*spatial_metric_tensor_der)[k][i][j] = 0.0;
    return;
  }

  struct tov_eval_bl bl = { 0 };
  gkyl_tov_eval_bl(gr_tov->tov, r, &bl);

  double comp = 1.0 - (1.0 / A);
  if (comp > 1.0 - 1.0e-12) {
    comp = 1.0 - 1.0e-12;
  }
  double m_prime = 4.0 * M_PI * r * r * bl.e;
  double comp_prime = (2.0 * m_prime / r) - (2.0 * bl.m / (r * r));

  double f = comp;
  double f_prime = comp_prime;
  if (!gr_tov->use_kerr_schild) {
    double one_minus_comp = 1.0 - comp;
    f = comp / one_minus_comp;
    f_prime = comp_prime / (one_minus_comp * one_minus_comp);
  }

  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        double dik = (i == k) ? 1.0 : 0.0;
        double djk = (j == k) ? 1.0 : 0.0;
        double dninj = ((dik - n[i] * n[k]) * n[j] + n[i] * (djk - n[j] * n[k])) / r;
        (*spatial_metric_tensor_der)[k][i][j] = f_prime * n[k] * n[i] * n[j] + f * dninj;
      }
    }
  }
}

static void
tov_spacetime_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_metric_tensor_der)
{
  gkyl_gr_spacetime_metric_tensor_diff(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_metric_tensor_der);
}

static void
tov_lapse_function_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** lapse_function_der)
{
  const struct gr_tov *gr_tov = container_of(spacetime, struct gr_tov, spacetime);

  double r, A, alpha_bl, n[3];
  tov_geometry(spacetime, x, y, z, &r, &A, &alpha_bl, n);
  if (r <= 1.0e-12) {
    for (int i = 0; i < 3; i++) {
      (*lapse_function_der)[i] = 0.0;
    }
    return;
  }

  struct tov_eval_bl bl = { 0 };
  gkyl_tov_eval_bl(gr_tov->tov, r, &bl);

  double comp = 1.0 - (1.0 / A);
  if (comp > 1.0 - 1.0e-12) {
    comp = 1.0 - 1.0e-12;
  }
  double m_prime = 4.0 * M_PI * r * r * bl.e;
  double comp_prime = (2.0 * m_prime / r) - (2.0 * bl.m / (r * r));
  double phi_prime = 0.0;
  if (r * (r - 2.0 * bl.m) > 0.0) {
    phi_prime = (bl.m + 4.0 * M_PI * r * r * r * bl.P) / (r * (r - 2.0 * bl.m));
  }

  double lapse = alpha_bl;
  double log_lapse_prime = phi_prime;
  if (gr_tov->use_kerr_schild) {
    lapse = alpha_bl / sqrt((1.0 - comp) * (1.0 + comp));
    log_lapse_prime += (comp * comp_prime) / ((1.0 - comp) * (1.0 + comp));
  }

  for (int i = 0; i < 3; i++) {
    (*lapse_function_der)[i] = lapse * log_lapse_prime * n[i];
  }
}

static void
tov_shift_vector_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** shift_vector_der)
{
  const struct gr_tov *gr_tov = container_of(spacetime, struct gr_tov, spacetime);

  if (!gr_tov->use_kerr_schild) {
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
        (*shift_vector_der)[j][i] = 0.0;
    return;
  }

  double r, A, alpha_bl, n[3];
  tov_geometry(spacetime, x, y, z, &r, &A, &alpha_bl, n);
  if (r <= 1.0e-12) {
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
        (*shift_vector_der)[j][i] = 0.0;
    return;
  }

  struct tov_eval_bl bl = { 0 };
  gkyl_tov_eval_bl(gr_tov->tov, r, &bl);

  double comp = 1.0 - (1.0 / A);
  if (comp <= 1.0e-14) {
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
        (*shift_vector_der)[j][i] = 0.0;
    return;
  }
  if (comp > 1.0 - 1.0e-12) {
    comp = 1.0 - 1.0e-12;
  }

  double m_prime = 4.0 * M_PI * r * r * bl.e;
  double comp_prime = (2.0 * m_prime / r) - (2.0 * bl.m / (r * r));
  double phi_prime = 0.0;
  if (r * (r - 2.0 * bl.m) > 0.0) {
    phi_prime = (bl.m + 4.0 * M_PI * r * r * r * bl.P) / (r * (r - 2.0 * bl.m));
  }

  double beta_rad = alpha_bl * comp / ((1.0 + comp) * sqrt(1.0 - comp));
  double beta_rad_prime = beta_rad *
    (phi_prime + (comp_prime / comp) - (comp_prime / (1.0 + comp)) + (0.5 * comp_prime / (1.0 - comp)));

  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      double dij = (i == j) ? 1.0 : 0.0;
      (*shift_vector_der)[j][i] = beta_rad_prime * n[j] * n[i] + beta_rad * (dij - n[i] * n[j]) / r;
    }
  }
}

static void
tov_spatial_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_christoffel)
{
  gkyl_gr_spatial_christoffel_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_christoffel);
}

static void
tov_spacetime_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_christoffel)
{
  gkyl_gr_spacetime_christoffel_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_christoffel);
}

static void
tov_spatial_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_riemann_tensor)
{
  gkyl_gr_spatial_riemann_tensor_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_riemann_tensor);
}

static void
tov_spacetime_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_riemann_tensor)
{
  gkyl_gr_spacetime_riemann_tensor_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_riemann_tensor);
}

static void
tov_spatial_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** spatial_ricci_tensor)
{
  gkyl_gr_spatial_ricci_tensor_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_ricci_tensor);
}

static void
tov_spacetime_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double*** spacetime_ricci_tensor)
{
  gkyl_gr_spacetime_ricci_tensor_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_ricci_tensor);
}

static void
tov_spatial_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double* spatial_ricci_scalar)
{
  gkyl_gr_spatial_ricci_scalar_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_ricci_scalar);
}

static void
tov_spacetime_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double* spacetime_ricci_scalar)
{
  gkyl_gr_spacetime_ricci_scalar_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_ricci_scalar);
}

static void
tov_spatial_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_weyl_tensor)
{
  gkyl_gr_spatial_weyl_tensor_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_weyl_tensor);
}

static void
tov_spacetime_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_weyl_tensor)
{
  gkyl_gr_spacetime_weyl_tensor_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_weyl_tensor);
}

static void
tov_extrinsic_curvature_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** extrinsic_curvature_tensor)
{
  const struct gr_tov *gr_tov = container_of(spacetime, struct gr_tov, spacetime);

  if (!gr_tov->use_kerr_schild) {
    // Areal (static) gauge: beta^i = 0, d_t gamma_ij = 0 => K_ij = 0.
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        (*extrinsic_curvature_tensor)[i][j] = 0.0;
      }
    }
    return;
  }

  // Kerr-Schild: stationary radial slicing with beta^r != 0. Use the analytic
  // spherical expression for K_ij instead of differentiating the Cartesian
  // n_i n_j representation with a tiny finite-difference step.
  double r, A, alpha_bl, n[3];
  tov_geometry(spacetime, x, y, z, &r, &A, &alpha_bl, n);
  if (r <= 1.0e-12) {
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        (*extrinsic_curvature_tensor)[i][j] = 0.0;
    return;
  }

  struct tov_eval_bl bl = { 0 };
  gkyl_tov_eval_bl(gr_tov->tov, r, &bl);

  double comp = 1.0 - (1.0 / A);
  if (comp <= 1.0e-14) {
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        (*extrinsic_curvature_tensor)[i][j] = 0.0;
    return;
  }
  if (comp > 1.0 - 1.0e-12) {
    comp = 1.0 - 1.0e-12;
  }

  double m_prime = 4.0 * M_PI * r * r * bl.e;
  double comp_prime = (2.0 * m_prime / r) - (2.0 * bl.m / (r * r));
  double phi_prime = 0.0;
  if (r * (r - 2.0 * bl.m) > 0.0) {
    phi_prime = (bl.m + 4.0 * M_PI * r * r * r * bl.P) / (r * (r - 2.0 * bl.m));
  }

  double lapse = alpha_bl / sqrt((1.0 - comp) * (1.0 + comp));
  double beta_rad = alpha_bl * comp / ((1.0 + comp) * sqrt(1.0 - comp));
  double beta_rad_prime = beta_rad *
    (phi_prime + (comp_prime / comp) - (comp_prime / (1.0 + comp)) + (0.5 * comp_prime / (1.0 - comp)));

  double B_rad = 1.0 + comp;
  double B_rad_prime = comp_prime;
  double K_rr = ((beta_rad * B_rad_prime) + (2.0 * B_rad * beta_rad_prime)) / (2.0 * lapse);
  double K_tan = beta_rad / (lapse * r);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double delta = (i == j) ? 1.0 : 0.0;
      (*extrinsic_curvature_tensor)[i][j] = K_tan * (delta - n[i] * n[j]) + K_rr * n[i] * n[j];
    }
  }
}

static void
tov_conformal_factor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* conformal_factor)
{
  double spatial_metric_det;
  tov_spatial_metric_det(spacetime, t, x, y, z, &spatial_metric_det);

  *conformal_factor = pow(spatial_metric_det, 1.0 / 12.0);
}

static void
tov_bssn_conformal_factor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* bssn_conformal_factor)
{
  double spatial_metric_det;
  tov_spatial_metric_det(spacetime, t, x, y, z, &spatial_metric_det);

  *bssn_conformal_factor = 1.0 / pow(spatial_metric_det, 1.0 / 6.0);
}

static void
tov_conformal_derivatives(const struct gkyl_gr_spacetime* spacetime,
  double x, double y, double z, bool bssn, double gradient[3], double hessian[3][3])
{
  const struct gr_tov *gr_tov = container_of(spacetime, struct gr_tov, spacetime);
  double r, A, alpha, n[3];
  tov_geometry(spacetime, x, y, z, &r, &A, &alpha, n);
  struct tov_eval_bl bl = {0};
  gkyl_tov_eval_bl(gr_tov->tov, r, &bl);
  double first = 0.0, second, first_over_r;
  if (r <= 1.0e-12) {
    // m = 4*pi*e_c*r^3/3 gives the regular Cartesian Hessian at the center.
    second = (bssn ? -8.0 : 4.0)*M_PI*bl.e/9.0;
    first_over_r = second;
  }
  else {
    double mp, mpp;
    if (!gkyl_tov_eval_mass_derivs(gr_tov->tov, r, &mp, &mpp)) {
      for (int i=0; i<3; ++i) {
        gradient[i] = NAN;
        for (int j=0; j<3; ++j) hessian[i][j] = NAN;
      }
      return;
    }
    double c = 1.0 - 1.0/A;
    double cp = 2.0*mp/r - 2.0*bl.m/(r*r);
    double cpp = 2.0*mpp/r - 4.0*mp/(r*r) + 4.0*bl.m/(r*r*r);
    // KS: det(gamma)=1+c. Areal: det(gamma)=1/(1-c).
    double sign = gr_tov->use_kerr_schild ? 1.0 : -1.0;
    double exponent = sign*(bssn ? -1.0/6.0 : 1.0/12.0);
    double f = 1.0 + sign*c, fp = sign*cp, fpp = sign*cpp;
    double factor = exponent*pow(f, exponent-1.0);
    first = factor*fp;
    second = factor*(fpp + (exponent-1.0)*fp*fp/f);
    first_over_r = first/r;
  }
  for (int i=0; i<3; ++i) {
    gradient[i] = first*n[i];
    for (int j=0; j<3; ++j)
      hessian[i][j] = (second-first_over_r)*n[i]*n[j] + (i==j ? first_over_r : 0.0);
  }
}

static void
tov_conformal_factor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** conformal_factor_der)
{
  double gradient[3], hessian[3][3];
  tov_conformal_derivatives(spacetime, x, y, z, false, gradient, hessian);
  for (int i=0; i<3; ++i) (*conformal_factor_der)[i] = gradient[i];
}

static void
tov_bssn_conformal_factor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** bssn_conformal_factor_der)
{
  double gradient[3], hessian[3][3];
  tov_conformal_derivatives(spacetime, x, y, z, true, gradient, hessian);
  for (int i=0; i<3; ++i) (*bssn_conformal_factor_der)[i] = gradient[i];
}

static void
tov_conformal_factor_der2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** conformal_factor_der2)
{
  double gradient[3], hessian[3][3];
  tov_conformal_derivatives(spacetime, x, y, z, false, gradient, hessian);
  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j)
    (*conformal_factor_der2)[i][j] = hessian[i][j];
}

static void
tov_bssn_conformal_factor_der2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** bssn_conformal_factor_der2)
{
  double gradient[3], hessian[3][3];
  tov_conformal_derivatives(spacetime, x, y, z, true, gradient, hessian);
  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j)
    (*bssn_conformal_factor_der2)[i][j] = hessian[i][j];
}

static void
tov_excision_region(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* in_excision_region)
{
  *in_excision_region = false;
}

void
gkyl_gr_tov_spacetime_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_gr_spacetime* base = container_of(ref, struct gkyl_gr_spacetime, ref_count);

  if (gkyl_gr_spacetime_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct gr_tov *gr_tov = container_of(base->on_dev, struct gr_tov, spacetime);
    gkyl_cu_free(gr_tov);
  }

  struct gr_tov *gr_tov = container_of(base, struct gr_tov, spacetime);
  gkyl_free(gr_tov);
}

struct gkyl_gr_spacetime*
gkyl_gr_tov_spacetime_new(bool use_gpu, const struct gkyl_tov *tov, double pos_x, double pos_y, double pos_z, bool use_kerr_schild)
{
  return gkyl_gr_tov_spacetime_inew(&(struct gkyl_gr_tov_inp) {
      .use_gpu = use_gpu,
      .tov = tov,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_kerr_schild = use_kerr_schild,
    }
  );
}

struct gkyl_gr_spacetime*
gkyl_gr_tov_spacetime_inew(const struct gkyl_gr_tov_inp* inp)
{
  struct gr_tov *gr_tov = gkyl_malloc(sizeof(struct gr_tov));

  gr_tov->tov = inp->tov;
  gr_tov->pos_x = inp->pos_x;
  gr_tov->pos_y = inp->pos_y;
  gr_tov->pos_z = inp->pos_z;
  gr_tov->use_kerr_schild = inp->use_kerr_schild;

  gr_tov->spacetime.spatial_metric_tensor_func = tov_spatial_metric_tensor;
  gr_tov->spacetime.spacetime_metric_tensor_func = tov_spacetime_metric_tensor;

  gr_tov->spacetime.spatial_inv_metric_tensor_func = tov_spatial_inv_metric_tensor;
  gr_tov->spacetime.spacetime_inv_metric_tensor_func = tov_spacetime_inv_metric_tensor;

  gr_tov->spacetime.spatial_metric_det_func = tov_spatial_metric_det;
  gr_tov->spacetime.spacetime_metric_det_func = tov_spacetime_metric_det;

  gr_tov->spacetime.spatial_metric_tensor_der_func = tov_spatial_metric_tensor_der;
  gr_tov->spacetime.spacetime_metric_tensor_der_func = tov_spacetime_metric_tensor_der;

  gr_tov->spacetime.lapse_function_func = tov_lapse_function;
  gr_tov->spacetime.shift_vector_func = tov_shift_vector;

  gr_tov->spacetime.lapse_function_der_func = tov_lapse_function_der;
  gr_tov->spacetime.shift_vector_der_func = tov_shift_vector_der;

  gr_tov->spacetime.spatial_christoffel_func = tov_spatial_christoffel;
  gr_tov->spacetime.spacetime_christoffel_func = tov_spacetime_christoffel;

  gr_tov->spacetime.spatial_riemann_tensor_func = tov_spatial_riemann_tensor;
  gr_tov->spacetime.spacetime_riemann_tensor_func = tov_spacetime_riemann_tensor;

  gr_tov->spacetime.spatial_ricci_tensor_func = tov_spatial_ricci_tensor;
  gr_tov->spacetime.spacetime_ricci_tensor_func = tov_spacetime_ricci_tensor;

  gr_tov->spacetime.spatial_ricci_scalar_func = tov_spatial_ricci_scalar;
  gr_tov->spacetime.spacetime_ricci_scalar_func = tov_spacetime_ricci_scalar;

  gr_tov->spacetime.spatial_weyl_tensor_func = tov_spatial_weyl_tensor;
  gr_tov->spacetime.spacetime_weyl_tensor_func = tov_spacetime_weyl_tensor;

  gr_tov->spacetime.extrinsic_curvature_tensor_func = tov_extrinsic_curvature_tensor;

  gr_tov->spacetime.conformal_factor_func = tov_conformal_factor;
  gr_tov->spacetime.bssn_conformal_factor_func = tov_bssn_conformal_factor;

  gr_tov->spacetime.conformal_factor_der_func = tov_conformal_factor_der;
  gr_tov->spacetime.bssn_conformal_factor_der_func = tov_bssn_conformal_factor_der;

  gr_tov->spacetime.conformal_factor_der2_func = tov_conformal_factor_der2;
  gr_tov->spacetime.bssn_conformal_factor_der2_func = tov_bssn_conformal_factor_der2;

  gr_tov->spacetime.excision_region_func = tov_excision_region;

  gr_tov->spacetime.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_tov->spacetime.flags);
  gr_tov->spacetime.ref_count = gkyl_ref_count_init(gkyl_gr_tov_spacetime_free);
  gr_tov->spacetime.on_dev = &gr_tov->spacetime; // On the CPU, the spacetime object points to itself.

  return &gr_tov->spacetime;
}
