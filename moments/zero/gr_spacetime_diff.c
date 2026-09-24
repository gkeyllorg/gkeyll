#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gr_spacetime_diff.h>

void
gkyl_gr_spatial_metric_tensor_diff(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
   const double dx, const double dy, const double dz, double**** spatial_metric_tensor_diff)
{
  double **spatial_metric_x_forward = gkyl_malloc(sizeof(double*[3]));
  double **spatial_metric_y_forward = gkyl_malloc(sizeof(double*[3]));
  double **spatial_metric_z_forward = gkyl_malloc(sizeof(double*[3]));

  double **spatial_metric_x_backward = gkyl_malloc(sizeof(double*[3]));
  double **spatial_metric_y_backward = gkyl_malloc(sizeof(double*[3]));
  double **spatial_metric_z_backward = gkyl_malloc(sizeof(double*[3]));

  for (int i = 0; i < 3; i++) {
    spatial_metric_x_forward[i] = gkyl_malloc(sizeof(double[3]));
    spatial_metric_y_forward[i] = gkyl_malloc(sizeof(double[3]));
    spatial_metric_z_forward[i] = gkyl_malloc(sizeof(double[3]));

    spatial_metric_x_backward[i] = gkyl_malloc(sizeof(double[3]));
    spatial_metric_y_backward[i] = gkyl_malloc(sizeof(double[3]));
    spatial_metric_z_backward[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->spatial_metric_tensor_func(spacetime, t, x + (0.5 * dx), y, z, &spatial_metric_x_forward);
  spacetime->spatial_metric_tensor_func(spacetime, t, x, y + (0.5 * dy), z, &spatial_metric_y_forward);
  spacetime->spatial_metric_tensor_func(spacetime, t, x, y, z + (0.5 * dz), &spatial_metric_z_forward);

  spacetime->spatial_metric_tensor_func(spacetime, t, x - (0.5 * dx), y, z, &spatial_metric_x_backward);
  spacetime->spatial_metric_tensor_func(spacetime, t, x, y - (0.5 * dy), z, &spatial_metric_y_backward);
  spacetime->spatial_metric_tensor_func(spacetime, t, x, y, z - (0.5 * dz), &spatial_metric_z_backward);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*spatial_metric_tensor_diff)[0][i][j] = (1.0 / dx) * (spatial_metric_x_forward[i][j] - spatial_metric_x_backward[i][j]);
      (*spatial_metric_tensor_diff)[1][i][j] = (1.0 / dy) * (spatial_metric_y_forward[i][j] - spatial_metric_y_backward[i][j]);
      (*spatial_metric_tensor_diff)[2][i][j] = (1.0 / dz) * (spatial_metric_z_forward[i][j] - spatial_metric_z_backward[i][j]);
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric_x_forward[i]);
    gkyl_free(spatial_metric_y_forward[i]);
    gkyl_free(spatial_metric_z_forward[i]);

    gkyl_free(spatial_metric_x_backward[i]);
    gkyl_free(spatial_metric_y_backward[i]);
    gkyl_free(spatial_metric_z_backward[i]);
  }
  gkyl_free(spatial_metric_x_forward);
  gkyl_free(spatial_metric_y_forward);
  gkyl_free(spatial_metric_z_forward);

  gkyl_free(spatial_metric_x_backward);
  gkyl_free(spatial_metric_y_backward);
  gkyl_free(spatial_metric_z_backward);
}

void
gkyl_gr_spacetime_metric_tensor_diff(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_metric_tensor_diff)
{
  double **spacetime_metric_t_forward = gkyl_malloc(sizeof(double*[4]));
  double **spacetime_metric_x_forward = gkyl_malloc(sizeof(double*[4]));
  double **spacetime_metric_y_forward = gkyl_malloc(sizeof(double*[4]));
  double **spacetime_metric_z_forward = gkyl_malloc(sizeof(double*[4]));

  double **spacetime_metric_t_backward = gkyl_malloc(sizeof(double*[4]));
  double **spacetime_metric_x_backward = gkyl_malloc(sizeof(double*[4]));
  double **spacetime_metric_y_backward = gkyl_malloc(sizeof(double*[4]));
  double **spacetime_metric_z_backward = gkyl_malloc(sizeof(double*[4]));

  for (int i = 0; i < 4; i++) {
    spacetime_metric_t_forward[i] = gkyl_malloc(sizeof(double[4]));
    spacetime_metric_x_forward[i] = gkyl_malloc(sizeof(double[4]));
    spacetime_metric_y_forward[i] = gkyl_malloc(sizeof(double[4]));
    spacetime_metric_z_forward[i] = gkyl_malloc(sizeof(double[4]));
    
    spacetime_metric_t_backward[i] = gkyl_malloc(sizeof(double[4]));
    spacetime_metric_x_backward[i] = gkyl_malloc(sizeof(double[4]));
    spacetime_metric_y_backward[i] = gkyl_malloc(sizeof(double[4]));
    spacetime_metric_z_backward[i] = gkyl_malloc(sizeof(double[4]));
  }

  spacetime->spacetime_metric_tensor_func(spacetime, t + (0.5 * dt), x, y, z, &spacetime_metric_t_forward);
  spacetime->spacetime_metric_tensor_func(spacetime, t, x + (0.5 * dx), y, z, &spacetime_metric_x_forward);
  spacetime->spacetime_metric_tensor_func(spacetime, t, x, y + (0.5 * dy), z, &spacetime_metric_y_forward);
  spacetime->spacetime_metric_tensor_func(spacetime, t, x, y, z + (0.5 * dz), &spacetime_metric_z_forward);

  spacetime->spacetime_metric_tensor_func(spacetime, t - (0.5 * dt), x, y, z, &spacetime_metric_t_backward);
  spacetime->spacetime_metric_tensor_func(spacetime, t, x - (0.5 * dx), y, z, &spacetime_metric_x_backward);
  spacetime->spacetime_metric_tensor_func(spacetime, t, x, y - (0.5 * dy), z, &spacetime_metric_y_backward);
  spacetime->spacetime_metric_tensor_func(spacetime, t, x, y, z - (0.5 * dz), &spacetime_metric_z_backward);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      (*spacetime_metric_tensor_diff)[0][i][j] = (1.0 / dt) * (spacetime_metric_t_forward[i][j] - spacetime_metric_t_backward[i][j]);
      (*spacetime_metric_tensor_diff)[1][i][j] = (1.0 / dx) * (spacetime_metric_x_forward[i][j] - spacetime_metric_x_backward[i][j]);
      (*spacetime_metric_tensor_diff)[2][i][j] = (1.0 / dy) * (spacetime_metric_y_forward[i][j] - spacetime_metric_y_backward[i][j]);
      (*spacetime_metric_tensor_diff)[3][i][j] = (1.0 / dz) * (spacetime_metric_z_forward[i][j] - spacetime_metric_z_backward[i][j]);
    }
  }

  for (int i = 0; i < 4; i++) {
    gkyl_free(spacetime_metric_t_forward[i]);
    gkyl_free(spacetime_metric_x_forward[i]);
    gkyl_free(spacetime_metric_y_forward[i]);
    gkyl_free(spacetime_metric_z_forward[i]);

    gkyl_free(spacetime_metric_t_backward[i]);
    gkyl_free(spacetime_metric_x_backward[i]);
    gkyl_free(spacetime_metric_y_backward[i]);
    gkyl_free(spacetime_metric_z_backward[i]);
  }
  gkyl_free(spacetime_metric_t_forward);
  gkyl_free(spacetime_metric_x_forward);
  gkyl_free(spacetime_metric_y_forward);
  gkyl_free(spacetime_metric_z_forward);

  gkyl_free(spacetime_metric_t_backward);
  gkyl_free(spacetime_metric_x_backward);
  gkyl_free(spacetime_metric_y_backward);
  gkyl_free(spacetime_metric_z_backward);
}

void
gkyl_gr_lapse_function_diff(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** lapse_function_diff)
{
  double lapse_function_x_forward;
  double lapse_function_y_forward;
  double lapse_function_z_forward;

  double lapse_function_x_backward;
  double lapse_function_y_backward;
  double lapse_function_z_backward;

  spacetime->lapse_function_func(spacetime, t, x + (0.5 * dx), y, z, &lapse_function_x_forward);
  spacetime->lapse_function_func(spacetime, t, x, y + (0.5 * dy), z, &lapse_function_y_forward);
  spacetime->lapse_function_func(spacetime, t, x, y, z + (0.5 * dz), &lapse_function_z_forward);

  spacetime->lapse_function_func(spacetime, t, x - (0.5 * dx), y, z, &lapse_function_x_backward);
  spacetime->lapse_function_func(spacetime, t, x, y - (0.5 * dy), z, &lapse_function_y_backward);
  spacetime->lapse_function_func(spacetime, t, x, y, z - (0.5 * dz), &lapse_function_z_backward);

  (*lapse_function_diff)[0] = (1.0 / dx) * (lapse_function_x_forward - lapse_function_x_backward);
  (*lapse_function_diff)[1] = (1.0 / dy) * (lapse_function_y_forward - lapse_function_y_backward);
  (*lapse_function_diff)[2] = (1.0 / dz) * (lapse_function_z_forward - lapse_function_z_backward);
}

void
gkyl_gr_shift_vector_diff(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** shift_vector_diff)
{
  double* shift_vector_x_forward = gkyl_malloc(sizeof(double[3]));
  double* shift_vector_y_forward = gkyl_malloc(sizeof(double[3]));
  double* shift_vector_z_forward = gkyl_malloc(sizeof(double[3]));
  
  double* shift_vector_x_backward = gkyl_malloc(sizeof(double[3]));
  double* shift_vector_y_backward = gkyl_malloc(sizeof(double[3]));
  double* shift_vector_z_backward = gkyl_malloc(sizeof(double[3]));

  spacetime->shift_vector_func(spacetime, t, x + (0.5 * dx), y, z, &shift_vector_x_forward);
  spacetime->shift_vector_func(spacetime, t, x, y + (0.5 * dy), z, &shift_vector_y_forward);
  spacetime->shift_vector_func(spacetime, t, x, y, z + (0.5 * dz), &shift_vector_z_forward);

  spacetime->shift_vector_func(spacetime, t, x - (0.5 * dx), y, z, &shift_vector_x_backward);
  spacetime->shift_vector_func(spacetime, t, x, y - (0.5 * dy), z, &shift_vector_y_backward);
  spacetime->shift_vector_func(spacetime, t, x, y, z - (0.5 * dz), &shift_vector_z_backward);
  
  for (int i = 0; i < 3; i++) {
    (*shift_vector_diff)[0][i] = (1.0 / dx) * (shift_vector_x_forward[i] - shift_vector_x_backward[i]);
    (*shift_vector_diff)[1][i] = (1.0 / dy) * (shift_vector_y_forward[i] - shift_vector_y_backward[i]);
    (*shift_vector_diff)[2][i] = (1.0 / dz) * (shift_vector_z_forward[i] - shift_vector_z_backward[i]);
  }

  gkyl_free(shift_vector_x_forward);
  gkyl_free(shift_vector_y_forward);
  gkyl_free(shift_vector_z_forward);

  gkyl_free(shift_vector_x_backward);
  gkyl_free(shift_vector_y_backward);
  gkyl_free(shift_vector_z_backward);
}

void
gkyl_gr_spatial_christoffel_fd(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_christoffel)
{
  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

    for (int j = 0; j < 3; j++) {
      spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
    }
  }

  spacetime->spatial_inv_metric_tensor_func(spacetime, t, x, y, z, &inv_spatial_metric);
  spacetime->spatial_metric_tensor_der_func(spacetime, t, x, y, z, dx, dy, dz, &spatial_metric_der);

  for (int i = 0; i < 3; i++) {
    for (int j = 0 ; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        (*spatial_christoffel)[i][j][k] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          (*spatial_christoffel)[i][j][k] += (0.5 * inv_spatial_metric[i][l]) * (spatial_metric_der[k][l][j] + spatial_metric_der[j][l][k]
            - spatial_metric_der[l][j][k]);
        }
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      gkyl_free(spatial_metric_der[i][j]);
    }

    gkyl_free(inv_spatial_metric[i]);
    gkyl_free(spatial_metric_der[i]);
  }
  gkyl_free(inv_spatial_metric);
  gkyl_free(spatial_metric_der);
}

void
gkyl_gr_spacetime_christoffel_fd(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_christoffel)
{
  double **inv_spacetime_metric = gkyl_malloc(sizeof(double*[4]));
  for (int i = 0; i < 4; i++) {
    inv_spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
  }

  double ***spacetime_metric_der = gkyl_malloc(sizeof(double**[4]));
  for (int i = 0; i < 4; i++) {
    spacetime_metric_der[i] = gkyl_malloc(sizeof(double*[4]));

    for (int j = 0; j < 4; j++) {
      spacetime_metric_der[i][j] = gkyl_malloc(sizeof(double[4]));
    }
  }

  spacetime->spacetime_inv_metric_tensor_func(spacetime, t, x, y, z, &inv_spacetime_metric);
  spacetime->spacetime_metric_tensor_der_func(spacetime, t, x, y, z, dt, dx, dy, dz, &spacetime_metric_der);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        (*spacetime_christoffel)[i][j][k] = 0.0;
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        for (int l = 0; l < 4; l++) {
          (*spacetime_christoffel)[i][j][k] += (0.5 * inv_spacetime_metric[i][l]) * (spacetime_metric_der[k][l][j] + spacetime_metric_der[j][l][k]
            - spacetime_metric_der[l][j][k]);
        }
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      gkyl_free(spacetime_metric_der[i][j]);
    }

    gkyl_free(inv_spacetime_metric[i]);
    gkyl_free(spacetime_metric_der[i]);
  }
  gkyl_free(inv_spacetime_metric);
  gkyl_free(spacetime_metric_der);
}

void
gkyl_gr_spatial_riemann_tensor_fd(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_riemann_tensor)
{
  double ***spatial_christoffel = gkyl_malloc(sizeof(double**[3]));
  double ***spatial_christoffel_x_forward = gkyl_malloc(sizeof(double**[3]));
  double ***spatial_christoffel_y_forward = gkyl_malloc(sizeof(double**[3]));
  double ***spatial_christoffel_z_forward = gkyl_malloc(sizeof(double**[3]));

  double ***spatial_christoffel_x_backward = gkyl_malloc(sizeof(double**[3]));
  double ***spatial_christoffel_y_backward = gkyl_malloc(sizeof(double**[3]));
  double ***spatial_christoffel_z_backward = gkyl_malloc(sizeof(double**[3]));

  for (int i = 0; i < 3; i++) {
    spatial_christoffel[i] = gkyl_malloc(sizeof(double*[3]));
    spatial_christoffel_x_forward[i] = gkyl_malloc(sizeof(double*[3]));
    spatial_christoffel_y_forward[i] = gkyl_malloc(sizeof(double*[3]));
    spatial_christoffel_z_forward[i] = gkyl_malloc(sizeof(double*[3]));

    spatial_christoffel_x_backward[i] = gkyl_malloc(sizeof(double*[3]));
    spatial_christoffel_y_backward[i] = gkyl_malloc(sizeof(double*[3]));
    spatial_christoffel_z_backward[i] = gkyl_malloc(sizeof(double*[3]));

    for (int j = 0; j < 3; j++) {
      spatial_christoffel[i][j] = gkyl_malloc(sizeof(double[3]));
      spatial_christoffel_x_forward[i][j] = gkyl_malloc(sizeof(double[3]));
      spatial_christoffel_y_forward[i][j] = gkyl_malloc(sizeof(double[3]));
      spatial_christoffel_z_forward[i][j] = gkyl_malloc(sizeof(double[3]));

      spatial_christoffel_x_backward[i][j] = gkyl_malloc(sizeof(double[3]));
      spatial_christoffel_y_backward[i][j] = gkyl_malloc(sizeof(double[3]));
      spatial_christoffel_z_backward[i][j] = gkyl_malloc(sizeof(double[3]));
    }
  }
  
  spacetime->spatial_christoffel_func(spacetime, t, x, y, z, dx, dy, dz, &spatial_christoffel);
  spacetime->spatial_christoffel_func(spacetime, t, x + (0.5 * dx), y, z, dx, dy, dz, &spatial_christoffel_x_forward);
  spacetime->spatial_christoffel_func(spacetime, t, x, y + (0.5 * dy), z, dx, dy, dz, &spatial_christoffel_y_forward);
  spacetime->spatial_christoffel_func(spacetime, t, x, y, z + (0.5 * dz), dx, dy, dz, &spatial_christoffel_z_forward);

  spacetime->spatial_christoffel_func(spacetime, t, x - (0.5 * dx), y, z, dx, dy, dz, &spatial_christoffel_x_backward);
  spacetime->spatial_christoffel_func(spacetime, t, x, y - (0.5 * dy), z, dx, dy, dz, &spatial_christoffel_y_backward);
  spacetime->spatial_christoffel_func(spacetime, t, x, y, z - (0.5 * dz), dx, dy, dz, &spatial_christoffel_z_backward);

  double ****spatial_christoffel_der = gkyl_malloc(sizeof(double***[3]));
  for (int i = 0; i < 3; i++) {
    spatial_christoffel_der[i] = gkyl_malloc(sizeof(double**[3]));

    for (int j = 0; j < 3; j++) {
      spatial_christoffel_der[i][j] = gkyl_malloc(sizeof(double*[3]));

      for (int k = 0; k < 3; k++) {
        spatial_christoffel_der[i][j][k] = gkyl_malloc(sizeof(double[3]));
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_christoffel_der[0][i][j][k] = (1.0 / dx) * (spatial_christoffel_x_forward[i][j][k] - spatial_christoffel_x_backward[i][j][k]);
        spatial_christoffel_der[1][i][j][k] = (1.0 / dy) * (spatial_christoffel_y_forward[i][j][k] - spatial_christoffel_y_backward[i][j][k]);
        spatial_christoffel_der[2][i][j][k] = (1.0 / dz) * (spatial_christoffel_z_forward[i][j][k] - spatial_christoffel_z_backward[i][j][k]);
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          (*spatial_riemann_tensor)[i][j][k][l] = spatial_christoffel_der[k][i][l][j] - spatial_christoffel_der[l][i][k][j];

          for (int m = 0; m < 3; m++) {
            (*spatial_riemann_tensor)[i][j][k][l] += (spatial_christoffel[i][k][m] * spatial_christoffel[m][l][j]);
            (*spatial_riemann_tensor)[i][j][k][l] -= (spatial_christoffel[i][l][m] * spatial_christoffel[m][k][j]);
          }
        }
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        gkyl_free(spatial_christoffel_der[i][j][k]);
      }

      gkyl_free(spatial_christoffel_der[i][j]);
      gkyl_free(spatial_christoffel[i][j]);
      gkyl_free(spatial_christoffel_x_forward[i][j]);
      gkyl_free(spatial_christoffel_y_forward[i][j]);
      gkyl_free(spatial_christoffel_z_forward[i][j]);
      gkyl_free(spatial_christoffel_x_backward[i][j]);
      gkyl_free(spatial_christoffel_y_backward[i][j]);
      gkyl_free(spatial_christoffel_z_backward[i][j]);
    }

    gkyl_free(spatial_christoffel_der[i]);
    gkyl_free(spatial_christoffel[i]);
    gkyl_free(spatial_christoffel_x_forward[i]);
    gkyl_free(spatial_christoffel_y_forward[i]);
    gkyl_free(spatial_christoffel_z_forward[i]);
    gkyl_free(spatial_christoffel_x_backward[i]);
    gkyl_free(spatial_christoffel_y_backward[i]);
    gkyl_free(spatial_christoffel_z_backward[i]);
  }

  gkyl_free(spatial_christoffel_der);
  gkyl_free(spatial_christoffel);
  gkyl_free(spatial_christoffel_x_forward);
  gkyl_free(spatial_christoffel_y_forward);
  gkyl_free(spatial_christoffel_z_forward);
  gkyl_free(spatial_christoffel_x_backward);
  gkyl_free(spatial_christoffel_y_backward);
  gkyl_free(spatial_christoffel_z_backward);
}

void
gkyl_gr_spacetime_riemann_tensor_fd(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_riemann_tensor)
{
  double ***spacetime_christoffel = gkyl_malloc(sizeof(double**[4]));
  double ***spacetime_christoffel_t_forward = gkyl_malloc(sizeof(double**[4]));
  double ***spacetime_christoffel_x_forward = gkyl_malloc(sizeof(double**[4]));
  double ***spacetime_christoffel_y_forward = gkyl_malloc(sizeof(double**[4]));
  double ***spacetime_christoffel_z_forward = gkyl_malloc(sizeof(double**[4]));

  double ***spacetime_christoffel_t_backward = gkyl_malloc(sizeof(double**[4]));
  double ***spacetime_christoffel_x_backward = gkyl_malloc(sizeof(double**[4]));
  double ***spacetime_christoffel_y_backward = gkyl_malloc(sizeof(double**[4]));
  double ***spacetime_christoffel_z_backward = gkyl_malloc(sizeof(double**[4]));

  for (int i = 0; i < 4; i++) {
    spacetime_christoffel[i] = gkyl_malloc(sizeof(double*[4]));
    spacetime_christoffel_t_forward[i] = gkyl_malloc(sizeof(double*[4]));
    spacetime_christoffel_x_forward[i] = gkyl_malloc(sizeof(double*[4]));
    spacetime_christoffel_y_forward[i] = gkyl_malloc(sizeof(double*[4]));
    spacetime_christoffel_z_forward[i] = gkyl_malloc(sizeof(double*[4]));

    spacetime_christoffel_t_backward[i] = gkyl_malloc(sizeof(double*[4]));
    spacetime_christoffel_x_backward[i] = gkyl_malloc(sizeof(double*[4]));
    spacetime_christoffel_y_backward[i] = gkyl_malloc(sizeof(double*[4]));
    spacetime_christoffel_z_backward[i] = gkyl_malloc(sizeof(double*[4]));

    for (int j = 0; j < 4; j++) {
      spacetime_christoffel[i][j] = gkyl_malloc(sizeof(double[4]));
      spacetime_christoffel_t_forward[i][j] = gkyl_malloc(sizeof(double[4]));
      spacetime_christoffel_x_forward[i][j] = gkyl_malloc(sizeof(double[4]));
      spacetime_christoffel_y_forward[i][j] = gkyl_malloc(sizeof(double[4]));
      spacetime_christoffel_z_forward[i][j] = gkyl_malloc(sizeof(double[4]));

      spacetime_christoffel_t_backward[i][j] = gkyl_malloc(sizeof(double[4]));
      spacetime_christoffel_x_backward[i][j] = gkyl_malloc(sizeof(double[4]));
      spacetime_christoffel_y_backward[i][j] = gkyl_malloc(sizeof(double[4]));
      spacetime_christoffel_z_backward[i][j] = gkyl_malloc(sizeof(double[4]));
    }
  }

  spacetime->spacetime_christoffel_func(spacetime, t, x, y, z, dt, dx, dy, dz, &spacetime_christoffel);
  spacetime->spacetime_christoffel_func(spacetime, t + (0.5 * dt), x, y, z, dt, dx, dy, dz, &spacetime_christoffel_t_forward);
  spacetime->spacetime_christoffel_func(spacetime, t, x + (0.5 * dx), y, z, dt, dx, dy, dz, &spacetime_christoffel_x_forward);
  spacetime->spacetime_christoffel_func(spacetime, t, x, y + (0.5 * dy), z, dt, dx, dy, dz, &spacetime_christoffel_y_forward);
  spacetime->spacetime_christoffel_func(spacetime, t, x, y, z + (0.5 * dz), dt, dx, dy, dz, &spacetime_christoffel_z_forward);

  spacetime->spacetime_christoffel_func(spacetime, t - (0.5 * dt), x, y, z, dt, dx, dy, dz, &spacetime_christoffel_t_backward);
  spacetime->spacetime_christoffel_func(spacetime, t, x - (0.5 * dx), y, z, dt, dx, dy, dz, &spacetime_christoffel_x_backward);
  spacetime->spacetime_christoffel_func(spacetime, t, x, y - (0.5 * dy), z, dt, dx, dy, dz, &spacetime_christoffel_y_backward);
  spacetime->spacetime_christoffel_func(spacetime, t, x, y, z - (0.5 * dz), dt, dx, dy, dz, &spacetime_christoffel_z_backward);

  double ****spacetime_christoffel_der = gkyl_malloc(sizeof(double***[4]));
  for (int i = 0; i < 4; i++) {
    spacetime_christoffel_der[i] = gkyl_malloc(sizeof(double**[4]));

    for (int j = 0; j < 4; j++) {
      spacetime_christoffel_der[i][j] = gkyl_malloc(sizeof(double*[4]));

      for (int k = 0; k < 4; k++) {
        spacetime_christoffel_der[i][j][k] = gkyl_malloc(sizeof(double[4]));
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        spacetime_christoffel_der[0][i][j][k] = (1.0 / dt) * (spacetime_christoffel_t_forward[i][j][k] - spacetime_christoffel_t_backward[i][j][k]);
        spacetime_christoffel_der[1][i][j][k] = (1.0 / dx) * (spacetime_christoffel_x_forward[i][j][k] - spacetime_christoffel_x_backward[i][j][k]);
        spacetime_christoffel_der[2][i][j][k] = (1.0 / dy) * (spacetime_christoffel_y_forward[i][j][k] - spacetime_christoffel_y_backward[i][j][k]);
        spacetime_christoffel_der[3][i][j][k] = (1.0 / dz) * (spacetime_christoffel_z_forward[i][j][k] - spacetime_christoffel_z_backward[i][j][k]);
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        for (int l = 0; l < 4; l++) {
          (*spacetime_riemann_tensor)[i][j][k][l] = spacetime_christoffel_der[k][i][l][j] - spacetime_christoffel_der[l][i][k][j];

          for (int m = 0; m < 4; m++) {
            (*spacetime_riemann_tensor)[i][j][k][l] += (spacetime_christoffel[i][k][m] * spacetime_christoffel[m][l][j]);
            (*spacetime_riemann_tensor)[i][j][k][l] -= (spacetime_christoffel[i][l][m] * spacetime_christoffel[m][k][j]);
          }
        }
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        gkyl_free(spacetime_christoffel_der[i][j][k]);
      }

      gkyl_free(spacetime_christoffel_der[i][j]);
      gkyl_free(spacetime_christoffel[i][j]);
      gkyl_free(spacetime_christoffel_t_forward[i][j]);
      gkyl_free(spacetime_christoffel_x_forward[i][j]);
      gkyl_free(spacetime_christoffel_y_forward[i][j]);
      gkyl_free(spacetime_christoffel_z_forward[i][j]);
      gkyl_free(spacetime_christoffel_t_backward[i][j]);
      gkyl_free(spacetime_christoffel_x_backward[i][j]);
      gkyl_free(spacetime_christoffel_y_backward[i][j]);
      gkyl_free(spacetime_christoffel_z_backward[i][j]);
    }

    gkyl_free(spacetime_christoffel_der[i]);
    gkyl_free(spacetime_christoffel[i]);
    gkyl_free(spacetime_christoffel_t_forward[i]);
    gkyl_free(spacetime_christoffel_x_forward[i]);
    gkyl_free(spacetime_christoffel_y_forward[i]);
    gkyl_free(spacetime_christoffel_z_forward[i]);
    gkyl_free(spacetime_christoffel_t_backward[i]);
    gkyl_free(spacetime_christoffel_x_backward[i]);
    gkyl_free(spacetime_christoffel_y_backward[i]);
    gkyl_free(spacetime_christoffel_z_backward[i]);
  }

  gkyl_free(spacetime_christoffel_der);
  gkyl_free(spacetime_christoffel);
  gkyl_free(spacetime_christoffel_t_forward);
  gkyl_free(spacetime_christoffel_x_forward);
  gkyl_free(spacetime_christoffel_y_forward);
  gkyl_free(spacetime_christoffel_z_forward);
  gkyl_free(spacetime_christoffel_t_backward);
  gkyl_free(spacetime_christoffel_x_backward);
  gkyl_free(spacetime_christoffel_y_backward);
  gkyl_free(spacetime_christoffel_z_backward);
}

void
gkyl_gr_spatial_ricci_tensor_fd(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** spatial_ricci_tensor)
{
  double ****spatial_riemann_tensor = gkyl_malloc(sizeof(double***[3]));
  for (int i = 0; i < 3; i++) {
    spatial_riemann_tensor[i] = gkyl_malloc(sizeof(double**[3]));

    for (int j = 0; j < 3; j++) {
      spatial_riemann_tensor[i][j] = gkyl_malloc(sizeof(double*[3]));

      for (int k = 0; k < 3; k++) {
        spatial_riemann_tensor[i][j][k] = gkyl_malloc(sizeof(double[3]));
      }
    }
  }

  spacetime->spatial_riemann_tensor_func(spacetime, t, x, y, z, dx, dy, dz, &spatial_riemann_tensor);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*spatial_ricci_tensor)[i][j] = 0.0;

      for (int k = 0; k < 3; k++) {
        (*spatial_ricci_tensor)[i][j] += spatial_riemann_tensor[k][i][k][j];
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        gkyl_free(spatial_riemann_tensor[i][j][k]);
      }
      
      gkyl_free(spatial_riemann_tensor[i][j]);
    }

    gkyl_free(spatial_riemann_tensor[i]);
  }

  gkyl_free(spatial_riemann_tensor);
}

void
gkyl_gr_spacetime_ricci_tensor_fd(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double*** spacetime_ricci_tensor)
{
  double ****spacetime_riemann_tensor = gkyl_malloc(sizeof(double***[4]));
  for (int i = 0; i < 4; i++) {
    spacetime_riemann_tensor[i] = gkyl_malloc(sizeof(double**[4]));

    for (int j = 0; j < 4; j++) {
      spacetime_riemann_tensor[i][j] = gkyl_malloc(sizeof(double*[4]));

      for (int k = 0; k < 4; k++) {
        spacetime_riemann_tensor[i][j][k] = gkyl_malloc(sizeof(double[4]));
      }
    }
  }

  spacetime->spacetime_riemann_tensor_func(spacetime, t, x, y, z, dt, dx, dy, dz, &spacetime_riemann_tensor);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      (*spacetime_ricci_tensor)[i][j] = 0.0;

      for (int k = 0; k < 4; k++) {
        (*spacetime_ricci_tensor)[i][j] += spacetime_riemann_tensor[k][i][k][j];
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        gkyl_free(spacetime_riemann_tensor[i][j][k]);
      }
      
      gkyl_free(spacetime_riemann_tensor[i][j]);
    }

    gkyl_free(spacetime_riemann_tensor[i]);
  }

  gkyl_free(spacetime_riemann_tensor);
}

void
gkyl_gr_spatial_ricci_scalar_fd(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double* spatial_ricci_scalar)
{
  double **spatial_inv_metric_tensor = gkyl_malloc(sizeof(double*[3]));
  double **spatial_ricci_tensor = gkyl_malloc(sizeof(double*[3]));
  
  for (int i = 0; i < 3; i++) {
    spatial_inv_metric_tensor[i] = gkyl_malloc(sizeof(double[3]));
    spatial_ricci_tensor[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->spatial_inv_metric_tensor_func(spacetime, t, x, y, z, &spatial_inv_metric_tensor);
  spacetime->spatial_ricci_tensor_func(spacetime, t, x, y, z, dx, dy, dz, &spatial_ricci_tensor);

  *spatial_ricci_scalar = 0.0;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      *spatial_ricci_scalar += spatial_inv_metric_tensor[i][j] * spatial_ricci_tensor[i][j];
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_inv_metric_tensor[i]);
    gkyl_free(spatial_ricci_tensor[i]);
  }

  gkyl_free(spatial_inv_metric_tensor);
  gkyl_free(spatial_ricci_tensor);
}

void
gkyl_gr_spacetime_ricci_scalar_fd(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double* spacetime_ricci_scalar)
{
  double **spacetime_inv_metric_tensor = gkyl_malloc(sizeof(double*[4]));
  double **spacetime_ricci_tensor = gkyl_malloc(sizeof(double*[4]));
  
  for (int i = 0; i < 4; i++) {
    spacetime_inv_metric_tensor[i] = gkyl_malloc(sizeof(double[4]));
    spacetime_ricci_tensor[i] = gkyl_malloc(sizeof(double[4]));
  }

  spacetime->spacetime_inv_metric_tensor_func(spacetime, t, x, y, z, &spacetime_inv_metric_tensor);
  spacetime->spacetime_ricci_tensor_func(spacetime, t, x, y, z, dt, dx, dy, dz, &spacetime_ricci_tensor);

  *spacetime_ricci_scalar = 0.0;

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      *spacetime_ricci_scalar += spacetime_inv_metric_tensor[i][j] * spacetime_ricci_tensor[i][j];
    }
  }

  for (int i = 0; i < 4; i++) {
    gkyl_free(spacetime_inv_metric_tensor[i]);
    gkyl_free(spacetime_ricci_tensor[i]);
  }

  gkyl_free(spacetime_inv_metric_tensor);
  gkyl_free(spacetime_ricci_tensor);
}

void
gkyl_gr_spatial_weyl_tensor_fd(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_weyl_tensor)
{
  double ****spatial_riemann_tensor = gkyl_malloc(sizeof(double***[3]));
  double ****covariant_spatial_riemann_tensor = gkyl_malloc(sizeof(double***[3]));

  for (int i = 0; i < 3; i++) {
    spatial_riemann_tensor[i] = gkyl_malloc(sizeof(double**[3]));
    covariant_spatial_riemann_tensor[i] = gkyl_malloc(sizeof(double**[3]));

    for (int j = 0; j < 3; j++) {
      spatial_riemann_tensor[i][j] = gkyl_malloc(sizeof(double*[3]));
      covariant_spatial_riemann_tensor[i][j] = gkyl_malloc(sizeof(double*[3]));

      for (int k = 0; k < 3; k++) {
        spatial_riemann_tensor[i][j][k] = gkyl_malloc(sizeof(double[3]));
        covariant_spatial_riemann_tensor[i][j][k] = gkyl_malloc(sizeof(double[3]));

        for (int l = 0; l < 3; l++) {
          covariant_spatial_riemann_tensor[i][j][k][l] = 0.0;
        }
      }
    }
  }

  double **spatial_metric_tensor = gkyl_malloc(sizeof(double*[3]));
  double **spatial_ricci_tensor = gkyl_malloc(sizeof(double*[3]));

  for (int i = 0; i < 3; i++) {
    spatial_metric_tensor[i] = gkyl_malloc(sizeof(double[3]));
    spatial_ricci_tensor[i] = gkyl_malloc(sizeof(double[3]));
  }

  double spatial_ricci_scalar;

  spacetime->spatial_riemann_tensor_func(spacetime, t, x, y, z, dx, dy, dz, &spatial_riemann_tensor);
  spacetime->spatial_metric_tensor_func(spacetime, t, x, y, z, &spatial_metric_tensor);
  spacetime->spatial_ricci_tensor_func(spacetime, t, x, y, z, dx, dy, dz, &spatial_ricci_tensor);
  spacetime->spatial_ricci_scalar_func(spacetime, t, x, y, z, dx, dy, dz, &spatial_ricci_scalar);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          for (int m = 0; m < 3; m++) {
            covariant_spatial_riemann_tensor[i][j][k][l] += spatial_metric_tensor[i][m] * spatial_riemann_tensor[m][j][k][l];
          }
        }
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          (*spatial_weyl_tensor)[i][j][k][l] = covariant_spatial_riemann_tensor[i][j][k][l] + ((spatial_ricci_tensor[i][l] * spatial_metric_tensor[j][k])
            - (spatial_ricci_tensor[i][k] * spatial_metric_tensor[j][l]) + (spatial_ricci_tensor[j][k] * spatial_metric_tensor[i][l])
            - (spatial_ricci_tensor[j][l] * spatial_metric_tensor[i][k])) + ((0.5 * spatial_ricci_scalar) * ((spatial_metric_tensor[i][k] * spatial_metric_tensor[j][l])
            - (spatial_metric_tensor[i][l] * spatial_metric_tensor[j][k])));
        }
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        gkyl_free(spatial_riemann_tensor[i][j][k]);
        gkyl_free(covariant_spatial_riemann_tensor[i][j][k]);
      }
      gkyl_free(spatial_riemann_tensor[i][j]);
      gkyl_free(covariant_spatial_riemann_tensor[i][j]);
    }
    gkyl_free(spatial_riemann_tensor[i]);
    gkyl_free(covariant_spatial_riemann_tensor[i]);
    gkyl_free(spatial_metric_tensor[i]);
    gkyl_free(spatial_ricci_tensor[i]);
  }
  gkyl_free(spatial_riemann_tensor);
  gkyl_free(covariant_spatial_riemann_tensor);
  gkyl_free(spatial_metric_tensor);
  gkyl_free(spatial_ricci_tensor);
}

void
gkyl_gr_spacetime_weyl_tensor_fd(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_weyl_tensor)
{
  double ****spacetime_riemann_tensor = gkyl_malloc(sizeof(double***[4]));
  double ****covariant_spacetime_riemann_tensor = gkyl_malloc(sizeof(double***[4]));

  for (int i = 0; i < 4; i++) {
    spacetime_riemann_tensor[i] = gkyl_malloc(sizeof(double**[4]));
    covariant_spacetime_riemann_tensor[i] = gkyl_malloc(sizeof(double**[4]));

    for (int j = 0; j < 4; j++) {
      spacetime_riemann_tensor[i][j] = gkyl_malloc(sizeof(double*[4]));
      covariant_spacetime_riemann_tensor[i][j] = gkyl_malloc(sizeof(double*[4]));

      for (int k = 0; k < 4; k++) {
        spacetime_riemann_tensor[i][j][k] = gkyl_malloc(sizeof(double[4]));
        covariant_spacetime_riemann_tensor[i][j][k] = gkyl_malloc(sizeof(double[4]));

        for (int l = 0; l < 4; l++) {
          covariant_spacetime_riemann_tensor[i][j][k][l] = 0.0;
        }
      }
    }
  }

  double **spacetime_metric_tensor = gkyl_malloc(sizeof(double*[4]));
  double **spacetime_ricci_tensor = gkyl_malloc(sizeof(double*[4]));

  for (int i = 0; i < 4; i++) {
    spacetime_metric_tensor[i] = gkyl_malloc(sizeof(double[4]));
    spacetime_ricci_tensor[i] = gkyl_malloc(sizeof(double[4]));
  }

  double spacetime_ricci_scalar;

  spacetime->spacetime_riemann_tensor_func(spacetime, t, x, y, z, dt, dx, dy, dz, &spacetime_riemann_tensor);
  spacetime->spacetime_metric_tensor_func(spacetime, t, x, y, z, &spacetime_metric_tensor);
  spacetime->spacetime_ricci_tensor_func(spacetime, t, x, y, z, dt, dx, dy, dz, &spacetime_ricci_tensor);
  spacetime->spacetime_ricci_scalar_func(spacetime, t, x, y, z, dt, dx, dy, dz, &spacetime_ricci_scalar);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        for (int l = 0; l < 4; l++) {
          for (int m = 0; m < 4; m++) {
            covariant_spacetime_riemann_tensor[i][j][k][l] += spacetime_metric_tensor[i][m] * spacetime_riemann_tensor[m][j][k][l];
          }
        }
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        for (int l = 0; l < 4; l++) {
          (*spacetime_weyl_tensor)[i][j][k][l] = covariant_spacetime_riemann_tensor[i][j][k][l] + (0.5 * ((spacetime_ricci_tensor[i][l] * spacetime_metric_tensor[j][k])
            - (spacetime_ricci_tensor[i][k] * spacetime_metric_tensor[j][l]) + (spacetime_ricci_tensor[j][k] * spacetime_metric_tensor[i][l])
            - (spacetime_ricci_tensor[j][l] * spacetime_metric_tensor[i][k]))) + ((spacetime_ricci_scalar / 6.0) * ((spacetime_metric_tensor[i][k] * spacetime_metric_tensor[j][l])
            - (spacetime_metric_tensor[i][l] * spacetime_metric_tensor[j][k])));
        }
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        gkyl_free(spacetime_riemann_tensor[i][j][k]);
        gkyl_free(covariant_spacetime_riemann_tensor[i][j][k]);
      }
      gkyl_free(spacetime_riemann_tensor[i][j]);
      gkyl_free(covariant_spacetime_riemann_tensor[i][j]);
    }
    gkyl_free(spacetime_riemann_tensor[i]);
    gkyl_free(covariant_spacetime_riemann_tensor[i]);
    gkyl_free(spacetime_metric_tensor[i]);
    gkyl_free(spacetime_ricci_tensor[i]);
  }
  gkyl_free(spacetime_riemann_tensor);
  gkyl_free(covariant_spacetime_riemann_tensor);
  gkyl_free(spacetime_metric_tensor);
  gkyl_free(spacetime_ricci_tensor);
}

void
gkyl_gr_conformal_factor_diff(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** conformal_factor_diff)
{
  double conformal_factor_x_forward;
  double conformal_factor_y_forward;
  double conformal_factor_z_forward;

  double conformal_factor_x_backward;
  double conformal_factor_y_backward;
  double conformal_factor_z_backward;

  spacetime->conformal_factor_func(spacetime, t, x + (0.5 * dx), y, z, &conformal_factor_x_forward);
  spacetime->conformal_factor_func(spacetime, t, x, y + (0.5 * dy), z, &conformal_factor_y_forward);
  spacetime->conformal_factor_func(spacetime, t, x, y, z + (0.5 * dz), &conformal_factor_z_forward);

  spacetime->conformal_factor_func(spacetime, t, x - (0.5 * dx), y, z, &conformal_factor_x_backward);
  spacetime->conformal_factor_func(spacetime, t, x, y - (0.5 * dy), z, &conformal_factor_y_backward);
  spacetime->conformal_factor_func(spacetime, t, x, y, z - (0.5 * dz), &conformal_factor_z_backward);

  (*conformal_factor_diff)[0] = (1.0 / dx) * (conformal_factor_x_forward - conformal_factor_x_backward);
  (*conformal_factor_diff)[1] = (1.0 / dy) * (conformal_factor_y_forward - conformal_factor_y_backward);
  (*conformal_factor_diff)[2] = (1.0 / dz) * (conformal_factor_z_forward - conformal_factor_z_backward);
}

void
gkyl_gr_bssn_conformal_factor_diff(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** bssn_conformal_factor_diff)
{
  double bssn_conformal_factor_x_forward;
  double bssn_conformal_factor_y_forward;
  double bssn_conformal_factor_z_forward;

  double bssn_conformal_factor_x_backward;
  double bssn_conformal_factor_y_backward;
  double bssn_conformal_factor_z_backward;

  spacetime->bssn_conformal_factor_func(spacetime, t, x + (0.5 * dx), y, z, &bssn_conformal_factor_x_forward);
  spacetime->bssn_conformal_factor_func(spacetime, t, x, y + (0.5 * dy), z, &bssn_conformal_factor_y_forward);
  spacetime->bssn_conformal_factor_func(spacetime, t, x, y, z + (0.5 * dz), &bssn_conformal_factor_z_forward);

  spacetime->bssn_conformal_factor_func(spacetime, t, x - (0.5 * dx), y, z, &bssn_conformal_factor_x_backward);
  spacetime->bssn_conformal_factor_func(spacetime, t, x, y - (0.5 * dy), z, &bssn_conformal_factor_y_backward);
  spacetime->bssn_conformal_factor_func(spacetime, t, x, y, z - (0.5 * dz), &bssn_conformal_factor_z_backward);

  (*bssn_conformal_factor_diff)[0] = (1.0 / dx) * (bssn_conformal_factor_x_forward - bssn_conformal_factor_x_backward);
  (*bssn_conformal_factor_diff)[1] = (1.0 / dy) * (bssn_conformal_factor_y_forward - bssn_conformal_factor_y_backward);
  (*bssn_conformal_factor_diff)[2] = (1.0 / dz) * (bssn_conformal_factor_z_forward - bssn_conformal_factor_z_backward);
}

void
gkyl_gr_conformal_factor_diff2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** conformal_factor_diff2)
{
  double* conformal_factor_der_x_forward = gkyl_malloc(sizeof(double[3]));
  double* conformal_factor_der_y_forward = gkyl_malloc(sizeof(double[3]));
  double* conformal_factor_der_z_forward = gkyl_malloc(sizeof(double[3]));

  double* conformal_factor_der_x_backward = gkyl_malloc(sizeof(double[3]));
  double* conformal_factor_der_y_backward = gkyl_malloc(sizeof(double[3]));
  double* conformal_factor_der_z_backward = gkyl_malloc(sizeof(double[3]));

  spacetime->conformal_factor_der_func(spacetime, t, x + (0.5 * dx), y, z, dx, dy, dz, &conformal_factor_der_x_forward);
  spacetime->conformal_factor_der_func(spacetime, t, x, y + (0.5 * dy), z, dx, dy, dz, &conformal_factor_der_y_forward);
  spacetime->conformal_factor_der_func(spacetime, t, x, y, z + (0.5 * dz), dx, dy, dz, &conformal_factor_der_z_forward);

  spacetime->conformal_factor_der_func(spacetime, t, x - (0.5 * dx), y, z, dx, dy, dz, &conformal_factor_der_x_backward);
  spacetime->conformal_factor_der_func(spacetime, t, x, y - (0.5 * dy), z, dx, dy, dz, &conformal_factor_der_y_backward);
  spacetime->conformal_factor_der_func(spacetime, t, x, y, z - (0.5 * dz), dx, dy, dz, &conformal_factor_der_z_backward);

  for (int i = 0; i < 3; i++) {
    (*conformal_factor_diff2)[0][i] = (1.0 / dx) * (conformal_factor_der_x_forward[i] - conformal_factor_der_x_backward[i]);
    (*conformal_factor_diff2)[1][i] = (1.0 / dy) * (conformal_factor_der_y_forward[i] - conformal_factor_der_y_backward[i]);
    (*conformal_factor_diff2)[2][i] = (1.0 / dz) * (conformal_factor_der_z_forward[i] - conformal_factor_der_z_backward[i]);
  }

  gkyl_free(conformal_factor_der_x_forward);
  gkyl_free(conformal_factor_der_y_forward);
  gkyl_free(conformal_factor_der_z_forward);

  gkyl_free(conformal_factor_der_x_backward);
  gkyl_free(conformal_factor_der_y_backward);
  gkyl_free(conformal_factor_der_z_backward);
}

void
gkyl_gr_bssn_conformal_factor_diff2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** bssn_conformal_factor_diff2)
{
  double* bssn_conformal_factor_der_x_forward = gkyl_malloc(sizeof(double[3]));
  double* bssn_conformal_factor_der_y_forward = gkyl_malloc(sizeof(double[3]));
  double* bssn_conformal_factor_der_z_forward = gkyl_malloc(sizeof(double[3]));

  double* bssn_conformal_factor_der_x_backward = gkyl_malloc(sizeof(double[3]));
  double* bssn_conformal_factor_der_y_backward = gkyl_malloc(sizeof(double[3]));
  double* bssn_conformal_factor_der_z_backward = gkyl_malloc(sizeof(double[3]));

  spacetime->bssn_conformal_factor_der_func(spacetime, t, x + (0.5 * dx), y, z, dx, dy, dz, &bssn_conformal_factor_der_x_forward);
  spacetime->bssn_conformal_factor_der_func(spacetime, t, x, y + (0.5 * dy), z, dx, dy, dz, &bssn_conformal_factor_der_y_forward);
  spacetime->bssn_conformal_factor_der_func(spacetime, t, x, y, z + (0.5 * dz), dx, dy, dz, &bssn_conformal_factor_der_z_forward);

  spacetime->bssn_conformal_factor_der_func(spacetime, t, x - (0.5 * dx), y, z, dx, dy, dz, &bssn_conformal_factor_der_x_backward);
  spacetime->bssn_conformal_factor_der_func(spacetime, t, x, y - (0.5 * dy), z, dx, dy, dz, &bssn_conformal_factor_der_y_backward);
  spacetime->bssn_conformal_factor_der_func(spacetime, t, x, y, z - (0.5 * dz), dx, dy, dz, &bssn_conformal_factor_der_z_backward);

  for (int i = 0; i < 3; i++) {
    (*bssn_conformal_factor_diff2)[0][i] = (1.0 / dx) * (bssn_conformal_factor_der_x_forward[i] - bssn_conformal_factor_der_x_backward[i]);
    (*bssn_conformal_factor_diff2)[1][i] = (1.0 / dy) * (bssn_conformal_factor_der_y_forward[i] - bssn_conformal_factor_der_y_backward[i]);
    (*bssn_conformal_factor_diff2)[2][i] = (1.0 / dz) * (bssn_conformal_factor_der_z_forward[i] - bssn_conformal_factor_der_z_backward[i]);
  }

  gkyl_free(bssn_conformal_factor_der_x_forward);
  gkyl_free(bssn_conformal_factor_der_y_forward);
  gkyl_free(bssn_conformal_factor_der_z_forward);

  gkyl_free(bssn_conformal_factor_der_x_backward);
  gkyl_free(bssn_conformal_factor_der_y_backward);
  gkyl_free(bssn_conformal_factor_der_z_backward);
}