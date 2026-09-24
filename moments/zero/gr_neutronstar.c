#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gr_neutronstar.h>
#include <gkyl_gr_spacetime_diff.h>

double
neutronstar_A_scalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double mass_quadrupole = neutronstar->mass_quadrupole;
  double mass_hexadecapole = neutronstar->mass_hexadecapole;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double ang_mom = spin * mass;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  return (8.0 * (rho * rho) * (z_cylindrical * z_cylindrical) * ((24.0 * (ang_mom * ang_mom) * mass) + (17.0 * (mass * mass) * mass_quadrupole) + (21.0 * mass_hexadecapole))) +
    ((rho * rho * rho * rho) * ((-10.0 * (ang_mom * ang_mom) * mass) + (7.0 * (mass * mass * mass * mass * mass)) + (32.0 * mass_quadrupole * (mass * mass)) - (21.0 * mass_hexadecapole))) +
    (8.0 * (z_cylindrical * z_cylindrical * z_cylindrical * z_cylindrical) * ((20.0 * (ang_mom * ang_mom) * mass) - (7.0 * (mass * mass * mass * mass * mass)) -
      (22.0 * mass_quadrupole * (mass * mass)) - (7.0 * mass_hexadecapole)));
}

double
neutronstar_B_scalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double mass_quadrupole = neutronstar->mass_quadrupole;
  double spin_octupole = neutronstar->spin_octupole;
  double mass_hexadecapole = neutronstar->mass_hexadecapole;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double ang_mom = spin * mass;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  return ((rho * rho * rho * rho) * ((10.0 * (ang_mom * ang_mom) * (mass * mass)) + (10.0 * mass_quadrupole * (mass * mass * mass)) + (21.0 * mass_hexadecapole * mass) +
      (7.0 * (mass_quadrupole * mass_quadrupole)))) +
    (4.0 * (z_cylindrical * z_cylindrical * z_cylindrical * z_cylindrical) * ((-40.0 * (ang_mom * ang_mom) * (mass * mass)) - (14.0 * ang_mom * spin_octupole) +
      (7.0 * (mass * mass * mass * mass * mass * mass)) + (30.0 * mass_quadrupole * (mass * mass * mass)) + (14.0 * mass_hexadecapole * mass) + (7.0 * (mass_quadrupole * mass_quadrupole)))) -
    (4.0 * (rho * rho) * (z_cylindrical * z_cylindrical) * ((27.0 * (ang_mom * ang_mom) * (mass * mass)) - (21.0 * ang_mom * spin_octupole) + (7.0 * (mass * mass * mass * mass * mass * mass)) +
      (48.0 * mass_quadrupole * (mass * mass * mass)) + (42.0 * mass_hexadecapole * mass) + (7.0 * (mass_quadrupole * mass_quadrupole))));
}

double
neutronstar_H_scalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double mass_quadrupole = neutronstar->mass_quadrupole;
  double spin_octupole = neutronstar->spin_octupole;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double ang_mom = spin * mass;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  return (4.0 * (rho * rho) * (z_cylindrical * z_cylindrical) * ((ang_mom * (mass_quadrupole - (2.0 * (mass * mass * mass)))) - (3.0 * mass * spin_octupole))) +
    ((rho * rho * rho * rho) * ((ang_mom * mass_quadrupole) + (3.0 * mass * spin_octupole)));
}

double
neutronstar_G_scalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double mass_quadrupole = neutronstar->mass_quadrupole;
  double spin_octupole = neutronstar->spin_octupole;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double ang_mom = spin * mass;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  return (rho * rho) * ((ang_mom * ang_mom * ang_mom) * (- ((rho * rho * rho * rho) + (8.0 * (z_cylindrical * z_cylindrical * z_cylindrical * z_cylindrical)) -
      (12.0 * (rho * rho) * (z_cylindrical * z_cylindrical)))) +
    ((ang_mom * mass) * ((((mass * mass * mass) + (2.0 * mass_quadrupole)) * (rho * rho * rho * rho)) - (8.0 * ((3.0 * (mass * mass * mass)) + (2.0 * mass_quadrupole)) *
      (z_cylindrical * z_cylindrical * z_cylindrical * z_cylindrical)) +
    (4.0 * (((mass * mass * mass) + (10.0 * mass_quadrupole)) * ((rho * rho) * (z_cylindrical * z_cylindrical)))))) +
    ((mass * mass) * spin_octupole) * ((3.0 * (rho * rho * rho * rho)) - (40.0 * (z_cylindrical * z_cylindrical * z_cylindrical * z_cylindrical)) +
      (12.0 * (rho * rho) * (z_cylindrical * z_cylindrical))));
}

double
neutronstar_F_scalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double spin_octupole = neutronstar->spin_octupole;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double ang_mom = spin * mass;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  return ((rho * rho * rho * rho) * (spin_octupole - (ang_mom * (mass * mass)))) - ((4.0 * (rho * rho) * (z_cylindrical * z_cylindrical)) * ((ang_mom * (mass * mass)) + spin_octupole));
}

double
neutronstar_f_function(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double mass_quadrupole = neutronstar->mass_quadrupole;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double ang_mom = spin * mass;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  double A_scalar = neutronstar_A_scalar(spacetime, x, y, z);
  double B_scalar = neutronstar_B_scalar(spacetime, x, y, z);

  return 1.0 - ((2.0 * mass) / sqrt((rho * rho) + (z_cylindrical * z_cylindrical))) + ((2.0 * (mass * mass)) / ((rho * rho) + (z_cylindrical * z_cylindrical))) +
    (((mass_quadrupole - (mass * mass * mass)) * (rho * rho)) - (2.0 * ((mass * mass * mass) + mass_quadrupole) * (z_cylindrical * z_cylindrical))) /
      pow((rho * rho) + (z_cylindrical * z_cylindrical), 5.0 / 2.0) +
    ((2.0 * (z_cylindrical * z_cylindrical) * (-(ang_mom * ang_mom) + (mass * mass * mass * mass) + (2.0 * mass_quadrupole * mass))) - (2.0 * mass * mass_quadrupole * (rho * rho))) /
      pow((rho * rho) + (z_cylindrical * z_cylindrical), 3.0) +
    (A_scalar / (28.0 * pow((rho * rho) + (z_cylindrical * z_cylindrical), 9.0 / 2.0))) + (B_scalar / (14.0 * pow((rho * rho) + (z_cylindrical * z_cylindrical), 5.0)));
}

double
neutronstar_omega_function(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double ang_mom = spin * mass;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  double H_scalar = neutronstar_H_scalar(spacetime, x, y, z);
  double G_scalar = neutronstar_G_scalar(spacetime, x, y, z);
  double F_scalar = neutronstar_F_scalar(spacetime, x, y, z);

  return (-(2.0 * ang_mom * (rho * rho)) / pow((rho * rho) + (z_cylindrical * z_cylindrical), 3.0 / 2.0)) - ((2.0 * ang_mom * mass * (rho * rho)) /
      pow((rho * rho) + (z_cylindrical * z_cylindrical), 2.0)) +
    (F_scalar / pow((rho * rho) + (z_cylindrical * z_cylindrical), 7.0 / 2.0)) + (H_scalar / (2.0 * pow((rho * rho) + (z_cylindrical * z_cylindrical), 4.0))) +
    (G_scalar / (4.0 * pow((rho * rho) + (z_cylindrical * z_cylindrical), 11.0 / 2.0)));
}

double
neutronstar_gamma_function(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double mass_quadrupole = neutronstar->mass_quadrupole;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double ang_mom = spin * mass;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  return (((rho * rho) * (((ang_mom * ang_mom) * ((rho * rho) - (8.0 * (z_cylindrical * z_cylindrical)))) + (mass * ((mass * mass * mass) +
    (3.0 * mass_quadrupole)) * ((rho * rho) - (4.0 * (z_cylindrical * z_cylindrical)))))) / (4.0 * pow((rho * rho) + (z_cylindrical * z_cylindrical), 4.0))) -
    (((mass * mass) * (rho * rho)) / (2.0 * pow((rho * rho) + (z_cylindrical * z_cylindrical), 2.0)));
}

double**
neutronstar_spatial_transformation_tensor(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double theta = atan2((y - pos_y), (x - pos_x));

  double **spatial_transformation_tensor = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_transformation_tensor[i] = gkyl_malloc(sizeof(double[3]));
  }

  spatial_transformation_tensor[0][0] = cos(theta);
  spatial_transformation_tensor[0][1] = sin(theta);
  spatial_transformation_tensor[0][2] = 0.0;

  spatial_transformation_tensor[1][0] = -sin(theta) / rho;
  spatial_transformation_tensor[1][1] = cos(theta) / rho;
  spatial_transformation_tensor[1][2] = 0.0;

  spatial_transformation_tensor[2][0] = 0.0;
  spatial_transformation_tensor[2][1] = 0.0;
  spatial_transformation_tensor[2][2] = 1.0;

  return spatial_transformation_tensor;
}

double**
neutronstar_spacetime_transformation_tensor(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double theta = atan2((y - pos_y), (x - pos_x));

  double **spacetime_transformation_tensor = gkyl_malloc(sizeof(double*[4]));
  for (int i = 0; i < 4; i++) {
    spacetime_transformation_tensor[i] = gkyl_malloc(sizeof(double[4]));
  }

  spacetime_transformation_tensor[0][0] = 1.0;
  spacetime_transformation_tensor[0][1] = 0.0;
  spacetime_transformation_tensor[0][2] = 0.0;
  spacetime_transformation_tensor[0][3] = 0.0;
  
  spacetime_transformation_tensor[1][0] = 0.0;
  spacetime_transformation_tensor[1][1] = cos(theta);
  spacetime_transformation_tensor[1][2] = sin(theta);
  spacetime_transformation_tensor[1][3] = 0.0;

  spacetime_transformation_tensor[2][0] = 0.0;
  spacetime_transformation_tensor[2][1] = -sin(theta) / rho;
  spacetime_transformation_tensor[2][2] = cos(theta) / rho;
  spacetime_transformation_tensor[2][3] = 0.0;

  spacetime_transformation_tensor[3][0] = 0.0;
  spacetime_transformation_tensor[3][1] = 0.0;
  spacetime_transformation_tensor[3][2] = 0.0;
  spacetime_transformation_tensor[3][3] = 1.0;

  return spacetime_transformation_tensor;
}

double**
neutronstar_spatial_inv_transformation_tensor(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double theta = atan2((y - pos_y), (x - pos_x));

  double **spatial_inv_transformation_tensor = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_inv_transformation_tensor[i] = gkyl_malloc(sizeof(double[3]));
  }

  spatial_inv_transformation_tensor[0][0] = cos(theta);
  spatial_inv_transformation_tensor[0][1] = -rho * sin(theta);
  spatial_inv_transformation_tensor[0][2] = 0.0;

  spatial_inv_transformation_tensor[1][0] = sin(theta);
  spatial_inv_transformation_tensor[1][1] = rho * cos(theta);
  spatial_inv_transformation_tensor[1][2] = 0.0;

  spatial_inv_transformation_tensor[2][0] = 0.0;
  spatial_inv_transformation_tensor[2][1] = 0.0;
  spatial_inv_transformation_tensor[2][2] = 1.0;

  return spatial_inv_transformation_tensor;
}

double**
neutronstar_spacetime_inv_transformation_tensor(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double theta = atan2((y - pos_y), (x - pos_x));

  double **spacetime_inv_transformation_tensor = gkyl_malloc(sizeof(double*[4]));
  for (int i = 0; i < 4; i++) {
    spacetime_inv_transformation_tensor[i] = gkyl_malloc(sizeof(double[4]));
  }

  spacetime_inv_transformation_tensor[0][0] = 1.0;
  spacetime_inv_transformation_tensor[0][1] = 0.0;
  spacetime_inv_transformation_tensor[0][2] = 0.0;
  spacetime_inv_transformation_tensor[0][3] = 0.0;
  
  spacetime_inv_transformation_tensor[1][0] = 0.0;
  spacetime_inv_transformation_tensor[1][1] = cos(theta);
  spacetime_inv_transformation_tensor[1][2] = -rho * sin(theta);
  spacetime_inv_transformation_tensor[1][3] = 0.0;

  spacetime_inv_transformation_tensor[2][0] = 0.0;
  spacetime_inv_transformation_tensor[2][1] = sin(theta);
  spacetime_inv_transformation_tensor[2][2] = rho * cos(theta);
  spacetime_inv_transformation_tensor[2][3] = 0.0;

  spacetime_inv_transformation_tensor[3][0] = 0.0;
  spacetime_inv_transformation_tensor[3][1] = 0.0;
  spacetime_inv_transformation_tensor[3][2] = 0.0;
  spacetime_inv_transformation_tensor[3][3] = 1.0;

  return spacetime_inv_transformation_tensor;
}

static void
neutronstar_spatial_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_metric_tensor)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));

  double **spatial_metric_tensor_cylindrical = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_tensor_cylindrical[i] = gkyl_malloc(sizeof(double[3]));
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      spatial_metric_tensor_cylindrical[i][j] = 0.0;
      (*spatial_metric_tensor)[i][j] = 0.0;
    }
  }

  double f_function = neutronstar_f_function(spacetime, x, y, z);
  double omega_function = neutronstar_omega_function(spacetime, x, y, z);
  double gamma_function = neutronstar_gamma_function(spacetime, x, y, z);

  spatial_metric_tensor_cylindrical[0][0] = exp(2.0 * gamma_function) / f_function;
  spatial_metric_tensor_cylindrical[1][1] = ((rho * rho) / f_function) - (f_function * (omega_function * omega_function));
  spatial_metric_tensor_cylindrical[2][2] = exp(2.0 * gamma_function) / f_function;

  double **spatial_transformation_tensor = neutronstar_spatial_transformation_tensor(spacetime, x, y, z);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          (*spatial_metric_tensor)[i][j] += spatial_transformation_tensor[k][i] * spatial_transformation_tensor[l][j] * spatial_metric_tensor_cylindrical[k][l];
        }
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric_tensor_cylindrical[i]);
    gkyl_free(spatial_transformation_tensor[i]);
  }
  gkyl_free(spatial_metric_tensor_cylindrical);
  gkyl_free(spatial_transformation_tensor);
}

static void
neutronstar_spacetime_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_metric_tensor)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));

  double **spacetime_metric_tensor_cylindrical = gkyl_malloc(sizeof(double*[4]));
  for (int i = 0; i < 4; i++) {
    spacetime_metric_tensor_cylindrical[i] = gkyl_malloc(sizeof(double[4]));
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      spacetime_metric_tensor_cylindrical[i][j] = 0.0;
      (*spacetime_metric_tensor)[i][j] = 0.0;
    }
  }

  double f_function = neutronstar_f_function(spacetime, x, y, z);
  double omega_function = neutronstar_omega_function(spacetime, x, y, z);
  double gamma_function = neutronstar_gamma_function(spacetime, x, y, z);

  spacetime_metric_tensor_cylindrical[0][0] = -f_function;
  spacetime_metric_tensor_cylindrical[1][1] = exp(2.0 * gamma_function) / f_function;
  spacetime_metric_tensor_cylindrical[2][2] = ((rho * rho) / f_function) - (f_function * (omega_function * omega_function));
  spacetime_metric_tensor_cylindrical[3][3] = exp(2.0 * gamma_function) / f_function;

  spacetime_metric_tensor_cylindrical[0][2] = f_function * omega_function;
  spacetime_metric_tensor_cylindrical[2][0] = f_function * omega_function;

  double **spacetime_transformation_tensor = neutronstar_spacetime_transformation_tensor(spacetime, x, y, z);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        for (int l = 0; l < 4; l++) {
          (*spacetime_metric_tensor)[i][j] += spacetime_transformation_tensor[k][i] * spacetime_transformation_tensor[l][j] * spacetime_metric_tensor_cylindrical[k][l];
        }
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    gkyl_free(spacetime_metric_tensor_cylindrical[i]);
    gkyl_free(spacetime_transformation_tensor[i]);
  }
  gkyl_free(spacetime_metric_tensor_cylindrical);
  gkyl_free(spacetime_transformation_tensor);
}

static void
neutronstar_spatial_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_inv_metric_tensor)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));

  double **spatial_inv_metric_tensor_cylindrical = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_inv_metric_tensor_cylindrical[i] = gkyl_malloc(sizeof(double[3]));
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      spatial_inv_metric_tensor_cylindrical[i][j] = 0.0;
      (*spatial_inv_metric_tensor)[i][j] = 0.0;
    }
  }

  double f_function = neutronstar_f_function(spacetime, x, y, z);
  double omega_function = neutronstar_omega_function(spacetime, x, y, z);
  double gamma_function = neutronstar_gamma_function(spacetime, x, y, z);
  
  spatial_inv_metric_tensor_cylindrical[0][0] = exp(-2.0 * gamma_function) * f_function;
  spatial_inv_metric_tensor_cylindrical[1][1] = f_function / ((rho * rho) - ((f_function * f_function) * (omega_function * omega_function)));
  spatial_inv_metric_tensor_cylindrical[2][2] = exp(-2.0 * gamma_function) * f_function;

  double **spatial_inv_transformation_tensor = neutronstar_spatial_inv_transformation_tensor(spacetime, x, y, z);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          (*spatial_inv_metric_tensor)[i][j] += spatial_inv_transformation_tensor[i][k] * spatial_inv_transformation_tensor[j][l] * spatial_inv_metric_tensor_cylindrical[k][l];
        }
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_inv_metric_tensor_cylindrical[i]);
    gkyl_free(spatial_inv_transformation_tensor[i]);
  }
  gkyl_free(spatial_inv_metric_tensor_cylindrical);
  gkyl_free(spatial_inv_transformation_tensor);
}

static void
neutronstar_spacetime_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_inv_metric_tensor)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));

  double **spacetime_inv_metric_tensor_cylindrical = gkyl_malloc(sizeof(double*[4]));
  for (int i = 0; i < 4; i++) {
    spacetime_inv_metric_tensor_cylindrical[i] = gkyl_malloc(sizeof(double[4]));
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      spacetime_inv_metric_tensor_cylindrical[i][j] = 0.0;
      (*spacetime_inv_metric_tensor)[i][j] = 0.0;
    }
  }

  double f_function = neutronstar_f_function(spacetime, x, y, z);
  double omega_function = neutronstar_omega_function(spacetime, x, y, z);
  double gamma_function = neutronstar_gamma_function(spacetime, x, y, z);

  spacetime_inv_metric_tensor_cylindrical[0][0] = -(1.0 / f_function) + ((f_function * (omega_function * omega_function)) / (rho * rho));
  spacetime_inv_metric_tensor_cylindrical[1][1] = exp(-2.0 * gamma_function) * f_function;
  spacetime_inv_metric_tensor_cylindrical[2][2] = f_function / (rho * rho);
  spacetime_inv_metric_tensor_cylindrical[3][3] = exp(-2.0 * gamma_function) * f_function;

  spacetime_inv_metric_tensor_cylindrical[0][2] = (f_function * omega_function) / (rho * rho);
  spacetime_inv_metric_tensor_cylindrical[2][0] = (f_function * omega_function) / (rho * rho);

  double **spacetime_inv_transformation_tensor = neutronstar_spacetime_inv_transformation_tensor(spacetime, x, y, z);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        for (int l = 0; l < 4; l++) {
          (*spacetime_inv_metric_tensor)[i][j] += spacetime_inv_transformation_tensor[i][k] * spacetime_inv_transformation_tensor[j][l] * spacetime_inv_metric_tensor_cylindrical[k][l];
        }
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    gkyl_free(spacetime_inv_metric_tensor_cylindrical[i]);
    gkyl_free(spacetime_inv_transformation_tensor[i]);
  }
  gkyl_free(spacetime_inv_metric_tensor_cylindrical);
  gkyl_free(spacetime_inv_transformation_tensor);
}

static void
neutronstar_spatial_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spatial_metric_det)
{
  double** spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  neutronstar_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);

  *spatial_metric_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
  }
  gkyl_free(spatial_metric);
}

static void
neutronstar_spacetime_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spacetime_metric_det)
{
  double spatial_metric_det;
  double lapse_function;
  neutronstar_spatial_metric_det(spacetime, t, x, y, z, &spatial_metric_det);
  neutronstar_lapse_function(spacetime, t, x, y, z, &lapse_function);

  *spacetime_metric_det = - (lapse_function * lapse_function) * spatial_metric_det;
}

static void
neutronstar_spatial_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
   const double dx, const double dy, const double dz, double**** spatial_metric_tensor_der)
{
  gkyl_gr_spatial_metric_tensor_diff(spacetime, t, x, y, z, dx, dy, dz, spatial_metric_tensor_der);
}

static void
neutronstar_spacetime_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_metric_tensor_der)
{
  gkyl_gr_spacetime_metric_tensor_diff(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_metric_tensor_der);
}

static void
neutronstar_lapse_function(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* lapse_function)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));

  double f_function = neutronstar_f_function(spacetime, x, y, z);
  double omega_function = neutronstar_omega_function(spacetime, x, y, z);

  *lapse_function = (sqrt(f_function) * rho) / sqrt((rho * rho) - ((f_function * f_function) * (omega_function * omega_function)));
}

static void
neutronstar_shift_vector(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double** shift_vector)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));

  double *shift_vector_cylindrical = gkyl_malloc(sizeof(double[3]));
  for (int i = 0; i < 3; i++) {
    (*shift_vector)[i] = 0.0;
  }

  double f_function = neutronstar_f_function(spacetime, x, y, z);
  double omega_function = neutronstar_omega_function(spacetime, x, y, z);

  shift_vector_cylindrical[0] = 0.0;
  shift_vector_cylindrical[1] = - ((f_function * f_function) * omega_function) / (-(rho * rho) + ((f_function * f_function) * (omega_function * omega_function)));
  shift_vector_cylindrical[2] = 0.0;

  double **spatial_inv_transformation_tensor = neutronstar_spatial_inv_transformation_tensor(spacetime, x, y, z);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*shift_vector)[i] += spatial_inv_transformation_tensor[i][j] * shift_vector_cylindrical[j];
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_inv_transformation_tensor[i]);
  }
  gkyl_free(shift_vector_cylindrical);
  gkyl_free(spatial_inv_transformation_tensor);
}

static void
neutronstar_lapse_function_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** lapse_function_der)
{
  gkyl_gr_lapse_function_diff(spacetime, t, x, y, z, dx, dy, dz, lapse_function_der);
}

static void
neutronstar_shift_vector_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** shift_vector_der)
{
  gkyl_gr_shift_vector_diff(spacetime, t, x, y, z, dx, dy, dz, shift_vector_der);
}

static void
neutronstar_spatial_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_christoffel)
{
  gkyl_gr_spatial_christoffel_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_christoffel);
}

static void
neutronstar_spacetime_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_christoffel)
{
  gkyl_gr_spacetime_christoffel_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_christoffel);
}

static void
neutronstar_spatial_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_riemann_tensor)
{
  gkyl_gr_spatial_riemann_tensor_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_riemann_tensor);
}

static void
neutronstar_spacetime_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_riemann_tensor)
{
  gkyl_gr_spacetime_riemann_tensor_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_riemann_tensor);
}

static void
neutronstar_spatial_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** spatial_ricci_tensor)
{
  gkyl_gr_spatial_ricci_tensor_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_ricci_tensor);
}

static void
neutronstar_spacetime_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double*** spacetime_ricci_tensor)
{
  gkyl_gr_spacetime_ricci_tensor_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_ricci_tensor);
}

static void
neutronstar_spatial_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double* spatial_ricci_scalar)
{
  gkyl_gr_spatial_ricci_scalar_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_ricci_scalar);
}

static void
neutronstar_spacetime_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double* spacetime_ricci_scalar)
{
  gkyl_gr_spacetime_ricci_scalar_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_ricci_scalar);
}

static void
neutronstar_spatial_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_weyl_tensor)
{
  gkyl_gr_spatial_weyl_tensor_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_weyl_tensor);
}

static void
neutronstar_spacetime_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_weyl_tensor)
{
  gkyl_gr_spacetime_weyl_tensor_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_weyl_tensor);
}

static void
neutronstar_extrinsic_curvature_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
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

  neutronstar_lapse_function(spacetime, t, x, y, z, &lapse_function);
  neutronstar_shift_vector(spacetime, t, x, y, z, &shift_vector);
  neutronstar_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);
  neutronstar_shift_vector_der(spacetime, t, x, y, z, dx, dy, dz, &shift_vector_der);
  neutronstar_spatial_christoffel(spacetime, t, x, y, z, dx, dy, dz, &spatial_christoffel);

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
neutronstar_excision_region(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* in_excision_region)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double r = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));

  //if (r <= (mass * (1.0 + sqrt(1.0 + (spin * spin))))) {
  if (r <= mass * 2.0) {
    *in_excision_region = true;
  }
  else {
    *in_excision_region = false;
  }
}

void
gkyl_gr_neutronstar_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_gr_spacetime* base = container_of(ref, struct gkyl_gr_spacetime, ref_count);

  if (gkyl_gr_spacetime_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct gr_neutronstar *gr_neutronstar = container_of(base->on_dev, struct gr_neutronstar, spacetime);
    gkyl_cu_free(gr_neutronstar);
  }

  struct gr_neutronstar *gr_neutronstar = container_of(base, struct gr_neutronstar, spacetime);
  gkyl_free(gr_neutronstar);
}

struct gkyl_gr_spacetime*
gkyl_gr_neutronstar_new(bool use_gpu, double mass, double spin, double mass_quadrupole, double spin_octupole, double mass_hexadecapole, double pos_x, double pos_y, double pos_z)
{
  return gkyl_gr_neutronstar_inew(&(struct gkyl_gr_neutronstar_inp) {
      .use_gpu = use_gpu,
      .mass = mass,
      .spin = spin,
      .mass_quadrupole = mass_quadrupole,
      .spin_octupole = spin_octupole,
      .mass_hexadecapole = mass_hexadecapole,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
    }
  );
}

struct gkyl_gr_spacetime*
gkyl_gr_neutronstar_inew(const struct gkyl_gr_neutronstar_inp* inp)
{
  struct gr_neutronstar *gr_neutronstar = gkyl_malloc(sizeof(struct gr_neutronstar));

  gr_neutronstar->mass = inp->mass;
  gr_neutronstar->spin = inp->spin;

  gr_neutronstar->mass_quadrupole = inp->mass_quadrupole;
  gr_neutronstar->spin_octupole = inp->spin_octupole;
  gr_neutronstar->mass_hexadecapole = inp->mass_hexadecapole;

  gr_neutronstar->pos_x = inp->pos_x;
  gr_neutronstar->pos_y = inp->pos_y;
  gr_neutronstar->pos_z = inp->pos_z;

  gr_neutronstar->spacetime.spatial_metric_tensor_func = neutronstar_spatial_metric_tensor;
  gr_neutronstar->spacetime.spacetime_metric_tensor_func = neutronstar_spacetime_metric_tensor;

  gr_neutronstar->spacetime.spatial_inv_metric_tensor_func = neutronstar_spatial_inv_metric_tensor;
  gr_neutronstar->spacetime.spacetime_inv_metric_tensor_func = neutronstar_spacetime_inv_metric_tensor;

  gr_neutronstar->spacetime.spatial_metric_det_func = neutronstar_spatial_metric_det;
  gr_neutronstar->spacetime.spacetime_metric_det_func = neutronstar_spacetime_metric_det;

  gr_neutronstar->spacetime.spatial_metric_tensor_der_func = neutronstar_spatial_metric_tensor_der;
  gr_neutronstar->spacetime.spacetime_metric_tensor_der_func = neutronstar_spacetime_metric_tensor_der;

  gr_neutronstar->spacetime.lapse_function_func = neutronstar_lapse_function;
  gr_neutronstar->spacetime.shift_vector_func = neutronstar_shift_vector;

  gr_neutronstar->spacetime.lapse_function_der_func = neutronstar_lapse_function_der;
  gr_neutronstar->spacetime.shift_vector_der_func = neutronstar_shift_vector_der;

  gr_neutronstar->spacetime.spatial_christoffel_func = neutronstar_spatial_christoffel;
  gr_neutronstar->spacetime.spacetime_christoffel_func = neutronstar_spacetime_christoffel;

  gr_neutronstar->spacetime.spatial_riemann_tensor_func = neutronstar_spatial_riemann_tensor;
  gr_neutronstar->spacetime.spacetime_riemann_tensor_func = neutronstar_spacetime_riemann_tensor;

  gr_neutronstar->spacetime.spatial_ricci_tensor_func = neutronstar_spatial_ricci_tensor;
  gr_neutronstar->spacetime.spacetime_ricci_tensor_func = neutronstar_spacetime_ricci_tensor;

  gr_neutronstar->spacetime.spatial_ricci_scalar_func = neutronstar_spatial_ricci_scalar;
  gr_neutronstar->spacetime.spacetime_ricci_scalar_func = neutronstar_spacetime_ricci_scalar;

  gr_neutronstar->spacetime.spatial_weyl_tensor_func = neutronstar_spatial_weyl_tensor;
  gr_neutronstar->spacetime.spacetime_weyl_tensor_func = neutronstar_spacetime_weyl_tensor;

  gr_neutronstar->spacetime.extrinsic_curvature_tensor_func = neutronstar_extrinsic_curvature_tensor;

  gr_neutronstar->spacetime.excision_region_func = neutronstar_excision_region;

  gr_neutronstar->spacetime.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_neutronstar->spacetime.flags);
  gr_neutronstar->spacetime.ref_count = gkyl_ref_count_init(gkyl_gr_neutronstar_free);
  gr_neutronstar->spacetime.on_dev = &gr_neutronstar->spacetime; // On the CPU, the spacetime object points to itself.

  return &gr_neutronstar->spacetime;
}