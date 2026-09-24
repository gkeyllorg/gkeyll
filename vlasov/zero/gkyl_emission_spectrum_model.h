#pragma once

#include <math.h>
#include <gkyl_alloc.h>
#include <gkyl_ref_count.h>

// Object type
struct gkyl_emission_spectrum_model;

typedef void (*emission_spectrum_dist_func_t)(double t, const double *xn, double *fout, void *ctx);
typedef void (*emission_spectrum_norm_func_t)(double *out, struct gkyl_emission_spectrum_model *spectrum,
  const double *flux, double effective_delta);

// Base model type
struct gkyl_emission_spectrum_model {
  int cdim;
  int vdim;
  double mass;
  double charge;
  emission_spectrum_dist_func_t distribution;
  emission_spectrum_norm_func_t normalization;
  
  uint32_t flags;
  struct gkyl_emission_spectrum_model *on_dev;
  struct gkyl_ref_count ref_count; // reference count
};

// Chung-Everhart model container
struct gkyl_emission_spectrum_chung_everhart {
  struct gkyl_emission_spectrum_model spectrum;
  double phi;
};

// Logarithmic Gaussian model container
struct gkyl_emission_spectrum_gaussian {
  struct gkyl_emission_spectrum_model spectrum;
  double E_0;
  double tau;
};

// Maxwellian model container
struct gkyl_emission_spectrum_maxwellian {
  struct gkyl_emission_spectrum_model spectrum;
  double vt;
};

// Free functions

/**
 * Check if model is on device.
 *
 * @param model Model to check
 * @return true if model on device, false otherwise
 */
bool
gkyl_emission_spectrum_model_is_cu_dev(const struct gkyl_emission_spectrum_model *model);

static void
gkyl_emission_spectrum_chung_everhart_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_emission_spectrum_model *spectrum =
    container_of(ref, struct gkyl_emission_spectrum_model, ref_count);

  if (gkyl_emission_spectrum_model_is_cu_dev(spectrum)) {
    struct gkyl_emission_spectrum_chung_everhart *model = container_of(spectrum->on_dev,
      struct gkyl_emission_spectrum_chung_everhart, spectrum);
    gkyl_cu_free(model);
  }

  struct gkyl_emission_spectrum_chung_everhart *model = container_of(spectrum,
    struct gkyl_emission_spectrum_chung_everhart, spectrum);
  gkyl_free(model);
}

static void
gkyl_emission_spectrum_gaussian_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_emission_spectrum_model *spectrum =
    container_of(ref, struct gkyl_emission_spectrum_model, ref_count);

  if (gkyl_emission_spectrum_model_is_cu_dev(spectrum)) {
    struct gkyl_emission_spectrum_gaussian *model = container_of(spectrum->on_dev,
      struct gkyl_emission_spectrum_gaussian, spectrum);
    gkyl_cu_free(model);
  }

  struct gkyl_emission_spectrum_gaussian *model = container_of(spectrum,
    struct gkyl_emission_spectrum_gaussian, spectrum);
  gkyl_free(model);
}

static void
gkyl_emission_spectrum_maxwellian_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_emission_spectrum_model *spectrum =
    container_of(ref, struct gkyl_emission_spectrum_model, ref_count);

  if (gkyl_emission_spectrum_model_is_cu_dev(spectrum)) {
    struct gkyl_emission_spectrum_maxwellian *model = container_of(spectrum->on_dev,
      struct gkyl_emission_spectrum_maxwellian, spectrum);
    gkyl_cu_free(model);
  }

  struct gkyl_emission_spectrum_maxwellian *model = container_of(spectrum,
    struct gkyl_emission_spectrum_maxwellian, spectrum);
  gkyl_free(model);
}

// Model distribution functions

GKYL_CU_D
static void
gkyl_emission_spectrum_chung_everhart_dist(double t, const double *xn, double *fout, void *ctx)
{
  struct gkyl_emission_spectrum_model *spectrum = (struct gkyl_emission_spectrum_model *) ctx;
  const struct gkyl_emission_spectrum_chung_everhart *model = container_of(spectrum,
    struct gkyl_emission_spectrum_chung_everhart, spectrum);
  
  int cdim = spectrum->cdim;
  int vdim = spectrum->vdim;
  double mass = spectrum->mass;
  double charge = spectrum->charge;
  double phi = model->phi;

  double E = 0.0;
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xn[cdim+d]*xn[cdim+d]/fabs(charge);
  }

  fout[0] = E/pow(E + phi, 4);
}

GKYL_CU_D
static void
gkyl_emission_spectrum_gaussian_dist(double t, const double *xn, double *fout, void *ctx)
{
  struct gkyl_emission_spectrum_model *spectrum = (struct gkyl_emission_spectrum_model *) ctx;
  const struct gkyl_emission_spectrum_gaussian *model = container_of(spectrum,
    struct gkyl_emission_spectrum_gaussian, spectrum);
  int cdim = spectrum->cdim;
  int vdim = spectrum->vdim;
  double mass = spectrum->mass;
  double charge = spectrum->charge;
  double E_0 = model->E_0;
  double tau = model->tau;

  double E = 0.0;
  double mu = 1.0; // currently hardcoded to normal, will add angular dependence later
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xn[cdim+d]*xn[cdim+d]/fabs(charge);
  }

  fout[0] = exp(-pow(log(E/E_0), 2)/(2.0*pow(tau, 2)));
}

GKYL_CU_D
static void
gkyl_emission_spectrum_maxwellian_dist(double t, const double *xn, double *fout, void *ctx)
{
  struct gkyl_emission_spectrum_model *spectrum = (struct gkyl_emission_spectrum_model *) ctx;
  const struct gkyl_emission_spectrum_maxwellian *model = container_of(spectrum,
    struct gkyl_emission_spectrum_maxwellian, spectrum);
  int cdim = spectrum->cdim;
  int vdim = spectrum->vdim;
  double mass = spectrum->mass;
  double charge = spectrum->charge;
  double vt = model->vt;

  double v_sq = 0.0;
  for (int d=0; d<vdim; d++) {
    v_sq += xn[cdim+d]*xn[cdim+d];
  }

  fout[0] = exp(-v_sq/(2.0*pow(vt, 2)));
}

// Chung-Everhart normalization factor
GKYL_CU_D
static void
gkyl_emission_spectrum_chung_everhart_norm(double *out,
  struct gkyl_emission_spectrum_model *spectrum, const double *flux,
  double effective_delta)
{
  const struct gkyl_emission_spectrum_chung_everhart *model = container_of(spectrum,
    struct gkyl_emission_spectrum_chung_everhart, spectrum);
  double mass = spectrum->mass;
  double charge = spectrum->charge;
  double phi = model->phi;
  
  out[0] = 6.0*effective_delta*flux[0]*phi*phi*mass/fabs(charge);
}

// Gaussian normalization factor
GKYL_CU_D
static void
gkyl_emission_spectrum_gaussian_norm(double *out,
  struct gkyl_emission_spectrum_model *spectrum, const double *flux,
  double effective_delta)
{
  const struct gkyl_emission_spectrum_gaussian *model = container_of(spectrum,
    struct gkyl_emission_spectrum_gaussian, spectrum);
  double mass = spectrum->mass;
  double charge = spectrum->charge;
  double E_0 = model->E_0;
  double tau = model->tau;

  out[0] = effective_delta*flux[0]*mass/(sqrt(2.0*M_PI)*E_0*tau*exp(tau*tau/2.0)*fabs(charge));
}

// Maxwellian normalization factor */
GKYL_CU_D
static void
gkyl_emission_spectrum_maxwellian_norm(double *out,
  struct gkyl_emission_spectrum_model *spectrum, const double *flux,
  double effective_delta)
{
  const struct gkyl_emission_spectrum_maxwellian *model = container_of(spectrum,
    struct gkyl_emission_spectrum_maxwellian, spectrum);
  double vt = model->vt;
  int vdim = spectrum->vdim;
  
  out[0] = effective_delta*flux[0]/(pow(2.0*M_PI, (vdim - 1)/2.0)*pow(vt, vdim + 1));
}

/**
 * Create the emission spectrum model using the Chung-Everhart distribution
 *
 * @param charge Elementary charge, used for eV units
 * @param phi Work function of the emitting material
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_spectrum_model*
gkyl_emission_spectrum_chung_everhart_new(double charge, double phi, bool use_gpu);

/**
 * Create the emission spectrum model using the logarithmic Gaussian distribution
 *
 * @param charge Elementary charge, used for eV units
 * @param E_0 Fitting parameter, energy location of distribution peak
 * @param tau Fitting parameter
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_spectrum_model*
gkyl_emission_spectrum_gaussian_new(double charge, double E_0, double tau, bool use_gpu);

/**
 * Create the emission spectrum model using the Maxwellian distribution
 *
 * @param charge Elementary charge, used for eV units
 * @param vt Thermal velocity of emitted distribution
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_spectrum_model*
gkyl_emission_spectrum_maxwellian_new(double charge, double vt, bool use_gpu);

/**
 * Acquire pointer to model object. Delete using the release()
 * method
 *
 * @param model Model object.
 * @return Acquired model obj pointer
 */
struct gkyl_emission_spectrum_model*
gkyl_emission_spectrum_model_acquire(const struct gkyl_emission_spectrum_model* model);

/**
 * Delete model object
 *
 * @param model Model object to delete.
 */
void
gkyl_emission_spectrum_model_release(const struct gkyl_emission_spectrum_model* model);

/**
 * Create the emission spectrum model using the Chung-Everhart distribution on NV-GPU
 *
 * @param charge Elementary charge, used for eV units
 * @param phi Work function of the emitting material
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_spectrum_model*
gkyl_emission_spectrum_chung_everhart_cu_dev_new(struct gkyl_emission_spectrum_chung_everhart *model,
  double charge, double phi);

/**
 * Create the emission spectrum model using the logarithmic Gaussian distribution on NV-GPU
 *
 * @param charge Elementary charge, used for eV units
 * @param E_0 Fitting parameter, energy location of distribution peak
 * @param tau Fitting parameter
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_spectrum_model*
gkyl_emission_spectrum_gaussian_cu_dev_new(struct gkyl_emission_spectrum_gaussian *model,
  double charge, double E_0, double tau);

/**
 * Create the emission spectrum model using the Maxwellian distribution on NV-GPU
 *
 * @param charge Elementary charge, used for eV units
 * @param vt Thermal velocity of emitted distribution
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_spectrum_model*
gkyl_emission_spectrum_maxwellian_cu_dev_new(struct gkyl_emission_spectrum_maxwellian *model,
  double charge, double vt);
