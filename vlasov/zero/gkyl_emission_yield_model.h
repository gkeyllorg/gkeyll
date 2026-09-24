#pragma once

#include <math.h>
#include <assert.h>
#include <gkyl_alloc.h>
#include <gkyl_ref_count.h>

// Object type
struct gkyl_emission_yield_model;

typedef void (*emission_yield_func_t)(double *out, struct gkyl_emission_yield_model *yield,
  double xc[GKYL_MAX_DIM]);

// Base model type
struct gkyl_emission_yield_model {
  int cdim;
  int vdim;
  double mass;
  double charge;
  emission_yield_func_t function;

  uint32_t flags;
  struct gkyl_emission_yield_model *on_dev;
  struct gkyl_ref_count ref_count; // reference count
};

// Furman-Pivi model container
struct gkyl_emission_yield_furman_pivi {
  struct gkyl_emission_yield_model yield;
  double deltahat_ts;
  double Ehat_ts;
  double t1;
  double t2;
  double t3;
  double t4;
  double s;
};

// Schou model container
struct gkyl_emission_yield_schou {
  struct gkyl_emission_yield_model yield;
  double int_wall;
  double a2;
  double a3;
  double a4;
  double a5;
  double nw;
};

// Schou model container (SRIM stopping power)
struct gkyl_emission_yield_schou_srim {
  struct gkyl_emission_yield_model yield;
  double int_wall;
  double lorentz_norm;
  double E0;
  double tau;
  double alpha;
  double beta;
  double gauss_norm;
  double gauss_E0;
  double gauss_tau;
};

// Constant yield model container
struct gkyl_emission_yield_constant {
  struct gkyl_emission_yield_model yield;
  double delta;
};

// Free functions

/**
 * Check if model is on device.
 *
 * @param model Model to check
 * @return true if model on device, false otherwise
 */
bool
gkyl_emission_yield_model_is_cu_dev(const struct gkyl_emission_yield_model *model);

static void
gkyl_emission_yield_furman_pivi_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_emission_yield_model *yield =
    container_of(ref, struct gkyl_emission_yield_model, ref_count);

  if (gkyl_emission_yield_model_is_cu_dev(yield)) {
    struct gkyl_emission_yield_furman_pivi *model = container_of(yield->on_dev,
      struct gkyl_emission_yield_furman_pivi, yield);
    gkyl_cu_free(model);
  }

  struct gkyl_emission_yield_furman_pivi *model = container_of(yield,
    struct gkyl_emission_yield_furman_pivi, yield);
  gkyl_free(model);
}

static void
gkyl_emission_yield_schou_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_emission_yield_model *yield =
    container_of(ref, struct gkyl_emission_yield_model, ref_count);

  if (gkyl_emission_yield_model_is_cu_dev(yield)) {
    struct gkyl_emission_yield_schou *model = container_of(yield->on_dev,
      struct gkyl_emission_yield_schou, yield);
    gkyl_cu_free(model);
  }

  struct gkyl_emission_yield_schou *model = container_of(yield,
    struct gkyl_emission_yield_schou, yield);
  gkyl_free(model);
}

// SRIM
static void
gkyl_emission_yield_schou_srim_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_emission_yield_model *yield =
    container_of(ref, struct gkyl_emission_yield_model, ref_count);

  if (gkyl_emission_yield_model_is_cu_dev(yield)) {
    struct gkyl_emission_yield_schou_srim *model = container_of(yield->on_dev,
      struct gkyl_emission_yield_schou_srim, yield);
    gkyl_cu_free(model);
  }

  struct gkyl_emission_yield_schou_srim *model = container_of(yield,
    struct gkyl_emission_yield_schou_srim, yield);
  gkyl_free(model);
}

static void
gkyl_emission_yield_constant_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_emission_yield_model *yield =
    container_of(ref, struct gkyl_emission_yield_model, ref_count);

  if (gkyl_emission_yield_model_is_cu_dev(yield)) {
    struct gkyl_emission_yield_constant *model = container_of(yield->on_dev,
      struct gkyl_emission_yield_constant, yield);
    gkyl_cu_free(model);
  }

  struct gkyl_emission_yield_constant *model = container_of(yield,
    struct gkyl_emission_yield_constant, yield);
  gkyl_free(model);
}

// Furman-Pivi SEY calculation */
GKYL_CU_D
static void
gkyl_emission_yield_furman_pivi_yield(double *out, struct gkyl_emission_yield_model *yield,
  double xc[GKYL_MAX_DIM])
// Electron impact model adapted from https://link.aps.org/doi/10.1103/PhysRevSTAB.5.124404
{
  const struct gkyl_emission_yield_furman_pivi *model = container_of(yield,
    struct gkyl_emission_yield_furman_pivi, yield);
  
  int cdim = yield->cdim;
  int vdim = yield->vdim;
  double mass = yield->mass;
  double charge = yield->charge;

  double deltahat_ts = model->deltahat_ts;
  double Ehat_ts = model->Ehat_ts;
  double t1 = model->t1;
  double t2 = model->t2;
  double t3 = model->t3;
  double t4 = model->t4;
  double s = model->s;
  
  double E = 0.0;
  double mu = 1.0; // currently hardcoded to normal, will add angular dependence later
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xc[cdim+d]*xc[cdim+d]/fabs(charge); // Calculate energy in eV
  }
  double deltahat = deltahat_ts*(1 + t1*(1 - pow(mu, t2)));
  double Ehat = Ehat_ts*(1 + t3*(1 - pow(mu, t4)));
  double x = E/Ehat;
  out[0] = deltahat*s*x/(s - 1 + pow(x, s));
}

// Schou SEY calculation
GKYL_CU_D
static void
gkyl_emission_yield_schou_yield(double *out, struct gkyl_emission_yield_model *yield,
  double xc[GKYL_MAX_DIM])
// Ion impact model adapted from https://doi.org/10.1103/PhysRevB.22.2141
{ // No angular dependence atm. Will have to add later
  const struct gkyl_emission_yield_schou *model = container_of(yield,
    struct gkyl_emission_yield_schou, yield);
  int cdim = yield->cdim;
  int vdim = yield->vdim;
  double mass = yield->mass;
  double charge = yield->charge;
  double int_wall = model->int_wall;
  double A2 = model->a2;  // Note: Starts at 2 to match notation from source: https://doi.org/10.1093/jicru_os25.2.18
  double A3 = model->a3;
  double A4 = model->a4;
  double A5 = model->a5;
  double nw = model->nw;   // Number density of wall material in m^-3

  double E = 0.0;
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xc[cdim+d]*xc[cdim+d]/fabs(charge)/1000;  // Calculate energy in keV
  }

  double Es = E*1.66053906660e-27/1.67262192369e-27;  // Scale energy by ratio of atomic mass unit to proton mass
  double eps_low = A2*pow(Es,0.45);
  double eps_high = 0.0;
  if (Es != 0) { // Divide by zero error catching. If Es is not 0, then do the normal calculation.
    eps_high = (A3/Es)*log(1+A4/Es+A5*Es);
  }
  double eps = eps_low*eps_high/(eps_low+eps_high);

  out[0] =  eps*nw*fabs(charge)*int_wall/1.0e19;
}

// Schou SEY calculation w/ SRIM
GKYL_CU_D
static void
gkyl_emission_yield_schou_srim_yield(double *out, struct gkyl_emission_yield_model *yield,
  double xc[GKYL_MAX_DIM])
// Ion impact model adapted from https://doi.org/10.1103/PhysRevB.22.2141
{ // No angular dependence atm. Will have to add later
  const struct gkyl_emission_yield_schou_srim *model = container_of(yield,
    struct gkyl_emission_yield_schou_srim, yield);
  int cdim = yield->cdim;
  int vdim = yield->vdim;
  double mass = yield->mass;
  double charge = yield->charge;
  double int_wall = model->int_wall;
  double lorentz_norm = model->lorentz_norm;
  double E0 = model->E0;
  double tau = model->tau;
  double alpha = model->alpha;
  double beta = model->beta;
  double gauss_norm = model->gauss_norm;
  double gauss_E0 = model->gauss_E0;
  double gauss_tau = model->gauss_tau;

  double E = 0.0;
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xc[cdim+d]*xc[cdim+d]/fabs(charge)/1000;  // Calculate energy in keV
  }

  double Si_e = 1 / (1 + (pow(log(E/E0), 2) / (2.0 * pow(tau, 2)))); // Electronic stopping power
  if (E <= E0) {
    Si_e = pow(Si_e, alpha);
  }
  else {
    Si_e = pow(Si_e, beta);
  }

  double Si_n = exp(-pow(log(E/gauss_E0), 2)/(2.0*pow(gauss_tau, 2))); // Nuclear stopping power

  out[0] = (Si_e*lorentz_norm + Si_n*gauss_norm)*int_wall;
}

// Fixed constant SEY
GKYL_CU_D
static void
gkyl_emission_yield_constant_yield(double *out, struct gkyl_emission_yield_model *yield,
  double xc[GKYL_MAX_DIM])
{
  const struct gkyl_emission_yield_constant *model = container_of(yield,
    struct gkyl_emission_yield_constant, yield);
  double delta = model->delta;

  out[0] = delta;
}

/**
 * Create the emission yield model using Furman-Pivi
 *
 * @param charge Elementary charge, used for eV units
 * @param deltahat_ts Fitting parameter
 * @param Ehat_ts Fitting parameter
 * @param t1 Fitting parameter
 * @param t2 Fitting parameter
 * @param t3 Fitting parameter
 * @param t4 Fitting parameter
 * @param s Fitting parameter
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_yield_model*
gkyl_emission_yield_furman_pivi_new(double charge, double deltahat_ts, double Ehat_ts, double t1,
  double t2, double t3, double t4, double s, bool use_gpu);

/**
 * Create the emission yield model using Schou
 *
 * @param charge Elementary charge, used for eV units
 * @param int_wall Fitting parameter
 * @param a2 Fitting parameter
 * @param a3 Fitting parameter
 * @param a4 Fitting parameter
 * @param a5 Fitting parameter
 * @param nw Fitting parameter
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_yield_model*
gkyl_emission_yield_schou_new(double charge, double int_wall, double a2, double a3, double a4,
  double a5, double nw, bool use_gpu);

/**
 * Create the emission yield model using Schou (SRIM ion stopping power)
 *
 * @param charge Elementary charge, used for eV units
 * @param int_wall Fitting parameter
 * @param lorentz_norm Normalization of electronic stopping power
 * @param E0 Fitting parameter (electronic stopping)
 * @param tau Fitting parameter (electronic stopping)
 * @param alpha Fitting parameter (electronic stopping)
 * @param beta Fitting parameter (electronic stopping)
 * @param gauss_norm Normaliation of nuclear stopping power
 * @param gauss_E0 Fitting parameter (nuclear stopping)
 * @param gauss_tau Fitting parameter (nuclear stopping)
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_yield_model*
gkyl_emission_yield_schou_srim_new(double charge, double int_wall, double lorentz_norm, double E0, double tau, double alpha, 
  double beta, double gauss_norm, double gauss_E0, double gauss_tau, bool use_gpu);

/**
 * Create the emission yield model using a constant yield
 *
 * @param charge Elementary charge, used for eV units
 * @param delta Yield value
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_yield_model*
gkyl_emission_yield_constant_new(double charge, double delta, bool use_gpu);

/**
 * Acquire pointer to model object. Delete using the release()
 * method
 *
 * @param model Model object.
 * @return Acquired model obj pointer
 */
struct gkyl_emission_yield_model*
gkyl_emission_yield_model_acquire(const struct gkyl_emission_yield_model* model);

/**
 * Delete model object
 *
 * @param model Model object to delete.
 */
void
gkyl_emission_yield_model_release(const struct gkyl_emission_yield_model* model);

/**
 * Create the emission yield model using Furman-Pivi on NV-GPU
 *
 * @param charge Elementary charge, used for eV units
 * @param deltahat_ts Fitting parameter
 * @param Ehat_ts Fitting parameter
 * @param t1 Fitting parameter
 * @param t2 Fitting parameter
 * @param t3 Fitting parameter
 * @param t4 Fitting parameter
 * @param s Fitting parameter
 * @return New model
 */
struct gkyl_emission_yield_model*
gkyl_emission_yield_furman_pivi_cu_dev_new(double charge, double deltahat_ts, double Ehat_ts,
  double t1, double t2, double t3, double t4, double s);


/**
 * Create the emission yield model using Schou on NV-GPU
 *
 * @param charge Elementary charge, used for eV units
 * @param int_wall Fitting parameter
 * @param a2 Fitting parameter
 * @param a3 Fitting parameter
 * @param a4 Fitting parameter
 * @param a5 Fitting parameter
 * @param nw Fitting parameter
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_yield_model*
gkyl_emission_yield_schou_cu_dev_new(double charge, double int_wall, double a2, double a3,
  double a4, double a5, double nw);

/**
 * Create the emission yield model using Schou (SRIM stopping power) on NV-GPU
 *
 * @param charge Elementary charge, used for eV units
 * @param int_wall Fitting parameter
 * @param lorentz_norm Normalization of electronic stopping power
 * @param E0 Fitting parameter (electronic stopping)
 * @param tau Fitting parameter (electronic stopping)
 * @param alpha Fitting parameter (electronic stopping)
 * @param beta Fitting parameter (electronic stopping)
 * @param gauss_norm Normalization of nuclear stopping power
 * @param gauss_E0 Fitting parameter (nuclear stopping)
 * @param gauss_tau Fitting parameter (nuclear stopping)
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_yield_model*
gkyl_emission_yield_schou_srim_cu_dev_new(double charge, double int_wall, double lorentz_norm, double E0, double tau,
  double alpha, double beta, double gauss_norm, double gauss_E0, double gauss_tau);

/**
 * Create the emission yield model using a constant yield on NV-GPU
 *
 * @param charge Elementary charge, used for eV units
 * @param delta Yield value
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_yield_model*
gkyl_emission_yield_constant_cu_dev_new(double charge, double delta);
