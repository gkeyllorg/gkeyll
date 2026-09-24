#pragma once

#include <math.h>
#include <gkyl_alloc.h>
#include <gkyl_ref_count.h>

typedef void (*emission_elastic_func_t)(double t, const double *xn, double *fout, void *ctx);

// Base model type
struct gkyl_emission_elastic_model {
  int cdim;
  int vdim;
  double mass;
  double charge;
  emission_elastic_func_t function;

  struct gkyl_emission_elastic_model *on_dev;
  struct gkyl_ref_count ref_count; // reference count
};

// Furman-Pivi model container
struct gkyl_emission_elastic_furman_pivi {
  struct gkyl_emission_elastic_model elastic;
  double P1_inf;
  double P1_hat;
  double E_hat;
  double W;
  double p;
};

// Cazaux model container
struct gkyl_emission_elastic_cazaux {
  struct gkyl_emission_elastic_model elastic;
  double E_f;
  double phi;
};

// Constant emission model container
struct gkyl_emission_elastic_constant {
  struct gkyl_emission_elastic_model elastic;
  double delta;
};

// Free functions

static void
gkyl_emission_elastic_furman_pivi_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_emission_elastic_model *elastic =
    container_of(ref, struct gkyl_emission_elastic_model, ref_count);
  struct gkyl_emission_elastic_furman_pivi *model = container_of(elastic,
    struct gkyl_emission_elastic_furman_pivi, elastic);
  gkyl_free(model);
}

static void
gkyl_emission_elastic_cazaux_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_emission_elastic_model *elastic =
    container_of(ref, struct gkyl_emission_elastic_model, ref_count);
  struct gkyl_emission_elastic_cazaux *model = container_of(elastic,
    struct gkyl_emission_elastic_cazaux, elastic);
  gkyl_free(model);
}

static void
gkyl_emission_elastic_constant_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_emission_elastic_model *elastic =
    container_of(ref, struct gkyl_emission_elastic_model, ref_count);
  struct gkyl_emission_elastic_constant *model = container_of(elastic,
    struct gkyl_emission_elastic_constant, elastic);
  gkyl_free(model);
}

// Furman-Pivi SEY calculation
GKYL_CU_D
static void
gkyl_emission_elastic_furman_pivi_yield(double t, const double *xn, double *fout, void *ctx)
// Electron impact model adapted from https://link.aps.org/doi/10.1103/PhysRevSTAB.5.124404
{
  struct gkyl_emission_elastic_model *elastic = (struct gkyl_emission_elastic_model *) ctx;
  const struct gkyl_emission_elastic_furman_pivi *model = container_of(elastic,
    struct gkyl_emission_elastic_furman_pivi, elastic);

  int cdim = elastic->cdim;
  int vdim = elastic->vdim;
  double mass = elastic->mass;
  double charge = elastic->charge;

  double P1_inf = model->P1_inf;
  double P1_hat = model->P1_hat;
  double E_hat = model->E_hat;
  double W = model->W;
  double p = model->p;

  double E = 0.0;
  double mu = 1.0; // currently hardcoded to normal, will add angular dependence later
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xn[cdim+d]*xn[cdim+d]/fabs(charge); // Calculate energy in eV
  }
  fout[0] = P1_inf + (P1_hat - P1_inf)*exp(pow(-fabs(E - E_hat)/W, p)/p);
}

// Cazaux backscattering */
GKYL_CU_D
static void
gkyl_emission_elastic_cazaux_yield(double t, const double *xn, double *fout, void *ctx)
// Low-energy backscattering model adapted from https://doi.org/10.1063/1.3691956
{
  struct gkyl_emission_elastic_model *elastic = (struct gkyl_emission_elastic_model *) ctx;
  const struct gkyl_emission_elastic_cazaux *model = container_of(elastic,
    struct gkyl_emission_elastic_cazaux, elastic);
  int cdim = elastic->cdim;
  int vdim = elastic->vdim;
  double mass = elastic->mass;
  double charge = elastic->charge;
  double E_f = model->E_f;
  double phi = model->phi;

  double E = 0.0;
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xn[cdim+d]*xn[cdim+d]/fabs(charge);  // Calculate energy in eV
  }
  double E_s = E + E_f + phi;
  double G = 1 + (E_s - E)/E;

  fout[0] = pow(1 - sqrt(G), 2)/pow(1 + sqrt(G), 2);
}

// Fixed constant reflection */
GKYL_CU_D
static void
gkyl_emission_elastic_constant_yield(double t, const double *xn, double *fout, void *ctx)
{
  struct gkyl_emission_elastic_model *elastic = (struct gkyl_emission_elastic_model *) ctx;
  const struct gkyl_emission_elastic_constant *model = container_of(elastic,
    struct gkyl_emission_elastic_constant, elastic);
  double delta = model->delta;

  fout[0] = delta;
}

/**
 * Create the elastic emission model using Furman-Pivi
 *
 * @param charge Elementary charge, used for eV units
 * @param P1_inf Fitting parameter
 * @param P1_hat Fitting parameter
 * @param E_hat Fitting parameter
 * @param W Fitting parameter
 * @param p Fitting parameter
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_elastic_model*
gkyl_emission_elastic_furman_pivi_new(double charge, double P1_inf, double P1_hat, double E_hat,
  double W, double p, bool use_gpu);

/**
 * Create the elastic emission model using Cazaux
 *
 * @param charge Elementary charge, used for eV units
 * @param E_f Fitting parameter
 * @param phi Fitting parameter
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_elastic_model*
gkyl_emission_elastic_cazaux_new(double charge, double E_f, double phi, bool use_gpu);

/**
 * Create the elastic emission model using constant yield
 *
 * @param charge Elementary charge, used for eV units
 * @param delta Yield value
 * @param use_gpu bool to determine if on GPU
 * @return New model
 */
struct gkyl_emission_elastic_model*
gkyl_emission_elastic_constant_new(double charge, double delta, bool use_gpu);

/**
 * Acquire pointer to model object. Delete using the release()
 * method
 *
 * @param model Model object.
 * @return Acquired model obj pointer
 */
struct gkyl_emission_elastic_model*
gkyl_emission_elastic_model_acquire(const struct gkyl_emission_elastic_model* model);

/**
 * Delete model object
 *
 * @param model Model object to delete.
 */
void
gkyl_emission_elastic_model_release(const struct gkyl_emission_elastic_model* model);
