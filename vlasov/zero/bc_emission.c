#include <string.h>
#include <gkyl_bc_emission.h>

// ctx with models specified by user
struct gkyl_bc_emission_ctx*
gkyl_bc_emission_new(int num_species, double t_bound, bool elastic,
  struct gkyl_emission_spectrum_model *spectrum_model[],
  struct gkyl_emission_yield_model *yield_model[],
  struct gkyl_emission_elastic_model *elastic_model, char in_species[][128])
{
  struct gkyl_bc_emission_ctx *ctx = gkyl_malloc(sizeof(struct gkyl_bc_emission_ctx));
  
  ctx->num_species = num_species;
  ctx->t_bound = t_bound;
  ctx->elastic = elastic;
  for (int i=0; i<num_species; ++i) {
    ctx->spectrum_model[i] = gkyl_emission_spectrum_model_acquire(spectrum_model[i]);
    ctx->yield_model[i] = gkyl_emission_yield_model_acquire(yield_model[i]);
    strcpy(ctx->in_species[i], in_species[i]);
  }
  if (elastic) ctx->elastic_model = gkyl_emission_elastic_model_acquire(elastic_model);

  return ctx;
}

// SEE copper preset
struct gkyl_bc_emission_ctx*
gkyl_bc_emission_secondary_electron_copper_new(int num_species, double t_bound,
  char in_species[][128], bool use_gpu)
{
  struct gkyl_bc_emission_ctx *ctx = gkyl_malloc(sizeof(struct gkyl_bc_emission_ctx));
  
  double q0 = 1.602e-19;
  double E_0 = 1.97;
  double tau = 0.88;

  double deltahat_ts = 1.885;
  double Ehat_ts = 276.8;
  double t1 = 0.66;
  double t2 = 0.8;
  double t3 = 0.7;
  double t4 = 1.0;
  double s = 1.54;

  double P1_inf = 0.02;
  double P1_hat = 0.496;
  double E_hat = 1.0e-6;
  double W = 60.86;
  double p = 1.0;

  ctx->num_species = num_species;
  ctx->t_bound = t_bound;
  ctx->elastic = true;

  for (int i=0; i<num_species; ++i) {
    ctx->spectrum_model[i] = gkyl_emission_spectrum_gaussian_new(q0, E_0, tau, use_gpu);
    ctx->yield_model[i] = gkyl_emission_yield_furman_pivi_new(q0, deltahat_ts, Ehat_ts, t1, t2, t3,
      t4, s, use_gpu);
    strcpy(ctx->in_species[i], in_species[i]);
  }
  ctx->elastic_model = gkyl_emission_elastic_furman_pivi_new(q0, P1_inf, P1_hat, E_hat,
    W, p, use_gpu);

  return ctx;
}

// SEE oxidized lithium preset
struct gkyl_bc_emission_ctx*
gkyl_bc_emission_secondary_electron_lithium_oxidized_new(int num_species, double t_bound,
  char in_species[][128], bool use_gpu)
{
  struct gkyl_bc_emission_ctx *ctx = gkyl_malloc(sizeof(struct gkyl_bc_emission_ctx));
  
  double q0 = 1.602e-19;
  double phi = 2.3;

  double deltahat_ts = 4.208;
  double Ehat_ts = 354.52;
  double t1 = 0.66;
  double t2 = 0.8;
  double t3 = 0.7;
  double t4 = 1.0;
  double s = 1.79;

  double E_f = 290.31;
  double phi_r = 144.9;

  ctx->num_species = num_species;
  ctx->t_bound = t_bound;
  ctx->elastic = true;

  for (int i=0; i<num_species; ++i) {
    ctx->spectrum_model[i] = gkyl_emission_spectrum_chung_everhart_new(q0, phi, use_gpu);
    ctx->yield_model[i] = gkyl_emission_yield_furman_pivi_new(q0, deltahat_ts, Ehat_ts, t1, t2, t3,
      t4, s, use_gpu);
    strcpy(ctx->in_species[i], in_species[i]);
  }
  ctx->elastic_model = gkyl_emission_elastic_cazaux_new(q0, E_f, phi_r, use_gpu);

  return ctx;
}

// SEE oxidized lithium preset
struct gkyl_bc_emission_ctx*
gkyl_bc_emission_secondary_electron_lithium_clean_new(int num_species, double t_bound,
  char in_species[][128], bool use_gpu)
{
  struct gkyl_bc_emission_ctx *ctx = gkyl_malloc(sizeof(struct gkyl_bc_emission_ctx));
  
  double q0 = 1.602e-19;
  double phi = 2.3;

  double deltahat_ts = 0.567;
  double Ehat_ts = 97.18;
  double t1 = 0.66;
  double t2 = 0.8;
  double t3 = 0.7;
  double t4 = 1.0;
  double s = 1.42;

  ctx->num_species = num_species;
  ctx->t_bound = t_bound;
  ctx->elastic = false;

  for (int i=0; i<num_species; ++i) {
    ctx->spectrum_model[i] = gkyl_emission_spectrum_chung_everhart_new(q0, phi, use_gpu);
    ctx->yield_model[i] = gkyl_emission_yield_furman_pivi_new(q0, deltahat_ts, Ehat_ts, t1, t2, t3,
      t4, s, use_gpu);
    strcpy(ctx->in_species[i], in_species[i]);
  }

  return ctx;
}

// Ion-impact SEE copper preset
struct gkyl_bc_emission_ctx*
gkyl_bc_emission_ion_impact_copper_new(int num_species, double t_bound,
  char in_species[][128], bool use_gpu)
{
  struct gkyl_bc_emission_ctx *ctx = gkyl_malloc(sizeof(struct gkyl_bc_emission_ctx));
  
  double q0 = 1.602e-19;
  double E_0 = 2.80670245635625;
  double tau = 1.21017574095341;

  double A2 = 4.194;  // a2-5 are empirical fits for ion stopping power
  double A3 = 4.649e3;  // Starting at a2 to match source and prevent confusion
  double A4 = 8.113e1;
  double A5 = 2.242e-2;
  double nw = 8.491231742083982e28;  // Number density of the wall material (m^-3). Calculated a priori using rho/m/u where rho is density in kg/m^3, m is mass in amu, and u is the unified atomic mass unit
  double int_wall = 6.33665090954913e7;  // Integration of wall term (J/m).  See the paper for more info on this.

  ctx->num_species = num_species;
  ctx->t_bound = t_bound;
  ctx->elastic = false;

  for (int i=0; i<num_species; ++i) {
    ctx->spectrum_model[i] = gkyl_emission_spectrum_gaussian_new(q0, E_0, tau, use_gpu);
    ctx->yield_model[i] = gkyl_emission_yield_schou_new(q0, int_wall, A2, A3, A4, A5, nw, use_gpu);
    strcpy(ctx->in_species[i], in_species[i]);
  }

  return ctx;
}

void gkyl_bc_emission_release(struct gkyl_bc_emission_ctx *ctx)
{
  for (int i=0; i<ctx->num_species; ++i) {
    gkyl_emission_spectrum_model_release(ctx->spectrum_model[i]);
    gkyl_emission_yield_model_release(ctx->yield_model[i]);
  }
  if (ctx->elastic) gkyl_emission_elastic_model_release(ctx->elastic_model);
  // Release ctx memory.
  gkyl_free(ctx);
}
