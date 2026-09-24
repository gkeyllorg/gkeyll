#include <gkyl_emission_yield_model.h>
#include <math.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>

struct gkyl_emission_yield_model*
gkyl_emission_yield_furman_pivi_new(double charge, double deltahat_ts, double Ehat_ts, double t1,
  double t2, double t3, double t4, double s, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_emission_yield_furman_pivi_cu_dev_new(charge, deltahat_ts, Ehat_ts,
      t1, t2, t3, t4, s);
  }
#endif
  struct gkyl_emission_yield_furman_pivi *model =
    gkyl_malloc(sizeof(struct gkyl_emission_yield_furman_pivi));
  
  model->deltahat_ts = deltahat_ts;
  model->Ehat_ts = Ehat_ts;
  model->t1 = t1;
  model->t2 = t2;
  model->t3 = t3;
  model->t4 = t4;
  model->s = s;
  model->yield.charge = charge;
  model->yield.function = gkyl_emission_yield_furman_pivi_yield;

  model->yield.flags = 0;
  GKYL_CLEAR_CU_ALLOC(model->yield.flags);
  model->yield.ref_count = gkyl_ref_count_init(gkyl_emission_yield_furman_pivi_free);

  return &model->yield;
}

struct gkyl_emission_yield_model*
gkyl_emission_yield_schou_new(double charge, double int_wall, double a2, double a3, double a4,
  double a5, double nw, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_emission_yield_schou_cu_dev_new(charge, int_wall, a2, a3, a4, a5, nw);
  }
#endif
  struct gkyl_emission_yield_schou *model = gkyl_malloc(sizeof(struct gkyl_emission_yield_schou));
  
  model->int_wall = int_wall;
  model->a2 = a2;
  model->a3 = a3;
  model->a4 = a4;
  model->a5 = a5;
  model->nw = nw;
  model->yield.charge = charge;
  model->yield.function = gkyl_emission_yield_schou_yield;

  model->yield.flags = 0;
  GKYL_CLEAR_CU_ALLOC(model->yield.flags);
  model->yield.ref_count = gkyl_ref_count_init(gkyl_emission_yield_schou_free);

  return &model->yield;
}

struct gkyl_emission_yield_model*
gkyl_emission_yield_schou_srim_new(double charge, double int_wall, double lorentz_norm, double E0,
  double tau, double alpha, double beta, double gauss_norm, double gauss_E0, double gauss_tau, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_emission_yield_schou_srim_cu_dev_new(charge, int_wall, lorentz_norm, E0, tau, alpha, beta,
    gauss_norm, gauss_E0, gauss_tau);
  }
#endif
  struct gkyl_emission_yield_schou_srim *model = gkyl_malloc(sizeof(struct gkyl_emission_yield_schou_srim));
  
  model->int_wall = int_wall;
  model->lorentz_norm = lorentz_norm;
  model->E0 = E0;
  model->tau = tau;
  model->alpha = alpha;
  model->beta = beta;
  model->gauss_norm = gauss_norm;
  model->gauss_E0 = gauss_E0;
  model->gauss_tau = gauss_tau;
  model->yield.charge = charge;
  model->yield.function = gkyl_emission_yield_schou_srim_yield;

  model->yield.flags = 0;
  GKYL_CLEAR_CU_ALLOC(model->yield.flags);
  model->yield.ref_count = gkyl_ref_count_init(gkyl_emission_yield_schou_srim_free);

  return &model->yield;
}

struct gkyl_emission_yield_model*
gkyl_emission_yield_constant_new(double charge, double delta, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_emission_yield_constant_cu_dev_new(charge, delta);
  }
#endif
  struct gkyl_emission_yield_constant *model =
    gkyl_malloc(sizeof(struct gkyl_emission_yield_constant));
  
  model->delta = delta;
  model->yield.charge = charge;
  model->yield.function = gkyl_emission_yield_constant_yield;

  model->yield.flags = 0;
  GKYL_CLEAR_CU_ALLOC(model->yield.flags);
  model->yield.ref_count = gkyl_ref_count_init(gkyl_emission_yield_constant_free);

  return &model->yield;
}

bool
gkyl_emission_yield_model_is_cu_dev(const struct gkyl_emission_yield_model *model)
{
  return GKYL_IS_CU_ALLOC(model->flags);
}

struct gkyl_emission_yield_model*
gkyl_emission_yield_model_acquire(const struct gkyl_emission_yield_model* model)
{
  gkyl_ref_count_inc(&model->ref_count);
  return (struct gkyl_emission_yield_model*) model;
}

void
gkyl_emission_yield_model_release(const struct gkyl_emission_yield_model* model)
{
  gkyl_ref_count_dec(&model->ref_count);
}
