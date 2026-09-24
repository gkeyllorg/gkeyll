#include <gkyl_emission_elastic_model.h>
#include <math.h>
#include <gkyl_alloc.h>

struct gkyl_emission_elastic_model*
gkyl_emission_elastic_furman_pivi_new(double charge, double P1_inf, double P1_hat,
  double E_hat, double W, double p, bool use_gpu)
{
  struct gkyl_emission_elastic_furman_pivi *model =
    gkyl_malloc(sizeof(struct gkyl_emission_elastic_furman_pivi));

  model->P1_inf = P1_inf;
  model->P1_hat = P1_hat;
  model->E_hat = E_hat;
  model->W = W;
  model->p = p;
  model->elastic.charge = charge;
  model->elastic.function = gkyl_emission_elastic_furman_pivi_yield;

  model->elastic.ref_count = gkyl_ref_count_init(gkyl_emission_elastic_furman_pivi_free);

  return &model->elastic;
}

struct gkyl_emission_elastic_model*
gkyl_emission_elastic_cazaux_new(double charge, double E_f, double phi, bool use_gpu)
{
  struct gkyl_emission_elastic_cazaux *model =
    gkyl_malloc(sizeof(struct gkyl_emission_elastic_cazaux));
  
  model->E_f = E_f;
  model->phi = phi;
  model->elastic.charge = charge;
  model->elastic.function = gkyl_emission_elastic_cazaux_yield;

  model->elastic.ref_count = gkyl_ref_count_init(gkyl_emission_elastic_cazaux_free);

  return &model->elastic;
}

struct gkyl_emission_elastic_model*
gkyl_emission_elastic_constant_new(double charge, double delta, bool use_gpu)
{
  struct gkyl_emission_elastic_constant *model =
   gkyl_malloc(sizeof(struct gkyl_emission_elastic_constant));
  
  model->delta = delta;
  model->elastic.charge = charge;
  model->elastic.function = gkyl_emission_elastic_constant_yield;

  model->elastic.ref_count = gkyl_ref_count_init(gkyl_emission_elastic_constant_free);

  return &model->elastic;
}

struct gkyl_emission_elastic_model*
gkyl_emission_elastic_model_acquire(const struct gkyl_emission_elastic_model* model)
{
  gkyl_ref_count_inc(&model->ref_count);
  return (struct gkyl_emission_elastic_model*) model;
}

void
gkyl_emission_elastic_model_release(const struct gkyl_emission_elastic_model* model)
{
  gkyl_ref_count_dec(&model->ref_count);
}
