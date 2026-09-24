/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_emission_yield_model.h> 
}

#include <cassert>  

__global__ static void
furman_pivi_set_cu_dev_ptrs(struct gkyl_emission_yield_furman_pivi *model)
{
  model->yield.function = gkyl_emission_yield_furman_pivi_yield;
}

__global__ static void
schou_set_cu_dev_ptrs(struct gkyl_emission_yield_schou *model)
{
  model->yield.function = gkyl_emission_yield_schou_yield;
}

__global__ static void
schou_srim_set_cu_dev_ptrs(struct gkyl_emission_yield_schou_srim *model)
{
  model->yield.function = gkyl_emission_yield_schou_srim_yield;
}

__global__ static void
constant_set_cu_dev_ptrs(struct gkyl_emission_yield_constant *model)
{
  model->yield.function = gkyl_emission_yield_constant_yield;
}

struct gkyl_emission_yield_model*
gkyl_emission_yield_furman_pivi_cu_dev_new(double charge, double deltahat_ts, double Ehat_ts,
  double t1, double t2, double t3, double t4, double s)
{
  struct gkyl_emission_yield_furman_pivi *model = (struct gkyl_emission_yield_furman_pivi*)
    gkyl_malloc(sizeof(struct gkyl_emission_yield_furman_pivi));

  model->deltahat_ts = deltahat_ts;
  model->Ehat_ts = Ehat_ts;
  model->t1 = t1;
  model->t2 = t2;
  model->t3 = t3;
  model->t4 = t4;
  model->s = s;
  model->yield.charge = charge;

  model->yield.flags = 0;
  GKYL_SET_CU_ALLOC(model->yield.flags);
  model->yield.ref_count = gkyl_ref_count_init(gkyl_emission_yield_furman_pivi_free);

  struct gkyl_emission_yield_furman_pivi *model_cu = (struct gkyl_emission_yield_furman_pivi*)
    gkyl_cu_malloc(sizeof(struct gkyl_emission_yield_furman_pivi));
  gkyl_cu_memcpy(model_cu, model, sizeof(struct gkyl_emission_yield_furman_pivi),
    GKYL_CU_MEMCPY_H2D);

  furman_pivi_set_cu_dev_ptrs<<<1,1>>>(model_cu);

  model->yield.on_dev = &model_cu->yield;

  return &model->yield;
}

struct gkyl_emission_yield_model*
gkyl_emission_yield_schou_cu_dev_new(double charge, double int_wall, double a2, double a3,
  double a4, double a5, double nw)
{
  struct gkyl_emission_yield_schou *model = (struct gkyl_emission_yield_schou*) 
    gkyl_malloc(sizeof(struct gkyl_emission_yield_schou));

  model->int_wall = int_wall;
  model->a2 = a2;
  model->a3 = a3;
  model->a4 = a4;
  model->a5 = a5;
  model->nw = nw;
  model->yield.charge = charge;

  model->yield.flags = 0;
  GKYL_SET_CU_ALLOC(model->yield.flags);
  model->yield.ref_count = gkyl_ref_count_init(gkyl_emission_yield_schou_free);

  struct gkyl_emission_yield_schou *model_cu = (struct gkyl_emission_yield_schou*)
    gkyl_cu_malloc(sizeof(struct gkyl_emission_yield_schou));
  gkyl_cu_memcpy(model_cu, model, sizeof(struct gkyl_emission_yield_schou), GKYL_CU_MEMCPY_H2D);

  schou_set_cu_dev_ptrs<<<1,1>>>(model_cu);

  model->yield.on_dev = &model_cu->yield;

  return &model->yield;
}

struct gkyl_emission_yield_model*
gkyl_emission_yield_schou_srim_cu_dev_new(double charge, double int_wall, double lorentz_norm, double E0, double tau,
  double alpha, double beta, double gauss_norm, double gauss_E0, double gauss_tau)
{
  struct gkyl_emission_yield_schou_srim *model = (struct gkyl_emission_yield_schou_srim*) 
    gkyl_malloc(sizeof(struct gkyl_emission_yield_schou_srim));

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

  model->yield.flags = 0;
  GKYL_SET_CU_ALLOC(model->yield.flags);
  model->yield.ref_count = gkyl_ref_count_init(gkyl_emission_yield_schou_srim_free);

  struct gkyl_emission_yield_schou_srim *model_cu = (struct gkyl_emission_yield_schou_srim*)
    gkyl_cu_malloc(sizeof(struct gkyl_emission_yield_schou_srim));
  gkyl_cu_memcpy(model_cu, model, sizeof(struct gkyl_emission_yield_schou_srim), GKYL_CU_MEMCPY_H2D);

  schou_srim_set_cu_dev_ptrs<<<1,1>>>(model_cu);

  model->yield.on_dev = &model_cu->yield;

  return &model->yield;
}

struct gkyl_emission_yield_model*
gkyl_emission_yield_constant_cu_dev_new(double charge, double delta)
{
  struct gkyl_emission_yield_constant *model = (struct gkyl_emission_yield_constant*)
    gkyl_malloc(sizeof(struct gkyl_emission_yield_constant));

  model->delta = delta;
  model->yield.charge = charge;

  model->yield.flags = 0;
  GKYL_SET_CU_ALLOC(model->yield.flags);
  model->yield.ref_count = gkyl_ref_count_init(gkyl_emission_yield_constant_free);

  struct gkyl_emission_yield_constant *model_cu = (struct gkyl_emission_yield_constant*)
    gkyl_cu_malloc(sizeof(struct gkyl_emission_yield_constant));
  gkyl_cu_memcpy(model_cu, model, sizeof(struct gkyl_emission_yield_constant), GKYL_CU_MEMCPY_H2D);

  constant_set_cu_dev_ptrs<<<1,1>>>(model_cu);

  model->yield.on_dev = &model_cu->yield;

  return &model->yield;
}
