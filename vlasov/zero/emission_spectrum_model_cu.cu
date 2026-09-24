/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_emission_spectrum_model.h> 
}

#include <cassert>  

__global__ static void
chung_everhart_set_cu_dev_ptrs(struct gkyl_emission_spectrum_chung_everhart *model)
{
  model->spectrum.normalization = gkyl_emission_spectrum_chung_everhart_norm;
}

__global__ static void
gaussian_set_cu_dev_ptrs(struct gkyl_emission_spectrum_gaussian *model)
{
  model->spectrum.normalization = gkyl_emission_spectrum_gaussian_norm;
}

__global__ static void
maxwellian_set_cu_dev_ptrs(struct gkyl_emission_spectrum_maxwellian *model)
{
  model->spectrum.normalization = gkyl_emission_spectrum_maxwellian_norm;
}

struct gkyl_emission_spectrum_model*
gkyl_emission_spectrum_chung_everhart_cu_dev_new(struct gkyl_emission_spectrum_chung_everhart *model,
  double charge, double phi)
{
  struct gkyl_emission_spectrum_chung_everhart *model_cu =
    (struct gkyl_emission_spectrum_chung_everhart*)
    gkyl_cu_malloc(sizeof(struct gkyl_emission_spectrum_chung_everhart));

  model->spectrum.flags = 0;
  GKYL_SET_CU_ALLOC(model->spectrum.flags);

  gkyl_cu_memcpy(model_cu, model, sizeof(struct gkyl_emission_spectrum_chung_everhart),
    GKYL_CU_MEMCPY_H2D);

  chung_everhart_set_cu_dev_ptrs<<<1,1>>>(model_cu);

  return &model_cu->spectrum;
}

struct gkyl_emission_spectrum_model*
gkyl_emission_spectrum_gaussian_cu_dev_new(struct gkyl_emission_spectrum_gaussian *model,
  double charge, double E_0, double tau)
{
  struct gkyl_emission_spectrum_gaussian *model_cu = (struct gkyl_emission_spectrum_gaussian*)
    gkyl_cu_malloc(sizeof(struct gkyl_emission_spectrum_gaussian));

  model->spectrum.flags = 0;
  GKYL_SET_CU_ALLOC(model->spectrum.flags);

  gkyl_cu_memcpy(model_cu, model, sizeof(struct gkyl_emission_spectrum_gaussian),
    GKYL_CU_MEMCPY_H2D);

  gaussian_set_cu_dev_ptrs<<<1,1>>>(model_cu);

  return &model_cu->spectrum;
}

struct gkyl_emission_spectrum_model*
gkyl_emission_spectrum_maxwellian_cu_dev_new(struct gkyl_emission_spectrum_maxwellian *model,
  double charge, double vt)
{
  struct gkyl_emission_spectrum_maxwellian *model_cu = (struct gkyl_emission_spectrum_maxwellian*)
    gkyl_cu_malloc(sizeof(struct gkyl_emission_spectrum_maxwellian));

  model->spectrum.flags = 0;
  GKYL_SET_CU_ALLOC(model->spectrum.flags);

  gkyl_cu_memcpy(model_cu, model, sizeof(struct gkyl_emission_spectrum_maxwellian),
    GKYL_CU_MEMCPY_H2D);

  maxwellian_set_cu_dev_ptrs<<<1,1>>>(model_cu);

  return &model_cu->spectrum;
}
