// Private header: not for direct use
#pragma once

#include "gkyl_vlasov_dist_moments.h"
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_dist_type.h>

// Forward declare
struct gkyl_vlasov_dist_moments;

/**
* Function pointer type to compute the density moment
*/
typedef void (*dist_density_moment_t)(struct gkyl_vlasov_dist_moments *dist_moms,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local,
  const struct gkyl_array *fin, struct gkyl_array *density_out);

/**
* Function pointer type to compute the moments required by the selected distribution
*/
typedef void (*dist_moments_t)(struct gkyl_vlasov_dist_moments *dist_moms,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local,
  const struct gkyl_array *fin, struct gkyl_array *moms_out);

/**
* Function pointer type to create a distribution moments updater
*/
typedef struct gkyl_vlasov_dist_moments* (*dist_moments_new_t)(
  const struct gkyl_vlasov_dist_moments_inp *inp);

/**
* Function pointer type to release a distribution moments updater.
*/
typedef void (*dist_moments_release_t)(struct gkyl_vlasov_dist_moments *dist_moms);



struct gkyl_vlasov_dist_moments
{
  struct gkyl_basis conf_basis; // Configuration-space basis
  struct gkyl_basis phase_basis; // Phase-space basis
  int num_conf_basis; // Number of configuration-space basis functions
  int vdim; // Number of velocity dimensions
  int num_mom; // Number of moments required by the selected distribution
  enum gkyl_distribution_projection dist_id;  // Selected distribution
  enum gkyl_model_id model_id; // Enum identifier for model type (e.g., SR, see gkyl_eqn_type.h)
  double mass; // Species mass
  bool use_gpu;
  dist_density_moment_t density_moment; // Density moment calculation
  dist_moments_t moments; // Distribution moments calculation
  dist_moments_release_t release; // Release distribution moments updater 
};

/**
* Create an LTE distribution moments updater.
*
* @param inp Input parameters.
* @return New generic distribution moments updater.
*/
struct gkyl_vlasov_dist_moments * gkyl_vlasov_lte_moments_new(const struct gkyl_vlasov_dist_moments_inp *inp);



// struct gkyl_vlasov_dist_moments * gkyl_vlasov_multitemp_max_moments_new(const struct gkyl_vlasov_dist_moments_inp *inp);
// struct gkyl_vlasov_dist_moments * gkyl_vlasov_kappa_moments_new(const struct gkyl_vlasov_dist_moments_inp *inp);
