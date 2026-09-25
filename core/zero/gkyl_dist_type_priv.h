#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_range.h>

#include <gkyl_dist_type.h>
#include <gkyl_eqn_type.h>


/**
 * Create an distribution projection object.
 *
 * The model_id determines the specific implementation:
 *
 * GKYL_MODEL_DEFAULT       -> Maxwellian
 * GKYL_MODEL_SR            -> Maxwell-Juttner
 * GKYL_MODEL_CANONICAL_PB  -> Canonical-PB
 *
 * @param cdim Configuration-space dimension
 * @param pdim Phase-space dimension
 * @param poly_order Polynomial order
 * @param num_config Number of configuration-space basis functions
 * @param num_phase Number of phase-space basis functions
 * @param model_id Vlasov model identifier
 * @return New generic distribution projection object
 */
struct gkyl_dist_proj_type*
gkyl_vlasov_lte_new(
  int cdim,
  int pdim,
  int poly_order,
  int num_config,
  int num_phase,
  enum gkyl_model_id model_id);
  
/**
 * Create a bi-Maxwellian distribution projection object.
 */
struct gkyl_dist_proj_type*
gkyl_vlasov_bimax_new(
  int cdim,
  int pdim,
  int poly_order,
  int num_config,
  int num_phase,
  enum gkyl_model_id model_id);


/**
 * Create a Kappa distribution projection object.
 */
struct gkyl_dist_proj_type*
gkyl_vlasov_kappa_new(
  int cdim,
  int pdim,
  int poly_order,
  int num_config,
  int num_phase,
  enum gkyl_model_id model_id);