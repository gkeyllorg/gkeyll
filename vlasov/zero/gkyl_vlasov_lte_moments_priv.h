// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_canonical_pb_vars.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>

struct gkyl_vlasov_lte_moments
{
  struct gkyl_basis conf_basis; // Configuration-space basis
  struct gkyl_basis phase_basis; // Phase-space basis
  int num_conf_basis; // Number of configuration-space basis functions
  int vdim; // Number of velocity dimensions
  enum gkyl_model_id model_id; // Enum identifier for model type (e.g., SR, see gkyl_eqn_type.h)
  double mass; // Species mass

  struct gkyl_array *M0; 
  struct gkyl_array *M1i;  
  struct gkyl_array *V_drift;
  struct gkyl_array *V_drift_dot_M1i;  
  struct gkyl_array *pressure;
  struct gkyl_array *temperature;
  struct gkyl_dg_bin_op_mem *mem;

  union {
    // special relativistic Vlasov-Maxwell model
    struct {
      struct gkyl_array *V_drift_sq;
      struct gkyl_array *GammaV;
      struct gkyl_array *GammaV_sq;
      struct gkyl_array *gamma;
      struct gkyl_array *gamma_inv;
      struct gkyl_dg_calc_sr_vars *sr_vars;
    };
    // canonical-pb model
    struct {
      struct gkyl_array *energy_moment;
      struct gkyl_array *h_ij;
      struct gkyl_array *h_ij_inv;
      struct gkyl_array *det_h;
      struct gkyl_dg_calc_canonical_pb_vars *can_pb_vars;
      struct gkyl_array *M1i_cov;  
      struct gkyl_array *V_drift_cov;
    };
  };

  struct gkyl_dg_updater_moment *M0_calc; 
  struct gkyl_dg_updater_moment *M1i_calc;
  struct gkyl_dg_updater_moment *Pcalc;
};
