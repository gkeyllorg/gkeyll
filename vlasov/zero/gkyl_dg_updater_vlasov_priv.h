#pragma once

#include <gkyl_eqn_type.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>

struct gkyl_dg_updater_vlasov {
  enum gkyl_model_id model_id; // Identifier for model (e.g., SR, PKPM, see gkyl_eqn_type.h)
  enum gkyl_field_id field_id; // Identifier for field type (e.g., Maxwell's, Poisson, see gkyl_eqn_type.h)
  bool use_gpu; // Boolean for if the update is performed on GPUs
  struct gkyl_dg_eqn *eqn_vlasov; // Equation object
  struct gkyl_hyper_dg *hdg_vlasov; // solvers for specific Vlasov equation

  double vlasov_tm; // total time spent in computing vlasov equation
};
