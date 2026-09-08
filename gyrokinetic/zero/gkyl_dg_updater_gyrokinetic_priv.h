#pragma once

#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_updater_gyrokinetic.h>

struct gkyl_dg_updater_gyrokinetic {
  bool use_gpu; // Boolean for if the update is performed on GPUs.
  struct gkyl_dg_eqn *eqn_gyrokinetic; // Equation object.
  struct gkyl_hyper_dg *up_gyrokinetic; // Solvers for specific gyrokinetic equation.

  struct gkyl_dg_eqn *eqn_add_apardot_gyrokinetic; // Equation object to only add Apardot contribution.
  struct gkyl_hyper_dg *up_add_apardot_gyrokinetic; // Solver to only add Apardot contribution.

  double gyrokinetic_tm; // Total time spent in computing gyrokinetic equation.
};
