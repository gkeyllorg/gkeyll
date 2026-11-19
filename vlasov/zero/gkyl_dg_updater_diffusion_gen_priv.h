#pragma once

#include <gkyl_eqn_type.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>

struct gkyl_dg_updater_diffusion_gen {
  struct gkyl_dg_eqn *dgeqn; // Equation object.
  struct gkyl_hyper_dg *hyperdg; // solvers for specific diffusion equation.
  int num_up_dirs; // Number of directions the diffusion acts in.
  int up_dirs[GKYL_MAX_DIM]; // Directions the diffusion acts in.
  bool use_gpu;

  double diffusion_tm; // total time spent in computing diffusion equation.
};
