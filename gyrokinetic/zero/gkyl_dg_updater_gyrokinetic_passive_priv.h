#pragma once

#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_updater_gyrokinetic_passive.h>

struct gkyl_dg_updater_gyrokinetic_passive {
  bool use_gpu;
  struct gkyl_dg_eqn *eqn_gyrokinetic_passive;     // Passive GK equation object.
  struct gkyl_hyper_dg *up_gyrokinetic_passive;     // hyper_dg solver (conf-space only).
  double gyrokinetic_passive_tm;                    // Total time in advance calls.
};
