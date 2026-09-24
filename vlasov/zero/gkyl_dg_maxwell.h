#pragma once

#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>

/**
 * Create a new Maxwell equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param lightSpeed Speed of light
 * @param elcErrorSpeedFactor Factor multiplying lightSpeed for div E correction
 * @param mgnErrorSpeedFactor Factor multiplying lightSpeed for div B correction
 * @return Pointer to Maxwell equation object
 */
struct gkyl_dg_eqn* gkyl_dg_maxwell_new(const struct gkyl_basis* cbasis,
  double lightSpeed, double elcErrorSpeedFactor, double mgnErrorSpeedFactor, bool use_gpu);

/*
 * Create a new Maxwell equation object that lives on NV-GPU.
 *
 * @param cbasis Configuration space basis functions
 * @param lightSpeed Speed of light
 * @param elcErrorSpeedFactor Factor multiplying lightSpeed for div E correction
 * @param mgnErrorSpeedFactor Factor multiplying lightSpeed for div B correction
 * @return Pointer to Maxwell equation object
 */
struct gkyl_dg_eqn* gkyl_dg_maxwell_cu_dev_new(const struct gkyl_basis* cbasis,
  double lightSpeed, double elcErrorSpeedFactor, double mgnErrorSpeedFactor);
