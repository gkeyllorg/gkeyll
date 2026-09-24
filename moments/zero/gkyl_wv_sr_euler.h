#pragma once

#include <gkyl_wv_eqn.h>

/**
 * Create a new SR Euler equation object.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @return Pointer to SR Euler equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_sr_euler_new(double gas_gamma);

/**
 * Get gas adiabatic constant.
 * 
 * @param wv SR Euler equation object
 * @return Get gas adiabatic constant
 */
double gkyl_wv_sr_euler_gas_gamma(const struct gkyl_wv_eqn* wv);
