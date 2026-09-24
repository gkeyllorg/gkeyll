#pragma once

#include <gkyl_wv_eqn.h>
#include <gkyl_gr_spacetime.h>

// Type of Riemann-solver to use:
enum gkyl_wv_vacuum_einstein_rp {
  WV_VACUUM_EINSTEIN_RP_HLL = 0, // Default (HLL fluxes).
  WV_VACUUM_EINSTEIN_RP_LAX,
};

// Input context, packaged as a struct.
struct gkyl_wv_vacuum_einstein_inp {
  double excision_threshold; // Excision threshold (lapse).
  enum gkyl_spacetime_slicing spacetime_slicing; // Spacetime slicing condition.
  enum gkyl_spacetime_evolution spacetime_evolution; // Spacetime evolution system.

  enum gkyl_wv_vacuum_einstein_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new vacuum Einstein equations object using the Bona-Masso formalism.
*
* @param excision_threshold Excision threshold (lapse).
* @param spacetime_slicing Spacetime slicing condition.
* @param spacetime_evolution Spacetime evolution system.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the vacuum Einstein equations object using the Bona-Masso formalism.
*/
struct gkyl_wv_eqn*
gkyl_wv_vacuum_einstein_new(double excision_threshold, enum gkyl_spacetime_slicing spacetime_slicing, enum gkyl_spacetime_evolution spacetime_evolution,
  bool use_gpu);

/**
* Create a new vacuum Einstein equations object using the Bona-Masso formalism, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the vacuum Einstein equations object using the Bona-Masso formalism.
*/
struct gkyl_wv_eqn*
gkyl_wv_vacuum_einstein_inew(const struct gkyl_wv_vacuum_einstein_inp* inp);

/**
* Get excision threshold.
*
* @param eqn Vacuum Einstein equations object using the Bona-Masso formalism.
* @return Excision threshold (lapse).
*/
double
gkyl_wv_vacuum_einstein_excision_threshold(const struct gkyl_wv_eqn* eqn);

/**
* Get spacetime slicing condition.
*
* @param eqn Vacuum Einstein equations object using the Bona-Masso formalism.
* @return Spacetime slicing condition.
*/
enum gkyl_spacetime_slicing
gkyl_wv_vacuum_einstein_spacetime_slicing(const struct gkyl_wv_eqn* eqn);

/**
* Get spacetime evolution system.
*
* @param eqn Vacuum Einstein equations object using the Bona-Masso formalism.
* @return Spacetime evolution system.
*/
enum gkyl_spacetime_evolution
gkyl_wv_vacuum_einstein_spacetime_evolution(const struct gkyl_wv_eqn* eqn);