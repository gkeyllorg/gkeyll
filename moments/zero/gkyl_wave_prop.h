#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_comm.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

// Limiters
enum gkyl_wave_limiter {
  GKYL_NO_LIMITER = 1, // to allow default to be 0
  GKYL_MIN_MOD,
  GKYL_SUPERBEE,
  GKYL_VAN_LEER,
  GKYL_MONOTONIZED_CENTERED,
  GKYL_BEAM_WARMING,
  GKYL_ZERO
};

// Wave-splitting to use
enum gkyl_wave_split_type {
  GKYL_WAVE_QWAVE = 0, // default: split jump in Q
  GKYL_WAVE_FWAVE // split jump in F
};

struct gkyl_wave_prop_status {
  int success; // 1 if step worked, 0 otherwise
  double dt_suggested; // suggested time-step
  double max_speed; // max wave speed due to sweep in one direction
};

// Object type for updater
typedef struct gkyl_wave_prop gkyl_wave_prop;

// Parameters for constructor
struct gkyl_wave_prop_inp {
  const struct gkyl_rect_grid *grid; // grid on which to solve equations
  const struct gkyl_wv_eqn *equation; // equation solver
  enum gkyl_wave_limiter limiter; // limiter to use
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  double cfl; // CFL number to use

  bool force_low_order_flux; // if true, only low-order flux is used
  bool check_inv_domain; // flag to indicate if invariant domains are checked

  enum gkyl_wave_split_type split_type; // type of splitting to use

  const struct gkyl_wave_geom *geom; // geometry
  const struct gkyl_comm *comm; // communcator
};

// Some statics from update calls
struct gkyl_wave_prop_stats {
  long n_calls; // number of calls to updater
  long n_bad_advance_calls; // number of calls in which positivity had to be fixed
  long n_bad_cells; // number  of cells fixed
  long n_max_bad_cells; // maximum number of cells fixed in any call
};

/**
 * Create new updater to update equations using wave-propagation
 * algorithm.
 *
 * @param winp Input for creating updater. See gkyl_wave_prop_inp above.
 */
gkyl_wave_prop* gkyl_wave_prop_new(const struct gkyl_wave_prop_inp *winp);

/**
 * Compute wave-propagation update. The update_rng MUST be a sub-range
 * of the range on which the array is defined. That is, it must be
 * either the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param wv Updater object
 * @param tm Current time
 * @param dt time-step
 * @param update_rng Range on which to compute.
 * @param qin Input to updater
 * @param qout Solution at tm+dt
 */
struct gkyl_wave_prop_status gkyl_wave_prop_advance(gkyl_wave_prop *wv,
  double tm, double dt, const struct gkyl_range *update_range,
  struct gkyl_array *phi, const struct gkyl_array *qin, struct gkyl_array *qout);

/**
 * Compute an estimate of maximum stable time-step for given input
 * state 'qin'
 *
 * @param wv Updater object
 * @param qin Input to compute dt for
 * @return maximum stable time-step
 */
double gkyl_wave_prop_max_dt(const gkyl_wave_prop *wv, const struct gkyl_range *update_range,
  const struct gkyl_array *qin);

/**
 * Fetch statics
 *
 * @param wv Updater object
 * @return statics from all calls to this updater
 */
struct gkyl_wave_prop_stats gkyl_wave_prop_stats(const gkyl_wave_prop *wv);
  
/**
 * Delete updater.
 *
 * @param wv Updater to delete.
 */
void gkyl_wave_prop_release(gkyl_wave_prop* wv);
