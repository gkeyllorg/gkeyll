#pragma once

#include <gkyl_block_geom.h>
#include <gkyl_moment.h>

typedef struct gkyl_moment_multib_app gkyl_moment_multib_app;

// Species input per-block
struct gkyl_moment_multib_species_pb {
  int block_id; // block ID
  
  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  void *app_accel_ctx; // context for applied acceleration function
  // pointer to applied acceleration function
  void (*app_accel)(double t, const double *xn, double *app_accel_out, void *ctx);
  bool app_accel_evolve; // set to true if applied acceleration function is time dependent

  void *nT_source_ctx; // context for nT source
  // pointer to user-defined number density and temperature sources
  void (*nT_source_func)(double t, const double *xn, double *fout, void *ctx);
  bool nT_source_set_only_once;
};

// Species input 
struct gkyl_moment_multib_species {
  char name[128]; // species name
  double charge, mass; // charge and mass

  bool is_static; // set to true if moment species does not change in time

  struct gkyl_wv_eqn *equation; // equation object
  enum gkyl_wave_limiter limiter; // limiter to use
  enum gkyl_wave_split_type split_type; // edge splitting to use

  bool force_low_order_flux; // should  we force low-order flux?

  bool duplicate_across_blocks; // set to true if all blocks are identical  
  // species inputs per-block: only one is needed is are_all_blocks_same = true
  const struct gkyl_moment_multib_species_pb *blocks;  

  int num_physical_bcs;
  const struct gkyl_block_physical_bcs *bcs;
};

// Field input per-block
struct gkyl_moment_multib_field_pb {
  int block_id; // block ID
  
  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  void *ext_em_ctx; // context for external electromagnetic fields function
  // pointer to external electromagnetic fields function
  void (*ext_em)(double t, const double *xn, double *ext_em_out, void *ctx);
  bool ext_em_evolve; // set to true if external electromagnetic field function is time dependent
  double t_ramp_E; // linear ramp for turning on external E field
  
  void *app_current_ctx; // context for external electromagnetic fields function
  // pointer to external electromagnetic fields function
  void (*app_current)(double t, const double *xn, double *app_current_out, void *ctx);
  bool app_current_evolve; // set to true if applied current function is time dependent
  double t_ramp_curr; // linear ramp for turning on applied currents

  bool use_explicit_em_coupling; // flag to indicate if using explicit em-coupling  
};

// Field input
struct gkyl_moment_multib_field {
  double epsilon0, mu0;
  double elc_error_speed_fact, mag_error_speed_fact;

  bool is_static; // set to true if field does not change in time

  enum gkyl_wave_limiter limiter; // limiter to use

  bool use_explicit_em_coupling; // flag to indicate if using explicit em-coupling

  bool duplicate_across_blocks; // set to true if all blocks are identical
  // field inputs per-block
  const struct gkyl_moment_multib_field_pb *blocks;

  int num_physical_bcs;
  const struct gkyl_block_physical_bcs *bcs;  
};

// Top-level app parameters: this
struct gkyl_moment_multib {
  char name[128]; // name of app
 // geometry and for blocks in simulation
  struct gkyl_block_geom *block_geom;

 // CFL fraction to use
  double cfl_frac;

  enum gkyl_moment_scheme scheme_type; // scheme to update fluid and moment eqns
  
  enum gkyl_mp_recon mp_recon; // reconstruction scheme to use
  bool skip_mp_limiter; // should MP limiter be skipped?
  bool use_hybrid_flux_kep; // should shock-hybrid scheme be used when using KEP?

  int num_skip_dirs; // number of directions to skip
  int skip_dirs[3]; // directions to skip

 // number of species  
  int num_species;
  // species inputs
  struct gkyl_moment_multib_species species[GKYL_MAX_SPECIES];
 
  // field inputs
  struct gkyl_moment_multib_field field;

 // communicator to used  
  struct gkyl_comm *comm;
};

/**
 * Construct a new moments multi-block app.
 *
 * @param mbinp Multi-block App inputs. See struct docs.
 * @return New multi-block moment app object.
 */
struct gkyl_moment_multib_app* gkyl_moment_multib_app_new(const struct gkyl_moment_multib *mbinp);
                                
/**
 * Compute maximum estimated stable dt wtih current app state. Call
 * after app initialized and after initial conditions set.
 *
 * @param app App object.
 * @retuen maximum estimated stable dt
 */
double gkyl_moment_multib_app_max_dt(gkyl_moment_multib_app* app);

/**
 * Initialize species and field.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void gkyl_moment_multib_app_apply_ic(gkyl_moment_multib_app* app, double t0);

/**
 * Initialize field.
 *
 * @param app App object.
 * @param t0 Time for initial conditions
 */
void gkyl_moment_multib_app_apply_ic_field(gkyl_moment_multib_app* app, double t0);

/**
 * Initialize species.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_moment_multib_app_apply_ic_species(gkyl_moment_multib_app* app, int sidx, double t0);

/**
 * Read field data from specified frame of previous simulation.
 *
 * @param app App object.
 * @param frame Frame number to read from
 * @return Status of read
 */
struct gkyl_app_restart_status gkyl_moment_multib_app_from_frame_field(gkyl_moment_multib_app *app,
  int frame);

/**
 * Read species data from specified frame of previous simulation.
 *
 * @param app App object.
 * @param sidx Index of species to read
 * @param frame Frame number to read from
 * @return Status of read
 */
struct gkyl_app_restart_status gkyl_moment_multib_app_from_frame_species(gkyl_moment_multib_app *app,
  int sidx, int frame);

/**
 * Write output to console: this is mainly for diagnostic messages the
 * driver code wants to write to console. It accounts for parallel
 * output by not messing up the console with messages from each rank.
 *
 * @param app App object
 * @param fp File pointer for open file for output
 * @param fmt Format string for console output
 * @param argp Objects to write
 */
void gkyl_moment_multib_app_cout(const gkyl_moment_multib_app* app, FILE *fp, const char *fmt, ...);

/**
 * Write block topology to file.
 * 
 * @param app App object.
 */
void gkyl_moment_multib_app_write_topo(const gkyl_moment_multib_app* app);

/**
 * Write field and species data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_moment_multib_app_write(const gkyl_moment_multib_app* app, double tm, int frame);

/**
 * Write field data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_moment_multib_app_write_field(const gkyl_moment_multib_app *app, double tm, int frame);

/**
 * Write species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_moment_multib_app_write_species(const gkyl_moment_multib_app* app, int sidx, double tm, int frame);

/**
 * Write field energy to file.
 *
 * @param app App object.
 */
void gkyl_moment_multib_app_write_field_energy(gkyl_moment_multib_app *app);

/**
 * Write integrated moments to file.
 *
 * @param app App object.
 */
void gkyl_moment_multib_app_write_integrated_mom(gkyl_moment_multib_app *app);

/**
 * Write stats to file. Data is written in json format.
 *
 * @param app App object.
 */
void gkyl_moment_multib_app_stat_write(const gkyl_moment_multib_app *app);

/**
 * Advance simulation by a suggested time-step 'dt'. The dt may be too
 * large in which case method will attempt to take a smaller time-step
 * and also return it as the 'dt_actual' field of the status
 * object. If the suggested time-step 'dt' is smaller than the largest
 * stable time-step the method will use the smaller value instead,
 * returning the larger time-step in the 'dt_suggested' field of the
 * status object. If the method fails to find any stable time-step
 * then the 'success' flag will be set to 0. At that point the calling
 * code must abort the simulation as this signals a catastrophic
 * failure and the simulation can't be safely continued.
 * 
 * @param app App object.
 * @param dt Suggested time-step to advance simulation
 * @return Status of update.
 */
struct gkyl_update_status gkyl_moment_multib_update(gkyl_moment_multib_app *app, double dt);

/**
 * Calculate integrated field energy. The "calc" method computes the
 * integrated moments and stores it. The "get" method returns the
 * values in the vals array, without storing it.
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_moment_multib_app_calc_field_energy(gkyl_moment_multib_app *app, double tm);
void gkyl_moment_multib_app_get_field_energy(gkyl_moment_multib_app *app, double *vals);

/**
 * Calculate integrated moments.
 *
 * @param app App object.
 * @param tm Time at which integrated diagnostic are to be computed
 */
void gkyl_moment_multib_app_calc_integrated_mom(gkyl_moment_multib_app *app, double tm);

/**
 * Return simulation statistics.
 * 
 * @return Return statistics.
 */
struct gkyl_moment_stat gkyl_moment_multib_app_stat(gkyl_moment_multib_app *app);

/**
 * Free moment app.
 *
 * @param app App to release.
 */
void gkyl_moment_multib_app_release(gkyl_moment_multib_app* app);  
