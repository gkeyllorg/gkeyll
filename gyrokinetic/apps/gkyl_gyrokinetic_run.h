#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_multib.h>
struct gkyl_gyrokinetic_time_stepping_inp {
  double t_end; // End time for the simulation
  int num_frames; // Number of output frames. Output every (t_end/num_frames) time units.
  double
    write_phase_freq; // Frequency (in multiples of num_frames) of writing phase-space data. e.g. 0.2 means write every 5 frames.
  int int_diag_calc_num; // Number of integrated diagnostic calculations to do. 100*num_frames means 100 calculations per frame.
  double dt_failure_tol; // Tolerance for small time-step failures. Typical value is 1e-4
  int num_failures_max; // Maximum number of consecutive small time-step failures before aborting simulation. Typical value is 20.
  // Provided by app_args
  bool is_restart; // Is this a restarted simulation?
  int restart_frame; // Frame to restart from.
  int num_steps; // Maximum number of time-steps to take.
};

// Settings for the text which is output into the terminal during the simulation
struct gkyl_gyrokinetic_run_verbosity_inp {
  bool enabled; // Is verbosity enabled? Prints information every time step. Defaults false
  double frequency; // Print information with given frequency. Defaults to 0.1 (every 10 steps)
  bool estimate_completion_time; // Estimate completion time based on current progress. Defaults false
  bool
    disable_timings; // Disable timing information output in the terminal. Defaults false, printing all timing information. This information is always output in the .json file
};

enum gkyl_gyrokinetic_run_app_type {
  GKYL_GK_SINGLEB, // Single-block simulation. Default
  GKYL_GK_MULTIB // Multi-block simulation
};

struct gkyl_gyrokinetic_run_inp {
  enum gkyl_gyrokinetic_run_app_type app_type; // Type of gyrokinetic application to run.
  union {
    struct gkyl_gk app_inp; // Single-block application input struct.
    struct gkyl_gyrokinetic_multib multib_app_inp; // Multi-block application input.
  };
  struct gkyl_gyrokinetic_time_stepping_inp time_stepping; // Timing parameters for the simulation.
  struct gkyl_gyrokinetic_run_verbosity_inp print_verbosity; // Verbosity settings for the simulation.
};

/**
 * Perform the complete gyrokinetic simulation, including initialization, the time loop,
 * and release of resources. If one wants to make modifications at the input file level,
 * copy this function to the input file and modify as needed
 *
 * @param inp Input parameters for the simulation. Includes application input structs.
 */
void gkyl_gyrokinetic_run_simulation(struct gkyl_gyrokinetic_run_inp *inp);
