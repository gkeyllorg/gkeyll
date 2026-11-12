#pragma once

/**
 *
 * SUNDIALS NVECTOR wrapper of Gkeyll data structure and methods.
 *
 */

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_comm.h>

// Time stepping method options.
enum gkyl_sundials_rk_method {
  GKYL_SUNDIALS_METHOD_NONE          = 0,
  GKYL_SUNDIALS_METHOD_RK_SSP_3_3    = 1 << 0, // Gkeyll's native 3rd order 3-stage SSP RK method.
  GKYL_SUNDIALS_METHOD_LSRK_SSP_S_2  = 1 << 1, // Optimal 2nd order SSP RK method s stages (>=2).
  GKYL_SUNDIALS_METHOD_LSRK_SSP_S_3  = 1 << 2, // Optimal 3rd order s-stage SSP RK method.
  GKYL_SUNDIALS_METHOD_LSRK_SSP_10_4 = 1 << 3, // Optimal 4th order 10-stage SSP RK method.
  GKYL_SUNDIALS_METHOD_LSRK_RKC_2    = 1 << 4, // 2nd order Runge-Kutta-Chebyshev (RKC).
  GKYL_SUNDIALS_METHOD_LSRK_RKL_2    = 1 << 5, // 2nd order Runge-Kutta-Legendre (RKL).
  GKYL_SUNDIALS_METHOD_MRI_EXTSTS    = 1 << 6, // MRI extended STS (combine with RKC_2 or RKL_2; default RKL_2).
};

// Operator splitting method, see
//   https://sundials.readthedocs.io/en/latest/arkode/Usage/SplittingStep/SplittingStepCoefficients.html#arkode-usage-splittingstep-splittingstepcoefficients
// for more info.
enum gkyl_sundials_opsplit_method {
  GKYL_SUNDIALS_OPSPLIT_METHOD_NONE = 0, // No operator splitting.
  GKYL_SUNDIALS_OPSPLIT_METHOD_LIE_TROTTER_1_1_2, // Default.
  GKYL_SUNDIALS_OPSPLIT_METHOD_STRANG_2_2_2,
  GKYL_SUNDIALS_OPSPLIT_METHOD_BEST_2_2_2,
  GKYL_SUNDIALS_OPSPLIT_METHOD_SUZUKI_3_3_2,
  GKYL_SUNDIALS_OPSPLIT_METHOD_RUTH_3_3_2,
  GKYL_SUNDIALS_OPSPLIT_METHOD_YOSHIDA_4_4_2,
  GKYL_SUNDIALS_OPSPLIT_METHOD_YOSHIDA_8_6_2,
};

// Context for functions that are app-specific and/or
// need to handle the app pointer as a void *.
struct gkyl_sundials_app_ctx {
  void *app_ptr; // Gkeyll app.
  void *fdot_args_ptr; // Arguments to df/dt calculation.
  void *arkode_mem_ssprk; // Memory for ARKODE SSP-RK stepper.
  void *arkode_mem_sts; // Memory for ARKODE STS stepper.
  void *arkode_mem_outer; // Memory for the outermost ARKODE stepper (e.g. operator split).
  // Function that computes weight for the error norm.
  int (*error_wgt_func)(void *app_gen, const struct gkyl_array *xarr, struct gkyl_array *wgt, struct gkyl_range *local_range);
  // Function that computes df/dt for all operators, or for SSP-RK or STS operators separately.
  double (*dfdt_func)(void *app_gen, double t_curr, void *fdot_args_gen);
  double (*dfdt_ssprk_func)(void *app_gen, double t_curr, void *fdot_args_gen);
  double (*dfdt_sts_func)(void *app_gen, double t_curr, void *fdot_args_gen);
  // Operations performed at the beginning/end of an outer step.
  void (*pre_process_step_opsplit_func)(void *app_gen, double tcurr, double dt, void *fdot_args_gen);
  void (*post_process_step_opsplit_func)(void *app_gen, double tcurr, double dt, void *fdot_args_gen);
  // Operations performed at the beginning/end of an SSP-RK step/stage or failed step.
  void (*pre_process_step_ssprk_func)(void *app_gen, double tcurr, double dt, void *fdot_args_gen);
  void (*pre_process_rk_stage_ssprk_func)(void *app_gen, double tcurr, double dt, void *fdot_args_gen, int stage_idx, int num_stages);
  void (*post_process_rk_stage_ssprk_func)(void* app_gen, double tcurr, double dt, void *fdot_args, int stage_idx, int num_stages);
  void (*post_process_step_ssprk_func)(void *app_gen, double tcurr, double dt, void *fdot_args_gen);
  // Operations performed at the beginning/end of an STS step/stage or failed step.
  void (*pre_process_step_sts_func)(void *app_gen, double tcurr, double dt, void *fdot_args_gen);
  void (*pre_process_rk_stage_sts_func)(void *app_gen, double tcurr, double dt, void *fdot_args_gen, int stage_idx, int num_stages);
  void (*post_process_rk_stage_sts_func)(void* app_gen, double tcurr, double dt, void *fdot_args, int stage_idx, int num_stages);
  void (*post_process_step_sts_func)(void *app_gen, double tcurr, double dt, void *fdot_args_gen);
  // Operations performed at the end of an ExtSTS outer stage.
  void (*post_process_stage_extsts_func)(void* app_gen, double tcurr, double dt, void *fdot_args, int stage_idx, int num_stages);
  // Function that reduces a local dt across MPI processes.
  double (*reduce_dt_func)(void *app_gen, double t_curr, double dt_local);
  // Objects below are private.
  double dt_local; // CFL constrained time step in local MPI process.
  double dt_local_ssprk, dt_local_sts; // dt_local for SSP-RK and STS components.
};

// Sundials inputs specified in input file.
struct gkyl_sundials_stepper_inp {
  double rel_tol; // Relative tolerance.
  double abs_tol; // Absolute tolerance.
  long max_steps; // Maximum number of steps (default: 1e5).
  unsigned int num_stages; // Number of stages in a step, for methods with fixed number of stages.
  unsigned int max_num_stages; // Maximum number of stages, for methods with adaptive number of stages (default: 200).
  enum gkyl_sundials_rk_method rk_method; // Time stepping method (default: GKYL_SUNDIALS_METHOD_RK_SSP_3_3).
  enum gkyl_sundials_opsplit_method opsplit_method; // Operator split type (default: GKYL_SUNDIALS_OPSPLIT_METHOD_LIE_TROTTER_1_1_2)
  int mri_coupling; // ARKODE_MRITableID of the outer ExtSTS method (default: 0 = SUNDIALS default).
  bool fixed_step; // Fixed time step set by the CFL limit of the SSP-RK terms, instead of adaptive steps (default: false; the operator split always uses fixed steps; not available for STS-only methods).
  struct gkyl_sundials_nvec *gsnv; // Input NVECTOR.
  double t_curr; // Current simulation time.
  struct gkyl_sundials_app_ctx *app_ctx; // Context with app-specific data and functions.
  bool print_arkode_stats; // Whether to print ARKODE statistics at the end of the run. 
  // Dominant eigenvalue estimator inputs.
  bool dee_by_gkeyll; // =true Gkeyll's DEE, =false SUNDIALS' (default: false).
  unsigned int dee_max_iter; // Maximum number of iterations (default: 1e3).
  double dee_rel_tol; // Relative tolerance (default: 0.01).
  unsigned int dee_num_init_warmups; // Number of initial warm ups (default: 100).
  unsigned int dee_num_succ_warmups; // Number of succeeding warm ups (default: 0).
  unsigned int dee_frequency; // Number of steps between DEEs (default: 10).
  double dee_safety_fac; // Safety factor setting the effective dominant eigenvalue (default: 1.01).
};

//
// SUNDIALS module methods.
//

typedef struct gkyl_sundials gkyl_sundials; // Sundials object.

/**
 * Create a new SUNDIALS object (e.g. a SUNDIALS context).
 *
 * @param use_gpu Whether data wrapped by SUNDIALS NVECTORS lives on the GPU.
 * @return A new SUNDIALS object.
 */
struct gkyl_sundials* gkyl_sundials_new(bool use_gpu);

/**
 * Check if rk_method equals:
 *   - opt_A if opt_B=0.
 *   - opt_B if opt_A=0.
 *   - opt_A | opt_B.
 *
 * @param rk_method Input method option(s).
 * @param opt_A First option to check for.
 * @param opt_B Second option to check for.
 * @return Whether input method meets desired condition.
 */
bool gkyl_sundials_check_rk_method(enum gkyl_sundials_rk_method rk_method,
  enum gkyl_sundials_rk_method opt_A, enum gkyl_sundials_rk_method opt_B);

/**
 * Initialize the SUNDIALS time stepper.
 *
 * @param gksun SUNDIALS object.
 * @param inp SUNDIALS inputs for time stepper.
 */
void gkyl_sundials_stepper_init(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp);

/**
 * Reset the ARKODE time stepper given a new time and state vector.
 *
 * @param gksun SUNDIALS object.
 * @param time New simulation time.
 * @param gsmanynv New state vectors (as ManyNVector).
 * @param gsmanynv_buff Buffer array, as big as state vector array.
 */
void gkyl_sundials_arkode_reset(struct gkyl_sundials *gksun, double time,
  struct gkyl_sundials_nvec *gsmanynv, struct gkyl_sundials_nvec *gsmanynv_buff);

/**
 * Disable adaptive time stepping (for outer stepper if using operator split)
 * and set the fixed time step.
 *
 * @param gksun SUNDIALS object.
 * @param dt Time step size.
 */
void gkyl_sundials_set_fixed_step(struct gkyl_sundials *gksun, double dt);

/**
 * Disable adaptive time stepping and set the fixed time step for SSP-RK.
 *
 * @param gksun SUNDIALS object.
 * @param dt Time step size.
 */
void gkyl_sundials_set_fixed_step_ssprk(struct gkyl_sundials *gksun, double dt);

/**
 * Disable adaptive time stepping and set the fixed time step for SSP-RK.
 *
 * @param gksun SUNDIALS object.
 * @param dt Time step size.
 */
void gkyl_sundials_set_fixed_step_sts(struct gkyl_sundials *gksun, double dt);

/**
 * Return the CFL stable time step size (dt) for the df/dt operator.
 *
 * @param gksun SUNDIALS object.
 * @return CFL stable time step size.
 */
double gkyl_sundials_get_cfl_dt(struct gkyl_sundials *gksun);

/**
 * Return the CFL stable time step size (dt) for the SSP-RK df/dt operator.
 *
 * @param gksun SUNDIALS object.
 * @return CFL stable time step size.
 */
double gkyl_sundials_get_cfl_dt_ssprk(struct gkyl_sundials *gksun);

/**
 * Return the CFL stable time step size (dt) for the STS df/dt operator.
 *
 * @param gksun SUNDIALS object.
 * @return CFL stable time step size.
 */
double gkyl_sundials_get_cfl_dt_sts(struct gkyl_sundials *gksun);

/**
 * Evolve the solution contained in a given Nvector from
 * the current to the new time.
 *
 * @param gksun SUNDIALS object.
 * @param t_new New time to evolve to.
 * @param gsnv Gkeyll wrapper of SUNDIALS NVECTOR.
 * @param t_curr Current time.
 */
int gkyl_sundials_evolve(struct gkyl_sundials *gksun, double t_new,
  struct gkyl_sundials_nvec *gsnv, double *t_curr);

/**
 * Fetch the number of error test_failures so far.
 *
 * @param gksun SUNDIALS object.
 * @return Number of failures.
 */
long gkyl_sundials_get_num_error_test_failures(struct gkyl_sundials *gksun);

/**
 * Fetch the number of evaluatons of the RHS (i.e. df/dt) function.
 *
 * @param gksun SUNDIALS object.
 * @return Number of RHS evaluations.
 */
long gkyl_sundials_get_num_rhs_evals(struct gkyl_sundials *gksun);

/**
 * Fetch size of the time step used in the last step (of the outer stepper if
 * using operator splitting).
 *
 * @param gksun SUNDIALS object.
 * @return Time step size.
 */
double gkyl_sundials_get_last_dt(struct gkyl_sundials *gksun);

/**
 * Fetch size of the time step used in the last step of the SSP-RK stepper.
 *
 * @param gksun SUNDIALS object.
 * @return Time step size.
 */
double gkyl_sundials_get_last_dt_ssprk(struct gkyl_sundials *gksun);

/**
 * Fetch size of the time step used in the last step of the STS stepper.
 *
 * @param gksun SUNDIALS object.
 * @return Time step size.
 */
double gkyl_sundials_get_last_dt_sts(struct gkyl_sundials *gksun);

/**
 * Fetch size of the current time step (of the outer stepper if
 * using operator splitting).
 *
 * @param gksun SUNDIALS object.
 * @return Time step size.
 */
double gkyl_sundials_get_current_dt(struct gkyl_sundials *gksun);

/**
 * Fetch size of the current time step in the SSP-RK stepper.
 *
 * @param gksun SUNDIALS object.
 * @return Time step size.
 */
double gkyl_sundials_get_current_dt_ssprk(struct gkyl_sundials *gksun);

/**
 * Fetch size of the current time step in the STS stepper.
 *
 * @param gksun SUNDIALS object.
 * @return Time step size.
 */
double gkyl_sundials_get_current_dt_sts(struct gkyl_sundials *gksun);

/**
 * Check if sundials is using an operator split approach (e.g. combining SSP-RK
 * and STS).
 *
 * @param gksun SUNDIALS object.
 * @return If using operator split.
 */
bool gkyl_sundials_use_operator_split(struct gkyl_sundials *gksun);

/**
 * Check if rk_method requires operator splitting.
 *
 * @param rk_method Method enum to be checked.
 * @return If using operator split.
 */
bool gkyl_sundials_operator_split_in_method(enum gkyl_sundials_rk_method rk_method);

/**
 * Free resources associates with a SUNDIALS object.
 *
 * @param gksun SUNDIALS object to be freed.
 */
void gkyl_sundials_release(struct gkyl_sundials *gksun);

//
// SUNDIALS NVECTOR methods.

typedef struct gkyl_sundials_nvec gkyl_sundials_nvec; // Sundials nvector.

/**
 * Create a new NVECTOR wrapping for a given gkyl_array.
 *
 * @param gksun SUNDIALS object.
 * @param arr Gkeyll array to wrap.
 * @param comm Gkeyll communicator object.
 * @param local_range Local range to loop over arr.
 * @param no_error_norm Whether to exclude from dt or error norm calculation.
 * @param not_stepped Whether to exclude from time stepping.
 * @return A new NVECTOR.
 */
struct gkyl_sundials_nvec* gkyl_sundials_nvec_new(struct gkyl_sundials *gksun,
  struct gkyl_array *arr, struct gkyl_comm *comm, struct gkyl_range *local_range,
  bool no_error_norm, bool not_stepped);

/**
 * Fetch the Gkeyll array wrapped by an NVECTOR.
 *
 * @param nvin Input NVECTOR.
 * @return Gkeyll array wrapped by NVECTOR.
 */
struct gkyl_array* gkyl_sundials_nvec_get_array(struct gkyl_sundials_nvec *nvin);

/**
 * Destroy SUNDIALS NVECTOR.
 *
 * @param nvin NVECTOR to be destroyed.
 */
void gkyl_sundials_nvec_release(struct gkyl_sundials_nvec *nvin);

/**
 * Create a new ManyNVECTOR wrapping of all our Nvectors (one for each
 * gkyl_array that is evolved via time integration of a PDE).
 *
 * @param gksun SUNDIALS object.
 * @param num_nvector Number of NVectors to wrap.
 * @param gsnv_arr Array of Gkeyll nvector objects.
 * @return A new NVECTOR, which holds the ManyVector.
 */
struct gkyl_sundials_nvec* gkyl_sundials_many_nvec_new(struct gkyl_sundials *gksun,
  int num_nvector, struct gkyl_sundials_nvec* gsnv_arr[]);

/**
 * Return the number of sub-vectors in a Sundials ManyNVector object.
 *
 * @param gksun SUNDIALS object.
 * @param gsmanynv A Gkeyll-Sundials ManyNVector.
 * @return Number of subvectors.
 */
int gkyl_sundials_many_nvec_get_num_subvec(struct gkyl_sundials *gksun,
  struct gkyl_sundials_nvec* gsmanynv);

/**
 * Destroy SUNDIALS ManyNVECTOR.
 *
 * @param gsmanynv ManyNVECTOR to be destroyed.
 */
void gkyl_sundials_many_nvec_release(struct gkyl_sundials_nvec *gsmanynv);

