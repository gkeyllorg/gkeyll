#pragma once

// Private header file for SUNDIALS nvector wrapper.

#include <gkyl_sundials.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector.h>
#include <nvector/nvector_manyvector.h>
#include <sundials/sundials_types.h> // Definition of type sunrealtype.
#include <arkode/arkode_lsrkstep.h> // Access to LSRKStep.
#include <arkode/arkode_erkstep.h> // Access to ERKStep.
#include <sundomeigest/sundomeigest_power.h> // Power iteration dominant eigen value estimator.
#include <arkode/arkode_splittingstep.h> // Operator splitting.
#include <arkode/arkode_mristep.h> // Extended STS (MRIStep).
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h> // Time step controller (ExtSTS).

// Give access to the Gkeyll vector from within the NVECTOR.
#define NV_CONTENT_GKZ(v) ((N_VectorContent_Gkeyll)(v->content))

struct gkyl_sundials_reduce_mem
{
  double *red_local, *red_global; // Result of local and global reductions.
  double *red_global_ho; // Result of global reduction on the host.
};

// Gkeyll implementation of N_Vector.
struct _N_VectorContent_Gkeyll
{
  struct gkyl_array *arr; // Array wrapped by this vector.
  sunbooleantype own_vector; // Ownership of the gkyl_array.
  sunbooleantype use_gpu; // Whether data is on GPU.
  struct gkyl_comm *comm; // Communicator.
  struct gkyl_range *local_range; // Local range.
  double global_volume; // Number of cells over all ranks (normalization of the WRMS norm).
  struct gkyl_sundials_reduce_mem *red_mem; // Memory for reductions.
  bool no_error_norm; // Whether to exclude it from the calculation of dt or error norms.
  bool not_stepped; // Whether to exclude it from time stepping.
};

typedef struct _N_VectorContent_Gkeyll* N_VectorContent_Gkeyll;

struct gkyl_sundials_nvec
{
  N_Vector nvec; // SUNDIALS NVector.
  struct gkyl_sundials_reduce_mem red_mem; // Memory for reductions.
};

struct gkyl_sundials
{
  struct gkyl_sundials_stepper_inp *stepper_inp; // SUNDIALS stepper inputs.
  SUNContext sunctx; // Sundials context.
  bool use_gpu; // Whether to run on GPU.
  bool has_ssprk; // Whether stepper has SSP-RK.
  bool has_sts; // Whether stepper has STS.
  void *arkode_mem_ssprk; // Memory for ARKODE SSP-RK stepper.
  void *arkode_mem_sts; // Memory for ARKODE STS stepper.
  struct gkyl_sundials_app_ctx *app_ctx; // App-specific context.
  SUNDomEigEstimator dom_eig_est; // Dominant eigenvalue estimator.
  bool is_opsplit; // Whether we are using an operator split approach.
  bool is_extsts; // Whether we are using an extended STS (MRI) method.
  SUNAdaptController adapt_controller; // Time step controller of adaptive ExtSTS (owned here; NULL if unused).
  SUNStepper stepper_ssprk, stepper_sts; // Steppers for operator splitting.
  void *arkode_mem_outer; // Memory for the outermost ARKODE stepper (e.g. operator split).
  // Methods assigned by specific apps.
  int (*snvec_efun_cell_norm_func)(N_Vector manyx, N_Vector manyw, void *ctx);

  int (*dfdt_func)(sunrealtype t_curr, N_Vector manynvec_y, N_Vector manynvec_ydot, void *ctx);
  int (*dfdt_ssprk_func)(sunrealtype t_curr, N_Vector manynvec_y, N_Vector manynvec_ydot, void *ctx);
  int (*dfdt_sts_func)(sunrealtype t_curr, N_Vector manynvec_y, N_Vector manynvec_ydot, void *ctx);

  int (*gk_dom_eig_func)(sunrealtype t_curr, N_Vector manynvec_y, N_Vector manynvec_ydot, sunrealtype* lambdaR,
    sunrealtype* lambdaI, void *ctx, N_Vector temp1, N_Vector temp2, N_Vector temp3);
  int (*gk_dom_eig_sts_func)(sunrealtype t_curr, N_Vector manynvec_y, N_Vector manynvec_ydot, sunrealtype* lambdaR,
    sunrealtype* lambdaI, void *ctx, N_Vector temp1, N_Vector temp2, N_Vector temp3);

  int (*cfl_stable_dt_func)(N_Vector nvec_y, sunrealtype t_curr, sunrealtype *dt_out, void *ctx);
  int (*cfl_stable_dt_ssprk_func)(N_Vector nvec_y, sunrealtype t_curr, sunrealtype *dt_out, void *ctx);
  int (*cfl_stable_dt_sts_func)(N_Vector nvec_y, sunrealtype t_curr, sunrealtype *dt_out, void *ctx);

  int (*pre_process_step_ssprk_func)(sunrealtype t_curr, N_Vector manynvec_y, long int step, int attempt, void* ctx);
  int (*post_process_step_ssprk_func)(sunrealtype t_curr, N_Vector manynvec_y, long int step, void* ctx);
  int (*pre_process_rk_stage_ssprk_func)(sunrealtype t_curr, N_Vector manynvec_y, void* ctx);
  int (*post_process_rk_stage_ssprk_func)(sunrealtype t_curr, N_Vector manynvec_y, void* ctx);

  int (*pre_process_step_sts_func)(sunrealtype t_curr, N_Vector manynvec_y, long int step, int attempt, void* ctx);
  int (*post_process_step_sts_func)(sunrealtype t_curr, N_Vector manynvec_y, long int step, void* ctx);
  int (*pre_process_rk_stage_sts_func)(sunrealtype t_curr, N_Vector manynvec_y, void* ctx);
  int (*post_process_rk_stage_sts_func)(sunrealtype t_curr, N_Vector manynvec_y, void* ctx);

  int (*pre_process_rk_stage_extsts_func)(sunrealtype t_curr, N_Vector manynvec_y, void* ctx);
  int (*post_process_stage_extsts_func)(sunrealtype t_curr, N_Vector manynvec_y, void* ctx);

  int (*pre_process_step_opsplit_func)(sunrealtype t_curr, N_Vector manynvec_y, long int step, int attempt, void* ctx);
  int (*post_process_step_opsplit_func)(sunrealtype t_curr, N_Vector manynvec_y, long int step, void* ctx);
};

/**
 * Check function return value.
 *   opt == 0 means function allocates memory so check if returned NULL pointer.
 *   opt == 1 means function returns a flag so check if flag >= 0.
 */
static void
sundials_check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  if (opt == 0) {
    if (flagvalue == NULL) {
      // Check if function returned NULL pointer - no memory allocated.
      fprintf(stderr, "\nError: %s() failed - returned NULL pointer\n", funcname);
      assert(false);
    }
  }
  else if (opt == 1) {
    // Check if flag != 0.
    errflag = (int*)flagvalue;
    if (*errflag != 0) {
      fprintf(stderr, "\nError %s() failed with flag = %d\n", funcname, *errflag);
      assert(false);
    }
  }
}

/**
 * Return the gkyl_array wrapped by an Nvector.
 *
 * @param nvin Input Nvector.
 * @return A pointer to the gkyl_array wrapped by this Nvector.
 */
struct gkyl_array* snvec_get_array(N_Vector nvin);

/**
 * Return the gkyl_array wrapped by an Nvector, which
 * is itself wrapped by a ManyNvector.
 *
 * @param manynvin Input ManyNvector.
 * @param nvidx Index of desired Nvector.
 * @return A pointer to the gkyl_array wrapped by this Nvector.
 */
struct gkyl_array* smanynvec_get_array(N_Vector manynvin, int nvidx);
