#ifdef GKYL_HAVE_SUNDIALS

//
// Gyrokinetic-specific SUNDIALS methods.
//

#include <gkyl_time_integ_gyrokinetic.h>
#include <gkyl_sundials_priv.h>
#include <gkyl_sundials_gyrokinetic.h>

static void
unpack_manynvec_gyrokinetic_distf_in(struct gkyl_gyrokinetic_fdot_args *fdot_args, N_Vector manynvec)
{
  for (int i=0; i<fdot_args->num_species; ++i) {
    fdot_args->fin[i] = smanynvec_get_array(manynvec, fdot_args->offset_distf_charged[i]);
  }
  for (int i=0; i<fdot_args->num_neut_species; ++i) {
    fdot_args->fin_neut[i] = smanynvec_get_array(manynvec, fdot_args->offset_distf_neut[i]);
  }
}

static void
unpack_manynvec_gyrokinetic_bflux_in(struct gkyl_gyrokinetic_fdot_args *fdot_args, N_Vector manynvec)
{
  for (int i=0; i<fdot_args->num_species; ++i) {
    struct gkyl_array **bflux_s = fdot_args->bflux_in[i];
    for (int j=0; j<fdot_args->num_arr_bflux_charged[i]; ++j) {
      int off = fdot_args->offset_bflux_charged[i];
      bflux_s[j] = smanynvec_get_array(manynvec, off+j);
    }
  }
  for (int i=0; i<fdot_args->num_neut_species; ++i) {
    struct gkyl_array **bflux_s = fdot_args->bflux_in_neut[i];
    for (int j=0; j<fdot_args->num_arr_bflux_neut[i]; ++j) {
      int off = fdot_args->offset_bflux_neut[i];
      bflux_s[j] = smanynvec_get_array(manynvec, off+j);
    }
  }
}

static void
unpack_manynvec_gyrokinetic_charged_in(struct gkyl_gyrokinetic_fdot_args *fdot_args, N_Vector manynvec)
{
  for (int i=0; i<fdot_args->num_species; ++i) {
    fdot_args->fin[i] = smanynvec_get_array(manynvec, fdot_args->offset_distf_charged[i]);

    struct gkyl_array **bflux_s = fdot_args->bflux_in[i];
    for (int j=0; j<fdot_args->num_arr_bflux_charged[i]; ++j) {
      int off = fdot_args->offset_bflux_charged[i];
      bflux_s[j] = smanynvec_get_array(manynvec, off+j);
    }
  }
}

static void
unpack_manynvec_gyrokinetic_neut_in(struct gkyl_gyrokinetic_fdot_args *fdot_args, N_Vector manynvec)
{
  for (int i=0; i<fdot_args->num_neut_species; ++i) {
    fdot_args->fin_neut[i] = smanynvec_get_array(manynvec, fdot_args->offset_distf_neut[i]);

    struct gkyl_array **bflux_s = fdot_args->bflux_in_neut[i];
    for (int j=0; j<fdot_args->num_arr_bflux_neut[i]; ++j) {
      int off = fdot_args->offset_bflux_neut[i];
      bflux_s[j] = smanynvec_get_array(manynvec, off+j);
    }
  }
}

static void
unpack_manynvec_gyrokinetic_field_in(struct gkyl_gyrokinetic_fdot_args *fdot_args, N_Vector manynvec)
{
  for (int i=0; i<fdot_args->num_fields; ++i)
    fdot_args->fieldin[i] = smanynvec_get_array(manynvec, fdot_args->offset_field[i]);
}

static void
unpack_manynvec_gyrokinetic_charged_out(struct gkyl_gyrokinetic_fdot_args *fdot_args, N_Vector manynvec)
{
  for (int i=0; i<fdot_args->num_species; ++i) {
    fdot_args->fout[i] = smanynvec_get_array(manynvec, fdot_args->offset_distf_charged[i]);

    struct gkyl_array **bflux_s = fdot_args->bflux_out[i];
    for (int j=0; j<fdot_args->num_arr_bflux_charged[i]; ++j) {
      int off = fdot_args->offset_bflux_charged[i];
      bflux_s[j] = smanynvec_get_array(manynvec, off+j);
    }
  }
}

static void
unpack_manynvec_gyrokinetic_neut_out(struct gkyl_gyrokinetic_fdot_args *fdot_args, N_Vector manynvec)
{
  for (int i=0; i<fdot_args->num_neut_species; ++i) {
    fdot_args->fout_neut[i] = smanynvec_get_array(manynvec, fdot_args->offset_distf_neut[i]);

    struct gkyl_array **bflux_s = fdot_args->bflux_out_neut[i];
    for (int j=0; j<fdot_args->num_arr_bflux_neut[i]; ++j) {
      int off = fdot_args->offset_bflux_neut[i];
      bflux_s[j] = smanynvec_get_array(manynvec, off+j);
    }
  }
}

static void
unpack_manynvec_gyrokinetic_field_out(struct gkyl_gyrokinetic_fdot_args *fdot_args, N_Vector manynvec)
{
  for (int i=0; i<fdot_args->num_fields; ++i)
    fdot_args->fieldout[i] = smanynvec_get_array(manynvec, fdot_args->offset_field[i]);
}

static void
unpack_manynvec_gyrokinetic_charged(struct gkyl_gyrokinetic_fdot_args *fdot_args,
  N_Vector manynvec_in, N_Vector manynvec_out)
{
  for (int i=0; i<fdot_args->num_species; ++i) {
    fdot_args->fin[i] = smanynvec_get_array(manynvec_in, fdot_args->offset_distf_charged[i]);
    fdot_args->fout[i] = smanynvec_get_array(manynvec_out, fdot_args->offset_distf_charged[i]);

    struct gkyl_array **bflux_in_s = fdot_args->bflux_in[i];
    struct gkyl_array **bflux_out_s = fdot_args->bflux_out[i];
    for (int j=0; j<fdot_args->num_arr_bflux_charged[i]; ++j) {
      int off = fdot_args->offset_bflux_charged[i];
      bflux_in_s[j] = smanynvec_get_array(manynvec_in, off+j);
      bflux_out_s[j] = smanynvec_get_array(manynvec_out, off+j);
    }
  }
}

static void
unpack_manynvec_gyrokinetic_neut(struct gkyl_gyrokinetic_fdot_args *fdot_args,
  N_Vector manynvec_in, N_Vector manynvec_out)
{
  for (int i=0; i<fdot_args->num_neut_species; ++i) {
    fdot_args->fin_neut[i] = smanynvec_get_array(manynvec_in, fdot_args->offset_distf_neut[i]);
    fdot_args->fout_neut[i] = smanynvec_get_array(manynvec_out, fdot_args->offset_distf_neut[i]);

    struct gkyl_array **bflux_in_s = fdot_args->bflux_in_neut[i];
    struct gkyl_array **bflux_out_s = fdot_args->bflux_out_neut[i];
    for (int j=0; j<fdot_args->num_arr_bflux_neut[i]; ++j) {
      int off = fdot_args->offset_bflux_neut[i];
      bflux_in_s[j] = smanynvec_get_array(manynvec_in, off+j);
      bflux_out_s[j] = smanynvec_get_array(manynvec_out, off+j);
    }
  }
}

static void
unpack_manynvec_gyrokinetic_field(struct gkyl_gyrokinetic_fdot_args *fdot_args,
  N_Vector manynvec_in, N_Vector manynvec_out)
{
  for (int i=0; i<fdot_args->num_fields; ++i) {
    fdot_args->fieldin[i] = smanynvec_get_array(manynvec_in, fdot_args->offset_field[i]);
    fdot_args->fieldout[i] = smanynvec_get_array(manynvec_out, fdot_args->offset_field[i]);
  }
}

static inline void
unpack_manynvec_gyrokinetic_in(struct gkyl_gyrokinetic_fdot_args *fdot_args,
  N_Vector manynvec_in)
{
  unpack_manynvec_gyrokinetic_charged_in(fdot_args, manynvec_in);
  unpack_manynvec_gyrokinetic_neut_in(fdot_args, manynvec_in);
  unpack_manynvec_gyrokinetic_field_in(fdot_args, manynvec_in);
}

static inline void
unpack_manynvec_gyrokinetic_out(struct gkyl_gyrokinetic_fdot_args *fdot_args,
  N_Vector manynvec_out)
{
  unpack_manynvec_gyrokinetic_charged_out(fdot_args, manynvec_out);
  unpack_manynvec_gyrokinetic_neut_out(fdot_args, manynvec_out);
  unpack_manynvec_gyrokinetic_field_out(fdot_args, manynvec_out);
}

static inline void
unpack_manynvec_gyrokinetic(struct gkyl_gyrokinetic_fdot_args *fdot_args,
  N_Vector manynvec_in, N_Vector manynvec_out)
{
  unpack_manynvec_gyrokinetic_charged(fdot_args, manynvec_in, manynvec_out);
  unpack_manynvec_gyrokinetic_neut(fdot_args, manynvec_in, manynvec_out);
  unpack_manynvec_gyrokinetic_field(fdot_args, manynvec_in, manynvec_out);
}

/**
 * Error weight function for error norm of y_{n-1}.
 *
 * @param x ManyNvector y_{n-1} whose norm appears in the weight.
 * @param w ManyNvector weight to be computed.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
snvec_efun_cell_norm_gyrokinetic(N_Vector manyx, N_Vector manyw, void *ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;

  int ns_charged = fdot_args->num_species;
  int ns_neut = fdot_args->num_neut_species;

  int flag = 0;
  for (int i=0; i<ns_charged; i++) {
    N_Vector x = N_VGetSubvector_ManyVector(manyx, fdot_args->offset_distf_charged[i]);
    N_Vector w = N_VGetSubvector_ManyVector(manyw, fdot_args->offset_distf_charged[i]);
    struct gkyl_array *x_arr = snvec_get_array(x);
    struct gkyl_array *w_arr = snvec_get_array(w);
    struct gkyl_range *local_range = NV_CONTENT_GKZ(w)->local_range;

    // Call the Gkeyll function that computes the error weight.
    flag = flag || app_ctx->error_wgt_func(app_ctx->app_ptr, x_arr, w_arr, local_range);
  }

  for (int i=0; i<ns_neut; i++) {
    N_Vector x = N_VGetSubvector_ManyVector(manyx, fdot_args->offset_distf_neut[i]);
    N_Vector w = N_VGetSubvector_ManyVector(manyw, fdot_args->offset_distf_neut[i]);
    struct gkyl_array *x_arr = snvec_get_array(x);
    struct gkyl_array *w_arr = snvec_get_array(w);
    struct gkyl_range *local_range = NV_CONTENT_GKZ(w)->local_range;

    // Call the Gkeyll function that computes the error weight.
    flag = flag || app_ctx->error_wgt_func(app_ctx->app_ptr, x_arr, w_arr, local_range);
  }

  return flag;
}

/**
 * Compute the RHS function df/dt.
 *
 * @param t_curr Current simulation time.
 * @param mnvec_y State vectors f_s.
 * @param mnvec_ydot Time rate of change df_s/dt of each state vector.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
dfdt_gyrokinetic(sunrealtype t_curr, N_Vector manynvec_y, N_Vector manynvec_ydot, void *ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;

  // Distribute state vector as Gkeyll expects.
  int stage_idx, num_stages;
  int flag = ARKodeGetStageIndex(app_ctx->arkode_mem_outer, &stage_idx, &num_stages);
  if (stage_idx == 0) {
    // Treat this stage different because SUNDIALS passes yin instead of ycur
    // in stage=0 to avoid a yin->ycur copy at the beginning of every step.
    N_Vector manynvec_ycur;
    flag = ARKodeGetCurrentState(app_ctx->arkode_mem_outer, &manynvec_ycur);
    unpack_manynvec_gyrokinetic(fdot_args, manynvec_ycur, manynvec_ydot);
  }
  else {
    unpack_manynvec_gyrokinetic(fdot_args, manynvec_y, manynvec_ydot);
  }

  // Call the Gkeyll function that computes df/dt. Store local CFL constrained
  // dt (may not be used, depends on stepping method used).
  app_ctx->dt_local = app_ctx->dfdt_func(app_ctx->app_ptr, t_curr, fdot_args);

  return 0; // Success.
}

/**
 * Compute the RHS function df/dt due to operators stepped with SSP-RK.
 *
 * @param t_curr Current simulation time.
 * @param mnvec_y State vectors f_s.
 * @param mnvec_ydot Time rate of change df_s/dt of each state vector.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
dfdt_ssprk_gyrokinetic(sunrealtype t_curr, N_Vector manynvec_y, N_Vector manynvec_ydot, void *ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;

  // Distribute state vector as Gkeyll expects.
  int stage_idx, num_stages;
  int flag = ARKodeGetStageIndex(app_ctx->arkode_mem_ssprk, &stage_idx, &num_stages);
  if (stage_idx == 0) {
    // Treat this stage different because SUNDIALS passes yin instead of ycur
    // in stage=0 to avoid a yin->ycur copy at the beginning of every step.
    N_Vector manynvec_ycur;
    flag = ARKodeGetCurrentState(app_ctx->arkode_mem_ssprk, &manynvec_ycur);
    unpack_manynvec_gyrokinetic(fdot_args, manynvec_ycur, manynvec_ydot);
  }
  else {
    unpack_manynvec_gyrokinetic(fdot_args, manynvec_y, manynvec_ydot);
  }

  // Call the Gkeyll function that computes df/dt. Store local CFL constrained
  // dt (may not be used, depends on stepping method used).
  app_ctx->dt_local_ssprk = app_ctx->dfdt_ssprk_func(app_ctx->app_ptr, t_curr, fdot_args);

  return 0; // Success.
}

/**
 * Compute the RHS function df/dt due to operators stepped with STS.
 *
 * @param t_curr Current simulation time.
 * @param mnvec_y State vectors f_s.
 * @param mnvec_ydot Time rate of change df_s/dt of each state vector.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
dfdt_sts_gyrokinetic(sunrealtype t_curr, N_Vector manynvec_y, N_Vector manynvec_ydot, void *ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;

  // Distribute state vector as Gkeyll expects.
  unpack_manynvec_gyrokinetic(fdot_args, manynvec_y, manynvec_ydot);

  // Call the Gkeyll function that computes df/dt. Store local CFL constrained
  // dt (may not be used, depends on stepping method used).
  app_ctx->dt_local_sts = app_ctx->dfdt_sts_func(app_ctx->app_ptr, t_curr, fdot_args);

  return 0; // Success.
}

/**
 * Pre-process a SSP-RK step.
 *
 * @param t_curr Current simulation time.
 * @param manynvec_y State vector.
 * @param step Step number.
 * @param attempt Step attempt number.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
pre_process_step_gyrokinetic_ssprk(sunrealtype t_curr, N_Vector manynvec_y, long int step, int attempt, void* ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  double dt;
  int flag = ARKodeGetCurrentStep(app_ctx->arkode_mem_ssprk, &dt);

  // Distribute state vector as Gkeyll expects.
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;
  unpack_manynvec_gyrokinetic_in(fdot_args, manynvec_y);

  app_ctx->pre_process_step_ssprk_func(app_ctx->app_ptr, t_curr, dt, fdot_args);

  return 0;
};

/**
 * Pre-process a SSP-RK stage.
 *
 * @param t_curr Current simulation time.
 * @param manynvec_y State vector.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
pre_process_rk_stage_gyrokinetic_ssprk(sunrealtype t_curr, N_Vector manynvec_y, void* ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  int stage_idx, num_stages;
  double dt;
  int flag = 0;
  flag = ARKodeGetStageIndex(app_ctx->arkode_mem_ssprk, &stage_idx, &num_stages);
  flag = ARKodeGetCurrentStep(app_ctx->arkode_mem_ssprk, &dt);

  // We only pre-process the first stage, which needs to zero out the boundary fluxes the step starts with.
  // Maybe the best thing to do is pass the current state of f, but the last step of boundary fluxes.
  N_Vector manynvec_yin;
  flag = ARKodeGetLastState(app_ctx->arkode_mem_ssprk, &manynvec_yin);

  // Distribute state vector as Gkeyll expects.
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;
  unpack_manynvec_gyrokinetic_in(fdot_args, manynvec_yin);

  app_ctx->pre_process_rk_stage_ssprk_func(app_ctx->app_ptr, t_curr, dt, fdot_args, stage_idx, num_stages);
  return 0;
};

/**
 * Post-process a SSP-RK stage.
 *
 * @param t_curr Current simulation time.
 * @param manynvec_ycur State vector.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
post_process_rk_stage_gyrokinetic_ssprk(sunrealtype t_curr, N_Vector manynvec_ycur, void* ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  int stage_idx, num_stages;
  double dt;
  int flag = 0;
  flag = ARKodeGetStageIndex(app_ctx->arkode_mem_ssprk, &stage_idx, &num_stages);
  flag = ARKodeGetCurrentStep(app_ctx->arkode_mem_ssprk, &dt);
  if (stage_idx > 0)
    stage_idx -= 1;

  // Distribute state vector as Gkeyll expects.
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;
  unpack_manynvec_gyrokinetic_out(fdot_args, manynvec_ycur);
  if (stage_idx == 0) {
    N_Vector manynvec_yn;
    flag = ARKodeGetLastState(app_ctx->arkode_mem_ssprk, &manynvec_yn);
    unpack_manynvec_gyrokinetic_distf_in(fdot_args, manynvec_yn);
    unpack_manynvec_gyrokinetic_field_in(fdot_args, manynvec_yn);

    unpack_manynvec_gyrokinetic_bflux_in(fdot_args, manynvec_ycur);
  }

  app_ctx->post_process_rk_stage_ssprk_func(app_ctx->app_ptr, t_curr, dt, fdot_args, stage_idx, num_stages);

  return 0;
};

/**
 * Post-process a SSP-RK step.
 *
 * @param t_curr Current simulation time.
 * @param manynvec_y State vector.
 * @param step Step number.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
post_process_step_gyrokinetic_ssprk(sunrealtype t_curr, N_Vector manynvec_y, long int step, void* ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  double dt;
  int flag = ARKodeGetCurrentStep(app_ctx->arkode_mem_ssprk, &dt);

  // Distribute state vector as Gkeyll expects.
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;
  unpack_manynvec_gyrokinetic_out(fdot_args, manynvec_y);

  app_ctx->post_process_step_ssprk_func(app_ctx->app_ptr, t_curr, dt, fdot_args);

  return 0;
};

/**
 * Pre-process a STS step.
 *
 * @param t_curr Current simulation time.
 * @param manynvec_y State vector.
 * @param step Step number.
 * @param attempt Step attempt number.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
pre_process_step_gyrokinetic_sts(sunrealtype t_curr, N_Vector manynvec_y, long int step, int attempt, void* ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  double dt;
  int flag = ARKodeGetCurrentStep(app_ctx->arkode_mem_sts, &dt);

  // Distribute state vector as Gkeyll expects.
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;
  unpack_manynvec_gyrokinetic_in(fdot_args, manynvec_y);

  app_ctx->pre_process_step_sts_func(app_ctx->app_ptr, t_curr, dt, fdot_args);

  return 0;
};

/**
 * Pre-process a STS stage.
 *
 * @param t_curr Current simulation time.
 * @param manynvec_y State vector.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
pre_process_rk_stage_gyrokinetic_sts(sunrealtype t_curr, N_Vector manynvec_y, void* ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  int stage_idx, num_stages;
  double dt;
  int flag = 0;
  flag = ARKodeGetStageIndex(app_ctx->arkode_mem_sts, &stage_idx, &num_stages);
  flag = ARKodeGetCurrentStep(app_ctx->arkode_mem_sts, &dt);

  // We only pre-process the first stage, which needs to zero out the boundary fluxes the step starts with.
  // Maybe the best thing to do is pass the current state of f, but the last step of boundary fluxes.
  N_Vector manynvec_yin;
  flag = ARKodeGetLastState(app_ctx->arkode_mem_sts, &manynvec_yin);

  // Distribute state vector as Gkeyll expects.
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;
  unpack_manynvec_gyrokinetic_in(fdot_args, manynvec_yin);

  app_ctx->pre_process_rk_stage_sts_func(app_ctx->app_ptr, t_curr, dt, fdot_args, stage_idx, num_stages);
  return 0;
};

/**
 * Post-process a STS stage.
 *
 * @param t_curr Current simulation time.
 * @param manynvec_y State vector.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
post_process_rk_stage_gyrokinetic_sts(sunrealtype t_curr, N_Vector manynvec_y, void* ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  int stage_idx, num_stages;
  double dt;
  int flag = 0;
  flag = ARKodeGetStageIndex(app_ctx->arkode_mem_sts, &stage_idx, &num_stages);
  flag = ARKodeGetCurrentStep(app_ctx->arkode_mem_sts, &dt);
  if (stage_idx > 0)
    stage_idx -= 1;

  // Distribute state vector as Gkeyll expects.
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;
  unpack_manynvec_gyrokinetic_out(fdot_args, manynvec_y);
  if (stage_idx == 0) {
    N_Vector manynvec_yin;
    flag = ARKodeGetLastState(app_ctx->arkode_mem_sts, &manynvec_yin);
    unpack_manynvec_gyrokinetic_distf_in(fdot_args, manynvec_yin);
    unpack_manynvec_gyrokinetic_field_in(fdot_args, manynvec_yin);

    flag = ARKodeGetCurrentState(app_ctx->arkode_mem_sts, &manynvec_yin);
    unpack_manynvec_gyrokinetic_bflux_in(fdot_args, manynvec_yin);
  }

  app_ctx->post_process_rk_stage_sts_func(app_ctx->app_ptr, t_curr, dt, fdot_args, stage_idx, num_stages);

  return 0;
};

/**
 * Post-process a STS step.
 *
 * @param t_curr Current simulation time.
 * @param manynvec_y State vector.
 * @param step Step number.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
post_process_step_gyrokinetic_sts(sunrealtype t_curr, N_Vector manynvec_y, long int step, void* ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  double dt;
  int flag = ARKodeGetCurrentStep(app_ctx->arkode_mem_sts, &dt);

  // Distribute state vector as Gkeyll expects.
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;
  unpack_manynvec_gyrokinetic_out(fdot_args, manynvec_y);

  app_ctx->post_process_step_sts_func(app_ctx->app_ptr, t_curr, dt, fdot_args);

  return 0;
};

/**
 * Pre-process the evaluation of the SSP-RK terms at an ExtSTS outer stage.
 *
 * Unlike SSP-RK methods, whose stages are built from the state at the start
 * of the step (yn), ExtSTS stages are built from the current state (ycur), so
 * the operations of the first stage (e.g. clearing boundary fluxes) are applied
 * to the input vector. They are also applied to yn, because if the step is
 * retried ycur is restored from yn without calling this function again (the
 * first-stage RHS is reused).
 *
 * @param t_curr Current simulation time.
 * @param manynvec_y State vector (ycur).
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
pre_process_rk_stage_gyrokinetic_extsts(sunrealtype t_curr, N_Vector manynvec_y, void* ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  int stage_idx, num_stages;
  double dt;
  int flag = 0;
  flag = ARKodeGetStageIndex(app_ctx->arkode_mem_outer, &stage_idx, &num_stages);
  flag = ARKodeGetCurrentStep(app_ctx->arkode_mem_outer, &dt);

  // Distribute state vector as Gkeyll expects.
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;
  unpack_manynvec_gyrokinetic_in(fdot_args, manynvec_y);

  app_ctx->pre_process_rk_stage_ssprk_func(app_ctx->app_ptr, t_curr, dt, fdot_args, stage_idx, num_stages);

  if (stage_idx == 0) {
    N_Vector manynvec_yn;
    flag = ARKodeGetLastState(app_ctx->arkode_mem_outer, &manynvec_yn);
    if (manynvec_yn != manynvec_y) {
      unpack_manynvec_gyrokinetic_in(fdot_args, manynvec_yn);
      app_ctx->pre_process_rk_stage_ssprk_func(app_ctx->app_ptr, t_curr, dt, fdot_args, stage_idx, num_stages);
    }
  }
  return 0;
};

/**
 * Post-process an ExtSTS outer stage.
 *
 * @param t_curr Current simulation time.
 * @param manynvec_ycur Stage solution.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
post_process_stage_gyrokinetic_extsts(sunrealtype t_curr, N_Vector manynvec_ycur, void* ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  int stage_idx, num_stages;
  double dt;
  int flag = 0;
  flag = ARKodeGetStageIndex(app_ctx->arkode_mem_outer, &stage_idx, &num_stages);
  flag = ARKodeGetCurrentStep(app_ctx->arkode_mem_outer, &dt);

  // Distribute state vector as Gkeyll expects.
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;
  unpack_manynvec_gyrokinetic_out(fdot_args, manynvec_ycur);
  if (stage_idx == 1) {
    // The first stage after the start of the step also needs the state at the
    // start of the step.
    N_Vector manynvec_yn;
    flag = ARKodeGetLastState(app_ctx->arkode_mem_outer, &manynvec_yn);
    unpack_manynvec_gyrokinetic_distf_in(fdot_args, manynvec_yn);
    unpack_manynvec_gyrokinetic_field_in(fdot_args, manynvec_yn);
  }

  app_ctx->post_process_stage_extsts_func(app_ctx->app_ptr, t_curr, dt, fdot_args, stage_idx, num_stages);

  return 0;
};

/**
 * Pre-process a operator split step.
 *
 * @param t_curr Current simulation time.
 * @param manynvec_y State vector.
 * @param step Step number.
 * @param attempt Step attempt number.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
pre_process_step_gyrokinetic_opsplit(sunrealtype t_curr, N_Vector manynvec_y, long int step, int attempt, void* ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  double dt;
  int flag = ARKodeGetCurrentStep(app_ctx->arkode_mem_outer, &dt);

  // Distribute state vector as Gkeyll expects.
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;
  unpack_manynvec_gyrokinetic_in(fdot_args, manynvec_y);

  app_ctx->pre_process_step_opsplit_func(app_ctx->app_ptr, t_curr, dt, fdot_args);

  return 0;
};

/**
 * Post-process a operator split step.
 *
 * @param t_curr Current simulation time.
 * @param manynvec_y State vector.
 * @param step Step number.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
post_process_step_gyrokinetic_opsplit(sunrealtype t_curr, N_Vector manynvec_y, long int step, void* ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  double dt;
  int flag = ARKodeGetCurrentStep(app_ctx->arkode_mem_outer, &dt);

  // Distribute state vector as Gkeyll expects.
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;
  unpack_manynvec_gyrokinetic_out(fdot_args, manynvec_y);

  app_ctx->post_process_step_opsplit_func(app_ctx->app_ptr, t_curr, dt, fdot_args);

  return 0;
};

/**
 * Gkeyll's estimate of the dominant eigenvalue of the df/dt operator.
 *
 * @param t_curr Current simulation time.
 * @param manynvec_y State vectors f_s.
 * @param manynvec_ydot Time rate of change df_s/dt of each state vector.
 * @param lambdaR Real part of the dominant eigenvalue.
 * @param lambdaI Imaginary part of the dominant eigenvalue.
 * @param ctx App-specific context.
 * @param temp1,temp2,temp3 Buffers that may be used locally if needed.
 * @return Sucess (=0) flag.
 */
static int
gk_dom_eig_gyrokinetic(sunrealtype t_curr, N_Vector manynvec_y, N_Vector manynvec_ydot,
  sunrealtype* lambdaR, sunrealtype* lambdaI, void *ctx,
  N_Vector temp1, N_Vector temp2, N_Vector temp3)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  // Distribute state vector as Gkeyll expects.
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;
  unpack_manynvec_gyrokinetic(fdot_args, manynvec_y, manynvec_ydot);

  // Call the Gkeyll function that computes df/dt due to the operator stepped
  // with STS and compute it slocal CFL constrained dt.
  double dt_local = app_ctx->dfdt_func(app_ctx->app_ptr, t_curr, fdot_args);

  double dt_global = app_ctx->reduce_dt_func(app_ctx->app_ptr, t_curr, dt_local);

  *lambdaR = -1.0/dt_global;
  *lambdaI = SUN_RCONST(0.0);

  return 0; // Success.
}

/**
 * Gkeyll's estimate of the dominant eigenvalue of the df/dt terms stepped with STS.
 *
 * @param t_curr Current simulation time.
 * @param manynvec_y State vectors f_s.
 * @param manynvec_ydot Time rate of change df_s/dt of each state vector.
 * @param lambdaR Real part of the dominant eigenvalue.
 * @param lambdaI Imaginary part of the dominant eigenvalue.
 * @param ctx App-specific context.
 * @param temp1,temp2,temp3 Buffers that may be used locally if needed.
 * @return Sucess (=0) flag.
 */
static int
gk_dom_eig_sts_gyrokinetic(sunrealtype t_curr, N_Vector manynvec_y, N_Vector manynvec_ydot,
  sunrealtype* lambdaR, sunrealtype* lambdaI, void *ctx,
  N_Vector temp1, N_Vector temp2, N_Vector temp3)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  // Distribute state vector as Gkeyll expects.
  struct gkyl_gyrokinetic_fdot_args *fdot_args = app_ctx->fdot_args_ptr;
  unpack_manynvec_gyrokinetic(fdot_args, manynvec_y, manynvec_ydot);

  // Call the Gkeyll function that computes df/dt due to the operator stepped
  // with STS and compute it slocal CFL constrained dt.
  double dt_local = app_ctx->dfdt_sts_func(app_ctx->app_ptr, t_curr, fdot_args);

  double dt_global = app_ctx->reduce_dt_func(app_ctx->app_ptr, t_curr, dt_local);

  *lambdaR = -1.0/dt_global;
  *lambdaI = SUN_RCONST(0.0);

  return 0; // Success.
}

/**
 * Return the time step dt required for stability according to Gkeyll's CFL
 * constraint.
 *
 * This takes the local dt computed when df/dt was called, and performs
 * a min reduction.
 *
 * @param nvec_y State vector.
 * @param t_curr Current simulation time.
 * @param dt_out CFL-stable dt.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
cfl_stable_dt_gyrokinetic(N_Vector nvec_y, sunrealtype t_curr, sunrealtype *dt_out, void *ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  double dt_local = app_ctx->dt_local;

  *dt_out = app_ctx->reduce_dt_func(app_ctx->app_ptr, t_curr, dt_local);

  return 0; // Success.
}

/**
 * Return the time step dt required for stability of terms stepped with SSP-RK
 * according to Gkeyll's CFL constraint.
 *
 * This takes the local dt computed when df/dt was called, and performs
 * a min reduction.
 *
 * @param nvec_y State vector.
 * @param t_curr Current simulation time.
 * @param dt_out CFL-stable dt.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
cfl_stable_dt_ssprk_gyrokinetic(N_Vector nvec_y, sunrealtype t_curr, sunrealtype *dt_out, void *ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  double dt_local = app_ctx->dt_local_ssprk;

  *dt_out = app_ctx->reduce_dt_func(app_ctx->app_ptr, t_curr, dt_local);

  return 0; // Success.
}

/**
 * Return the time step dt required for stability of terms stepped with STS
 * according to Gkeyll's CFL constraint.
 *
 * This takes the local dt computed when df/dt was called, and performs
 * a min reduction.
 *
 * @param nvec_y State vector.
 * @param t_curr Current simulation time.
 * @param dt_out CFL-stable dt.
 * @param ctx App-specific context.
 * @return Sucess (=0) flag.
 */
static int
cfl_stable_dt_sts_gyrokinetic(N_Vector nvec_y, sunrealtype t_curr, sunrealtype *dt_out, void *ctx)
{
  struct gkyl_sundials_app_ctx *app_ctx = ctx;

  double dt_local = app_ctx->dt_local_sts;

  *dt_out = app_ctx->reduce_dt_func(app_ctx->app_ptr, t_curr, dt_local);

  return 0; // Success.
}

void
gkyl_sundials_gyrokinetic_assign_methods(struct gkyl_sundials *gksun)
{
  gksun->dfdt_func = dfdt_gyrokinetic;
  gksun->dfdt_ssprk_func = dfdt_ssprk_gyrokinetic;
  gksun->dfdt_sts_func = dfdt_sts_gyrokinetic;

  gksun->gk_dom_eig_func = gk_dom_eig_gyrokinetic;
  gksun->gk_dom_eig_sts_func = gk_dom_eig_sts_gyrokinetic;

  gksun->cfl_stable_dt_func = cfl_stable_dt_gyrokinetic;
  gksun->cfl_stable_dt_ssprk_func = cfl_stable_dt_ssprk_gyrokinetic;
  gksun->cfl_stable_dt_sts_func = cfl_stable_dt_sts_gyrokinetic;

  gksun->snvec_efun_cell_norm_func = snvec_efun_cell_norm_gyrokinetic;

  gksun->pre_process_step_ssprk_func = pre_process_step_gyrokinetic_ssprk;
  gksun->pre_process_rk_stage_ssprk_func = pre_process_rk_stage_gyrokinetic_ssprk;
  gksun->post_process_rk_stage_ssprk_func = post_process_rk_stage_gyrokinetic_ssprk;
  gksun->post_process_step_ssprk_func = post_process_step_gyrokinetic_ssprk;

  gksun->pre_process_step_sts_func = pre_process_step_gyrokinetic_sts;
  gksun->pre_process_rk_stage_sts_func = pre_process_rk_stage_gyrokinetic_sts;
  gksun->post_process_rk_stage_sts_func = post_process_rk_stage_gyrokinetic_sts;
  gksun->post_process_step_sts_func = post_process_step_gyrokinetic_sts;

  gksun->pre_process_rk_stage_extsts_func = pre_process_rk_stage_gyrokinetic_extsts;
  gksun->post_process_stage_extsts_func = post_process_stage_gyrokinetic_extsts;

  gksun->pre_process_step_opsplit_func = pre_process_step_gyrokinetic_opsplit;
  gksun->post_process_step_opsplit_func = post_process_step_gyrokinetic_opsplit;
}

#endif
