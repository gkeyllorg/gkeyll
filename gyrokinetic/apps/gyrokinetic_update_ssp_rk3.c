#include <gkyl_gyrokinetic_priv.h>

static void
gyrokinetic_forward_euler(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args, struct gkyl_update_status *st)
{
  struct timespec wst_fe = gkyl_wall_clock();
  // Take a forward Euler step with the suggested time-step dt. This may
  // not be the actual time-step taken. However, the function will never
  // take a time-step larger than dt even if it is allowed by
  // stability. The actual time-step and dt_suggested are returned in
  // the status object.
  app->stat.nfeuler += 1;

  // Compute the time rate of change of the distributions, df/dt.
  gyrokinetic_rhs(app, tcurr, dt, fdot_args, st);

  // Complete update of distribution functions.
  struct timespec wst = gkyl_wall_clock();
  const struct gkyl_array **fin = fdot_args->fin;
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array ***bflux_in = fdot_args->bflux_in;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;

  const struct gkyl_array **fin_neut = fdot_args->fin_neut;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***bflux_in_neut = fdot_args->bflux_in_neut;
  struct gkyl_array ***bflux_out_neut = fdot_args->bflux_out_neut;

  double dta = st->dt_actual;

  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_step_f(gks, fout[i], dta, fin[i]);
    gk_species_bflux_step_f(app, &gks->bflux, bflux_out[i], dta, bflux_in[i]);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_step_f(gkns, fout_neut[i], dta, fin_neut[i]);
    gk_neut_species_bflux_step_f(app, &gkns->bflux, bflux_out_neut[i], dta, bflux_in_neut[i]);
  }
  app->stat.fwd_euler_step_f_tm += gkyl_time_diff_now_sec(wst);

  app->stat.fwd_euler_tm += gkyl_time_diff_now_sec(wst_fe);
}

static void
gyrokinetic_bflux_scale(gkyl_gyrokinetic_app* app, struct gkyl_array ***bflux,
  struct gkyl_array ***bflux_neut, double fac)
{
  // Scale boundary fluxes by a factor.
  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_bflux_scale(app, &gks->bflux, bflux[i], fac);
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_bflux_scale(app, &gkns->bflux, bflux_neut[i], fac);
  }
}

void
gyrokinetic_pre_process_step_ssprk(gkyl_gyrokinetic_app* app, double tcurr,
  double dt, struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the beginning of a SSP-RK step.
  // NOTE: do not use dt here because its size is not yet determined.
  struct gkyl_array ***bflux_in = fdot_args->bflux_in;
  struct gkyl_array ***bflux_in_neut = fdot_args->bflux_in_neut;

  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;

  // Adapt sources.
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_source_adapt(app, gks, &gks->src, gks->lte.f_lte, bflux_in, tcurr);
  }
}

static void
gyrokinetic_pre_process_rk_stage_initial_ssprk(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the beginning of the first SSP-RK stage.
  // NOTE: do not use dt here because its size is not yet determined.
  struct gkyl_array ***bflux_in = fdot_args->bflux_in;
  struct gkyl_array ***bflux_in_neut = fdot_args->bflux_in_neut;

  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;

  // Clear input boundary fluxes so we compute fluxes through the boundary in one step.
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_bflux_clear(app, &gks->bflux, bflux_in[i], 0.0);
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_bflux_clear(app, &gkns->bflux, bflux_in_neut[i], 0.0);
  }
}

void
gyrokinetic_pre_process_rk_stage_ssprk(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args, int stage_idx, int num_stages)
{
  // Operations at the end of an SSP-RK stage.
  if (stage_idx == 0) {
    // First RK stage.
    gyrokinetic_pre_process_rk_stage_initial_ssprk(app, tcurr, dt, fdot_args);
  }
}

static void
gyrokinetic_post_process_rk_stage_initial_ssprk(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the end of the first SSP-RK stage.
  const struct gkyl_array **fin = fdot_args->fin;
  const struct gkyl_array **fieldin = fdot_args->fieldin;
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;
  struct gkyl_array ***bflux_out_neut = fdot_args->bflux_out_neut;
  struct gkyl_array **fieldout = fdot_args->fieldout;

  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;

  // Divide new bflux by 1/dt (since forward Euler multiplied by dt).
  struct timespec wst = gkyl_wall_clock();
  gyrokinetic_bflux_scale(app, bflux_out, bflux_out_neut, 1.0/dt);
  app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    // Compute moment of f_old to later compute moment of df/dt.
    gk_species_calc_int_mom_dt(app, gks, (struct gkyl_array *) fin[i],
      (struct gkyl_array **) fieldin, dt, gks->fdot_mom_old);
  }

  // Compute field energy divided by dt for energy balance diagnostics.
  gk_field_calc_energy_dt(app, app->field, dt, (struct gkyl_array **) fieldin, app->field->em_energy_red_old);

  // Compute the fields and apply BCs.
  gyrokinetic_calc_field(app, tcurr, fout, bflux_out, fieldout);
  gyrokinetic_apply_bc(app, tcurr, fieldout, fout, fout_neut);

  // Multiply new bflux by dt (so forward Euler multiplied can add dt*bflux).
  wst = gkyl_wall_clock();
  gyrokinetic_bflux_scale(app, bflux_out, bflux_out_neut, dt);
  app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);
}

static void
gyrokinetic_post_process_rk_stage_intermediate_ssprk(gkyl_gyrokinetic_app* app, double tcurr,
  double dt, struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the end of intermediate SSP-RK stages.
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;
  struct gkyl_array ***bflux_out_neut = fdot_args->bflux_out_neut;
  struct gkyl_array **fieldout = fdot_args->fieldout;

  // Divide new bflux by 1/dt (since forward Euler multiplied by dt).
  struct timespec wst = gkyl_wall_clock();
  gyrokinetic_bflux_scale(app, bflux_out, bflux_out_neut, 1.0/dt);
  app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

  // Update fields and apply BCs.
  gyrokinetic_calc_field(app, tcurr, fout, bflux_out, fieldout);
  gyrokinetic_apply_bc(app, tcurr, fieldout, fout, fout_neut);

  // Multiply new bflux by dt (so forward Euler multiplied can add dt*bflux).
  wst = gkyl_wall_clock();
  gyrokinetic_bflux_scale(app, bflux_out, bflux_out_neut, dt);
  app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);
}

void
gyrokinetic_post_process_rk_stage_ssprk(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args, int stage_idx, int num_stages)
{
  // Operations at the end of an SSP-RK stage.
  if (stage_idx == 0) {
    // First RK stage.
    gyrokinetic_post_process_rk_stage_initial_ssprk(app, tcurr, dt, fdot_args);
  }
  else if (stage_idx < num_stages-1) {
    // Intermediate RK stages.
    gyrokinetic_post_process_rk_stage_intermediate_ssprk(app, tcurr, dt, fdot_args);
  }
}

void
gyrokinetic_post_process_stage_extsts(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args, int stage_idx, int num_stages)
{
  // Operations at the end of an ExtSTS outer stage. Stage 0 is the state at
  // the start of the step and the last stage is handled at the end of the
  // step, so stages 1 to num_stages-2 are post-processed: the first of these
  // like the first SSP-RK stage, the others like intermediate SSP-RK stages.
  if (stage_idx == 1)
    gyrokinetic_post_process_rk_stage_initial_ssprk(app, tcurr, dt, fdot_args);
  else
    gyrokinetic_post_process_rk_stage_intermediate_ssprk(app, tcurr, dt, fdot_args);
}

void
gyrokinetic_post_process_step_ssprk(gkyl_gyrokinetic_app* app, double tcurr,
  double dt, struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the end of the SSP-RK step.
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;
  struct gkyl_array ***bflux_out_neut = fdot_args->bflux_out_neut;
  struct gkyl_array **fieldout = fdot_args->fieldout;

  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;

  // Divide new bflux by 1/dt (since forward Euler multiplied by dt).
  struct timespec wst = gkyl_wall_clock();
  gyrokinetic_bflux_scale(app, bflux_out, bflux_out_neut, 1.0/dt);
  app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

  // Compute volume- and time-integrated boundary fluxes.
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_bflux_calc_voltime_integrated_mom(app, gks, &gks->bflux, tcurr);
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_bflux_calc_voltime_integrated_mom(app, gkns, &gkns->bflux, tcurr);
  }

  // Apply positivity shift if requested.
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_positivity_apply(app, gks, &gks->positivity, gks->fnew, fout[i]);
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_positivity_apply(app, gkns, &gkns->positivity, gkns->fnew, fout_neut[i]);
  }

  // Enforce quasineutrality of the positivity shifts.
  gyrokinetic_post_positivity_quasineut(app, fout);

  // Compute the fields and apply BCs
  gyrokinetic_calc_field(app, tcurr, fout, bflux_out, fieldout);
  gyrokinetic_apply_bc(app, tcurr, fieldout, fout, fout_neut);

  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    // Compute moment of f_new to compute moment of df/dt.
    // Need to do it after the fields are updated.
    gk_species_calc_int_mom_dt(app, gks, fout[i], fieldout, dt, gks->fdot_mom_new);
  }

  // Scale species according to some criteria.
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_scaling_apply(app, gks, &gks->sca, gks->f, fieldout, bflux_out);
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_scaling_apply(app, gkns, &gkns->sca, gkns->f, fieldout, bflux_out);
  }

  // Compute field energy divided by dt for energy balance diagnostics.
  gk_field_calc_energy_dt(app, app->field, dt, fieldout, app->field->em_energy_red_new);
}

struct gkyl_update_status
gyrokinetic_update_ssp_rk3(gkyl_gyrokinetic_app* app, double dt0)
{
  // Take time-step using the RK3 method. Also sets the status object
  // which has the actual and suggested dts used. These can be different
  // from the actual time-step.

  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;
  int nfields = app->field->num_fields;

  struct gkyl_gyrokinetic_fdot_args *fdot_args = &app->fdot_args;
  const struct gkyl_array **fin = fdot_args->fin;
  const struct gkyl_array **fin_neut = fdot_args->fin_neut;
  const struct gkyl_array **fieldin = fdot_args->fieldin;
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array **fieldout = fdot_args->fieldout;
  struct gkyl_array ***bflux_in = fdot_args->bflux_in;
  struct gkyl_array ***bflux_in_neut = fdot_args->bflux_in_neut;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;
  struct gkyl_array ***bflux_out_neut = fdot_args->bflux_out_neut;

  struct gkyl_update_status st = { .success = true };
 
  enum {
    SSPRK_3_3_STAGE_0 = 0,
    SSPRK_3_3_STAGE_1,
    SSPRK_3_3_STAGE_2,
    SSPRK_3_3_STAGE_COMPLETE,
  } state = SSPRK_3_3_STAGE_0; // Time-stepper state.

  double tcurr = app->tcurr, dt = dt0;

  // Pre-process step.
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    fin[i] = gks->f;
    bflux_in[i] = gks->bflux.f;
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    fin_neut[i] = gkns->f;
    bflux_in_neut[i] = gkns->bflux.f;
  }
  for (int i=0; i<nfields; ++i) {
    struct gk_field *gkf = app->field;
    fieldin[i] = gkf->f[i];
  }
  gyrokinetic_pre_process_step_ssprk(app, tcurr, dt, fdot_args);

  while (state != SSPRK_3_3_STAGE_COMPLETE) {
    switch (state) {
      case SSPRK_3_3_STAGE_0:
        for (int i=0; i<ns_charged; ++i) {
          struct gk_species *gks = &app->species[i];
          fin[i] = gks->f;
          fout[i] = gks->f1;
          bflux_in[i] = gks->bflux.f;
          bflux_out[i] = gks->bflux.f1;
        }
        for (int i=0; i<ns_neut; ++i) {
          struct gk_neut_species *gkns = &app->neut_species[i];
          fin_neut[i] = gkns->f;
          fout_neut[i] = gkns->f1;
          bflux_in_neut[i] = gkns->bflux.f;
          bflux_out_neut[i] = gkns->bflux.f1;
        }
        for (int i=0; i<nfields; ++i) {
          struct gk_field *gkf = app->field;
          fieldin[i] = gkf->f[i];
          fieldout[i] = gkf->f1[i];
        }

        // Pre-process stage.
        gyrokinetic_pre_process_rk_stage_initial_ssprk(app, tcurr, dt, fdot_args);

        // Step solution.
        gyrokinetic_forward_euler(app, tcurr, dt, fdot_args, &st);
        dt = st.dt_actual;

        // Post process stage.
        gyrokinetic_post_process_rk_stage_initial_ssprk(app, tcurr, dt, fdot_args);

        state = SSPRK_3_3_STAGE_1;
        break;

      case SSPRK_3_3_STAGE_1:
        for (int i=0; i<ns_charged; ++i) {
          struct gk_species *gks = &app->species[i];
          fin[i] = gks->f1;
          fout[i] = gks->fnew;
          bflux_in[i] = gks->bflux.f1;
          bflux_out[i] = gks->bflux.fnew;
        }
        for (int i=0; i<ns_neut; ++i) {
          struct gk_neut_species *gkns = &app->neut_species[i];
          fin_neut[i] = gkns->f1;
          fout_neut[i] = gkns->fnew;
          bflux_in_neut[i] = gkns->bflux.f1;
          bflux_out_neut[i] = gkns->bflux.fnew;
        }
        for (int i=0; i<nfields; ++i) {
          struct gk_field *gkf = app->field;
          fieldin[i] = gkf->f1[i];
          fieldout[i] = gkf->fnew[i];
        }

        // Step solution.
        gyrokinetic_forward_euler(app, tcurr+dt, dt, fdot_args, &st);

        if (st.dt_actual < dt) {

          // Failed step. Collect stats.
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_2_fail += 1;

          dt = st.dt_actual;
          state = SSPRK_3_3_STAGE_0; // Restart from stage 0.

        } 
        else {
          struct timespec wst = gkyl_wall_clock();
          for (int i=0; i<ns_charged; ++i) {
            struct gk_species *gks = &app->species[i];
            gk_species_combine(gks, gks->f1, 3.0/4.0, gks->f, 1.0/4.0, gks->fnew, &gks->local_ext);
            gk_species_bflux_combine(app, &gks->bflux, gks->bflux.f1, 3.0/4.0, gks->bflux.f, 1.0/4.0, gks->bflux.fnew);
          }
          for (int i=0; i<ns_neut; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            gk_neut_species_combine(gkns, gkns->f1, 3.0/4.0, gkns->f, 1.0/4.0, gkns->fnew, &gkns->local_ext);
            gk_neut_species_bflux_combine(app, &gkns->bflux, gkns->bflux.f1, 3.0/4.0, gkns->bflux.f, 1.0/4.0, gkns->bflux.fnew);
          }
          app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

          // Post process stage.
          for (int i=0; i<ns_charged; ++i) {
            struct gk_species *gks = &app->species[i];
            fout[i] = gks->f1;
            bflux_out[i] = gks->bflux.f1;
          }
          for (int i=0; i<ns_neut; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            fout_neut[i] = gkns->f1;
            bflux_out_neut[i] = gkns->bflux.f1;
          }
          for (int i=0; i<nfields; ++i) {
            struct gk_field *gkf = app->field;
            fieldout[i] = gkf->f1[i];
          }
          gyrokinetic_post_process_rk_stage_intermediate_ssprk(app, tcurr, dt, fdot_args);

          state = SSPRK_3_3_STAGE_2;
        }
        break;

      case SSPRK_3_3_STAGE_2:
        for (int i=0; i<ns_charged; ++i) {
          struct gk_species *gks = &app->species[i];
          fin[i] = gks->f1;
          fout[i] = gks->fnew;
          bflux_in[i] = gks->bflux.f1;
          bflux_out[i] = gks->bflux.fnew;
        }
        for (int i=0; i<ns_neut; ++i) {
          struct gk_neut_species *gkns = &app->neut_species[i];
          fin_neut[i] = gkns->f1;
          fout_neut[i] = gkns->fnew;
          bflux_in_neut[i] = gkns->bflux.f1;
          bflux_out_neut[i] = gkns->bflux.fnew;
        }
        for (int i=0; i<nfields; ++i) {
          struct gk_field *gkf = app->field;
          fieldin[i] = gkf->f1[i];
          fieldout[i] = gkf->fnew[i];
        }

        // Step solution.
        gyrokinetic_forward_euler(app, tcurr+dt/2, dt, fdot_args, &st);

        if (st.dt_actual < dt) {

          // Failed step. Collect stats.
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_3_fail += 1;

          dt = st.dt_actual;
          state = SSPRK_3_3_STAGE_0; // Restart from stage 0.

          app->stat.nstage_2_fail += 1;

        }
        else {
          struct timespec wst = gkyl_wall_clock();
          for (int i=0; i<ns_charged; ++i) {
            struct gk_species *gks = &app->species[i];
            // Step f.
            gk_species_combine(gks, gks->f, 1.0/3.0, gks->f, 2.0/3.0, gks->fnew, &gks->local_ext);
            // Step boundary fluxes.
            gk_species_bflux_combine(app, &gks->bflux, gks->bflux.f, 1.0/3.0, gks->bflux.f, 2.0/3.0, gks->bflux.fnew);
          }

          for (int i=0; i<ns_neut; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            // Step f.
            gk_neut_species_combine(gkns, gkns->f, 1.0/3.0, gkns->f, 2.0/3.0, gkns->fnew, &gkns->local_ext);
            // Step boundary fluxes.
            gk_neut_species_bflux_combine(app, &gkns->bflux, gkns->bflux.f, 1.0/3.0, gkns->bflux.f, 2.0/3.0, gkns->bflux.fnew);
          }
          app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

          state = SSPRK_3_3_STAGE_COMPLETE;
        }
        break;

      case SSPRK_3_3_STAGE_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }

  // Post-process step.
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    fout[i] = gks->f;
    bflux_out[i] = gks->bflux.f;
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    fout_neut[i] = gkns->f;
    bflux_out_neut[i] = gkns->bflux.f;
  }
  for (int i=0; i<nfields; ++i) {
    struct gk_field *gkf = app->field;
    fieldout[i] = gkf->f[i];
  }
  gyrokinetic_post_process_step_ssprk(app, tcurr, dt, fdot_args);

  return st;
}

void
gyrokinetic_pre_process_step_sts(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the beginning of an STS step.
}

void
gyrokinetic_pre_process_rk_stage_sts(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args, int stage_idx, int num_stages)
{
  // Operations at the end of an RK stage.
}

void
gyrokinetic_post_process_rk_stage_sts(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args, int stage_idx, int num_stages)
{
  // Operations at the end of an STS stage.
}

void
gyrokinetic_post_process_step_sts(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the end of the STS step.
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;
  struct gkyl_array ***bflux_out_neut = fdot_args->bflux_out_neut;
  struct gkyl_array **fieldout = fdot_args->fieldout;

  // Compute the fields and apply BCs
  gyrokinetic_calc_field(app, tcurr, fout, bflux_out, fieldout);
  gyrokinetic_apply_bc(app, tcurr, fieldout, fout, fout_neut);
}

void
gyrokinetic_pre_process_step_opsplit(gkyl_gyrokinetic_app* app, double tcurr,
  double dt, struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the beginning of an operator split step.
}

void
gyrokinetic_post_process_step_opsplit(gkyl_gyrokinetic_app* app, double tcurr,
  double dt, struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the end of the operator split step.
}

