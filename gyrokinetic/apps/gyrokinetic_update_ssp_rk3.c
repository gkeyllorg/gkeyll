#include <gkyl_gyrokinetic_priv.h>

void
gyrokinetic_pre_process_step(
  gkyl_gyrokinetic_app *app, double tcurr, const struct gkyl_gyrokinetic_fdot_args *args
)
{
  // Adapt sources.
  for (int i = 0; i < app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_source_adapt(app, gks, &gks->src, gks->lte.f_lte, args->smom_in, tcurr);
  }
}

// Capture old-state diagnostics while the app still holds the initial field.
void
gyrokinetic_save_old_step_diagnostics(
  gkyl_gyrokinetic_app *app, double dt, const struct gkyl_gyrokinetic_fdot_args *args
)
{
  for (int i = 0; i < app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_calc_integ_mom_dt(app, gks, args->fin[i], dt, gks->fdot_mom_old);
  }
  gk_field_calc_energy_dt(app, app->field, dt, true);
}

// Complete single-block bookkeeping after accepted-state field/BC preparation.
void
gyrokinetic_post_process_step(
  gkyl_gyrokinetic_app *app, double t_end, double dt, struct gkyl_gyrokinetic_fdot_args *args
)
{
  for (int i = 0; i < app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_calc_integ_mom_dt(app, gks, args->fout[i], dt, gks->fdot_mom_new);
  }

  // Preserve the existing order: new moments, species scaling, then field energy.
  for (int i = 0; i < app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_scaling_apply(app, gks, &gks->sca, args->fout[i], args->smom_out);
  }
  for (int i = 0; i < app->num_neut_species; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_scaling_apply(app, gkns, &gkns->sca, args->fout_neut[i], args->smom_out);
  }
  gk_field_calc_energy_dt(app, app->field, dt, false);
}

static void
gyrokinetic_forward_euler(
  gkyl_gyrokinetic_app *app, double tcurr, double dt, struct gkyl_gyrokinetic_fdot_args *fdot_args,
  struct gkyl_update_status *st
)
{
  const struct gkyl_array **fin = fdot_args->fin;
  const struct gkyl_array **fin_neut = fdot_args->fin_neut;
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***smom_in = fdot_args->smom_in;
  struct gkyl_array ***smom_in_neut = fdot_args->smom_in_neut;
  struct gkyl_array ***smom_out = fdot_args->smom_out;
  struct gkyl_array ***smom_out_neut = fdot_args->smom_out_neut;

  struct timespec wst_fe = gkyl_wall_clock();
  // Take a forward Euler step with the suggested time-step dt. This may
  // not be the actual time-step taken. However, the function will never
  // take a time-step larger than dt even if it is allowed by
  // stability. The actual time-step and dt_suggested are returned in
  // the status object.
  app->stat.nfeuler += 1;

  // Compute the time rate of change of the distributions, df/dt.
  gyrokinetic_rhs(app, tcurr, dt, fdot_args, st);

  struct timespec wst = gkyl_wall_clock();
  // Complete update of distribution functions.
  double dta = st->dt_actual;
  for (int i = 0; i < app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_step_f(gks, fout[i], dta, fin[i]);
    gk_stepped_moms_accumulate(&gks->step_moms, smom_out[i], 1.0, smom_in[i]);
  }
  for (int i = 0; i < app->num_neut_species; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_step_f(gkns, fout_neut[i], dta, fin_neut[i]);
    gk_stepped_moms_accumulate(&gkns->step_moms, smom_out_neut[i], 1.0, smom_in_neut[i]);
  }
  app->stat.fwd_euler_step_f_tm += gkyl_time_diff_now_sec(wst);
  app->stat.fwd_euler_tm += gkyl_time_diff_now_sec(wst_fe);
}

struct gkyl_update_status
gyrokinetic_update_ssp_rk3(gkyl_gyrokinetic_app *app, double dt0)
{
  // Take time-step using the RK3 method. Also sets the status object
  // which has the actual and suggested dts used. These can be different
  // from the actual time-step.
  struct gkyl_gyrokinetic_fdot_args *fdot_args = &app->fdot_args;
  const struct gkyl_array **fin = fdot_args->fin;
  const struct gkyl_array **fin_neut = fdot_args->fin_neut;
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***smom_in = fdot_args->smom_in;
  struct gkyl_array ***smom_in_neut = fdot_args->smom_in_neut;
  struct gkyl_array ***smom_out = fdot_args->smom_out;
  struct gkyl_array ***smom_out_neut = fdot_args->smom_out_neut;

  struct gkyl_update_status st = {.success = true};

  enum {
    SSPRK_3_3_STAGE_0 = 0,
    SSPRK_3_3_STAGE_1,
    SSPRK_3_3_STAGE_2,
    SSPRK_3_3_STAGE_COMPLETE
  } state = SSPRK_3_3_STAGE_0; // Time-stepper state.

  double tcurr = app->tcurr, dt = dt0;

  // Perform operations before RK stages begin.
  for (int i = 0; i < app->num_species; ++i) {
    fin[i] = app->species[i].f;
    smom_in[i] = app->species[i].step_moms.f;
  }
  for (int i = 0; i < app->num_neut_species; ++i) {
    fin_neut[i] = app->neut_species[i].f;
    smom_in_neut[i] = app->neut_species[i].step_moms.f;
  }
  gyrokinetic_pre_process_step(app, tcurr, fdot_args);

  while (state != SSPRK_3_3_STAGE_COMPLETE) {
    switch (state) {
      case SSPRK_3_3_STAGE_0:
        for (int i = 0; i < app->num_species; ++i) {
          struct gk_species *gks = &app->species[i];
          fin[i] = gks->f;
          fout[i] = gks->f1;
          // Stepped moments.
          smom_in[i] = gks->step_moms.f;
          smom_out[i] = gks->step_moms.f1;
        }
        for (int i = 0; i < app->num_neut_species; ++i) {
          struct gk_neut_species *gkns = &app->neut_species[i];
          fin_neut[i] = gkns->f;
          fout_neut[i] = gkns->f1;
          // Stepped moments.
          smom_in_neut[i] = gkns->step_moms.f;
          smom_out_neut[i] = gkns->step_moms.f1;
        }

        gyrokinetic_forward_euler(app, tcurr, dt, fdot_args, &st);
        dt = st.dt_actual;

        // Subtract the old moments from f1 to accumulate only the current
        // timestep's contribution. Keep the old moments in f for RK-stage retries.
        struct timespec wst = gkyl_wall_clock();
        for (int i = 0; i < app->num_species; ++i) {
          struct gk_species *gks = &app->species[i];
          gk_stepped_moms_accumulate(&gks->step_moms, smom_out[i], -1.0, smom_in[i]);
        }
        for (int i = 0; i < app->num_neut_species; ++i) {
          struct gk_neut_species *gkns = &app->neut_species[i];
          gk_stepped_moms_accumulate(&gkns->step_moms, smom_out_neut[i], -1.0, smom_in_neut[i]);
        }
        app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

        gyrokinetic_save_old_step_diagnostics(app, dt, fdot_args);

        // Compute the fields and apply BCs.
        gyrokinetic_calc_field_and_apply_bc(app, tcurr + dt, fout, smom_out, fout_neut);

        state = SSPRK_3_3_STAGE_1;
        break;

      case SSPRK_3_3_STAGE_1:
        for (int i = 0; i < app->num_species; ++i) {
          struct gk_species *gks = &app->species[i];
          fin[i] = gks->f1;
          fout[i] = gks->fnew;
          // Stepped moments.
          smom_in[i] = gks->step_moms.f1;
          smom_out[i] = gks->step_moms.fnew;
        }
        for (int i = 0; i < app->num_neut_species; ++i) {
          struct gk_neut_species *gkns = &app->neut_species[i];
          fin_neut[i] = gkns->f1;
          fout_neut[i] = gkns->fnew;
          // Stepped moments.
          smom_in_neut[i] = gkns->step_moms.f1;
          smom_out_neut[i] = gkns->step_moms.fnew;
        }

        gyrokinetic_forward_euler(app, tcurr + dt, dt, fdot_args, &st);

        if (st.dt_actual < dt) {
          // Restore the initial field; the saved state already has valid ghost cells.
          for (int i = 0; i < app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            fin[i] = gks->f;
            smom_in[i] = gks->step_moms.f;
          }
          gyrokinetic_calc_field(app, tcurr, fin, smom_in);

          // Collect stats.
          double dt_rel_diff = (dt - st.dt_actual) / st.dt_actual;
          app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0], dt_rel_diff);
          app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1], dt_rel_diff);
          app->stat.nstage_2_fail += 1;

          dt = st.dt_actual;
          state = SSPRK_3_3_STAGE_0; // Restart from stage 0.

        } else {
          struct timespec wst = gkyl_wall_clock();
          for (int i = 0; i < app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            gk_species_combine(
              gks, gks->f1, 3.0 / 4.0, gks->f, 1.0 / 4.0, gks->fnew, &gks->local_ext
            );
            gk_stepped_moms_set(&gks->step_moms, gks->step_moms.f1, 1.0 / 4.0, gks->step_moms.fnew);
          }
          for (int i = 0; i < app->num_neut_species; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            gk_neut_species_combine(
              gkns, gkns->f1, 3.0 / 4.0, gkns->f, 1.0 / 4.0, gkns->fnew, &gkns->local_ext
            );
            gk_stepped_moms_set(
              &gkns->step_moms, gkns->step_moms.f1, 1.0 / 4.0, gkns->step_moms.fnew
            );
          }
          app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

          // Compute the fields and apply BCs.
          for (int i = 0; i < app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            fout[i] = gks->f1;
            smom_out[i] = gks->step_moms.f1;
          }
          for (int i = 0; i < app->num_neut_species; ++i) {
            fout_neut[i] = app->neut_species[i].f1;
          }
          gyrokinetic_calc_field_and_apply_bc(app, tcurr + dt / 2, fout, smom_out, fout_neut);

          state = SSPRK_3_3_STAGE_2;
        }
        break;

      case SSPRK_3_3_STAGE_2:
        for (int i = 0; i < app->num_species; ++i) {
          struct gk_species *gks = &app->species[i];
          fin[i] = gks->f1;
          fout[i] = gks->fnew;
          // Stepped moments.
          smom_in[i] = gks->step_moms.f1;
          smom_out[i] = gks->step_moms.fnew;
        }
        for (int i = 0; i < app->num_neut_species; ++i) {
          struct gk_neut_species *gkns = &app->neut_species[i];
          fin_neut[i] = gkns->f1;
          fout_neut[i] = gkns->fnew;
          // Stepped moments.
          smom_in_neut[i] = gkns->step_moms.f1;
          smom_out_neut[i] = gkns->step_moms.fnew;
        }

        gyrokinetic_forward_euler(app, tcurr + dt / 2, dt, fdot_args, &st);

        if (st.dt_actual < dt) {
          // Restore the initial field; the saved state already has valid ghost cells.
          for (int i = 0; i < app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            fin[i] = gks->f;
            smom_in[i] = gks->step_moms.f;
          }
          gyrokinetic_calc_field(app, tcurr, fin, smom_in);

          // Collect stats.
          double dt_rel_diff = (dt - st.dt_actual) / st.dt_actual;
          app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0], dt_rel_diff);
          app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1], dt_rel_diff);
          app->stat.nstage_3_fail += 1;

          dt = st.dt_actual;
          state = SSPRK_3_3_STAGE_0; // Restart from stage 0.

          app->stat.nstage_2_fail += 1;
        } else {
          struct timespec wst = gkyl_wall_clock();
          for (int i = 0; i < app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            // Step f.
            gk_species_combine(
              gks, gks->f1, 1.0 / 3.0, gks->f, 2.0 / 3.0, gks->fnew, &gks->local_ext
            );
            gk_species_copy_range(gks, gks->f, gks->f1, &gks->local_ext);
            // Update stepped moments.
            gk_stepped_moms_set(&gks->step_moms, gks->step_moms.f, 2.0 / 3.0, gks->step_moms.fnew);
          }

          for (int i = 0; i < app->num_neut_species; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            gk_neut_species_combine(
              gkns, gkns->f1, 1.0 / 3.0, gkns->f, 2.0 / 3.0, gkns->fnew, &gkns->local_ext
            );
            gk_neut_species_copy_range(gkns, gkns->f, gkns->f1, &gkns->local_ext);
            // Update stepped moments.
            gk_stepped_moms_set(
              &gkns->step_moms, gkns->step_moms.f, 2.0 / 3.0, gkns->step_moms.fnew
            );
          }
          app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

          // Apply positivity shift if requested.
          for (int i = 0; i < app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            gk_species_positivity_apply(app, gks, &gks->positivity, gks->fnew, gks->f);
          }

          for (int i = 0; i < app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            fout[i] = gks->f;
            smom_out[i] = gks->step_moms.f;
          }
          for (int i = 0; i < app->num_neut_species; ++i) {
            fout_neut[i] = app->neut_species[i].f;
          }

          // Enforce quasineutrality of the positivity shifts.
          gyrokinetic_post_positivity_quasineut(app, fout);

          // Compute the fields and apply BCs
          gyrokinetic_calc_field_and_apply_bc(app, tcurr + dt, fout, smom_out, fout_neut);

          state = SSPRK_3_3_STAGE_COMPLETE;
        }
        break;

      case SSPRK_3_3_STAGE_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }

  // Perform operations after RK stages end.
  gyrokinetic_post_process_step(app, tcurr + dt, dt, fdot_args);

  return st;
}
