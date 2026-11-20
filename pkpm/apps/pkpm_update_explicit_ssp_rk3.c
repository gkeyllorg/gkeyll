#include <gkyl_pkpm_priv.h>

// Take a time-step with the SSP-RK3 method. The status object holds the step
// taken and the step suggested for the next update. With retake, stage 1 takes
// the step its CFL bound allows and a smaller bound in stage 2 or 3 restarts
// from stage 1; without, the function returns at the first stage whose bound is
// below dt0, with only the stage arrays written and dt_actual < dt0 in the
// status.
static inline struct gkyl_update_status
update_explicit_ssp_rk3(gkyl_pkpm_app *app, double dt0, bool retake)
{
  int ns = app->num_species;

  struct gkyl_array *fin[ns];
  struct gkyl_array *fout[ns];
  const struct gkyl_array *fluidin[ns];
  struct gkyl_array *fluidout[ns];
  struct gkyl_update_status st = {.success = true};

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        do {
          struct timespec rk3_s1_tm = gkyl_wall_clock();

          for (int i = 0; i < ns; ++i) {
            fin[i] = app->species[i].kinetic->f;
            fluidin[i] = app->species[i].fluid->fluid;
            fout[i] = app->species[i].kinetic->f1;
            fluidout[i] = app->species[i].fluid->fluid1;
          }
          pkpm_forward_euler(
            app, 1, tcurr, dt, fin, fluidin, app->field->em, fout, fluidout, app->field->em1, &st
          );

          if (!retake && st.dt_actual < dt) {
            app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s1_tm);
            return st;
          }

          dt = st.dt_actual;
          pkpm_coll_stage_advance(app, 1, dt);
          state = RK_STAGE_2;

          app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s1_tm);
        } while (0);
        break;

      case RK_STAGE_2:
        do {
          struct timespec rk3_s2_tm = gkyl_wall_clock();

          for (int i = 0; i < ns; ++i) {
            fin[i] = app->species[i].kinetic->f1;
            fluidin[i] = app->species[i].fluid->fluid1;
            fout[i] = app->species[i].kinetic->fnew;
            fluidout[i] = app->species[i].fluid->fluidnew;
          }
          pkpm_forward_euler(
            app, 2, tcurr + dt, dt, fin, fluidin, app->field->em1, fout, fluidout,
            app->field->emnew, &st
          );

          if (st.dt_actual < dt) {
            // collect stats
            double dt_rel_diff = (dt - st.dt_actual) / st.dt_actual;
            app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0], dt_rel_diff);
            app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1], dt_rel_diff);
            app->stat.nstage_2_fail += 1;
            if (!retake) {
              app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s2_tm);
              return st;
            }

            dt = st.dt_actual;
            state = RK_STAGE_1; // restart from stage 1
          } else {
            for (int i = 0; i < ns; ++i) {
              array_combine(
                app->species[i].kinetic->f1, 3.0 / 4.0, app->species[i].kinetic->f, 1.0 / 4.0,
                app->species[i].kinetic->fnew, &app->species[i].kinetic->local_ext
              );
            }
            for (int i = 0; i < ns; ++i) {
              array_combine(
                app->species[i].fluid->fluid1, 3.0 / 4.0, app->species[i].fluid->fluid, 1.0 / 4.0,
                app->species[i].fluid->fluidnew, &app->local_ext
              );
            }
            array_combine(
              app->field->em1, 3.0 / 4.0, app->field->em, 1.0 / 4.0, app->field->emnew,
              &app->local_ext
            );
            pkpm_coll_stage_advance(app, 2, dt);

            state = RK_STAGE_3;
          }

          app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s2_tm);
        } while (0);
        break;

      case RK_STAGE_3:
        do {
          struct timespec rk3_s3_tm = gkyl_wall_clock();

          for (int i = 0; i < ns; ++i) {
            fin[i] = app->species[i].kinetic->f1;
            fluidin[i] = app->species[i].fluid->fluid1;
            fout[i] = app->species[i].kinetic->fnew;
            fluidout[i] = app->species[i].fluid->fluidnew;
          }
          pkpm_forward_euler(
            app, 3, tcurr + dt / 2, dt, fin, fluidin, app->field->em1, fout, fluidout,
            app->field->emnew, &st
          );

          if (st.dt_actual < dt) {
            // collect stats
            double dt_rel_diff = (dt - st.dt_actual) / st.dt_actual;
            app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0], dt_rel_diff);
            app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1], dt_rel_diff);
            app->stat.nstage_3_fail += 1;
            if (!retake) {
              app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s3_tm);
              return st;
            }

            dt = st.dt_actual;
            state = RK_STAGE_1; // restart from stage 1

            app->stat.nstage_2_fail += 1;
          } else {
            for (int i = 0; i < ns; ++i) {
              array_combine(
                app->species[i].kinetic->f1, 1.0 / 3.0, app->species[i].kinetic->f, 2.0 / 3.0,
                app->species[i].kinetic->fnew, &app->species[i].kinetic->local_ext
              );
            }
            for (int i = 0; i < ns; ++i) {
              array_combine(
                app->species[i].fluid->fluid1, 1.0 / 3.0, app->species[i].fluid->fluid, 2.0 / 3.0,
                app->species[i].fluid->fluidnew, &app->local_ext
              );
            }
            pkpm_coll_stage_advance(app, 3, dt);
            for (int i = 0; i < ns; ++i) {
              pkpm_species_complete_step(app, &app->species[i]);
              gkyl_array_copy_range(
                app->species[i].kinetic->f, app->species[i].kinetic->f1,
                &app->species[i].kinetic->local_ext
              );
              gkyl_array_copy_range(
                app->species[i].fluid->fluid, app->species[i].fluid->fluid1, &app->local_ext
              );
            }
            array_combine(
              app->field->em1, 1.0 / 3.0, app->field->em, 2.0 / 3.0, app->field->emnew,
              &app->local_ext
            );
            gkyl_array_copy_range(app->field->em, app->field->em1, &app->local_ext);

            state = RK_COMPLETE;
          }

          app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s3_tm);
        } while (0);
        break;

      case RK_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }

  return st;
}

struct gkyl_update_status
pkpm_update_explicit_ssp_rk3(gkyl_pkpm_app *app, double dt0)
{
  return update_explicit_ssp_rk3(app, dt0, true);
}

struct gkyl_update_status
pkpm_update_explicit_ssp_rk3_attempt(gkyl_pkpm_app *app, double dt0)
{
  return update_explicit_ssp_rk3(app, dt0, false);
}
