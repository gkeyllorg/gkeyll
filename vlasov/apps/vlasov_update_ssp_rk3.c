#include <gkyl_vlasov_priv.h>

// Take a time-step with SSP-RK3; the status holds the step taken and the
// step suggested, the smallest of the bounds of the three stages of the step
// taken: each stage evaluates its bound on its own state, and a step suggested
// from the last stage alone is retaken wherever an earlier stage allows less
// (the state after the forward-Euler step of stage 1 is the farthest from the
// solution). With retake, stage 1 takes the step its CFL bound allows
// and a smaller bound in stage 2 or 3 restarts from stage 1; without, the
// function returns at the first stage whose bound is below dt0, with only
// the stage arrays written and dt_actual < dt0 in the status. What the
// species record over the stages (vlasov_species_discard_step) is discarded
// on a restart, and left to the caller on that return.
static inline struct gkyl_update_status
update_ssp_rk3(gkyl_vlasov_app *app, double dt0, bool retake)
{
  int num_species = app->num_species;

  struct gkyl_array *fin[num_species];
  struct gkyl_array *fout[num_species];
  const struct gkyl_array *fluidin[num_species];
  struct gkyl_array *fluidout[num_species];
  struct gkyl_update_status st = {.success = true};

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  double dt_suggested = DBL_MAX; // smallest bound of the stages of the step being taken
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        do {
          struct timespec rk3_s1_tm = gkyl_wall_clock();

          vlasov_species_gather_rk_state(app, VM_RK_F, VM_RK_F1, fin, fout, fluidin, fluidout);
          vlasov_forward_euler(
            app, 1, tcurr, dt, fin, fluidin, app->field->em, fout, fluidout, app->field->em1, &st
          );

          if (!retake && st.dt_actual < dt) {
            app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s1_tm);
            return st;
          }

          dt_suggested = st.dt_suggested;
          vm_apply_bc(app, tcurr, fout, fluidout, app->field->em1);

          // Limit fluid and EM solutions if desired (done after update as post-hoc fix)
          for (int i = 0; i < num_species; ++i) {
            vlasov_species_limiter(app, &app->species[i], fluidout[i]);
          }
          vlasov_field_limiter(app, app->field->em1); // no-op for the null field
          dt = st.dt_actual;
          vlasov_stage_implicit(app, 1, dt);
          state = RK_STAGE_2;

          app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s1_tm);
        } while (0);
        break;

      case RK_STAGE_2:
        do {
          struct timespec rk3_s2_tm = gkyl_wall_clock();

          vlasov_species_gather_rk_state(app, VM_RK_F1, VM_RK_FNEW, fin, fout, fluidin, fluidout);
          vlasov_forward_euler(
            app, 2, tcurr + dt, dt, fin, fluidin, app->field->em1, fout, fluidout,
            app->field->emnew, &st
          );

          dt_suggested = fmin(dt_suggested, st.dt_suggested);
          vm_apply_bc(app, tcurr, fout, fluidout, app->field->emnew);

          // Limit fluid and EM solutions if desired (done after update as post-hoc fix)
          for (int i = 0; i < num_species; ++i) {
            vlasov_species_limiter(app, &app->species[i], fluidout[i]);
          }
          vlasov_field_limiter(app, app->field->emnew); // no-op for the null field
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
            for (int i = 0; i < num_species; ++i) {
              vlasov_species_discard_step(app, &app->species[i]);
            }
          } else {
            for (int i = 0; i < num_species; ++i) {
              vlasov_species_combine(
                app, &app->species[i], VM_RK_F1, 3.0 / 4.0, VM_RK_F, 1.0 / 4.0, VM_RK_FNEW
              );
            }
            vlasov_field_combine(
              app, app->field->em1, 3.0 / 4.0, app->field->em, 1.0 / 4.0, app->field->emnew
            ); // no-op for null field
            vlasov_stage_implicit(app, 2, dt);

            state = RK_STAGE_3;
          }

          app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s2_tm);
        } while (0);
        break;

      case RK_STAGE_3:
        do {
          struct timespec rk3_s3_tm = gkyl_wall_clock();

          vlasov_species_gather_rk_state(app, VM_RK_F1, VM_RK_FNEW, fin, fout, fluidin, fluidout);
          vlasov_forward_euler(
            app, 3, tcurr + dt / 2, dt, fin, fluidin, app->field->em1, fout, fluidout,
            app->field->emnew, &st
          );

          dt_suggested = fmin(dt_suggested, st.dt_suggested);
          vm_apply_bc(app, tcurr, fout, fluidout, app->field->emnew);

          // Limit fluid and EM solutions if desired (done after update as post-hoc fix)
          for (int i = 0; i < num_species; ++i) {
            vlasov_species_limiter(app, &app->species[i], fluidout[i]);
          }
          vlasov_field_limiter(app, app->field->emnew); // no-op for the null field
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
            for (int i = 0; i < num_species; ++i) {
              vlasov_species_discard_step(app, &app->species[i]);
            }

            app->stat.nstage_2_fail += 1;
          } else {
            for (int i = 0; i < num_species; ++i) {
              vlasov_species_combine(
                app, &app->species[i], VM_RK_F1, 1.0 / 3.0, VM_RK_F, 2.0 / 3.0, VM_RK_FNEW
              );
            }
            vlasov_field_combine(
              app, app->field->em1, 1.0 / 3.0, app->field->em, 2.0 / 3.0, app->field->emnew
            ); // no-op for null field
            // The implicit step of the stage acts on the combined state, which
            // then becomes the solution.
            vlasov_stage_implicit(app, 3, dt);
            for (int i = 0; i < num_species; ++i) {
              vlasov_species_complete_step(app, &app->species[i]);
              vlasov_species_copy_range(app, &app->species[i], VM_RK_F, VM_RK_F1);
            }
            vlasov_field_copy_range(app, app->field->em, app->field->em1); // no-op for null field

            st.dt_suggested = dt_suggested;
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
vlasov_update_ssp_rk3(gkyl_vlasov_app *app, double dt0)
{
  return update_ssp_rk3(app, dt0, true);
}

struct gkyl_update_status
vlasov_update_ssp_rk3_attempt(gkyl_vlasov_app *app, double dt0)
{
  return update_ssp_rk3(app, dt0, false);
}
