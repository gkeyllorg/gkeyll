#include <gkyl_pkpm_priv.h>

// Take time-step using the RK3 method. Also sets the status object
// which has the actual and suggested dts used. These can be different
// from the actual time-step.
struct gkyl_update_status
pkpm_update_explicit_ssp_rk3(gkyl_pkpm_app* app, double dt0)
{
  int ns = app->num_species;  

  const struct gkyl_array *fin[ns];
  struct gkyl_array *fout[ns];
  const struct gkyl_array *fluidin[ns];
  struct gkyl_array *fluidout[ns];
  struct gkyl_update_status st = { .success = true };

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        do {
          struct timespec rk3_s1_tm = gkyl_wall_clock();

          for (int i=0; i<ns; ++i) {
            fin[i] = app->species[i].f;
            fluidin[i] = app->species[i].fluid;
            fout[i] = app->species[i].f1;
            fluidout[i] = app->species[i].fluid1;
          }
          pkpm_forward_euler(app, tcurr, dt, fin, fluidin, app->field->em,
            fout, fluidout, app->field->em1,
            &st
          );
          // Limit fluid and EM solutions if desired (done after update as post-hoc fix)
          for (int i=0; i<ns; ++i) {
            pkpm_fluid_species_limiter(app, &app->species[i], fout[i], fluidout[i]);
          }
          pkpm_field_limiter(app, app->field, app->field->emnew);

          dt = st.dt_actual;
          state = RK_STAGE_2;

          app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s1_tm);
        } while(0);
        break;

      case RK_STAGE_2:
        do {
          struct timespec rk3_s2_tm = gkyl_wall_clock();

          for (int i=0; i<ns; ++i) {
            fin[i] = app->species[i].f1;
            fluidin[i] = app->species[i].fluid1;
            fout[i] = app->species[i].fnew;
            fluidout[i] = app->species[i].fluidnew;
          }
          pkpm_forward_euler(app, tcurr+dt, dt, fin, fluidin, app->field->em1,
            fout, fluidout, app->field->emnew,
            &st
          );
          // Limit fluid and EM solutions if desired (done after update as post-hoc fix)
          for (int i=0; i<ns; ++i) {
            pkpm_fluid_species_limiter(app, &app->species[i], fout[i], fluidout[i]);
          }
          pkpm_field_limiter(app, app->field, app->field->emnew);

          if (st.dt_actual < dt) {
            // collect stats
            double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
            app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0],
              dt_rel_diff);
            app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1],
              dt_rel_diff);
            app->stat.nstage_2_fail += 1;

            dt = st.dt_actual;
            state = RK_STAGE_1; // restart from stage 1
          } 
          else {
            for (int i=0; i<ns; ++i) {
              array_combine(app->species[i].f1,
                3.0/4.0, app->species[i].f, 1.0/4.0, app->species[i].fnew, &app->species[i].local_ext);
            }
            for (int i=0; i<ns; ++i) {
              array_combine(app->species[i].fluid1,
                3.0/4.0, app->species[i].fluid, 1.0/4.0, app->species[i].fluidnew, &app->local_ext);
            }
            array_combine(app->field->em1,
              3.0/4.0, app->field->em, 1.0/4.0, app->field->emnew, &app->local_ext);

            state = RK_STAGE_3;
          }

          app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s2_tm);
        } while(0);
        break;

      case RK_STAGE_3:
        do {
          struct timespec rk3_s3_tm = gkyl_wall_clock();

          for (int i=0; i<ns; ++i) {
            fin[i] = app->species[i].f1;
            fluidin[i] = app->species[i].fluid1;
            fout[i] = app->species[i].fnew;
            fluidout[i] = app->species[i].fluidnew;
          }
          pkpm_forward_euler(app, tcurr+dt/2, dt, fin, fluidin, app->field->em1,
            fout, fluidout, app->field->emnew,
            &st
          );
          // Limit fluid and EM solutions if desired (done after update as post-hoc fix)
          for (int i=0; i<ns; ++i) {
            pkpm_fluid_species_limiter(app, &app->species[i], fout[i], fluidout[i]);
          }
          pkpm_field_limiter(app, app->field, app->field->emnew);

          if (st.dt_actual < dt) {
            // collect stats
            double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
            app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0],
              dt_rel_diff);
            app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1],
              dt_rel_diff);
            app->stat.nstage_3_fail += 1;

            dt = st.dt_actual;
            state = RK_STAGE_1; // restart from stage 1

            app->stat.nstage_2_fail += 1;
          }
          else {
            for (int i=0; i<ns; ++i) {
              array_combine(app->species[i].f1,
                1.0/3.0, app->species[i].f, 2.0/3.0, app->species[i].fnew, &app->species[i].local_ext);
              gkyl_array_copy_range(app->species[i].f, app->species[i].f1, &app->species[i].local_ext);
            }
            for (int i=0; i<ns; ++i) {
              array_combine(app->species[i].fluid1,
                1.0/3.0, app->species[i].fluid, 2.0/3.0, app->species[i].fluidnew, &app->local_ext);
              gkyl_array_copy_range(app->species[i].fluid, app->species[i].fluid1, &app->local_ext);
            }
            array_combine(app->field->em1,
              1.0/3.0, app->field->em, 2.0/3.0, app->field->emnew, &app->local_ext);
            gkyl_array_copy_range(app->field->em, app->field->em1, &app->local_ext);

            state = RK_COMPLETE;
          }
        
          app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s3_tm);
        } while(0);
        break;

      case RK_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }

  return st;
}