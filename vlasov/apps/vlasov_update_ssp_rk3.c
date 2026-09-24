#include <gkyl_vlasov_priv.h>

// Take time-step using the RK3 method. Also sets the status object
// which has the actual and suggested dts used. These can be different
// from the actual time-step.
struct gkyl_update_status
vlasov_update_ssp_rk3(gkyl_vlasov_app* app, double dt0)
{
  int num_species = app->num_species;

  const struct gkyl_array *fin[num_species];
  struct gkyl_array *fout[num_species];
  const struct gkyl_array *fluidin[num_species];
  struct gkyl_array *fluidout[num_species];
  struct gkyl_update_status st = { .success = true };

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        do {
          struct timespec rk3_s1_tm = gkyl_wall_clock();

          vlasov_species_gather_rk_state(app, VM_RK_F, VM_RK_F1, fin, fout, fluidin, fluidout);
          vlasov_forward_euler(app, tcurr, dt, fin, fluidin, app->field->em,
            fout, fluidout, app->field->em1,
            &st
          );

          vm_apply_bc(app, tcurr, fout, fluidout, app->field->em1);

          // Limit fluid and EM solutions if desired (done after update as post-hoc fix)
          for (int i=0; i<num_species; ++i)
            vlasov_species_limiter(app, &app->species[i], fluidout[i]);
          vlasov_field_limiter(app, app->field->em1); // no-op for the null field
          dt = st.dt_actual;
          state = RK_STAGE_2;

          app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s1_tm);
        } while(0);
        break;

      case RK_STAGE_2:
        do {
          struct timespec rk3_s2_tm = gkyl_wall_clock();

          vlasov_species_gather_rk_state(app, VM_RK_F1, VM_RK_FNEW, fin, fout, fluidin, fluidout);
          vlasov_forward_euler(app, tcurr+dt, dt, fin, fluidin, app->field->em1,
            fout, fluidout, app->field->emnew,
            &st
          );

          vm_apply_bc(app, tcurr, fout, fluidout, app->field->emnew);

          // Limit fluid and EM solutions if desired (done after update as post-hoc fix)
          for (int i=0; i<num_species; ++i)
            vlasov_species_limiter(app, &app->species[i], fluidout[i]);
          vlasov_field_limiter(app, app->field->emnew); // no-op for the null field
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
            for (int i=0; i<num_species; ++i)
              vlasov_species_combine(app, &app->species[i], VM_RK_F1, 3.0/4.0, VM_RK_F, 1.0/4.0, VM_RK_FNEW);
            vlasov_field_combine(app, app->field->em1,
              3.0/4.0, app->field->em, 1.0/4.0, app->field->emnew); // no-op for null field

            state = RK_STAGE_3;
          }

          app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s2_tm);
        } while(0);
        break;

      case RK_STAGE_3:
        do {
          struct timespec rk3_s3_tm = gkyl_wall_clock();

          vlasov_species_gather_rk_state(app, VM_RK_F1, VM_RK_FNEW, fin, fout, fluidin, fluidout);
          vlasov_forward_euler(app, tcurr+dt/2, dt, fin, fluidin, app->field->em1,
            fout, fluidout, app->field->emnew,
            &st
          );

          vm_apply_bc(app, tcurr, fout, fluidout, app->field->emnew);

          // Limit fluid and EM solutions if desired (done after update as post-hoc fix)
          for (int i=0; i<num_species; ++i)
            vlasov_species_limiter(app, &app->species[i], fluidout[i]);
          vlasov_field_limiter(app, app->field->emnew); // no-op for the null field
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
            for (int i=0; i<num_species; ++i) {
              vlasov_species_combine(app, &app->species[i], VM_RK_F1, 1.0/3.0, VM_RK_F, 2.0/3.0, VM_RK_FNEW);
              vlasov_species_copy_range(app, &app->species[i], VM_RK_F, VM_RK_F1);
            }
            // no-ops for the null field
            vlasov_field_combine(app, app->field->em1,
              1.0/3.0, app->field->em, 2.0/3.0, app->field->emnew);
            vlasov_field_copy_range(app, app->field->em, app->field->em1);

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
