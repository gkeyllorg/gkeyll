#include <gkyl_vlasov_priv.h>

// Take time-step using the RK3 method. Also sets the status object
// which has the actual and suggested dts used. These can be different
// from the actual time-step.
struct gkyl_update_status
vlasov_poisson_update_ssp_rk3(gkyl_vlasov_app* app, double dt0)
{
  int ns = app->num_species;
  int nfs = app->num_fluid_species;

  const struct gkyl_array *fin[ns];
  const struct gkyl_array *fluidin[nfs];
  struct gkyl_array *fout[ns];
  struct gkyl_array *fluidout[nfs];
  struct gkyl_update_status st = { .success = true };

  // Time-stepper state.
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        do {
          struct timespec rk3_s1_tm = gkyl_wall_clock();

          for (int i=0; i<ns; ++i) {
            fin[i] = app->species[i].f;
            fout[i] = app->species[i].f1;
          }
  
          vlasov_forward_euler(app, tcurr, dt, fin, fluidin, 0, fout, fluidout, 0, &st);

          // Compute the fields and apply BCs.
          vp_calc_field_and_apply_bc(app, tcurr, fout);
  
          dt = st.dt_actual;
          state = RK_STAGE_2;

          app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s1_tm);
        } while (0);
        break;

      case RK_STAGE_2:
        do {
          struct timespec rk3_s2_tm = gkyl_wall_clock();

          for (int i=0; i<ns; ++i) {
            fin[i] = app->species[i].f1;
            fout[i] = app->species[i].fnew;
          }
  
          vlasov_forward_euler(app, tcurr+dt, dt, fin, fluidin, 0, fout, fluidout, 0, &st);
  
          if (st.dt_actual < dt) {
  
            // Recalculate the field.
            for (int i=0; i<ns; ++i)
              fin[i] = app->species[i].f;
            vp_calc_field(app, tcurr, fin);
  
            // collect stats
            double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
            app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0],
              dt_rel_diff);
            app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1],
              dt_rel_diff);
            app->stat.nstage_2_fail += 1;
  
            dt = st.dt_actual;
            state = RK_STAGE_1; // restart from stage 1
  
          } else {
            for (int i=0; i<ns; ++i)
              array_combine(app->species[i].f1,
                3.0/4.0, app->species[i].f, 1.0/4.0, app->species[i].fnew, &app->species[i].local_ext);
  
            // Compute the fields and apply BCs.
            for (int i=0; i<ns; ++i) {
              fout[i] = app->species[i].f1;
            }
            vp_calc_field_and_apply_bc(app, tcurr, fout);
  
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
            fout[i] = app->species[i].fnew;
          }
  
          vlasov_forward_euler(app, tcurr+dt/2, dt, fin, fluidin, 0, fout, fluidout, 0, &st);
  
          if (st.dt_actual < dt) {
            // Recalculate the field.
            for (int i=0; i<ns; ++i)
              fin[i] = app->species[i].f;
            vp_calc_field(app, tcurr, fin);
  
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
  
            // Compute the fields and apply BCs
            for (int i=0; i<ns; ++i) {
              fout[i] = app->species[i].f;
            }
            vp_calc_field_and_apply_bc(app, tcurr, fout);
  
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
