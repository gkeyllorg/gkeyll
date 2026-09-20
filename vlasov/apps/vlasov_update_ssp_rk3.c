#include <gkyl_vlasov_priv.h>

// Set the per-stage RK input/output array pointers for every species, indexed
// over the overall species count. 'in'/'out' select which of the kinetic RK
// buffers (f/f1/fnew) and fluid RK buffers (fluid/fluid1/fluidnew) this stage
// reads and writes; entries are NULL where a species lacks that aspect.
enum vm_rk_buf { VM_RK_F, VM_RK_F1, VM_RK_FNEW };

static struct gkyl_array*
dist_buf(struct vm_species *d, enum vm_rk_buf b)
{
  return b == VM_RK_F ? d->f : (b == VM_RK_F1 ? d->f1 : d->fnew);
}
static struct gkyl_array*
fluid_buf(struct vm_fluid_species *f, enum vm_rk_buf b)
{
  return b == VM_RK_F ? f->fluid : (b == VM_RK_F1 ? f->fluid1 : f->fluidnew);
}

static void
set_rk_arrays(gkyl_vlasov_app *app, enum vm_rk_buf in, enum vm_rk_buf out,
  const struct gkyl_array *fin[], struct gkyl_array *fout[],
  const struct gkyl_array *fluidin[], struct gkyl_array *fluidout[])
{
  int num_species = app->num_species;
  for (int i=0; i<num_species; ++i) {
    struct vlasov_species *sp = &app->species[i];
    fin[i]      = sp->kinetic  ? dist_buf(sp->kinetic, in)   : 0;
    fout[i]     = sp->kinetic  ? dist_buf(sp->kinetic, out)  : 0;
    fluidin[i]  = sp->fluid ? fluid_buf(sp->fluid, in)  : 0;
    fluidout[i] = sp->fluid ? fluid_buf(sp->fluid, out) : 0;
  }
}

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

          set_rk_arrays(app, VM_RK_F, VM_RK_F1, fin, fout, fluidin, fluidout);
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

          set_rk_arrays(app, VM_RK_F1, VM_RK_FNEW, fin, fout, fluidin, fluidout);
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
              vlasov_species_combine(app, &app->species[i], 3.0/4.0, 1.0/4.0);
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

          set_rk_arrays(app, VM_RK_F1, VM_RK_FNEW, fin, fout, fluidin, fluidout);
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
              vlasov_species_combine(app, &app->species[i], 1.0/3.0, 2.0/3.0);
              vlasov_species_copy_range(app, &app->species[i]);
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
