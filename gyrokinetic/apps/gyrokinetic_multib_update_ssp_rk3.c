#include <gkyl_gyrokinetic_multib_priv.h>

static void
gyrokinetic_multib_forward_euler(
  struct gkyl_gyrokinetic_multib_app *app, double tcurr, double dt,
  struct gkyl_gyrokinetic_multib_fdot_args *fdot_args, struct gkyl_update_status *st
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

  double dtmin = DBL_MAX;

  // Compute the time rate of change of the distributions, df/dt.
  for (int b = 0; b < app->num_local_blocks; ++b) {
    int li_charged = b * app->num_species;
    int li_neut = b * app->num_neut_species;
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    struct gkyl_gyrokinetic_fdot_args *sb_args = &sbapp->fdot_args;
    for (int i = 0; i < app->num_species; ++i) {
      sb_args->fin[i] = fin[li_charged + i];
      sb_args->fout[i] = fout[li_charged + i];
      sb_args->smom_in[i] = smom_in[li_charged + i];
      sb_args->smom_out[i] = smom_out[li_charged + i];
    }
    for (int i = 0; i < app->num_neut_species; ++i) {
      sb_args->fin_neut[i] = fin_neut[li_neut + i];
      sb_args->fout_neut[i] = fout_neut[li_neut + i];
      sb_args->smom_in_neut[i] = smom_in_neut[li_neut + i];
      sb_args->smom_out_neut[i] = smom_out_neut[li_neut + i];
    }
    // A_parallel is not evolved in multiblock.
    sb_args->aparin = NULL;
    sb_args->aparout = NULL;
    gyrokinetic_rhs(sbapp, tcurr, dt, sb_args, st);
    dtmin = fmin(dtmin, st->dt_actual);
  }

  struct timespec wtm = gkyl_wall_clock();
  // Compute minimum time-step across all processors.
  double dtmin_local = dtmin, dtmin_global;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &dtmin_local, &dtmin_global);
  st->dt_actual = dtmin_global;
  app->stat.dfdt_dt_reduce_tm += gkyl_time_diff_now_sec(wtm);

  struct timespec wst = gkyl_wall_clock();
  // Complete update of distribution functions.
  double dta = st->dt_actual;
  for (int b = 0; b < app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int li_charged = b * app->num_species;
    int li_neut = b * app->num_neut_species;
    for (int i = 0; i < app->num_species; ++i) {
      struct gk_species *gks = &sbapp->species[i];
      gk_species_step_f(gks, fout[li_charged + i], dta, fin[li_charged + i]);
      gk_stepped_moms_step_f(
        &gks->step_moms, smom_out[li_charged + i], 1.0, smom_in[li_charged + i]
      );
    }
    for (int i = 0; i < app->num_neut_species; ++i) {
      struct gk_neut_species *gkns = &sbapp->neut_species[i];
      gk_neut_species_step_f(gkns, fout_neut[li_neut + i], dta, fin_neut[li_neut + i]);
      gk_stepped_moms_step_f(
        &gkns->step_moms, smom_out_neut[li_neut + i], 1.0, smom_in_neut[li_neut + i]
      );
    }
  }

  app->stat.fwd_euler_step_f_tm += gkyl_time_diff_now_sec(wst);
  app->stat.fwd_euler_tm += gkyl_time_diff_now_sec(wst_fe);
}

struct gkyl_update_status
gyrokinetic_multib_update_ssp_rk3(struct gkyl_gyrokinetic_multib_app *app, double dt0)
{
  // Take time-step using the RK3 method. Also sets the status object
  // which has the actual and suggested dts used. These can be different
  // from the actual time-step.
  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;
  int nblocks_local = app->num_local_blocks;

  struct gkyl_gyrokinetic_multib_fdot_args *fdot_args = &app->fdot_args;
  const struct gkyl_array **fin = fdot_args->fin;
  const struct gkyl_array **fin_neut = fdot_args->fin_neut;
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***smom_in = fdot_args->smom_in;
  struct gkyl_array ***smom_in_neut = fdot_args->smom_in_neut;
  struct gkyl_array ***smom_out = fdot_args->smom_out;
  struct gkyl_array ***smom_out_neut = fdot_args->smom_out_neut;

  struct gkyl_update_status st = {.success = true};

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        for (int b = 0; b < nblocks_local; ++b) {
          struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
          int li_charged = b * ns_charged;
          int li_neut = b * ns_neut;
          for (int i = 0; i < ns_charged; ++i) {
            struct gk_species *gks = &sbapp->species[i];
            fin[li_charged + i] = gks->f;
            fout[li_charged + i] = gks->f1;
            // Stepped moments.
            smom_in[li_charged + i] = gks->step_moms.f;
            smom_out[li_charged + i] = gks->step_moms.f1;
          }
          for (int i = 0; i < ns_neut; ++i) {
            struct gk_neut_species *gkns = &sbapp->neut_species[i];
            fin_neut[li_neut + i] = gkns->f;
            fout_neut[li_neut + i] = gkns->f1;
            // Stepped moments.
            smom_in_neut[li_neut + i] = gkns->step_moms.f;
            smom_out_neut[li_neut + i] = gkns->step_moms.f1;
          }
        }

        gyrokinetic_multib_forward_euler(app, tcurr, dt, fdot_args, &st);
        dt = st.dt_actual;

        // Subtract the old moments from f1 to accumulate only the current
        // timestep's contribution. Keep the old moments in f for RK-stage retries.
        struct timespec wst = gkyl_wall_clock();
        for (int b = 0; b < nblocks_local; ++b) {
          struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
          int li_charged = b * ns_charged;
          int li_neut = b * ns_neut;
          for (int i = 0; i < ns_charged; ++i) {
            struct gk_species *gks = &sbapp->species[i];
            gk_stepped_moms_accumulate(
              &gks->step_moms, smom_out[li_charged + i], -1.0, smom_in[li_charged + i]
            );
          }
          for (int i = 0; i < ns_neut; ++i) {
            struct gk_neut_species *gkns = &sbapp->neut_species[i];
            gk_stepped_moms_accumulate(
              &gkns->step_moms, smom_out_neut[li_neut + i], -1.0, smom_in_neut[li_neut + i]
            );
          }
        }
        app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

        for (int b = 0; b < nblocks_local; ++b) {
          struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
          for (int i = 0; i < ns_charged; ++i) {
            struct gk_species *gks = &sbapp->species[i];
            // Compute moment of f_old to later compute moment of df/dt.
            // Do it before the fields are updated, but after dt is calculated.
            gk_species_calc_integ_mom_dt(sbapp, gks, dt, gks->fdot_mom_old);
          }

          // Store the field energies divided by dt for the energy rate and balance diagnostics.
          gk_field_calc_energy_dt(sbapp, sbapp->field, dt, true);
        }

        // Compute the fields and apply BCs.
        gyrokinetic_multib_calc_field_and_apply_bc(app, tcurr, fout, smom_out, fout_neut);

        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int b = 0; b < nblocks_local; ++b) {
          struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
          int li_charged = b * ns_charged;
          int li_neut = b * ns_neut;
          for (int i = 0; i < ns_charged; ++i) {
            struct gk_species *gks = &sbapp->species[i];
            fin[li_charged + i] = gks->f1;
            fout[li_charged + i] = gks->fnew;
            // Stepped moments.
            smom_in[li_charged + i] = gks->step_moms.f1;
            smom_out[li_charged + i] = gks->step_moms.fnew;
          }
          for (int i = 0; i < ns_neut; ++i) {
            struct gk_neut_species *gkns = &sbapp->neut_species[i];
            fin_neut[li_neut + i] = gkns->f1;
            fout_neut[li_neut + i] = gkns->fnew;
            // Stepped moments.
            smom_in_neut[li_neut + i] = gkns->step_moms.f1;
            smom_out_neut[li_neut + i] = gkns->step_moms.fnew;
          }
        }

        gyrokinetic_multib_forward_euler(app, tcurr + dt, dt, fdot_args, &st);

        if (st.dt_actual < dt) {
          // Recalculate the field.
          for (int b = 0; b < nblocks_local; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            int li_charged = b * ns_charged;
            for (int i = 0; i < ns_charged; ++i) {
              struct gk_species *gks = &sbapp->species[i];
              fin[li_charged + i] = gks->f;
              smom_in[li_charged + i] = gks->step_moms.f;
            }
          }
          gyrokinetic_multib_calc_field(app, tcurr, fin, smom_in);

          // Collect stats.
          double dt_rel_diff = (dt - st.dt_actual) / st.dt_actual;
          app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0], dt_rel_diff);
          app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1], dt_rel_diff);
          app->stat.nstage_2_fail += 1;

          dt = st.dt_actual;
          state = RK_STAGE_1; // Restart from stage 1.

        } else {
          struct timespec wst = gkyl_wall_clock();
          for (int b = 0; b < nblocks_local; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            for (int i = 0; i < ns_charged; ++i) {
              struct gk_species *gks = &sbapp->species[i];
              gk_species_combine(
                gks, gks->f1, 3.0 / 4.0, gks->f, 1.0 / 4.0, gks->fnew, &gks->local_ext
              );
              gk_stepped_moms_set(
                &gks->step_moms, gks->step_moms.f1, 1.0 / 4.0, gks->step_moms.fnew
              );
            }
            for (int i = 0; i < ns_neut; ++i) {
              struct gk_neut_species *gkns = &sbapp->neut_species[i];
              gk_neut_species_combine(
                gkns, gkns->f1, 3.0 / 4.0, gkns->f, 1.0 / 4.0, gkns->fnew, &gkns->local_ext
              );
              gk_stepped_moms_set(
                &gkns->step_moms, gkns->step_moms.f1, 1.0 / 4.0, gkns->step_moms.fnew
              );
            }
          }
          app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

          for (int b = 0; b < nblocks_local; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            int li_charged = b * ns_charged;
            int li_neut = b * ns_neut;
            // Compute the fields and apply BCs.
            for (int i = 0; i < ns_charged; ++i) {
              struct gk_species *gks = &sbapp->species[i];
              fout[li_charged + i] = gks->f1;
              smom_out[li_charged + i] = gks->step_moms.f1;
            }
            for (int i = 0; i < ns_neut; ++i) {
              fout_neut[li_neut + i] = sbapp->neut_species[i].f1;
            }
          }
          gyrokinetic_multib_calc_field_and_apply_bc(app, tcurr, fout, smom_out, fout_neut);

          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int b = 0; b < nblocks_local; ++b) {
          struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
          int li_charged = b * ns_charged;
          int li_neut = b * ns_neut;
          for (int i = 0; i < ns_charged; ++i) {
            struct gk_species *gks = &sbapp->species[i];
            fin[li_charged + i] = gks->f1;
            fout[li_charged + i] = gks->fnew;
            // Stepped moments.
            smom_in[li_charged + i] = gks->step_moms.f1;
            smom_out[li_charged + i] = gks->step_moms.fnew;
          }
          for (int i = 0; i < ns_neut; ++i) {
            struct gk_neut_species *gkns = &sbapp->neut_species[i];
            fin_neut[li_neut + i] = sbapp->neut_species[i].f1;
            fout_neut[li_neut + i] = sbapp->neut_species[i].fnew;
            // Stepped moments.
            smom_in_neut[li_neut + i] = gkns->step_moms.f1;
            smom_out_neut[li_neut + i] = gkns->step_moms.fnew;
          }
        }

        gyrokinetic_multib_forward_euler(app, tcurr + dt / 2, dt, fdot_args, &st);

        if (st.dt_actual < dt) {
          // Recalculate the field.
          for (int b = 0; b < nblocks_local; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            int li_charged = b * ns_charged;
            for (int i = 0; i < ns_charged; ++i) {
              struct gk_species *gks = &sbapp->species[i];
              fin[li_charged + i] = gks->f;
              smom_in[li_charged + i] = gks->step_moms.f;
            }
          }
          gyrokinetic_multib_calc_field(app, tcurr, fin, smom_in);

          // Collect stats.
          double dt_rel_diff = (dt - st.dt_actual) / st.dt_actual;
          app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0], dt_rel_diff);
          app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1], dt_rel_diff);
          app->stat.nstage_3_fail += 1;

          dt = st.dt_actual;
          state = RK_STAGE_1; // Restart from stage 1.

          app->stat.nstage_2_fail += 1;
        } else {
          struct timespec wst = gkyl_wall_clock();
          for (int b = 0; b < nblocks_local; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            for (int i = 0; i < ns_charged; ++i) {
              struct gk_species *gks = &sbapp->species[i];
              // Step f.
              gk_species_combine(
                gks, gks->f1, 1.0 / 3.0, gks->f, 2.0 / 3.0, gks->fnew, &gks->local_ext
              );
              gk_species_copy_range(gks, gks->f, gks->f1, &gks->local_ext);
              // Update stepped moments.
              gk_stepped_moms_set(&gks->step_moms, gks->step_moms.f, 2.0 / 3.0, gks->step_moms.fnew);
            }
            for (int i = 0; i < ns_neut; ++i) {
              struct gk_neut_species *gkns = &sbapp->neut_species[i];
              gk_neut_species_combine(
                gkns, gkns->f1, 1.0 / 3.0, gkns->f, 2.0 / 3.0, gkns->fnew, &gkns->local_ext
              );
              gk_neut_species_copy_range(gkns, gkns->f, gkns->f1, &gkns->local_ext);
              // Update stepped moments.
              gk_stepped_moms_set(
                &gkns->step_moms, gkns->step_moms.f, 2.0 / 3.0, gkns->step_moms.fnew
              );
            }
          }
          app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

          // Apply positivity shift if requested.
          for (int b = 0; b < nblocks_local; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            for (int i = 0; i < ns_charged; ++i) {
              struct gk_species *gks = &sbapp->species[i];
              gk_species_positivity_apply(sbapp, gks, &gks->positivity, gks->fnew, gks->f);
            }
          }

          for (int b = 0; b < nblocks_local; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            int li_charged = b * ns_charged;
            int li_neut = b * ns_neut;
            for (int i = 0; i < ns_charged; ++i) {
              struct gk_species *gks = &sbapp->species[i];
              fout[li_charged + i] = gks->f;
              smom_out[li_charged + i] = gks->step_moms.f;
            }
            for (int i = 0; i < ns_neut; ++i) {
              fout_neut[li_neut + i] = sbapp->neut_species[i].f;
            }
          }

          //          // Enforce quasineutrality of the positivity shifts.
          //          gyrokinetic_multib_post_positivity_quasineut(app, fout);

          // Compute the fields and apply BCs
          gyrokinetic_multib_calc_field_and_apply_bc(app, tcurr, fout, smom_out, fout_neut);

          for (int b = 0; b < nblocks_local; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            for (int i = 0; i < ns_charged; ++i) {
              struct gk_species *gks = &sbapp->species[i];
              // Compute moment of f_new to compute moment of df/dt.
              // Need to do it after the fields are updated.
              gk_species_calc_integ_mom_dt(sbapp, gks, dt, gks->fdot_mom_new);
            }

            // Store the field energies divided by dt for the energy rate and balance diagnostics.
            gk_field_calc_energy_dt(sbapp, sbapp->field, dt, false);
          }

          state = RK_COMPLETE;
        }
        break;

      case RK_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }

  return st;
}
