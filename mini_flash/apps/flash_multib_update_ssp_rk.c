#include <gkyl_app_priv.h>
#include <gkyl_flash_multib_priv.h>

// Arrays passed between the multib time-stepper routines hold one entry per
// (local block, species) pair, indexed as b*num_species + i.

// Sync interpolated fluid momentum density at interblock boundaries.
static void
flash_multib_sync_fluid_rhou(struct gkyl_flash_multib_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    struct gkyl_array *rhou_blocks[app->num_local_blocks];
    for (int b=0; b<app->num_local_blocks; ++b)
      rhou_blocks[b] = app->block_apps[b]->species[i].fluid_rhou;

    gkyl_multib_comm_conn_array_transfer(app->comm, app->num_local_blocks, app->local_blocks,
      app->mbcc_sync_conf->send, app->mbcc_sync_conf->recv, rhou_blocks, rhou_blocks);
  }
}

// Take a single forward Euler step of the hybrid Euler-PKPM system on all
// blocks with suggested timestep dt0. The step taken is the same on every
// block and rank, is never larger than dt0, and is stored in th status object.
static void
flash_multib_forward_euler(struct gkyl_flash_multib_app *app, double tcurr, double dt0,
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[],
  struct gkyl_array *fout[], struct gkyl_array *fluidout[],
  struct gkyl_update_status *st)
{
  app->stat.nfeuler += 1;
  struct timespec wst = gkyl_wall_clock();

  int ns = app->num_species;

  // Compute field variables and interpolate bulk velocity onto kinetic grid.
  // Note that we interpolate the velocity, but the kinetic update looks for 
  // *momentum density*. Therefore, we need to compute the density moment first
  // to scale the interpolated velocity.
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_flash_app *sbapp = app->block_apps[b];
    flash_field_calc_bvar(sbapp, sbapp->field, sbapp->field->em);

    for (int i=0; i<ns; ++i) {
      flash_species_calc_pkpm_moms(sbapp, &sbapp->species[i], fin[b*ns+i]);
      flash_species_interp_rhou_to_dg(sbapp, &sbapp->species[i], fluidin[b*ns+i]);
    }
  }
  // Sync fluid momentum density across connected blocks.
  flash_multib_sync_fluid_rhou(app);

  // Compute moments for the PKPM update and coupling fluxes for the kinetic equation,
  // then compute the RHS of the kinetic update.
  double dt_kin = DBL_MAX;
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_flash_app *sbapp = app->block_apps[b];
    for (int i=0; i<ns; ++i) {
      struct flash_species *s = &sbapp->species[i];
      flash_species_calc_pkpm_vars(sbapp, s, fin[b*ns+i]);
      flash_species_calc_pkpm_update_vars(sbapp, s, fin[b*ns+i]);
      flash_species_calc_coupling_fluxes(sbapp, s, fin[b*ns+i]);
      dt_kin = fmin(dt_kin, flash_species_rhs(sbapp, s, fin[b*ns+i], fout[b*ns+i]));
    }
  }

  // Advance the Euler equations with the suggested timestep,
  // storing fluid timestep to determine if CFL was violated.
  double dt_flu = DBL_MAX;
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_flash_app *sbapp = app->block_apps[b];
    for (int i=0; i<ns; ++i) {
      dt_flu = fmin(dt_flu, flash_fluid_species_update(sbapp, tcurr, dt0,
        &sbapp->species[i], &sbapp->fluid_species[i], fluidin[b*ns+i], fluidout[b*ns+i]));
    }
  }

  // Determine limiting timesteps across all blocks and ranks.
  double dt_local[2] = { dt_kin, dt_flu }, dt_global[2];
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 2, dt_local, dt_global);
  dt_kin = dt_global[0];
  dt_flu = dt_global[1];

  double dtmin = fmin(dt_kin, dt_flu);
  // If dtmin is only slightly smaller than dt0, use dt0 (avoids retaking steps
  // for very small dt changes). Not applied if the fluid update failed.
  double dt_max_rel_diff = 0.01;
  double dt_rel_diff = (dt0-dtmin)/dt0;
  if (dt_flu >= dt0 && dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt0;

  // Never take a timestep larger than dt0.
  double dta = st->dt_actual = fmin(dt0, dtmin);
  st->dt_suggested = dtmin;

  // If wave propagation update failed earlier due to violating CFL, redo with the proper timestep.
  // Redo the fluid update with the timestep actually taken.
  if (dta < dt0) {
    for (int b=0; b<app->num_local_blocks; ++b) {
      struct gkyl_flash_app *sbapp = app->block_apps[b];
      for (int i=0; i<ns; ++i) {
        flash_fluid_species_update(sbapp, tcurr, dta, &sbapp->species[i], &sbapp->fluid_species[i],
          fluidin[b*ns+i], fluidout[b*ns+i]);
      }
    }
  }

  // Complete the forward Euler step of the distribution functions and 
  // apply BCs to distribution functions and Euler variables.
  for (int b=0; b<app->num_local_blocks; ++b) {
    for (int i=0; i<ns; ++i)
      gkyl_array_accumulate(gkyl_array_scale(fout[b*ns+i], dta), 1.0, fin[b*ns+i]);
  }

  // Physical and inter-block boundary conditions.
  gkyl_flash_multib_app_apply_bc(app, tcurr+dta, fout, fluidout);

  app->stat.fwd_euler_tm += gkyl_time_diff_now_sec(wst);
}

// Take a single time step with the 3-stage SSP-RK method with suggested
// timestep dt0. This function will never take a timestep larger than dt0,
// but may take one smaller. The actual dt is returned in the status struct.
struct gkyl_update_status
gkyl_flash_multib_update_ssp_rk3(struct gkyl_flash_multib_app *app, double dt0)
{
  int ns = app->num_species;
  int nsb = ns*app->num_local_blocks;

  const struct gkyl_array *fin[nsb], *fluidin[nsb];
  struct gkyl_array *fout[nsb], *fluidout[nsb];
  struct gkyl_update_status st = { .success = true };

  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;

  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        for (int b=0; b<app->num_local_blocks; ++b) {
          struct gkyl_flash_app *sbapp = app->block_apps[b];
          int li = b*ns;
          for (int i=0; i<ns; ++i) {
            fin[li+i] = sbapp->species[i].f;
            fout[li+i] = sbapp->species[i].f1;
            fluidin[li+i] = sbapp->fluid_species[i].fluid;
            fluidout[li+i] = sbapp->fluid_species[i].fluid1;
          }
        }
        flash_multib_forward_euler(app, tcurr, dt, fin, fluidin, fout, fluidout, &st);
        dt = st.dt_actual;

        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int b=0; b<app->num_local_blocks; ++b) {
          struct gkyl_flash_app *sbapp = app->block_apps[b];
          int li = b*ns;
          for (int i=0; i<ns; ++i) {
            fin[li+i] = sbapp->species[i].f1;
            fout[li+i] = sbapp->species[i].fnew;
            fluidin[li+i] = sbapp->fluid_species[i].fluid1;
            fluidout[li+i] = sbapp->fluid_species[i].fluidnew;
          }
        }
        flash_multib_forward_euler(app, tcurr+dt, dt, fin, fluidin, fout, fluidout, &st);

        if (st.dt_actual < dt) {
          // Timestep inconsistent with stage 1: restart with the smaller one.
          app->stat.nfail += 1;
          dt = st.dt_actual;
          state = RK_STAGE_1;
        }
        else {
          for (int b=0; b<app->num_local_blocks; ++b) {
            struct gkyl_flash_app *sbapp = app->block_apps[b];
            for (int i=0; i<ns; ++i) {
              struct flash_species *s = &sbapp->species[i];
              struct flash_fluid_species *fs = &sbapp->fluid_species[i];
              array_combine(s->f1, 3.0/4.0, s->f, 1.0/4.0, s->fnew, &s->local_ext);
              array_combine(fs->fluid1, 3.0/4.0, fs->fluid, 1.0/4.0, fs->fluidnew, &fs->local_ext);
            }
          }
          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int b=0; b<app->num_local_blocks; ++b) {
          struct gkyl_flash_app *sbapp = app->block_apps[b];
          int li = b*ns;
          for (int i=0; i<ns; ++i) {
            fin[li+i] = sbapp->species[i].f1;
            fout[li+i] = sbapp->species[i].fnew;
            fluidin[li+i] = sbapp->fluid_species[i].fluid1;
            fluidout[li+i] = sbapp->fluid_species[i].fluidnew;
          }
        }
        flash_multib_forward_euler(app, tcurr+dt/2, dt, fin, fluidin, fout, fluidout, &st);

        if (st.dt_actual < dt) {
          // Timestep inconsistent with earlier stages: restart with the smaller one.
          app->stat.nfail += 1;
          dt = st.dt_actual;
          state = RK_STAGE_1;
        }
        else {
          for (int b=0; b<app->num_local_blocks; ++b) {
            struct gkyl_flash_app *sbapp = app->block_apps[b];
            for (int i=0; i<ns; ++i) {
              struct flash_species *s = &sbapp->species[i];
              struct flash_fluid_species *fs = &sbapp->fluid_species[i];
              array_combine(s->f1, 1.0/3.0, s->f, 2.0/3.0, s->fnew, &s->local_ext);
              gkyl_array_copy_range(s->f, s->f1, &s->local_ext);
              array_combine(fs->fluid1, 1.0/3.0, fs->fluid, 2.0/3.0, fs->fluidnew, &fs->local_ext);
              gkyl_array_copy_range(fs->fluid, fs->fluid1, &fs->local_ext);
            }
          }
          state = RK_COMPLETE;
        }
        break;

      case RK_COMPLETE:
        break;
    }
  }

  return st;
}
