#include <gkyl_moment_priv.h>

// internal function that takes a single time-step using a single-step
// Strang-split scheme
struct gkyl_update_status
moment_update_one_step(gkyl_moment_app* app, double dt0)
{
  int ns = app->num_species, ndim = app->ndim;
  bool have_nans_occured = false;
  
  double dt_suggested = DBL_MAX;
  
  // time-stepper states
  enum {
    UPDATE_DONE = 0,
    PRE_UPDATE,
    POST_UPDATE,
    FIRST_COUPLING_UPDATE,
    FIELD_UPDATE,
    SPECIES_UPDATE,
    SECOND_COUPLING_UPDATE,
    UPDATE_REDO,
  } state = PRE_UPDATE;

  double tcurr = app->tcurr, dt = dt0;
  while (state != UPDATE_DONE) {
    switch (state) {
      case PRE_UPDATE:
        state = FIRST_COUPLING_UPDATE; // next state
          
        // copy old solution in case we need to redo this step
        for (int i=0; i<ns; ++i)
          gkyl_array_copy(app->species[i].fdup, app->species[i].f[0]);
        if (app->has_field)
          gkyl_array_copy(app->field.fdup, app->field.f[0]);

        break;
          
      case FIRST_COUPLING_UPDATE:
        state = FIELD_UPDATE; // next state

        if (app->update_sources) {
          struct timespec src1_tm = gkyl_wall_clock();
          struct gkyl_update_status s = moment_coupling_update(app, &app->sources,
            0, tcurr, dt/2);
          if (!s.success) {
            app->stat.nfail += 1;
            dt = s.dt_suggested;
            state = UPDATE_REDO;
            break;
          }
          dt_suggested = fmin(dt_suggested, s.dt_suggested);
          app->stat.sources_tm += gkyl_time_diff_now_sec(src1_tm);
        }
        if (app->update_mhd_source) {
          struct timespec src1_tm = gkyl_wall_clock();
          mhd_src_update(app, &app->mhd_source, 0, tcurr, dt/2);
          app->stat.sources_tm += gkyl_time_diff_now_sec(src1_tm);
        }

        break;

      case FIELD_UPDATE:
        state = SPECIES_UPDATE; // next state

        if (app->has_field) {
          struct timespec fl_tm = gkyl_wall_clock();
          struct gkyl_update_status s = moment_field_update(app, &app->field, tcurr, dt);
          if (!s.success) {
            app->stat.nfail += 1;
            dt = s.dt_suggested;
            state = UPDATE_REDO;
            break;
          }
            
          dt_suggested = fmin(dt_suggested, s.dt_suggested);
          app->stat.field_tm += gkyl_time_diff_now_sec(fl_tm);
        }
          
        break;

      case SPECIES_UPDATE:
        state = SECOND_COUPLING_UPDATE; // next state

        struct timespec sp_tm = gkyl_wall_clock();
        for (int i=0; i<ns; ++i) {         
          struct gkyl_update_status s =
            moment_species_update(app, &app->species[i], tcurr, dt);

          if (!s.success) {
            app->stat.nfail += 1;
            dt = s.dt_suggested;
            state = UPDATE_REDO;
            break;
          }
          dt_suggested = fmin(dt_suggested, s.dt_suggested);
        }
        app->stat.species_tm += gkyl_time_diff_now_sec(sp_tm);
         
        break;

      case SECOND_COUPLING_UPDATE:
        state = POST_UPDATE; // next state

        if (app->update_sources) {
          struct timespec src2_tm = gkyl_wall_clock();
          moment_coupling_update(app, &app->sources, 1, tcurr, dt/2);
          app->stat.sources_tm += gkyl_time_diff_now_sec(src2_tm);
        }
        if (app->update_mhd_source) {
          struct timespec src2_tm = gkyl_wall_clock();
          mhd_src_update(app, &app->mhd_source, 1, tcurr, dt/2);
          app->stat.sources_tm += gkyl_time_diff_now_sec(src2_tm);
        }

        break;

      case POST_UPDATE:
        state = UPDATE_DONE;

        // copy solution in prep for next time-step
        for (int i=0; i<ns; ++i) {
          // check for nans before copying
          if (check_for_nans(app->species[i].f[ndim], app->local))
            have_nans_occured = true;
          else // only copy in case no nans, so old solution can be written out
            gkyl_array_copy(app->species[i].f[0], app->species[i].f[ndim]);
        }
        
        if (app->has_field)
          gkyl_array_copy(app->field.f[0], app->field.f[ndim]);
          
        break;

      case UPDATE_REDO:
        state = PRE_UPDATE; // start all-over again
          
        // restore solution and retake step
        for (int i=0; i<ns; ++i)
          gkyl_array_copy(app->species[i].f[0], app->species[i].fdup);
        if (app->has_field)
          gkyl_array_copy(app->field.f[0], app->field.fdup);
          
        break;

      case UPDATE_DONE: // unreachable code! (suppresses warning)
        break;
    }
  }

  return (struct gkyl_update_status) {
    .success = have_nans_occured ? false : true,
    .dt_actual = dt,
    .dt_suggested = dt_suggested,
  };
}
