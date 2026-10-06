#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>

void
gk_field_biased_wall_new(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gkf->phi_wall_lo = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkf->has_phi_wall_lo = false;
  gkf->phi_wall_lo_evolve = false;
  if (gkf->info.phi_wall_lo) {
    gkf->has_phi_wall_lo = true;
    if (gkf->info.phi_wall_lo_evolve) {
      gkf->phi_wall_lo_evolve = gkf->info.phi_wall_lo_evolve;
    }

    gkf->phi_wall_lo_host = gkf->phi_wall_lo;
    if (app->use_gpu) {
      gkf->phi_wall_lo_host = mkarr(false, gkf->phi_wall_lo->ncomp, gkf->phi_wall_lo->size);
    }

    gkf->phi_wall_lo_proj = gkyl_eval_on_nodes_new(
      &app->grid, &app->basis, 1, gkf->info.phi_wall_lo, gkf->info.phi_wall_lo_ctx
    );

    // Compute phi_wall_lo at t = 0
    gkyl_eval_on_nodes_advance(gkf->phi_wall_lo_proj, 0.0, &app->local_ext, gkf->phi_wall_lo_host);
    if (app->use_gpu) { // note: phi_wall_lo_host is same as phi_wall_lo when not on GPUs
      gkyl_array_copy(gkf->phi_wall_lo, gkf->phi_wall_lo_host);
    }
  }

  // Set up biased upper wall (same size as electrostatic potential), by default is 0.0
  gkf->phi_wall_up = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkf->has_phi_wall_up = false;
  gkf->phi_wall_up_evolve = false;
  if (gkf->info.phi_wall_up) {
    gkf->has_phi_wall_up = true;
    if (gkf->info.phi_wall_up_evolve) {
      gkf->phi_wall_up_evolve = gkf->info.phi_wall_up_evolve;
    }

    gkf->phi_wall_up_host = gkf->phi_wall_up;
    if (app->use_gpu) {
      gkf->phi_wall_up_host = mkarr(false, gkf->phi_wall_up->ncomp, gkf->phi_wall_up->size);
    }

    gkf->phi_wall_up_proj = gkyl_eval_on_nodes_new(
      &app->grid, &app->basis, 1, gkf->info.phi_wall_up, gkf->info.phi_wall_up_ctx
    );

    // Compute phi_wall_up at t = 0.
    gkyl_eval_on_nodes_advance(gkf->phi_wall_up_proj, 0.0, &app->local_ext, gkf->phi_wall_up_host);
    if (app->use_gpu) { // Note: phi_wall_up_host is same as phi_wall_up when not on GPUs.
      gkyl_array_copy(gkf->phi_wall_up, gkf->phi_wall_up_host);
    }
  }
}

void
gk_field_calc_phi_wall(gkyl_gyrokinetic_app *app, struct gk_field *field, double tm)
{
  if (field->has_phi_wall_lo && field->phi_wall_lo_evolve) {
    gkyl_eval_on_nodes_advance(
      field->phi_wall_lo_proj, tm, &app->local_ext, field->phi_wall_lo_host
    );
    if (app->use_gpu) {
      gkyl_array_copy(field->phi_wall_lo, field->phi_wall_lo_host);
    }
  }
  if (field->has_phi_wall_up && field->phi_wall_up_evolve) {
    gkyl_eval_on_nodes_advance(
      field->phi_wall_up_proj, tm, &app->local_ext, field->phi_wall_up_host
    );
    if (app->use_gpu) {
      gkyl_array_copy(field->phi_wall_up, field->phi_wall_up_host);
    }
  }
}

void
gk_field_biased_wall_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gkyl_array_release(gkf->phi_wall_lo);
  if (gkf->has_phi_wall_lo) {
    gkyl_eval_on_nodes_release(gkf->phi_wall_lo_proj);
    if (app->use_gpu) {
      gkyl_array_release(gkf->phi_wall_lo_host);
    }
  }

  gkyl_array_release(gkf->phi_wall_up);
  if (gkf->has_phi_wall_up) {
    gkyl_eval_on_nodes_release(gkf->phi_wall_up_proj);
    if (app->use_gpu) {
      gkyl_array_release(gkf->phi_wall_up_host);
    }
  }
}
