#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>

void
gk_field_biased_wall_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  f->phi_wall_lo = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->has_phi_wall_lo = false;
  f->phi_wall_lo_evolve = false;
  if (f->info.phi_wall_lo)
  {
    f->has_phi_wall_lo = true;
    if (f->info.phi_wall_lo_evolve) {
      f->phi_wall_lo_evolve = f->info.phi_wall_lo_evolve;
    }

    f->phi_wall_lo_host = f->phi_wall_lo;
    if (app->use_gpu) {
      f->phi_wall_lo_host = mkarr(false, f->phi_wall_lo->ncomp, f->phi_wall_lo->size);
    }

    f->phi_wall_lo_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis,
                                                 1, f->info.phi_wall_lo, f->info.phi_wall_lo_ctx);

    // Compute phi_wall_lo at t = 0
    gkyl_eval_on_nodes_advance(f->phi_wall_lo_proj, 0.0, &app->local_ext, f->phi_wall_lo_host);
    if (app->use_gpu) { // note: phi_wall_lo_host is same as phi_wall_lo when not on GPUs
      gkyl_array_copy(f->phi_wall_lo, f->phi_wall_lo_host);
    }
  }

  // Set up biased upper wall (same size as electrostatic potential), by default is 0.0
  f->phi_wall_up = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->has_phi_wall_up = false;
  f->phi_wall_up_evolve = false;
  if (f->info.phi_wall_up)
  {
    f->has_phi_wall_up = true;
    if (f->info.phi_wall_up_evolve) {
      f->phi_wall_up_evolve = f->info.phi_wall_up_evolve;
    }

    f->phi_wall_up_host = f->phi_wall_up;
    if (app->use_gpu) {
      f->phi_wall_up_host = mkarr(false, f->phi_wall_up->ncomp, f->phi_wall_up->size);
    }

    f->phi_wall_up_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis,
                                                 1, f->info.phi_wall_up, f->info.phi_wall_up_ctx);

    // Compute phi_wall_up at t = 0.
    gkyl_eval_on_nodes_advance(f->phi_wall_up_proj, 0.0, &app->local_ext, f->phi_wall_up_host);
    if (app->use_gpu) { // Note: phi_wall_up_host is same as phi_wall_up when not on GPUs.
      gkyl_array_copy(f->phi_wall_up, f->phi_wall_up_host);
    }
  }
}

void
gk_field_calc_phi_wall(gkyl_gyrokinetic_app *app, struct gk_field *field, double tm)
{
  if (field->has_phi_wall_lo && field->phi_wall_lo_evolve) {
    gkyl_eval_on_nodes_advance(field->phi_wall_lo_proj, tm, &app->local_ext, field->phi_wall_lo_host);
    if (app->use_gpu) {
      gkyl_array_copy(field->phi_wall_lo, field->phi_wall_lo_host);
    }
  }
  if (field->has_phi_wall_up && field->phi_wall_up_evolve) {
    gkyl_eval_on_nodes_advance(field->phi_wall_up_proj, tm, &app->local_ext, field->phi_wall_up_host);
    if (app->use_gpu) {
      gkyl_array_copy(field->phi_wall_up, field->phi_wall_up_host);
    }
  }
}

void
gk_field_biased_wall_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_array_release(f->phi_wall_lo);
  if (f->has_phi_wall_lo) {
    gkyl_eval_on_nodes_release(f->phi_wall_lo_proj);
    if (app->use_gpu) {
      gkyl_array_release(f->phi_wall_lo_host);
    }
  }

  gkyl_array_release(f->phi_wall_up);
  if (f->has_phi_wall_up) { 
    gkyl_eval_on_nodes_release(f->phi_wall_up_proj);
    if (app->use_gpu) {
      gkyl_array_release(f->phi_wall_up_host);
    }
  }
}
