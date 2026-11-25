#include <gkyl_dg_bin_ops.h>
#include <gkyl_vlasov_priv.h>

// Held background of a fluid aspect (First_Fourier_Harmonic_Notes/background_density_gradient.tex):
// an Euler or isothermal Euler fluid species, or the fluid block of a PKPM species. Its force, the
// momentum rows of the right-hand side on the background divided by the density, joins the applied
// acceleration (rho a in the momentum and rho u.a in the energy, in every stage: a body force does
// work and scales with the density that evolves, which a fixed source would not); the remainder, the
// right-hand side with that force, is subtracted in every stage. The capture is in
// gkyl_vlasov_app_apply_ic (vlasov.c).

// Stage methods: nothing (no background, or before the remainder is stored), or its subtraction.
static void
vm_fluid_species_background_rhs_none(const struct vm_background *bg, struct gkyl_array *rhs)
{
}

static void
vm_fluid_species_background_rhs_remainder(const struct vm_background *bg, struct gkyl_array *rhs)
{
  gkyl_array_accumulate(rhs, -1.0, bg->rhs);
}

void
vm_fluid_species_background_init_none(struct vm_background *bg)
{
  *bg = (struct vm_background){.rhs_func = vm_fluid_species_background_rhs_none};
}

void
vm_fluid_species_background_init(
  struct gkyl_vlasov_app *app, struct vm_fluid_species *f, struct vm_background *bg
)
{
  vm_fluid_species_background_init_none(bg);
  if (!f->info.background) {
    return;
  }
  if (f->eqn_type != GKYL_EQN_EULER && f->eqn_type != GKYL_EQN_ISO_EULER) {
    gkyl_exit(
      "vm_fluid_species_background: a held background needs an Euler or isothermal Euler fluid."
    );
  }
  vm_fluid_species_background_alloc(app, f, bg);
}

void
vm_fluid_species_background_alloc(
  struct gkyl_vlasov_app *app, const struct vm_fluid_species *f, struct vm_background *bg
)
{
  bg->is_held = true;
  bg->accel = mkarr(app->use_gpu, 3 * app->basis.num_basis, app->local_ext.volume);
  gkyl_array_clear(bg->accel, 0.0);
  bg->rhs = mkarr(app->use_gpu, f->fluid->ncomp, f->fluid->size);
  gkyl_array_clear(bg->rhs, 0.0);
}

void
vm_fluid_species_background_reset(
  gkyl_vlasov_app *app, struct vm_fluid_species *f, struct vm_background *bg
)
{
  bg->rhs_func = vm_fluid_species_background_rhs_none;
  gkyl_array_clear(bg->accel, 0.0);
  vm_fluid_species_calc_app_accel(app, f, app->tcurr);
}

void
vm_fluid_species_background_project(gkyl_vlasov_app *app, struct vm_fluid_species *f, double t0)
{
  vm_fluid_species_project_state(app, f, f->info.background, f->info.background_ctx, t0);
  // The background is what the limiter leaves of it with the boundary conditions, as the initial
  // condition is.
  vm_fluid_species_apply_bc(app, f, f->fluid);
  vm_fluid_species_limiter(app, f, f->fluid);
}

void
vm_fluid_species_background_calc_force(
  gkyl_vlasov_app *app, struct vm_fluid_species *f, struct vm_background *bg,
  const struct gkyl_array *fluidin, const struct gkyl_array *rhs
)
{
  gkyl_dg_bin_op_mem *mem =
    app->use_gpu ? gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis) :
                   gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);
  for (int d = 0; d < 3; ++d) {
    gkyl_dg_div_op_range(mem, &app->basis, d, bg->accel, d + 1, rhs, 0, fluidin, &app->local);
  }
  gkyl_array_scale(bg->accel, -1.0);
  gkyl_dg_bin_op_mem_release(mem);
  vm_fluid_species_calc_app_accel(app, f, app->tcurr);
}

void
vm_fluid_species_background_write(
  gkyl_vlasov_app *app, const struct vm_fluid_species *f, const struct vm_background *bg
)
{
  struct gkyl_array *accel_ho = mkarr(false, bg->accel->ncomp, bg->accel->size);
  gkyl_array_copy(accel_ho, bg->accel);
  const char *fmt = "%s-%s_background_accel.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, f->name);
  char fname[sz + 1];
  snprintf(fname, sizeof fname, fmt, app->name, f->name);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 0, accel_ho, fname);
  gkyl_array_release(accel_ho);
}

void
vm_fluid_species_background_store(struct vm_background *bg, const struct gkyl_array *rhs)
{
  gkyl_array_copy(bg->rhs, rhs);
  bg->rhs_func = vm_fluid_species_background_rhs_remainder;
}

void
vm_fluid_species_background_rhs(const struct vm_background *bg, struct gkyl_array *rhs)
{
  bg->rhs_func(bg, rhs);
}

void
vm_fluid_species_background_release(
  const struct gkyl_vlasov_app *app, const struct vm_background *bg
)
{
  if (bg->is_held) {
    gkyl_array_release(bg->accel);
    gkyl_array_release(bg->rhs);
  }
}
