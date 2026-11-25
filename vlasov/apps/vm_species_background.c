#include <gkyl_dg_bin_ops.h>
#include <gkyl_vlasov_priv.h>

// Held background of a kinetic species (First_Fourier_Harmonic_Notes/background_density_gradient.tex).
// The background is a projection in the format of those of the initial condition. Its force, the
// first velocity moment of the right-hand side on the background divided by the density, joins the
// applied acceleration of the collisionless terms (it acts on everything that evolves about the
// background); the remainder, the right-hand side with that force, is subtracted in every stage.
// The capture is in gkyl_vlasov_app_apply_ic (vlasov.c).

// Stage methods: nothing (no background, or before the remainder is stored), or its subtraction.
static void
vm_species_background_rhs_none(const struct vm_background *bg, struct gkyl_array *rhs)
{
}

static void
vm_species_background_rhs_remainder(const struct vm_background *bg, struct gkyl_array *rhs)
{
  gkyl_array_accumulate(rhs, -1.0, bg->rhs);
}

void
vm_species_background_init_none(struct vm_background *bg)
{
  *bg = (struct vm_background){.rhs_func = vm_species_background_rhs_none};
}

void
vm_species_background_init(
  struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_background *bg
)
{
  vm_species_background_init_none(bg);
  const struct gkyl_vlasov_projection *inp = &s->info.background;
  if (!(inp->proj_id == GKYL_PROJ_FUNC ? inp->func != 0 : inp->density != 0)) {
    return;
  }
  if (s->model_id != GKYL_MODEL_DEFAULT || !s->vel_map->is_identity || !s->pos_map->is_identity ||
      s->info.is_static) {
    gkyl_exit(
      "vm_species_background: a held background needs a non-relativistic, evolving species without velocity or position maps."
    );
  }
  vm_species_projection_init(app, s, *inp, &bg->proj);
  bg->accel = mkarr(app->use_gpu, 3 * app->basis.num_basis, app->local_ext.volume);
  gkyl_array_clear(bg->accel, 0.0);
  vm_species_background_alloc_rhs(app, s, bg);
}

void
vm_species_background_alloc_rhs(
  struct gkyl_vlasov_app *app, const struct vm_species *s, struct vm_background *bg
)
{
  bg->is_held = true;
  bg->rhs = mkarr(app->use_gpu, s->f->ncomp, s->f->size);
  gkyl_array_clear(bg->rhs, 0.0);
}

void
vm_species_background_reset(gkyl_vlasov_app *app, struct vm_species *s, struct vm_background *bg)
{
  bg->rhs_func = vm_species_background_rhs_none;
  if (bg->accel) {
    gkyl_array_clear(bg->accel, 0.0);
    vm_species_collisionless_app_accel(app, &s->collisionless, app->tcurr);
  }
}

void
vm_species_background_project(
  gkyl_vlasov_app *app, struct vm_species *s, struct vm_background *bg, double t0
)
{
  vm_species_projection_calc(app, s, &bg->proj, s->f, t0);
}

void
vm_species_background_calc_force(
  gkyl_vlasov_app *app, struct vm_species *s, struct vm_background *bg,
  const struct gkyl_array *fin, const struct gkyl_array *rhs
)
{
  vm_species_moment_calc(&s->m1i, s->local, app->local, rhs);
  vm_species_moment_calc(&s->m0, s->local, app->local, fin);
  gkyl_dg_bin_op_mem *mem =
    app->use_gpu ? gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis) :
                   gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);
  gkyl_array_clear(bg->accel, 0.0);
  for (int d = 0; d < app->vdim; ++d) {
    gkyl_dg_div_op_range(mem, &app->basis, d, bg->accel, d, s->m1i.marr, 0, s->m0.marr, &app->local);
  }
  gkyl_array_scale(bg->accel, -1.0);
  gkyl_dg_bin_op_mem_release(mem);
  vm_species_collisionless_app_accel(app, &s->collisionless, app->tcurr);
}

void
vm_species_background_write(
  gkyl_vlasov_app *app, const struct vm_species *s, const struct vm_background *bg
)
{
  struct gkyl_array *accel_ho = mkarr(false, bg->accel->ncomp, bg->accel->size);
  gkyl_array_copy(accel_ho, bg->accel);
  const char *fmt = "%s-%s_background_accel.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, s->name);
  char fname[sz + 1];
  snprintf(fname, sizeof fname, fmt, app->name, s->name);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 0, accel_ho, fname);
  gkyl_array_release(accel_ho);
}

void
vm_species_background_store(struct vm_background *bg, const struct gkyl_array *rhs)
{
  gkyl_array_copy(bg->rhs, rhs);
  bg->rhs_func = vm_species_background_rhs_remainder;
}

void
vm_species_background_rhs(const struct vm_background *bg, struct gkyl_array *rhs)
{
  bg->rhs_func(bg, rhs);
}

void
vm_species_background_release(const struct gkyl_vlasov_app *app, const struct vm_background *bg)
{
  if (!bg->is_held) {
    return;
  }
  // A kinetic species has the projector and the force; the kinetic aspect of a PKPM species neither.
  if (bg->accel) {
    vm_species_projection_release(app, &bg->proj);
    gkyl_array_release(bg->accel);
  }
  gkyl_array_release(bg->rhs);
}
