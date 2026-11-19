#include <assert.h>
#include <gkyl_vlasov_priv.h>

// The temperature source S = nu_Q(x) (f_M(n, u, T_Q(x)) - f) of a kinetic species: a BGK
// operator whose target has the density and the drift of the species and the prescribed
// temperature T_Q. Density and momentum are untouched and every other moment relaxes to that
// of the target at the rate nu_Q, so d_t T = nu_Q (T_Q - T). The explicit form is independent
// of the collision operator of the species. The implicit form is one backward-Euler step
// together with the implicit BGK collisions of the species, if it has them:
//   f' = (f + dt nu_c f_M(n, u, T') + dt nu_Q f_M(n, u, T_Q))/(1 + dt (nu_c + nu_Q)),
// with cell-average rates in the denominator. The self-collision target is f_M(n, u, T'),
// T' = (T + dt nu_Q T_Q)/(1 + dt nu_Q); the cross-collision targets are built from T' too.

// Moments of the target: the LTE moments (n, V_drift, T/m) of the input distribution with the
// temperature replaced by T_Q/m. The rate is zeroed in cells where the average density or T_Q is
// not positive.
static void
source_bgk_target_moms(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_source_bgk *src,
  const struct gkyl_array *fin
)
{
  vm_species_lte_moms(app, vms, &vms->lte, fin);
  gkyl_array_copy_range(src->moms, vms->lte.moms.marr, &app->local);
  gkyl_array_set_offset_range(
    src->moms, 1.0, src->vtsq, (app->vdim + 1) * app->basis.num_basis, &app->local
  );
  gkyl_bgk_collisions_correct_nu(src->up_bgk, &app->local, src->moms, src->ref_nu, src->nu);
}

// The rate times the target Maxwellian, projected from the moments of the target.
static void
source_bgk_target(gkyl_vlasov_app *app, struct vm_species *vms, struct vm_source_bgk *src)
{
  vm_species_lte_from_moms(app, vms, &vms->lte, src->moms);
  gkyl_dg_mul_conf_phase_op_range(
    &app->basis, &vms->basis, src->nu_f_lte, src->nu, vms->lte.f_lte, &app->local, &vms->local
  );
}

static void
source_bgk_moms_disabled(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_source_bgk *src,
  const struct gkyl_array *fin
)
{
}

static void
source_bgk_moms_enabled(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_source_bgk *src,
  const struct gkyl_array *fin
)
{
  struct timespec wst = gkyl_wall_clock();
  source_bgk_target_moms(app, vms, src, fin);
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
}

static void
source_bgk_rhs_disabled(
  gkyl_vlasov_app *app, struct vm_species *vms, struct vm_source_bgk *src,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
}

static void
source_bgk_rhs_enabled(
  gkyl_vlasov_app *app, struct vm_species *vms, struct vm_source_bgk *src,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  struct timespec wst = gkyl_wall_clock();
  source_bgk_target(app, vms, src);
  gkyl_bgk_collisions_advance(
    src->up_bgk, &app->local, &vms->local, src->nu, src->nu_f_lte, fin, false, 0.0, rhs,
    vms->cflrate
  );
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

static void
source_bgk_moms_implicit_disabled(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_source_bgk *src,
  const struct gkyl_array *fin, double dt
)
{
}

// Moments of the implicit step: those of the target, and the temperature of the species after
// the step, T' = (T + dt nu_Q T_Q)/(1 + dt nu_Q), written into the LTE moments of the species,
// where the implicit collisions of the species and of its partners read their targets.
static void
source_bgk_moms_implicit_enabled(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_source_bgk *src,
  const struct gkyl_array *fin, double dt
)
{
  struct timespec wst = gkyl_wall_clock();
  source_bgk_target_moms(app, vms, src, fin);

  int temp_off = (app->vdim + 1) * app->basis.num_basis;
  gkyl_array_set_offset_range(src->vtsq_new, 1.0, vms->lte.moms.marr, temp_off, &app->local);
  // numerator T/m + dt nu_Q T_Q/m and denominator 1 + dt nu_Q
  gkyl_dg_mul_op_range(&app->basis, 0, src->div_num, 0, src->nu, 0, src->vtsq, &app->local);
  gkyl_array_scale_range(src->div_num, dt, &app->local);
  gkyl_array_accumulate_range(src->div_num, 1.0, src->vtsq_new, &app->local);
  gkyl_array_set_range(src->div_den, dt, src->nu, &app->local);
  gkyl_array_shiftc_range(src->div_den, pow(sqrt(2.0), app->cdim), 0, &app->local);
  gkyl_dg_div_op_range(
    src->div_mem, &app->basis, 0, src->vtsq_new, 0, src->div_num, 0, src->div_den, &app->local
  );
  gkyl_array_set_offset_range(vms->lte.moms.marr, 1.0, src->vtsq_new, temp_off, &app->local);

  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
}

static void
source_bgk_implicit_disabled(
  gkyl_vlasov_app *app, struct vm_species *vms, struct vm_source_bgk *src,
  const struct gkyl_array *fin, double dt, struct gkyl_array *out
)
{
}

// Implicit step of the source alone: (nu_Q f_M - nu_Q f)/(1 + dt nu_Q) onto out.
static void
source_bgk_implicit_alone(
  gkyl_vlasov_app *app, struct vm_species *vms, struct vm_source_bgk *src,
  const struct gkyl_array *fin, double dt, struct gkyl_array *out
)
{
  struct timespec wst = gkyl_wall_clock();
  source_bgk_target(app, vms, src);
  gkyl_bgk_collisions_advance(
    src->up_bgk, &app->local, &vms->local, src->nu, src->nu_f_lte, fin, true, dt, out, vms->cflrate
  );
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

// Implicit step of the source and the implicit BGK collisions together: the targets and the
// frequencies of the two are summed and the operator is applied once.
static void
source_bgk_implicit_with_collisions(
  gkyl_vlasov_app *app, struct vm_species *vms, struct vm_source_bgk *src,
  const struct gkyl_array *fin, double dt, struct gkyl_array *out
)
{
  struct timespec wst = gkyl_wall_clock();
  // Collision and source targets and rates, summed.
  vm_species_bgk_target_implicit(app, vms, &vms->bgk, fin);
  source_bgk_target(app, vms, src);
  gkyl_array_accumulate(src->nu_f_lte, 1.0, vms->bgk.nu_f_lte);
  gkyl_array_set(src->nu_tot, 1.0, src->nu);
  gkyl_array_accumulate(src->nu_tot, 1.0, vms->bgk.nu_sum);
  gkyl_bgk_collisions_advance(
    src->up_bgk, &app->local, &vms->local, src->nu_tot, src->nu_f_lte, fin, true, dt, out,
    vms->cflrate
  );
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void
vm_species_source_bgk_init(
  struct gkyl_vlasov_app *app, struct vm_species *vms, struct vm_source_bgk *src
)
{
  const struct gkyl_vlasov_source_bgk *info = &vms->info.source_bgk;
  src->source_bgk_id = info->source_bgk_id;
  src->is_implicit = false;
  src->moms_func = source_bgk_moms_disabled;
  src->rhs_func = source_bgk_rhs_disabled;
  src->moms_implicit_func = source_bgk_moms_implicit_disabled;
  src->implicit_func = source_bgk_implicit_disabled;
  if (src->source_bgk_id == GKYL_NO_SOURCE_BGK) {
    return;
  }
  // The prescribed-temperature source is the only model.
  assert(src->source_bgk_id == GKYL_SOURCE_BGK_HEATING);
  assert(info->rate_profile && info->temp_shape);

  // Source rate nu_Q and source temperature over the mass, T_Q/m, projected once.
  src->nu = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  src->ref_nu = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  src->vtsq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  struct gkyl_array *proj_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  vm_species_proj_conf_func(app, vms, info->rate_profile, info->rate_profile_ctx, proj_ho);
  gkyl_array_copy(src->ref_nu, proj_ho);
  gkyl_array_copy(src->nu, proj_ho);
  vm_species_proj_conf_func(app, vms, info->temp_shape, info->temp_shape_ctx, proj_ho);
  gkyl_array_copy(src->vtsq, proj_ho);
  gkyl_array_scale(src->vtsq, 1.0 / vms->mass);
  gkyl_array_release(proj_ho);

  // Moments of the target, and the rate times the target Maxwellian.
  src->moms = mkarr(app->use_gpu, vms->lte.moms.marr->ncomp, app->local_ext.volume);
  src->nu_f_lte = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume);

  src->up_bgk = gkyl_bgk_collisions_new(&app->basis, &vms->basis, app->use_gpu);

  // A static species does not take its source.
  src->is_implicit = info->is_implicit && !vms->info.is_static;
  src->has_implicit_coll = false;
  if (src->is_implicit) {
    // Scratch for the temperature after the step, and the sum of the rates when the species
    // also has implicit BGK collisions (its collision object is initialized before the source).
    src->vtsq_new = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    src->div_num = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    src->div_den = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    src->div_mem = app->use_gpu ?
                     gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis) :
                     gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);
    src->moms_implicit_func = source_bgk_moms_implicit_enabled;
    src->has_implicit_coll = vms->bgk.collision_id == GKYL_BGK_COLLISIONS && vms->bgk.implicit_step;
    if (src->has_implicit_coll) {
      src->nu_tot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      src->implicit_func = source_bgk_implicit_with_collisions;
    } else {
      src->implicit_func = source_bgk_implicit_alone;
    }
  } else if (!vms->info.is_static) {
    src->moms_func = source_bgk_moms_enabled;
    src->rhs_func = source_bgk_rhs_enabled;
  }
}

void
vm_species_source_bgk_moms(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_source_bgk *src,
  const struct gkyl_array *fin
)
{
  src->moms_func(app, vms, src, fin);
}

void
vm_species_source_bgk_rhs(
  gkyl_vlasov_app *app, struct vm_species *vms, struct vm_source_bgk *src,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  src->rhs_func(app, vms, src, fin, rhs);
}

void
vm_species_source_bgk_moms_implicit(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_source_bgk *src,
  const struct gkyl_array *fin, double dt
)
{
  src->moms_implicit_func(app, vms, src, fin, dt);
}

void
vm_species_source_bgk_implicit_advance(
  gkyl_vlasov_app *app, struct vm_species *vms, struct vm_source_bgk *src,
  const struct gkyl_array *fin, double dt, struct gkyl_array *out
)
{
  src->implicit_func(app, vms, src, fin, dt, out);
}

void
vm_species_source_bgk_release(const struct gkyl_vlasov_app *app, const struct vm_source_bgk *src)
{
  if (src->source_bgk_id == GKYL_NO_SOURCE_BGK) {
    return;
  }
  gkyl_array_release(src->nu);
  gkyl_array_release(src->ref_nu);
  gkyl_array_release(src->vtsq);
  gkyl_array_release(src->moms);
  gkyl_array_release(src->nu_f_lte);
  gkyl_bgk_collisions_release(src->up_bgk);
  if (src->is_implicit) {
    gkyl_array_release(src->vtsq_new);
    gkyl_array_release(src->div_num);
    gkyl_array_release(src->div_den);
    gkyl_dg_bin_op_mem_release(src->div_mem);
    if (src->has_implicit_coll) {
      gkyl_array_release(src->nu_tot);
    }
  }
}
