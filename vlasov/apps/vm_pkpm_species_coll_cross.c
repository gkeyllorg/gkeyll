#include <assert.h>
#include <gkyl_vlasov_priv.h>

// Collision frequencies computed from the local moments and the kinetic part of the cross-species
// collisions of a PKPM species, for the LBO and BGK operators alike, at the configuration-space
// nodes. The self frequency is the profile of the input function, or, without one, the Spitzer
// frequency nu_ss = nu_frac alpha_E n_s / (2 vtsq_s)^(3/2), with the constant part alpha_E of
// Morse's energy-exchange coefficient at the reference parameters (no magnetic field in the
// Coulomb logarithm).
//
// Cross-species collisions have two parts that meet only through a rate (gkyl_pkpm_coll_cross).
// The momentum and energy exchange enters the fluid states of the species with a fluid aspect
// (vm_coll_exchange.c). The kinetic operator of each PKPM species, the self operator in the
// normalized variables, takes a rate that includes its partners; it is part of the LBO or BGK
// update of the species and is formed here.

// Node values of a fluid state of the species and its self-collision frequency from them.
static void
vm_pkpm_species_coll_cross_prim_self_nu(
  gkyl_vlasov_app *app, struct vlasov_species *s, const struct gkyl_array *fluid
)
{
  struct vm_coll_cross *cross = &s->kinetic->pkpm.cross;
  vm_coll_cross_prim_nodal(app, s, fluid);
  gkyl_pkpm_coll_nu_advance(
    cross->up_self_nu, &app->local, cross->self_nu_ref, cross->prim_fluid, cross->prim_fluid,
    cross->nu_kin_self
  );
}

// The rate of the kinetic operator of the stage, rebuilt from the self-collision frequency and
// the rate each partner adds.
static void
vm_pkpm_species_coll_cross_moms_nodal(gkyl_vlasov_app *app, struct vlasov_species *s)
{
  struct vm_coll_cross *cross = &s->kinetic->pkpm.cross;
  gkyl_array_copy_range(cross->nu_kin, cross->nu_kin_self, &app->local);
  for (int i = 0; i < cross->num_cross_collisions; ++i) {
    const struct vlasov_species *other = cross->collide_with[i];
    vm_coll_cross_pair_nu(app, s, i);
    gkyl_pkpm_coll_cross_kinetic_rate(
      cross->up, &app->local, s->mass, other->mass, cross->nu_self, cross->nu_other,
      cross->prim_fluid, other->coll_cross->prim_fluid, cross->nu_kin
    );
  }
  // Kinetic partners: the same rate, from the node values of their moments.
  for (int k = 0; k < cross->num_kin_partners; ++k) {
    gkyl_pkpm_coll_cross_kinetic_rate(
      cross->up, &app->local, s->mass, cross->kin_partner[k]->mass, cross->kin_nu_self[k],
      cross->kin_nu_other[k], cross->prim_fluid, cross->kin_prim[k], cross->nu_kin
    );
  }
  gkyl_pkpm_coll_cross_floor(cross->up, &app->local, cross->nu_kin);
}

void
vm_pkpm_species_coll_cross_init(gkyl_vlasov_app *app, struct vlasov_species *s)
{
  struct vm_coll_cross *cross = &s->kinetic->pkpm.cross;
  const struct gkyl_vlasov_collisions *coll = &s->kinetic->info.collisions;
  bool is_lbo = s->kinetic->pkpm.collision_id == GKYL_LBO_COLLISIONS;
  bool is_bgk = s->kinetic->pkpm.collision_id == GKYL_BGK_COLLISIONS;
  bool has_self_nu = (is_lbo || is_bgk) && !coll->self_nu; // computed self frequency

  int num_partners = (is_lbo || is_bgk) ? coll->num_cross_collisions : 0;
  cross->num_cross_collisions = 0;
  cross->num_kin_partners = 0;
  s->fluid->pkpm.num_kin_exch = 0;
  cross->is_active = has_self_nu || num_partners > 0;
  cross->calc_prim = vm_coll_cross_prim_none;
  cross->calc_step_prim = vm_coll_cross_prim_none;
  cross->calc_moms = vm_coll_cross_moms_none;
  if (!cross->is_active) {
    return;
  }
  // Neither is available in the expansion-based form of the operators (GPU runs).
  assert(!app->use_gpu);

  // The partners, their frequencies and the node values of the fluid state.
  vm_coll_cross_partners_init(app, s, num_partners);

  // The kinetic operator in the normalized variables is the self operator at the rate
  // nu_sr vtsq_sr / vtsq_s (LBO), or the relaxation to the unit target at the rate
  // nu_sr + d ln(n T_s)/dt (BGK).
  if (is_lbo) {
    cross->nu_kin_self = gkyl_array_acquire(s->kinetic->pkpm.lbo.nu_self);
    cross->nu_kin = gkyl_array_acquire(s->kinetic->pkpm.lbo.nu);
    cross->up = gkyl_pkpm_coll_cross_inew(&(struct gkyl_pkpm_coll_cross_inp){
      .conf_basis = &app->basis,
      .heat_fac = 1.0,
      .use_gpu = false,
    });
  } else {
    cross->nu_kin_self = gkyl_array_acquire(s->kinetic->pkpm.bgk.rate_self);
    cross->nu_kin = gkyl_array_acquire(s->kinetic->pkpm.bgk.rate);
    cross->up = gkyl_pkpm_coll_cross_inew(&(struct gkyl_pkpm_coll_cross_inp){
      .conf_basis = &app->basis,
      .heat_fac = 2.0,
      .use_gpu = false,
    });
  }

  // The computed self frequency, with delta_ss (1 + beta) = 2.
  int num_nodes = app->basis.num_basis;
  double nu_frac = coll->nu_frac > 0.0 ? coll->nu_frac : 1.0;
  double ref_self[4];
  vm_coll_cross_ref_params(coll, ref_self);
  double alpha_E_self = gkyl_calc_Morse_alpha_E_const(
    ref_self[0], ref_self[0], s->mass, s->mass, s->charge, s->charge, ref_self[1], ref_self[1], 0.0,
    ref_self[3], ref_self[2]
  );
  cross->self_nu_ref = mkarr(false, num_nodes, app->local_ext.volume);
  gkyl_array_clear(cross->self_nu_ref, nu_frac * alpha_E_self);
  cross->up_self_nu = gkyl_pkpm_coll_nu_inew(&(struct gkyl_pkpm_coll_nu_inp){
    .conf_basis = &app->basis,
    .spitzer = true,
    .vtsq_min = 0.0,
    .use_gpu = false,
  });
  vm_coll_cross_prim_t calc_prim = has_self_nu ? vm_pkpm_species_coll_cross_prim_self_nu :
                                                 vm_coll_cross_prim_nodal;

  // The exchange reads the node values of the stage state after its forward-Euler step. The
  // kinetic operator reads them (and the self frequency) before that step if it is explicit, with
  // the exchange if it is implicit.
  bool is_implicit = is_bgk && coll->is_implicit;
  bool has_exchange = cross->num_cross_collisions > 0;
  if (is_implicit) {
    cross->calc_step_prim = calc_prim;
  } else {
    cross->calc_prim = calc_prim;
    if (has_exchange) {
      cross->calc_step_prim = vm_coll_cross_prim_nodal;
    }
    if (num_partners > 0) {
      cross->calc_moms = vm_pkpm_species_coll_cross_moms_nodal;
    }
  }
}

void
vm_pkpm_species_coll_cross_prim(
  gkyl_vlasov_app *app, struct vlasov_species *s, const struct gkyl_array *fluidin
)
{
  struct timespec wst = gkyl_wall_clock();
  s->kinetic->pkpm.cross.calc_prim(app, s, fluidin);
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
}

void
vm_pkpm_species_coll_cross_moms(gkyl_vlasov_app *app, struct vlasov_species *s)
{
  struct timespec wst = gkyl_wall_clock();
  s->kinetic->pkpm.cross.calc_moms(app, s);
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
}

void
vm_pkpm_species_coll_cross_shape(
  gkyl_vlasov_app *app, struct vlasov_species *s, double dt, struct gkyl_array *shape
)
{
  gkyl_pkpm_coll_cross_shape(
    s->kinetic->pkpm.cross.up, &app->local, dt, s->kinetic->pkpm.cross.nu_kin,
    s->kinetic->pkpm.cross.temp_ratio, shape
  );
}

void
vm_pkpm_species_coll_cross_release(const gkyl_vlasov_app *app, const struct vm_coll_cross *cross)
{
  if (!cross->is_active) {
    return;
  }
  vm_coll_cross_partners_release(cross);
  gkyl_array_release(cross->self_nu_ref);
  gkyl_pkpm_coll_nu_release(cross->up_self_nu);
  gkyl_array_release(cross->nu_kin_self);
  gkyl_array_release(cross->nu_kin);
  gkyl_pkpm_coll_cross_release(cross->up);
}
