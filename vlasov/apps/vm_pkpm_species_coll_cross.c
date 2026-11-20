#include <assert.h>
#include <gkyl_vlasov_priv.h>

// Collision frequencies computed from the local moments and cross-species collisions, for the
// LBO and BGK operators alike, at the configuration-space nodes. A frequency (self, or with a
// partner) is the profile of the input function, or, without one, the Spitzer frequency
//   nu_sr = nu_frac alpha_E (m_s + m_r) / (2 m_s) n_r / (vtsq_s + vtsq_r)^(3/2)
// (nu_ss = nu_frac alpha_E n_s / (2 vtsq_s)^(3/2)), with the constant part alpha_E of Morse's
// energy-exchange coefficient at the reference parameters (no magnetic field in the Coulomb
// logarithm), so m_s n_s nu_sr = m_r n_r nu_rs.
//
// Cross-species collisions have two parts that meet only through a rate (gkyl_pkpm_coll_cross).
// The kinetic operator of each species, the self operator in the normalized variables, takes a
// rate that includes its partners; it is part of the LBO or BGK update of the species. The
// momentum and energy exchange enters the fluid states: one backward-Euler step over all species
// in each SSP-RK3 stage (vlasov_update_implicit.c), which conserves the momentum and energy of the
// species at every node and does not depend on the operators or on their time integration.

// Operations of a species that needs none.
static void
vm_pkpm_species_coll_cross_prim_none(
  gkyl_vlasov_app *app, struct vlasov_species *s, const struct gkyl_array *fluid
)
{
}

static void
vm_pkpm_species_coll_cross_moms_none(gkyl_vlasov_app *app, struct vlasov_species *s)
{
}

// Node values [n, u (3), vtsq] of a fluid state of the species.
static void
vm_pkpm_species_coll_cross_prim_nodal(
  gkyl_vlasov_app *app, struct vlasov_species *s, const struct gkyl_array *fluid
)
{
  gkyl_pkpm_coll_cross_prim(
    s->kinetic->pkpm.cross.up, &app->local, s->mass, app->pos_map->jacob_pos_gauss, fluid,
    s->kinetic->pkpm.cross.prim_fluid
  );
}

// Node values of a fluid state of the species and its self-collision frequency from them.
static void
vm_pkpm_species_coll_cross_prim_self_nu(
  gkyl_vlasov_app *app, struct vlasov_species *s, const struct gkyl_array *fluid
)
{
  struct vm_pkpm_coll_cross *cross = &s->kinetic->pkpm.cross;
  vm_pkpm_species_coll_cross_prim_nodal(app, s, fluid);
  gkyl_pkpm_coll_nu_advance(
    cross->up_self_nu, &app->local, cross->self_nu_ref, cross->prim_fluid, cross->prim_fluid,
    cross->nu_kin_self
  );
}

// The two frequencies of the pair (species, partner i) at the nodes, into cross->nu_self and
// cross->nu_other. The partner forms the same two from the same node values.
static void
vm_pkpm_species_coll_cross_pair_nu(gkyl_vlasov_app *app, struct vlasov_species *s, int i)
{
  struct vm_pkpm_coll_cross *cross = &s->kinetic->pkpm.cross;
  const struct vm_pkpm_coll_cross *other = &cross->collide_with[i]->kinetic->pkpm.cross;
  gkyl_pkpm_coll_nu_advance(
    cross->up_cross_nu, &app->local, cross->nu_ref[i], cross->prim_fluid, other->prim_fluid,
    cross->nu_self
  );
  gkyl_pkpm_coll_nu_advance(
    other->up_cross_nu, &app->local, other->nu_ref[cross->other_idx[i]], other->prim_fluid,
    cross->prim_fluid, cross->nu_other
  );
}

// The rate of the kinetic operator of the stage, rebuilt from the self-collision frequency and
// the rate each partner adds.
static void
vm_pkpm_species_coll_cross_moms_nodal(gkyl_vlasov_app *app, struct vlasov_species *s)
{
  struct vm_pkpm_coll_cross *cross = &s->kinetic->pkpm.cross;
  gkyl_array_copy_range(cross->nu_kin, cross->nu_kin_self, &app->local);
  for (int i = 0; i < cross->num_cross_collisions; ++i) {
    const struct vlasov_species *other = cross->collide_with[i];
    vm_pkpm_species_coll_cross_pair_nu(app, s, i);
    gkyl_pkpm_coll_cross_kinetic_rate(
      cross->up, &app->local, s->mass, other->mass, cross->nu_self, cross->nu_other,
      cross->prim_fluid, other->kinetic->pkpm.cross.prim_fluid, cross->nu_kin
    );
  }
  gkyl_pkpm_coll_cross_floor(cross->up, &app->local, cross->nu_kin);
}

// Reference density, temperature, Planck constant and permittivity of the Coulomb logarithm of
// a computed frequency (1 where the input does not set them).
static void
vm_pkpm_species_coll_cross_ref_params(const struct gkyl_vlasov_collisions *coll, double ref[4])
{
  ref[0] = coll->den_ref ? coll->den_ref : 1.0;
  ref[1] = coll->temp_ref ? coll->temp_ref : 1.0;
  ref[2] = coll->hbar ? coll->hbar : 1.0;
  ref[3] = coll->eps0 ? coll->eps0 : 1.0;
}

// Every partner must list this species, and the cross frequencies of a pair are either both
// given (cross_nu) or both computed.
void
vm_pkpm_species_coll_cross_init(gkyl_vlasov_app *app, struct vlasov_species *s)
{
  struct vm_pkpm_coll_cross *cross = &s->kinetic->pkpm.cross;
  const struct gkyl_vlasov_collisions *coll = &s->kinetic->info.collisions;
  bool is_lbo = s->kinetic->pkpm.collision_id == GKYL_LBO_COLLISIONS;
  bool is_bgk = s->kinetic->pkpm.collision_id == GKYL_BGK_COLLISIONS;
  bool has_self_nu = (is_lbo || is_bgk) && !coll->self_nu; // computed self frequency

  cross->num_cross_collisions = (is_lbo || is_bgk) ? coll->num_cross_collisions : 0;
  cross->is_active = has_self_nu || cross->num_cross_collisions > 0;
  cross->calc_prim = vm_pkpm_species_coll_cross_prim_none;
  cross->calc_step_prim = vm_pkpm_species_coll_cross_prim_none;
  cross->calc_moms = vm_pkpm_species_coll_cross_moms_none;
  if (!cross->is_active) {
    return;
  }
  // Neither is available in the expansion-based form of the operators (GPU runs).
  assert(!app->use_gpu);

  int num_nodes = app->basis.num_basis;
  double nu_frac = coll->nu_frac > 0.0 ? coll->nu_frac : 1.0;
  bool is_computed = !coll->cross_nu[0];
  double ref_self[4];
  vm_pkpm_species_coll_cross_ref_params(coll, ref_self);
  double mass_self = s->mass;

  for (int i = 0; i < cross->num_cross_collisions; ++i) {
    struct vlasov_species *other = vlasov_find_species(app, coll->collide_with[i]);
    const struct gkyl_vlasov_collisions *other_coll = &other->kinetic->info.collisions;
    cross->collide_with[i] = other;
    cross->other_idx[i] = -1;
    for (int j = 0; j < other_coll->num_cross_collisions; ++j) {
      if (0 == strcmp(other_coll->collide_with[j], s->name)) {
        cross->other_idx[i] = j;
      }
    }
    assert(cross->other_idx[i] >= 0);
    assert(!coll->cross_nu[i] == is_computed);
    assert(!other_coll->cross_nu[cross->other_idx[i]] == is_computed);
    if (is_computed) {
      double ref_other[4];
      vm_pkpm_species_coll_cross_ref_params(other_coll, ref_other);
      double mass_other = other->mass;
      double alpha_E = gkyl_calc_Morse_alpha_E_const(
        ref_self[0], ref_other[0], mass_self, mass_other, s->charge, other->charge, ref_self[1],
        ref_other[1], 0.0, ref_self[3], ref_self[2]
      );
      cross->nu_ref[i] = mkarr(false, num_nodes, app->local_ext.volume);
      gkyl_array_clear(
        cross->nu_ref[i], nu_frac * alpha_E * (mass_self + mass_other) / (2.0 * mass_self)
      );
    } else {
      cross->nu_ref[i] = vm_pkpm_species_coll_profile_new(
        app, coll->cross_nu[i], coll->cross_nu_ctx[i], nu_frac, true
      );
    }
    cross->coef[i] = mkarr(false, num_nodes, app->local_ext.volume);
  }
  cross->up_cross_nu = gkyl_pkpm_coll_nu_inew(&(struct gkyl_pkpm_coll_nu_inp){
    .conf_basis = &app->basis,
    .spitzer = is_computed,
    .vtsq_min = 0.0,
    .use_gpu = false,
  });

  cross->prim_fluid = mkarr(false, 5 * num_nodes, app->local_ext.volume);
  cross->nu_self = mkarr(false, num_nodes, app->local_ext.volume);
  cross->nu_other = mkarr(false, num_nodes, app->local_ext.volume);
  cross->temp_ratio = mkarr(false, num_nodes, app->local_ext.volume);

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
  double alpha_E_self = gkyl_calc_Morse_alpha_E_const(
    ref_self[0], ref_self[0], mass_self, mass_self, s->charge, s->charge, ref_self[1], ref_self[1],
    0.0, ref_self[3], ref_self[2]
  );
  cross->self_nu_ref = mkarr(false, num_nodes, app->local_ext.volume);
  gkyl_array_clear(cross->self_nu_ref, nu_frac * alpha_E_self);
  cross->up_self_nu = gkyl_pkpm_coll_nu_inew(&(struct gkyl_pkpm_coll_nu_inp){
    .conf_basis = &app->basis,
    .spitzer = true,
    .vtsq_min = 0.0,
    .use_gpu = false,
  });
  vm_pkpm_coll_cross_prim_t calc_prim = has_self_nu ? vm_pkpm_species_coll_cross_prim_self_nu :
                                                      vm_pkpm_species_coll_cross_prim_nodal;

  // The exchange reads the node values of the stage state after its forward-Euler step. The
  // kinetic operator reads them (and the self frequency) before that step if it is explicit, with
  // the exchange if it is implicit.
  bool is_implicit = is_bgk && coll->is_implicit;
  bool has_partners = cross->num_cross_collisions > 0;
  if (is_implicit) {
    cross->calc_step_prim = calc_prim;
  } else {
    cross->calc_prim = calc_prim;
    if (has_partners) {
      cross->calc_step_prim = vm_pkpm_species_coll_cross_prim_nodal;
      cross->calc_moms = vm_pkpm_species_coll_cross_moms_nodal;
    }
  }
}

// Exchange coefficient of every pair, and the total frequency of each species (self plus
// partners; read by the implicit BGK step, rebuilt in the next stage otherwise), from the node
// values of the species.
static void
vm_pkpm_coll_cross_pair_coefs(gkyl_vlasov_app *app)
{
  struct vm_pkpm_coll_exchange *exch = &app->coll_exchange;
  for (int e = 0; e < exch->num_species; ++e) {
    struct vlasov_species *s = exch->species[e];
    struct vm_pkpm_coll_cross *cross = &s->kinetic->pkpm.cross;
    gkyl_array_copy_range(cross->nu_kin, cross->nu_kin_self, &app->local);
    for (int i = 0; i < cross->num_cross_collisions; ++i) {
      const struct vlasov_species *other = cross->collide_with[i];
      vm_pkpm_species_coll_cross_pair_nu(app, s, i);
      gkyl_pkpm_coll_cross_coef(
        cross->up, &app->local, s->mass, other->mass, cross->nu_self, cross->nu_other,
        cross->prim_fluid, other->kinetic->pkpm.cross.prim_fluid, cross->coef[i]
      );
      gkyl_array_accumulate_range(cross->nu_kin, 1.0, cross->nu_self, &app->local);
    }
  }
}

// No species has partners.
static void
vm_pkpm_coll_cross_exchange_none(gkyl_vlasov_app *app, double dt, struct gkyl_array *fluid[])
{
}

// Backward-Euler step of the exchange over the species with partners, on the given fluid states.
static void
vm_pkpm_coll_cross_exchange_nodal(gkyl_vlasov_app *app, double dt, struct gkyl_array *fluid[])
{
  struct vm_pkpm_coll_exchange *exch = &app->coll_exchange;
  vm_pkpm_coll_cross_pair_coefs(app);
  for (int e = 0; e < exch->num_species; ++e) {
    exch->inp[e].fluid = fluid[exch->species[e] - app->species];
  }
  gkyl_pkpm_coll_cross_advance(
    exch->up, &app->local, dt, app->pos_map->jacob_pos_gauss, exch->num_species, exch->inp
  );

  // The step keeps the temperatures positive at the nodes, not at the vertices, so the vertex
  // limiter and the fluid boundary conditions follow it.
  for (int e = 0; e < exch->num_species; ++e) {
    struct vlasov_species *s = exch->species[e];
    vm_pkpm_species_apply_positivity(app, s->fluid, exch->inp[e].fluid);
    vm_pkpm_fluid_species_apply_bc(app, s->fluid, exch->inp[e].fluid);
  }
}

// The species with partners and their data in the exchange.
void
vm_pkpm_coll_cross_init(gkyl_vlasov_app *app)
{
  struct vm_pkpm_coll_exchange *exch = &app->coll_exchange;
  exch->up = gkyl_pkpm_coll_cross_inew(&(struct gkyl_pkpm_coll_cross_inp){
    .conf_basis = &app->basis,
    .heat_fac = 1.0,
    .use_gpu = false,
  });
  exch->num_species = 0;
  int pos[GKYL_MAX_SPECIES]; // position of each species of the app in the exchange
  for (int i = 0; i < app->num_species; ++i) {
    pos[i] = -1;
    if (app->species[i].type == GKYL_SPECIES_PKPM &&
        app->species[i].kinetic->pkpm.cross.num_cross_collisions > 0) {
      pos[i] = exch->num_species;
      exch->species[exch->num_species++] = &app->species[i];
    }
  }
  for (int e = 0; e < exch->num_species; ++e) {
    struct vlasov_species *s = exch->species[e];
    struct gkyl_pkpm_coll_cross_species *inp = &exch->inp[e];
    inp->mass = s->mass;
    inp->prim = s->kinetic->pkpm.cross.prim_fluid;
    inp->fluid = s->fluid->fluid;
    inp->temp_ratio = s->kinetic->pkpm.cross.temp_ratio;
    inp->num_partners = s->kinetic->pkpm.cross.num_cross_collisions;
    for (int i = 0; i < inp->num_partners; ++i) {
      inp->partner[i] = pos[s->kinetic->pkpm.cross.collide_with[i] - app->species];
      inp->coef[i] = s->kinetic->pkpm.cross.coef[i];
    }
  }
  exch->exchange = exch->num_species > 0 ? vm_pkpm_coll_cross_exchange_nodal :
                                           vm_pkpm_coll_cross_exchange_none;
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

// Node values of the given fluid states (and the computed self frequencies of the implicit
// operators), of every species.
static void
vm_pkpm_coll_cross_step_prim(gkyl_vlasov_app *app, struct gkyl_array *fluid[])
{
  for (int i = 0; i < app->num_species; ++i) {
    struct vlasov_species *s = &app->species[i];
    s->kinetic->pkpm.cross.calc_step_prim(app, s, fluid[i]);
  }
}

void
vm_pkpm_coll_cross_advance(gkyl_vlasov_app *app, double dt, struct gkyl_array *fluid[])
{
  struct timespec wst = gkyl_wall_clock();
  // node values of every species before any exchange
  vm_pkpm_coll_cross_step_prim(app, fluid);
  app->coll_exchange.exchange(app, dt, fluid);
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
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
vm_pkpm_species_coll_cross_release(
  const gkyl_vlasov_app *app, const struct vm_pkpm_coll_cross *cross
)
{
  if (!cross->is_active) {
    return;
  }
  for (int i = 0; i < cross->num_cross_collisions; ++i) {
    gkyl_array_release(cross->nu_ref[i]);
    gkyl_array_release(cross->coef[i]);
  }
  gkyl_array_release(cross->self_nu_ref);
  gkyl_pkpm_coll_nu_release(cross->up_self_nu);
  gkyl_pkpm_coll_nu_release(cross->up_cross_nu);
  gkyl_array_release(cross->prim_fluid);
  gkyl_array_release(cross->nu_self);
  gkyl_array_release(cross->nu_other);
  gkyl_array_release(cross->temp_ratio);
  gkyl_array_release(cross->nu_kin_self);
  gkyl_array_release(cross->nu_kin);
  gkyl_pkpm_coll_cross_release(cross->up);
}

void
vm_pkpm_coll_cross_release(const gkyl_vlasov_app *app)
{
  gkyl_pkpm_coll_cross_release(app->coll_exchange.up);
}
