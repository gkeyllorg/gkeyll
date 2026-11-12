#include <assert.h>
#include <gkyl_pkpm_priv.h>

// Collision frequencies computed from the local moments and cross-species collisions, for the
// LBO and BGK operators alike, at the configuration-space nodes. The input follows the Vlasov
// app: a frequency (self, or with a partner) is the profile of the input function, or, without
// one, the Spitzer frequency
//   nu_sr = nuFrac alpha_E (m_s + m_r) / (2 m_s) n_r / (vtsq_s + vtsq_r)^(3/2)
// (nu_ss = nuFrac alpha_E n_s / (2 vtsq_s)^(3/2)), with the constant part alpha_E of Morse's
// energy-exchange coefficient at the reference parameters (no magnetic field in the Coulomb
// logarithm), so m_s n_s nu_sr = m_r n_r nu_rs.
//
// Cross-species collisions have two parts that meet only through a rate (gkyl_pkpm_coll_cross).
// The kinetic operator of each species, the self operator in the normalized variables, takes a
// rate that includes its partners; it is part of the LBO or BGK update of the species. The
// momentum and energy exchange enters the fluid states: one backward-Euler step over all species
// after each SSP-RK3 step, which conserves the momentum and energy of the species at every node
// and does not depend on the operators or on their time integration.

// Operations of a species that needs none.
static void
pkpm_species_coll_cross_prim_none(
  gkyl_pkpm_app *app, struct pkpm_species *s, const struct gkyl_array *fluid
)
{
}

static void
pkpm_species_coll_cross_moms_none(gkyl_pkpm_app *app, struct pkpm_species *s)
{
}

// Node values [n, u (3), vtsq] of a fluid state of the species.
static void
pkpm_species_coll_cross_prim_nodal(
  gkyl_pkpm_app *app, struct pkpm_species *s, const struct gkyl_array *fluid
)
{
  gkyl_pkpm_coll_cross_prim(
    s->cross.up, &app->local, s->info.mass, app->pos_map->jacob_pos_gauss, fluid,
    s->cross.prim_fluid
  );
}

// Node values of a fluid state of the species and its self-collision frequency from them.
static void
pkpm_species_coll_cross_prim_self_nu(
  gkyl_pkpm_app *app, struct pkpm_species *s, const struct gkyl_array *fluid
)
{
  struct pkpm_coll_cross *cross = &s->cross;
  pkpm_species_coll_cross_prim_nodal(app, s, fluid);
  gkyl_pkpm_coll_nu_advance(
    cross->up_self_nu, &app->local, cross->self_nu_ref, cross->prim_fluid, cross->prim_fluid,
    cross->nu_kin_self
  );
}

// The two frequencies of the pair (species, partner i) at the nodes, into cross->nu_self and
// cross->nu_other. The partner forms the same two from the same node values.
static void
pkpm_species_coll_cross_pair_nu(gkyl_pkpm_app *app, struct pkpm_species *s, int i)
{
  struct pkpm_coll_cross *cross = &s->cross;
  const struct pkpm_coll_cross *other = &cross->collide_with[i]->cross;
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
pkpm_species_coll_cross_moms_nodal(gkyl_pkpm_app *app, struct pkpm_species *s)
{
  struct pkpm_coll_cross *cross = &s->cross;
  gkyl_array_copy_range(cross->nu_kin, cross->nu_kin_self, &app->local);
  for (int i = 0; i < cross->num_cross_collisions; ++i) {
    const struct pkpm_species *other = cross->collide_with[i];
    pkpm_species_coll_cross_pair_nu(app, s, i);
    gkyl_pkpm_coll_cross_kinetic_rate(
      cross->up, &app->local, s->info.mass, other->info.mass, cross->nu_self, cross->nu_other,
      cross->prim_fluid, other->cross.prim_fluid, cross->nu_kin
    );
  }
  gkyl_pkpm_coll_cross_floor(cross->up, &app->local, cross->nu_kin);
}

// Reference density, temperature, Planck constant and permittivity of the Coulomb logarithm of
// a computed frequency (1 where the input does not set them).
static void
pkpm_species_coll_cross_ref_params(const struct gkyl_pkpm_collisions *coll, double ref[4])
{
  ref[0] = coll->den_ref ? coll->den_ref : 1.0;
  ref[1] = coll->temp_ref ? coll->temp_ref : 1.0;
  ref[2] = coll->hbar ? coll->hbar : 1.0;
  ref[3] = coll->eps0 ? coll->eps0 : 1.0;
}

// As in the Vlasov app, every partner must list this species, and the cross frequencies of a
// pair are either both given (cross_nu) or both computed.
void
pkpm_species_coll_cross_init(gkyl_pkpm_app *app, struct pkpm_species *s)
{
  struct pkpm_coll_cross *cross = &s->cross;
  const struct gkyl_pkpm_collisions *coll = &s->info.collisions;
  bool is_lbo = s->collision_id == GKYL_LBO_COLLISIONS;
  bool is_bgk = s->collision_id == GKYL_BGK_COLLISIONS;
  bool has_self_nu = (is_lbo || is_bgk) && !coll->self_nu; // computed self frequency

  cross->num_cross_collisions = (is_lbo || is_bgk) ? coll->num_cross_collisions : 0;
  cross->is_active = has_self_nu || cross->num_cross_collisions > 0;
  cross->calc_prim = pkpm_species_coll_cross_prim_none;
  cross->calc_step_prim = pkpm_species_coll_cross_prim_none;
  cross->calc_moms = pkpm_species_coll_cross_moms_none;
  if (!cross->is_active) {
    return;
  }
  // Neither is available in the expansion-based form of the operators (GPU runs).
  assert(!app->use_gpu);

  int num_nodes = app->confBasis.num_basis;
  double nu_frac = coll->nuFrac > 0.0 ? coll->nuFrac : 1.0;
  bool is_computed = !coll->cross_nu[0];
  double ref_self[4];
  pkpm_species_coll_cross_ref_params(coll, ref_self);
  double mass_self = s->info.mass;

  for (int i = 0; i < cross->num_cross_collisions; ++i) {
    struct pkpm_species *other = pkpm_find_species(app, coll->collide_with[i]);
    const struct gkyl_pkpm_collisions *other_coll = &other->info.collisions;
    cross->collide_with[i] = other;
    cross->other_idx[i] = -1;
    for (int j = 0; j < other_coll->num_cross_collisions; ++j) {
      if (0 == strcmp(other_coll->collide_with[j], s->info.name)) {
        cross->other_idx[i] = j;
      }
    }
    assert(cross->other_idx[i] >= 0);
    assert(!coll->cross_nu[i] == is_computed);
    assert(!other_coll->cross_nu[cross->other_idx[i]] == is_computed);
    if (is_computed) {
      double ref_other[4];
      pkpm_species_coll_cross_ref_params(other_coll, ref_other);
      double mass_other = other->info.mass;
      double alpha_E = gkyl_calc_Morse_alpha_E_const(
        ref_self[0], ref_other[0], mass_self, mass_other, s->info.charge, other->info.charge,
        ref_self[1], ref_other[1], 0.0, ref_self[3], ref_self[2]
      );
      cross->nu_ref[i] = mkarr(false, num_nodes, app->local_ext.volume);
      gkyl_array_clear(
        cross->nu_ref[i], nu_frac * alpha_E * (mass_self + mass_other) / (2.0 * mass_self)
      );
    } else {
      cross->nu_ref[i] =
        pkpm_species_coll_profile_new(app, coll->cross_nu[i], coll->cross_nu_ctx[i], nu_frac, true);
    }
    cross->coef[i] = mkarr(false, num_nodes, app->local_ext.volume);
  }
  cross->up_cross_nu = gkyl_pkpm_coll_nu_new(&app->confBasis, is_computed, 0.0, false);

  cross->prim_fluid = mkarr(false, 5 * num_nodes, app->local_ext.volume);
  cross->nu_self = mkarr(false, num_nodes, app->local_ext.volume);
  cross->nu_other = mkarr(false, num_nodes, app->local_ext.volume);
  cross->temp_ratio = mkarr(false, num_nodes, app->local_ext.volume);

  // The kinetic operator in the normalized variables is the self operator at the rate
  // nu_sr vtsq_sr / vtsq_s (LBO), or the relaxation to the unit target at the rate
  // nu_sr + d ln(n T_s)/dt (BGK).
  if (is_lbo) {
    cross->nu_kin_self = gkyl_array_acquire(s->lbo.nu_self);
    cross->nu_kin = gkyl_array_acquire(s->lbo.nu);
    cross->up = gkyl_pkpm_coll_cross_new(&app->confBasis, 1.0, false);
  } else {
    cross->nu_kin_self = gkyl_array_acquire(s->bgk.rate_self);
    cross->nu_kin = gkyl_array_acquire(s->bgk.rate);
    cross->up = gkyl_pkpm_coll_cross_new(&app->confBasis, 2.0, false);
  }

  // The computed self frequency, with delta_ss (1 + beta) = 2 as in the Vlasov app.
  double alpha_E_self = gkyl_calc_Morse_alpha_E_const(
    ref_self[0], ref_self[0], mass_self, mass_self, s->info.charge, s->info.charge, ref_self[1],
    ref_self[1], 0.0, ref_self[3], ref_self[2]
  );
  cross->self_nu_ref = mkarr(false, num_nodes, app->local_ext.volume);
  gkyl_array_clear(cross->self_nu_ref, nu_frac * alpha_E_self);
  cross->up_self_nu = gkyl_pkpm_coll_nu_new(&app->confBasis, true, 0.0, false);
  pkpm_coll_cross_prim_t calc_prim = has_self_nu ? pkpm_species_coll_cross_prim_self_nu :
                                                   pkpm_species_coll_cross_prim_nodal;

  // The exchange reads the node values after each step. The kinetic operator reads them (and
  // the self frequency) in the stages, or after the step if it is implicit.
  bool is_implicit = is_bgk && coll->is_implicit;
  bool has_partners = cross->num_cross_collisions > 0;
  if (is_implicit) {
    cross->calc_step_prim = calc_prim;
  } else {
    cross->calc_prim = calc_prim;
    if (has_partners) {
      cross->calc_step_prim = pkpm_species_coll_cross_prim_nodal;
      cross->calc_moms = pkpm_species_coll_cross_moms_nodal;
    }
  }
}

// No species has partners.
static void
pkpm_coll_cross_exchange_none(gkyl_pkpm_app *app, double dt)
{
}

// Backward-Euler step of the exchange over the species with partners.
static void
pkpm_coll_cross_exchange_nodal(gkyl_pkpm_app *app, double dt)
{
  struct pkpm_coll_exchange *exch = &app->coll_exchange;

  // Exchange coefficient of every pair, and the total frequency of each species (self plus
  // partners; read by the implicit BGK step, rebuilt in the next stage otherwise).
  for (int e = 0; e < exch->num_species; ++e) {
    struct pkpm_species *s = exch->species[e];
    struct pkpm_coll_cross *cross = &s->cross;
    gkyl_array_copy_range(cross->nu_kin, cross->nu_kin_self, &app->local);
    for (int i = 0; i < cross->num_cross_collisions; ++i) {
      const struct pkpm_species *other = cross->collide_with[i];
      pkpm_species_coll_cross_pair_nu(app, s, i);
      gkyl_pkpm_coll_cross_coef(
        cross->up, &app->local, s->info.mass, other->info.mass, cross->nu_self, cross->nu_other,
        cross->prim_fluid, other->cross.prim_fluid, cross->coef[i]
      );
      gkyl_array_accumulate_range(cross->nu_kin, 1.0, cross->nu_self, &app->local);
    }
  }

  gkyl_pkpm_coll_cross_advance(
    exch->up, &app->local, dt, app->pos_map->jacob_pos_gauss, exch->num_species, exch->inp
  );

  // The step keeps the temperatures positive at the nodes, not at the vertices, so the vertex
  // limiter and the fluid boundary conditions follow it.
  for (int e = 0; e < exch->num_species; ++e) {
    struct pkpm_species *s = exch->species[e];
    pkpm_species_apply_positivity(app, s, s->fluid);
    pkpm_fluid_species_apply_bc(app, s, s->fluid);
  }
}

// The species with partners and their data in the exchange.
void
pkpm_coll_cross_init(gkyl_pkpm_app *app)
{
  struct pkpm_coll_exchange *exch = &app->coll_exchange;
  exch->up = gkyl_pkpm_coll_cross_new(&app->confBasis, 1.0, false);
  exch->num_species = 0;
  int pos[GKYL_MAX_SPECIES]; // position of each species of the app in the exchange
  for (int i = 0; i < app->num_species; ++i) {
    pos[i] = -1;
    if (app->species[i].cross.num_cross_collisions > 0) {
      pos[i] = exch->num_species;
      exch->species[exch->num_species++] = &app->species[i];
    }
  }
  for (int e = 0; e < exch->num_species; ++e) {
    struct pkpm_species *s = exch->species[e];
    struct gkyl_pkpm_coll_cross_species *inp = &exch->inp[e];
    inp->mass = s->info.mass;
    inp->prim = s->cross.prim_fluid;
    inp->fluid = s->fluid;
    inp->temp_ratio = s->cross.temp_ratio;
    inp->num_partners = s->cross.num_cross_collisions;
    for (int i = 0; i < inp->num_partners; ++i) {
      inp->partner[i] = pos[s->cross.collide_with[i] - app->species];
      inp->coef[i] = s->cross.coef[i];
    }
  }
  exch->exchange = exch->num_species > 0 ? pkpm_coll_cross_exchange_nodal :
                                           pkpm_coll_cross_exchange_none;
}

void
pkpm_species_coll_cross_prim(
  gkyl_pkpm_app *app, struct pkpm_species *s, const struct gkyl_array *fluidin
)
{
  struct timespec wst = gkyl_wall_clock();
  s->cross.calc_prim(app, s, fluidin);
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
}

void
pkpm_species_coll_cross_moms(gkyl_pkpm_app *app, struct pkpm_species *s)
{
  struct timespec wst = gkyl_wall_clock();
  s->cross.calc_moms(app, s);
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
}

void
pkpm_coll_cross_advance(gkyl_pkpm_app *app, double dt)
{
  struct timespec wst = gkyl_wall_clock();
  // Node values of the state after the SSP-RK3 step (and the computed self frequencies of the
  // implicit operators), of every species before any exchange.
  for (int i = 0; i < app->num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    s->cross.calc_step_prim(app, s, s->fluid);
  }
  app->coll_exchange.exchange(app, dt);
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void
pkpm_species_coll_cross_shape(
  gkyl_pkpm_app *app, struct pkpm_species *s, double dt, struct gkyl_array *shape
)
{
  gkyl_pkpm_coll_cross_shape(
    s->cross.up, &app->local, dt, s->cross.nu_kin, s->cross.temp_ratio, shape
  );
}

void
pkpm_species_coll_cross_release(const gkyl_pkpm_app *app, const struct pkpm_coll_cross *cross)
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
pkpm_coll_cross_release(const gkyl_pkpm_app *app)
{
  gkyl_pkpm_coll_cross_release(app->coll_exchange.up);
}
