#include <assert.h>
#include <gkyl_vlasov_priv.h>

// Cross-species collisions of the species with a fluid aspect, PKPM species and DG Euler fluids
// with gas_gamma = 5/3: the momentum and energy exchange between their fluid states at the
// configuration-space nodes. A frequency with a partner is the profile of the input function, or,
// without one, the Spitzer frequency
//   nu_sr = nu_frac alpha_E (m_s + m_r) / (2 m_s) n_r / (vtsq_s + vtsq_r)^(3/2),
// with the constant part alpha_E of Morse's energy-exchange coefficient at the reference
// parameters (no magnetic field in the Coulomb logarithm), so m_s n_s nu_sr = m_r n_r nu_rs.
//
// The exchange is one backward-Euler step over all species in each SSP-RK3 stage
// (vlasov_update_implicit.c), which conserves the momentum and energy of the species at every
// node and does not depend on the collision operators or on their time integration. A PKPM
// species adds the rate its partners give its kinetic operator (vm_pkpm_species_coll_cross.c); a
// fluid has no operator of its own.

// The collisions input of a species with a fluid aspect.
static const struct gkyl_vlasov_collisions *
vm_coll_cross_inp(const struct vlasov_species *s)
{
  return s->type == GKYL_SPECIES_PKPM ? &s->kinetic->info.collisions : &s->fluid->info.collisions;
}

void
vm_coll_cross_prim_none(
  gkyl_vlasov_app *app, struct vlasov_species *s, const struct gkyl_array *fluid
)
{
}

void
vm_coll_cross_moms_none(gkyl_vlasov_app *app, struct vlasov_species *s)
{
}

void
vm_coll_cross_prim_nodal(
  gkyl_vlasov_app *app, struct vlasov_species *s, const struct gkyl_array *fluid
)
{
  gkyl_pkpm_coll_cross_prim(
    s->coll_cross->up, &app->local, s->mass, app->pos_map->jacob_pos_gauss, fluid,
    s->coll_cross->prim_fluid
  );
}

void
vm_coll_cross_pair_nu(gkyl_vlasov_app *app, struct vlasov_species *s, int i)
{
  struct vm_coll_cross *cross = s->coll_cross;
  const struct vm_coll_cross *other = cross->collide_with[i]->coll_cross;
  gkyl_pkpm_coll_nu_advance(
    cross->up_cross_nu, &app->local, cross->nu_ref[i], cross->prim_fluid, other->prim_fluid,
    cross->nu_self
  );
  gkyl_pkpm_coll_nu_advance(
    other->up_cross_nu, &app->local, other->nu_ref[cross->other_idx[i]], other->prim_fluid,
    cross->prim_fluid, cross->nu_other
  );
}

void
vm_coll_cross_ref_params(const struct gkyl_vlasov_collisions *coll, double ref[4])
{
  ref[0] = coll->den_ref ? coll->den_ref : 1.0;
  ref[1] = coll->temp_ref ? coll->temp_ref : 1.0;
  ref[2] = coll->hbar ? coll->hbar : 1.0;
  ref[3] = coll->eps0 ? coll->eps0 : 1.0;
}

// Position of the species named nm among the partners of the input coll: among all of them, or
// (in_exchange) among those with a fluid aspect, which are the ones in the exchange. -1 if absent.
static int
vm_coll_cross_partner_pos(
  const gkyl_vlasov_app *app, const struct gkyl_vlasov_collisions *coll, const char *nm,
  bool in_exchange
)
{
  int pos = 0;
  for (int j = 0; j < coll->num_cross_collisions; ++j) {
    if (0 == strcmp(coll->collide_with[j], nm)) {
      return pos;
    }
    bool is_kinetic = vlasov_find_species(app, coll->collide_with[j])->type == GKYL_SPECIES_VLASOV;
    pos += (in_exchange && is_kinetic) ? 0 : 1;
  }
  return -1;
}

// Every partner must list this species, and the cross frequencies of a pair are either both
// given (cross_nu) or both computed.
void
vm_coll_cross_partners_init(gkyl_vlasov_app *app, struct vlasov_species *s, int num_partners)
{
  struct vm_coll_cross *cross = s->coll_cross;
  const struct gkyl_vlasov_collisions *coll = vm_coll_cross_inp(s);

  int num_nodes = app->basis.num_basis;
  double nu_frac = coll->nu_frac > 0.0 ? coll->nu_frac : 1.0;
  bool is_computed = !coll->cross_nu[0];
  double ref_self[4];
  vm_coll_cross_ref_params(coll, ref_self);
  double mass_self = s->mass;

  cross->num_cross_collisions = 0;
  cross->num_kin_partners = 0;
  for (int j = 0; j < num_partners; ++j) {
    struct vlasov_species *other = vlasov_find_species(app, coll->collide_with[j]);
    if (other->type == GKYL_SPECIES_VLASOV) {
      // A kinetic partner (frequencies given on both sides): the node values of its moments and
      // the two frequencies of the pair.
      const struct gkyl_vlasov_collisions *other_coll = &other->kinetic->info.collisions;
      int pos = vm_coll_cross_partner_pos(app, other_coll, s->name, false);
      assert(pos >= 0 && coll->cross_nu[j] && other_coll->cross_nu[pos]);
      double other_nu_frac = other_coll->nu_frac > 0.0 ? other_coll->nu_frac : 1.0;
      int k = cross->num_kin_partners++;
      cross->kin_partner[k] = other;
      cross->kin_nu_self[k] = vm_pkpm_species_coll_profile_new(
        app, coll->cross_nu[j], coll->cross_nu_ctx[j], nu_frac, true
      );
      cross->kin_nu_other[k] = vm_pkpm_species_coll_profile_new(
        app, other_coll->cross_nu[pos], other_coll->cross_nu_ctx[pos], other_nu_frac, true
      );
      cross->kin_prim[k] = other->kinetic->lbo.prim_nodal;
      s->fluid->pkpm.kin_exch[k] = other->kinetic->lbo.exch[pos];
      continue;
    }

    const struct gkyl_vlasov_collisions *other_coll = vm_coll_cross_inp(other);
    int i = cross->num_cross_collisions++;
    cross->collide_with[i] = other;
    cross->other_idx[i] = vm_coll_cross_partner_pos(app, other_coll, s->name, true);
    assert(cross->other_idx[i] >= 0);
    assert(!coll->cross_nu[j] == is_computed);
    assert(
      !other_coll->cross_nu[vm_coll_cross_partner_pos(app, other_coll, s->name, false)] ==
      is_computed
    );
    if (is_computed) {
      double ref_other[4];
      vm_coll_cross_ref_params(other_coll, ref_other);
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
        app, coll->cross_nu[j], coll->cross_nu_ctx[j], nu_frac, true
      );
    }
    cross->coef[i] = mkarr(false, num_nodes, app->local_ext.volume);
  }
  if (s->type == GKYL_SPECIES_PKPM) {
    s->fluid->pkpm.num_kin_exch = cross->num_kin_partners;
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
}

void
vm_coll_cross_partners_release(const struct vm_coll_cross *cross)
{
  for (int i = 0; i < cross->num_cross_collisions; ++i) {
    gkyl_array_release(cross->nu_ref[i]);
    gkyl_array_release(cross->coef[i]);
  }
  for (int k = 0; k < cross->num_kin_partners; ++k) {
    gkyl_array_release(cross->kin_nu_self[k]);
    gkyl_array_release(cross->kin_nu_other[k]);
  }
  gkyl_pkpm_coll_nu_release(cross->up_cross_nu);
  gkyl_array_release(cross->prim_fluid);
  gkyl_array_release(cross->nu_self);
  gkyl_array_release(cross->nu_other);
  gkyl_array_release(cross->temp_ratio);
}

// A fluid species in the exchange: the node values of its state and its partners.
void
vm_fluid_species_coll_cross_init(gkyl_vlasov_app *app, struct vlasov_species *s)
{
  struct vm_coll_cross *cross = s->coll_cross;
  int num_partners = s->fluid->info.collisions.num_cross_collisions;
  cross->num_cross_collisions = 0;
  cross->num_kin_partners = 0;
  cross->is_active = num_partners > 0;
  cross->calc_prim = vm_coll_cross_prim_none;
  cross->calc_step_prim = vm_coll_cross_prim_none;
  cross->calc_moms = vm_coll_cross_moms_none;
  if (!cross->is_active) {
    return;
  }
  // The exchange is at the nodes on the host.
  assert(!app->use_gpu);

  vm_coll_cross_partners_init(app, s, num_partners);
  cross->up = gkyl_pkpm_coll_cross_inew(&(struct gkyl_pkpm_coll_cross_inp){
    .conf_basis = &app->basis,
    .heat_fac = 1.0,
    .use_gpu = false,
  });
  // A fluid has no collision operator: its rate holds the sum of the frequencies with the
  // partners, which the exchange forms for every species.
  int num_nodes = app->basis.num_basis;
  cross->nu_kin_self = mkarr(false, num_nodes, app->local_ext.volume);
  gkyl_array_clear(cross->nu_kin_self, 0.0);
  cross->nu_kin = mkarr(false, num_nodes, app->local_ext.volume);

  // The node values of the stage state for the kinetic operators of PKPM partners, and those of
  // the state after its forward-Euler step for the exchange.
  cross->calc_prim = vm_coll_cross_prim_nodal;
  cross->calc_step_prim = vm_coll_cross_prim_nodal;
}

void
vm_fluid_species_coll_cross_release(const gkyl_vlasov_app *app, const struct vm_coll_cross *cross)
{
  if (!cross->is_active) {
    return;
  }
  vm_coll_cross_partners_release(cross);
  gkyl_pkpm_coll_cross_release(cross->up);
  gkyl_array_release(cross->nu_kin_self);
  gkyl_array_release(cross->nu_kin);
}

// Exchange coefficient of every pair, and the total frequency of each species (self plus
// partners; read by the implicit BGK step of a PKPM species, rebuilt in the next stage otherwise),
// from the node values of the species.
static void
vm_coll_exchange_pair_coefs(gkyl_vlasov_app *app)
{
  struct vm_coll_exchange *exch = &app->coll_exchange;
  for (int e = 0; e < exch->num_species; ++e) {
    struct vlasov_species *s = exch->species[e];
    struct vm_coll_cross *cross = s->coll_cross;
    gkyl_array_copy_range(cross->nu_kin, cross->nu_kin_self, &app->local);
    for (int i = 0; i < cross->num_cross_collisions; ++i) {
      const struct vlasov_species *other = cross->collide_with[i];
      vm_coll_cross_pair_nu(app, s, i);
      gkyl_pkpm_coll_cross_coef(
        cross->up, &app->local, s->mass, other->mass, cross->nu_self, cross->nu_other,
        cross->prim_fluid, other->coll_cross->prim_fluid, cross->coef[i]
      );
      gkyl_array_accumulate_range(cross->nu_kin, 1.0, cross->nu_self, &app->local);
    }
  }
}

// No species has partners.
static void
vm_coll_exchange_none(gkyl_vlasov_app *app, double dt, struct gkyl_array *fluid[])
{
}

// Backward-Euler step of the exchange over the species with partners, on the given fluid states.
static void
vm_coll_exchange_nodal(gkyl_vlasov_app *app, double dt, struct gkyl_array *fluid[])
{
  struct vm_coll_exchange *exch = &app->coll_exchange;
  vm_coll_exchange_pair_coefs(app);
  for (int e = 0; e < exch->num_species; ++e) {
    exch->inp[e].fluid = fluid[exch->species[e] - app->species];
  }
  gkyl_pkpm_coll_cross_advance(
    exch->up, &app->local, dt, app->pos_map->jacob_pos_gauss, exch->num_species, exch->inp
  );

  // The step keeps the temperatures positive at the nodes, not at the vertices, so the vertex
  // limiter and the fluid boundary conditions of each species follow it.
  for (int e = 0; e < exch->num_species; ++e) {
    struct vlasov_species *s = exch->species[e];
    vm_fluid_species_apply_positivity(app, s->fluid, exch->inp[e].fluid);
    vm_fluid_species_apply_bc(app, s->fluid, exch->inp[e].fluid);
  }
}

// The species with collision data, those with partners, and the data of the latter in the
// exchange.
void
vm_coll_exchange_init(gkyl_vlasov_app *app)
{
  struct vm_coll_exchange *exch = &app->coll_exchange;
  exch->num_active = 0;
  exch->num_species = 0;
  int pos[GKYL_MAX_SPECIES]; // position of each species of the app in the exchange
  for (int i = 0; i < app->num_species; ++i) {
    struct vlasov_species *s = &app->species[i];
    pos[i] = -1;
    if (s->coll_cross && s->coll_cross->is_active) {
      exch->active[exch->num_active++] = s;
      if (s->coll_cross->num_cross_collisions > 0) {
        pos[i] = exch->num_species;
        exch->species[exch->num_species++] = s;
      }
    }
  }
  for (int e = 0; e < exch->num_species; ++e) {
    struct vlasov_species *s = exch->species[e];
    struct gkyl_pkpm_coll_cross_species *inp = &exch->inp[e];
    inp->mass = s->mass;
    inp->prim = s->coll_cross->prim_fluid;
    inp->fluid = s->fluid->fluid;
    inp->temp_ratio = s->coll_cross->temp_ratio;
    inp->num_partners = s->coll_cross->num_cross_collisions;
    for (int i = 0; i < inp->num_partners; ++i) {
      inp->partner[i] = pos[s->coll_cross->collide_with[i] - app->species];
      inp->coef[i] = s->coll_cross->coef[i];
    }
  }
  exch->exchange = vm_coll_exchange_none;
  if (exch->num_species > 0) {
    exch->up = gkyl_pkpm_coll_cross_inew(&(struct gkyl_pkpm_coll_cross_inp){
      .conf_basis = &app->basis,
      .heat_fac = 1.0,
      .use_gpu = false,
    });
    exch->exchange = vm_coll_exchange_nodal;
  }
}

void
vm_coll_exchange_advance(gkyl_vlasov_app *app, double dt, struct gkyl_array *fluid[])
{
  struct timespec wst = gkyl_wall_clock();
  struct vm_coll_exchange *exch = &app->coll_exchange;
  // Node values of the given fluid states (and the computed self frequencies of the implicit
  // operators) of every species, before any exchange.
  for (int i = 0; i < exch->num_active; ++i) {
    struct vlasov_species *s = exch->active[i];
    s->coll_cross->calc_step_prim(app, s, fluid[s - app->species]);
  }
  exch->exchange(app, dt, fluid);
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void
vm_coll_exchange_release(const gkyl_vlasov_app *app)
{
  if (app->coll_exchange.num_species > 0) {
    gkyl_pkpm_coll_cross_release(app->coll_exchange.up);
  }
}
