#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm_coll_cross.h>
#include <gkyl_pkpm_coll_cross_priv.h>

struct gkyl_pkpm_coll_cross *
gkyl_pkpm_coll_cross_new(const struct gkyl_basis *conf_basis, double heat_fac, bool use_gpu)
{
  assert(!use_gpu);
  assert(conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY && conf_basis->poly_order == 1);

  struct gkyl_pkpm_coll_cross *up = gkyl_malloc(sizeof(*up));
  up->num_nodes = conf_basis->num_basis;
  up->heat_fac = heat_fac;
  up->conf_to_nodal = pkpm_coll_cross_conf_to_nodal_list[conf_basis->ndim - 1];
  up->nodal_to_conf = pkpm_coll_cross_nodal_to_conf_list[conf_basis->ndim - 1];
  return up;
}

void
gkyl_pkpm_coll_cross_prim(
  const struct gkyl_pkpm_coll_cross *up, const struct gkyl_range *conf_range, double mass,
  const struct gkyl_array *jacob, const struct gkyl_array *fluid, struct gkyl_array *prim
)
{
  int num_nodes = up->num_nodes;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    const double *jacob_d = gkyl_array_cfetch(jacob, loc);
    const double *fluid_d = gkyl_array_cfetch(fluid, loc);
    double *prim_d = gkyl_array_fetch(prim, loc);

    double fluid_n[PKPM_FLUID_NUM_COMP][PKPM_MAX_CONF_BASIS];
    for (int c = 0; c < PKPM_FLUID_NUM_COMP; ++c) {
      up->conf_to_nodal(&fluid_d[c * num_nodes], fluid_n[c]);
    }
    for (int q = 0; q < num_nodes; ++q) {
      double rho = fluid_n[PKPM_FLUID_RHO][q], mom_sq = 0.0;
      prim_d[q] = rho / (jacob_d[0] * mass);
      for (int k = 0; k < 3; ++k) {
        prim_d[(PKPM_PRIM_U + k) * num_nodes + q] = fluid_n[PKPM_FLUID_MOM + k][q] / rho;
        mom_sq += fluid_n[PKPM_FLUID_MOM + k][q] * fluid_n[PKPM_FLUID_MOM + k][q];
      }
      // thermal energy (3/2) rho vtsq = E - |M|^2 / (2 rho)
      prim_d[PKPM_PRIM_VTSQ * num_nodes + q] =
        (2.0 * fluid_n[PKPM_FLUID_ENERGY][q] - mom_sq / rho) / (3.0 * rho);
    }
  }
}

void
gkyl_pkpm_coll_cross_coef(
  const struct gkyl_pkpm_coll_cross *up, const struct gkyl_range *conf_range, double mass_self,
  double mass_other, const struct gkyl_array *nu_self, const struct gkyl_array *nu_other,
  const struct gkyl_array *prim_self, const struct gkyl_array *prim_other, struct gkyl_array *coef
)
{
  int num_nodes = up->num_nodes;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    const double *nu_self_d = gkyl_array_cfetch(nu_self, loc);
    const double *nu_other_d = gkyl_array_cfetch(nu_other, loc);
    const double *self_d = gkyl_array_cfetch(prim_self, loc);
    const double *other_d = gkyl_array_cfetch(prim_other, loc);
    double *coef_d = gkyl_array_fetch(coef, loc);
    for (int q = 0; q < num_nodes; ++q) {
      // A_s = rho_s nu_sr, A_r = rho_r nu_rs
      double a_self = mass_self * self_d[q] * nu_self_d[q];
      double a_other = mass_other * other_d[q] * nu_other_d[q];
      double a_sum = a_self + a_other;
      coef_d[q] = a_sum > 0.0 ? 2.0 * a_self * a_other / a_sum : 0.0;
    }
  }
}

void
gkyl_pkpm_coll_cross_kinetic_rate(
  const struct gkyl_pkpm_coll_cross *up, const struct gkyl_range *conf_range, double mass_self,
  double mass_other, const struct gkyl_array *nu_self, const struct gkyl_array *nu_other,
  const struct gkyl_array *prim_self, const struct gkyl_array *prim_other, struct gkyl_array *nu_kin
)
{
  int num_nodes = up->num_nodes;
  double mass_sum = mass_self + mass_other, heat_fac = up->heat_fac;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    const double *nu_self_d = gkyl_array_cfetch(nu_self, loc);
    const double *nu_other_d = gkyl_array_cfetch(nu_other, loc);
    const double *self_d = gkyl_array_cfetch(prim_self, loc);
    const double *other_d = gkyl_array_cfetch(prim_other, loc);
    double *nu_kin_d = gkyl_array_fetch(nu_kin, loc);

    for (int q = 0; q < num_nodes; ++q) {
      // f_s = 2 A_r / (A_s + A_r)
      double a_self = mass_self * self_d[q] * nu_self_d[q];
      double a_other = mass_other * other_d[q] * nu_other_d[q];
      double a_sum = a_self + a_other;
      double fac_self = a_sum > 0.0 ? 2.0 * a_other / a_sum : 0.0;

      double du_sq = 0.0;
      for (int k = 0; k < 3; ++k) {
        double du =
          other_d[(PKPM_PRIM_U + k) * num_nodes + q] - self_d[(PKPM_PRIM_U + k) * num_nodes + q];
        du_sq += du * du;
      }
      double vtsq_self = self_d[PKPM_PRIM_VTSQ * num_nodes + q];
      double vtsq_other = other_d[PKPM_PRIM_VTSQ * num_nodes + q];
      // vtsq_sr - vtsq_s = f_s [(T_r - T_s) + m_r |du|^2 / 3] / (m_s + m_r)
      double dvtsq =
        (mass_other * vtsq_other - mass_self * vtsq_self + mass_other * du_sq / 3.0) / mass_sum;
      nu_kin_d[q] += nu_self_d[q] * (1.0 + heat_fac * fac_self * dvtsq / vtsq_self);
    }
  }
}

void
gkyl_pkpm_coll_cross_floor(
  const struct gkyl_pkpm_coll_cross *up, const struct gkyl_range *conf_range,
  struct gkyl_array *nu_kin
)
{
  int num_nodes = up->num_nodes;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    double *nu_kin_d = gkyl_array_fetch(nu_kin, gkyl_range_idx(conf_range, iter.idx));
    for (int q = 0; q < num_nodes; ++q) {
      nu_kin_d[q] = fmax(nu_kin_d[q], 0.0);
    }
  }
}

void
gkyl_pkpm_coll_cross_advance(
  const struct gkyl_pkpm_coll_cross *up, const struct gkyl_range *conf_range, double dt,
  const struct gkyl_array *jacob, int num_species,
  const struct gkyl_pkpm_coll_cross_species species[]
)
{
  int num_nodes = up->num_nodes;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    const double *jacob_d = gkyl_array_cfetch(jacob, loc);

    // increments of [M (3), E] of every species at the nodes
    double incr[GKYL_MAX_SPECIES][PKPM_COLL_CROSS_NUM_INCR][PKPM_MAX_CONF_BASIS];
    for (int q = 0; q < num_nodes; ++q) {
      double mom_mat[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES];
      double temp_mat[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES];
      double dens[GKYL_MAX_SPECIES], rho[GKYL_MAX_SPECIES], temp[GKYL_MAX_SPECIES];
      double flow[3][GKYL_MAX_SPECIES], flow_new[3][GKYL_MAX_SPECIES];

      // (rho_s + dt sum_r K_sr) u_s' - dt sum_r K_sr u_r' = rho_s u_s, and the same structure
      // for the temperatures with (3/2) n_s and 3 K_sr / (m_s + m_r)
      for (int s = 0; s < num_species; ++s) {
        const double *prim_d = gkyl_array_cfetch(species[s].prim, loc);
        dens[s] = prim_d[q];
        rho[s] = species[s].mass * dens[s];
        temp[s] = species[s].mass * prim_d[PKPM_PRIM_VTSQ * num_nodes + q];
        for (int k = 0; k < 3; ++k) {
          flow[k][s] = prim_d[(PKPM_PRIM_U + k) * num_nodes + q];
          flow_new[k][s] = rho[s] * flow[k][s];
        }
        for (int r = 0; r < num_species; ++r) {
          mom_mat[s][r] = 0.0;
          temp_mat[s][r] = 0.0;
        }
        mom_mat[s][s] = rho[s];
        temp_mat[s][s] = 1.5 * dens[s];
        for (int i = 0; i < species[s].num_partners; ++i) {
          int r = species[s].partner[i];
          const double *coef_d = gkyl_array_cfetch(species[s].coef[i], loc);
          double mom_rate = dt * coef_d[q];
          double temp_rate = 3.0 * mom_rate / (species[s].mass + species[r].mass);
          mom_mat[s][s] += mom_rate;
          mom_mat[s][r] -= mom_rate;
          temp_mat[s][s] += temp_rate;
          temp_mat[s][r] -= temp_rate;
        }
      }
      pkpm_coll_cross_factor(num_species, mom_mat);
      for (int k = 0; k < 3; ++k) {
        pkpm_coll_cross_solve(num_species, mom_mat, flow_new[k]);
      }

      // thermal energy: the kinetic energy the step of each species gives up, rho |du|^2 / 2,
      // and the frictional heating of each pair at the flows after the step
      double temp_new[GKYL_MAX_SPECIES], kin_incr[GKYL_MAX_SPECIES];
      for (int s = 0; s < num_species; ++s) {
        double step_sq = 0.0;
        kin_incr[s] = 0.0;
        for (int k = 0; k < 3; ++k) {
          double step = flow_new[k][s] - flow[k][s];
          step_sq += step * step;
          kin_incr[s] += rho[s] * step * (flow[k][s] + 0.5 * step);
          incr[s][k][q] = jacob_d[0] * rho[s] * step;
        }
        double heat = 0.5 * rho[s] * step_sq;
        for (int i = 0; i < species[s].num_partners; ++i) {
          int r = species[s].partner[i];
          const double *coef_d = gkyl_array_cfetch(species[s].coef[i], loc);
          double du_sq = 0.0;
          for (int k = 0; k < 3; ++k) {
            double du = flow_new[k][r] - flow_new[k][s];
            du_sq += du * du;
          }
          heat += dt * coef_d[q] * species[r].mass * du_sq / (species[s].mass + species[r].mass);
        }
        temp_new[s] = 1.5 * dens[s] * temp[s] + heat;
      }
      pkpm_coll_cross_factor(num_species, temp_mat);
      pkpm_coll_cross_solve(num_species, temp_mat, temp_new);

      for (int s = 0; s < num_species; ++s) {
        double *temp_ratio_d = gkyl_array_fetch(species[s].temp_ratio, loc);
        incr[s][PKPM_FLUID_ENERGY - PKPM_FLUID_MOM][q] =
          jacob_d[0] * (kin_incr[s] + 1.5 * dens[s] * (temp_new[s] - temp[s]));
        temp_ratio_d[q] = temp[s] / temp_new[s];
      }
    }

    for (int s = 0; s < num_species; ++s) {
      double *fluid_d = gkyl_array_fetch(species[s].fluid, loc);
      for (int c = 0; c < PKPM_COLL_CROSS_NUM_INCR; ++c) {
        double incr_modal[PKPM_MAX_CONF_BASIS];
        up->nodal_to_conf(incr[s][c], incr_modal);
        for (int i = 0; i < num_nodes; ++i) {
          fluid_d[(PKPM_FLUID_MOM + c) * num_nodes + i] += incr_modal[i];
        }
      }
    }
  }
}

void
gkyl_pkpm_coll_cross_shape(
  const struct gkyl_pkpm_coll_cross *up, const struct gkyl_range *conf_range, double dt,
  const struct gkyl_array *nu_tot, const struct gkyl_array *temp_ratio, struct gkyl_array *shape
)
{
  int num_nodes = up->num_nodes;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    const double *nu_tot_d = gkyl_array_cfetch(nu_tot, loc);
    const double *temp_ratio_d = gkyl_array_cfetch(temp_ratio, loc);
    double *shape_d = gkyl_array_fetch(shape, loc);
    for (int q = 0; q < num_nodes; ++q) {
      shape_d[q] = fmin(1.0, temp_ratio_d[q] / (1.0 + nu_tot_d[q] * dt));
    }
  }
}

void
gkyl_pkpm_coll_cross_release(struct gkyl_pkpm_coll_cross *up)
{
  gkyl_free(up);
}
