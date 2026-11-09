// Semi-discrete stability check of the coupled factored, thermal-velocity normalized
// PKPM update (kinetic two-pass update of the normalized distributions on the w grid
// plus the fluid system [rho, M, E]) on a 1x1v periodic entropy wave: a density pulse
// advected by a uniform flow with uniform pressure. The exact solution keeps u and p
// uniform; the test integrates the coupled zero-layer operators with SSP-RK3 and checks
// that the flow stays bounded and that the integrated fluid energy is conserved (the
// energy equation is in conservative form; the kinetic closures only redistribute it).
// (The legacy split-form momentum flux with a Lax mass flux was unstable here at a rate
// ~ u/dx.) The kinetic and fluid halves can be frozen separately to localize an
// instability (diagnostics, printed with -v).
#include <math.h>
#include <stdio.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_pkpm_conf_flux_surf.h>
#include <gkyl_dg_pkpm_fluid.h>
#include <gkyl_dg_pkpm_fluid_flux_surf.h>
#include <gkyl_dg_pkpm_harm_vol.h>
#include <gkyl_dg_pkpm_kinetic_conf.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_dg_pkpm_kinetic_vel.h>
#include <gkyl_dg_pkpm_vel_flux_surf.h>
#include <gkyl_dg_updater_moment_pkpm.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// The weak gradient: the in-cell gradient with the whole jump of every face lifted into the cell.
static void
weak_grad(
  struct gkyl_dg_calc_pkpm_vars *vars, const struct gkyl_range *range, const struct gkyl_array *q,
  struct gkyl_array *grad
)
{
  struct gkyl_array *blend = gkyl_array_new(GKYL_DOUBLE, 3, grad->size);
  gkyl_array_clear(blend, 1.0);
  gkyl_dg_calc_pkpm_vars_in_cell_grad(vars, range, blend, q, grad);
  gkyl_array_release(blend);
}

// The flux both sides of every face receive, from the single-valued label-paired flux:
// [dir][mode] -> [dir][side][mode], in place.
static void
expand_two_sided(struct gkyl_array *flux, int cdim, int nflux)
{
  for (long i = 0; i < flux->size; ++i) {
    double *f = gkyl_array_fetch(flux, i);
    for (int d = cdim - 1; d >= 0; --d) {
      for (int m = nflux - 1; m >= 0; --m) {
        double val = f[d * nflux + m];
        f[2 * d * nflux + m] = val;
        f[(2 * d + 1) * nflux + m] = val;
      }
    }
  }
}

struct wave_ctx {
  double u0, alpha, vt;
  double b_theta, b_perp; // static field: B = (sin theta, -cos theta, -b_perp cos(pi x))
  double uy_amp; // sheared perpendicular drift u_y = uy_amp cos(pi x)
  double ux_amp; // compressive flow u_x = u0 + ux_amp cos(pi x)
};

static struct gkyl_array *
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// Identity configuration-space position map (the updaters always take one).
static struct gkyl_vlasov_position_map *
identity_pos_map(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *clocal, const struct gkyl_range *clocal_ext,
  const struct gkyl_basis *cbasis
)
{
  struct gkyl_vlasov_position_map_inp inp[GKYL_MAX_CDIM] = {0};
  return gkyl_vlasov_position_map_new(cgrid, clocal, clocal_ext, cbasis, inp, false);
}

static void
eval_em(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct wave_ctx *c = ctx;
  for (int i = 0; i < 8; ++i) fout[i] = 0.0;
  if (c->b_perp == 0.0 && c->b_theta == 0.0) {
    fout[3] = 1.0; // B = x-hat
    return;
  }
  fout[3] = sin(c->b_theta);
  fout[4] = -cos(c->b_theta);
  fout[5] = -c->b_perp * cos(M_PI * xn[0]);
}

// Normalized distributions on the w grid of a Maxwellian with isotropic temperature:
// F = exp(-w^2/2)/sqrt(2 pi), G = (T_perp/T) F = F, independent of x.
static void
eval_f(double t, const double *xn, double *restrict fout, void *ctx)
{
  double w = xn[1];
  double F = exp(-0.5 * w * w) / sqrt(2.0 * M_PI);
  fout[0] = F;
  fout[1] = F;
}

// Fluid state [rho, M, E] (mass 1): n = 1 + alpha sin(pi x), uniform pressure p = vt^2,
// E = rho |u|^2/2 + 3 p/2.
static void
eval_fluid(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct wave_ctx *c = ctx;
  double n = 1.0 + c->alpha * sin(M_PI * xn[0]);
  double ux = c->u0 + c->ux_amp * cos(M_PI * xn[0]), uy = c->uy_amp * cos(M_PI * xn[0]);
  fout[0] = n;
  fout[1] = n * ux;
  fout[2] = n * uy;
  fout[3] = 0.0;
  fout[4] = 0.5 * n * (ux * ux + uy * uy) + 1.5 * c->vt * c->vt;
}

// Isothermal fluid state: n = 1 + alpha sin(pi x), uniform temperature T = vt^2 (p = n vt^2),
// so that xi is uniform and the uniform Maxwellian shape is a steady state of the kinetic
// equation (the fluid is not: it accelerates down the pressure gradient).
static void
eval_fluid_isothermal(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct wave_ctx *c = ctx;
  double x = xn[0], n = 1.0 + c->alpha * sin(M_PI * x);
  fout[0] = n;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = 1.5 * n * c->vt * c->vt;
}

struct run_opts {
  double u0, alpha, t_end;
  int Nx, Nv;
  bool freeze_fluid, freeze_kinetic;
  double b_theta, b_perp; // oblique, x-dependent static field (0, 0: b = x-hat)
  double uy_amp; // sheared perpendicular drift amplitude
  double ux_amp; // compressive flow amplitude
};

struct sim {
  struct gkyl_vlasov_position_map *pos_map; // identity configuration-space map
  struct gkyl_rect_grid cgrid, pgrid;
  struct gkyl_basis cbasis, pbasis;
  struct gkyl_range clocal, clocal_ext, plocal, plocal_ext, plocal_ext_conf, vlocal, vlocal_ext;
  struct gkyl_range cskin[2], cghost[2], pskin[2], pghost[2];
  int nc, np, nkc, nkf;
  struct gkyl_array *bvar, *grad_b, *u, *grad_u, *theta, *nu_sum, *nu_vthsq_sum, *alpha_0, *c_h, *c_0;
  struct gkyl_array *a_perp, *em, *dtb; // perpendicular force closure, static field, d_t B/|B| = 0
  struct gkyl_array *xi, *grad_xi, *scratch;
  struct gkyl_array *nmoms, *moms, *rhs_moms, *rho, *rho3, *grad_rho, *grad_ln_n;
  struct gkyl_array *cflux, *vflux, *fflux, *cflrate_f, *cflrate_m, *face_coef;
  struct gkyl_dg_calc_pkpm_bvar *calc_b;
  struct gkyl_dg_calc_pkpm_vars *vars;
  struct gkyl_dg_updater_moment *mcalc;
  struct gkyl_dg_pkpm_conf_flux_surf *calc_cflux;
  struct gkyl_dg_pkpm_vel_flux_surf *calc_vflux;
  struct gkyl_dg_pkpm_fluid_flux_surf *calc_fflux;
  struct gkyl_dg_pkpm_harm_vol *harm;
  struct gkyl_array *hcoef; // stored first-harmonic operator of the stage
  struct gkyl_array *hwn; // raw w-flux values of the stage at the w-face nodes
  struct gkyl_dg_eqn *eqn_a, *eqn_b, *eqn_m;
  gkyl_hyper_dg *slvr_a, *slvr_b, *slvr_m;
};

// Phase-space ranges with configuration ghosts on both sides and ONE velocity ghost cell
// above w_max (the empty state beyond the characteristic velocity edge, whose flux slot
// holds the upper-edge flux); plocal_ext_conf is extended in configuration space only.
static void
phase_ranges(
  const struct gkyl_rect_grid *pgrid, int cdim, struct gkyl_range *plocal_ext, struct gkyl_range *plocal,
  struct gkyl_range *plocal_ext_conf
)
{
  int pdim = pgrid->ndim, lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int elo[GKYL_MAX_DIM] = {0}, eup[GKYL_MAX_DIM] = {0};
  for (int d = 0; d < pdim; ++d) {
    lower[d] = 1;
    upper[d] = pgrid->cells[d];
  }
  for (int d = 0; d < cdim; ++d) {
    elo[d] = eup[d] = 1;
  }
  eup[cdim] = 1;
  struct gkyl_range local;
  gkyl_range_init(&local, pdim, lower, upper);
  gkyl_range_extend(plocal_ext, &local, elo, eup);
  gkyl_sub_range_init(plocal, plocal_ext, lower, upper);
  int lo_c[GKYL_MAX_DIM], up_c[GKYL_MAX_DIM];
  for (int d = 0; d < pdim; ++d) {
    lo_c[d] = d < cdim ? plocal_ext->lower[d] : lower[d];
    up_c[d] = d < cdim ? plocal_ext->upper[d] : upper[d];
  }
  gkyl_sub_range_init(plocal_ext_conf, plocal_ext, lo_c, up_c);
}

static void
periodic(struct gkyl_array *arr, const struct gkyl_range *skin, const struct gkyl_range *ghost)
{
  // lower ghost <- upper skin, upper ghost <- lower skin
  gkyl_array_copy_range_to_range(arr, arr, &ghost[0], &skin[1]);
  gkyl_array_copy_range_to_range(arr, arr, &ghost[1], &skin[0]);
}

static void
sim_init(struct sim *s, const struct run_opts *o)
{
  struct wave_ctx bctx = {.b_theta = o->b_theta, .b_perp = o->b_perp};
  double lower[2] = {-1.0, -6.0}, upper[2] = {1.0, 6.0};
  int cells[2] = {o->Nx, o->Nv};
  int cghost[2] = {1, 0}, vghost[1] = {0};
  gkyl_rect_grid_init(&s->cgrid, 1, lower, upper, cells);
  gkyl_rect_grid_init(&s->pgrid, 2, lower, upper, cells);
  gkyl_cart_modal_serendip(&s->cbasis, 1, 1);
  gkyl_cart_modal_hybrid(&s->pbasis, 1, 1);
  gkyl_create_grid_ranges(&s->cgrid, cghost, &s->clocal_ext, &s->clocal);
  s->pos_map = identity_pos_map(&s->cgrid, &s->clocal, &s->clocal_ext, &s->cbasis);
  phase_ranges(&s->pgrid, 1, &s->plocal_ext, &s->plocal, &s->plocal_ext_conf);
  struct gkyl_rect_grid vgrid;
  gkyl_rect_grid_init(&vgrid, 1, &lower[1], &upper[1], &cells[1]);
  gkyl_create_grid_ranges(&vgrid, vghost, &s->vlocal_ext, &s->vlocal);
  gkyl_skin_ghost_ranges(&s->cskin[0], &s->cghost[0], 0, GKYL_LOWER_EDGE, &s->clocal_ext, cghost);
  gkyl_skin_ghost_ranges(&s->cskin[1], &s->cghost[1], 0, GKYL_UPPER_EDGE, &s->clocal_ext, cghost);
  gkyl_skin_ghost_ranges(&s->pskin[0], &s->pghost[0], 0, GKYL_LOWER_EDGE, &s->plocal_ext_conf, cghost);
  gkyl_skin_ghost_ranges(&s->pskin[1], &s->pghost[1], 0, GKYL_UPPER_EDGE, &s->plocal_ext_conf, cghost);

  int nc = s->nc = s->cbasis.num_basis, np = s->np = s->pbasis.num_basis;
  long cvol = s->clocal_ext.volume, pvol = s->plocal_ext.volume;
  s->bvar = mkarr(3 * nc, cvol);
  s->u = mkarr(3 * nc, cvol);
  s->grad_u = mkarr(3 * nc, cvol);
  gkyl_array_clear(s->grad_u, 0.0);
  s->theta = mkarr(nc, cvol);
  s->nu_sum = mkarr(nc, cvol);
  s->nu_vthsq_sum = mkarr(nc, cvol);
  s->alpha_0 = mkarr(nc, cvol);
  s->c_h = mkarr(nc, cvol);
  gkyl_array_clear(s->c_h, 0.0);
  s->c_0 = mkarr(nc, cvol);
  gkyl_array_clear(s->c_0, 0.0);
  s->a_perp = mkarr(3 * nc, cvol);
  s->dtb = mkarr(3 * nc, cvol);
  gkyl_array_clear(s->dtb, 0.0);
  s->xi = mkarr(3 * nc, cvol);
  gkyl_array_clear(s->xi, 0.0);
  s->grad_xi = mkarr(3 * nc, cvol);
  gkyl_array_clear(s->grad_xi, 0.0);
  s->scratch = mkarr(nc, cvol);
  s->nmoms = mkarr(18 * nc, cvol);
  s->moms = mkarr(10 * nc, cvol);
  s->rhs_moms = mkarr(18 * nc, cvol);
  s->rho = mkarr(nc, cvol);
  s->rho3 = mkarr(3 * nc, cvol);
  gkyl_array_clear(s->rho3, 0.0);
  s->grad_rho = mkarr(3 * nc, cvol);
  gkyl_array_clear(s->grad_rho, 0.0);
  s->grad_ln_n = mkarr(nc, cvol);
  gkyl_array_clear(s->grad_ln_n, 0.0);
  s->cflrate_f = mkarr(1, pvol);
  s->cflrate_m = mkarr(1, cvol);
  gkyl_array_clear(s->nu_sum, 0.0);
  gkyl_array_clear(s->nu_vthsq_sum, 0.0);
  gkyl_array_clear(s->alpha_0, 0.0);

  s->em = mkarr(8 * nc, cvol);
  struct gkyl_array *em = s->em;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 8, eval_em, &bctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, em);
  gkyl_proj_on_basis_release(proj);
  s->calc_b = gkyl_dg_calc_pkpm_bvar_new(&s->cgrid, &s->cbasis, s->pos_map, false);
  gkyl_dg_calc_pkpm_bvar_advance(s->calc_b, &s->clocal_ext, em, s->bvar);
  s->grad_b = mkarr(3 * nc, cvol);
  gkyl_array_clear(s->grad_b, 0.0);
  gkyl_dg_calc_pkpm_bvar_grad(s->calc_b, &s->clocal, s->bvar, s->grad_b);

  s->vars =
    gkyl_dg_calc_pkpm_vars_new(&s->cgrid, &s->cbasis, s->pos_map, false);
  s->mcalc = gkyl_dg_updater_moment_pkpm_new(
    &s->pgrid, &s->cbasis, &s->pbasis, &s->clocal_ext, &s->vlocal, 1.0, false, false
  );

  struct gkyl_dg_pkpm_conf_flux_surf_inp cinp = {
    .phase_grid = &s->pgrid, .conf_basis = &s->cbasis, .phase_basis = &s->pbasis, .pos_map = s->pos_map, .use_gpu = false};
  s->calc_cflux = gkyl_dg_pkpm_conf_flux_surf_inew(&cinp);
  s->nkc = gkyl_dg_pkpm_conf_flux_surf_num_surf_basis(s->calc_cflux);
  s->cflux = mkarr(2 * 8 * s->nkc, pvol);
  s->face_coef = mkarr(gkyl_dg_pkpm_conf_flux_surf_face_coef_ncomp(s->calc_cflux), s->clocal_ext.volume);
  struct gkyl_dg_pkpm_kinetic_conf_inp ainp = {
    .conf_basis = &s->cbasis,
    .phase_basis = &s->pbasis,
    .conf_range = &s->clocal,
    .phase_range = &s->plocal,
    .auxfields =
      {
        .bvar = s->bvar,
        .u = s->u,
        .xi = s->xi,
        .conf_flux_surf = s->cflux,
      },
    .pos_map = s->pos_map, .use_gpu = false,
  };
  s->eqn_a = gkyl_dg_pkpm_kinetic_conf_inew(&ainp);
  // First-harmonic stage-A terms of a magnetized species: Omega = 20 |B| (|B| ~ 1), gyroradius
  // vt/Omega = 0.05 below the cell size, Omega dt < 0.4 for the explicit gyration. Without
  // gyration (the unmagnetized limit, outside the model's ordering) the sources drive M0 and N
  // secularly in the oblique cases until they blow up.
  struct gkyl_dg_pkpm_harm_vol_inp hinp = {
    .phase_grid = &s->pgrid, .conf_basis = &s->cbasis, .phase_basis = &s->pbasis, .qbym = 20.0, .pos_map = s->pos_map, .use_gpu = false};
  s->harm = gkyl_dg_pkpm_harm_vol_inew(&hinp);
  s->hcoef = mkarr(gkyl_dg_pkpm_harm_vol_num_coef(s->harm), s->clocal_ext.volume);
  s->hwn = mkarr(gkyl_dg_pkpm_harm_vol_num_wnode(s->harm), s->clocal_ext.volume);
  int up_dirs_a[GKYL_MAX_DIM] = {0}, zero_flux_a[2 * GKYL_MAX_DIM] = {0};
  s->slvr_a = gkyl_hyper_dg_new(&s->pgrid, &s->pbasis, s->eqn_a, 1, up_dirs_a, zero_flux_a, 1, false);

  struct gkyl_dg_pkpm_vel_flux_surf_inp vinp = {
    .phase_grid = &s->pgrid, .conf_basis = &s->cbasis, .phase_basis = &s->pbasis, .use_gpu = false};
  s->calc_vflux = gkyl_dg_pkpm_vel_flux_surf_inew(&vinp);
  s->vflux = mkarr(8 * nc, pvol);
  struct gkyl_dg_pkpm_kinetic_vel_inp binp = {
    .conf_basis = &s->cbasis,
    .phase_basis = &s->pbasis,
    .conf_range = &s->clocal,
    .phase_range = &s->plocal,
    .auxfields = {.c_0 = s->c_0, .alpha_0 = s->alpha_0, .c_h = s->c_h, .a_perp = s->a_perp, .vel_flux_surf = s->vflux},
    .use_gpu = false,
  };
  s->eqn_b = gkyl_dg_pkpm_kinetic_vel_inew(&binp);
  int up_dirs_b[GKYL_MAX_DIM] = {1}, zero_flux_b[2 * GKYL_MAX_DIM] = {0};
  zero_flux_b[1] = 1; // no ghost cell below w = -w_max; the upper edge is a regular face against the ghost cell
  s->slvr_b = gkyl_hyper_dg_new(&s->pgrid, &s->pbasis, s->eqn_b, 1, up_dirs_b, zero_flux_b, 1, false);

  struct gkyl_dg_pkpm_fluid_flux_surf_inp finp = {
    .conf_grid = &s->cgrid,
    .conf_basis = &s->cbasis,
    .phase_grid = &s->pgrid,
    .pos_map = s->pos_map, .use_gpu = false,
  };
  s->calc_fflux = gkyl_dg_pkpm_fluid_flux_surf_inew(&finp);
  s->nkf = gkyl_dg_pkpm_fluid_flux_surf_num_surf_basis(s->calc_fflux);
  s->fflux = mkarr(5 * s->nkf, cvol);
  struct gkyl_dg_pkpm_fluid_inp minp = {
    .conf_basis = &s->cbasis,
    .conf_range = &s->clocal,
    .auxfields = {.vlasov_pkpm_moms = s->moms, .bvar = s->bvar, .u = s->u, .fluid_flux_surf = s->fflux},
    .pos_map = s->pos_map, .use_gpu = false,
  };
  s->eqn_m = gkyl_dg_pkpm_fluid_inew(&minp);
  int up_dirs_m[GKYL_MAX_DIM] = {0}, zero_flux_m[2 * GKYL_MAX_DIM] = {0};
  s->slvr_m = gkyl_hyper_dg_new(&s->cgrid, &s->cbasis, s->eqn_m, 1, up_dirs_m, zero_flux_m, 1, false);
}

static void
sim_release(struct sim *s)
{
  gkyl_vlasov_position_map_release(s->pos_map);
  gkyl_hyper_dg_release(s->slvr_m);
  gkyl_dg_eqn_release(s->eqn_m);
  gkyl_array_release(s->fflux);
  gkyl_dg_pkpm_fluid_flux_surf_release(s->calc_fflux);
  gkyl_hyper_dg_release(s->slvr_b);
  gkyl_dg_eqn_release(s->eqn_b);
  gkyl_array_release(s->vflux);
  gkyl_dg_pkpm_vel_flux_surf_release(s->calc_vflux);
  gkyl_hyper_dg_release(s->slvr_a);
  gkyl_dg_eqn_release(s->eqn_a);
  gkyl_dg_pkpm_harm_vol_release(s->harm);
  gkyl_array_release(s->hcoef);
  gkyl_array_release(s->hwn);
  gkyl_array_release(s->cflux);
  gkyl_array_release(s->face_coef);
  gkyl_dg_pkpm_conf_flux_surf_release(s->calc_cflux);
  gkyl_dg_updater_moment_pkpm_release(s->mcalc);
  gkyl_dg_calc_pkpm_vars_release(s->vars);
  gkyl_dg_calc_pkpm_bvar_release(s->calc_b);
  gkyl_array_release(s->cflrate_m);
  gkyl_array_release(s->cflrate_f);
  gkyl_array_release(s->rhs_moms);
  gkyl_array_release(s->moms);
  gkyl_array_release(s->nmoms);
  gkyl_array_release(s->rho);
  gkyl_array_release(s->rho3);
  gkyl_array_release(s->grad_rho);
  gkyl_array_release(s->grad_ln_n);
  gkyl_array_release(s->c_0);
  gkyl_array_release(s->a_perp);
  gkyl_array_release(s->dtb);
  gkyl_array_release(s->em);
  gkyl_array_release(s->c_h);
  gkyl_array_release(s->alpha_0);
  gkyl_array_release(s->xi);
  gkyl_array_release(s->grad_xi);
  gkyl_array_release(s->scratch);
  gkyl_array_release(s->nu_vthsq_sum);
  gkyl_array_release(s->nu_sum);
  gkyl_array_release(s->theta);
  gkyl_array_release(s->grad_u);
  gkyl_array_release(s->u);
  gkyl_array_release(s->grad_b);
  gkyl_array_release(s->bvar);
}

// Normalized moments of [F, G], the nodal primitive variables u and [xi, xi^2, p] of the
// fluid state and the kinetic coupling moments [p_par, p_perp, q_par, q_perp] = [p M2,
// p M0_G, p xi M3, p xi M1_G] of the stage (rho is the leading block of the fluid array).
static void
sim_moms(
  struct sim *s, const struct gkyl_array *f, const struct gkyl_array *M,
  const struct gkyl_range *prng, const struct gkyl_range *crng
)
{
  int nc = s->nc;
  // the coupling moments of the app (pkpm_species_coupling_moms): p [M2, M0_G, xi M3, xi M1_G, m1,
  // xi (m2/2 + n0)]
  static const int raw[7] = {1, 2, 4, 5, 9, 10, 11};
  gkyl_array_clear(s->nmoms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(s->mcalc, prng, crng, f, s->nmoms);
  gkyl_dg_calc_pkpm_vars_prim(s->vars, crng, M, s->u, s->xi);
  gkyl_array_set_offset(s->rho, 1.0, M, 0);
  for (int c = 0; c < 10; ++c) {
    if (c < 7) {
      gkyl_dg_mul_op_range(&s->cbasis, c, s->moms, 2, s->xi, raw[c], s->nmoms, crng);
    } else {
      gkyl_array_set_offset(s->scratch, 0.5, s->nmoms, (12 + c - 7) * nc);
      gkyl_array_accumulate_offset(s->scratch, 1.0, s->nmoms, (15 + c - 7) * nc);
      gkyl_dg_mul_op_range(&s->cbasis, c, s->moms, 2, s->xi, 0, s->scratch, crng);
    }
    if (c == 2 || c == 3 || c >= 7) {
      gkyl_dg_mul_op_range(&s->cbasis, 0, s->scratch, c, s->moms, 0, s->xi, crng);
      gkyl_array_set_offset_range(s->moms, 1.0, s->scratch, c * nc, crng);
    }
  }
}

// Transversality projection of the stage state (as at the start of every app stage), primitive
// variables of the stage (ghost cells included), then the right-hand sides. Returns the maximum
// CFL rate.
static double
sim_rhs(
  struct sim *s, const struct run_opts *o, struct gkyl_array *f, const struct gkyl_array *M,
  struct gkyl_array *rhs_f, struct gkyl_array *rhs_M
)
{
  int nc = s->nc;
  gkyl_dg_calc_pkpm_vars_transverse(s->vars, &s->clocal_ext, &s->plocal_ext_conf, s->bvar, f);
  sim_moms(s, f, M, &s->plocal_ext_conf, &s->clocal_ext);
  gkyl_array_set_offset(s->theta, 1.0, s->nmoms, 2 * nc);
  weak_grad(s->vars, &s->clocal, s->u, s->grad_u);
  gkyl_array_set_offset(s->rho3, 1.0, s->rho, 0);
  weak_grad(s->vars, &s->clocal, s->rho3, s->grad_rho);
  gkyl_dg_calc_pkpm_vars_grad_ln_n(s->vars, &s->clocal, s->rho, s->grad_rho, s->grad_ln_n);
  weak_grad(s->vars, &s->clocal, s->xi, s->grad_xi);
  gkyl_dg_pkpm_conf_flux_surf_face_coef(s->calc_cflux, &s->clocal, s->bvar, s->u, s->xi, s->theta, s->face_coef);

  gkyl_array_clear(s->cflrate_f, 0.0);
  gkyl_array_clear(s->cflrate_m, 0.0);
  gkyl_array_clear(rhs_f, 0.0);
  gkyl_array_clear(rhs_M, 0.0);

  // Fluid face fluxes [F_rho, F_x, F_y, F_z, F_E] of the stage.
  gkyl_array_clear(s->fflux, 0.0);
  gkyl_dg_pkpm_fluid_flux_surf_advance(
    s->calc_fflux, &s->clocal, &s->clocal_ext, &s->plocal, &s->plocal_ext, s->bvar, s->u, s->xi, s->face_coef,
    M, f, s->cflrate_m, s->fflux
  );
  if (!o->freeze_kinetic) {
    gkyl_array_clear(s->cflux, 0.0);
    gkyl_dg_pkpm_conf_flux_surf_advance(
      s->calc_cflux, &s->clocal, &s->plocal, &s->plocal_ext, s->face_coef, f, s->cflrate_f, s->cflux
    );
    expand_two_sided(s->cflux, 1, 8 * s->nkc);
    gkyl_hyper_dg_advance(s->slvr_a, &s->plocal, f, s->cflrate_f, rhs_f);
    gkyl_dg_pkpm_harm_vol_coef(s->harm, &s->clocal, s->bvar, s->grad_b, s->u, s->grad_u, s->xi, s->grad_xi,
      s->grad_ln_n, s->theta, s->dtb, s->em, s->hcoef, s->hwn);
    gkyl_dg_pkpm_harm_vol_advance(s->harm, &s->clocal, &s->plocal, s->hcoef, f, s->cflrate_f, rhs_f);
    gkyl_array_clear(s->rhs_moms, 0.0);
    gkyl_dg_updater_moment_pkpm_advance(s->mcalc, &s->plocal, &s->clocal, rhs_f, s->rhs_moms);
    // closures with the characteristic velocity edges (the switch from the previous stage)
    gkyl_dg_calc_pkpm_vars_closure(
      s->vars, &s->clocal, &s->plocal, &s->pgrid, s->nmoms, s->rhs_moms, s->bvar, s->hwn, s->theta, f, s->c_0,
      s->alpha_0, s->c_h, s->a_perp, s->vflux, s->cflrate_f
    );
    gkyl_dg_pkpm_vel_flux_surf_advance(
      s->calc_vflux, &s->clocal, &s->plocal, s->hwn, s->theta, s->bvar, s->alpha_0, s->c_h, f,
      s->vflux
    );
    gkyl_hyper_dg_advance(s->slvr_b, &s->plocal, f, s->cflrate_f, rhs_f);
  }
  if (!o->freeze_fluid) {
    gkyl_hyper_dg_advance(s->slvr_m, &s->clocal, M, s->cflrate_m, rhs_M);
  }
  double wf = 0.0, wm = 0.0;
  gkyl_array_reduce_range(&wf, s->cflrate_f, GKYL_MAX, &s->plocal);
  gkyl_array_reduce_range(&wm, s->cflrate_m, GKYL_MAX, &s->clocal);
  return fmax(wf, wm);
}

struct diag {
  double mass, mom, ke, ppar, pperp, etot, du_max, dp_max;
};

struct run_result {
  double dev_max; // max |u - u0| (|p - 1| with the fluid frozen) over the run
  double denergy_max; // max |E(t) - E(0)| of the integrated fluid energy E
  double dinternal_max; // max drift of the kinetic internal energy p_par/2 + p_perp
};

static struct diag
sim_diag(struct sim *s, const struct run_opts *o, const struct gkyl_array *f, const struct gkyl_array *M)
{
  int nc = s->nc;
  double dx = s->cgrid.dx[0];
  struct diag d = {0};
  sim_moms(s, f, M, &s->plocal, &s->clocal);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s->clocal);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&s->clocal_ext, iter.idx);
    const double *m = gkyl_array_cfetch(s->moms, cidx);
    const double *Mx = gkyl_array_cfetch(M, cidx);
    double rho = Mx[0] / sqrt(2.0), ppar = m[0] / sqrt(2.0), pperp = m[nc] / sqrt(2.0);
    double mx = Mx[nc] / sqrt(2.0), my = Mx[2 * nc] / sqrt(2.0), mz = Mx[3 * nc] / sqrt(2.0);
    d.etot += Mx[4 * nc] / sqrt(2.0) * dx;
    d.mass += rho * dx;
    d.mom += mx * dx;
    d.ke += (mx * mx + my * my + mz * mz) / rho * dx;
    d.ppar += ppar * dx;
    d.pperp += pperp * dx;
    double du = fabs(mx / rho - o->u0), dp = fabs(ppar - 1.0);
    d.du_max = du > d.du_max ? du : d.du_max;
    d.dp_max = dp > d.dp_max ? dp : d.dp_max;
  }
  return d;
}

// Runs the coupled system and returns the maximum deviation of the cell-average flow
// from u0 over the run (of the parallel pressure from 1 when the fluid is frozen) and
// the maximum drift of the integrated energy; the deviation is 1e300 if the run blew up.
static struct run_result
run(const struct run_opts *o, bool verbose)
{
  struct sim s;
  sim_init(&s, o);
  int np = s.np, nc = s.nc;
  long pvol = s.plocal_ext.volume, cvol = s.clocal_ext.volume;
  struct wave_ctx ctx = {
    .u0 = o->u0, .alpha = o->alpha, .vt = 1.0, .uy_amp = o->uy_amp, .ux_amp = o->ux_amp};

  struct gkyl_array *f = mkarr(8 * np, pvol), *f1 = mkarr(8 * np, pvol), *fnew = mkarr(8 * np, pvol);
  struct gkyl_array *M = mkarr(5 * nc, cvol), *M1 = mkarr(5 * nc, cvol), *Mnew = mkarr(5 * nc, cvol);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s.pgrid, &s.pbasis, 8, 2, eval_f, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s.plocal_ext_conf, f);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_clear(M, 0.0);
  proj = gkyl_proj_on_basis_new(&s.cgrid, &s.cbasis, 8, 5, eval_fluid, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s.clocal_ext, M);
  gkyl_proj_on_basis_release(proj);
  // Constraints M0(F) = 1 and M2(F) + 2 M0(G) = 3 exactly from the projected moments.
  gkyl_array_clear(s.nmoms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(s.mcalc, &s.plocal_ext_conf, &s.clocal_ext, f, s.nmoms);
  gkyl_array_copy(f1, f);
  gkyl_dg_calc_pkpm_vars_normalize_dist(s.vars, &s.clocal_ext, &s.plocal_ext_conf, s.nmoms, f1, f);
  periodic(f, s.pskin, s.pghost);
  periodic(M, s.cskin, s.cghost);

  double t = 0.0, du_run = 0.0, de_run = 0.0, di_run = 0.0, t_print = 0.0;
  int step = 0;
  struct diag d0 = sim_diag(&s, o, f, M);
  if (verbose) {
    printf("  Nx=%d Nv=%d u0=%g alpha=%g freeze_fluid=%d freeze_kinetic=%d\n", o->Nx, o->Nv, o->u0,
      o->alpha, o->freeze_fluid, o->freeze_kinetic);
    printf("  %8s %6s %10s %10s %10s %10s %10s %10s\n", "t", "step", "dmass", "dKE", "dppar", "dpperp",
      "max|u-u0|", "max|p-1|");
  }
  while (t < o->t_end) {
    // SSP-RK3.
    double omega = sim_rhs(&s, o, f, M, fnew, Mnew);
    double dt = fmin(1.0 / omega, o->t_end - t);
    gkyl_array_accumulate(gkyl_array_scale(fnew, dt), 1.0, f);
    gkyl_array_accumulate(gkyl_array_scale(Mnew, dt), 1.0, M);
    periodic(fnew, s.pskin, s.pghost);
    periodic(Mnew, s.cskin, s.cghost);
    gkyl_array_copy(f1, fnew);
    gkyl_array_copy(M1, Mnew);

    sim_rhs(&s, o, f1, M1, fnew, Mnew);
    gkyl_array_accumulate(gkyl_array_scale(fnew, dt), 1.0, f1);
    gkyl_array_accumulate(gkyl_array_scale(Mnew, dt), 1.0, M1);
    gkyl_array_accumulate(gkyl_array_scale(fnew, 0.25), 0.75, f);
    gkyl_array_accumulate(gkyl_array_scale(Mnew, 0.25), 0.75, M);
    periodic(fnew, s.pskin, s.pghost);
    periodic(Mnew, s.cskin, s.cghost);
    gkyl_array_copy(f1, fnew);
    gkyl_array_copy(M1, Mnew);

    sim_rhs(&s, o, f1, M1, fnew, Mnew);
    gkyl_array_accumulate(gkyl_array_scale(fnew, dt), 1.0, f1);
    gkyl_array_accumulate(gkyl_array_scale(Mnew, dt), 1.0, M1);
    gkyl_array_accumulate(gkyl_array_scale(fnew, 2.0 / 3.0), 1.0 / 3.0, f);
    gkyl_array_accumulate(gkyl_array_scale(Mnew, 2.0 / 3.0), 1.0 / 3.0, M);
    periodic(fnew, s.pskin, s.pghost);
    periodic(Mnew, s.cskin, s.cghost);
    gkyl_array_copy(f, fnew);
    gkyl_array_copy(M, Mnew);

    t += dt;
    step += 1;
    struct diag d = sim_diag(&s, o, f, M);
    // With the fluid frozen the flow is fixed, so the pressure deviation is the measure.
    double dev = o->freeze_fluid ? d.dp_max : d.du_max;
    du_run = dev > du_run ? dev : du_run;
    double de = fabs(d.etot - d0.etot);
    de_run = de > de_run ? de : de_run;
    double di = fabs((0.5 * d.ppar + d.pperp) - (0.5 * d0.ppar + d0.pperp));
    di_run = di > di_run ? di : di_run;
    if (verbose && (t >= t_print || t >= o->t_end)) {
      printf("  %8.4f %6d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n", t, step, d.mass - d0.mass,
        d.ke - d0.ke, d.ppar - d0.ppar, d.pperp - d0.pperp, d.du_max, d.dp_max);
      t_print += o->t_end / 10.0;
    }
    if (!(du_run < 10.0)) {
      du_run = 1e300;
      break;
    }
  }
  if (verbose) {
    // Cell averages and slopes of rho and u_x at the end: a slope much larger than the
    // neighbor difference / 2 sqrt(3/2) is a sub-cell sawtooth.
    sim_moms(&s, f, M, &s.plocal, &s.clocal);
    printf("  final rho avg / slope, ux avg / slope per cell:\n");
    struct gkyl_range_iter it;
    gkyl_range_iter_init(&it, &s.clocal);
    while (gkyl_range_iter_next(&it)) {
      long cidx = gkyl_range_idx(&s.clocal_ext, it.idx);
      const double *m = gkyl_array_cfetch(s.moms, cidx), *u = gkyl_array_cfetch(s.u, cidx);
      printf("   %3d  %9.5f %9.5f   %9.5f %9.5f\n", it.idx[0], m[0] / sqrt(2.0), m[1], u[0] / sqrt(2.0), u[1]);
    }
  }

  gkyl_array_release(Mnew);
  gkyl_array_release(M1);
  gkyl_array_release(M);
  gkyl_array_release(fnew);
  gkyl_array_release(f1);
  gkyl_array_release(f);
  sim_release(&s);
  return (struct run_result){.dev_max = du_run, .denergy_max = de_run, .dinternal_max = di_run};
}

static bool
verbose(void)
{
  return test_verbose_level_ >= 3;
}

// Coupled runs: the flow stays bounded (the collisionless temperature perturbation phase
// mixes into a sound wave of modest amplitude) and the energy drift is small.
static void
check_coupled(const struct run_opts *o)
{
  struct run_result r = run(o, verbose());
  TEST_CHECK(r.dev_max < 0.2);
  TEST_MSG("max |u - u0| over the run: %g", r.dev_max);
  TEST_CHECK(r.denergy_max < 1e-10);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}

void
test_coupled_u1()
{
  struct run_opts o = {.u0 = 1.0, .alpha = 0.2, .t_end = 0.6, .Nx = 16, .Nv = 16};
  check_coupled(&o);
}

void
test_coupled_u0()
{
  struct run_opts o = {.u0 = 0.0, .alpha = 0.2, .t_end = 0.6, .Nx = 16, .Nv = 16};
  check_coupled(&o);
}

// Oblique, x-dependent static field (div b != 0): the mirror force exchanges parallel
// and perpendicular pressure in a plasma initially at rest; the total energy must be
// conserved up to the (small) dissipation of the fluxes.
void
test_mirror_oblique()
{
  struct run_opts o = {.u0 = 0.0, .alpha = 0.0, .t_end = 2.0, .Nx = 16, .Nv = 16,
    .b_theta = 0.3, .b_perp = 0.5};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.dev_max < 1.0);
  TEST_MSG("max |u - u0| over the run: %g", r.dev_max);
  TEST_CHECK(r.denergy_max < 1e-10);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}

// Oblique varying field with a flowing density pulse: streaming along b, mirror
// forces and the K = bb:grad(u) work all act at once.
void
test_mirror_oblique_flow()
{
  struct run_opts o = {.u0 = 0.5, .alpha = 0.2, .t_end = 2.0, .Nx = 16, .Nv = 16,
    .b_theta = 0.3, .b_perp = 0.5};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.dev_max < 1.0);
  TEST_MSG("max |u - u0| over the run: %g", r.dev_max);
  TEST_CHECK(r.denergy_max < 1e-10);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}

// Same field, kinetic update alone (u = u0 exactly) with a uniform density: without
// compression work the internal energy p_par/2 + p_perp must be conserved by the
// mirror-force pairing. (In the factored update the density lives in the frozen fluid
// array, so a non-uniform frozen density with a flow would not be a consistent state:
// the normalized distribution would stream through a density it cannot change.)
void
test_mirror_oblique_frozen_fluid()
{
  struct run_opts o = {.u0 = 0.5, .alpha = 0.0, .t_end = 2.0, .Nx = 16, .Nv = 16,
    .b_theta = 0.3, .b_perp = 0.5, .freeze_fluid = true};
  struct run_result r = run(&o, verbose());
  // The momentum is not evolved here, so only the internal energy is meaningful.
  TEST_CHECK(r.dinternal_max < 1e-5);
  TEST_MSG("max internal energy drift over the run: %g", r.dinternal_max);
}

// Sheared perpendicular drift through an oblique uniform field: K = b_x b_y du_y/dx
// exchanges perpendicular kinetic energy with the pressures through P_xy; the total
// energy must be conserved.
static void
check_shear(int Nx)
{
  struct run_opts o = {.u0 = 0.0, .alpha = 0.0, .t_end = 2.0, .Nx = Nx, .Nv = 16,
    .b_theta = 0.3, .b_perp = 0.0, .uy_amp = 0.5};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.denergy_max < 1e-10);
  TEST_MSG("Nx=%d: max energy drift over the run: %g", Nx, r.denergy_max);
}
void test_shear_nx16() { check_shear(16); }

// Same shear with a nearly perpendicular field: the streaming speed along x is
// small for every particle, so the slope modes get little upwind damping.
// Compressive flow varying along x across a nearly perpendicular field: the upwind
// switch w* = -u_x/b_x of the face flux sits inside different velocity cells at the
// two faces of a cell, so its quadrature error acts on the slope modes.
static void
check_compress(double theta, double amp, bool b_along_x)
{
  struct run_opts o = {.u0 = 0.0, .alpha = 0.0, .t_end = 2.0, .Nx = 16, .Nv = 16,
    .b_theta = b_along_x ? 0.0 : theta, .b_perp = 0.0, .ux_amp = amp};
  if (b_along_x) {
    o.b_theta = 0.0;
  }
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.denergy_max < 1e-10);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}
void test_compress_quasi_perp() { check_compress(0.0871, 0.1, false); }
void test_compress_parallel() { check_compress(0.0, 0.1, true); }

void
test_shear_quasi_perp()
{
  struct run_opts o = {.u0 = 0.0, .alpha = 0.0, .t_end = 2.0, .Nx = 16, .Nv = 16,
    .b_theta = 0.0871, .b_perp = 0.0, .uy_amp = 0.5};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.denergy_max < 1e-10);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}

// Nearly perpendicular field (b_x = 0.087, as in the kinetic Alfven wave input) with a
// slow flow comparable to the parallel streaming speed along x: a slope-mode check.
void
test_quasi_perp()
{
  struct run_opts o = {.u0 = 0.05, .alpha = 0.2, .t_end = 20.0, .Nx = 16, .Nv = 16,
    .b_theta = 0.0871, .b_perp = 0.0};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.dev_max < 1.0);
  TEST_MSG("max |u - u0| over the run: %g", r.dev_max);
  TEST_CHECK(r.denergy_max < 1e-10);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}
void test_shear_nx32() { check_shear(32); }

void
test_mirror_oblique_flow_nx32()
{
  struct run_opts o = {.u0 = 0.5, .alpha = 0.2, .t_end = 2.0, .Nx = 32, .Nv = 16,
    .b_theta = 0.3, .b_perp = 0.5};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.denergy_max < 1e-10);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}

void
test_mirror_oblique_flow_nx64()
{
  struct run_opts o = {.u0 = 0.5, .alpha = 0.2, .t_end = 2.0, .Nx = 64, .Nv = 16,
    .b_theta = 0.3, .b_perp = 0.5};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.denergy_max < 1e-10);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}

// Finer grid: the inconsistent coupling grew at a rate ~ u/dx, so this run is the
// sharper instability check.
void
test_coupled_u1_nx64()
{
  struct run_opts o = {.u0 = 1.0, .alpha = 0.2, .t_end = 0.4, .Nx = 64, .Nv = 16};
  check_coupled(&o);
}

void
test_frozen_fluid_u1()
{
  struct run_opts o = {
    .u0 = 1.0, .alpha = 0.2, .t_end = 0.6, .Nx = 16, .Nv = 16, .freeze_fluid = true};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.dev_max < 1.0);
  TEST_MSG("max |p - 1| over the run: %g", r.dev_max);
}

void
test_frozen_kinetic_u1()
{
  struct run_opts o = {
    .u0 = 1.0, .alpha = 0.2, .t_end = 0.6, .Nx = 16, .Nv = 16, .freeze_kinetic = true};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.dev_max < 1.0);
  TEST_MSG("max |u - u0| over the run: %g", r.dev_max);
}

// Discrete steady state on a steep density profile (plan section 11.11, step 2): a uniform
// Maxwellian [F, G] with u = 0 and uniform xi on n = 1 + alpha sin(pi x) is a steady state of
// the normalized kinetic equation (the density-gradient source is cancelled pointwise in w by
// the closure force alpha_0 = xi b.grad(ln n)). The discrete right-hand side of one evaluation
// is measured by mode (cell average, slope) and by moment; the constrained moments vanish by
// construction, the residual of the other coefficients is the source/closure mismatch.
struct steady_res {
  double avg, slope, m3_total;
};

static struct steady_res
steady_gradient_residual(int Nx, int Nv, double alpha, bool verbose)
{
  struct run_opts o = {.u0 = 0.0, .alpha = alpha, .t_end = 0.0, .Nx = Nx, .Nv = Nv};
  struct sim s;
  sim_init(&s, &o);
  int np = s.np, nc = s.nc;
  long pvol = s.plocal_ext.volume, cvol = s.clocal_ext.volume;
  struct wave_ctx ctx = {.u0 = 0.0, .alpha = alpha, .vt = 1.0};
  struct gkyl_array *f = mkarr(8 * np, pvol), *f1 = mkarr(8 * np, pvol), *rhs = mkarr(8 * np, pvol);
  struct gkyl_array *M = mkarr(5 * nc, cvol), *rhs_M = mkarr(5 * nc, cvol);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s.pgrid, &s.pbasis, 8, 2, eval_f, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s.plocal_ext_conf, f);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_clear(M, 0.0);
  proj = gkyl_proj_on_basis_new(&s.cgrid, &s.cbasis, 8, 5, eval_fluid_isothermal, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s.clocal_ext, M);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_clear(s.nmoms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(s.mcalc, &s.plocal_ext_conf, &s.clocal_ext, f, s.nmoms);
  gkyl_array_copy(f1, f);
  gkyl_dg_calc_pkpm_vars_normalize_dist(s.vars, &s.clocal_ext, &s.plocal_ext_conf, s.nmoms, f1, f);
  periodic(f, s.pskin, s.pghost);
  periodic(M, s.cskin, s.cghost);

  sim_rhs(&s, &o, f, M, rhs, rhs_M);
  // sim_rhs leaves the moments of the configuration pass (stage A) in rhs_moms
  struct gkyl_array *momsA = mkarr(18 * nc, cvol);
  gkyl_array_copy(momsA, s.rhs_moms);
  gkyl_array_clear(s.rhs_moms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(s.mcalc, &s.plocal, &s.clocal, rhs, s.rhs_moms);
  double m1A = 0.0, m3A = 0.0, m3T = 0.0, a0 = 0.0, gln = 0.0, ratioA = 0.0;
  {
    struct gkyl_range_iter it;
    gkyl_range_iter_init(&it, &s.clocal);
    while (gkyl_range_iter_next(&it)) {
      long loc = gkyl_range_idx(&s.clocal_ext, it.idx);
      const double *mA = gkyl_array_cfetch(momsA, loc), *mT = gkyl_array_cfetch(s.rhs_moms, loc);
      const double *al = gkyl_array_cfetch(s.alpha_0, loc), *g = gkyl_array_cfetch(s.grad_ln_n, loc);
      if (fabs(mA[3 * nc]) > m1A) {
        m1A = fabs(mA[3 * nc]);
        ratioA = mA[4 * nc] / mA[3 * nc];
      }
      m3A = fmax(m3A, fabs(mA[4 * nc]));
      m3T = fmax(m3T, fabs(mT[4 * nc]));
      a0 = fmax(a0, fabs(al[0]));
      gln = fmax(gln, fabs(g[0]));
    }
  }
  gkyl_array_release(momsA);
  // Decomposition of the configuration pass: volume kernel alone versus the full pass
  // (volume + stored face fluxes) for the same state.
  {
    struct gkyl_array *vol = mkarr(8 * np, pvol), *rhsA = mkarr(8 * np, pvol), *mv = mkarr(18 * nc, cvol), *ma = mkarr(18 * nc, cvol);
    gkyl_array_clear(vol, 0.0);
    gkyl_array_clear(rhsA, 0.0);
    gkyl_array_clear(s.cflrate_f, 0.0);
    gkyl_array_clear(s.cflux, 0.0);
    gkyl_dg_pkpm_conf_flux_surf_face_coef(s.calc_cflux, &s.clocal, s.bvar, s.u, s.xi, s.theta, s.face_coef);
    gkyl_dg_pkpm_conf_flux_surf_advance(s.calc_cflux, &s.clocal, &s.plocal, &s.plocal_ext, s.face_coef, f, s.cflrate_f, s.cflux);
    expand_two_sided(s.cflux, 1, 8 * s.nkc);
    gkyl_hyper_dg_advance(s.slvr_a, &s.plocal, f, s.cflrate_f, rhsA);
    struct gkyl_range_iter it;
    gkyl_range_iter_init(&it, &s.plocal);
    while (gkyl_range_iter_next(&it)) {
      long pidx = gkyl_range_idx(&s.plocal_ext, it.idx), cidx = gkyl_range_idx(&s.clocal_ext, it.idx);
      double xc[2];
      gkyl_rect_grid_cell_center(&s.pgrid, it.idx, xc);
      pkpm_kinetic_conf_vol_1x1v_ser_p1(xc, s.pgrid.dx, gkyl_array_cfetch(s.pos_map->jacob_pos, cidx),
        gkyl_array_cfetch(s.bvar, cidx), gkyl_array_cfetch(s.u, cidx), gkyl_array_cfetch(s.xi, cidx),
        gkyl_array_cfetch(f, pidx), gkyl_array_fetch(vol, pidx));
    }
    gkyl_array_clear(mv, 0.0);
    gkyl_array_clear(ma, 0.0);
    gkyl_dg_updater_moment_pkpm_advance(s.mcalc, &s.plocal, &s.clocal, vol, mv);
    gkyl_dg_updater_moment_pkpm_advance(s.mcalc, &s.plocal, &s.clocal, rhsA, ma);
    double m1v = 0.0, m3v = 0.0, m1a = 0.0, m3a = 0.0, m1s = 0.0, m3s = 0.0;
    gkyl_range_iter_init(&it, &s.clocal);
    while (gkyl_range_iter_next(&it)) {
      long loc = gkyl_range_idx(&s.clocal_ext, it.idx);
      const double *a = gkyl_array_cfetch(mv, loc), *b = gkyl_array_cfetch(ma, loc);
      m1v = fmax(m1v, fabs(a[3 * nc])); m3v = fmax(m3v, fabs(a[4 * nc]));
      m1a = fmax(m1a, fabs(b[3 * nc])); m3a = fmax(m3a, fabs(b[4 * nc]));
      m1s = fmax(m1s, fabs(b[3 * nc] - a[3 * nc])); m3s = fmax(m3s, fabs(b[4 * nc] - a[4 * nc]));
    }
    if (verbose) {
      printf("      volume kernel: M1 %.3e M3 %.3e | surface part: M1 %.3e M3 %.3e | pass: M1 %.3e M3 %.3e\n", m1v, m3v, m1s, m3s, m1a, m3a);
    }
    gkyl_array_release(ma); gkyl_array_release(mv); gkyl_array_release(rhsA); gkyl_array_release(vol);
  }

  // Scale: the size of the density-gradient source alone, xi w |d ln n/dx| F at w = 2.
  double dx = s.cgrid.dx[0], gmax = alpha * M_PI / (1.0 - alpha);
  double scale = 2.0 * gmax * exp(-2.0) / sqrt(2.0 * M_PI);
  double res_avg = 0.0, res_slope = 0.0, mom_cons = 0.0, mom_other = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *r = gkyl_array_cfetch(rhs, gkyl_range_idx(&s.plocal_ext, iter.idx));
    for (int c = 0; c < 2; ++c) {
      for (int k = 0; k < np; ++k) {
        // hybrid 1x1v p1: basis 1 and 3 (and 5) carry the x-slope
        bool slope = k == 1 || k == 3 || k == 5;
        double v = fabs(r[c * np + k]) / sqrt(2.0); // phase-space psi_0 = 1/2, conf psi_0 = 1/sqrt(2)
        if (slope) res_slope = fmax(res_slope, v);
        else res_avg = fmax(res_avg, v);
      }
    }
  }
  gkyl_range_iter_init(&iter, &s.clocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *m = gkyl_array_cfetch(s.rhs_moms, gkyl_range_idx(&s.clocal_ext, iter.idx));
    for (int k = 0; k < nc; ++k) {
      mom_cons = fmax(mom_cons, fmax(fabs(m[k]), fmax(fabs(m[3 * nc + k]), fabs(m[nc + k] + 2.0 * m[2 * nc + k]))));
      mom_other = fmax(mom_other, fmax(fabs(m[4 * nc + k]), fabs(m[5 * nc + k])));
    }
  }
  if (verbose) {
    printf("  Nx=%d Nv=%d alpha=%g (max |dlnn/dx| %.2f): residual/scale avg %.3e slope %.3e | stage A: M1 %.3e M3 %.3e (M3/M1 %.3f, Maxwellian 3) | total M3 %.3e | alpha_0 %.3e vs xi grad ln n %.3e\n",
      Nx, Nv, alpha, gmax, res_avg / scale, res_slope / scale, m1A, m3A, ratioA, m3T, a0, gln);
  }
  TEST_CHECK(mom_cons < 1e-11);
  TEST_MSG("Nx=%d Nv=%d: constrained moments of the right-hand side %g", Nx, Nv, mom_cons);
  struct steady_res r = {.avg = res_avg / scale, .slope = res_slope / scale, .m3_total = m3T};
  gkyl_array_release(rhs_M);
  gkyl_array_release(M);
  gkyl_array_release(rhs);
  gkyl_array_release(f1);
  gkyl_array_release(f);
  sim_release(&s);
  return r;
}

static void
test_steady_gradient()
{
  bool verbose = test_verbose_level_ > 2;
  // The residual of the exact steady state is truncation error only. The closure force
  // alpha_0 = xi grad(ln n) cancels the density-gradient source up to the error of the w
  // discretization, so the residual is grad(ln n)(x) E(w) with E = O(dv^2): the cell-average
  // residual converges with dv^2, the slope residual is the in-cell slope of grad(ln n)
  // (O(dx)) times E, and the unconstrained third moment is small at every resolution.
  struct steady_res r16 = steady_gradient_residual(16, 16, 0.6, verbose);
  struct steady_res r32 = steady_gradient_residual(16, 32, 0.6, verbose);
  struct steady_res r64 = steady_gradient_residual(16, 64, 0.6, verbose);
  TEST_CHECK(r16.avg / r32.avg > 3.0 && r32.avg / r64.avg > 3.0);
  TEST_MSG("average residual vs Nv: %g %g %g", r16.avg, r32.avg, r64.avg);
  TEST_CHECK(r64.avg < 5e-3);
  TEST_CHECK(r16.slope / r32.slope > 3.0 && r32.slope / r64.slope > 3.0);
  TEST_MSG("slope residual vs Nv: %g %g %g", r16.slope, r32.slope, r64.slope);
  struct steady_res x8 = steady_gradient_residual(8, 32, 0.6, verbose);
  struct steady_res x32 = steady_gradient_residual(32, 32, 0.6, verbose);
  TEST_CHECK(x8.slope / r32.slope > 1.7 && r32.slope / x32.slope > 1.7);
  TEST_MSG("slope residual vs Nx: %g %g %g", x8.slope, r32.slope, x32.slope);
  TEST_CHECK(r32.m3_total < 1e-4 && x32.m3_total < 1e-4);
  TEST_MSG("third moment of the total right-hand side: %g %g", r32.m3_total, x32.m3_total);
}

TEST_LIST = {
  {"steady_gradient", test_steady_gradient},
  {"coupled_u1", test_coupled_u1},
  {"coupled_u0", test_coupled_u0},
  {"coupled_u1_nx64", test_coupled_u1_nx64},
  {"mirror_oblique", test_mirror_oblique},
  {"mirror_oblique_flow", test_mirror_oblique_flow},
  {"mirror_oblique_frozen_fluid", test_mirror_oblique_frozen_fluid},
  {"shear_nx16", test_shear_nx16},
  {"quasi_perp", test_quasi_perp},
  {"shear_quasi_perp", test_shear_quasi_perp},
  {"compress_quasi_perp", test_compress_quasi_perp},
  {"compress_parallel", test_compress_parallel},
  {"shear_nx32", test_shear_nx32},
  {"mirror_oblique_flow_nx32", test_mirror_oblique_flow_nx32},
  {"mirror_oblique_flow_nx64", test_mirror_oblique_flow_nx64},
  {"frozen_fluid_u1", test_frozen_fluid_u1},
  {"frozen_kinetic_u1", test_frozen_kinetic_u1},
  {NULL, NULL},
};
