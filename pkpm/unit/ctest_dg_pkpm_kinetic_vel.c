// Tests of the velocity-space pass (stage B) of the refactored PKPM kinetic update of the eight
// blocks [F, G, M0_perp, N_perp]: the v_par flux (dg_pkpm_vel_flux_surf, HLL of the coupled w
// flux) against an independent nodal evaluation, the stage-B volume kernel against quadrature,
// the invariants of the complete two-pass right-hand side (M0, M1 and M2(F) + 2 M0(G) of
// rhs_A + rhs_B vanish per cell and mode, and I_perp int rhs(M0_perp) dw at the closure nodes, also
// off the constraint manifold, int M0_perp dw != 0), the characteristic velocity edges, and the
// transversality projection of the first-harmonic blocks.
#include <math.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_pkpm_conf_flux_surf.h>
#include <gkyl_dg_pkpm_harm_vol.h>
#include <gkyl_dg_pkpm_kinetic_conf.h>
#include <gkyl_dg_pkpm_kinetic_vel.h>
#include <gkyl_dg_pkpm_vel_flux_surf.h>
#include <gkyl_dg_updater_moment_pkpm.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_pkpm_harm_coeffs.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_ten_moment.h>

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
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  for (int i = 0; i < 8; ++i) {
    fout[i] = 0.0;
  }
  fout[3] = 1.0 + 0.4 * sin(x) * cos(y);
  fout[4] = 0.6 * cos(x) * sin(z) + 0.5 * sin(y);
  fout[5] = -0.5 + 0.3 * sin(y) * cos(z) + 0.2 * sin(x);
}

static void
eval_u(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  fout[0] = 0.8 * sin(x) + 0.3 * cos(y);
  fout[1] = -0.7 * cos(y) * sin(z) + 0.2 * sin(x);
  fout[2] = 0.5 * sin(z) - 0.4 * cos(x) * cos(y);
}

static void
eval_conf_scalars(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  fout[0] = 0.7 + 0.2 * cos(x) * sin(y) + 0.1 * z; // theta
  fout[1] = 0.3 + 0.1 * sin(x); // nu_sum
  fout[2] = 0.2 + 0.05 * cos(y) + 0.03 * sin(z); // nu_vthsq_sum
  fout[3] = 0.4 * sin(x) - 0.3 * cos(y) * sin(z); // a stand-in alpha_0 = a_P/xi
  fout[4] = 0.2 * cos(x) + 0.1 * sin(y) * cos(z); // a stand-in heating rate c
  fout[5] = 0.3 * sin(x) * sin(y) - 0.1 * cos(z); // a stand-in c_0
  fout[6] = 0.2 * cos(y) - 0.1 * sin(x); // a stand-in a_perp (x, y, z)
  fout[7] = 0.15 * sin(x) * cos(z);
  fout[8] = -0.25 + 0.1 * cos(x + y);
}

// d_t B/|B| of a time-dependent field.
static void
eval_dtb(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  fout[0] = 0.3 * cos(x) * sin(y);
  fout[1] = -0.2 + 0.4 * sin(x - z);
  fout[2] = 0.25 * cos(y) + 0.1 * sin(x);
}

static void
eval_gln(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  for (int d = 0; d < cdim; ++d) {
    fout[d] = 0.3 * cos(x + d) * sin(y) + 0.1 * (d + 1) * cos(z);
  }
}

static const double qbym_test = 0.7; // charge to mass ratio of the gyration term

// Thermal speed [xi, xi^2, p] (positive; p is not used by the kinetic kernels).
static void
eval_xi(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  double xi = 1.0 + 0.2 * cos(x) * sin(y) + 0.05 * z;
  fout[0] = xi;
  fout[1] = xi * xi;
  fout[2] = 0.9 * xi * xi;
}

static void
eval_f(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  double v = xn[cdim];
  double n = 1.0 + 0.3 * sin(x) * cos(y) + 0.1 * cos(z);
  double vt = 1.0 + 0.2 * cos(x);
  double F0 = n / sqrt(2.0 * M_PI * vt * vt) * exp(-v * v / (2.0 * vt * vt));
  fout[0] = F0 * (1.0 + 0.1 * v); // a deliberately nonzero first moment
  fout[1] = 0.5 * vt * vt * F0 * (1.0 - 0.1 * v);
  // generic first-harmonic blocks (int M0_perp dw != 0: off the constraint manifold); the setup
  // projects them perpendicular to b as the app does
  fout[2] = F0 * (0.3 * sin(x) + 0.1 * v);
  fout[3] = F0 * (-0.2 * cos(y) + 0.15 * v * v - 0.1);
  fout[4] = F0 * (0.25 * cos(z) - 0.05 * v);
  fout[5] = 0.5 * vt * vt * F0 * (0.2 * cos(x) - 0.1 * v);
  fout[6] = 0.5 * vt * vt * F0 * (0.1 + 0.2 * sin(y) * v);
  fout[7] = 0.5 * vt * vt * F0 * (-0.3 * sin(x + z) + 0.05 * v * v);
}

// The v_par flux samples the configuration face at 5 Gauss-Legendre points per direction (exact
// for its degree-9 integrand).
static const double gl5_x[5] = {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831,
  0.9061798459386640};
static const double gl5_w[5] = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
  0.4786286704993665, 0.2369268850561891};
static const double gl4_x[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
  0.8611363115940526};
static const double gl4_w[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
  0.3478548451374538};

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

struct test_setup {
  struct gkyl_vlasov_position_map *pos_map; // identity configuration-space map
  int cdim, pdim;
  struct gkyl_rect_grid cgrid, pgrid;
  struct gkyl_basis cbasis, pbasis;
  struct gkyl_range clocal, clocal_ext, plocal, plocal_ext, plocal_ext_conf, vlocal, vlocal_ext;
  struct gkyl_array *bvar, *grad_b, *u, *grad_u, *grad_ln_n, *theta, *nu_sum, *nu_vthsq_sum, *alpha_0, *c_h, *c_0, *f;
  struct gkyl_array *xi, *grad_xi, *a_perp, *em, *dtb;
  struct gkyl_array *hcoef; // stored first-harmonic operator of the fields
  struct gkyl_array *hwn; // raw w-flux values at the w-face nodes
  struct gkyl_dg_calc_pkpm_bvar *calc_b;
};

// The stored first-harmonic operator and the raw w-flux node values of the current fields
// (gkyl_dg_pkpm_harm_vol_coef).
static void
build_coef(struct test_setup *s)
{
  struct gkyl_dg_pkpm_harm_vol_inp hinp = {.phase_grid = &s->pgrid, .conf_basis = &s->cbasis,
    .phase_basis = &s->pbasis, .qbym = qbym_test, .pos_map = s->pos_map, .use_gpu = false};
  struct gkyl_dg_pkpm_harm_vol *harm = gkyl_dg_pkpm_harm_vol_inew(&hinp);
  if (!s->hcoef) {
    s->hcoef = mkarr(gkyl_dg_pkpm_harm_vol_num_coef(harm), s->clocal_ext.volume);
    s->hwn = mkarr(gkyl_dg_pkpm_harm_vol_num_wnode(harm), s->clocal_ext.volume);
  }
  gkyl_dg_pkpm_harm_vol_coef(harm, &s->clocal, s->bvar, s->grad_b, s->u, s->grad_u, s->xi, s->grad_xi,
    s->grad_ln_n, s->theta, s->dtb, s->em, s->hcoef, s->hwn);
  gkyl_dg_pkpm_harm_vol_release(harm);
}

// Weak gradients of u and of [xi, xi^2, p] on the interior cells (both projected on the
// extended range); only the first block of grad_xi is used by the kernels.
static void
compute_grad_u(struct test_setup *s)
{
  int cdim = s->cdim, nc = s->cbasis.num_basis;
  s->grad_u = mkarr(cdim * 3 * nc, s->clocal_ext.volume);
  s->grad_xi = mkarr(cdim * 3 * nc, s->clocal_ext.volume);
  gkyl_array_clear(s->grad_u, 0.0);
  gkyl_array_clear(s->grad_xi, 0.0);
  struct gkyl_dg_calc_pkpm_vars *vars =
    gkyl_dg_calc_pkpm_vars_new(&s->cgrid, &s->cbasis, s->pos_map, false);
  weak_grad(vars, &s->clocal, s->u, s->grad_u);
  weak_grad(vars, &s->clocal, s->xi, s->grad_xi);
  gkyl_dg_calc_pkpm_vars_release(vars);
}

static void
setup(struct test_setup *s, int cdim)
{
  int pdim = cdim + 1;
  s->cdim = cdim;
  s->pdim = pdim;
  double lower[4] = {-M_PI, -M_PI, -M_PI, -4.0}, upper[4] = {M_PI, M_PI, M_PI, 4.0};
  int cells[4] = {6, 4, 3, 8};
  int ccells[3] = {6, 4, 3};
  int cghost[3] = {1, 1, 1}, vghost[1] = {0};
  if (cdim < 3) {
    lower[cdim] = -4.0;
    upper[cdim] = 4.0;
    cells[cdim] = 8;
  }
  gkyl_rect_grid_init(&s->cgrid, cdim, lower, upper, ccells);
  gkyl_rect_grid_init(&s->pgrid, pdim, lower, upper, cells);
  gkyl_cart_modal_serendip(&s->cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&s->pbasis, cdim, 1);
  gkyl_create_grid_ranges(&s->cgrid, cghost, &s->clocal_ext, &s->clocal);
  s->pos_map = identity_pos_map(&s->cgrid, &s->clocal, &s->clocal_ext, &s->cbasis);
  phase_ranges(&s->pgrid, cdim, &s->plocal_ext, &s->plocal, &s->plocal_ext_conf);
  struct gkyl_rect_grid vgrid;
  gkyl_rect_grid_init(&vgrid, 1, &lower[cdim], &upper[cdim], &cells[cdim]);
  gkyl_create_grid_ranges(&vgrid, vghost, &s->vlocal_ext, &s->vlocal);

  int nc = s->cbasis.num_basis, np = s->pbasis.num_basis;
  s->em = mkarr(8 * nc, s->clocal_ext.volume);
  struct gkyl_array *em = s->em;
  struct gkyl_array *scal = mkarr(9 * nc, s->clocal_ext.volume);
  s->a_perp = mkarr(3 * nc, s->clocal_ext.volume);
  s->dtb = mkarr(3 * nc, s->clocal_ext.volume);
  s->bvar = mkarr(3 * nc, s->clocal_ext.volume);
  s->u = mkarr(3 * nc, s->clocal_ext.volume);
  s->xi = mkarr(3 * nc, s->clocal_ext.volume);
  s->theta = mkarr(nc, s->clocal_ext.volume);
  s->nu_sum = mkarr(nc, s->clocal_ext.volume);
  s->nu_vthsq_sum = mkarr(nc, s->clocal_ext.volume);
  s->alpha_0 = mkarr(nc, s->clocal_ext.volume);
  s->c_h = mkarr(nc, s->clocal_ext.volume);
  s->c_0 = mkarr(nc, s->clocal_ext.volume);
  s->grad_ln_n = mkarr(cdim * nc, s->clocal_ext.volume);
  s->f = mkarr(8 * np, s->plocal_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 8, eval_em, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, em);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 3, eval_u, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->u);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 9, eval_conf_scalars, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, scal);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_set_offset(s->theta, 1.0, scal, 0 * nc);
  gkyl_array_set_offset(s->nu_sum, 1.0, scal, 1 * nc);
  gkyl_array_set_offset(s->nu_vthsq_sum, 1.0, scal, 2 * nc);
  gkyl_array_set_offset(s->alpha_0, 1.0, scal, 3 * nc);
  gkyl_array_set_offset(s->c_h, 1.0, scal, 4 * nc);
  gkyl_array_set_offset(s->c_0, 1.0, scal, 5 * nc);
  for (int a = 0; a < 3; ++a) {
    gkyl_array_set_offset_range(s->a_perp, 1.0, scal, (6 + a) * nc, &s->clocal_ext);
  }
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 3, eval_dtb, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->dtb);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, cdim, eval_gln, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->grad_ln_n);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 3, eval_xi, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->xi);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->pgrid, &s->pbasis, 8, 8, eval_f, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->plocal_ext_conf, s->f);
  gkyl_proj_on_basis_release(proj);

  s->calc_b = gkyl_dg_calc_pkpm_bvar_new(&s->cgrid, &s->cbasis, s->pos_map, false);
  gkyl_dg_calc_pkpm_bvar_advance(s->calc_b, &s->clocal_ext, em, s->bvar);
  s->grad_b = mkarr(cdim * 3 * nc, s->clocal_ext.volume);
  gkyl_array_clear(s->grad_b, 0.0);
  gkyl_dg_calc_pkpm_bvar_grad(s->calc_b, &s->clocal, s->bvar, s->grad_b);
  gkyl_array_release(scal);
  compute_grad_u(s);
  // the stage state is transverse at the nodes of b (as at the start of every app stage)
  struct gkyl_dg_calc_pkpm_vars *vars =
    gkyl_dg_calc_pkpm_vars_new(&s->cgrid, &s->cbasis, s->pos_map, false);
  gkyl_dg_calc_pkpm_vars_transverse(vars, &s->clocal_ext, &s->plocal_ext_conf, s->bvar, s->f);
  gkyl_dg_calc_pkpm_vars_release(vars);
  s->hcoef = s->hwn = 0;
  build_coef(s);
}

static void
teardown(struct test_setup *s)
{
  gkyl_vlasov_position_map_release(s->pos_map);
  gkyl_dg_calc_pkpm_bvar_release(s->calc_b);
  gkyl_array_release(s->grad_b);
  gkyl_array_release(s->grad_u);
  gkyl_array_release(s->grad_xi);
  gkyl_array_release(s->xi);
  gkyl_array_release(s->f);
  gkyl_array_release(s->alpha_0);
  gkyl_array_release(s->c_h);
  gkyl_array_release(s->c_0);
  gkyl_array_release(s->nu_vthsq_sum);
  gkyl_array_release(s->nu_sum);
  gkyl_array_release(s->theta);
  gkyl_array_release(s->u);
  gkyl_array_release(s->bvar);
  gkyl_array_release(s->a_perp);
  gkyl_array_release(s->em);
  gkyl_array_release(s->dtb);
  gkyl_array_release(s->grad_ln_n);
  gkyl_array_release(s->hcoef);
  gkyl_array_release(s->hwn);
}

// Pointwise first-harmonic coefficients at reference point z of cell cidx (the definitions of
// gkyl_pkpm_harm_coeffs.h from the fields, not the stored projections).
static void
conf_coeffs(struct test_setup *s, long cidx, const double *z, struct pkpm_harm_coeffs *hc)
{
  int cdim = s->cdim, nc = s->cbasis.num_basis;
  const double *b = gkyl_array_cfetch(s->bvar, cidx), *u = gkyl_array_cfetch(s->u, cidx);
  const double *gu = gkyl_array_cfetch(s->grad_u, cidx), *gb = gkyl_array_cfetch(s->grad_b, cidx);
  const double *gxi = gkyl_array_cfetch(s->grad_xi, cidx), *dtb = gkyl_array_cfetch(s->dtb, cidx);
  double bq[3], uq[3], dq[3], guq[9], gbq[9], gxq[3];
  for (int i = 0; i < 3; ++i) {
    bq[i] = s->cbasis.eval_expand(z, &b[i * nc]);
    uq[i] = s->cbasis.eval_expand(z, &u[i * nc]);
    dq[i] = s->cbasis.eval_expand(z, &dtb[i * nc]);
  }
  for (int d = 0; d < cdim; ++d) {
    gxq[d] = s->cbasis.eval_expand(z, &gxi[(d * 3 + 0) * nc]);
    for (int j = 0; j < 3; ++j) {
      guq[d * 3 + j] = s->cbasis.eval_expand(z, &gu[(d * 3 + j) * nc]);
      gbq[d * 3 + j] = s->cbasis.eval_expand(z, &gb[(d * 3 + j) * nc]);
    }
  }
  pkpm_harm_coeffs_set(hc, cdim, bq, uq, s->cbasis.eval_expand(z, gkyl_array_cfetch(s->xi, cidx)),
    s->cbasis.eval_expand(z, gkyl_array_cfetch(s->theta, cidx)), dq, guq, gbq, gxq, 0, 0.0, 0);
}

// Node values of the w flux at reference point z of cell cidx (PKPM_HS_NODE_NV): the w-flux
// coefficients from the pointwise definitions of gkyl_pkpm_harm_coeffs.h (conf_coeffs; the build
// writes these raw values at its nodes), theta, b, the closures and the spectrum part of the bounds.
static void
node_values(
  struct test_setup *s, long cidx, const double *z, const struct gkyl_array *alpha_0,
  const struct gkyl_array *c_h, double *nv
)
{
  int nc = s->cbasis.num_basis;
  struct pkpm_harm_coeffs hc;
  conf_coeffs(s, cidx, z, &hc);
  double th = hc.theta, th2 = th * th, t[2][3], R[3][3];
  pkpm_harm_t(&hc, 0.0, t[0]);
  for (int a = 0; a < 3; ++a) {
    t[1][a] = hc.t1[a];
  }
  for (int c = 0; c < 3; ++c) {
    double e[3] = {0.0}, Re[3];
    e[c] = 1.0;
    pkpm_harm_rb1(&hc, e, Re);
    for (int a = 0; a < 3; ++a) {
      R[a][c] = Re[a];
    }
  }
  nv[PKPM_HS_KW] = -hc.K;
  nv[PKPM_HS_BG] = -hc.bg;
  for (int a = 0; a < 3; ++a) {
    for (int j = 0; j < 2; ++j) {
      nv[PKPM_HS_TF(a, j)] = t[j][a];
      nv[PKPM_HS_TG(a, j)] = t[j][a];
      nv[PKPM_HS_TM(a, j)] = t[j][a];
      nv[PKPM_HS_TNG(a, j)] = 4.0 * th * t[j][a];
      nv[PKPM_HS_TNF(a, j)] = -2.0 * th2 * t[j][a];
    }
    for (int c = 0; c < 3; ++c) {
      nv[PKPM_HS_RMN(a, c)] = R[a][c];
      nv[PKPM_HS_RNN(a, c)] = 6.0 * th * R[a][c];
      nv[PKPM_HS_RNM(a, c)] = -6.0 * th2 * R[a][c];
    }
  }
  nv[PKPM_HS_BF] = hc.beta;
  nv[PKPM_HS_BGG] = 4.0 * th * hc.beta;
  nv[PKPM_HS_BGF] = -2.0 * th2 * hc.beta;
  const double *b = gkyl_array_cfetch(s->bvar, cidx);
  nv[PKPM_HS_NODE_THETA] = th;
  for (int a = 0; a < 3; ++a) {
    nv[PKPM_HS_NODE_B + a] = s->cbasis.eval_expand(z, &b[a * nc]);
  }
  nv[PKPM_HS_NODE_A0] = s->cbasis.eval_expand(z, gkyl_array_cfetch(alpha_0, cidx));
  nv[PKPM_HS_NODE_CH] = s->cbasis.eval_expand(z, gkyl_array_cfetch(c_h, cidx));
  pkpm_harm_node_spec(nv);
}

// Independent reference of the stored v_par flux on the lower face of cell idx (idx_l is the cell
// below): the HLL flux of the pointwise node values (node_values, pkpm_harm_node_wflux) at the 5-point Gauss-Legendre nodes of every configuration direction,
// projected with that rule. Returns the CFL speed estimate at the nodes.
static double
reference_vel_face(
  struct test_setup *s, const struct gkyl_array *alpha_0, const struct gkyl_array *c_h,
  const int *idx_l, const int *idx, double *flux_ref, double *ratio
)
{
  int cdim = s->cdim;
  int nc = s->cbasis.num_basis, np = s->pbasis.num_basis;
  long cidx = gkyl_range_idx(&s->clocal_ext, idx);
  const double *f_l = gkyl_array_cfetch(s->f, gkyl_range_idx(&s->plocal_ext, idx_l));
  const double *f_r = gkyl_array_cfetch(s->f, gkyl_range_idx(&s->plocal_ext, idx));
  double xc[4];
  gkyl_rect_grid_cell_center(&s->pgrid, idx, xc);
  double v_face = xc[cdim] - 0.5 * s->pgrid.dx[cdim];

  for (int k = 0; k < 8 * nc; ++k) {
    flux_ref[k] = 0.0;
  }
  double cfl_max = 0.0;
  // the column bound of the w-flux speeds (absolute coefficients at the larger |w| edge, plan 17.27)
  double wm = fmax(fabs(s->pgrid.lower[cdim]), fabs(s->pgrid.upper[cdim]));
  *ratio = 0.0;
  int nnodes = 1;
  for (int d = 0; d < cdim; ++d) {
    nnodes *= 5;
  }
  for (int n = 0; n < nnodes; ++n) {
    double zs[4], zl[4], zr[4], w = 1.0;
    int m = n;
    for (int d = cdim - 1; d >= 0; --d) {
      zs[d] = gl5_x[m % 5];
      w *= gl5_w[m % 5];
      m /= 5;
    }
    for (int d = 0; d < cdim; ++d) {
      zl[d] = zr[d] = zs[d];
    }
    zl[cdim] = 1.0;
    zr[cdim] = -1.0;
    double nv[PKPM_HS_NODE_NV], Ul[8], Ur[8], ghat[8];
    node_values(s, cidx, zs, alpha_0, c_h, nv);
    for (int c = 0; c < 8; ++c) {
      Ul[c] = s->pbasis.eval_expand(zl, &f_l[c * np]);
      Ur[c] = s->pbasis.eval_expand(zr, &f_r[c * np]);
    }
    double speed = pkpm_harm_node_wflux(nv, v_face, Ul, Ur, ghat);
    cfl_max = fmax(cfl_max, speed);
    double bound = fabs(nv[PKPM_HS_NODE_A0]) + fabs(nv[PKPM_HS_NODE_CH] + nv[PKPM_HS_KW]) * wm +
      pkpm_harm_node_speed_rest(nv, wm);
    *ratio = fmax(*ratio, speed / bound);
    for (int c = 0; c < 8; ++c) {
      for (int k = 0; k < nc; ++k) {
        double unit[8] = {0.0};
        unit[k] = 1.0;
        double znorm[3] = {0.4, 0.4, 0.4};
        double chi_k = s->cbasis.eval_expand(zs, unit) / s->cbasis.eval_expand(znorm, unit);
        flux_ref[c * nc + k] += w * chi_k * ghat[c];
      }
    }
  }
  return cfl_max;
}

// The raw w-flux values the build writes (gkyl_dg_pkpm_harm_vol_coef) at the 5^cdim w-face nodes
// (node i0 + 5 i1 + 25 i2) are the pointwise definitions of gkyl_pkpm_harm_coeffs.h (node_values)
// at those points.
static void
test_build_nodes(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nn = 1;
  for (int d = 0; d < cdim; ++d) {
    nn *= 5;
  }
  double max_err = 0.0, max_ref = 0.0;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &s.clocal);
  while (gkyl_range_iter_next(&it)) {
    long cidx = gkyl_range_idx(&s.clocal_ext, it.idx);
    const double *wn = gkyl_array_cfetch(s.hwn, cidx);
    for (int n = 0; n < nn; ++n) {
      double z[3], nv[PKPM_HS_NODE_NV];
      for (int d = 0, q = n; d < cdim; ++d, q /= 5) {
        z[d] = gl5_x[q % 5];
      }
      node_values(&s, cidx, z, s.alpha_0, s.c_h, nv);
      for (int sl = 0; sl < PKPM_HS_W_END; ++sl) {
        max_err = fmax(max_err, fabs(wn[n * PKPM_HS_W_END + sl] - nv[sl]));
        max_ref = fmax(max_ref, fabs(nv[sl]));
      }
    }
  }
  TEST_CHECK(max_ref > 0.0);
  TEST_CHECK(max_err < 1e-13 * max_ref);
  TEST_MSG("cdim=%d: raw node values against the pointwise definitions, max error %g (max %g)", cdim,
    max_err, max_ref);
  teardown(&s);
}

static void
test_vel_flux(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nc = s.cbasis.num_basis;

  struct gkyl_dg_pkpm_vel_flux_surf_inp inp = {
    .phase_grid = &s.pgrid,
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .use_gpu = false,
  };
  struct gkyl_dg_pkpm_vel_flux_surf *up = gkyl_dg_pkpm_vel_flux_surf_inew(&inp);
  int nk = gkyl_dg_pkpm_vel_flux_surf_num_surf_basis(up);
  TEST_CHECK(nk == nc);
  struct gkyl_array *flux = mkarr(8 * nk, s.plocal_ext.volume);
  gkyl_array_clear(flux, 1.0e30); // every interior face must be overwritten; the edge slot is left alone
  gkyl_dg_pkpm_vel_flux_surf_advance(
    up, &s.clocal, &s.plocal, s.hwn, s.theta, s.bvar, s.alpha_0, s.c_h, s.f, flux
  );

  double flux_ref[64];
  double max_err = 0.0, max_flux = 0.0, max_ratio = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    long pidx = gkyl_range_idx(&s.plocal_ext, iter.idx);
    const double *fl = gkyl_array_cfetch(flux, pidx);
    // The lowest cell's slot is the velocity domain edge, written by the closure updater.
    if (iter.idx[cdim] == s.plocal.lower[cdim]) {
      for (int k = 0; k < 8 * nk; ++k) {
        TEST_CHECK(fl[k] == 1.0e30);
      }
      continue;
    }
    int idx_l[4];
    gkyl_copy_int_arr(s.pdim, iter.idx, idx_l);
    idx_l[cdim] -= 1;
    double ratio;
    reference_vel_face(&s, s.alpha_0, s.c_h, idx_l, iter.idx, flux_ref, &ratio);
    max_ratio = fmax(max_ratio, ratio);
    for (int k = 0; k < 8 * nk; ++k) {
      double err = fabs(fl[k] - flux_ref[k]);
      max_err = err > max_err ? err : max_err;
      max_flux = fabs(flux_ref[k]) > max_flux ? fabs(flux_ref[k]) : max_flux;
    }
  }
  TEST_CHECK(max_flux > 0.0);
  TEST_CHECK(max_err < 1e-12 * max_flux);
  TEST_MSG("cdim=%d: max flux error %g (max flux %g)", cdim, max_err, max_flux);
  // no interior face speed exceeds the column bound at the edges (plan 17.27)
  TEST_CHECK(max_ratio <= 1.0 + 1e-12);
  TEST_MSG("cdim=%d: largest face speed / column bound %g", cdim, max_ratio);

  gkyl_array_release(flux);
  gkyl_dg_pkpm_vel_flux_surf_release(up);
  teardown(&s);
}

// Exactness of the w faces: where the flux is upwind over the whole face (alpha_0 large, every
// speed positive), the stored flux equals the exact Galerkin projection of the POINTWISE flux
// alpha U_l + T(U_l) (pointwise coefficients of gkyl_pkpm_harm_coeffs.h, integrated with 6
// Gauss-Legendre points per configuration direction): the face rule integrates the raw node
// values exactly (degree 1 + 7 + 1 <= 9 per direction with 5 points).
static void
test_vel_flux_exact(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nc = s.cbasis.num_basis, np = s.pbasis.num_basis;
  double c0 = pow(sqrt(2.0), cdim);
  gkyl_array_clear(s.alpha_0, 0.0);
  gkyl_array_shiftc(s.alpha_0, 30.0 * c0, 0);
  gkyl_array_clear(s.c_h, 0.0);

  struct gkyl_dg_pkpm_vel_flux_surf_inp inp = {
    .phase_grid = &s.pgrid, .conf_basis = &s.cbasis, .phase_basis = &s.pbasis, .use_gpu = false};
  struct gkyl_dg_pkpm_vel_flux_surf *up = gkyl_dg_pkpm_vel_flux_surf_inew(&inp);
  struct gkyl_array *flux = mkarr(8 * nc, s.plocal_ext.volume), *cflrate = mkarr(1, s.plocal_ext.volume);
  gkyl_array_clear(flux, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_dg_pkpm_vel_flux_surf_advance(up, &s.clocal, &s.plocal, s.hwn, s.theta, s.bvar, s.alpha_0, s.c_h,
    s.f, flux);

  static const double gl6_x[6] = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
    0.2386191860831969, 0.6612093864662645, 0.9324695142031521};
  static const double gl6_w[6] = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910,
    0.4679139345726910, 0.3607615730481386, 0.1713244923791704};
  double max_err = 0.0, max_ref = 0.0, max_ref_h = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    if (iter.idx[cdim] == s.plocal.lower[cdim]) {
      continue;
    }
    int idx_l[4];
    gkyl_copy_int_arr(s.pdim, iter.idx, idx_l);
    idx_l[cdim] -= 1;
    long cidx = gkyl_range_idx(&s.clocal_ext, iter.idx);
    const double *f_l = gkyl_array_cfetch(s.f, gkyl_range_idx(&s.plocal_ext, idx_l));
    double xc[4];
    gkyl_rect_grid_cell_center(&s.pgrid, iter.idx, xc);
    double v_face = xc[cdim] - 0.5 * s.pgrid.dx[cdim], ref[8 * 8] = {0.0};
    int nnodes = 1;
    for (int d = 0; d < cdim; ++d) {
      nnodes *= 6;
    }
    for (int n = 0; n < nnodes; ++n) {
      double z[4], zl[4], w = 1.0;
      for (int d = cdim - 1, m = n; d >= 0; --d, m /= 6) {
        z[d] = zl[d] = gl6_x[m % 6];
        w *= gl6_w[m % 6];
      }
      zl[cdim] = 1.0;
      struct pkpm_harm_coeffs hc;
      conf_coeffs(&s, cidx, z, &hc);
      double alpha = s.cbasis.eval_expand(z, gkyl_array_cfetch(s.alpha_0, cidx)) - hc.K * v_face -
        hc.bg * v_face * v_face;
      TEST_CHECK(alpha > 0.0);
      double Ul[8], T[8];
      for (int c = 0; c < 8; ++c) {
        Ul[c] = s.pbasis.eval_expand(zl, &f_l[c * np]);
      }
      pkpm_harm_coupling_flux(&hc, v_face, Ul, T);
      for (int c = 0; c < 8; ++c) {
        for (int k = 0; k < nc; ++k) {
          double unit[8] = {0.0};
          unit[k] = 1.0;
          double znorm[3] = {0.4, 0.4, 0.4};
          double chi_k = s.cbasis.eval_expand(z, unit) / s.cbasis.eval_expand(znorm, unit);
          ref[c * nc + k] += w * chi_k * (alpha * Ul[c] + T[c]);
        }
      }
    }
    const double *fl = gkyl_array_cfetch(flux, gkyl_range_idx(&s.plocal_ext, iter.idx));
    for (int k = 0; k < 8 * nc; ++k) {
      max_err = fmax(max_err, fabs(fl[k] - ref[k]));
      max_ref = fmax(max_ref, fabs(ref[k]));
      if (k >= 2 * nc) {
        max_ref_h = fmax(max_ref_h, fabs(ref[k]));
      }
    }
  }
  TEST_CHECK(max_ref_h > 0.0);
  TEST_CHECK(max_err < 1e-12 * max_ref);
  TEST_MSG("cdim=%d: upwind w faces against the exact projection of the pointwise flux: max error %g (max %g)",
    cdim, max_err, max_ref);

  gkyl_array_release(cflrate);
  gkyl_array_release(flux);
  gkyl_dg_pkpm_vel_flux_surf_release(up);
  teardown(&s);
}

// The velocity-domain edges are constructed as the interior w faces: with previous closures that
// make the upper edge pure outflow and the lower edge inflow at every node (alpha_0 large), the
// stored upper edge flux is the exact Galerkin projection of the POINTWISE flux alpha U + T(U) at
// w_max with the NEW closures (6 Gauss-Legendre points per configuration direction), and the lower
// edge carries no flux.
static void
test_edge_flux_exact(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nc = s.cbasis.num_basis, np = s.pbasis.num_basis;
  double c0 = pow(sqrt(2.0), cdim);
  gkyl_array_clear(s.alpha_0, 0.0);
  gkyl_array_shiftc(s.alpha_0, 30.0 * c0, 0);
  gkyl_array_clear(s.c_h, 0.0);
  gkyl_array_clear(s.c_0, 0.0);

  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_pkpm_new(
    &s.pgrid, &s.cbasis, &s.pbasis, &s.clocal, &s.vlocal, 1.0, false, false
  );
  struct gkyl_array *moms = mkarr(18 * nc, s.clocal_ext.volume), *rhs_moms = mkarr(18 * nc, s.clocal_ext.volume);
  gkyl_array_clear(moms, 0.0);
  gkyl_array_clear(rhs_moms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(mcalc, &s.plocal, &s.clocal, s.f, moms);
  struct gkyl_dg_calc_pkpm_vars *vars = gkyl_dg_calc_pkpm_vars_new(&s.cgrid, &s.cbasis, s.pos_map, false);
  struct gkyl_array *vflux = mkarr(8 * nc, s.plocal_ext.volume), *cflrate = mkarr(1, s.plocal_ext.volume);
  gkyl_array_clear(vflux, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  int num_null = gkyl_dg_calc_pkpm_vars_closure(vars, &s.clocal, &s.plocal, &s.pgrid, moms, rhs_moms, s.bvar,
    s.hwn, s.theta, s.f, s.c_0, s.alpha_0, s.c_h, s.a_perp, vflux, cflrate);
  TEST_CHECK(num_null == 0);

  static const double gl6_x[6] = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
    0.2386191860831969, 0.6612093864662645, 0.9324695142031521};
  static const double gl6_w[6] = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910,
    0.4679139345726910, 0.3607615730481386, 0.1713244923791704};
  double w_up = s.pgrid.upper[cdim], max_err = 0.0, max_ref = 0.0, max_lo = 0.0;
  int nnodes = 1;
  for (int d = 0; d < cdim; ++d) {
    nnodes *= 6;
  }
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &s.clocal);
  while (gkyl_range_iter_next(&it)) {
    long cidx = gkyl_range_idx(&s.clocal_ext, it.idx);
    int idx[4];
    gkyl_copy_int_arr(cdim, it.idx, idx);
    idx[cdim] = s.plocal.upper[cdim];
    const double *f_up = gkyl_array_cfetch(s.f, gkyl_range_idx(&s.plocal_ext, idx));
    idx[cdim] = s.plocal.upper[cdim] + 1; // the velocity ghost cell's slot: the upper edge
    const double *fu = gkyl_array_cfetch(vflux, gkyl_range_idx(&s.plocal_ext, idx));
    idx[cdim] = s.plocal.lower[cdim]; // the lowest cell's slot: the lower edge
    const double *fl = gkyl_array_cfetch(vflux, gkyl_range_idx(&s.plocal_ext, idx));
    double ref[8 * 8] = {0.0};
    for (int n = 0; n < nnodes; ++n) {
      double z[4], zu[4], w = 1.0;
      for (int d = cdim - 1, m = n; d >= 0; --d, m /= 6) {
        z[d] = zu[d] = gl6_x[m % 6];
        w *= gl6_w[m % 6];
      }
      zu[cdim] = 1.0;
      struct pkpm_harm_coeffs hc;
      conf_coeffs(&s, cidx, z, &hc);
      double alpha = s.cbasis.eval_expand(z, gkyl_array_cfetch(s.alpha_0, cidx)) +
        s.cbasis.eval_expand(z, gkyl_array_cfetch(s.c_h, cidx)) * w_up - hc.K * w_up - hc.bg * w_up * w_up;
      double U[8], T[8];
      for (int c = 0; c < 8; ++c) {
        U[c] = s.pbasis.eval_expand(zu, &f_up[c * np]);
      }
      pkpm_harm_coupling_flux(&hc, w_up, U, T);
      for (int c = 0; c < 8; ++c) {
        for (int k = 0; k < nc; ++k) {
          double unit[8] = {0.0};
          unit[k] = 1.0;
          double znorm[3] = {0.4, 0.4, 0.4};
          double chi_k = s.cbasis.eval_expand(z, unit) / s.cbasis.eval_expand(znorm, unit);
          ref[c * nc + k] += w * chi_k * (alpha * U[c] + T[c]);
        }
      }
    }
    for (int k = 0; k < 8 * nc; ++k) {
      max_err = fmax(max_err, fabs(fu[k] - ref[k]));
      max_ref = fmax(max_ref, fabs(ref[k]));
      max_lo = fmax(max_lo, fabs(fl[k]));
    }
  }
  TEST_CHECK(max_ref > 0.0);
  TEST_CHECK(max_err < 1e-12 * max_ref);
  TEST_MSG("cdim=%d: upper edge against the exact projection of the pointwise flux: max error %g (max %g)",
    cdim, max_err, max_ref);
  TEST_CHECK(max_lo == 0.0);
  TEST_MSG("cdim=%d: inflow lower edge carries a flux (%g)", cdim, max_lo);

  gkyl_array_release(cflrate);
  gkyl_array_release(vflux);
  gkyl_dg_calc_pkpm_vars_release(vars);
  gkyl_array_release(rhs_moms);
  gkyl_array_release(moms);
  gkyl_dg_updater_moment_pkpm_release(mcalc);
  teardown(&s);
}

static void
test_vel_vol(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int np = s.pbasis.num_basis;

  struct gkyl_dg_pkpm_kinetic_vel_inp inp = {
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .conf_range = &s.clocal,
    .phase_range = &s.plocal,
    .auxfields = {.c_0 = s.c_0, .alpha_0 = s.alpha_0, .c_h = s.c_h, .a_perp = s.a_perp, .vel_flux_surf = NULL},
    .use_gpu = false,
  };
  struct gkyl_dg_eqn *eqn = gkyl_dg_pkpm_kinetic_vel_inew(&inp);

  double max_err = 0.0, max_ref = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    long pidx = gkyl_range_idx(&s.plocal_ext, iter.idx);
    long cidx = gkyl_range_idx(&s.clocal_ext, iter.idx);
    const double *f = gkyl_array_cfetch(s.f, pidx);
    const double *aP = gkyl_array_cfetch(s.alpha_0, cidx), *ch = gkyl_array_cfetch(s.c_h, cidx);
    const double *c0 = gkyl_array_cfetch(s.c_0, cidx), *ap = gkyl_array_cfetch(s.a_perp, cidx);
    int nc = s.cbasis.num_basis;
    double xc[4], out[8 * 24] = {0.0}, ref[8 * 24] = {0.0};
    gkyl_rect_grid_cell_center(&s.pgrid, iter.idx, xc);
    eqn->vol_term(eqn, xc, s.pgrid.dx, iter.idx, f, out);

    int nnodes = 1;
    for (int d = 0; d < s.pdim; ++d) {
      nnodes *= 4;
    }
    double dv1 = 2.0 / s.pgrid.dx[cdim];
    for (int n = 0; n < nnodes; ++n) {
      double z[4], w = 1.0;
      int m = n;
      for (int d = s.pdim - 1; d >= 0; --d) {
        z[d] = gl4_x[m % 4];
        w *= gl4_w[m % 4];
        m /= 4;
      }
      double aPq = s.cbasis.eval_expand(z, aP), chq = s.cbasis.eval_expand(z, ch);
      double c0q = s.cbasis.eval_expand(z, c0);
      double v = xc[cdim] + 0.5 * s.pgrid.dx[cdim] * z[cdim];
      double Uq[8], apq[3];
      for (int c = 0; c < 8; ++c) {
        Uq[c] = s.pbasis.eval_expand(z, &f[c * np]);
      }
      for (int a = 0; a < 3; ++a) {
        apq[a] = s.cbasis.eval_expand(z, &ap[a * nc]);
      }
      double Fq = Uq[0], Gq = Uq[1];
      // closure terms: force (alpha_0 + c w) U, sources c_0 U + k c U (k = 0, 2, 1, 3), and the
      // perpendicular force a_perp.M0 (G), a_perp F (M0), 2 G a_perp (N)
      double src[8], kc[8] = {0.0, 2.0, 1.0, 1.0, 1.0, 3.0, 3.0, 3.0};
      for (int c = 0; c < 8; ++c) {
        src[c] = (c0q + kc[c] * chq) * Uq[c];
      }
      src[1] += apq[0] * Uq[2] + apq[1] * Uq[3] + apq[2] * Uq[4];
      for (int a = 0; a < 3; ++a) {
        src[2 + a] += apq[a] * Fq;
        src[5 + a] += 2.0 * Gq * apq[a];
      }
      for (int k = 0; k < np; ++k) {
        double unit[64] = {0.0};
        unit[k] = 1.0;
        double phi = s.pbasis.eval_expand(z, unit);
        double dphi_v = dv1 * s.pbasis.eval_grad_expand(cdim, z, unit);
        for (int c = 0; c < 8; ++c) {
          ref[c * np + k] += w * (dphi_v * (aPq + chq * v) * Uq[c] + phi * src[c]);
        }
      }
    }
    for (int k = 0; k < 8 * np; ++k) {
      double err = fabs(out[k] - ref[k]);
      max_err = err > max_err ? err : max_err;
      max_ref = fabs(ref[k]) > max_ref ? fabs(ref[k]) : max_ref;
    }
  }
  TEST_CHECK(max_ref > 0.0);
  TEST_CHECK(max_err < 1e-12 * max_ref);
  TEST_MSG("cdim=%d: max volume error %g (max reference %g)", cdim, max_err, max_ref);

  gkyl_dg_eqn_release(eqn);
  teardown(&s);
}

// Largest perpendicular part of int M0_perp dw (moments mom, components 6-8) at the 2^cdim
// Gauss-Legendre nodes of cell cidx, where b is a unit vector.
static double
perp_moment_max(struct test_setup *s, long cidx, const struct gkyl_array *mom)
{
  int cdim = s->cdim, nc = s->cbasis.num_basis, nnodes = 1 << cdim;
  const double *m = gkyl_array_cfetch(mom, cidx), *b = gkyl_array_cfetch(s->bvar, cidx);
  double mx = 0.0;
  for (int n = 0; n < nnodes; ++n) {
    double z[3], bq[3], mq[3], mp[3];
    for (int d = 0; d < cdim; ++d) {
      z[d] = ((n >> (cdim - 1 - d)) & 1) ? 1.0 / sqrt(3.0) : -1.0 / sqrt(3.0);
    }
    for (int a = 0; a < 3; ++a) {
      bq[a] = s->cbasis.eval_expand(z, &b[a * nc]);
      mq[a] = s->cbasis.eval_expand(z, &m[(6 + a) * nc]);
    }
    pkpm_harm_perp(bq, mq, mp);
    mx = fmax(mx, sqrt(pkpm_harm_dot(mp, mp)));
  }
  return mx;
}

// The complete two-pass right-hand side: configuration fluxes, stage A, closures, w
// fluxes, stage B. The zeroth and first moments of the total and its scalar second
// moment M2(F) + 2 M0(G) must vanish in every cell and configuration mode, to round-off
// relative to the cancelling stage-A moments.
static void
test_m1_invariant(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nc = s.cbasis.num_basis, np = s.pbasis.num_basis;

  // Stage A.
  struct gkyl_dg_pkpm_conf_flux_surf_inp cinp = {
    .phase_grid = &s.pgrid, .conf_basis = &s.cbasis, .phase_basis = &s.pbasis, .pos_map = s.pos_map, .use_gpu = false};
  struct gkyl_dg_pkpm_conf_flux_surf *calc_cflux = gkyl_dg_pkpm_conf_flux_surf_inew(&cinp);
  int nkc = gkyl_dg_pkpm_conf_flux_surf_num_surf_basis(calc_cflux);
  struct gkyl_array *cflux = mkarr(2 * 8 * cdim * nkc, s.plocal_ext.volume);
  struct gkyl_array *cflrate = mkarr(1, s.plocal_ext.volume);
  struct gkyl_array *rhs = mkarr(8 * np, s.plocal_ext.volume);
  gkyl_array_clear(cflux, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);
  struct gkyl_array *face_coef = mkarr(gkyl_dg_pkpm_conf_flux_surf_face_coef_ncomp(calc_cflux), s.clocal_ext.volume);
  gkyl_dg_pkpm_conf_flux_surf_face_coef(calc_cflux, &s.clocal, s.bvar, s.u, s.xi, s.theta, face_coef);
  gkyl_dg_pkpm_conf_flux_surf_advance(
    calc_cflux, &s.clocal, &s.plocal, &s.plocal_ext, face_coef, s.f, cflrate, cflux
  );
  expand_two_sided(cflux, cdim, 8 * nkc);
  struct gkyl_dg_pkpm_kinetic_conf_inp ainp = {
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .conf_range = &s.clocal,
    .phase_range = &s.plocal,
    .auxfields =
      {
        .bvar = s.bvar,
        .u = s.u,
        .xi = s.xi,
        .conf_flux_surf = cflux,
      },
    .pos_map = s.pos_map, .use_gpu = false,
  };
  struct gkyl_dg_eqn *eqn_a = gkyl_dg_pkpm_kinetic_conf_inew(&ainp);
  int up_dirs_a[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_a[2 * GKYL_MAX_DIM] = {0};
  gkyl_hyper_dg *slvr_a =
    gkyl_hyper_dg_new(&s.pgrid, &s.pbasis, eqn_a, cdim, up_dirs_a, zero_flux_a, 1, false);
  gkyl_hyper_dg_advance(slvr_a, &s.plocal, s.f, cflrate, rhs);
  struct gkyl_dg_pkpm_harm_vol_inp hinp = {
    .phase_grid = &s.pgrid, .conf_basis = &s.cbasis, .phase_basis = &s.pbasis, .qbym = qbym_test, .pos_map = s.pos_map, .use_gpu = false};
  struct gkyl_dg_pkpm_harm_vol *harm = gkyl_dg_pkpm_harm_vol_inew(&hinp);
  {
    gkyl_dg_pkpm_harm_vol_advance(harm, &s.clocal, &s.plocal, s.hcoef, s.f, cflrate, rhs);
  }
  gkyl_dg_pkpm_harm_vol_release(harm);

  // Closure.
  double mass = 1.0; // the normalized distributions carry mass 1 (the edge moments assume it)
  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_pkpm_new(
    &s.pgrid, &s.cbasis, &s.pbasis, &s.clocal, &s.vlocal, mass, false, false
  );
  struct gkyl_array *moms = mkarr(18 * nc, s.clocal_ext.volume);
  struct gkyl_array *rhs_moms = mkarr(18 * nc, s.clocal_ext.volume);
  gkyl_array_clear(moms, 0.0);
  gkyl_array_clear(rhs_moms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(mcalc, &s.plocal, &s.clocal, s.f, moms);
  gkyl_dg_updater_moment_pkpm_advance(mcalc, &s.plocal, &s.clocal, rhs, rhs_moms);
  struct gkyl_dg_pkpm_vel_flux_surf_inp vinp = {
    .phase_grid = &s.pgrid, .conf_basis = &s.cbasis, .phase_basis = &s.pbasis, .use_gpu = false};
  struct gkyl_dg_pkpm_vel_flux_surf *calc_vflux = gkyl_dg_pkpm_vel_flux_surf_inew(&vinp);
  struct gkyl_array *vflux = mkarr(8 * nc, s.plocal_ext.volume);
  gkyl_array_clear(vflux, 0.0);

  struct gkyl_dg_calc_pkpm_vars *vars =
    gkyl_dg_calc_pkpm_vars_new(&s.cgrid, &s.cbasis, s.pos_map, false);
  struct gkyl_array *a_P = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *c_h = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *c_0 = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *a_perp = mkarr(3 * nc, s.clocal_ext.volume);
  // closures of the previous stage (the upwind switch of the edges): the setup fields
  gkyl_array_copy(a_P, s.alpha_0);
  gkyl_array_copy(c_h, s.c_h);
  gkyl_array_clear(c_0, 0.0);
  int num_null = gkyl_dg_calc_pkpm_vars_closure(
    vars, &s.clocal, &s.plocal, &s.pgrid, moms, rhs_moms, s.bvar, s.hwn, s.theta, s.f, c_0, a_P, c_h, a_perp,
    vflux, cflrate
  );
  TEST_CHECK(num_null == 0);

  double max_m1_a = 0.0, max_m0_a = 0.0, max_m2_a = 0.0, max_mp_a = 0.0, max_mu = 0.0;
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, &s.clocal);
  while (gkyl_range_iter_next(&citer)) {
    long cidx = gkyl_range_idx(&s.clocal_ext, citer.idx);
    const double *rm = gkyl_array_cfetch(rhs_moms, cidx);
    for (int k = 0; k < nc; ++k) {
      max_m1_a = fabs(rm[3 * nc + k]) > max_m1_a ? fabs(rm[3 * nc + k]) : max_m1_a;
      max_m0_a = fabs(rm[k]) > max_m0_a ? fabs(rm[k]) : max_m0_a;
      double m2 = fabs(rm[nc + k] + 2.0 * rm[2 * nc + k]);
      max_m2_a = m2 > max_m2_a ? m2 : max_m2_a;
    }
    max_mp_a = fmax(max_mp_a, perp_moment_max(&s, cidx, rhs_moms));
    max_mu = fmax(max_mu, perp_moment_max(&s, cidx, moms));
  }

  // Stage B (interior faces with the closures; the edge faces were set above).
  gkyl_dg_pkpm_vel_flux_surf_advance(
    calc_vflux, &s.clocal, &s.plocal, s.hwn, s.theta, s.bvar, a_P, c_h, s.f, vflux
  );
  struct gkyl_dg_pkpm_kinetic_vel_inp binp = {
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .conf_range = &s.clocal,
    .phase_range = &s.plocal,
    .auxfields = {.c_0 = c_0, .alpha_0 = a_P, .c_h = c_h, .a_perp = a_perp, .vel_flux_surf = vflux},
    .use_gpu = false,
  };
  struct gkyl_dg_eqn *eqn_b = gkyl_dg_pkpm_kinetic_vel_inew(&binp);
  int up_dirs_b[GKYL_MAX_DIM] = {cdim}, zero_flux_b[2 * GKYL_MAX_DIM] = {0};
  zero_flux_b[cdim] = 1; // no ghost cell below w = -w_max; the upper edge is a regular face against the ghost cell
  gkyl_hyper_dg *slvr_b =
    gkyl_hyper_dg_new(&s.pgrid, &s.pbasis, eqn_b, 1, up_dirs_b, zero_flux_b, 1, false);
  gkyl_hyper_dg_advance(slvr_b, &s.plocal, s.f, cflrate, rhs);

  // Zeroth and first moments of the total right-hand side.
  gkyl_array_clear(rhs_moms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(mcalc, &s.plocal, &s.clocal, rhs, rhs_moms);
  double max_m1 = 0.0, max_m0 = 0.0, max_m2 = 0.0, max_mp = 0.0;
  gkyl_range_iter_init(&citer, &s.clocal);
  while (gkyl_range_iter_next(&citer)) {
    long cidx = gkyl_range_idx(&s.clocal_ext, citer.idx);
    const double *rm = gkyl_array_cfetch(rhs_moms, cidx);
    for (int k = 0; k < nc; ++k) {
      max_m1 = fabs(rm[3 * nc + k]) > max_m1 ? fabs(rm[3 * nc + k]) : max_m1;
      max_m0 = fabs(rm[k]) > max_m0 ? fabs(rm[k]) : max_m0;
      double m2 = fabs(rm[nc + k] + 2.0 * rm[2 * nc + k]);
      max_m2 = m2 > max_m2 ? m2 : max_m2;
    }
    max_mp = fmax(max_mp, perp_moment_max(&s, cidx, rhs_moms));
  }
  TEST_CHECK(max_m1_a > 0.0);
  TEST_CHECK(max_m1 < 1e-12 * max_m1_a);
  TEST_MSG("cdim=%d: max |M1(rhs)| %g after closure (stage-A moment %g)", cdim, max_m1, max_m1_a);
  TEST_CHECK(max_m0_a > 0.0);
  TEST_CHECK(max_m0 < 1e-12 * max_m0_a);
  TEST_MSG("cdim=%d: max |M0(rhs)| %g after closure (stage-A moment %g)", cdim, max_m0, max_m0_a);
  TEST_CHECK(max_m2_a > 0.0);
  TEST_CHECK(max_m2 < 1e-12 * max_m2_a);
  TEST_MSG("cdim=%d: max |M2(rhs_F) + 2 M0(rhs_G)| %g after closure (stage-A moment %g)", cdim, max_m2, max_m2_a);
  TEST_CHECK(max_mu > 0.0); // off the constraint manifold: the coupled 5x5 is exercised
  TEST_CHECK(max_mp_a > 0.0);
  TEST_CHECK(max_mp < 1e-12 * max_mp_a);
  TEST_MSG("cdim=%d: max |I_perp int rhs(M0_perp) dw| %g at the nodes after closure (stage-A moment %g, |mu| %g)",
    cdim, max_mp, max_mp_a, max_mu);
  gkyl_array_release(a_perp);

  gkyl_hyper_dg_release(slvr_b);
  gkyl_dg_eqn_release(eqn_b);
  gkyl_array_release(vflux);
  gkyl_dg_pkpm_vel_flux_surf_release(calc_vflux);
  gkyl_array_release(a_P);
  gkyl_array_release(c_h);
  gkyl_array_release(c_0);
  gkyl_dg_calc_pkpm_vars_release(vars);
  gkyl_array_release(rhs_moms);
  gkyl_array_release(moms);
  gkyl_dg_updater_moment_pkpm_release(mcalc);
  gkyl_hyper_dg_release(slvr_a);
  gkyl_dg_eqn_release(eqn_a);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_array_release(cflux);
  gkyl_array_release(face_coef);
  gkyl_dg_pkpm_conf_flux_surf_release(calc_cflux);
  teardown(&s);
}

void test_build_nodes_1x() { test_build_nodes(1); }
void test_build_nodes_2x() { test_build_nodes(2); }
void test_build_nodes_3x() { test_build_nodes(3); }
void test_vel_flux_1x() { test_vel_flux(1); }
void test_vel_flux_exact_1x() { test_vel_flux_exact(1); }
void test_edge_flux_exact_1x() { test_edge_flux_exact(1); }
void test_edge_flux_exact_2x() { test_edge_flux_exact(2); }
void test_edge_flux_exact_3x() { test_edge_flux_exact(3); }
void test_vel_flux_exact_2x() { test_vel_flux_exact(2); }
void test_vel_flux_exact_3x() { test_vel_flux_exact(3); }
void test_vel_flux_2x() { test_vel_flux(2); }
void test_vel_flux_3x() { test_vel_flux(3); }
void test_vel_vol_1x() { test_vel_vol(1); }
void test_vel_vol_2x() { test_vel_vol(2); }
void test_vel_vol_3x() { test_vel_vol(3); }
void test_m1_invariant_1x() { test_m1_invariant(1); }
void test_m1_invariant_2x() { test_m1_invariant(2); }
void test_m1_invariant_3x() { test_m1_invariant(3); }

// Characteristic condition at the velocity domain edges inside the closure: a uniform state
// whose distribution is concentrated near the upper edge w = 4, a prescribed stage-A
// energy residual (heating, so the new alpha(w_max) is positive), and the previous
// closures pointing outward (alpha_0 = 0.7) or inward (-0.7) at the upper edge. In both
// cases the 3x3 solve balances the stage-A residual with c_0, alpha_0, c such that the
// complete stage-B right-hand side (closure terms + edge flux) cancels it in the
// constrained moments. Outward: the stored edge flux equals alpha(w_max) F_edge of the
// NEW closures. Inward: nothing leaves (the lower edge holds no distribution).
static void
eval_edge_f(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double w = xn[cdim], F = exp(-0.5 * (w - 3.0) * (w - 3.0) / 0.25);
  fout[0] = F;
  fout[1] = F;
}

static void
test_edge_closure(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nc = s.cbasis.num_basis, np = s.pbasis.num_basis;
  double c0 = pow(sqrt(2.0), cdim);
  gkyl_array_clear(s.grad_b, 0.0);
  gkyl_array_clear(s.grad_u, 0.0);
  gkyl_array_clear(s.grad_xi, 0.0);
  gkyl_array_clear(s.u, 0.0);
  gkyl_array_clear(s.bvar, 0.0);
  gkyl_array_shiftc(s.bvar, c0, 0);
  gkyl_array_clear(s.xi, 0.0);
  gkyl_array_shiftc(s.xi, c0, 0);
  gkyl_array_shiftc(s.xi, c0, nc);
  gkyl_array_clear(s.theta, 0.0);
  gkyl_array_shiftc(s.theta, c0, 0);
  gkyl_array_clear(s.dtb, 0.0); // a static field: no first-harmonic coupling, plain upwind edges
  gkyl_array_clear(s.grad_ln_n, 0.0);
  gkyl_array_clear(s.f, 0.0);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s.pgrid, &s.pbasis, 8, 2, eval_edge_f, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s.plocal_ext_conf, s.f);
  gkyl_proj_on_basis_release(proj);
  build_coef(&s);
  int idx_top[GKYL_MAX_DIM];
  gkyl_copy_int_arr(s.pdim, s.plocal.lower, idx_top);
  idx_top[cdim] = s.plocal.upper[cdim];
  double z[4] = {0.0, 0.0, 0.0, 0.0};
  z[cdim] = 1.0;
  double F_edge = s.pbasis.eval_expand(z, gkyl_array_cfetch(s.f, gkyl_range_idx(&s.plocal_ext, idx_top)));
  TEST_CHECK(F_edge > 0.05);
  double w_up = s.pgrid.upper[cdim];

  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_pkpm_new(
    &s.pgrid, &s.cbasis, &s.pbasis, &s.clocal, &s.vlocal, 1.0, false, false
  );
  struct gkyl_array *moms = mkarr(18 * nc, s.clocal_ext.volume), *rhs_moms = mkarr(18 * nc, s.clocal_ext.volume);
  struct gkyl_array *b_moms = mkarr(18 * nc, s.clocal_ext.volume);
  gkyl_array_clear(moms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(mcalc, &s.plocal, &s.clocal, s.f, moms);
  struct gkyl_dg_calc_pkpm_vars *vars = gkyl_dg_calc_pkpm_vars_new(&s.cgrid, &s.cbasis, s.pos_map, false);
  struct gkyl_array *rhs = mkarr(8 * np, s.plocal_ext.volume), *cflrate = mkarr(1, s.plocal_ext.volume);
  struct gkyl_dg_pkpm_vel_flux_surf_inp vinp = {
    .phase_grid = &s.pgrid, .conf_basis = &s.cbasis, .phase_basis = &s.pbasis, .use_gpu = false};
  struct gkyl_dg_pkpm_vel_flux_surf *calc_vflux = gkyl_dg_pkpm_vel_flux_surf_inew(&vinp);
  struct gkyl_array *vflux = mkarr(8 * nc, s.plocal_ext.volume);
  struct gkyl_dg_pkpm_kinetic_vel_inp binp = {
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .conf_range = &s.clocal,
    .phase_range = &s.plocal,
    .auxfields = {.c_0 = s.c_0, .alpha_0 = s.alpha_0, .c_h = s.c_h, .a_perp = s.a_perp, .vel_flux_surf = vflux},
    .use_gpu = false,
  };
  struct gkyl_dg_eqn *eqn_b = gkyl_dg_pkpm_kinetic_vel_inew(&binp);
  int up_dirs_b[GKYL_MAX_DIM] = {cdim}, zero_flux_b[2 * GKYL_MAX_DIM] = {0};
  zero_flux_b[cdim] = 1; // no ghost cell below w = -w_max; the upper edge is a regular face against the ghost cell
  gkyl_hyper_dg *slvr_b = gkyl_hyper_dg_new(&s.pgrid, &s.pbasis, eqn_b, 1, up_dirs_b, zero_flux_b, 1, false);

  double alphas[2] = {0.7, -0.7};
  for (int a = 0; a < 2; ++a) {
    // previous-stage closures (the switch), no configuration pass
    gkyl_array_clear(s.alpha_0, 0.0);
    gkyl_array_shiftc(s.alpha_0, alphas[a] * c0, 0);
    gkyl_array_clear(s.c_h, 0.0);
    gkyl_array_clear(s.c_0, 0.0);
    gkyl_array_clear(rhs_moms, 0.0);
    gkyl_array_shiftc(rhs_moms, -4.0 * c0, nc); // stage-A rate of M2(F): the energy residual
    gkyl_array_clear(vflux, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    int num_null = gkyl_dg_calc_pkpm_vars_closure(
      vars, &s.clocal, &s.plocal, &s.pgrid, moms, rhs_moms, s.bvar, s.hwn, s.theta, s.f, s.c_0, s.alpha_0,
      s.c_h, s.a_perp, vflux, cflrate
    );
    TEST_CHECK(num_null == 0);
    gkyl_dg_pkpm_vel_flux_surf_advance(
      calc_vflux, &s.clocal, &s.plocal, s.hwn, s.theta, s.bvar, s.alpha_0, s.c_h, s.f, vflux
    );
    gkyl_array_clear(rhs, 0.0);
    gkyl_hyper_dg_advance(slvr_b, &s.plocal, s.f, cflrate, rhs);
    gkyl_array_clear(b_moms, 0.0);
    gkyl_dg_updater_moment_pkpm_advance(mcalc, &s.plocal, &s.clocal, rhs, b_moms);
    double max_cons = 0.0, max_a0 = 0.0, max_flux_err = 0.0, max_flux = 0.0, max_stored = 0.0;
    struct gkyl_range_iter citer;
    gkyl_range_iter_init(&citer, &s.clocal);
    while (gkyl_range_iter_next(&citer)) {
      long loc = gkyl_range_idx(&s.clocal_ext, citer.idx);
      const double *ra = gkyl_array_cfetch(rhs_moms, loc), *rb = gkyl_array_cfetch(b_moms, loc);
      const double *a0 = gkyl_array_cfetch(s.alpha_0, loc), *ch = gkyl_array_cfetch(s.c_h, loc);
      int idx_gh[GKYL_MAX_DIM];
      gkyl_copy_int_arr(cdim, citer.idx, idx_gh);
      idx_gh[cdim] = s.plocal.upper[cdim] + 1; // the velocity ghost cell's slot: the upper edge
      const double *fu = gkyl_array_cfetch(vflux, gkyl_range_idx(&s.plocal_ext, idx_gh));
      for (int k = 0; k < nc; ++k) {
        double m0 = ra[k] + rb[k], m1 = ra[3 * nc + k] + rb[3 * nc + k];
        double m2 = ra[nc + k] + 2.0 * ra[2 * nc + k] + rb[nc + k] + 2.0 * rb[2 * nc + k];
        max_cons = fmax(max_cons, fmax(fabs(m0), fmax(fabs(m1), fabs(m2))));
      }
      max_a0 = fmax(max_a0, fabs(a0[0]) / c0);
      // stored upper flux (uniform state: mode 0 only) vs max(alpha(w_max), 0) F_edge of the new closures,
      // in the stored format (orthonormal coefficient times sqrt(2)^cdim)
      double alpha_new = (a0[0] + ch[0] * w_up) / c0;
      double expect = fmax(alpha_new, 0.0) * F_edge * c0 * c0;
      max_flux_err = fmax(max_flux_err, fabs(fu[0] - expect));
      max_flux = fmax(max_flux, fabs(expect));
      for (int k = 0; k < 8 * nc; ++k) max_stored = fmax(max_stored, fabs(fu[k]));
    }
    TEST_CHECK(max_cons < 1e-11);
    TEST_MSG("cdim=%d alpha_prev=%g: constrained moments of stage A + stage B %g", cdim, alphas[a], max_cons);
    TEST_CHECK(max_a0 > 0.1);
    TEST_MSG("cdim=%d alpha_prev=%g: the forcing produced no closure (|alpha_0| %g)", cdim, alphas[a], max_a0);
    if (alphas[a] > 0.0) {
      TEST_CHECK(max_flux > 0.0 && max_flux_err < 1e-11 * max_flux);
      TEST_MSG("cdim=%d: stored edge flux vs closure-consistent outflow, error %g (flux %g)", cdim, max_flux_err, max_flux);
    } else {
      TEST_CHECK(max_stored < 1e-13);
      TEST_MSG("cdim=%d: inward case produced a flux (%g)", cdim, max_stored);
    }
  }

  gkyl_hyper_dg_release(slvr_b);
  gkyl_dg_eqn_release(eqn_b);
  gkyl_array_release(vflux);
  gkyl_dg_pkpm_vel_flux_surf_release(calc_vflux);
  gkyl_array_release(cflrate);
  gkyl_array_release(rhs);
  gkyl_dg_calc_pkpm_vars_release(vars);
  gkyl_array_release(b_moms);
  gkyl_array_release(rhs_moms);
  gkyl_array_release(moms);
  gkyl_dg_updater_moment_pkpm_release(mcalc);
  teardown(&s);
}
void test_edge_closure_1x() { test_edge_closure(1); }
void test_edge_closure_2x() { test_edge_closure(2); }

// Transversality projection: at the Gauss-Legendre nodes of b, b.M0_perp = b.N_perp = 0 for every
// w after the projection (it was not before), F and G are untouched, the perpendicular part is
// unchanged, and a second projection changes nothing.
static void
test_transverse(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nc = s.cbasis.num_basis, np = s.pbasis.num_basis, nnodes = 1 << cdim;
  // a fresh, non-transverse state
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s.pgrid, &s.pbasis, 8, 8, eval_f, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s.plocal_ext_conf, s.f);
  gkyl_proj_on_basis_release(proj);
  struct gkyl_array *f0 = mkarr(8 * np, s.plocal_ext.volume), *f1 = mkarr(8 * np, s.plocal_ext.volume);
  gkyl_array_copy(f0, s.f);
  struct gkyl_dg_calc_pkpm_vars *vars =
    gkyl_dg_calc_pkpm_vars_new(&s.cgrid, &s.cbasis, s.pos_map, false);
  gkyl_dg_calc_pkpm_vars_transverse(vars, &s.clocal_ext, &s.plocal_ext_conf, s.bvar, s.f);
  gkyl_array_copy(f1, s.f);
  gkyl_dg_calc_pkpm_vars_transverse(vars, &s.clocal_ext, &s.plocal_ext_conf, s.bvar, f1);

  double max_par0 = 0.0, max_par = 0.0, max_perr = 0.0, max_fg = 0.0, max_idem = 0.0, max_val = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal_ext_conf);
  while (gkyl_range_iter_next(&iter)) {
    long pidx = gkyl_range_idx(&s.plocal_ext, iter.idx), cidx = gkyl_range_idx(&s.clocal_ext, iter.idx);
    const double *fa = gkyl_array_cfetch(f0, pidx), *fb = gkyl_array_cfetch(s.f, pidx);
    const double *fc = gkyl_array_cfetch(f1, pidx), *b = gkyl_array_cfetch(s.bvar, cidx);
    for (int k = 0; k < 2 * np; ++k) {
      max_fg = fmax(max_fg, fabs(fb[k] - fa[k]));
    }
    for (int k = 0; k < 8 * np; ++k) {
      max_idem = fmax(max_idem, fabs(fc[k] - fb[k]));
      max_val = fmax(max_val, fabs(fb[k]));
    }
    for (int n = 0; n < nnodes; ++n) {
      for (int j = 0; j < 5; ++j) {
        double z[4], bq[3];
        for (int d = 0; d < cdim; ++d) {
          z[d] = ((n >> (cdim - 1 - d)) & 1) ? 1.0 / sqrt(3.0) : -1.0 / sqrt(3.0);
        }
        z[cdim] = -1.0 + 0.5 * j;
        for (int a = 0; a < 3; ++a) {
          bq[a] = s.cbasis.eval_expand(z, &b[a * nc]);
        }
        for (int vec = 0; vec < 2; ++vec) {
          double va[3], vb[3], pa[3];
          for (int a = 0; a < 3; ++a) {
            va[a] = s.pbasis.eval_expand(z, &fa[(2 + 3 * vec + a) * np]);
            vb[a] = s.pbasis.eval_expand(z, &fb[(2 + 3 * vec + a) * np]);
          }
          pkpm_harm_perp(bq, va, pa);
          max_par0 = fmax(max_par0, fabs(pkpm_harm_dot(bq, va)));
          max_par = fmax(max_par, fabs(pkpm_harm_dot(bq, vb)));
          for (int a = 0; a < 3; ++a) {
            max_perr = fmax(max_perr, fabs(vb[a] - pa[a]));
          }
        }
      }
    }
  }
  TEST_CHECK(max_par0 > 1e-3);
  TEST_CHECK(max_par < 1e-14 * max_val);
  TEST_MSG("cdim=%d: b.M at the nodes %g after the projection (%g before)", cdim, max_par, max_par0);
  TEST_CHECK(max_perr < 1e-13 * max_val);
  TEST_MSG("cdim=%d: perpendicular part changed by %g", cdim, max_perr);
  TEST_CHECK(max_fg == 0.0);
  TEST_CHECK(max_idem < 1e-13 * max_val);
  TEST_MSG("cdim=%d: second projection changed the state by %g (max %g)", cdim, max_idem, max_val);

  gkyl_dg_calc_pkpm_vars_release(vars);
  gkyl_array_release(f0);
  gkyl_array_release(f1);
  teardown(&s);
}
void test_transverse_1x() { test_transverse(1); }
void test_transverse_2x() { test_transverse(2); }
void test_transverse_3x() { test_transverse(3); }


TEST_LIST = {
  {"test_edge_closure_1x", test_edge_closure_1x},
  {"test_edge_closure_2x", test_edge_closure_2x},
  {"test_build_nodes_1x", test_build_nodes_1x},
  {"test_build_nodes_2x", test_build_nodes_2x},
  {"test_build_nodes_3x", test_build_nodes_3x},
  {"test_vel_flux_1x", test_vel_flux_1x},
  {"test_vel_flux_2x", test_vel_flux_2x},
  {"test_vel_flux_3x", test_vel_flux_3x},
  {"test_vel_flux_exact_1x", test_vel_flux_exact_1x},
  {"test_vel_flux_exact_2x", test_vel_flux_exact_2x},
  {"test_vel_flux_exact_3x", test_vel_flux_exact_3x},
  {"test_edge_flux_exact_1x", test_edge_flux_exact_1x},
  {"test_edge_flux_exact_2x", test_edge_flux_exact_2x},
  {"test_edge_flux_exact_3x", test_edge_flux_exact_3x},
  {"test_vel_vol_1x", test_vel_vol_1x},
  {"test_vel_vol_2x", test_vel_vol_2x},
  {"test_vel_vol_3x", test_vel_vol_3x},
  {"test_m1_invariant_1x", test_m1_invariant_1x},
  {"test_m1_invariant_2x", test_m1_invariant_2x},
  {"test_m1_invariant_3x", test_m1_invariant_3x},
  {"test_transverse_1x", test_transverse_1x},
  {"test_transverse_2x", test_transverse_2x},
  {"test_transverse_3x", test_transverse_3x},
  {NULL, NULL},
};
