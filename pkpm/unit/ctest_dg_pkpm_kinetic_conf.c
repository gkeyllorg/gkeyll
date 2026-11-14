// Tests of the configuration-space pass (stage A) of the refactored PKPM kinetic
// update of the eight blocks [F, G, M0_perp, N_perp]: the exact modal volume kernel of the
// two-field terms against an independent Gauss-Legendre quadrature of the same terms, the modal
// first-harmonic volume integrals (dg_pkpm_harm_vol) against a quadrature exact for them of the
// pointwise coefficients of gkyl_pkpm_harm_coeffs.h, and the closures (weak residuals of
// M0(f) c_0 + M0(rhs_A) = 0 and M0(f) a_P + M1(rhs_A) + c_0 M1(f) = 0 in every mode, and the
// perpendicular rows M0 a_perp + I_perp int rhs_A(M0_perp) dw = 0 at the closure nodes).
#include <math.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_pkpm_remap_flux_surf.h>
#include <gkyl_dg_pkpm_harm_vol.h>
#include <gkyl_dg_pkpm_harm_vol_priv.h>
#include <gkyl_dg_pkpm_kinetic_conf.h>
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
  fout[0] = F0 * (1.0 + 0.1 * v);
  fout[1] = 0.5 * vt * vt * F0 * (1.0 - 0.1 * v);
  // generic (not transverse) first-harmonic blocks; the setup projects them as the app does
  fout[2] = F0 * (0.3 * sin(x) + 0.1 * v);
  fout[3] = F0 * (-0.2 * cos(y) + 0.15 * v * v - 0.1);
  fout[4] = F0 * (0.25 * cos(z) - 0.05 * v);
  fout[5] = 0.5 * vt * vt * F0 * (0.2 * cos(x) - 0.1 * v);
  fout[6] = 0.5 * vt * vt * F0 * (0.1 + 0.2 * sin(y) * v);
  fout[7] = 0.5 * vt * vt * F0 * (-0.3 * sin(x + z) + 0.05 * v * v);
}

// Fluid state [rho, rho u, E]: the face flux reads the mass density only.
static void
eval_fluid(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  fout[0] = 2.0 * (1.0 + 0.3 * sin(x) * cos(y) + 0.1 * cos(z));
  for (int i = 1; i < 5; ++i) {
    fout[i] = 0.0;
  }
}

// d_t B/|B| of a time-dependent field, and grad(ln n).
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

// 6-point Gauss-Legendre rule: exact through degree 11 (the first-harmonic integrands reach
// degree 9 per configuration direction and 6 in w).
static const double gl6_x[6] = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
  0.2386191860831969, 0.6612093864662645, 0.9324695142031521};
static const double gl6_w[6] = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910,
  0.4679139345726910, 0.3607615730481386, 0.1713244923791704};
static const double qbym_test = 0.7; // charge to mass ratio of the gyration term

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
  struct gkyl_array *bvar, *grad_b, *u, *grad_u, *grad_ln_n, *theta, *f;
  struct gkyl_array *xi, *grad_xi, *em, *dtb;
  struct gkyl_dg_calc_pkpm_bvar *calc_b;
};

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

// Weak gradient of u on the interior cells (u is projected on the extended range).
static void
compute_grad_u(struct test_setup *s)
{
  int cdim = s->cdim, nc = s->cbasis.num_basis;
  s->grad_u = mkarr(cdim * 3 * nc, s->clocal_ext.volume);
  s->grad_xi = mkarr(cdim * 3 * nc, s->clocal_ext.volume);
  s->grad_ln_n = mkarr(cdim * nc, s->clocal_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, cdim, eval_gln, &s->cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->grad_ln_n);
  gkyl_proj_on_basis_release(proj);
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
  s->dtb = mkarr(3 * nc, s->clocal_ext.volume);
  struct gkyl_array *scal = mkarr(nc, s->clocal_ext.volume);
  s->bvar = mkarr(3 * nc, s->clocal_ext.volume);
  s->u = mkarr(3 * nc, s->clocal_ext.volume);
  s->xi = mkarr(3 * nc, s->clocal_ext.volume);
  s->theta = mkarr(nc, s->clocal_ext.volume);
  s->f = mkarr(8 * np, s->plocal_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 8, eval_em, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, em);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 3, eval_dtb, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->dtb);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 3, eval_u, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->u);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 1, eval_conf_scalars, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, scal);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_set_offset(s->theta, 1.0, scal, 0 * nc);
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
}

static void
teardown(struct test_setup *s)
{
  gkyl_vlasov_position_map_release(s->pos_map);
  gkyl_dg_calc_pkpm_bvar_release(s->calc_b);
  gkyl_array_release(s->grad_b);
  gkyl_array_release(s->grad_u);
  gkyl_array_release(s->f);
  gkyl_array_release(s->theta);
  gkyl_array_release(s->xi);
  gkyl_array_release(s->grad_xi);
  gkyl_array_release(s->grad_ln_n);
  gkyl_array_release(s->em);
  gkyl_array_release(s->dtb);
  gkyl_array_release(s->u);
  gkyl_array_release(s->bvar);
}

static struct gkyl_dg_eqn *
make_eqn(struct test_setup *s, const struct gkyl_array *flux)
{
  struct gkyl_dg_pkpm_kinetic_conf_inp inp = {
    .conf_basis = &s->cbasis,
    .phase_basis = &s->pbasis,
    .conf_range = &s->clocal,
    .phase_range = &s->plocal,
    .auxfields =
      {
        .bvar = s->bvar,
        .u = s->u,
        .xi = s->xi,
        .conf_flux_surf = flux,
      },
    .pos_map = s->pos_map, .use_gpu = false,
  };
  return gkyl_dg_pkpm_kinetic_conf_inew(&inp);
}

// Independent evaluation of the first-harmonic volume terms of one phase cell: the 6-point rule
// in every phase direction (exact for them) applied to the pointwise coefficients of
// gkyl_pkpm_harm_coeffs.h (whose definitions ctest_pkpm_harm_coeffs checks) from the fields
// evaluated at the nodes.
static void
reference_harm_vol(struct test_setup *s, const int *idx, double *out_ref)
{
  int cdim = s->cdim, pdim = s->pdim;
  int nc = s->cbasis.num_basis, np = s->pbasis.num_basis;
  long cidx = gkyl_range_idx(&s->clocal_ext, idx), pidx = gkyl_range_idx(&s->plocal_ext, idx);
  const double *b = gkyl_array_cfetch(s->bvar, cidx), *u = gkyl_array_cfetch(s->u, cidx);
  const double *gu = gkyl_array_cfetch(s->grad_u, cidx), *gb = gkyl_array_cfetch(s->grad_b, cidx);
  const double *gln = gkyl_array_cfetch(s->grad_ln_n, cidx), *gxi = gkyl_array_cfetch(s->grad_xi, cidx);
  const double *xi = gkyl_array_cfetch(s->xi, cidx), *th = gkyl_array_cfetch(s->theta, cidx);
  const double *dtb = gkyl_array_cfetch(s->dtb, cidx), *em = gkyl_array_cfetch(s->em, cidx);
  const double *f = gkyl_array_cfetch(s->f, pidx);
  double xc[4];
  gkyl_rect_grid_cell_center(&s->pgrid, idx, xc);
  double dv1 = 2.0 / s->pgrid.dx[cdim];
  int nnodes = 1;
  for (int d = 0; d < pdim; ++d) {
    nnodes *= 6;
  }
  for (int n = 0; n < nnodes; ++n) {
    double z[4], wq = 1.0;
    int m = n;
    for (int d = pdim - 1; d >= 0; --d) {
      z[d] = gl6_x[m % 6];
      wq *= gl6_w[m % 6];
      m /= 6;
    }
    double bq[3], uq[3], dq[3], Bq[3], guq[9], gbq[9], gxq[3], glq[3];
    for (int i = 0; i < 3; ++i) {
      bq[i] = s->cbasis.eval_expand(z, &b[i * nc]);
      uq[i] = s->cbasis.eval_expand(z, &u[i * nc]);
      dq[i] = s->cbasis.eval_expand(z, &dtb[i * nc]);
      Bq[i] = s->cbasis.eval_expand(z, &em[(3 + i) * nc]);
    }
    for (int d = 0; d < cdim; ++d) {
      gxq[d] = s->cbasis.eval_expand(z, &gxi[(d * 3) * nc]);
      glq[d] = s->cbasis.eval_expand(z, &gln[d * nc]);
      for (int j = 0; j < 3; ++j) {
        guq[d * 3 + j] = s->cbasis.eval_expand(z, &gu[(d * 3 + j) * nc]);
        gbq[d * 3 + j] = s->cbasis.eval_expand(z, &gb[(d * 3 + j) * nc]);
      }
    }
    struct pkpm_harm_coeffs hc;
    pkpm_harm_coeffs_set(&hc, cdim, bq, uq, s->cbasis.eval_expand(z, xi), s->cbasis.eval_expand(z, th),
      dq, guq, gbq, gxq, glq, qbym_test, Bq);
    double v = xc[cdim] + 0.5 * s->pgrid.dx[cdim] * z[cdim], Uq[8], flux[8], src[8];
    for (int c = 0; c < 8; ++c) {
      Uq[c] = s->pbasis.eval_expand(z, &f[c * np]);
    }
    pkpm_harm_flux_known(&hc, v, Uq, flux);
    pkpm_harm_source(&hc, v, Uq, src);
    // configuration-space flux couplings in direction d (written out): xi M0_d, xi N_d,
    // xi G (e_d - b b_d), xi (4 theta G - 2 theta^2 F) (e_d - b b_d)
    double xiq = s->cbasis.eval_expand(z, xi), thq = s->cbasis.eval_expand(z, th), T[3][8];
    for (int d = 0; d < cdim; ++d) {
      T[d][0] = xiq * Uq[2 + d];
      T[d][1] = xiq * Uq[5 + d];
      for (int a = 0; a < 3; ++a) {
        double pn = (a == d ? 1.0 : 0.0) - bq[a] * bq[d];
        T[d][2 + a] = xiq * Uq[1] * pn;
        T[d][5 + a] = xiq * (4.0 * thq * Uq[1] - 2.0 * thq * thq * Uq[0]) * pn;
      }
    }
    for (int k = 0; k < np; ++k) {
      double unit[64] = {0.0};
      unit[k] = 1.0;
      double phi = s->pbasis.eval_expand(z, unit);
      double dphi_v = dv1 * s->pbasis.eval_grad_expand(cdim, z, unit);
      for (int c = 0; c < 8; ++c) {
        out_ref[c * np + k] += wq * (phi * src[c] + dphi_v * flux[c]);
        for (int d = 0; d < cdim; ++d) {
          out_ref[c * np + k] +=
            wq * 2.0 / s->pgrid.dx[d] * s->pbasis.eval_grad_expand(d, z, unit) * T[d][c];
        }
      }
    }
  }
}

// Independent evaluation of the stage-A volume terms of one cell, the streaming of every block, by
// 4-point Gauss-Legendre quadrature in every phase-space direction (exact for the products).
static void
reference_vol(struct test_setup *s, const int *idx, double *out_ref)
{
  int cdim = s->cdim, pdim = s->pdim;
  int nc = s->cbasis.num_basis, np = s->pbasis.num_basis;
  long cidx = gkyl_range_idx(&s->clocal_ext, idx), pidx = gkyl_range_idx(&s->plocal_ext, idx);
  const double *b = gkyl_array_cfetch(s->bvar, cidx), *u = gkyl_array_cfetch(s->u, cidx);
  const double *xi = gkyl_array_cfetch(s->xi, cidx);
  const double *f = gkyl_array_cfetch(s->f, pidx);
  double xc[4];
  gkyl_rect_grid_cell_center(&s->pgrid, idx, xc);
  double wv = xc[cdim], dv = s->pgrid.dx[cdim];
  double dx1[4];
  for (int d = 0; d < pdim; ++d) {
    dx1[d] = 2.0 / s->pgrid.dx[d];
  }
  for (int k = 0; k < 8 * np; ++k) {
    out_ref[k] = 0.0;
  }
  int nnodes = 1;
  for (int d = 0; d < pdim; ++d) {
    nnodes *= 4;
  }
  for (int n = 0; n < nnodes; ++n) {
    double z[4], w = 1.0;
    int m = n;
    for (int d = pdim - 1; d >= 0; --d) {
      z[d] = gl4_x[m % 4];
      w *= gl4_w[m % 4];
      m /= 4;
    }
    double v = wv + 0.5 * dv * z[cdim];
    double bq[3], uq[3];
    for (int i = 0; i < 3; ++i) {
      bq[i] = s->cbasis.eval_expand(z, &b[i * nc]);
      uq[i] = s->cbasis.eval_expand(z, &u[i * nc]);
    }
    double xiq = s->cbasis.eval_expand(z, xi);
    for (int k = 0; k < np; ++k) {
      double unit[64] = {0.0};
      unit[k] = 1.0;
      for (int c = 0; c < 8; ++c) {
        double Xq = s->pbasis.eval_expand(z, &f[c * np]);
        for (int d = 0; d < cdim; ++d) {
          out_ref[c * np + k] += w * dx1[d] * s->pbasis.eval_grad_expand(d, z, unit) * (xiq * v * bq[d] + uq[d]) * Xq;
        }
      }
    }
  }
}

static void
test_vol(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int np = s.pbasis.num_basis;

  struct gkyl_dg_eqn *eqn = make_eqn(&s, NULL);
  struct gkyl_array *rhs = mkarr(8 * np, s.plocal_ext.volume);
  gkyl_array_clear(rhs, 0.0);

  double out_ref[8 * 24];
  double max_err = 0.0, max_ref = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    long pidx = gkyl_range_idx(&s.plocal_ext, iter.idx);
    double xc[4];
    gkyl_rect_grid_cell_center(&s.pgrid, iter.idx, xc);
    double *out = gkyl_array_fetch(rhs, pidx);
    double cfl = eqn->vol_term(eqn, xc, s.pgrid.dx, iter.idx, gkyl_array_cfetch(s.f, pidx), out);
    TEST_CHECK(cfl >= 0.0);
    reference_vol(&s, iter.idx, out_ref);
    for (int k = 0; k < 8 * np; ++k) {
      double err = fabs(out[k] - out_ref[k]);
      max_err = err > max_err ? err : max_err;
      max_ref = fabs(out_ref[k]) > max_ref ? fabs(out_ref[k]) : max_ref;
    }
  }
  TEST_CHECK(max_ref > 0.0);
  TEST_CHECK(max_err < 1e-11 * max_ref);
  TEST_MSG("cdim=%d: max volume error %g (max reference %g)", cdim, max_err, max_ref);

  gkyl_array_release(rhs);
  gkyl_dg_eqn_release(eqn);
  teardown(&s);
}

// Full stage-A right-hand side (fluxes + volume) through hyper_dg, then the closures of
// the factored update: the weak residuals <psi_k, M0(f) c_0 + M0(rhs_A)> and
// <psi_k, M0(f) a_P + M1(rhs_A) + c_0 M1(f)> vanish in every mode (f need not be normalized
// for this: the closures are weak divisions by its zeroth moment).
// The modal first-harmonic volume integrals against the exact quadrature of the pointwise
// coefficients, every block of every phase cell.
static void
test_harm_vol(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int np = s.pbasis.num_basis;
  struct gkyl_dg_pkpm_harm_vol_inp hinp = {
    .phase_grid = &s.pgrid, .conf_basis = &s.cbasis, .phase_basis = &s.pbasis, .qbym = qbym_test, .pos_map = s.pos_map, .use_gpu = false};
  struct gkyl_dg_pkpm_harm_vol *harm = gkyl_dg_pkpm_harm_vol_inew(&hinp);
  struct gkyl_array *rhs = mkarr(8 * np, s.plocal_ext.volume), *cflrate = mkarr(1, s.plocal_ext.volume);
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  {
    struct gkyl_array *hcoef = mkarr(gkyl_dg_pkpm_harm_vol_num_coef(harm), s.clocal_ext.volume);
    struct gkyl_array *hwn = mkarr(gkyl_dg_pkpm_harm_vol_num_wnode(harm), s.clocal_ext.volume);
    gkyl_dg_pkpm_harm_vol_coef(harm, &s.clocal, s.bvar, s.grad_b, s.u, s.grad_u, s.xi,
      s.grad_xi, s.grad_ln_n, s.theta, s.dtb, s.em, hcoef, hwn);
    gkyl_dg_pkpm_harm_vol_advance(harm, &s.clocal, &s.plocal, hcoef, s.f, cflrate, rhs);
    // the column source bound (plan 17.27, 17.44) against the Gershgorin row sums of the stored source
    // entries plus |Omega| from the gyration entries, sampled over the cell (corners included) and
    // over |w| <= W
    double W = harm->wmax, max_ratio = 0.0, om_max = 0.0;
    int ns = 1;
    for (int d = 0; d < cdim; ++d) {
      ns *= 5;
    }
    struct gkyl_range_iter cit;
    gkyl_range_iter_init(&cit, &s.clocal);
    while (gkyl_range_iter_next(&cit)) {
      const double *ho = gkyl_array_cfetch(hcoef, gkyl_range_idx(&s.clocal_ext, cit.idx));
      double bound = pkpm_hm_src_bound(harm, ho);
      for (int m = 0; m < ns; ++m) {
        // the degree <= 2 modes at the sample point (1D orthonormal Legendre per direction)
        double ph[3][3], mode[27];
        for (int d = 0, mm = m; d < cdim; ++d, mm /= 5) {
          double z = -1.0 + 0.5 * (mm % 5);
          ph[d][0] = 0.7071067811865475;
          ph[d][1] = 1.224744871391589 * z;
          ph[d][2] = 0.7905694150420949 * (3.0 * z * z - 1.0);
        }
        for (int k = 0; k < harm->n2; ++k) {
          mode[k] = 1.0;
          for (int d = 0, kk = k; d < cdim; ++d, kk /= 3) {
            mode[k] *= ph[d][kk % 3];
          }
        }
        double om2 = 0.0;
        for (int e = 0; e < 3; ++e) {
          double v = 0.0;
          for (int k = 0; k < harm->n2; ++k) {
            v += ho[PKPM_HO_GYR(e) * harm->n2 + k] * mode[k];
          }
          om2 += v * v;
        }
        om_max = fmax(om_max, sqrt(om2));
        for (int iw = 0; iw <= 8; ++iw) {
          double w = -W + 0.25 * W * iw, row[8] = {0.0};
          for (int p = 0; p < harm->nsrc; ++p) {
            double v = 0.0, wj = 1.0;
            for (int j = 0; j < 3; ++j, wj *= w) {
              if (harm->src_ent[p][j] >= 0) {
                const double *c = &ho[harm->src_ent[p][j] * harm->n2];
                for (int k = 0; k < harm->n2; ++k) {
                  v += wj * c[k] * mode[k];
                }
              }
            }
            row[harm->src_o[p]] += fabs(v);
          }
          for (int o = 0; o < 8; ++o) {
            max_ratio = fmax(max_ratio, (row[o] + sqrt(om2)) / bound);
          }
        }
      }
    }
    TEST_CHECK(max_ratio <= 1.0 + 1e-12);
    TEST_MSG("cdim=%d: largest sampled source row sum + |Omega| / column bound %g", cdim, max_ratio);
    TEST_CHECK(om_max > 0.0);
    gkyl_array_release(hcoef);
    gkyl_array_release(hwn);
  }

  double out_ref[8 * 24], max_err[8] = {0.0}, max_ref[8] = {0.0}, min_cfl = 1e300;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    long pidx = gkyl_range_idx(&s.plocal_ext, iter.idx);
    const double *out = gkyl_array_cfetch(rhs, pidx);
    for (int k = 0; k < 8 * np; ++k) {
      out_ref[k] = 0.0;
    }
    reference_harm_vol(&s, iter.idx, out_ref);
    for (int c = 0; c < 8; ++c) {
      for (int k = 0; k < np; ++k) {
        max_err[c] = fmax(max_err[c], fabs(out[c * np + k] - out_ref[c * np + k]));
        max_ref[c] = fmax(max_ref[c], fabs(out_ref[c * np + k]));
      }
    }
    min_cfl = fmin(min_cfl, ((const double *)gkyl_array_cfetch(cflrate, pidx))[0]);
  }
  double ref_all = 0.0;
  for (int c = 0; c < 8; ++c) {
    ref_all = fmax(ref_all, max_ref[c]);
  }
  for (int c = 0; c < 8; ++c) {
    TEST_CHECK(max_ref[c] > 0.0);
    TEST_CHECK(max_err[c] < 1e-13 * ref_all);
    TEST_MSG("cdim=%d block %d: max error %g (max reference %g, overall %g)", cdim, c, max_err[c], max_ref[c], ref_all);
  }
  TEST_CHECK(min_cfl > 0.0);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_dg_pkpm_harm_vol_release(harm);
  teardown(&s);
}

static void
test_closure(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nc = s.cbasis.num_basis, np = s.pbasis.num_basis;
  // No mirror force, and no first harmonic in the state (int M0_perp dw = 0, the constraint
  // manifold): the rows checked below are then the edge-free, uncoupled ones. The perpendicular
  // rows are still forced, by the known first-harmonic sources of the M0_perp blocks.
  gkyl_array_clear(s.grad_b, 0.0);
  gkyl_array_clear_range(s.f, 0.0, &s.plocal_ext);
  {
    gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s.pgrid, &s.pbasis, 8, 2, eval_f, &cdim);
    gkyl_proj_on_basis_advance(proj, 0.0, &s.plocal_ext_conf, s.f);
    gkyl_proj_on_basis_release(proj);
  }

  // Face fluxes of the eight blocks (the fluid flux and the blend weight are not used here).
  struct gkyl_dg_pkpm_remap_flux_surf_inp finp = {
    .phase_grid = &s.pgrid,
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .pos_map = s.pos_map,
    .num_face_nodes = 3,
    .use_gpu = false,
  };
  struct gkyl_dg_pkpm_remap_flux_surf *calc_flux = gkyl_dg_pkpm_remap_flux_surf_inew(&finp);
  int nk = gkyl_dg_pkpm_remap_flux_surf_num_surf_basis(calc_flux);
  int nk_fluid = gkyl_dg_pkpm_remap_flux_surf_num_surf_basis_fluid(calc_flux);
  struct gkyl_array *fluid = mkarr(5 * nc, s.clocal_ext.volume);
  {
    gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s.cgrid, &s.cbasis, 8, 5, eval_fluid, &cdim);
    gkyl_proj_on_basis_advance(proj, 0.0, &s.clocal_ext, fluid);
    gkyl_proj_on_basis_release(proj);
  }
  struct gkyl_array *flux = mkarr(2 * 8 * cdim * nk, s.plocal_ext.volume);
  struct gkyl_array *fluid_flux = mkarr(5 * cdim * nk_fluid, s.clocal_ext.volume);
  struct gkyl_array *cflrate = mkarr(1, s.plocal_ext.volume);
  struct gkyl_array *cflrate_fluid = mkarr(1, s.clocal_ext.volume);
  struct gkyl_array *blend = mkarr(cdim, s.clocal_ext.volume);
  gkyl_array_clear(flux, 0.0);
  gkyl_array_clear(fluid_flux, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_array_clear(cflrate_fluid, 0.0);
  gkyl_array_clear(blend, 0.0);
  gkyl_dg_pkpm_remap_flux_surf_advance(
    calc_flux, &s.clocal, &s.clocal_ext, &s.plocal, &s.plocal_ext, s.bvar, fluid, s.u, s.xi,
    s.theta, s.f, flux, fluid_flux, cflrate, cflrate_fluid, blend
  );

  struct gkyl_dg_eqn *eqn = make_eqn(&s, flux);
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux[2 * GKYL_MAX_DIM] = {0};
  gkyl_hyper_dg *slvr = gkyl_hyper_dg_new(&s.pgrid, &s.pbasis, eqn, cdim, up_dirs, zero_flux, 1, false);
  struct gkyl_array *rhs = mkarr(8 * np, s.plocal_ext.volume);
  gkyl_array_clear(rhs, 0.0);
  gkyl_hyper_dg_advance(slvr, &s.plocal, s.f, cflrate, rhs);
  struct gkyl_dg_pkpm_harm_vol_inp hinp = {
    .phase_grid = &s.pgrid, .conf_basis = &s.cbasis, .phase_basis = &s.pbasis, .qbym = qbym_test, .pos_map = s.pos_map, .use_gpu = false};
  struct gkyl_dg_pkpm_harm_vol *harm = gkyl_dg_pkpm_harm_vol_inew(&hinp);
  struct gkyl_array *hcoef = mkarr(gkyl_dg_pkpm_harm_vol_num_coef(harm), s.clocal_ext.volume);
  struct gkyl_array *hwn = mkarr(gkyl_dg_pkpm_harm_vol_num_wnode(harm), s.clocal_ext.volume);
  gkyl_dg_pkpm_harm_vol_coef(harm, &s.clocal, s.bvar, s.grad_b, s.u, s.grad_u, s.xi, s.grad_xi, s.grad_ln_n,
    s.theta, s.dtb, s.em, hcoef, hwn);
  gkyl_dg_pkpm_harm_vol_advance(harm, &s.clocal, &s.plocal, hcoef, s.f, cflrate, rhs);
  gkyl_dg_pkpm_harm_vol_release(harm);

  // Moments of the stage state and of the partial right-hand side.
  double mass = 2.0;
  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_pkpm_new(
    &s.pgrid, &s.cbasis, &s.pbasis, &s.clocal, &s.vlocal, mass, false, false
  );
  struct gkyl_array *moms = mkarr(18 * nc, s.clocal_ext.volume);
  struct gkyl_array *rhs_moms = mkarr(18 * nc, s.clocal_ext.volume);
  gkyl_array_clear(moms, 0.0);
  gkyl_array_clear(rhs_moms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(mcalc, &s.plocal, &s.clocal, s.f, moms);
  gkyl_dg_updater_moment_pkpm_advance(mcalc, &s.plocal, &s.clocal, rhs, rhs_moms);

  struct gkyl_dg_calc_pkpm_vars *vars =
    gkyl_dg_calc_pkpm_vars_new(&s.cgrid, &s.cbasis, s.pos_map, false);
  struct gkyl_array *a_P = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *c_h = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *c_0 = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *a_perp = mkarr(3 * nc, s.clocal_ext.volume);
  // Previous closures alpha_0 = 0, c = -100 make both velocity edges inflow at every node
  // (the known forces are far smaller than 100 w_max), so the edge fluxes vanish and the
  // closure is the edge-free solve whose weak residuals are checked below.
  gkyl_array_clear(a_P, 0.0);
  gkyl_array_clear(c_h, 0.0);
  gkyl_array_shiftc(c_h, -100.0 * pow(sqrt(2.0), cdim), 0);
  gkyl_array_clear(c_0, 0.0);
  struct gkyl_array *vflux = mkarr(8 * nc, s.plocal_ext.volume);
  gkyl_array_clear(vflux, 0.0);
  int num_null = gkyl_dg_calc_pkpm_vars_closure(
    vars, &s.clocal, &s.plocal, &s.pgrid, moms, rhs_moms, s.bvar, hwn, s.theta, s.f, c_0, a_P, c_h, a_perp,
    vflux, cflrate
  );
  gkyl_array_release(hcoef);
  gkyl_array_release(hwn);
  // both edge slots (lowest cell, velocity ghost cell) must hold no flux
  double edge_max = 0.0;
  struct gkyl_range_iter eit;
  gkyl_range_iter_init(&eit, &s.clocal);
  while (gkyl_range_iter_next(&eit)) {
    int idx[GKYL_MAX_DIM];
    gkyl_copy_int_arr(cdim, eit.idx, idx);
    for (int e = 0; e < 2; ++e) {
      idx[cdim] = e == 0 ? s.plocal.lower[cdim] : s.plocal.upper[cdim] + 1;
      const double *fe = gkyl_array_cfetch(vflux, gkyl_range_idx(&s.plocal_ext, idx));
      for (int k = 0; k < 8 * nc; ++k) edge_max = fmax(edge_max, fabs(fe[k]));
    }
  }
  TEST_CHECK(edge_max == 0.0);
  gkyl_array_release(vflux);
  TEST_CHECK(num_null == 0);

  // Weak residuals in every mode: M0 c_0 + r0, M0 alpha_0 + M1 c + r1 + c_0 M1 and
  // 2 M1 alpha_0 + 2 (M2 + 2 M0G) c + r2 + 2 r0G + c_0 (M2 + 2 M0G).
  struct gkyl_array *n_c0 = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *n_aP = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *c0_m1 = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *ch_m1 = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *m1_aP = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *s2 = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *ch_s2 = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *c0_m2 = mkarr(nc, s.clocal_ext.volume);
  gkyl_array_clear(n_c0, 0.0);
  gkyl_array_clear(n_aP, 0.0);
  gkyl_array_clear(c0_m1, 0.0);
  gkyl_array_clear(ch_m1, 0.0);
  gkyl_array_clear(m1_aP, 0.0);
  gkyl_array_clear(ch_s2, 0.0);
  gkyl_array_clear(c0_m2, 0.0);
  gkyl_array_set_offset(s2, 1.0, moms, nc);
  gkyl_array_accumulate_offset(s2, 2.0, moms, 2 * nc);
  gkyl_dg_mul_op_range(&s.cbasis, 0, n_c0, 0, moms, 0, c_0, &s.clocal);
  gkyl_dg_mul_op_range(&s.cbasis, 0, n_aP, 0, moms, 0, a_P, &s.clocal);
  gkyl_dg_mul_op_range(&s.cbasis, 0, c0_m1, 0, c_0, 3, moms, &s.clocal);
  gkyl_dg_mul_op_range(&s.cbasis, 0, ch_m1, 0, c_h, 3, moms, &s.clocal);
  gkyl_dg_mul_op_range(&s.cbasis, 0, m1_aP, 0, a_P, 3, moms, &s.clocal);
  gkyl_dg_mul_op_range(&s.cbasis, 0, ch_s2, 0, c_h, 0, s2, &s.clocal);
  gkyl_dg_mul_op_range(&s.cbasis, 0, c0_m2, 0, c_0, 0, s2, &s.clocal);
  double max_res = 0.0, max_m1 = 0.0, max_aP = 0.0, max_res0 = 0.0, max_m0 = 0.0;
  double max_res2 = 0.0, max_m2 = 0.0, max_ch = 0.0;
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, &s.clocal);
  while (gkyl_range_iter_next(&citer)) {
    long cidx = gkyl_range_idx(&s.clocal_ext, citer.idx);
    const double *na = gkyl_array_cfetch(n_aP, cidx);
    const double *nc0 = gkyl_array_cfetch(n_c0, cidx);
    const double *cm1 = gkyl_array_cfetch(c0_m1, cidx), *chm1 = gkyl_array_cfetch(ch_m1, cidx);
    const double *m1a = gkyl_array_cfetch(m1_aP, cidx), *chs = gkyl_array_cfetch(ch_s2, cidx);
    const double *cm2 = gkyl_array_cfetch(c0_m2, cidx);
    const double *rm = gkyl_array_cfetch(rhs_moms, cidx);
    const double *ap = gkyl_array_cfetch(a_P, cidx), *ch = gkyl_array_cfetch(c_h, cidx);
    for (int k = 0; k < nc; ++k) {
      double res0 = fabs(nc0[k] + rm[k]);
      max_res0 = res0 > max_res0 ? res0 : max_res0;
      max_m0 = fabs(rm[k]) > max_m0 ? fabs(rm[k]) : max_m0;
      double res = fabs(na[k] + chm1[k] + rm[3 * nc + k] + cm1[k]);
      max_res = res > max_res ? res : max_res;
      max_m1 = fabs(rm[3 * nc + k]) > max_m1 ? fabs(rm[3 * nc + k]) : max_m1;
      max_aP = fabs(ap[k]) > max_aP ? fabs(ap[k]) : max_aP;
      double r2 = rm[nc + k] + 2.0 * rm[2 * nc + k];
      double res2 = fabs(2.0 * m1a[k] + 2.0 * chs[k] + r2 + cm2[k]);
      max_res2 = res2 > max_res2 ? res2 : max_res2;
      max_m2 = fabs(r2) > max_m2 ? fabs(r2) : max_m2;
      max_ch = fabs(ch[k]) > max_ch ? fabs(ch[k]) : max_ch;
    }
  }
  TEST_CHECK(max_m0 > 0.0);
  TEST_CHECK(max_res0 < 1e-12 * max_m0);
  TEST_MSG("cdim=%d: max normalization residual %g (max M0(rhs) %g)", cdim, max_res0, max_m0);
  gkyl_array_release(n_c0);
  gkyl_array_release(n_aP);
  gkyl_array_release(c0_m1);
  gkyl_array_release(ch_m1);
  gkyl_array_release(m1_aP);
  gkyl_array_release(s2);
  gkyl_array_release(ch_s2);
  gkyl_array_release(c0_m2);
  gkyl_array_release(c_0);
  TEST_CHECK(max_m1 > 0.0);
  TEST_CHECK(max_aP > 0.0);
  TEST_CHECK(max_res < 1e-12 * max_m1);
  TEST_MSG("cdim=%d: first-moment residual %g (max |M1(rhs)| %g, max |alpha_0| %g)", cdim, max_res, max_m1, max_aP);
  TEST_CHECK(max_m2 > 0.0);
  TEST_CHECK(max_ch > 0.0);
  TEST_CHECK(max_res2 < 1e-12 * max_m2);
  TEST_MSG("cdim=%d: second-moment residual %g (max |M2 + 2 M0G (rhs)| %g, max |c| %g)", cdim, max_res2, max_m2, max_ch);

  // Perpendicular rows at the closure nodes (where b is a unit vector): M0 a_perp + I_perp rho_M = 0,
  // rho_M = int rhs_A(M0_perp) dw.
  double max_res_p = 0.0, max_rho_p = 0.0;
  int nnodes = 1;
  for (int d = 0; d < cdim; ++d) {
    nnodes *= 2;
  }
  gkyl_range_iter_init(&citer, &s.clocal);
  while (gkyl_range_iter_next(&citer)) {
    long cidx = gkyl_range_idx(&s.clocal_ext, citer.idx);
    const double *m = gkyl_array_cfetch(moms, cidx), *rm = gkyl_array_cfetch(rhs_moms, cidx);
    const double *ap = gkyl_array_cfetch(a_perp, cidx), *b = gkyl_array_cfetch(s.bvar, cidx);
    for (int n = 0; n < nnodes; ++n) {
      double z[3];
      for (int d = 0; d < cdim; ++d) {
        z[d] = ((n >> (cdim - 1 - d)) & 1) ? 1.0 / sqrt(3.0) : -1.0 / sqrt(3.0);
      }
      double bq[3], rho[3], rho_p[3];
      for (int a = 0; a < 3; ++a) {
        bq[a] = s.cbasis.eval_expand(z, &b[a * nc]);
        rho[a] = s.cbasis.eval_expand(z, &rm[(6 + a) * nc]);
      }
      pkpm_harm_perp(bq, rho, rho_p);
      double m0 = s.cbasis.eval_expand(z, &m[0]);
      for (int a = 0; a < 3; ++a) {
        max_res_p = fmax(max_res_p, fabs(m0 * s.cbasis.eval_expand(z, &ap[a * nc]) + rho_p[a]));
        max_rho_p = fmax(max_rho_p, fabs(rho_p[a]));
      }
    }
  }
  TEST_CHECK(max_rho_p > 0.0);
  TEST_CHECK(max_res_p < 1e-12 * max_rho_p);
  TEST_MSG("cdim=%d: perpendicular-row residual %g (max |I_perp int rhs(M0_perp) dw| %g)", cdim, max_res_p, max_rho_p);
  gkyl_array_release(a_perp);

  gkyl_array_release(a_P);
  gkyl_array_release(c_h);
  gkyl_dg_calc_pkpm_vars_release(vars);
  gkyl_array_release(rhs_moms);
  gkyl_array_release(moms);
  gkyl_dg_updater_moment_pkpm_release(mcalc);
  gkyl_array_release(rhs);
  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);
  gkyl_array_release(cflrate);
  gkyl_array_release(cflrate_fluid);
  gkyl_array_release(flux);
  gkyl_array_release(fluid_flux);
  gkyl_array_release(fluid);
  gkyl_array_release(blend);
  gkyl_dg_pkpm_remap_flux_surf_release(calc_flux);
  teardown(&s);
}

// The weak gradient kernels against an independent quadrature of the weak form
//   <psi_k, d_i u_j> = dx1_i [ <psi_k(xi_i=+1), ubar_r> - <psi_k(xi_i=-1), ubar_l> - <d_xi_i psi_k, u_c> ]
// with ubar the average of the two traces at the face.
static void
test_grad_u(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nc = s.cbasis.num_basis;
  double max_err = 0.0, max_val = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.clocal);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&s.clocal_ext, iter.idx);
    const double *gu = gkyl_array_cfetch(s.grad_u, cidx), *uc = gkyl_array_cfetch(s.u, cidx);
    for (int dir = 0; dir < cdim; ++dir) {
      int idx_l[3], idx_r[3];
      gkyl_copy_int_arr(cdim, iter.idx, idx_l);
      gkyl_copy_int_arr(cdim, iter.idx, idx_r);
      idx_l[dir] -= 1;
      idx_r[dir] += 1;
      const double *ul = gkyl_array_cfetch(s.u, gkyl_range_idx(&s.clocal_ext, idx_l));
      const double *ur = gkyl_array_cfetch(s.u, gkyl_range_idx(&s.clocal_ext, idx_r));
      double dx1 = 2.0 / s.cgrid.dx[dir];
      for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < nc; ++k) {
          double unit[8] = {0.0};
          unit[k] = 1.0;
          double ref = 0.0;
          // Volume term: 4-point Gauss-Legendre in every direction.
          int nv = 1;
          for (int d = 0; d < cdim; ++d) {
            nv *= 4;
          }
          for (int n = 0; n < nv; ++n) {
            double z[3], w = 1.0;
            int m = n;
            for (int d = cdim - 1; d >= 0; --d) {
              z[d] = gl4_x[m % 4];
              w *= gl4_w[m % 4];
              m /= 4;
            }
            ref -= w * dx1 * s.cbasis.eval_grad_expand(dir, z, unit) * s.cbasis.eval_expand(z, &uc[j * nc]);
          }
          // Face terms: 4-point rule in the transverse directions.
          int nf = 1;
          for (int d = 0; d < cdim - 1; ++d) {
            nf *= 4;
          }
          for (int n = 0; n < nf; ++n) {
            double zp[3], zm[3], w = 1.0;
            int m = n, t = 0;
            double zt[3];
            for (int d = cdim - 2; d >= 0; --d) {
              zt[d] = gl4_x[m % 4];
              w *= gl4_w[m % 4];
              m /= 4;
            }
            for (int d = 0; d < cdim; ++d) {
              if (d == dir) {
                zp[d] = 1.0;
                zm[d] = -1.0;
              } else {
                zp[d] = zm[d] = zt[t++];
              }
            }
            // Upper face of this cell: average of this cell at +1 and the right cell at -1.
            double ubar_r = 0.5 * (s.cbasis.eval_expand(zp, &uc[j * nc]) + s.cbasis.eval_expand(zm, &ur[j * nc]));
            double ubar_l = 0.5 * (s.cbasis.eval_expand(zp, &ul[j * nc]) + s.cbasis.eval_expand(zm, &uc[j * nc]));
            ref += w * dx1 * (s.cbasis.eval_expand(zp, unit) * ubar_r - s.cbasis.eval_expand(zm, unit) * ubar_l);
          }
          double val = gu[(dir * 3 + j) * nc + k];
          max_err = fmax(max_err, fabs(val - ref));
          max_val = fmax(max_val, fabs(ref));
        }
      }
    }
  }
  TEST_CHECK(max_val > 0.0);
  TEST_CHECK(max_err < 1e-11 * max_val);
  TEST_MSG("cdim=%d: weak gradient max error %g (max value %g)", cdim, max_err, max_val);
  teardown(&s);
}

void test_grad_u_1x() { test_grad_u(1); }
void test_grad_u_2x() { test_grad_u(2); }
void test_grad_u_3x() { test_grad_u(3); }

void test_vol_1x() { test_vol(1); }
void test_vol_2x() { test_vol(2); }
void test_vol_3x() { test_vol(3); }
void test_harm_vol_1x() { test_harm_vol(1); }
void test_harm_vol_2x() { test_harm_vol(2); }
void test_harm_vol_3x() { test_harm_vol(3); }
void test_closure_1x() { test_closure(1); }
void test_closure_2x() { test_closure(2); }
void test_closure_3x() { test_closure(3); }

TEST_LIST = {
  {"test_grad_u_1x", test_grad_u_1x},
  {"test_grad_u_2x", test_grad_u_2x},
  {"test_grad_u_3x", test_grad_u_3x},
  {"test_vol_1x", test_vol_1x},
  {"test_vol_2x", test_vol_2x},
  {"test_vol_3x", test_vol_3x},
  {"test_harm_vol_1x", test_harm_vol_1x},
  {"test_harm_vol_2x", test_harm_vol_2x},
  {"test_harm_vol_3x", test_harm_vol_3x},
  {"test_closure_1x", test_closure_1x},
  {"test_closure_2x", test_closure_2x},
  {"test_closure_3x", test_closure_3x},
  {NULL, NULL},
};
