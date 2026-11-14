// The Jacobian-weighted (position-mapped) form of the PKPM configuration-space operators.
// With the C^0 linear map x_d = x_d(xi_d) of gkyl_vlasov_position_map the Jacobian
// J = prod_d J_d is constant per cell and the solver evolves J [F, G] and J [rho, M, E]:
//   d_t (J q) + sum_d d_xi_d [ (J/J_d) F_d(q) ] = J S.
// Two checks of the configuration-space operators (weak gradients, HLL face fluxes,
// absorbing-wall fluxes, kinetic configuration pass and fluid right-hand sides):
//   1. A constant stretch x_d = s_d xi_d (different in every direction) is the uniform
//      physical grid in disguise: the modal coefficients of every field agree with those of
//      the uniform grid (times J for the extensive ones), so the weak gradients must be equal,
//      every stored face flux J/J_d times the physical one, the cfl rates equal, and the
//      right-hand sides J times the physical ones.
//   2. A smooth periodic map x = xi + a sin(xi) with J varying from cell to cell: the face
//      factor J/J_d is single-valued, so the cell-averaged right-hand sides of the conserved
//      components (F without its normalization sources, rho, M, E) sum to zero on the
//      periodic domain.
#include <math.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_pkpm_remap_flux_surf.h>
#include <gkyl_dg_pkpm_fluid.h>
#include <gkyl_dg_pkpm_harm_vol.h>
#include <gkyl_dg_pkpm_kinetic_conf.h>
#include <gkyl_evalf_def.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>

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

// Smooth periodic fields on the physical domain [-pi, pi]^cdim (xn are physical
// coordinates; the projection maps the quadrature nodes through the position map).
static void
eval_em(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.4, z = cdim > 2 ? xn[2] : -0.3;
  fout[0] = 0.1 * sin(x);
  fout[1] = 0.2 * cos(y);
  fout[2] = 0.05;
  fout[3] = 1.0 + 0.3 * cos(x) * sin(y);
  fout[4] = 0.4 * sin(x);
  fout[5] = 0.5 + 0.2 * cos(z);
  fout[6] = fout[7] = 0.0;
}

static void
eval_u(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.4, z = cdim > 2 ? xn[2] : -0.3;
  fout[0] = 0.3 * sin(x) + 0.1 * cos(y);
  fout[1] = 0.2 * cos(x);
  fout[2] = -0.1 + 0.1 * sin(z);
}

static double
xi_phys(const double *xn, int cdim)
{
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.4, z = cdim > 2 ? xn[2] : -0.3;
  return 1.0 + 0.2 * cos(x) * sin(y) + 0.05 * cos(z);
}

static double
rho_phys(const double *xn, int cdim)
{
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.4;
  return 1.0 + 0.3 * cos(x) + 0.1 * sin(y);
}

// Thermal speed [xi, xi^2, p].
static void
eval_xi(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double xi = xi_phys(xn, cdim);
  fout[0] = xi;
  fout[1] = xi * xi;
  fout[2] = 0.9 * rho_phys(xn, cdim) * xi * xi;
}

static void
eval_theta(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.4;
  fout[0] = 0.8 + 0.2 * sin(x) * cos(y);
}

// Physical [F, G].
static void
eval_f(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.4, w = xn[cdim];
  double F = (1.0 + 0.2 * sin(x)) * (1.0 + 0.1 * cos(y)) * (1.0 + 0.1 * w) * exp(-0.5 * w * w)
    / sqrt(2.0 * M_PI);
  fout[0] = F;
  fout[1] = 0.7 * F;
}

// Physical fluid state [rho, M, E] with p_par = rho xi^2, p_perp = 0.8 rho xi^2.
static void
eval_fluid(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double rho = rho_phys(xn, cdim), xi = xi_phys(xn, cdim), u[3];
  eval_u(t, xn, u, ctx);
  double uu = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
  fout[0] = rho;
  for (int i = 0; i < 3; ++i) {
    fout[1 + i] = rho * u[i];
  }
  fout[4] = 0.5 * rho * uu + 0.5 * rho * xi * xi + 0.8 * rho * xi * xi;
}

// Physical kinetic coupling moments [p_par, p_perp, q_par, q_perp].
static void
eval_cmoms(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double rho = rho_phys(xn, cdim), xi = xi_phys(xn, cdim);
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.4;
  fout[0] = rho * xi * xi;
  fout[1] = 0.8 * rho * xi * xi;
  fout[2] = 0.1 * sin(x);
  fout[3] = 0.05 * cos(y);
  // first-harmonic stress s and heat flux r
  fout[4] = 0.05 * cos(x);
  fout[5] = -0.03 * sin(y);
  fout[6] = 0.02 * sin(x + y);
  fout[7] = 0.04 * cos(x - y);
  fout[8] = 0.03 * sin(x);
  fout[9] = -0.02 * cos(y);
}

// Position maps: a constant stretch x = s xi (s from the context) and the sinusoid
// x = xi + a sin(xi), which fixes the end points of [-pi, pi] and has a periodic Jacobian.
static void
pmap_stretch(double t, const double *xc, double *xp, void *ctx)
{
  xp[0] = *(double *)ctx * xc[0];
}

static void
pmap_sine(double t, const double *xc, double *xp, void *ctx)
{
  xp[0] = xc[0] + 0.4 * sin(xc[0]);
}

static void
pos_c2p(const double *xcomp, double *xphys, void *ctx)
{
  gkyl_vlasov_position_map_eval_mc2p(ctx, xcomp, xphys);
}

enum map_kind { MAP_IDENTITY, MAP_STRETCH, MAP_SINE };

struct sim {
  int cdim, pdim;
  double scale[3];
  struct gkyl_rect_grid cgrid, pgrid;
  struct gkyl_basis cbasis, pbasis;
  struct gkyl_range clocal, clocal_ext, plocal, plocal_ext, plocal_ext_conf;
  struct gkyl_vlasov_position_map *pos_map;
  // intensive fields
  struct gkyl_array *em, *dtb, *bvar, *grad_b, *u, *grad_u, *xi, *grad_xi, *rho, *rho3, *grad_rho, *grad_ln_n, *theta;
  // physical fluid state and the J-weighted evolved quantities
  struct gkyl_array *fluid_phys, *fluid, *cmoms, *f;
  struct gkyl_dg_calc_pkpm_bvar *calc_b;
  struct gkyl_dg_calc_pkpm_vars *vars;
  struct gkyl_dg_pkpm_remap_flux_surf *flux_up;
  struct gkyl_array *kin_flux, *kin_cfl, *kin_rhs, *fl_flux, *fl_cfl, *fl_rhs, *blend;
};

// Phase-space ranges with the app's layout: conf ghost cells and one velocity ghost cell
// above the upper velocity edge.
static void
phase_ranges(struct sim *s)
{
  int pdim = s->pdim, cdim = s->cdim, lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int elo[GKYL_MAX_DIM] = {0}, eup[GKYL_MAX_DIM] = {0};
  for (int d = 0; d < pdim; ++d) {
    lower[d] = 1;
    upper[d] = s->pgrid.cells[d];
  }
  for (int d = 0; d < cdim; ++d) {
    elo[d] = eup[d] = 1;
  }
  eup[cdim] = 1;
  struct gkyl_range local;
  gkyl_range_init(&local, pdim, lower, upper);
  gkyl_range_extend(&s->plocal_ext, &local, elo, eup);
  gkyl_sub_range_init(&s->plocal, &s->plocal_ext, lower, upper);
  int lo_c[GKYL_MAX_DIM], up_c[GKYL_MAX_DIM];
  for (int d = 0; d < pdim; ++d) {
    lo_c[d] = d < cdim ? s->plocal_ext.lower[d] : lower[d];
    up_c[d] = d < cdim ? s->plocal_ext.upper[d] : upper[d];
  }
  gkyl_sub_range_init(&s->plocal_ext_conf, &s->plocal_ext, lo_c, up_c);
}

// Periodic ghost cells of a configuration or phase-space array with one ghost layer in
// every configuration direction: the skin and ghost layers of direction d span the full
// extended range in the other directions, so sequential sweeps fill the corners.
static void
apply_periodic(struct gkyl_array *arr, const struct gkyl_range *ext, int cdim)
{
  for (int d = 0; d < cdim; ++d) {
    int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
    struct gkyl_range ghost_lo, ghost_up, skin_lo, skin_up;
    gkyl_copy_int_arr(ext->ndim, ext->lower, lo);
    gkyl_copy_int_arr(ext->ndim, ext->upper, up);
    up[d] = ext->lower[d];
    gkyl_sub_range_init(&ghost_lo, ext, lo, up);
    lo[d] = up[d] = ext->lower[d] + 1;
    gkyl_sub_range_init(&skin_lo, ext, lo, up);
    lo[d] = up[d] = ext->upper[d] - 1;
    gkyl_sub_range_init(&skin_up, ext, lo, up);
    lo[d] = up[d] = ext->upper[d];
    gkyl_sub_range_init(&ghost_up, ext, lo, up);
    gkyl_array_copy_range_to_range(arr, arr, &ghost_lo, &skin_up);
    gkyl_array_copy_range_to_range(arr, arr, &ghost_up, &skin_lo);
  }
}

static void
project_conf(struct sim *s, int nret, evalf_t eval, struct gkyl_array *out)
{
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &s->cgrid,
    .basis = &s->cbasis,
    .num_quad = 8,
    .num_ret_vals = nret,
    .eval = eval,
    .ctx = &s->cdim,
    .c2p_func = pos_c2p,
    .c2p_func_ctx = s->pos_map,
  });
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal, out);
  gkyl_proj_on_basis_release(proj);
  apply_periodic(out, &s->clocal_ext, s->cdim);
}

static void
setup(struct sim *s, int cdim, enum map_kind kind)
{
  int pdim = cdim + 1;
  s->cdim = cdim;
  s->pdim = pdim;
  double stretch[3] = {1.5, 0.75, 2.5};
  double lower[4] = {-M_PI, -M_PI, -M_PI, -4.0}, upper[4] = {M_PI, M_PI, M_PI, 4.0};
  int cells[4] = {6, 4, 3, 8};
  int cghost[3] = {1, 1, 1};
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = {0};
  for (int d = 0; d < cdim; ++d) {
    s->scale[d] = 1.0;
    if (kind == MAP_STRETCH) {
      s->scale[d] = stretch[d];
      lower[d] /= stretch[d];
      upper[d] /= stretch[d];
      inp_pmap[d].eval_pmap = pmap_stretch;
      inp_pmap[d].ctx = &s->scale[d];
    } else if (kind == MAP_SINE) {
      inp_pmap[d].eval_pmap = pmap_sine;
    }
  }
  lower[cdim] = -4.0;
  upper[cdim] = 4.0;
  cells[cdim] = 8;
  gkyl_rect_grid_init(&s->cgrid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&s->pgrid, pdim, lower, upper, cells);
  gkyl_cart_modal_serendip(&s->cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&s->pbasis, cdim, 1);
  gkyl_create_grid_ranges(&s->cgrid, cghost, &s->clocal_ext, &s->clocal);
  phase_ranges(s);
  s->pos_map = gkyl_vlasov_position_map_new(&s->cgrid, &s->clocal, &s->clocal_ext, &s->cbasis, inp_pmap, false);

  int nc = s->cbasis.num_basis, np = s->pbasis.num_basis;
  long cvol = s->clocal_ext.volume, pvol = s->plocal_ext.volume;
  s->em = mkarr(8 * nc, cvol);
  s->bvar = mkarr(3 * nc, cvol);
  s->grad_b = mkarr(cdim * 3 * nc, cvol);
  s->u = mkarr(3 * nc, cvol);
  s->grad_u = mkarr(cdim * 3 * nc, cvol);
  s->xi = mkarr(3 * nc, cvol);
  s->grad_xi = mkarr(cdim * 3 * nc, cvol);
  s->rho = mkarr(nc, cvol);
  s->rho3 = mkarr(3 * nc, cvol);
  s->grad_rho = mkarr(cdim * 3 * nc, cvol);
  s->grad_ln_n = mkarr(cdim * nc, cvol);
  s->theta = mkarr(nc, cvol);
  s->fluid_phys = mkarr(5 * nc, cvol);
  s->fluid = mkarr(5 * nc, cvol);
  s->cmoms = mkarr(10 * nc, cvol);
  s->dtb = mkarr(3 * nc, cvol); // a static field
  gkyl_array_clear(s->dtb, 0.0);
  s->f = mkarr(8 * np, pvol); // [F, G, M0_perp, N_perp]; the first harmonic is zero here
  gkyl_array_clear(s->grad_b, 0.0);
  gkyl_array_clear(s->grad_u, 0.0);
  gkyl_array_clear(s->grad_xi, 0.0);
  gkyl_array_clear(s->rho3, 0.0);
  gkyl_array_clear(s->grad_rho, 0.0);
  gkyl_array_clear(s->grad_ln_n, 0.0);
  gkyl_array_clear(s->f, 0.0);

  project_conf(s, 8, eval_em, s->em);
  project_conf(s, 3, eval_u, s->u);
  project_conf(s, 3, eval_xi, s->xi);
  project_conf(s, 1, eval_theta, s->theta);
  project_conf(s, 5, eval_fluid, s->fluid_phys);
  project_conf(s, 10, eval_cmoms, s->cmoms);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &s->pgrid,
    .basis = &s->pbasis,
    .num_quad = 8,
    .num_ret_vals = 2,
    .eval = eval_f,
    .ctx = &s->cdim,
    .c2p_func = pos_c2p,
    .c2p_func_ctx = s->pos_map,
  });
  gkyl_proj_on_basis_advance(proj, 0.0, &s->plocal, s->f);
  gkyl_proj_on_basis_release(proj);
  apply_periodic(s->f, &s->plocal_ext, cdim);

  // The evolved quantities carry the Jacobian.
  gkyl_vlasov_position_map_rescale_jacobpos_conf(s->pos_map, &s->clocal_ext, s->fluid_phys, s->fluid);
  gkyl_vlasov_position_map_rescale_jacobpos_conf(s->pos_map, &s->clocal_ext, s->cmoms, s->cmoms);
  gkyl_vlasov_position_map_rescale_jacobpos(s->pos_map, &s->pbasis, &s->plocal_ext_conf, s->f, s->f);

  // Intensive derived fields: b and the weak gradients of b, u, xi and ln n.
  s->calc_b = gkyl_dg_calc_pkpm_bvar_new(&s->cgrid, &s->cbasis, s->pos_map, false);
  gkyl_dg_calc_pkpm_bvar_advance(s->calc_b, &s->clocal_ext, s->em, s->bvar);
  gkyl_dg_calc_pkpm_bvar_grad(s->calc_b, &s->clocal, s->bvar, s->grad_b);
  s->vars = gkyl_dg_calc_pkpm_vars_new(&s->cgrid, &s->cbasis, s->pos_map, false);
  weak_grad(s->vars, &s->clocal, s->u, s->grad_u);
  weak_grad(s->vars, &s->clocal, s->xi, s->grad_xi);
  gkyl_array_set_offset(s->rho, 1.0, s->fluid_phys, 0);
  gkyl_array_set_offset(s->rho3, 1.0, s->rho, 0);
  weak_grad(s->vars, &s->clocal, s->rho3, s->grad_rho);
  int null_rho = gkyl_dg_calc_pkpm_vars_grad_ln_n(s->vars, &s->clocal, s->rho, s->grad_rho, s->grad_ln_n);
  TEST_CHECK(null_rho == 0);

}

// The configuration-space operators on the current fields: the kinetic configuration
// pass (stored face fluxes, then hyper_dg: volume + scatter), the fluid (stored HLL
// face fluxes, hyper_dg) and the absorbing-wall fluxes of the x walls.
static void
advance(struct sim *s)
{
  int cdim = s->cdim, nc = s->cbasis.num_basis, np = s->pbasis.num_basis;
  long cvol = s->clocal_ext.volume, pvol = s->plocal_ext.volume;
  int cghost[3] = {1, 1, 1};

  // Face fluxes of the eight blocks and of the fluid, which take the physical fluid state and the
  // weighted distribution.
  struct gkyl_dg_pkpm_remap_flux_surf_inp flux_inp = {
    .phase_grid = &s->pgrid,
    .conf_basis = &s->cbasis,
    .phase_basis = &s->pbasis,
    .pos_map = s->pos_map,
    .num_face_nodes = 3,
    .use_gpu = false,
  };
  s->flux_up = gkyl_dg_pkpm_remap_flux_surf_inew(&flux_inp);
  int nkc = gkyl_dg_pkpm_remap_flux_surf_num_surf_basis(s->flux_up);
  int nkf = gkyl_dg_pkpm_remap_flux_surf_num_surf_basis_fluid(s->flux_up);
  s->kin_flux = mkarr(2 * 8 * cdim * nkc, pvol);
  s->kin_cfl = mkarr(1, pvol);
  s->kin_rhs = mkarr(8 * np, pvol);
  s->fl_flux = mkarr(5 * cdim * nkf, cvol);
  s->fl_cfl = mkarr(1, cvol);
  s->fl_rhs = mkarr(5 * nc, cvol);
  s->blend = mkarr(cdim, cvol);
  gkyl_array_clear(s->kin_flux, 0.0);
  gkyl_array_clear(s->kin_cfl, 0.0);
  gkyl_array_clear(s->kin_rhs, 0.0);
  gkyl_array_clear(s->fl_flux, 0.0);
  gkyl_array_clear(s->fl_cfl, 0.0);
  gkyl_array_clear(s->fl_rhs, 0.0);
  gkyl_array_clear(s->blend, 0.0);
  gkyl_dg_pkpm_remap_flux_surf_advance(
    s->flux_up, &s->clocal, &s->clocal_ext, &s->plocal, &s->plocal_ext, s->bvar, s->fluid_phys,
    s->u, s->xi, s->theta, s->f, s->kin_flux, s->fl_flux, s->kin_cfl, s->fl_cfl, s->blend
  );
  struct gkyl_dg_pkpm_kinetic_conf_inp kin_eqn_inp = {
    .conf_basis = &s->cbasis,
    .phase_basis = &s->pbasis,
    .conf_range = &s->clocal,
    .phase_range = &s->plocal,
    .auxfields =
      {
        .bvar = s->bvar,
        .u = s->u,
        .xi = s->xi,
        .conf_flux_surf = s->kin_flux,
      },
    .pos_map = s->pos_map,
    .use_gpu = false,
  };
  struct gkyl_dg_eqn *kin_eqn = gkyl_dg_pkpm_kinetic_conf_inew(&kin_eqn_inp);
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux[2 * GKYL_MAX_DIM] = {0};
  gkyl_hyper_dg *kin_slvr = gkyl_hyper_dg_new(&s->pgrid, &s->pbasis, kin_eqn, cdim, up_dirs, zero_flux, 1, false);
  gkyl_hyper_dg_advance(kin_slvr, &s->plocal, s->f, s->kin_cfl, s->kin_rhs);
  gkyl_hyper_dg_release(kin_slvr);
  gkyl_dg_eqn_release(kin_eqn);
  // the first-harmonic terms of the pass (their configuration-flux couplings carry 1/J_d)
  struct gkyl_dg_pkpm_harm_vol *harm = gkyl_dg_pkpm_harm_vol_inew(&(struct gkyl_dg_pkpm_harm_vol_inp){
    .phase_grid = &s->pgrid, .conf_basis = &s->cbasis, .phase_basis = &s->pbasis, .qbym = 1.0,
    .pos_map = s->pos_map, .use_gpu = false});
  {
    struct gkyl_array *hcoef = mkarr(gkyl_dg_pkpm_harm_vol_num_coef(harm), s->clocal_ext.volume);
    struct gkyl_array *hwn = mkarr(gkyl_dg_pkpm_harm_vol_num_wnode(harm), s->clocal_ext.volume);
    gkyl_dg_pkpm_harm_vol_coef(harm, &s->clocal, s->bvar, s->grad_b, s->u, s->grad_u, s->xi,
      s->grad_xi, s->grad_ln_n, s->theta, s->dtb, s->em, hcoef, hwn);
    gkyl_dg_pkpm_harm_vol_advance(harm, &s->clocal, &s->plocal, hcoef, s->f, s->kin_cfl, s->kin_rhs);
    gkyl_array_release(hcoef);
    gkyl_array_release(hwn);
  }
  gkyl_dg_pkpm_harm_vol_release(harm);

  // Fluid.
  struct gkyl_dg_pkpm_fluid_inp fl_eqn_inp = {
    .conf_basis = &s->cbasis,
    .conf_range = &s->clocal,
    .auxfields = {.vlasov_pkpm_moms = s->cmoms, .bvar = s->bvar, .u = s->u, .fluid_flux_surf = s->fl_flux},
    .pos_map = s->pos_map,
    .use_gpu = false,
  };
  struct gkyl_dg_eqn *fl_eqn = gkyl_dg_pkpm_fluid_inew(&fl_eqn_inp);
  gkyl_hyper_dg *fl_slvr = gkyl_hyper_dg_new(&s->cgrid, &s->cbasis, fl_eqn, cdim, up_dirs, zero_flux, 1, false);
  gkyl_hyper_dg_advance(fl_slvr, &s->clocal, s->fluid, s->fl_cfl, s->fl_rhs);
  gkyl_hyper_dg_release(fl_slvr);
  gkyl_dg_eqn_release(fl_eqn);
}

static void
teardown(struct sim *s)
{
  gkyl_dg_pkpm_remap_flux_surf_release(s->flux_up);
  gkyl_dg_calc_pkpm_vars_release(s->vars);
  gkyl_dg_calc_pkpm_bvar_release(s->calc_b);
  gkyl_vlasov_position_map_release(s->pos_map);
  struct gkyl_array *arrs[] = {s->em, s->dtb, s->bvar, s->grad_b, s->u, s->grad_u, s->xi, s->grad_xi, s->rho, s->rho3,
    s->grad_rho, s->grad_ln_n, s->theta, s->fluid_phys, s->fluid, s->cmoms, s->f, s->kin_flux, s->kin_cfl,
    s->kin_rhs, s->fl_flux, s->fl_cfl, s->fl_rhs, s->blend};
  for (int i = 0; i < sizeof(arrs) / sizeof(arrs[0]); ++i) {
    gkyl_array_release(arrs[i]);
  }
}

// Total Jacobian J and the per-direction J_d of a configuration cell.
static double
jac_tot(const struct sim *s, long cidx)
{
  return ((const double *)gkyl_array_cfetch(s->pos_map->jacob_pos_gauss, cidx))[0];
}

static double
jac_dir(const struct sim *s, long cidx, int dir)
{
  return ((const double *)gkyl_array_cfetch(s->pos_map->jacob_pos, cidx))[2 * dir];
}

// max |b - w * a| and max |a| over range for arrays with the same layout, w a per-cell
// weight: w = 1 (equal), J (extensive right-hand sides), or J/J_d per direction block
// (face fluxes, nblk blocks of blk components).
enum weight_kind { W_ONE, W_JAC, W_JPERP };

static void
compare(
  const struct sim *a, const struct sim *b, const struct gkyl_array *arr_a, const struct gkyl_array *arr_b,
  const struct gkyl_range *range, const struct gkyl_range *range_ext, enum weight_kind kind, int nblk,
  double *max_err, double *max_ref
)
{
  int cdim = a->cdim, ncomp = arr_a->ncomp, blk = ncomp / nblk;
  *max_err = 0.0;
  *max_ref = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(range_ext, iter.idx);
    long cidx = gkyl_range_idx(&b->clocal_ext, iter.idx); // leading cdim indices
    const double *va = gkyl_array_cfetch(arr_a, idx), *vb = gkyl_array_cfetch(arr_b, idx);
    double J = jac_tot(b, cidx);
    for (int k = 0; k < ncomp; ++k) {
      double w = 1.0;
      if (kind == W_JAC) {
        w = J;
      } else if (kind == W_JPERP) {
        w = J / jac_dir(b, cidx, (k / blk) % cdim);
      }
      double err = fabs(vb[k] - w * va[k]);
      *max_err = err > *max_err ? err : *max_err;
      *max_ref = fabs(va[k]) > *max_ref ? fabs(va[k]) : *max_ref;
    }
  }
}

static void
test_stretch(int cdim)
{
  struct sim a, b;
  setup(&a, cdim, MAP_IDENTITY);
  setup(&b, cdim, MAP_STRETCH);
  advance(&a);
  advance(&b);
  TEST_CHECK(a.pos_map->is_identity);
  TEST_CHECK(!b.pos_map->is_identity);
  // The map is what it claims: J_d = s_d in every cell (including ghosts), J their product.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &b.clocal_ext);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&b.clocal_ext, iter.idx);
    double J = 1.0;
    for (int d = 0; d < cdim; ++d) {
      TEST_CHECK(gkyl_compare(jac_dir(&b, cidx, d), b.scale[d], 1e-14));
      J *= b.scale[d];
    }
    TEST_CHECK(gkyl_compare(jac_tot(&b, cidx), J, 1e-14));
  }

  double tol = 1e-10, err, ref;
  struct {
    const char *name;
    const struct gkyl_array *arr_a, *arr_b;
    const struct gkyl_range *range, *range_ext;
    enum weight_kind kind;
    int nblk;
  } checks[] = {
    // intensive inputs agree (the quadrature nodes of the stretched grid map to the physical ones)
    {"u", a.u, b.u, &a.clocal_ext, &a.clocal_ext, W_ONE, 1},
    {"bvar", a.bvar, b.bvar, &a.clocal_ext, &a.clocal_ext, W_ONE, 1},
    // extensive inputs carry J
    {"fluid", a.fluid, b.fluid, &a.clocal_ext, &a.clocal_ext, W_JAC, 1},
    {"f", a.f, b.f, &a.plocal_ext_conf, &a.plocal_ext, W_JAC, 1},
    // weak gradients are intensive
    {"grad_b", a.grad_b, b.grad_b, &a.clocal, &a.clocal_ext, W_ONE, 1},
    {"grad_u", a.grad_u, b.grad_u, &a.clocal, &a.clocal_ext, W_ONE, 1},
    {"grad_xi", a.grad_xi, b.grad_xi, &a.clocal, &a.clocal_ext, W_ONE, 1},
    {"grad_ln_n", a.grad_ln_n, b.grad_ln_n, &a.clocal, &a.clocal_ext, W_ONE, 1},
    // stored face fluxes carry J/J_d of their direction, cfl rates are physical
    {"kinetic flux", a.kin_flux, b.kin_flux, &a.plocal_ext_conf, &a.plocal_ext, W_JPERP, cdim},
    {"kinetic cfl", a.kin_cfl, b.kin_cfl, &a.plocal, &a.plocal_ext, W_ONE, 1},
    {"fluid flux", a.fl_flux, b.fl_flux, &a.clocal_ext, &a.clocal_ext, W_JPERP, cdim},
    {"fluid cfl", a.fl_cfl, b.fl_cfl, &a.clocal, &a.clocal_ext, W_ONE, 1},
    // right-hand sides of the evolved J q
    {"kinetic rhs", a.kin_rhs, b.kin_rhs, &a.plocal, &a.plocal_ext, W_JAC, 1},
    {"fluid rhs", a.fl_rhs, b.fl_rhs, &a.clocal, &a.clocal_ext, W_JAC, 1},
  };
  for (int i = 0; i < sizeof(checks) / sizeof(checks[0]); ++i) {
    compare(&a, &b, checks[i].arr_a, checks[i].arr_b, checks[i].range, checks[i].range_ext, checks[i].kind,
      checks[i].nblk, &err, &ref);
    TEST_CHECK(ref > 0.0);
    TEST_MSG("cdim=%d %s: reference vanishes", cdim, checks[i].name);
    TEST_CHECK(err < tol * ref);
    TEST_MSG("cdim=%d %s: max error %g (max reference %g)", cdim, checks[i].name, err, ref);
  }

  teardown(&b);
  teardown(&a);
}

// Sum of the cell-averaged right-hand side of component comp over the interior cells
// against the sum of its magnitudes.
static void
sum_rhs(const struct gkyl_array *rhs, const struct gkyl_range *range, const struct gkyl_range *range_ext,
  int comp, double *sum, double *mag)
{
  *sum = 0.0;
  *mag = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    const double *r = gkyl_array_cfetch(rhs, gkyl_range_idx(range_ext, iter.idx));
    *sum += r[comp];
    *mag += fabs(r[comp]);
  }
}

// Set every component of a configuration-space array to a constant on the extended range.
static void
set_uniform(struct sim *s, struct gkyl_array *arr, int ncomp, const double *vals)
{
  int nc = s->cbasis.num_basis;
  double c0 = pow(sqrt(2.0), s->cdim);
  gkyl_array_clear(arr, 0.0);
  for (long i = 0; i < arr->size; ++i) {
    double *d = gkyl_array_fetch(arr, i);
    for (int c = 0; c < ncomp; ++c) {
      d[c * nc] = c0 * vals[c];
    }
  }
}

static void
test_sine(int cdim)
{
  struct sim s;
  setup(&s, cdim, MAP_SINE);
  // F has volume sources proportional to the gradients of ln n, u, b and xi (the
  // normalization by n and xi), so for the conservation check only the streaming and
  // the face fluxes (which do depend on b, u and xi) are kept.
  gkyl_array_clear(s.grad_b, 0.0);
  gkyl_array_clear(s.grad_u, 0.0);
  gkyl_array_clear(s.grad_xi, 0.0);
  gkyl_array_clear(s.grad_ln_n, 0.0);
  gkyl_array_clear(s.theta, 0.0);
  advance(&s);
  int nc = s.cbasis.num_basis;
  // J varies from cell to cell.
  double jmin = 1e300, jmax = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.clocal);
  while (gkyl_range_iter_next(&iter)) {
    double J = jac_tot(&s, gkyl_range_idx(&s.clocal_ext, iter.idx));
    jmin = J < jmin ? J : jmin;
    jmax = J > jmax ? J : jmax;
  }
  TEST_CHECK(jmax > 1.5 * jmin);
  TEST_MSG("cdim=%d: J in [%g, %g]", cdim, jmin, jmax);

  // The fluid flux is single-valued on a face, so the fluid state is conserved.
  double sum, mag;
  for (int c = 0; c < 5; ++c) {
    sum_rhs(s.fl_rhs, &s.clocal, &s.clocal_ext, c * nc, &sum, &mag);
    TEST_CHECK(mag > 0.0);
    TEST_CHECK(fabs(sum) < 1e-12 * mag);
    TEST_MSG("cdim=%d fluid component %d: sum %g, |sum| %g", cdim, c, sum, mag);
  }
  teardown(&s);

  // The kinetic flux is the flux each side receives in its own frame (rho, u, xi), so the sum of
  // the normalized F is conserved only where the two sides of every face share a frame: a uniform
  // frame, with F still varying in space.
  setup(&s, cdim, MAP_SINE);
  gkyl_array_clear(s.grad_b, 0.0);
  gkyl_array_clear(s.grad_u, 0.0);
  gkyl_array_clear(s.grad_xi, 0.0);
  gkyl_array_clear(s.grad_ln_n, 0.0);
  gkyl_array_clear(s.theta, 0.0);
  const double b0[3] = {0.6, 0.64, 0.48}, u0[3] = {0.3, -0.2, 0.1};
  const double xi0[3] = {1.1, 1.1 * 1.1, 0.9 * 1.3 * 1.1 * 1.1};
  const double fluid0[5] = {1.3, 1.3 * u0[0], 1.3 * u0[1], 1.3 * u0[2], 2.0};
  set_uniform(&s, s.bvar, 3, b0);
  set_uniform(&s, s.u, 3, u0);
  set_uniform(&s, s.xi, 3, xi0);
  set_uniform(&s, s.fluid_phys, 5, fluid0);
  advance(&s);
  sum_rhs(s.kin_rhs, &s.plocal, &s.plocal_ext, 0, &sum, &mag);
  TEST_CHECK(mag > 0.0);
  TEST_CHECK(fabs(sum) < 1e-12 * mag);
  TEST_MSG("cdim=%d F: sum %g, |sum| %g", cdim, sum, mag);
  teardown(&s);
}

void test_stretch_1x() { test_stretch(1); }
void test_stretch_2x() { test_stretch(2); }
void test_stretch_3x() { test_stretch(3); }
void test_sine_1x() { test_sine(1); }
void test_sine_2x() { test_sine(2); }
void test_sine_3x() { test_sine(3); }

TEST_LIST = {
  {"test_stretch_1x", test_stretch_1x},
  {"test_stretch_2x", test_stretch_2x},
  {"test_stretch_3x", test_stretch_3x},
  {"test_sine_1x", test_sine_1x},
  {"test_sine_2x", test_sine_2x},
  {"test_sine_3x", test_sine_3x},
  {NULL, NULL},
};
