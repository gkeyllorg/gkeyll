// Tests of the momentum diffusion of the PKPM fluid on the position map (dg_pkpm_mom_diffusion):
// on the identity map it is the operator of gkyl_dg_updater_diffusion_fluid; on a sinusoidal map
// it conserves the momentum and converges to the physical d^2/dx^2 and -d^4/dx^4.
#include <math.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_pkpm_mom_diffusion.h>
#include <gkyl_dg_updater_diffusion_fluid.h>
#include <gkyl_evalf_def.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>

static struct gkyl_array *
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// x = xi + 0.4 sin(xi) on [-pi, pi]: fixed end points, periodic Jacobian 1 + 0.4 cos(xi).
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

// A generic smooth periodic fluid state [rho, M_x, M_y, M_z, E] (physical).
static void
eval_fluid(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.4, z = cdim > 2 ? xn[2] : -0.3;
  fout[0] = 1.0 + 0.3 * cos(x);
  fout[1] = sin(x) + 0.3 * cos(2.0 * y) + 0.2 * sin(z + 0.5);
  fout[2] = cos(x) * sin(y) + 0.1 * cos(3.0 * x);
  fout[3] = 0.5 * sin(2.0 * x + y) - 0.2 * cos(z);
  fout[4] = 2.0 + 0.1 * sin(y);
}

// A separable momentum M = sin(x) + cos(y) + sin(z) in every component, and its exact diffusion.
static void
eval_sep(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double m = sin(xn[0]) + (cdim > 1 ? cos(xn[1]) : 0.0) + (cdim > 2 ? sin(xn[2]) : 0.0);
  fout[0] = 1.0;
  fout[1] = fout[2] = fout[3] = m;
  fout[4] = 1.0;
}

// d^2 M = -M and d^4 M = M per direction: order 2 gives D sum_d d_d^2 M = -D M, order 4 gives
// -D sum_d d_d^4 M = -D M.
static void
eval_sep_exact(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double m = sin(xn[0]) + (cdim > 1 ? cos(xn[1]) : 0.0) + (cdim > 2 ? sin(xn[2]) : 0.0);
  fout[0] = fout[4] = 0.0;
  fout[1] = fout[2] = fout[3] = -m;
}

struct sim {
  int cdim;
  struct gkyl_rect_grid grid;
  struct gkyl_basis basis;
  struct gkyl_range local, local_ext;
  struct gkyl_vlasov_position_map *pos_map;
  struct gkyl_array *fluid, *rhs, *cfl;
};

// Periodic ghost cells of a configuration array with one ghost layer in every direction.
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

// The physical field projected on the (mapped) grid, J-weighted, with periodic ghost cells.
static void
project_weighted(struct sim *s, evalf_t eval, struct gkyl_array *out)
{
  struct gkyl_array *phys = mkarr(out->ncomp, out->size);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &s->grid,
    .basis = &s->basis,
    .num_quad = 8,
    .num_ret_vals = 5,
    .eval = eval,
    .ctx = &s->cdim,
    .c2p_func = pos_c2p,
    .c2p_func_ctx = s->pos_map,
  });
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local, phys);
  gkyl_proj_on_basis_release(proj);
  gkyl_vlasov_position_map_rescale_jacobpos_conf(s->pos_map, &s->local, phys, out);
  apply_periodic(out, &s->local_ext, s->cdim);
  gkyl_array_release(phys);
}

static void
setup(struct sim *s, int cdim, int ncell, bool mapped)
{
  s->cdim = cdim;
  double lower[3] = {-M_PI, -M_PI, -M_PI}, upper[3] = {M_PI, M_PI, M_PI};
  int cells[3] = {ncell, ncell, ncell}, ghost[3] = {1, 1, 1};
  gkyl_rect_grid_init(&s->grid, cdim, lower, upper, cells);
  gkyl_cart_modal_serendip(&s->basis, cdim, 1);
  gkyl_create_grid_ranges(&s->grid, ghost, &s->local_ext, &s->local);
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = {0};
  for (int d = 0; d < cdim; ++d) {
    if (mapped) {
      inp_pmap[d].eval_pmap = pmap_sine;
    }
  }
  s->pos_map = gkyl_vlasov_position_map_new(&s->grid, &s->local, &s->local_ext, &s->basis, inp_pmap, false);
  int nc = s->basis.num_basis;
  s->fluid = mkarr(5 * nc, s->local_ext.volume);
  s->rhs = mkarr(5 * nc, s->local_ext.volume);
  s->cfl = mkarr(1, s->local_ext.volume);
}

static void
teardown(struct sim *s)
{
  gkyl_array_release(s->fluid);
  gkyl_array_release(s->rhs);
  gkyl_array_release(s->cfl);
  gkyl_vlasov_position_map_release(s->pos_map);
}

static struct gkyl_dg_pkpm_mom_diffusion *
new_diff(struct sim *s, int order, double D)
{
  struct gkyl_dg_pkpm_mom_diffusion_inp inp = {
    .conf_grid = &s->grid, .conf_basis = &s->basis, .pos_map = s->pos_map, .order = order, .use_gpu = false};
  for (int d = 0; d < s->cdim; ++d) {
    inp.D[d] = D;
  }
  return gkyl_dg_pkpm_mom_diffusion_inew(&inp);
}

// Identity map: the operator and the rate of gkyl_dg_updater_diffusion_fluid on the momentum.
static void
test_uniform(int cdim, int order)
{
  struct sim s;
  setup(&s, cdim, 8, false);
  int nc = s.basis.num_basis;
  double D = 0.37;
  project_weighted(&s, eval_fluid, s.fluid);
  gkyl_array_clear(s.rhs, 0.0);
  gkyl_array_clear(s.cfl, 0.0);
  struct gkyl_dg_pkpm_mom_diffusion *up = new_diff(&s, order, D);
  gkyl_dg_pkpm_mom_diffusion_advance(up, &s.local, s.fluid, s.cfl, s.rhs);

  struct gkyl_array *mom = mkarr(3 * nc, s.local_ext.volume), *rhs3 = mkarr(3 * nc, s.local_ext.volume);
  struct gkyl_array *cfl3 = mkarr(1, s.local_ext.volume), *diffD = mkarr(cdim, 1);
  gkyl_array_set_offset(mom, 1.0, s.fluid, nc);
  gkyl_array_clear(rhs3, 0.0);
  gkyl_array_clear(cfl3, 0.0);
  gkyl_array_clear(diffD, 0.0);
  for (int d = 0; d < cdim; ++d) {
    gkyl_array_shiftc(diffD, D, d);
  }
  bool is_zero_flux[GKYL_MAX_CDIM] = {false};
  struct gkyl_dg_updater_diffusion_fluid *ref =
    gkyl_dg_updater_diffusion_fluid_new(&s.grid, &s.basis, 3, NULL, order, is_zero_flux, false);
  gkyl_dg_updater_diffusion_fluid_advance(ref, &s.local, diffD, mom, cfl3, rhs3);

  double max_err = 0.0, max_ref = 0.0, max_cfl_err = 0.0;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &s.local);
  while (gkyl_range_iter_next(&it)) {
    long loc = gkyl_range_idx(&s.local, it.idx);
    const double *r = gkyl_array_cfetch(s.rhs, loc), *r3 = gkyl_array_cfetch(rhs3, loc);
    for (int k = 0; k < 3 * nc; ++k) {
      max_err = fmax(max_err, fabs(r[nc + k] - r3[k]));
      max_ref = fmax(max_ref, fabs(r3[k]));
    }
    TEST_CHECK(r[0] == 0.0); // rho and E untouched
    max_cfl_err = fmax(max_cfl_err, fabs(((const double *)gkyl_array_cfetch(s.cfl, loc))[0] -
      ((const double *)gkyl_array_cfetch(cfl3, loc))[0]));
  }
  TEST_CHECK(max_ref > 0.0);
  TEST_CHECK(max_err < 1e-12 * max_ref);
  TEST_MSG("cdim=%d order=%d: against dg_diffusion_fluid, max error %g (max %g)", cdim, order, max_err, max_ref);
  TEST_CHECK(max_cfl_err < 1e-12 * ((const double *)gkyl_array_cfetch(cfl3, gkyl_range_idx(&s.local, s.local.lower)))[0]);
  TEST_MSG("cdim=%d order=%d: rate error %g", cdim, order, max_cfl_err);

  gkyl_dg_updater_diffusion_fluid_release(ref);
  gkyl_array_release(mom);
  gkyl_array_release(rhs3);
  gkyl_array_release(cfl3);
  gkyl_array_release(diffD);
  gkyl_dg_pkpm_mom_diffusion_release(up);
  teardown(&s);
}

// Sinusoidal map, periodic: the cell averages of the J-weighted momentum change by face fluxes only.
static void
test_conservation(int cdim, int order)
{
  struct sim s;
  setup(&s, cdim, 8, true);
  int nc = s.basis.num_basis;
  project_weighted(&s, eval_fluid, s.fluid);
  gkyl_array_clear(s.rhs, 0.0);
  gkyl_array_clear(s.cfl, 0.0);
  struct gkyl_dg_pkpm_mom_diffusion *up = new_diff(&s, order, 0.37);
  gkyl_dg_pkpm_mom_diffusion_advance(up, &s.local, s.fluid, s.cfl, s.rhs);
  double sum[3] = {0.0}, abs_sum = 0.0;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &s.local);
  while (gkyl_range_iter_next(&it)) {
    const double *r = gkyl_array_cfetch(s.rhs, gkyl_range_idx(&s.local, it.idx));
    for (int c = 0; c < 3; ++c) {
      sum[c] += r[(c + 1) * nc];
      abs_sum += fabs(r[(c + 1) * nc]);
    }
  }
  TEST_CHECK(abs_sum > 0.0);
  for (int c = 0; c < 3; ++c) {
    TEST_CHECK(fabs(sum[c]) < 1e-13 * abs_sum);
    TEST_MSG("cdim=%d order=%d: momentum component %d changes by %g (sum of |changes| %g)", cdim, order, c,
      sum[c], abs_sum);
  }
  gkyl_dg_pkpm_mom_diffusion_release(up);
  teardown(&s);
}

// Sinusoidal map: the cell averages of the physical result against the exact diffusion of a
// separable momentum, at two resolutions; returns the two errors.
static void
accuracy(int cdim, int order, int ncell, double *err)
{
  struct sim s;
  setup(&s, cdim, ncell, true);
  int nc = s.basis.num_basis;
  double D = 1.0;
  project_weighted(&s, eval_sep, s.fluid);
  struct gkyl_array *exact = mkarr(5 * nc, s.local_ext.volume), *rhs_phys = mkarr(5 * nc, s.local_ext.volume);
  project_weighted(&s, eval_sep_exact, exact);
  gkyl_array_clear(s.rhs, 0.0);
  gkyl_array_clear(s.cfl, 0.0);
  struct gkyl_dg_pkpm_mom_diffusion *up = new_diff(&s, order, D);
  gkyl_dg_pkpm_mom_diffusion_advance(up, &s.local, s.fluid, s.cfl, s.rhs);
  // both J-weighted: compare the cell averages of J (op M) with those of J (exact)
  double e = 0.0, m = 0.0;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &s.local);
  while (gkyl_range_iter_next(&it)) {
    long loc = gkyl_range_idx(&s.local, it.idx);
    const double *r = gkyl_array_cfetch(s.rhs, loc), *x = gkyl_array_cfetch(exact, loc);
    for (int c = 1; c < 4; ++c) {
      e = fmax(e, fabs(r[c * nc] - x[c * nc]));
      m = fmax(m, fabs(x[c * nc]));
    }
  }
  *err = e / m;
  gkyl_array_release(exact);
  gkyl_array_release(rhs_phys);
  gkyl_dg_pkpm_mom_diffusion_release(up);
  teardown(&s);
}

static void
test_accuracy(int cdim, int order)
{
  double e1, e2;
  int n1 = cdim == 1 ? 32 : 16;
  accuracy(cdim, order, n1, &e1);
  accuracy(cdim, order, 2 * n1, &e2);
  double rate = log(e1 / e2) / log(2.0);
  TEST_CHECK(e2 < e1);
  TEST_CHECK(rate > 1.8);
  TEST_MSG("cdim=%d order=%d: relative error of the cell averages %g (N = %d), %g (N = %d), order %g", cdim,
    order, e1, n1, e2, 2 * n1, rate);
}

void test_uniform_1x_o2() { test_uniform(1, 2); }
void test_uniform_1x_o4() { test_uniform(1, 4); }
void test_uniform_2x_o2() { test_uniform(2, 2); }
void test_uniform_2x_o4() { test_uniform(2, 4); }
void test_uniform_3x_o4() { test_uniform(3, 4); }
void test_conservation_1x_o4() { test_conservation(1, 4); }
void test_conservation_2x_o2() { test_conservation(2, 2); }
void test_conservation_2x_o4() { test_conservation(2, 4); }
void test_conservation_3x_o4() { test_conservation(3, 4); }
void test_accuracy_1x_o2() { test_accuracy(1, 2); }
void test_accuracy_1x_o4() { test_accuracy(1, 4); }
void test_accuracy_2x_o4() { test_accuracy(2, 4); }

TEST_LIST = {
  {"test_uniform_1x_o2", test_uniform_1x_o2},
  {"test_uniform_1x_o4", test_uniform_1x_o4},
  {"test_uniform_2x_o2", test_uniform_2x_o2},
  {"test_uniform_2x_o4", test_uniform_2x_o4},
  {"test_uniform_3x_o4", test_uniform_3x_o4},
  {"test_conservation_1x_o4", test_conservation_1x_o4},
  {"test_conservation_2x_o2", test_conservation_2x_o2},
  {"test_conservation_2x_o4", test_conservation_2x_o4},
  {"test_conservation_3x_o4", test_conservation_3x_o4},
  {"test_accuracy_1x_o2", test_accuracy_1x_o2},
  {"test_accuracy_1x_o4", test_accuracy_1x_o4},
  {"test_accuracy_2x_o4", test_accuracy_2x_o4},
  {NULL, NULL},
};
