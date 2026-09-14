// Test the dg_lowpass_filter updater.
#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_lowpass_filter.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

#include <math.h>

#define FILTER_DIR 0 // Direction filtered in every test below.
#define FOUT_FILL 7.0 // Value fout is cleared to before each application.

static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a = on_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static double
filter_gain(int M, double fc, double freq)
{
  // Frequency response of the kernel at freq cycles/cell.
  double wsum = 0.0, gain = 0.0;
  for (int k=-M; k<M+1; k++) {
    double hk = k == 0? 2.0*fc : sin(2.0*M_PI*fc*k)/(M_PI*k);
    double wk = 0.42 + 0.5*cos(M_PI*k/M) + 0.08*cos(2.0*M_PI*k/M);
    wsum += hk*wk;
    gain += hk*wk*cos(2.0*M_PI*freq*k);
  }
  return gain/wsum;
}

struct profile_ctx {
  int ndim; // Dimensionality of the field.
  double mode_num; // Mode number along the filtered direction, domain is [0,1].
  double x0, w; // Center and width of the bump along the filtered direction.
};

void eval_const(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = 3.0;
}

void eval_linear(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = 1.5 + 0.5*xn[FILTER_DIR];
}

static double
transverse_mod(const double *xn, int ndim)
{
  // Smooth modulation across the directions that are not filtered.
  double mod = 1.0;
  for (int d=0; d<ndim; d++) {
    if (d != FILTER_DIR)
      mod *= 1.0 + 0.3*cos(2.0*M_PI*xn[d]);
  }
  return mod;
}

void eval_mode(double t, const double *xn, double *fout, void *ctx)
{
  struct profile_ctx *pctx = ctx;
  fout[0] = cos(2.0*M_PI*pctx->mode_num*xn[FILTER_DIR])*transverse_mod(xn, pctx->ndim);
}

void eval_bump(double t, const double *xn, double *fout, void *ctx)
{
  // Off-center bump along the filtered direction.
  struct profile_ctx *pctx = ctx;
  fout[0] = (0.1 + exp(-pow((xn[FILTER_DIR]-pctx->x0)/pctx->w, 2)))*transverse_mod(xn, pctx->ndim);
}

struct filter_env {
  bool use_gpu;
  struct gkyl_rect_grid grid;
  struct gkyl_basis basis;
  struct gkyl_range local, local_ext;
  struct gkyl_array *fin, *fout; // What the updater sees.
  struct gkyl_array *fin_ho, *fout_ho; // Host copies the checks read.
};

static void
filter_env_new(struct filter_env *env, bool use_gpu, int ndim, const int *cells, int poly_order)
{
  // Grid, basis, ranges and arrays for one test, on the unit cube.
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int nghost[GKYL_MAX_DIM];
  for (int d=0; d<ndim; d++) {
    lower[d] = 0.0;
    upper[d] = 1.0;
    nghost[d] = 1;
  }

  env->use_gpu = use_gpu;
  gkyl_rect_grid_init(&env->grid, ndim, lower, upper, cells);
  gkyl_cart_modal_serendip(&env->basis, ndim, poly_order);
  gkyl_create_grid_ranges(&env->grid, nghost, &env->local_ext, &env->local);

  env->fin = mkarr(use_gpu, env->basis.num_basis, env->local_ext.volume);
  env->fout = mkarr(use_gpu, env->basis.num_basis, env->local_ext.volume);
  env->fin_ho = use_gpu? mkarr(false, env->fin->ncomp, env->fin->size) : gkyl_array_acquire(env->fin);
  env->fout_ho = use_gpu? mkarr(false, env->fout->ncomp, env->fout->size) : gkyl_array_acquire(env->fout);
}

static void
filter_env_release(struct filter_env *env)
{
  gkyl_array_release(env->fin);
  gkyl_array_release(env->fout);
  gkyl_array_release(env->fin_ho);
  gkyl_array_release(env->fout_ho);
}

static void
filter_apply(struct filter_env *env, const struct gkyl_range *sub, int M,
  double cutoff_wavelength, evalf_t func, void *func_ctx)
{
  // Project func and filter it over sub, leaving both fields on the host.
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&env->grid, &env->basis,
    env->basis.poly_order+1, 1, func, func_ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &env->local, env->fin_ho);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(env->fin, env->fin_ho);

  gkyl_array_clear(env->fout, FOUT_FILL);
  struct gkyl_dg_lowpass_filter *lpf = gkyl_dg_lowpass_filter_new(FILTER_DIR, M,
    cutoff_wavelength, &env->basis, &env->grid, sub, env->use_gpu);
  gkyl_dg_lowpass_filter_advance(lpf, env->fin, env->fout);
  gkyl_dg_lowpass_filter_release(lpf);
  gkyl_array_copy(env->fout_ho, env->fout);
}

static double
filter_integral_check(struct filter_env *env, const struct gkyl_range *sub,
  double tol, const char *what)
{
  // Change of the integral over sub, as a fraction of the mass in the range.
  double tot_in = 0.0, tot_out = 0.0, mass = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, sub);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(sub, iter.idx);
    double in = ((const double *) gkyl_array_cfetch(env->fin_ho, linidx))[0];
    tot_in += in;
    mass += fabs(in);
    tot_out += ((const double *) gkyl_array_cfetch(env->fout_ho, linidx))[0];
  }

  double rel = (tot_out-tot_in)/mass;
  TEST_CHECK( fabs(rel) < tol );
  TEST_MSG("%s: the integral changed by %.3e of the mass in the range", what, rel);
  return rel;
}

static void
filter_gain_check(struct filter_env *env, const struct gkyl_range *sub, int edge_skip,
  double gain, double tol, const char *what)
{
  // Check fout = gain*fin, skipping edge_skip cells at each edge of sub.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, sub);
  while (gkyl_range_iter_next(&iter)) {
    if (iter.idx[FILTER_DIR] < sub->lower[FILTER_DIR]+edge_skip ||
        iter.idx[FILTER_DIR] > sub->upper[FILTER_DIR]-edge_skip)
      continue;

    long linidx = gkyl_range_idx(sub, iter.idx);
    const double *in_c = gkyl_array_cfetch(env->fin_ho, linidx);
    const double *out_c = gkyl_array_cfetch(env->fout_ho, linidx);
    for (int c=0; c<env->basis.num_basis; c++) {
      TEST_CHECK( fabs(out_c[c] - gain*in_c[c]) < tol );
      TEST_MSG("%s, cell %d coeff %d: expected %.13e, got %.13e",
        what, iter.idx[FILTER_DIR], c, gain*in_c[c], out_c[c]);
    }
  }
}

static void
test_response(bool use_gpu, int ndim, const int *cells)
{
  // Away from the edges the filter scales each field by the kernel response.
  int M = 8;
  double fc = 0.3;

  struct filter_env env;
  filter_env_new(&env, use_gpu, ndim, cells, 1);
  double cutoff = env.grid.dx[FILTER_DIR]/fc;

  double mod_max = 1.0; // Peak of the transverse modulation.
  for (int d=0; d<ndim; d++)
    mod_max *= d == FILTER_DIR? 1.0 : 1.3;

  // A constant passes through everywhere, and fout's ghost cells are untouched.
  filter_apply(&env, &env.local, M, cutoff, eval_const, NULL);
  filter_gain_check(&env, &env.local, 0, 1.0, 1e-12, "a constant");

  int idx_ghost[GKYL_MAX_DIM];
  for (int d=0; d<ndim; d++)
    idx_ghost[d] = env.local.lower[d];
  idx_ghost[FILTER_DIR] -= 1;
  const double *out_g = gkyl_array_cfetch(env.fout_ho, gkyl_range_idx(&env.local_ext, idx_ghost));
  for (int c=0; c<env.basis.num_basis; c++)
    TEST_CHECK( out_g[c] == FOUT_FILL );

  // A ramp passes through in the interior and opens no jump at any interface.
  filter_apply(&env, &env.local, M, cutoff, eval_linear, NULL);
  filter_gain_check(&env, &env.local, M, 1.0, 1e-12, "a ramp");

  double eta_lo[GKYL_MAX_DIM] = {0.0}, eta_up[GKYL_MAX_DIM] = {0.0};
  eta_lo[FILTER_DIR] = -1.0;
  eta_up[FILTER_DIR] =  1.0;
  double b_lo[32], b_up[32];
  env.basis.eval(eta_lo, b_lo);
  env.basis.eval(eta_up, b_up);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &env.local);
  while (gkyl_range_iter_next(&iter)) {
    if (iter.idx[FILTER_DIR] == env.local.upper[FILTER_DIR])
      continue;

    int idx_up[GKYL_MAX_DIM];
    for (int d=0; d<ndim; d++)
      idx_up[d] = iter.idx[d];
    idx_up[FILTER_DIR] += 1;

    const double *left = gkyl_array_cfetch(env.fout_ho, gkyl_range_idx(&env.local, iter.idx));
    const double *right = gkyl_array_cfetch(env.fout_ho, gkyl_range_idx(&env.local, idx_up));
    double jump = 0.0;
    for (int c=0; c<env.basis.num_basis; c++)
      jump += right[c]*b_lo[c] - left[c]*b_up[c];
    TEST_CHECK( fabs(jump) < 1.0e-12 );
    TEST_MSG("interface %d: filtering a ramp opened a jump of %.3e", iter.idx[FILTER_DIR], jump);
  }

  // A single Fourier mode is scaled by the response, the Nyquist one killed.
  double mode_nums[] = {2.0, 8.0, cells[FILTER_DIR]/2.0};
  for (int im=0; im<3; im++) {
    struct profile_ctx pctx = { .ndim = ndim, .mode_num = mode_nums[im] };
    double gain = filter_gain(M, fc, pctx.mode_num/cells[FILTER_DIR]);

    filter_apply(&env, &env.local, M, cutoff, eval_mode, &pctx);
    filter_gain_check(&env, &env.local, M, gain, 1e-12, "a Fourier mode");

    if (im == 2) {
      TEST_CHECK( fabs(gain) < 1e-3 );
      filter_gain_check(&env, &env.local, M, 0.0, 1e-3*mod_max, "the Nyquist mode");
    }
  }

  filter_env_release(&env);
}

static void
test_conservation(bool use_gpu, int ndim, const int *cells)
{
  // The integral survives a reflected stencil, and a truncated one only for a constant.
  int M = 8;

  struct filter_env env;
  filter_env_new(&env, use_gpu, ndim, cells, 1);
  double cutoff = env.grid.dx[FILTER_DIR]/0.25;

  // Whole grid: the stencil is reflected at both ends.
  // Sub-range ending inside the grid: the stencil is truncated at the upper end,
  // as the twist-shift's is at the LCFS.
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  for (int d=0; d<ndim; d++) {
    lo[d] = env.local.lower[d];
    up[d] = env.local.upper[d];
  }
  up[FILTER_DIR] = 3*cells[FILTER_DIR]/4;
  struct gkyl_range cut;
  gkyl_sub_range_init(&cut, &env.local, lo, up);

  struct gkyl_range *ranges[] = {&env.local, &cut};
  const char *const_names[] = {"reflected stencil, constant", "truncated stencil, constant"};
  // The truncated stencil with non constant field does not conserve the particle number.
  const char *bump_names[] = {"reflected stencil, bump"};//, "truncated stencil, bump"};

  // A bump next to the truncated edge, so the profile overlaps the cut rows.
  struct profile_ctx bump = { .ndim = ndim, .x0 = 0.68, .w = 0.05 };

  for (int r=0; r<2; r++) {
    filter_apply(&env, ranges[r], M, cutoff, eval_const, NULL);
    filter_integral_check(&env, ranges[r], 1e-12, const_names[r]);
  }
  filter_apply(&env, ranges[0], M, cutoff, eval_bump, &bump);
  filter_integral_check(&env, ranges[0], 1e-12, bump_names[0]);

  filter_env_release(&env);
}

static void test_1x_response(bool use_gpu) { test_response(use_gpu, 1, (int[]) {32}); }
static void test_2x_response(bool use_gpu) { test_response(use_gpu, 2, (int[]) {32, 8}); }
static void test_3x_response(bool use_gpu) { test_response(use_gpu, 3, (int[]) {32, 8, 6}); }

static void test_1x_conservation(bool use_gpu) { test_conservation(use_gpu, 1, (int[]) {64}); }
static void test_2x_conservation(bool use_gpu) { test_conservation(use_gpu, 2, (int[]) {64, 8}); }
static void test_3x_conservation(bool use_gpu) { test_conservation(use_gpu, 3, (int[]) {64, 8, 6}); }

void test_1x_response_ho() { test_1x_response(false); }
void test_2x_response_ho() { test_2x_response(false); }
void test_3x_response_ho() { test_3x_response(false); }
void test_1x_conservation_ho() { test_1x_conservation(false); }
void test_2x_conservation_ho() { test_2x_conservation(false); }
void test_3x_conservation_ho() { test_3x_conservation(false); }

#ifdef GKYL_HAVE_CUDA
void test_1x_response_cu() { test_1x_response(true); }
void test_2x_response_cu() { test_2x_response(true); }
void test_3x_response_cu() { test_3x_response(true); }
void test_1x_conservation_cu() { test_1x_conservation(true); }
void test_2x_conservation_cu() { test_2x_conservation(true); }
void test_3x_conservation_cu() { test_3x_conservation(true); }
#endif

TEST_LIST = {
  { "test_1x_response", test_1x_response_ho },
  { "test_2x_response", test_2x_response_ho },
  { "test_3x_response", test_3x_response_ho },
  { "test_1x_conservation", test_1x_conservation_ho },
  { "test_2x_conservation", test_2x_conservation_ho },
  { "test_3x_conservation", test_3x_conservation_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_1x_response_cu", test_1x_response_cu },
  { "test_2x_response_cu", test_2x_response_cu },
  { "test_3x_response_cu", test_3x_response_cu },
  { "test_1x_conservation_cu", test_1x_conservation_cu },
  { "test_2x_conservation_cu", test_2x_conservation_cu },
  { "test_3x_conservation_cu", test_3x_conservation_cu },
#endif
  { NULL, NULL },
};
