#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_flash_proj_pkpm_lte_on_basis.h>
#include <gkyl_util.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct flash_ctx {
  double mass;
  double n;
  double ux, uy, uz;
  double temp;

  double L;
  double vmax;
  int NC;
};

struct flash_ctx
init_ctx()
{
  double mass = 1.0;
  double n = 1.0;
  double ux = 0.5;
  double uy = 0.0;
  double uz = 0.0;
  double temp = 3.0;
  double vth = sqrt(temp/mass);
  double vmax = 5.0*vth;

  double L = 1.0;

  return (struct flash_ctx) {
    .mass = mass,
    .n = n,
    .ux = ux,
    .uy = uy,
    .uz = uz,
    .temp = temp,

    .L = L,
    .vmax = vmax,
    .NC = 16,
  };
}

void init_fluid(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct flash_ctx *app = ctx;
  fout[0] = app->n;
  fout[1] = app->ux;
  fout[2] = app->uy;
  fout[3] = app->uz;
  fout[4] = app->temp/app->mass;
}

// Indices of the phase-space cells checked against reference values, all in
// the first configuration-space cell:
// 1. Lower edge in velocity space
// 2. Interior in velocity space
#define NUM_CHECK_CELLS 2
// Maximum number of checked coefficients per cell: [F0, G] on the 1x1v hybrid basis.
#define MAX_CHECK_COEFFS 12
static void
check_cells(int cdim, int vdim, int num_cells_v, int idx[NUM_CHECK_CELLS][GKYL_MAX_DIM])
{
  for (int d=0; d<cdim; ++d) {
    idx[0][d] = 1;
    idx[1][d] = 1;
  }
  for (int d=cdim; d<cdim+vdim; ++d) {
    idx[0][d] = 1;
    idx[1][d] = num_cells_v/2;
  }
}

void
test_proj_lte(int cdim, int vdim, int NV, const double f_ref[NUM_CHECK_CELLS][MAX_CHECK_COEFFS])
{
  struct flash_ctx ctx = init_ctx();
  int poly_order = 1; // Only ever using hybrid phase basis.

  int pdim = cdim+vdim;

  int cells[] = {ctx.NC, NV, NV, NV};
  int ghost[] = {0, 0, 0, 0};

  double lower[] = {0.0, -ctx.vmax, -ctx.vmax, -ctx.vmax};
  double upper[] = {ctx.L, ctx.vmax, ctx.vmax, ctx.vmax};

  struct gkyl_rect_grid conf_grid;
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  struct gkyl_rect_grid phase_grid;
  struct gkyl_range phase_range, phase_range_ext;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phase_grid, ghost, &phase_range_ext, &phase_range);

  // initialize basis
  struct gkyl_basis phase_basis, conf_basis;
  gkyl_cart_modal_hybrid(&phase_basis, cdim, vdim);
  gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);

  // Initialize arrays
  struct gkyl_array *fluidin, *distf;
  fluidin = mkarr(5*conf_basis.num_basis, conf_range_ext.volume);
  distf = mkarr(2*phase_basis.num_basis, phase_range_ext.volume);

  // Project initial fluid moments
  int num_quad = poly_order+1;
  gkyl_proj_on_basis *proj_fluid = gkyl_proj_on_basis_new(&conf_grid,
    &conf_basis,num_quad, 5, init_fluid, &ctx);
  
  gkyl_proj_on_basis_advance(proj_fluid, 0.0, &conf_range, fluidin);
  gkyl_proj_on_basis_release(proj_fluid);

  // Initialize updater
  struct gkyl_flash_proj_pkpm_lte_on_basis_inp inp = {
    .phase_grid = &phase_grid,
    .conf_basis = &conf_basis,
    .phase_basis = &phase_basis,
    .use_gpu = false,
    .quad_type = 0
  };
  gkyl_flash_proj_pkpm_lte_on_basis *proj = gkyl_flash_proj_pkpm_lte_on_basis_new(&inp);

  // Project distribution function
  gkyl_flash_proj_pkpm_lte_on_basis_advance(proj, &phase_range, &conf_range, fluidin, distf);

  // Check the [F0, G] coefficients in each representative cell.
  int num_coeffs = 2*phase_basis.num_basis;
  TEST_ASSERT( num_coeffs <= MAX_CHECK_COEFFS );

  int idx[NUM_CHECK_CELLS][GKYL_MAX_DIM];
  check_cells(cdim, vdim, NV, idx);
  for (int i=0; i<NUM_CHECK_CELLS; ++i) {
    const double *f_d = gkyl_array_cfetch(distf, gkyl_range_idx(&phase_range, idx[i]));
    for (int k=0; k<num_coeffs; ++k) {
      TEST_CHECK( gkyl_compare_double(f_ref[i][k], f_d[k], 1e-12) );
      TEST_MSG("Cell %d, coefficient %d: expected %.16e, got %.16e", i, k, f_ref[i][k], f_d[k]);
    }
  }

  // Write input moments and distribution functions.
  struct gkyl_msgpack_map_elem io_meta_conf[] = {
    {.key = "polyOrder", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = conf_basis.poly_order},
    {.key = "basisType", .elem_type = GKYL_MP_STRING, .cval = conf_basis.id}
  };
  struct gkyl_msgpack_data *mt_conf = gkyl_msgpack_create(
    sizeof(io_meta_conf)/sizeof(io_meta_conf[0]), io_meta_conf);

  struct gkyl_msgpack_map_elem io_meta_phase[] = {
    {.key = "polyOrder", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = phase_basis.poly_order},
    {.key = "basisType", .elem_type = GKYL_MP_STRING, .cval = phase_basis.id},
    {.key = "num_cdim", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = cdim},
    {.key = "num_vdim", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = vdim}
  };
  struct gkyl_msgpack_data *mt_phase = gkyl_msgpack_create(
    sizeof(io_meta_phase)/sizeof(io_meta_phase[0]), io_meta_phase);

  char fname[256];
  snprintf(fname, sizeof(fname), "ctest_flash_proj_pkpm_lte_moms_%dx_p%d_%d.gkyl", cdim, poly_order, NV);
  gkyl_grid_sub_array_write(&conf_grid, &conf_range, mt_conf, fluidin, fname);
  snprintf(fname, sizeof(fname), "ctest_flash_proj_pkpm_lte_%dx%dv_p%d_%d.gkyl", cdim, vdim, poly_order, NV);
  gkyl_grid_sub_array_write(&phase_grid, &phase_range, mt_phase, distf, fname);

  gkyl_msgpack_data_release(mt_conf);
  gkyl_msgpack_data_release(mt_phase);

  gkyl_flash_proj_pkpm_lte_on_basis_release(proj);
  gkyl_array_release(fluidin);
  gkyl_array_release(distf);
}

// NOTE: Cell number is for velocity space. Reference values are the [F0, G]
// coefficients in the cells given by check_cells.

void test_1x1v_p1_4()
{
  const double ref[NUM_CHECK_CELLS][MAX_CHECK_COEFFS] = {
    {  2.8679956368919002e-03,  0.0000000000000000e+00,  3.6855348348247478e-03,  0.0000000000000000e+00,  2.6353058158845046e-03,  0.0000000000000000e+00,
       8.6039869106757007e-03,  0.0000000000000000e+00,  1.1056604504474246e-02,  0.0000000000000000e+00,  7.9059174476535168e-03,  0.0000000000000000e+00 },
    {  2.2807140026223627e-01,  0.0000000000000000e+00,  1.5092914984723757e-01,  0.0000000000000000e+00,  1.0723858747704259e-02,  0.0000000000000000e+00,
       6.8421420078670903e-01,  0.0000000000000000e+00,  4.5278744954171274e-01,  0.0000000000000000e+00,  3.2171576243112780e-02,  0.0000000000000000e+00 },
  };
  test_proj_lte(1, 1, 4, ref);
}

void test_1x1v_p1_8()
{
  const double ref[NUM_CHECK_CELLS][MAX_CHECK_COEFFS] = {
    {  8.1411606287512921e-05,  0.0000000000000000e+00,  8.8224205287364664e-05,  0.0000000000000000e+00,  4.8601910087969888e-05,  0.0000000000000000e+00,
       2.4423481886253872e-04,  0.0000000000000000e+00,  2.6467261586209410e-04,  0.0000000000000000e+00,  1.4580573026390969e-04,  0.0000000000000000e+00 },
    {  3.6428513436224480e-01,  0.0000000000000000e+00,  7.7249771675131917e-02,  0.0000000000000000e+00, -1.2946790186796578e-02,  0.0000000000000000e+00,
       1.0928554030867350e+00,  0.0000000000000000e+00,  2.3174931502539584e-01,  0.0000000000000000e+00, -3.8840370560389685e-02,  0.0000000000000000e+00 },
  };
  test_proj_lte(1, 1, 8, ref);
}

void test_1x1v_p1_16()
{
  const double ref[NUM_CHECK_CELLS][MAX_CHECK_COEFFS] = {
    {  1.0687857082806368e-05,  0.0000000000000000e+00,  7.8821695698877775e-06,  0.0000000000000000e+00,  2.6711213362193507e-06,  0.0000000000000000e+00,
       3.2063571248419113e-05,  0.0000000000000000e+00,  2.3646508709663336e-05,  0.0000000000000000e+00,  8.0133640086580541e-06,  0.0000000000000000e+00 },
    {  4.3234661695439447e-01,  0.0000000000000000e+00,  2.4045673667392367e-02,  0.0000000000000000e+00, -5.6561964531881923e-03,  0.0000000000000000e+00,
       1.2970398508631835e+00,  0.0000000000000000e+00,  7.2137021002177129e-02,  0.0000000000000000e+00, -1.6968589359564686e-02,  0.0000000000000000e+00 },
  };
  test_proj_lte(1, 1, 16, ref);
}

TEST_LIST = {
  { "test_1x1v_p1_4", test_1x1v_p1_4 },
  { "test_1x1v_p1_8", test_1x1v_p1_8 },
  { "test_1x1v_p1_16", test_1x1v_p1_16 },
  { NULL, NULL }
};
