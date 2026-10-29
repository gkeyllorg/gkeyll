#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_fv_interpolate.h>
#include <gkyl_fv_proj.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>

struct interp_ctx {
  double ndim;
};

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Initialize test function to interpolate.
void
proj_fn(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  struct interp_ctx *app = (struct interp_ctx *)ctx;

  if (app->ndim == 1) {
    fout[0] = sin(M_PI*x);
  }
  else if (app->ndim == 2) {
    fout[0] = sin(M_PI*x)*cos(M_PI*y);
  }
  else if (app->ndim == 3) {
    fout[0] = sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*z+M_PI/3.0);
  }
  else {
    fout[0] = 0.0;
  }
}

// Indices of the DG cells checked against reference values: 
// 1. Lower corner
// 2. Interior
// 3. upper corner
#define NUM_CHECK_CELLS 3
static void
check_cells(int ndim, int num_cells_dg, int idx[NUM_CHECK_CELLS][GKYL_MAX_DIM])
{
  for (int d=0; d<ndim; ++d) {
    idx[0][d] = 1;
    idx[1][d] = 2;
    idx[2][d] = num_cells_dg;
  }
}

void
test_interp_modal_to_fv(int ndim, int num_cells_fv, const double fv_ref[NUM_CHECK_CELLS][8])
{
  struct interp_ctx ctx = {.ndim = ndim};

  // Initialize grids.
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM] = {0.0};
  int cells_dg[GKYL_MAX_DIM], cells_fv[GKYL_MAX_DIM] = {0};
  int ghost[GKYL_MAX_DIM] = {0};

  double L = 1.0;
  for (int d=0; d<ndim; ++d) {
    lower[d] = -L;
    upper[d] = L;
    cells_fv[d] = num_cells_fv;
    cells_dg[d] = (int)(num_cells_fv/2);
  }

  struct gkyl_rect_grid grid_dg, grid_fv;
  gkyl_rect_grid_init(&grid_dg, ndim, lower, upper, cells_dg);
  gkyl_rect_grid_init(&grid_fv, ndim, lower, upper, cells_fv);

  // Construct grid ranges.
  struct gkyl_range range_dg, range_dg_ext, range_fv, range_fv_ext;
  gkyl_create_grid_ranges(&grid_dg, ghost, &range_dg_ext, &range_dg);
  gkyl_create_grid_ranges(&grid_fv, ghost, &range_fv_ext, &range_fv);

  // Initialize basis, always p=1.
  struct gkyl_basis basis;
  int poly_order = 1;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Initialize arrays for DG and FV data.
  struct gkyl_array *dg_data = mkarr(basis.num_basis, range_dg.volume);
  struct gkyl_array *fv_data = mkarr(1, range_fv.volume);
  struct gkyl_array *fv_data_exact = mkarr(1, range_fv.volume);

  // Project initial condition onto DG grid.
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid_dg, &basis, 5, 1, proj_fn, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &range_dg, dg_data);
  gkyl_proj_on_basis_release(proj);

  // Project initial condition onto finite volume grid for comparison.
  gkyl_fv_proj *fv_proj = gkyl_fv_proj_new(&grid_fv, 2, 1, proj_fn, &ctx);
  gkyl_fv_proj_advance(fv_proj, 0.0, &range_fv, fv_data_exact);
  gkyl_fv_proj_release(fv_proj);

  // Interpolate modal DG data onto finite volume grid.
  struct gkyl_dg_fv_interpolate *interp = gkyl_dg_fv_interpolate_new(&basis, &grid_dg, &grid_fv);

  gkyl_dg_fv_interpolate_modal_to_fv_advance(interp, &range_dg, &range_fv, dg_data, 0, fv_data, 0, 1);

  // Check the interpolated finite volume cells.
  int idx_dg[NUM_CHECK_CELLS][GKYL_MAX_DIM];
  check_cells(ndim, cells_dg[0], idx_dg);

  // Number of finite volume cells per DG cell.
  int num_fv = 1;
  for (int d=0; d<ndim; ++d) num_fv *= 2;
  int offset = 0;

  for (int i=0; i<NUM_CHECK_CELLS; ++i) {
    for (int k=0; k<num_fv; ++k) {
      int idx_fv[GKYL_MAX_DIM];

      // Fetch index for each finite volume cell overlapping the DG cell.
      for (int d=0; d<ndim; ++d) {
        offset = (k / (int)pow(2, ndim-1-d)) % 2;
        idx_fv[d] = 2*idx_dg[i][d] - 1 + offset;
      }

      long lin_fv = gkyl_range_idx(&range_fv, idx_fv);
      const double *fv_d = gkyl_array_cfetch(fv_data, lin_fv);
      TEST_CHECK( gkyl_compare_double(fv_ref[i][k], fv_d[0], 1e-12) );
      TEST_MSG("DG cell %d, FV cell %d: expected %.16e, got %.16e", i, k, fv_ref[i][k], fv_d[0]);
    }
  }

  // Write interpolated and exact FV data.
  char fname[256];
  snprintf(fname, sizeof(fname), "modal_to_fv_%dx_%d.gkyl", ndim, num_cells_fv);
  gkyl_grid_sub_array_write(&grid_fv, &range_fv, 0, fv_data, fname);
  snprintf(fname, sizeof(fname), "modal_to_fv_exact_%dx_%d.gkyl", ndim, num_cells_fv);
  gkyl_grid_sub_array_write(&grid_fv, &range_fv, 0, fv_data_exact, fname);

  gkyl_dg_fv_interpolate_release(interp);
  gkyl_array_release(dg_data);
  gkyl_array_release(fv_data);
  gkyl_array_release(fv_data_exact);
}

void
test_interp_fv_to_modal(int ndim, int num_cells_fv, const double dg_ref[NUM_CHECK_CELLS][8])
{
  struct interp_ctx ctx = {.ndim = ndim};

  // Initialize grids.
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM] = {0.0};
  int cells_dg[GKYL_MAX_DIM], cells_fv[GKYL_MAX_DIM] = {0};
  int ghost[GKYL_MAX_DIM] = {0};

  double L = 1.0;
  for (int d=0; d<ndim; ++d) {
    lower[d] = -L;
    upper[d] = L;
    cells_fv[d] = num_cells_fv;
    cells_dg[d] = (int)(num_cells_fv/2);
  }

  struct gkyl_rect_grid grid_dg, grid_fv;
  gkyl_rect_grid_init(&grid_dg, ndim, lower, upper, cells_dg);
  gkyl_rect_grid_init(&grid_fv, ndim, lower, upper, cells_fv);

  // Construct grid ranges.
  struct gkyl_range range_dg, range_dg_ext, range_fv, range_fv_ext;
  gkyl_create_grid_ranges(&grid_dg, ghost, &range_dg_ext, &range_dg);
  gkyl_create_grid_ranges(&grid_fv, ghost, &range_fv_ext, &range_fv);

  // Initialize basis, always p=1.
  struct gkyl_basis basis;
  int poly_order = 1;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Initialize arrays for FV and DG data.
  struct gkyl_array *fv_data = mkarr(1, range_fv.volume);
  struct gkyl_array *dg_data = mkarr(basis.num_basis, range_dg.volume);
  struct gkyl_array *dg_data_exact = mkarr(basis.num_basis, range_dg.volume);

  // Project initial condition onto finite volume grid.
  gkyl_fv_proj *fv_proj = gkyl_fv_proj_new(&grid_fv, 2, 1, proj_fn, &ctx);
  gkyl_fv_proj_advance(fv_proj, 0.0, &range_fv, fv_data);
  gkyl_fv_proj_release(fv_proj);

  // Project initial condition onto DG grid for comparison.
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid_dg, &basis, 5, 1, proj_fn, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &range_dg, dg_data_exact);
  gkyl_proj_on_basis_release(proj);

  // Interpolate finite volume data onto modal DG grid.
  struct gkyl_dg_fv_interpolate *interp = gkyl_dg_fv_interpolate_new(&basis, &grid_dg, &grid_fv);

  gkyl_dg_fv_interpolate_fv_to_modal_advance(interp, &range_dg, &range_fv, fv_data, 0, dg_data, 0, 1);

  // Check the DG coefficients in three representative cells.
  int idx_dg[NUM_CHECK_CELLS][GKYL_MAX_DIM];
  check_cells(ndim, cells_dg[0], idx_dg);

  for (int i=0; i<NUM_CHECK_CELLS; ++i) {
    const double *dg_d = gkyl_array_cfetch(dg_data, gkyl_range_idx(&range_dg, idx_dg[i]));

    // Check each expansion coefficient.
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare_double(dg_ref[i][k], dg_d[k], 1e-12) );
      TEST_MSG("DG cell %d, coefficient %d: expected %.16e, got %.16e", i, k, dg_ref[i][k], dg_d[k]);
    }
  }

  // Write interpolated and exact DG data.
  struct gkyl_msgpack_map_elem io_meta[] = {
    {.key = "polyOrder", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = basis.poly_order},
    {.key = "basisType", .elem_type = GKYL_MP_STRING, .cval = basis.id}
  };
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create(sizeof(io_meta)/sizeof(io_meta[0]), io_meta);

  char fname[256];
  snprintf(fname, sizeof(fname), "fv_to_modal_%dx_%d.gkyl", ndim, num_cells_fv);
  gkyl_grid_sub_array_write(&grid_dg, &range_dg, mt, dg_data, fname);
  snprintf(fname, sizeof(fname), "fv_to_modal_exact_%dx_%d.gkyl", ndim, num_cells_fv);
  gkyl_grid_sub_array_write(&grid_dg, &range_dg, mt, dg_data_exact, fname);

  gkyl_msgpack_data_release(mt);

  gkyl_dg_fv_interpolate_release(interp);
  gkyl_array_release(fv_data);
  gkyl_array_release(dg_data);
  gkyl_array_release(dg_data_exact);
}

// NOTE: Cell number is for the refined grid. Reference values are listed for
// the cells given by check_cells: FV cell averages in row-major order for
// modal_to_fv, DG coefficients for fv_to_modal.

void test_modal_to_fv_1x_p1_4()
{
  const double ref[NUM_CHECK_CELLS][8] = {
    { -6.3661980747221925e-01, -6.3661980747221902e-01 },
    {  6.3661980747221891e-01,  6.3661980747221913e-01 },
    {  6.3661980747221902e-01,  6.3661980747221913e-01 },
  };
  test_interp_modal_to_fv(1, 4, ref);
}

void test_modal_to_fv_1x_p1_8()
{
  const double ref[NUM_CHECK_CELLS][8] = {
    { -3.7569522771591340e-01, -8.9754431706962534e-01 },
    { -8.9754431706962512e-01, -3.7569522771591313e-01 },
    {  8.9754431706962512e-01,  3.7569522771591329e-01 },
  };
  test_interp_modal_to_fv(1, 8, ref);
}

void test_modal_to_fv_1x_p1_16()
{
  const double ref[NUM_CHECK_CELLS][8] = {
    { -1.9430201041317435e-01, -5.5154444674296599e-01 },
    { -8.2632898506563535e-01, -9.7430364724864182e-01 },
    {  5.5154444674296599e-01,  1.9430201041317435e-01 },
  };
  test_interp_modal_to_fv(1, 16, ref);
}

void test_modal_to_fv_1x_p1_32()
{
  const double ref[NUM_CHECK_CELLS][8] = {
    { -9.7921520152287819e-02, -2.8975683732914048e-01 },
    { -4.7069221385132060e-01, -6.3332234297947809e-01 },
    {  2.8975683732914048e-01,  9.7921520152287819e-02 },
  };
  test_interp_modal_to_fv(1, 32, ref);
}

void test_modal_to_fv_2x_p1_32()
{
  const double ref[NUM_CHECK_CELLS][8] = {
    {  9.7292331980176741e-02,  9.3555801772444463e-02,  2.8789502416945639e-01,  2.7683836191693512e-01 },
    {  4.1443049731820625e-01,  3.6328226748517151e-01,  5.5762148988218352e-01,  4.8880089798825732e-01 },
    { -2.7683836191693512e-01, -2.8789502416945639e-01, -9.3555801772444477e-02, -9.7292331980176797e-02 },
  };
  test_interp_modal_to_fv(2, 32, ref);
}

void test_modal_to_fv_3x_p1_32()
{
  const double ref[NUM_CHECK_CELLS][8] = {
    { -8.8479745843895022e-02, -9.4596660696240301e-02, -8.5081664655075856e-02, -9.0963658248379153e-02, -2.6181794648961237e-01, -2.7991833850835096e-01, -2.5176277928296958e-01, -2.6916802236064269e-01 },
    { -4.1354184310119418e-01, -4.0824026870858149e-01, -3.6250328929448467e-01, -3.5785602520785809e-01, -5.5642579436342487e-01, -5.4929245878445543e-01, -4.8775277295381281e-01, -4.8149982018940535e-01 },
    {  1.8895221416021171e-01,  2.2465434604759424e-01,  1.9649878682221295e-01,  2.3362682807865973e-01,  6.3855225013003752e-02,  7.5920538326452083e-02,  6.6405542285287625e-02,  7.8952732797235192e-02 },
  };
  test_interp_modal_to_fv(3, 32, ref);
}

void test_fv_to_modal_1x_p1_4()
{
  const double ref[NUM_CHECK_CELLS][8] = {
    { -8.9894118508375098e-01,  0.0000000000000000e+00 },
    {  8.9894118508375087e-01, 0.0000000000000000e+00 },
    {  8.9894118508375087e-01,  0.0000000000000000e+00 },
  };
  test_interp_fv_to_modal(1, 4, ref);
}

void test_fv_to_modal_1x_p1_8()
{
  const double ref[NUM_CHECK_CELLS][8] = {
    { -9.0023543783430482e-01, -4.3057596930801562e-01 },
    { -9.0023543783430471e-01,  4.3057596930801578e-01 },
    {  9.0023543783430482e-01, -4.3057596930801562e-01 },
  };
  test_interp_fv_to_modal(1, 8, ref);
}

void test_fv_to_modal_1x_p1_16()
{
  const double ref[NUM_CHECK_CELLS][8] = {
    { -5.2739017001841149e-01, -2.9244141097468407e-01 },
    { -1.2732325011207011e+00, -1.2113319862523829e-01 },
    {  5.2739017001841149e-01, -2.9244141097468407e-01 },
  };
  test_interp_fv_to_modal(1, 16, ref);
}

void test_fv_to_modal_1x_p1_32()
{
  const double ref[NUM_CHECK_CELLS][8] = {
    { -2.7412990106130614e-01, -1.5673367505835104e-01 },
    { -7.8065591074107543e-01, -1.3287239382501956e-01 },
    {  2.7412990106130602e-01, -1.5673367505835101e-01 },
  };
  test_interp_fv_to_modal(1, 32, ref);
}

void test_fv_to_modal_2x_p1_32()
{
  const double ref[NUM_CHECK_CELLS][8] = {
    {  3.7779049963469713e-01,  2.1600158603870973e-01, -8.5463468110748980e-03, -4.8863708002558157e-03 },
    {  9.1206694795379328e-01,  1.5523935326161298e-01, -6.9308579588171784e-02, -1.1796742656761523e-02 },
    { -3.7779049963469713e-01,  2.1600158603870973e-01, -8.5463468110748980e-03,  4.8863708002558157e-03 },
  };
  test_interp_fv_to_modal(2, 32, ref);
}

void test_fv_to_modal_3x_p1_32()
{
  const double ref[NUM_CHECK_CELLS][8] = {
    { -5.0267773714261343e-01, -2.8740582040613599e-01,  1.1371536023222332e-02, -1.9406115284593798e-02,  6.5016717442270737e-03, -1.1095439626924190e-02,  4.3900360553855320e-04,  2.5100015793067855e-04 },
    { -1.2788409865182453e+00, -2.1766652998107891e-01,  9.7179984971022343e-02,  9.5326494590851057e-03,  1.6540625719110101e-02,  1.6225150360047865e-03, -7.2439243106378864e-04, -1.2329600667832703e-04 },
    {  3.9911406485585255e-01, -2.2819332699625935e-01,  9.0287268174667860e-03,  3.9806378125282829e-02, -5.1621714002041052e-03, -2.2759282771355357e-02,  9.0049673848446735e-04, -5.1485869529118810e-04 },
  };
  test_interp_fv_to_modal(3, 32, ref);
}

TEST_LIST = {
  { "test_modal_to_fv_1x_p1_4", test_modal_to_fv_1x_p1_4 },
  { "test_modal_to_fv_1x_p1_8", test_modal_to_fv_1x_p1_8 },
  { "test_modal_to_fv_1x_p1_16", test_modal_to_fv_1x_p1_16 },
  { "test_modal_to_fv_1x_p1_32", test_modal_to_fv_1x_p1_32 },
  { "test_modal_to_fv_2x_p1_32", test_modal_to_fv_2x_p1_32 },
  { "test_modal_to_fv_3x_p1_32", test_modal_to_fv_3x_p1_32 },

  { "test_fv_to_modal_1x_p1_4", test_fv_to_modal_1x_p1_4 },
  { "test_fv_to_modal_1x_p1_8", test_fv_to_modal_1x_p1_8 },
  { "test_fv_to_modal_1x_p1_16", test_fv_to_modal_1x_p1_16 },
  { "test_fv_to_modal_1x_p1_32", test_fv_to_modal_1x_p1_32 },
  { "test_fv_to_modal_2x_p1_32", test_fv_to_modal_2x_p1_32 },
  { "test_fv_to_modal_3x_p1_32", test_fv_to_modal_3x_p1_32 },
  { NULL, NULL }
};
