// Test creation and deallocation of updater that applies the
// twist shift BCs.
//
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_basis.h>
#include <gkyl_proj_on_basis.h>
#include <mpack.h>
#include <gkyl_array_rio.h>
#include <gkyl_twistshift_dg.h>
#include <gkyl_velocity_map.h>
#include <gkyl_position_map.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>
#include <gkyl_math.h>
#include <gkyl_const.h>

// Compare two DG coefficients with a mixed relative/absolute tolerance.
// Needed to pass the CBC 3x test on stellar nvcc installation.
static bool
ts_compare_coeff(double ref, double val)
{
  return gkyl_compare(ref, val, 1e-13) || fabs(ref-val) < 1e-14;
}

// Meta-data for IO
struct test_bc_twistshift_output_meta {
  int poly_order; // polynomial order
  const char *basis_type; // name of basis functions
};

// returned gkyl_array_meta must be freed using gyrokinetic_array_meta_release
static struct gkyl_msgpack_data*
test_bc_twistshift_array_meta_new(struct test_bc_twistshift_output_meta meta)
{
  struct gkyl_msgpack_data *mt = gkyl_malloc(sizeof(*mt));

  mt->meta_sz = 0;
  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mt->meta, &mt->meta_sz);

  // add some data to mpack
  mpack_build_map(&writer);

  mpack_write_cstr(&writer, "polyOrder");
  mpack_write_i64(&writer, meta.poly_order);

  mpack_write_cstr(&writer, "basisType");
  mpack_write_cstr(&writer, meta.basis_type);

  mpack_complete_map(&writer);

  int status = mpack_writer_destroy(&writer);

  if (status != mpack_ok) {
    free(mt->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mt);
    mt = 0;
  }

  return mt;
}

static void
test_bc_twistshift_array_meta_release(struct gkyl_msgpack_data *mt)
{
  if (!mt) return;
  MPACK_FREE(mt->meta);
  gkyl_free(mt);
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}
// Apply periodic BCs along parallel direction
void
apply_periodic_bc(struct gkyl_array *buff, struct gkyl_array *fld, const int dir, const struct skin_ghost_ranges sgr)
{
  gkyl_array_copy_to_buffer(buff->data, fld, &(sgr.lower_skin[dir]));
  gkyl_array_copy_from_buffer(fld, buff->data, &(sgr.upper_ghost[dir]));

  gkyl_array_copy_to_buffer(buff->data, fld, &(sgr.upper_skin[dir]));
  gkyl_array_copy_from_buffer(fld, buff->data, &(sgr.lower_ghost[dir]));
}

static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct test_bc_twistshift_ctx {
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM];
  double B0;
  double vt;
  double mass;
  enum gkyl_edge_loc edge;
};

void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void eval_bfield_3x(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

  struct test_bc_twistshift_ctx *pars = ctx;
  double B0 = pars->B0;

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = B0;
}

void
shift1_fig6(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct test_bc_twistshift_ctx *pars = ctx;
  double Lx[2] = {pars->upper[0]-pars->lower[0], pars->upper[1]-pars->lower[1]};
  double dx[2] = {Lx[0]/pars->cells[0], Lx[1]/pars->cells[1]};

  fout[0] = 4.0*dx[1];
}

void
shift1m_fig6(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  shift1_fig6(t, xn, fout, ctx);
  fout[0] *= -1.0;
}

void
shift2_fig6(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1.1;
}

void
shift2m_fig6(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  shift2_fig6(t, xn, fout, ctx);
  fout[0] *= -1.0;
}

void
init_donor_fig6(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double y = xn[1];

  double mu = 0.0;
  double sigma = 0.3;

  fout[0] = ( 1.0/sqrt(2.0*M_PI*pow(sigma,2)) ) * exp( -pow(y-mu,2)/(2.0*pow(sigma,2)) );
}

void
shift_fig9(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct test_bc_twistshift_ctx *pars = ctx;
  double Lx[2] = {pars->upper[0]-pars->lower[0], pars->upper[1]-pars->lower[1]};
  double dx[2] = {Lx[0]/pars->cells[0], Lx[1]/pars->cells[1]};

  fout[0] = dx[1]/2.0;
}

void
shiftm_fig9(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  shift_fig9(t, xn, fout, ctx);
  fout[0] *= -1.0;
}

void
init_donor_fig9(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double y = xn[1];

  struct test_bc_twistshift_ctx *pars = ctx;
  double ymid = 0.5*(pars->upper[1]+pars->lower[1]);
  double dy = (pars->upper[1]-pars->lower[1])/pars->cells[1];

  fout[0] = 0.;
  if (ymid < y && y < ymid+dy)
    fout[0] = 1.;
}

void
test_bc_twistshift_3x_fig6_wcells(const int *cells, enum gkyl_edge_loc edge,
  bool check_distf, bool use_gpu, bool write_f)
{
  double vt = 1.0; // Thermal speed.
  double mass = 1.0;
  double B0 = 1.0; // Magnetic field magnitude.
  int bc_dir = 2; // Direction in which to apply TS.

  int poly_order = 1;
  const double lower[] = {-2.0, -1.50, -3.0};
  const double upper[] = { 2.0,  1.50,  3.0};
  int vdim = 0;
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int cdim = ndim - vdim;

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d=0; d<cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = cells[d];
  }

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  // Ranges.
  int ghost_conf[cdim];
  for (int d=0; d<cdim; d++) ghost_conf[d] = 1;
  struct gkyl_range local_conf, local_ext_conf; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);

  int ghost[ndim];
  for (int d=0; d<cdim; d++) ghost[d] = ghost_conf[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);
  struct skin_ghost_ranges skin_ghost_conf; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost_conf, &local_ext_conf, ghost_conf);

  // Pick skin and ghost ranges based on 'edge'.
  struct gkyl_range skin_rng, ghost_rng;
  struct gkyl_range skin_rng_conf, ghost_rng_conf;
  if (edge == GKYL_LOWER_EDGE) {
    skin_rng = skin_ghost.upper_skin[bc_dir];
    ghost_rng = skin_ghost.lower_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.upper_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.lower_ghost[bc_dir];
  }
  else {
    skin_rng = skin_ghost.lower_skin[bc_dir];
    ghost_rng = skin_ghost.upper_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.lower_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.upper_ghost[bc_dir];
  }

  struct test_bc_twistshift_ctx proj_ctx = {
    .lower = {lower[0], lower[1], lower[2]},
    .upper = {upper[0], upper[1], upper[2]},
    .cells = {cells[0], cells[1], cells[2]},
    .B0 = B0,
    .vt = vt,
    .mass = mass,
  };

  // Initialize the distribution
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, basis.num_basis, local_ext.volume) : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .num_ret_vals = 1,
      .eval = init_donor_fig6,
//      .eval = init_donor_fig9,
      .ctx = &proj_ctx,
    }
  );
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);
  struct gkyl_msgpack_data *mt = test_bc_twistshift_array_meta_new( (struct test_bc_twistshift_output_meta) {
      .poly_order = poly_order,
      .basis_type = basis.id
    }
  );
  if (write_f)
    gkyl_grid_sub_array_write(&grid, &local, mt, distf_ho, "ctest_bc_twistshift_3x_fig6_do.gkyl");

  // Create a range only extended in bc_dir.
  struct gkyl_range update_rng;
  int lower_bcdir_ext[ndim], upper_bcdir_ext[ndim];
  for (int d=0; d<ndim; d++) {
    lower_bcdir_ext[d] = local.lower[d];
    upper_bcdir_ext[d] = local.upper[d];
  }
  lower_bcdir_ext[bc_dir] = local.lower[bc_dir] - ghost[bc_dir];
  upper_bcdir_ext[bc_dir] = local.upper[bc_dir] + ghost[bc_dir];
  gkyl_sub_range_init(&update_rng, &local_ext, lower_bcdir_ext, upper_bcdir_ext);

  // Create the twist-shift updater and shift the donor field.
  struct gkyl_twistshift_dg_inp tsinp = {
    .bc_dir = bc_dir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = edge,
    .cdim = cdim,
    .bcdir_ext_update_r = &update_rng,
    .num_ghost = ghost,
    .basis = &basis,
    .grid = &grid,
    .shift_func = shift1_fig6,
//    .shift_func = shift_fig9,
    .shift_func_ctx = &proj_ctx,
    .use_gpu = use_gpu,
  };

  struct gkyl_twistshift_dg *tsup = gkyl_twistshift_dg_inew(&tsinp);

  // First apply periodicity in z.
  struct gkyl_array *buff_per = mkarr(use_gpu, basis.num_basis, skin_rng.volume);
  apply_periodic_bc(buff_per, distf, bc_dir, skin_ghost);

  gkyl_twistshift_dg_advance(tsup, distf, distf);
  gkyl_array_copy(distf_ho, distf);

  // Write out the target in the extended range.
  if (write_f) {
    double lower_ext[ndim], upper_ext[ndim];
    int cells_ext[ndim];
    for (int d=0; d<ndim; d++) {
      double dx = (upper[d]-lower[d])/cells[d];
      lower_ext[d] = lower[d]-dx*ghost[d];
      upper_ext[d] = upper[d]+dx*ghost[d];
      cells_ext[d] = cells[d]+2*ghost[d];
    }
    struct gkyl_rect_grid grid_ext;
    gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);
    gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "ctest_bc_twistshift_3x_fig6_tar.gkyl");
  }

  if (check_distf) {
    // Check 0th and 2nd DG coeffs.
    const double f0[] =
    {
      2.5656054446469541e+00, 4.0323377043325592e-01, 2.4549376970832818e-02,
      5.6930416832388616e-04, 5.6930416832388561e-04, 2.4549376970832790e-02,
      4.0323377043325587e-01, 2.5656054446469554e+00, 6.4341275514494081e+00,
      6.4341275514494081e+00,
    };
    const double f2[] =
    {
      -1.0463462350335195e+00, -2.4917980000416121e-01, -1.8802664978950782e-02,
      -4.9044149652606374e-04,  4.9044149652607144e-04,  1.8802664978951094e-02,
       2.4917980000416726e-01,  1.0463462350335588e+00,  9.2229040193425293e-01,
      -9.2229040193415335e-01,
    };

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &ghost_rng);
    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(&ghost_rng, iter.idx);
      double *f_c = gkyl_array_fetch(distf_ho, linidx);
      int refidx = (iter.idx[1]-1)*cells[0] + iter.idx[0]-1;
      TEST_CHECK( gkyl_compare(f0[refidx]/pow(sqrt(2.0),2), f_c[0], 1e-13) );
      TEST_CHECK( gkyl_compare(f2[refidx]/pow(sqrt(2.0),2), f_c[2], 1e-13) );
    }
  }

  // Copy the ghost cell back to the skin cell, and apply the negated shift to
  // see if we recover the donor field approximately.
  gkyl_array_copy_range_to_range(distf, distf, &skin_rng, &ghost_rng);
  tsinp.shift_func = shift1m_fig6;
//  tsinp.shift_func = shiftm_fig9;
  struct gkyl_twistshift_dg *tsup_m = gkyl_twistshift_dg_inew(&tsinp);
  gkyl_twistshift_dg_advance(tsup_m, distf, distf);
  gkyl_array_copy(distf_ho, distf);

  if (write_f) {
    double lower_ext[ndim], upper_ext[ndim];
    int cells_ext[ndim];
    for (int d=0; d<ndim; d++) {
      double dx = (upper[d]-lower[d])/cells[d];
      lower_ext[d] = lower[d]-dx*ghost[d];
      upper_ext[d] = upper[d]+dx*ghost[d];
      cells_ext[d] = cells[d]+2*ghost[d];
    }
    struct gkyl_rect_grid grid_ext;
    gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);
    gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "ctest_bc_twistshift_3x_fig6_tar_shifted.gkyl");
  }

  if (check_distf) {
    // Check 0th and 2nd DG coeffs.
    const double f0[] =
    {
      5.6930416832388551e-04, 2.4549376970832783e-02, 4.0323377043325570e-01,
      2.5656054446469541e+00, 6.4341275514494018e+00, 6.4341275514494001e+00,
      2.5656054446469501e+00, 4.0323377043325520e-01, 2.4549376970832770e-02,
      5.6930416832388496e-04,
    };
    const double f2[] =
    {
      4.9044149652607047e-04,  1.8802664978951059e-02,  2.4917980000416681e-01,
      1.0463462350335573e+00,  9.2229040193425282e-01, -9.2229040193414913e-01,
     -1.0463462350335162e+00, -2.4917980000416048e-01, -1.8802664978950730e-02,
     -4.9044149652606244e-04,
    };

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &ghost_rng);
    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(&ghost_rng, iter.idx);
      double *f_c = gkyl_array_fetch(distf_ho, linidx);
      int refidx = (iter.idx[1]-1)*cells[0] + iter.idx[0]-1;
      TEST_CHECK( gkyl_compare(f0[refidx]/pow(sqrt(2.0),2), f_c[0], 1e-13) );
      TEST_CHECK( gkyl_compare(f2[refidx]/pow(sqrt(2.0),2), f_c[2], 1e-13) );
    }
  }

  gkyl_array_release(buff_per);
  test_bc_twistshift_array_meta_release(mt);
  gkyl_twistshift_dg_release(tsup);
  gkyl_twistshift_dg_release(tsup_m);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(distf);

}
void
test_bc_twistshift_3x2v_fig6_wcells(const int *cells, enum gkyl_edge_loc edge,
  bool check_distf, bool use_gpu, bool write_f)
{
  double vt = 1.0; // Thermal speed.
  double mass = 1.0;
  double B0 = 1.0; // Magnetic field magnitude.
  int bc_dir = 2; // Direction in which to apply TS.

  int poly_order = 1;
  const double lower[] = {-2.0, -1.50, -3.0, -5.0*vt, 0.};
  const double upper[] = { 2.0,  1.50,  3.0,  5.0*vt, mass*(pow(5.0*vt,2))/(2.0*B0)};
  int vdim = 2;
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int cdim = ndim - vdim;

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d=0; d<cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = cells[d];
  }
  double lower_vel[vdim], upper_vel[vdim];
  int cells_vel[vdim];
  for (int d=0; d<vdim; d++) {
    lower_vel[d] = lower[cdim+d];
    upper_vel[d] = upper[cdim+d];
    cells_vel[d] = cells[cdim+d];
  }

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  struct gkyl_rect_grid grid_vel;
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

  // Basis functions.
  struct gkyl_basis basis;
  if (poly_order == 1) 
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  // Ranges.
  int ghost_conf[cdim];
  for (int d=0; d<cdim; d++) ghost_conf[d] = 1;
  struct gkyl_range local_conf, local_ext_conf; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);

  int ghost_vel[vdim];
  for (int d=0; d<vdim; d++) ghost_vel[d] = 0;
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_ext_vel, &local_vel);

  int ghost[ndim];
  for (int d=0; d<cdim; d++) ghost[d] = ghost_conf[d];
  for (int d=cdim; d<ndim; d++) ghost[d] = ghost_vel[d-cdim];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);
  struct skin_ghost_ranges skin_ghost_conf; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost_conf, &local_ext_conf, ghost_conf);

  // Pick skin and ghost ranges based on 'edge'.
  struct gkyl_range skin_rng, ghost_rng;
  struct gkyl_range skin_rng_conf, ghost_rng_conf;
  if (edge == GKYL_LOWER_EDGE) {
    skin_rng = skin_ghost.upper_skin[bc_dir];
    ghost_rng = skin_ghost.lower_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.upper_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.lower_ghost[bc_dir];
  }
  else {
    skin_rng = skin_ghost.lower_skin[bc_dir];
    ghost_rng = skin_ghost.upper_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.lower_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.upper_ghost[bc_dir];
  }

  struct test_bc_twistshift_ctx proj_ctx = {
    .lower = {lower[0], lower[1], lower[2], lower[3], lower[4]},
    .upper = {upper[0], upper[1], upper[2], upper[3], upper[4]},
    .cells = {cells[0], cells[1], cells[2], cells[3], cells[4]},
    .B0 = B0,
    .vt = vt,
    .mass = mass,
  };

  // Initialize the distribution
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, basis.num_basis, local_ext.volume) : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .num_ret_vals = 1,
      .eval = init_donor_fig6,
//      .eval = init_donor_fig9,
      .ctx = &proj_ctx,
    }
  );
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);
  struct gkyl_msgpack_data *mt = test_bc_twistshift_array_meta_new( (struct test_bc_twistshift_output_meta) {
      .poly_order = poly_order,
      .basis_type = basis.id
    }
  );
  if (write_f)
    gkyl_grid_sub_array_write(&grid, &local, mt, distf_ho, "ctest_bc_twistshift_3x2v_fig6_do.gkyl");

  // Create a range only extended in bc_dir.
  struct gkyl_range update_rng;
  int lower_bcdir_ext[ndim], upper_bcdir_ext[ndim];
  for (int d=0; d<ndim; d++) {
    lower_bcdir_ext[d] = local.lower[d];
    upper_bcdir_ext[d] = local.upper[d];
  }
  lower_bcdir_ext[bc_dir] = local.lower[bc_dir] - ghost[bc_dir];
  upper_bcdir_ext[bc_dir] = local.upper[bc_dir] + ghost[bc_dir];
  gkyl_sub_range_init(&update_rng, &local_ext, lower_bcdir_ext, upper_bcdir_ext);

  // Create the twist-shift updater and shift the donor field.
  struct gkyl_twistshift_dg_inp tsinp = {
    .bc_dir = bc_dir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = edge,
    .cdim = cdim,
    .bcdir_ext_update_r = &update_rng,
    .num_ghost = ghost,
    .basis = &basis,
    .grid = &grid,
    .shift_func = shift1_fig6,
//    .shift_func = shift_fig9,
    .shift_func_ctx = &proj_ctx,
    .use_gpu = use_gpu,
  };

  struct gkyl_twistshift_dg *tsup = gkyl_twistshift_dg_inew(&tsinp);

  // First apply periodicity in z.
  struct gkyl_array *buff_per = mkarr(use_gpu, basis.num_basis, skin_rng.volume);
  apply_periodic_bc(buff_per, distf, bc_dir, skin_ghost);

  gkyl_twistshift_dg_advance(tsup, distf, distf);
  gkyl_array_copy(distf_ho, distf);

  // Write out the target in the extended range.
  if (write_f) {
    double lower_ext[ndim], upper_ext[ndim];
    int cells_ext[ndim];
    for (int d=0; d<ndim; d++) {
      double dx = (upper[d]-lower[d])/cells[d];
      lower_ext[d] = lower[d]-dx*ghost[d];
      upper_ext[d] = upper[d]+dx*ghost[d];
      cells_ext[d] = cells[d]+2*ghost[d];
    }
    struct gkyl_rect_grid grid_ext;
    gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);
    gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "ctest_bc_twistshift_3x2v_fig6_tar.gkyl");
  }

  // Compute the integrated moments of the skin cell and the ghost cell.
  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_ext_vel, use_gpu);
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .c2p_ctx = 0,
    .mapc2p = mapc2p,
    .bfield_ctx = &proj_ctx,
    .bfield_func = eval_bfield_3x,
    .position_map = pmap,
    .grid = grid_conf,
    .local = local_conf,
    .local_ext = local_ext_conf,
    .global = local_conf,
    .global_ext = local_ext_conf,
    .basis = basis_conf,
    .geo_grid = grid_conf,
    .geo_local = local_conf,
    .geo_local_ext = local_ext_conf,
    .geo_global = local_conf,
    .geo_global_ext = local_ext_conf,
    .geo_basis = basis_conf,
  };
  struct gk_geometry* gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_inp);
  struct gk_geometry* gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);
  gkyl_gk_geometry_release(gk_geom_3d); // release temporary 3d geometry
  if (use_gpu) {  // If we are on the gpu, copy from host
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_inp, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }
  // Need the magnetic field to be initialized in the ghost cell in order to
  // compute integrated moments in the ghost cell (for checking moment
  // conservation).
  gkyl_array_clear(gk_geom->geo_corn.bmag, 0.0);
  gkyl_array_shiftc(gk_geom->geo_corn.bmag, B0*pow(sqrt(2.0),cdim), 0);

  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_gyrokinetic_new(&grid, &basis_conf,
    &basis, &local_conf, mass, 0, gvm, gk_geom, NULL, GKYL_F_MOMENT_M0M1M2, true, use_gpu);
  int num_mom = gkyl_dg_updater_moment_gyrokinetic_num_mom(mcalc);

  struct gkyl_array *marr = mkarr(use_gpu, num_mom, local_ext_conf.volume);
  double *red_integ_mom_skin, *red_integ_mom_ghost;
  if (use_gpu) {
    red_integ_mom_skin = gkyl_cu_malloc(sizeof(double[num_mom]));
    red_integ_mom_ghost = gkyl_cu_malloc(sizeof(double[num_mom]));
  }
  else {
    red_integ_mom_skin = gkyl_malloc(sizeof(double[num_mom]));
    red_integ_mom_ghost = gkyl_malloc(sizeof(double[num_mom]));
  }
  double *red_integ_mom_skin_ho = gkyl_malloc(sizeof(double[num_mom]));
  double *red_integ_mom_ghost_ho = gkyl_malloc(sizeof(double[num_mom]));

  gkyl_dg_updater_moment_gyrokinetic_advance(mcalc,
      &skin_rng, &skin_rng_conf, distf, marr);
  gkyl_array_reduce_range(red_integ_mom_skin, marr, GKYL_SUM, &skin_rng_conf);

  gkyl_dg_updater_moment_gyrokinetic_advance(mcalc,
      &ghost_rng, &ghost_rng_conf, distf, marr);
  gkyl_array_reduce_range(red_integ_mom_ghost, marr, GKYL_SUM, &ghost_rng_conf);

  if (use_gpu) {
    gkyl_cu_memcpy(red_integ_mom_skin_ho, red_integ_mom_skin, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(red_integ_mom_ghost_ho, red_integ_mom_ghost, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    memcpy(red_integ_mom_skin_ho, red_integ_mom_skin, sizeof(double[num_mom]));
    memcpy(red_integ_mom_ghost_ho, red_integ_mom_ghost, sizeof(double[num_mom]));
  }

  for (int k=0; k<num_mom; k++) {
    TEST_CHECK( gkyl_compare(red_integ_mom_skin_ho[k], red_integ_mom_ghost_ho[k], 1e-12));
    TEST_MSG( "integ_mom %d | Expected: %.14e | Got: %.14e\n",k,red_integ_mom_skin_ho[k],red_integ_mom_ghost_ho[k]);
  }

  if (check_distf) {
    // Check 0th and 2nd DG coeffs.
    const double f0[] =
    {
      2.5656054446469541e+00, 4.0323377043325592e-01, 2.4549376970832818e-02,
      5.6930416832388616e-04, 5.6930416832388561e-04, 2.4549376970832790e-02,
      4.0323377043325587e-01, 2.5656054446469554e+00, 6.4341275514494081e+00,
      6.4341275514494081e+00,
    };
    const double f2[] =
    {
      -1.0463462350335195e+00, -2.4917980000416121e-01, -1.8802664978950782e-02,
      -4.9044149652606374e-04,  4.9044149652607144e-04,  1.8802664978951094e-02,
       2.4917980000416726e-01,  1.0463462350335588e+00,  9.2229040193425293e-01,
      -9.2229040193415335e-01,
    };

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &ghost_rng);
    while (gkyl_range_iter_next(&iter)) {
      if (iter.idx[3]==1 && iter.idx[4]==1) {
        long linidx = gkyl_range_idx(&ghost_rng, iter.idx);
        double *f_c = gkyl_array_fetch(distf_ho, linidx);
        int refidx = (iter.idx[1]-1)*cells[0] + iter.idx[0]-1;
        TEST_CHECK( gkyl_compare(f0[refidx], f_c[0], 1e-13) );
        TEST_CHECK( gkyl_compare(f2[refidx], f_c[2], 1e-13) );
      }
    }
  }

  // Copy the ghost cell back to the skin cell, and apply the negated shift to
  // see if we recover the donor field approximately.
  gkyl_array_copy_range_to_range(distf, distf, &skin_rng, &ghost_rng);
  tsinp.shift_func = shift1m_fig6;
//  tsinp.shift_func = shiftm_fig9;
  struct gkyl_twistshift_dg *tsup_m = gkyl_twistshift_dg_inew(&tsinp);
  gkyl_twistshift_dg_advance(tsup_m, distf, distf);
  gkyl_array_copy(distf_ho, distf);

  if (write_f) {
    double lower_ext[ndim], upper_ext[ndim];
    int cells_ext[ndim];
    for (int d=0; d<ndim; d++) {
      double dx = (upper[d]-lower[d])/cells[d];
      lower_ext[d] = lower[d]-dx*ghost[d];
      upper_ext[d] = upper[d]+dx*ghost[d];
      cells_ext[d] = cells[d]+2*ghost[d];
    }
    struct gkyl_rect_grid grid_ext;
    gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);
    gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "ctest_bc_twistshift_3x2v_fig6_tar_shifted.gkyl");
  }

  gkyl_dg_updater_moment_gyrokinetic_advance(mcalc,
      &ghost_rng, &ghost_rng_conf, distf, marr);
  gkyl_array_reduce_range(red_integ_mom_ghost, marr, GKYL_SUM, &ghost_rng_conf);

  if (use_gpu)
    gkyl_cu_memcpy(red_integ_mom_ghost_ho, red_integ_mom_ghost, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(red_integ_mom_ghost_ho, red_integ_mom_ghost, sizeof(double[num_mom]));

  for (int k=0; k<num_mom; k++) {
    TEST_CHECK( gkyl_compare(red_integ_mom_skin_ho[k], red_integ_mom_ghost_ho[k], 1e-12));
    TEST_MSG( "integ_mom %d | Expected: %.14e | Got: %.14e\n",k,red_integ_mom_skin_ho[k],red_integ_mom_ghost_ho[k]);
  }

  if (check_distf) {
    // Check 0th and 2nd DG coeffs.
    const double f0[] =
    {
      5.6930416832388551e-04, 2.4549376970832783e-02, 4.0323377043325570e-01,
      2.5656054446469541e+00, 6.4341275514494018e+00, 6.4341275514494001e+00,
      2.5656054446469501e+00, 4.0323377043325520e-01, 2.4549376970832770e-02,
      5.6930416832388496e-04,
    };
    const double f2[] =
    {
      4.9044149652607047e-04,  1.8802664978951059e-02,  2.4917980000416681e-01,
      1.0463462350335573e+00,  9.2229040193425282e-01, -9.2229040193414913e-01,
     -1.0463462350335162e+00, -2.4917980000416048e-01, -1.8802664978950730e-02,
     -4.9044149652606244e-04,
    };

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &ghost_rng);
    while (gkyl_range_iter_next(&iter)) {
      if (iter.idx[3]==1 && iter.idx[4]==1) {
        long linidx = gkyl_range_idx(&ghost_rng, iter.idx);
        double *f_c = gkyl_array_fetch(distf_ho, linidx);
        int refidx = (iter.idx[1]-1)*cells[0] + iter.idx[0]-1;
        TEST_CHECK( gkyl_compare(f0[refidx], f_c[0], 1e-13) );
        TEST_CHECK( gkyl_compare(f2[refidx], f_c[2], 1e-13) );
      }
    }
  }

  gkyl_free(red_integ_mom_skin_ho);
  gkyl_free(red_integ_mom_ghost_ho);
  if (use_gpu) {
    gkyl_cu_free(red_integ_mom_skin);
    gkyl_cu_free(red_integ_mom_ghost);
  }
  else {
    gkyl_free(red_integ_mom_skin);
    gkyl_free(red_integ_mom_ghost);
  }
  gkyl_dg_updater_moment_gyrokinetic_release(mcalc);
  gkyl_array_release(marr);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_position_map_release(pmap);
  gkyl_velocity_map_release(gvm);
  gkyl_array_release(buff_per);
  test_bc_twistshift_array_meta_release(mt);
  gkyl_twistshift_dg_release(tsup);
  gkyl_twistshift_dg_release(tsup_m);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(distf);

}

void
shift_fig11(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];

  struct test_bc_twistshift_ctx *pars = ctx;
  double Lx[2] = {pars->upper[0]-pars->lower[0], pars->upper[1]-pars->lower[1]};
  double dx[2] = {Lx[0]/pars->cells[0], Lx[1]/pars->cells[1]};

  fout[0] = 0.6*x+1.8;
}

void
init_donor_3x_fig11(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

  struct test_bc_twistshift_ctx *pars = ctx;
  double Lx[2] = {pars->upper[0]-pars->lower[0], pars->upper[1]-pars->lower[1]};
  double B0 = pars->B0;
  double vt = pars->vt;
  double mass = pars->mass;
  double vtsq = vt*vt;

  double beta[2] = {0.0, 0.0};
  double sigma[2] = {0.6, 0.2};

  fout[0] = ( 1.0/pow(sqrt(2.0*M_PI*vtsq),3) )
    * exp( -pow(x-beta[0],2)/(2.0*pow(sigma[0],2)) -pow(y-beta[1],2)/(2.0*pow(sigma[1],2)) );  
}

void
init_donor_3x2v_fig11(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2], vpar = xn[3], mu = xn[4];

  struct test_bc_twistshift_ctx *pars = ctx;
  double B0 = pars->B0;
  double vt = pars->vt;
  double mass = pars->mass;
  double vtsq = vt*vt;

  init_donor_3x_fig11(t, xn, fout, ctx);
  fout[0] *= exp( -(pow(vpar,2)+2.0*mu*B0/mass)/(2.0*vtsq) );  
}

void
test_bc_twistshift_3x_fig11_wcells(const int *cells, enum gkyl_edge_loc edge,
  int apply_in_half_x, bool check_distf, bool use_gpu, bool write_f)
{
  double vt = 1.0; // Thermal speed.
  double mass = 1.0;
  double B0 = 1.0; // Magnetic field magnitude.
  int bc_dir = 2; // Direction in which to apply TS.

  int poly_order = 1;
  const double lower[] = {-2.0, -1.50, -3.0};
  const double upper[] = { 2.0,  1.50,  3.0};
  int vdim = 0;
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int cdim = ndim - vdim;

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d=0; d<cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = cells[d];
  }

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  // Ranges.
  int ghost_conf[cdim];
  for (int d=0; d<cdim; d++) ghost_conf[d] = 1;
  struct gkyl_range local_conf, local_ext_conf; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);

  int ghost[ndim];
  for (int d=0; d<cdim; d++) ghost[d] = ghost_conf[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);
  struct skin_ghost_ranges skin_ghost_conf; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost_conf, &local_ext_conf, ghost_conf);

  // Pick skin and ghost ranges based on 'edge'.
  struct gkyl_range skin_rng, ghost_rng;
  struct gkyl_range skin_rng_conf, ghost_rng_conf;
  if (edge == GKYL_LOWER_EDGE) {
    skin_rng = skin_ghost.upper_skin[bc_dir];
    ghost_rng = skin_ghost.lower_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.upper_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.lower_ghost[bc_dir];
  }
  else {
    skin_rng = skin_ghost.lower_skin[bc_dir];
    ghost_rng = skin_ghost.upper_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.lower_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.upper_ghost[bc_dir];
  }

  struct test_bc_twistshift_ctx proj_ctx = {
    .lower = {lower[0], lower[1], lower[2]},
    .upper = {upper[0], upper[1], upper[2]},
    .cells = {cells[0], cells[1], cells[2]},
    .B0 = B0,
    .vt = vt,
    .mass = mass,
  };

  // Initialize the distribution
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, basis.num_basis, local_ext.volume) : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .num_ret_vals = 1,
      .eval = init_donor_3x_fig11,
      .ctx = &proj_ctx,
    }
  );
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);
  struct gkyl_msgpack_data *mt = test_bc_twistshift_array_meta_new( (struct test_bc_twistshift_output_meta) {
      .poly_order = poly_order,
      .basis_type = basis.id
    }
  );
  if (write_f)
    gkyl_grid_sub_array_write(&grid, &local, mt, distf_ho, "ctest_bc_twistshift_3x_fig11_do.gkyl");

  // Create a range only extended in bc_dir.
  struct gkyl_range update_rng;
  int lower_bcdir_ext[ndim], upper_bcdir_ext[ndim];
  for (int d=0; d<ndim; d++) {
    lower_bcdir_ext[d] = local.lower[d];
    upper_bcdir_ext[d] = local.upper[d];
  }
  lower_bcdir_ext[bc_dir] = local.lower[bc_dir] - ghost[bc_dir];
  upper_bcdir_ext[bc_dir] = local.upper[bc_dir] + ghost[bc_dir];
  gkyl_sub_range_init(&update_rng, &local_ext, lower_bcdir_ext, upper_bcdir_ext);

  if (apply_in_half_x < 0) {
    // Apply the BC only on the lower half of the domain.
    int x_half_len = (update_rng.upper[0] - update_rng.lower[0] + 1)/2;
    gkyl_range_shorten_from_above(&update_rng, &update_rng, 0, x_half_len);
  }
  else if (apply_in_half_x > 0) {
    // Apply the BC only on the upper half of the domain.
    int x_half_len = (update_rng.upper[0] - update_rng.lower[0] + 1)/2;
    gkyl_range_shorten_from_below(&update_rng, &update_rng, 0, x_half_len);
  }

  // Create the twist-shift updater and shift the donor field.
  struct gkyl_twistshift_dg_inp tsinp = {
    .bc_dir = bc_dir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = edge,
    .cdim = cdim,
    .bcdir_ext_update_r = &update_rng,
    .num_ghost = ghost,
    .basis = &basis,
    .grid = &grid,
    .shift_func = shift_fig11,
    .shift_func_ctx = &proj_ctx,
    .use_gpu = use_gpu,
  };

  struct gkyl_twistshift_dg *tsup = gkyl_twistshift_dg_inew(&tsinp);

  // First apply periodicity in z.
  struct gkyl_array *buff_per = mkarr(use_gpu, basis.num_basis, skin_rng.volume);
  apply_periodic_bc(buff_per, distf, bc_dir, skin_ghost);

  gkyl_twistshift_dg_advance(tsup, distf, distf);
  gkyl_array_copy(distf_ho, distf);

  // Write out the target in the extended range.
  if (write_f) {
    double lower_ext[ndim], upper_ext[ndim];
    int cells_ext[ndim];
    for (int d=0; d<ndim; d++) {
      double dx = (upper[d]-lower[d])/cells[d];
      lower_ext[d] = lower[d]-dx*ghost[d];
      upper_ext[d] = upper[d]+dx*ghost[d];
      cells_ext[d] = cells[d]+2*ghost[d];
    }
    struct gkyl_rect_grid grid_ext;
    gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);
    gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "ctest_bc_twistshift_3x_fig11_tar.gkyl");
  }

  if (check_distf) {
    // Check 0th, 1st, 2nd and 6th DG coeffs.
    const double f0[] =
    {
     -6.2583195868812408e-06,  7.9567311608137231e-10,  5.0509906078140035e-05,
      1.2343352957621956e-03,  5.3472486075120380e-04,
      1.8256130045163409e-04,  3.0009163112130789e-09, -1.3326694681034118e-05,
      3.7575576376431747e-03,  6.1716772244852208e-03,
      4.4292093668301044e-03, -7.7229340458416127e-05, -3.8275760062292437e-05,
      3.0421687066165188e-03,  2.9194555762249566e-02,
      3.7056203209553537e-02,  1.4705657602094287e-04,  2.0550274812940050e-08,
      1.0240094716951157e-03,  4.7863752496844179e-02,
      9.6491415121706434e-02,  4.8597593964605536e-03,  8.0126693914598593e-08,
     -4.3972043519208224e-04,  3.1180735654073199e-02,
      9.6491415121706420e-02,  3.1180735654073172e-02, -4.3972043519208387e-04,
      8.0126693914598646e-08,  4.8597593964605458e-03,
      3.7056203209553502e-02,  4.7863752496844159e-02,  1.0240094716951155e-03,
      2.0550274812939984e-08,  1.4705657602094392e-04,
      4.4292093668300931e-03,  2.9194555762249559e-02,  3.0421687066165266e-03,
     -3.8275760062292789e-05, -7.7229340458415897e-05,
      1.8256130045163365e-04,  6.1716772244852129e-03,  3.7575576376431765e-03,
     -1.3326694681033820e-05,  3.0009163112130839e-09,
     -6.2583195868812094e-06,  5.3472486075120294e-04,  1.2343352957621958e-03,
      5.0509906078140123e-05,  7.9567311608137107e-10,
    };
    const double f1[] =
    {
     -2.4097311206019044e-06,  7.6833036063548360e-11,  7.2344842311530572e-06,
      5.2914450025189400e-04,  4.1028939039880680e-04,
      1.6599629313756654e-04,  2.2461980415131899e-09, -4.1694729744125476e-05,
      7.8438224724348157e-04,  3.3439621031217184e-03,
      3.6040320392418893e-03, -8.7133864207591700e-05,  3.8447881347591303e-05,
     -8.6762333862857460e-04,  8.6496739552592294e-03,
      1.4945209682648789e-02,  4.5852653062261286e-04, -3.3068454717620671e-09,
     -6.2529388927217496e-04,  1.5882293198893030e-03,
      1.7368783554643934e-02,  2.3314728455433925e-03,  4.3253655108451250e-08,
      9.4609992614394520e-06, -1.1247642251296977e-02,
     -1.7368783554643944e-02,  1.1247642251296978e-02, -9.4609992614413358e-06,
     -4.3253655108451283e-08, -2.3314728455433894e-03,
     -1.4945209682648786e-02, -1.5882293198893067e-03,  6.2529388927217593e-04,
      3.3068454717620489e-09, -4.5852653062261221e-04,
     -3.6040320392418802e-03, -8.6496739552592190e-03,  8.6762333862857905e-04,
     -3.8447881347591479e-05,  8.7133864207591442e-05,
     -1.6599629313756624e-04, -3.3439621031217180e-03, -7.8438224724348212e-04,
      4.1694729744125585e-05, -2.2461980415131949e-09,
      2.4097311206018891e-06, -4.1028939039880621e-04, -5.2914450025189443e-04,
     -7.2344842311531317e-06, -7.6833036063549213e-11,
    };
    const double f2[] =
    {
      5.3070807682913396e-06,  1.3751826183769400e-09,  6.8926303153258719e-05,
      4.1882549426144460e-04, -4.9306025336561313e-04,
     -2.8467596971097078e-04, -5.0150358149508259e-09,  1.3960594126194340e-05,
      3.0887163471536865e-03, -2.8179959565330977e-03,
     -4.7580701893096058e-03,  1.1300292822992238e-04, -5.7287804480022883e-05,
      3.3338137165750199e-03,  1.3685413489847044e-03,
     -2.7243387375233196e-02, -5.6694093500983371e-04,  2.2529555439069578e-08,
      1.7831031494975383e-03,  2.6027202631190245e-02,
     -2.3565987663701497e-02, -5.9773227818891196e-03, -1.1962791428901063e-07,
     -4.8028129954664068e-04,  3.0023711373051848e-02,
      2.3565987663701948e-02, -3.0023711373051674e-02,  4.8028129954667072e-04,
      1.1962791428901214e-07,  5.9773227818891586e-03,
      2.7243387375233262e-02, -2.6027202631190141e-02, -1.7831031494975104e-03,
     -2.2529555439068645e-08,  5.6694093500985290e-04,
      4.7580701893096223e-03, -1.3685413489846329e-03, -3.3338137165750282e-03,
      5.7287804480022280e-05, -1.1300292822991837e-04,
      2.8467596971097305e-04,  2.8179959565330925e-03, -3.0887163471536604e-03,
     -1.3960594126190875e-05,  5.0150358149509500e-09,
     -5.3070807682908602e-06,  4.9306025336561432e-04, -4.1882549426144140e-04,
     -6.8926303153258244e-05, -1.3751826183769166e-09,
    };
    const double f6[] =
    {
      2.0154146681198872e-24, -3.5688851136225351e-27,  8.0173640710989959e-24,
     -3.7540344263165092e-21, -8.5473941658947867e-21,
      2.1098399411304440e-21, -3.2987430986791354e-26, -1.2786887990890647e-22,
      5.5774262155530719e-21, -5.2422242910241063e-20,
     -1.3146241530559257e-19, -2.7562828879312281e-22, -1.0734417010552807e-22,
     -4.3870259352886555e-20,  2.4194632944932309e-19,
      2.1137135832652420e-19,  1.9015676044847855e-21, -8.5272613545944742e-25,
     -4.3773733626837943e-20, -9.3076782806108862e-19,
     -2.2823024255915727e-19,  1.2271105194682399e-19, -1.3400884207613767e-24,
     -1.0649198357855089e-21,  3.4926449572105291e-20,
     -3.0499841362447923e-19, -1.8965105476591979e-20, -2.6902054933988774e-21,
      7.0876511626018948e-26,  2.9290027937202072e-20,
      8.7582274320145349e-19, -6.0508072579716399e-19, -2.2396479411516622e-20,
     -2.8336961187045548e-26,  1.0484792769983417e-20,
     -1.2817391312059265e-19, -1.2999883788145454e-19,  5.2404683918956681e-22,
      7.0852375629219591e-22, -1.6110515602546699e-24,
     -4.7728594396620859e-21, -5.5496727075010419e-20,  7.7181639308670794e-21,
     -3.5828835842015954e-22, -2.4162654496516419e-26,
     -1.2872700269810765e-22,  3.0405747283361938e-21, -1.5767712674352063e-20,
     -9.1281108005001305e-22,  6.4377759843551990e-27,
    };

    struct gkyl_range check_ghost_rng, check_other_ghost_rng;
    if (apply_in_half_x < 0) {
      // Applied the BC only on the lower half of the domain.
      int x_half_len = (ghost_rng.upper[0] - ghost_rng.lower[0] + 1)/2;
      gkyl_range_shorten_from_above(&check_ghost_rng, &ghost_rng, 0, x_half_len);
      gkyl_range_shorten_from_below(&check_other_ghost_rng, &ghost_rng, 0, x_half_len);
    }
    else if (apply_in_half_x > 0) {
      // Applied the BC only on the upper half of the domain.
      int x_half_len = (ghost_rng.upper[0] - ghost_rng.lower[0] + 1)/2;
      gkyl_range_shorten_from_below(&check_ghost_rng, &ghost_rng, 0, x_half_len);
      gkyl_range_shorten_from_above(&check_other_ghost_rng, &ghost_rng, 0, x_half_len);
    }
    else
      check_ghost_rng = ghost_rng;

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &check_ghost_rng);
    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(&check_ghost_rng, iter.idx);
      double *f_c = gkyl_array_fetch(distf_ho, linidx);
      int refidx = (iter.idx[0]-1)*cells[1] + iter.idx[1]-1;
      TEST_CHECK( gkyl_compare(f0[refidx], f_c[0], 1e-13) );
      TEST_CHECK( gkyl_compare(f1[refidx], f_c[1], 1e-13) );
      TEST_CHECK( gkyl_compare(f2[refidx], f_c[2], 1e-12) );
      TEST_CHECK( gkyl_compare(f6[refidx], f_c[6], 1e-12) );
    }

    if (apply_in_half_x != 0) {
      // Check that the other half is untouched. 
      int skin_idx[GKYL_MAX_DIM];
      gkyl_range_iter_init(&iter, &check_other_ghost_rng);
      while (gkyl_range_iter_next(&iter)) {
        long linidx = gkyl_range_idx(&check_other_ghost_rng, iter.idx);
        double *f_c = gkyl_array_fetch(distf_ho, linidx);
  
        for (int d=0; d<check_other_ghost_rng.ndim; d++)
          skin_idx[d] = iter.idx[d];
        if (edge == GKYL_LOWER_EDGE)
          skin_idx[bc_dir] = local.upper[bc_dir];
        else
          skin_idx[bc_dir] = local.lower[bc_dir];
  
        linidx = gkyl_range_idx(&local_ext, skin_idx);
        double *fskin_c = gkyl_array_fetch(distf_ho, linidx);
  
        for (int k=0; k<distf_ho->ncomp; k++)
          TEST_CHECK( gkyl_compare(fskin_c[k], f_c[k], 1e-15) );
      }
    }
  }

  gkyl_array_release(buff_per);
  test_bc_twistshift_array_meta_release(mt);
  gkyl_twistshift_dg_release(tsup);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(distf);

}

void
init_donor_3x_fig14(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

  struct test_bc_twistshift_ctx *pars = ctx;
  double Lx[3] = {
    pars->upper[0]-pars->lower[0],
    pars->upper[1]-pars->lower[1],
    pars->upper[2]-pars->lower[2],
  };

  double f_amplitude = 1.0;
  double f_floor = 1.0e-10;

  // Cube
  double rx2 = pow(x-Lx[0]/2,2);
  double ry2 = pow(y-Lx[1]/2,2);
  double rz2 = pow(z-Lx[2]/2,2);

  if (rx2 < pow(Lx[0]/4,2) && ry2 < pow(Lx[1]/4,2))
    fout[0] = f_amplitude;
  else
    fout[0] = f_floor;
}

void
shift_fig14(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];

  struct test_bc_twistshift_ctx *pars = ctx;
  double Lx[3] = {
    pars->upper[0]-pars->lower[0],
    pars->upper[1]-pars->lower[1],
    pars->upper[2]-pars->lower[2],
  };
  double dx[3] = {
    Lx[0]/pars->cells[0],
    Lx[1]/pars->cells[1],
    Lx[2]/pars->cells[2],
  };
  enum gkyl_edge_loc edge = pars->edge;

  fout[0] = -(x-0.5);

  if (edge == GKYL_UPPER_EDGE)
    fout[0] *= -1.0;
}

void
test_bc_twistshift_3x_fig14_wcells(const int *cells, enum gkyl_edge_loc edge,
  int apply_in_half_x, bool check_distf, bool use_gpu, bool write_f)
{
  double vt = 1.0; // Thermal speed.
  double mass = 1.0;
  double B0 = 1.0; // Magnetic field magnitude.
  int bc_dir = 2; // Direction in which to apply TS.

  int poly_order = 1;
  const double lower[] = {0.0, 0.0, 0.0};
  const double upper[] = {1.0, 1.0, 1.0};
  int vdim = 0;
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int cdim = ndim - vdim;

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d=0; d<cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = cells[d];
  }

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  // Ranges.
  int ghost_conf[cdim];
  for (int d=0; d<cdim; d++) ghost_conf[d] = 1;
  struct gkyl_range local_conf, local_ext_conf; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);

  int ghost[ndim];
  for (int d=0; d<cdim; d++) ghost[d] = ghost_conf[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);
  struct skin_ghost_ranges skin_ghost_conf; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost_conf, &local_ext_conf, ghost_conf);

  // Pick skin and ghost ranges based on 'edge'.
  struct gkyl_range skin_rng, ghost_rng;
  struct gkyl_range skin_rng_conf, ghost_rng_conf;
  if (edge == GKYL_LOWER_EDGE) {
    skin_rng = skin_ghost.upper_skin[bc_dir];
    ghost_rng = skin_ghost.lower_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.upper_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.lower_ghost[bc_dir];
  }
  else {
    skin_rng = skin_ghost.lower_skin[bc_dir];
    ghost_rng = skin_ghost.upper_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.lower_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.upper_ghost[bc_dir];
  }

  struct test_bc_twistshift_ctx proj_ctx = {
    .lower = {lower[0], lower[1], lower[2]},
    .upper = {upper[0], upper[1], upper[2]},
    .cells = {cells[0], cells[1], cells[2]},
    .B0 = B0,
    .vt = vt,
    .mass = mass,
    .edge = edge,
  };

  // Initialize the distribution
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, basis.num_basis, local_ext.volume) : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .num_ret_vals = 1,
      .eval = init_donor_3x_fig14,
      .ctx = &proj_ctx,
    }
  );
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);
  struct gkyl_msgpack_data *mt = test_bc_twistshift_array_meta_new( (struct test_bc_twistshift_output_meta) {
      .poly_order = poly_order,
      .basis_type = basis.id
    }
  );
  if (write_f)
    gkyl_grid_sub_array_write(&grid, &local, mt, distf_ho, "ctest_bc_twistshift_3x_fig14_do.gkyl");

  // Create a range only extended in bc_dir.
  struct gkyl_range update_rng;
  int lower_bcdir_ext[ndim], upper_bcdir_ext[ndim];
  for (int d=0; d<ndim; d++) {
    lower_bcdir_ext[d] = local.lower[d];
    upper_bcdir_ext[d] = local.upper[d];
  }
  lower_bcdir_ext[bc_dir] = local.lower[bc_dir] - ghost[bc_dir];
  upper_bcdir_ext[bc_dir] = local.upper[bc_dir] + ghost[bc_dir];
  gkyl_sub_range_init(&update_rng, &local_ext, lower_bcdir_ext, upper_bcdir_ext);

  if (apply_in_half_x < 0) {
    // Apply the BC only on the lower half of the domain.
    int x_half_len = (update_rng.upper[0] - update_rng.lower[0] + 1)/2;
    gkyl_range_shorten_from_above(&update_rng, &update_rng, 0, x_half_len);
  }
  else if (apply_in_half_x > 0) {
    // Apply the BC only on the upper half of the domain.
    int x_half_len = (update_rng.upper[0] - update_rng.lower[0] + 1)/2;
    gkyl_range_shorten_from_below(&update_rng, &update_rng, 0, x_half_len);
  }

  // Create the twist-shift updater and shift the donor field.
  struct gkyl_twistshift_dg_inp tsinp = {
    .bc_dir = bc_dir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = edge,
    .cdim = cdim,
    .bcdir_ext_update_r = &update_rng,
    .num_ghost = ghost,
    .basis = &basis,
    .grid = &grid,
    .shift_func = shift_fig14,
    .shift_func_ctx = &proj_ctx,
    .use_gpu = use_gpu,
  };

  struct gkyl_twistshift_dg *tsup = gkyl_twistshift_dg_inew(&tsinp);

  // First apply periodicity in z.
  struct gkyl_array *buff_per = mkarr(use_gpu, basis.num_basis, skin_rng.volume);
  apply_periodic_bc(buff_per, distf, bc_dir, skin_ghost);

  gkyl_twistshift_dg_advance(tsup, distf, distf);
  gkyl_array_copy(distf_ho, distf);

  // Write out the target in the extended range.
  if (write_f) {
    double lower_ext[ndim], upper_ext[ndim];
    int cells_ext[ndim];
    for (int d=0; d<ndim; d++) {
      double dx = (upper[d]-lower[d])/cells[d];
      lower_ext[d] = lower[d]-dx*ghost[d];
      upper_ext[d] = upper[d]+dx*ghost[d];
      cells_ext[d] = cells[d]+2*ghost[d];
    }
    struct gkyl_rect_grid grid_ext;
    gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);
    gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "ctest_bc_twistshift_3x_fig14_tar.gkyl");
  }

  if (check_distf) {
    // Check 0th, 1st, 2nd and 6th DG coeffs.
    double f0[64], f1[64], f2[64], f6[64];

    if (edge == GKYL_UPPER_EDGE) {
      const double f0m[] =
      {
        2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        1.4142135625145160e+00, 2.8284271247461894e+00,
        2.8284271247461894e+00, 2.8284271247461894e+00, 1.4142135625145160e+00,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 1.4142135625145160e+00,
        2.8284271247461894e+00, 2.8284271247461894e+00, 2.8284271247461894e+00,
        1.4142135625145160e+00, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10,
        1.4142135625145160e+00, 2.8284271247461894e+00, 2.8284271247461894e+00,
        2.8284271247461894e+00, 1.4142135625145160e+00, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 1.4142135625145160e+00, 2.8284271247461894e+00,
        2.8284271247461894e+00, 2.8284271247461894e+00, 1.4142135625145160e+00,
        2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
      };
      const double f1m[] =
      {
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26,
        -8.1649658084607613e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,
         0.0000000000000000e+00,  8.1649658084607613e-01, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -8.1649658084607613e-01,  0.0000000000000000e+00,
         0.0000000000000000e+00,  0.0000000000000000e+00,  8.1649658084607613e-01,
        -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -8.1649658084607613e-01,
         0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
         8.1649658084607613e-01, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -8.1649658084607613e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,
         0.0000000000000000e+00,  8.1649658084607613e-01,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26,
      };
      const double f2m[] =
      {
        -5.2991257991678335e-25, -5.2991257991678335e-25, -5.2991257991678335e-25,
        -5.2991257991678335e-25, -5.2991257991678335e-25, -5.2991257991678335e-25,
        -5.2991257991678335e-25, -5.2991257991678335e-25,
        -3.4896682092080855e-25, -3.4896682092080855e-25, -3.4896682092080855e-25,
        -3.4896682092080855e-25, -3.4896682092080855e-25, -3.4896682092080855e-25,
        -3.4896682092080855e-25, -3.4896682092080855e-25,
         8.1649658084607457e-01, -2.5535129566378600e-15, -2.5535129566378600e-15,
        -2.5535129566378600e-15, -8.1649658084607712e-01, -2.5849394142282115e-25,
        -2.5849394142282115e-25, -2.5849394142282115e-25,
        -6.4623485355705287e-26,  8.1649658084607535e-01, -6.6613381477509392e-16,
        -6.6613381477509392e-16, -6.6613381477509392e-16, -8.1649658084607601e-01,
        -6.4623485355705287e-26, -6.4623485355705287e-26,
         6.4623485355705287e-26,  6.4623485355705287e-26,  8.1649658084607601e-01,
         6.6613381477509392e-16,  6.6613381477509392e-16,  6.6613381477509392e-16,
        -8.1649658084607535e-01,  6.4623485355705287e-26,
         2.5849394142282115e-25,  2.5849394142282115e-25,  2.5849394142282115e-25,
         8.1649658084607712e-01,  2.5535129566378600e-15,  2.5535129566378600e-15,
         2.5535129566378600e-15, -8.1649658084607457e-01,
         5.2991257991678335e-25,  5.2991257991678335e-25,  5.2991257991678335e-25,
         5.2991257991678335e-25,  5.2991257991678335e-25,  5.2991257991678335e-25,
         5.2991257991678335e-25,  5.2991257991678335e-25,
         5.2991257991678335e-25,  5.2991257991678335e-25,  5.2991257991678335e-25,
         5.2991257991678335e-25,  5.2991257991678335e-25,  5.2991257991678335e-25,
         5.2991257991678335e-25,  5.2991257991678335e-25,
      };
      const double f6m[] =
      {
         1.2712579668354740e-41,  1.2712579668354740e-41,  1.2712579668354740e-41,
         1.2712579668354740e-41,  1.2712579668354740e-41,  1.2712579668354740e-41,
         1.2712579668354740e-41,  1.2712579668354740e-41,
        -3.1389085600875902e-43, -3.1389085600875902e-43, -3.1389085600875902e-43,
        -3.1389085600875902e-43, -3.1389085600875902e-43, -3.1389085600875902e-43,
        -3.1389085600875902e-43, -3.1389085600875902e-43,
        -1.1538233466087235e-28,  0.0000000000000000e+00,  0.0000000000000000e+00,
         0.0000000000000000e+00,  1.1538233466087551e-28,  3.1613293355167873e-42,
         3.1613293355167873e-42,  3.1613293355167873e-42,
         3.2061708863751815e-42, -1.1538233466087102e-28,  0.0000000000000000e+00,
         0.0000000000000000e+00,  0.0000000000000000e+00,  1.1538233466087423e-28,
         3.2061708863751815e-42,  3.2061708863751815e-42,
         1.9730282377693424e-42,  1.9730282377693424e-42, -1.1538233466087044e-28,
         0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
         1.1538233466087241e-28,  1.9730282377693424e-42,
         9.8651411888467122e-43,  9.8651411888467122e-43,  9.8651411888467122e-43,
        -1.1538233466086849e-28,  0.0000000000000000e+00,  0.0000000000000000e+00,
         0.0000000000000000e+00,  1.1538233466086948e-28,
         1.9057659114817512e-42,  1.9057659114817512e-42,  1.9057659114817512e-42,
         1.9057659114817512e-42,  1.9057659114817512e-42,  1.9057659114817512e-42,
         1.9057659114817512e-42,  1.9057659114817512e-42,
        -1.5470335046145980e-42, -1.5470335046145980e-42, -1.5470335046145980e-42,
        -1.5470335046145980e-42, -1.5470335046145980e-42, -1.5470335046145980e-42,
        -1.5470335046145980e-42, -1.5470335046145980e-42,
      };
      int nvals = sizeof(f0m)/sizeof(f0m[0]);
      for (int i=0; i<nvals; i++) {
        f0[i] = f0m[i];  f1[i] = f1m[i];  f2[i] = f2m[i]; f6[i] = f6m[i];
      }
    }
    else {
      const double f0m[] =
      {
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        1.4142135625145160e+00, 2.8284271247461894e+00, 2.8284271247461894e+00,
        2.8284271247461894e+00, 1.4142135625145160e+00,
        2.8284271247461891e-10, 2.8284271247461891e-10, 1.4142135625145160e+00,
        2.8284271247461894e+00, 2.8284271247461894e+00, 2.8284271247461894e+00,
        1.4142135625145160e+00, 2.8284271247461891e-10,
        2.8284271247461891e-10, 1.4142135625145160e+00, 2.8284271247461894e+00,
        2.8284271247461894e+00, 2.8284271247461894e+00, 1.4142135625145160e+00,
        2.8284271247461891e-10, 2.8284271247461891e-10,
        1.4142135625145160e+00, 2.8284271247461894e+00, 2.8284271247461894e+00,
        2.8284271247461894e+00, 1.4142135625145160e+00, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10, 2.8284271247461891e-10,
        2.8284271247461891e-10, 2.8284271247461891e-10,
      };
      const double f1m[] =
      {
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
         8.1649658084607613e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,
         0.0000000000000000e+00, -8.1649658084607613e-01,
        -2.5849394142282115e-26, -2.5849394142282115e-26,  8.1649658084607613e-01,
         0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
        -8.1649658084607613e-01, -2.5849394142282115e-26,
        -2.5849394142282115e-26,  8.1649658084607613e-01,  0.0000000000000000e+00,
         0.0000000000000000e+00,  0.0000000000000000e+00, -8.1649658084607613e-01,
        -2.5849394142282115e-26, -2.5849394142282115e-26,
         8.1649658084607613e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,
         0.0000000000000000e+00, -8.1649658084607613e-01, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26, -2.5849394142282115e-26,
        -2.5849394142282115e-26, -2.5849394142282115e-26,
      };
      const double f2m[] =
      {
         5.2991257991678335e-25,  5.2991257991678335e-25,  5.2991257991678335e-25,
         5.2991257991678335e-25,  5.2991257991678335e-25,  5.2991257991678335e-25,
         5.2991257991678335e-25,  5.2991257991678335e-25,
         3.4896682092080855e-25,  3.4896682092080855e-25,  3.4896682092080855e-25,
         3.4896682092080855e-25,  3.4896682092080855e-25,  3.4896682092080855e-25,
         3.4896682092080855e-25,  3.4896682092080855e-25,
         2.5849394142282115e-25,  2.5849394142282115e-25,  2.5849394142282115e-25,
         8.1649658084607712e-01,  2.5535129566378600e-15,  2.5535129566378600e-15,
         2.5535129566378600e-15, -8.1649658084607457e-01,
         6.4623485355705287e-26,  6.4623485355705287e-26,  8.1649658084607601e-01,
         6.6613381477509392e-16,  6.6613381477509392e-16,  6.6613381477509392e-16,
        -8.1649658084607535e-01,  6.4623485355705287e-26,
        -6.4623485355705287e-26,  8.1649658084607535e-01, -6.6613381477509392e-16,
        -6.6613381477509392e-16, -6.6613381477509392e-16, -8.1649658084607601e-01,
        -6.4623485355705287e-26, -6.4623485355705287e-26,
         8.1649658084607457e-01, -2.5535129566378600e-15, -2.5535129566378600e-15,
        -2.5535129566378600e-15, -8.1649658084607712e-01, -2.5849394142282115e-25,
        -2.5849394142282115e-25, -2.5849394142282115e-25,
        -5.2991257991678335e-25, -5.2991257991678335e-25, -5.2991257991678335e-25,
        -5.2991257991678335e-25, -5.2991257991678335e-25, -5.2991257991678335e-25,
        -5.2991257991678335e-25, -5.2991257991678335e-25,
        -5.2991257991678335e-25, -5.2991257991678335e-25, -5.2991257991678335e-25,
        -5.2991257991678335e-25, -5.2991257991678335e-25, -5.2991257991678335e-25,
        -5.2991257991678335e-25, -5.2991257991678335e-25,
      };
      const double f6m[] =
      {
        -5.3809861030072976e-43, -5.3809861030072976e-43, -5.3809861030072976e-43,
        -5.3809861030072976e-43, -5.3809861030072976e-43, -5.3809861030072976e-43,
        -5.3809861030072976e-43, -5.3809861030072976e-43,
        -1.0941338409448172e-41, -1.0941338409448172e-41, -1.0941338409448172e-41,
        -1.0941338409448172e-41, -1.0941338409448172e-41, -1.0941338409448172e-41,
        -1.0941338409448172e-41, -1.0941338409448172e-41,
        -6.0984509167416039e-42, -6.0984509167416039e-42, -6.0984509167416039e-42,
         1.2692056812695856e-27,  0.0000000000000000e+00,  0.0000000000000000e+00,
         0.0000000000000000e+00, -1.2692056812695917e-27,
        -1.7936620343357659e-42, -1.7936620343357659e-42,  1.2692056812695897e-27,
         0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
        -1.2692056812695915e-27, -1.7936620343357659e-42,
        -1.2555634240350361e-42,  1.2692056812695930e-27,  0.0000000000000000e+00,
         0.0000000000000000e+00,  0.0000000000000000e+00, -1.2692056812695942e-27,
        -1.2555634240350361e-42, -1.2555634240350361e-42,
         1.2692056812695946e-27,  0.0000000000000000e+00,  0.0000000000000000e+00,
         0.0000000000000000e+00, -1.2692056812695955e-27, -8.9683101716788293e-43,
        -8.9683101716788293e-43, -8.9683101716788293e-43,
         3.7666902721051083e-42,  3.7666902721051083e-42,  3.7666902721051083e-42,
         3.7666902721051083e-42,  3.7666902721051083e-42,  3.7666902721051083e-42,
         3.7666902721051083e-42,  3.7666902721051083e-42,
         5.3809861030072976e-43,  5.3809861030072976e-43,  5.3809861030072976e-43,
         5.3809861030072976e-43,  5.3809861030072976e-43,  5.3809861030072976e-43,
         5.3809861030072976e-43,  5.3809861030072976e-43,
      };
      int nvals = sizeof(f0m)/sizeof(f0m[0]);
      for (int i=0; i<nvals; i++) {
        f0[i] = f0m[i];  f1[i] = f1m[i];  f2[i] = f2m[i]; f6[i] = f6m[i];
      }
    }

    struct gkyl_range check_ghost_rng, check_other_ghost_rng;
    if (apply_in_half_x < 0) {
      // Applied the BC only on the lower half of the domain.
      int x_half_len = (ghost_rng.upper[0] - ghost_rng.lower[0] + 1)/2;
      gkyl_range_shorten_from_above(&check_ghost_rng, &ghost_rng, 0, x_half_len);
      gkyl_range_shorten_from_below(&check_other_ghost_rng, &ghost_rng, 0, x_half_len);
    }
    else if (apply_in_half_x > 0) {
      // Applied the BC only on the upper half of the domain.
      int x_half_len = (ghost_rng.upper[0] - ghost_rng.lower[0] + 1)/2;
      gkyl_range_shorten_from_below(&check_ghost_rng, &ghost_rng, 0, x_half_len);
      gkyl_range_shorten_from_above(&check_other_ghost_rng, &ghost_rng, 0, x_half_len);
    }
    else
      check_ghost_rng = ghost_rng;

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &check_ghost_rng);
    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(&check_ghost_rng, iter.idx);
      double *f_c = gkyl_array_fetch(distf_ho, linidx);
      int refidx = (iter.idx[0]-1)*cells[1] + iter.idx[1]-1;
      TEST_CHECK( gkyl_compare(f0[refidx], f_c[0], 1e-13) );
      TEST_CHECK( gkyl_compare(f1[refidx], f_c[1], 1e-13) );
      TEST_CHECK( gkyl_compare(f2[refidx], f_c[2], 1e-12) );
      TEST_CHECK( gkyl_compare(f6[refidx], f_c[6], 1e-12) );
    }

  }

  gkyl_array_release(buff_per);
  test_bc_twistshift_array_meta_release(mt);
  gkyl_twistshift_dg_release(tsup);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(distf);

}

void
test_bc_twistshift_3x2v_fig11_wcells(const int *cells, enum gkyl_edge_loc edge,
  int apply_in_half_x, bool check_distf, bool use_gpu, bool write_f)
{
  double vt = 1.0; // Thermal speed.
  double mass = 1.0;
  double B0 = 1.0; // Magnetic field magnitude.
  int bc_dir = 2; // Direction in which to apply TS.

  int poly_order = 1;
  const double lower[] = {-2.0, -1.50, -3.0, -5.0*vt, 0.};
  const double upper[] = { 2.0,  1.50,  3.0,  5.0*vt, mass*(pow(5.0*vt,2))/(2.0*B0)};
  int vdim = 2;
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int cdim = ndim - vdim;

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d=0; d<cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = cells[d];
  }
  double lower_vel[vdim], upper_vel[vdim];
  int cells_vel[vdim];
  for (int d=0; d<vdim; d++) {
    lower_vel[d] = lower[cdim+d];
    upper_vel[d] = upper[cdim+d];
    cells_vel[d] = cells[cdim+d];
  }

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  struct gkyl_rect_grid grid_vel;
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

  // Basis functions.
  struct gkyl_basis basis;
  if (poly_order == 1) 
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  // Ranges.
  int ghost_conf[cdim];
  for (int d=0; d<cdim; d++) ghost_conf[d] = 1;
  struct gkyl_range local_conf, local_ext_conf; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);

  int ghost_vel[vdim];
  for (int d=0; d<vdim; d++) ghost_vel[d] = 0;
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_ext_vel, &local_vel);

  int ghost[ndim];
  for (int d=0; d<cdim; d++) ghost[d] = ghost_conf[d];
  for (int d=cdim; d<ndim; d++) ghost[d] = ghost_vel[d-cdim];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);
  struct skin_ghost_ranges skin_ghost_conf; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost_conf, &local_ext_conf, ghost_conf);

  // Pick skin and ghost ranges based on 'edge'.
  struct gkyl_range skin_rng, ghost_rng;
  struct gkyl_range skin_rng_conf, ghost_rng_conf;
  if (edge == GKYL_LOWER_EDGE) {
    skin_rng = skin_ghost.upper_skin[bc_dir];
    ghost_rng = skin_ghost.lower_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.upper_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.lower_ghost[bc_dir];
  }
  else {
    skin_rng = skin_ghost.lower_skin[bc_dir];
    ghost_rng = skin_ghost.upper_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.lower_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.upper_ghost[bc_dir];
  }

  struct test_bc_twistshift_ctx proj_ctx = {
    .lower = {lower[0], lower[1], lower[2], lower[3], lower[4]},
    .upper = {upper[0], upper[1], upper[2], upper[3], upper[4]},
    .cells = {cells[0], cells[1], cells[2], cells[3], cells[4]},
    .B0 = B0,
    .vt = vt,
    .mass = mass,
  };

  // Initialize the distribution
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, basis.num_basis, local_ext.volume) : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .num_ret_vals = 1,
      .eval = init_donor_3x2v_fig11,
      .ctx = &proj_ctx,
    }
  );
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);
  struct gkyl_msgpack_data *mt = test_bc_twistshift_array_meta_new( (struct test_bc_twistshift_output_meta) {
      .poly_order = poly_order,
      .basis_type = basis.id
    }
  );
  if (write_f)
    gkyl_grid_sub_array_write(&grid, &local, mt, distf_ho, "ctest_bc_twistshift_3x2v_fig11_do.gkyl");

  // Create a range only extended in bc_dir.
  struct gkyl_range update_rng;
  int lower_bcdir_ext[ndim], upper_bcdir_ext[ndim];
  for (int d=0; d<ndim; d++) {
    lower_bcdir_ext[d] = local.lower[d];
    upper_bcdir_ext[d] = local.upper[d];
  }
  lower_bcdir_ext[bc_dir] = local.lower[bc_dir] - ghost[bc_dir];
  upper_bcdir_ext[bc_dir] = local.upper[bc_dir] + ghost[bc_dir];
  gkyl_sub_range_init(&update_rng, &local_ext, lower_bcdir_ext, upper_bcdir_ext);

  if (apply_in_half_x < 0) {
    // Apply the BC only on the lower half of the domain.
    int x_half_len = (update_rng.upper[0] - update_rng.lower[0] + 1)/2;
    gkyl_range_shorten_from_above(&update_rng, &update_rng, 0, x_half_len);
  }
  else if (apply_in_half_x > 0) {
    // Apply the BC only on the upper half of the domain.
    int x_half_len = (update_rng.upper[0] - update_rng.lower[0] + 1)/2;
    gkyl_range_shorten_from_below(&update_rng, &update_rng, 0, x_half_len);
  }

  // Create the twist-shift updater and shift the donor field.
  struct gkyl_twistshift_dg_inp tsinp = {
    .bc_dir = bc_dir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = edge,
    .cdim = cdim,
    .bcdir_ext_update_r = &update_rng,
    .num_ghost = ghost,
    .basis = &basis,
    .grid = &grid,
    .shift_func = shift_fig11,
    .shift_func_ctx = &proj_ctx,
    .use_gpu = use_gpu,
  };

  struct gkyl_twistshift_dg *tsup = gkyl_twistshift_dg_inew(&tsinp);

  // First apply periodicity in z.
  struct gkyl_array *buff_per = mkarr(use_gpu, basis.num_basis, skin_rng.volume);
  apply_periodic_bc(buff_per, distf, bc_dir, skin_ghost);

  gkyl_twistshift_dg_advance(tsup, distf, distf);
  gkyl_array_copy(distf_ho, distf);

  // Write out the target in the extended range.
  if (write_f) {
    double lower_ext[ndim], upper_ext[ndim];
    int cells_ext[ndim];
    for (int d=0; d<ndim; d++) {
      double dx = (upper[d]-lower[d])/cells[d];
      lower_ext[d] = lower[d]-dx*ghost[d];
      upper_ext[d] = upper[d]+dx*ghost[d];
      cells_ext[d] = cells[d]+2*ghost[d];
    }
    struct gkyl_rect_grid grid_ext;
    gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);
    gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "ctest_bc_twistshift_3x2v_fig11_tar.gkyl");
  }

  // Compute the integrated moments of the skin cell and the ghost cell.
  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_ext_vel, use_gpu);
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .c2p_ctx = 0,
    .mapc2p = mapc2p,
    .bfield_ctx = &proj_ctx,
    .bfield_func = eval_bfield_3x,
    .position_map = pmap,
    .grid = grid_conf,
    .local = local_conf,
    .local_ext = local_ext_conf,
    .global = local_conf,
    .global_ext = local_ext_conf,
    .basis = basis_conf,
    .geo_grid = grid_conf,
    .geo_local = local_conf,
    .geo_local_ext = local_ext_conf,
    .geo_global = local_conf,
    .geo_global_ext = local_ext_conf,
    .geo_basis = basis_conf,
  };
  struct gk_geometry* gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_inp);
  struct gk_geometry* gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);
  gkyl_gk_geometry_release(gk_geom_3d); // release temporary 3d geometry
  if (use_gpu) {  // If we are on the gpu, copy from host
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_inp, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }
  // Need the magnetic field to be initialized in the ghost cell in order to
  // compute integrated moments in the ghost cell (for checking moment
  // conservation).
  gkyl_array_clear(gk_geom->geo_corn.bmag, 0.0);
  gkyl_array_shiftc(gk_geom->geo_corn.bmag, B0*pow(sqrt(2.0),cdim), 0);

  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_gyrokinetic_new(&grid, &basis_conf,
    &basis, &local_conf, mass, 0, gvm, gk_geom, NULL, GKYL_F_MOMENT_M0M1M2, true, use_gpu);
  int num_mom = gkyl_dg_updater_moment_gyrokinetic_num_mom(mcalc);

  struct gkyl_array *marr = mkarr(use_gpu, num_mom, local_ext_conf.volume);
  double *red_integ_mom_skin, *red_integ_mom_ghost;
  if (use_gpu) {
    red_integ_mom_skin = gkyl_cu_malloc(sizeof(double[num_mom]));
    red_integ_mom_ghost = gkyl_cu_malloc(sizeof(double[num_mom]));
  }
  else {
    red_integ_mom_skin = gkyl_malloc(sizeof(double[num_mom]));
    red_integ_mom_ghost = gkyl_malloc(sizeof(double[num_mom]));
  }
  double *red_integ_mom_skin_ho = gkyl_malloc(sizeof(double[num_mom]));
  double *red_integ_mom_ghost_ho = gkyl_malloc(sizeof(double[num_mom]));

  gkyl_dg_updater_moment_gyrokinetic_advance(mcalc,
      &skin_rng, &skin_rng_conf, distf, marr);
  gkyl_array_reduce_range(red_integ_mom_skin, marr, GKYL_SUM, &skin_rng_conf);

  gkyl_dg_updater_moment_gyrokinetic_advance(mcalc,
      &ghost_rng, &ghost_rng_conf, distf, marr);
  gkyl_array_reduce_range(red_integ_mom_ghost, marr, GKYL_SUM, &ghost_rng_conf);

  if (use_gpu) {
    gkyl_cu_memcpy(red_integ_mom_skin_ho, red_integ_mom_skin, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(red_integ_mom_ghost_ho, red_integ_mom_ghost, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    memcpy(red_integ_mom_skin_ho, red_integ_mom_skin, sizeof(double[num_mom]));
    memcpy(red_integ_mom_ghost_ho, red_integ_mom_ghost, sizeof(double[num_mom]));
  }

  for (int k=0; k<num_mom; k++) {
    TEST_CHECK( gkyl_compare(red_integ_mom_skin_ho[k], red_integ_mom_ghost_ho[k], 1e-12));
    TEST_MSG( "integ_mom %d | Expected: %.14e | Got: %.14e\n",k,red_integ_mom_skin_ho[k],red_integ_mom_ghost_ho[k]);
  }

  if (check_distf) {
    // Check 0th, 1st, 2nd and 6th DG coeffs.
    const double f0[] =
    {
      -1.2776583204321330e-07,  1.6243951159613727e-11,  1.0311777925221165e-06,
       2.5199396402501986e-05,  1.0916599224379850e-05,
       3.7270542239467049e-06,  6.1264779478169708e-11, -2.7206923690464824e-07,
       7.6711882696143613e-05,  1.2599699723571620e-04,
       9.0423892898165817e-05, -1.5766645989926764e-06, -7.8141332726062177e-07,
       6.2106961880222910e-05,  5.9601729446907966e-04,
       7.5651563805650783e-04,  3.0022125798956283e-06,  4.1954120810574047e-10,
       2.0905519501674566e-05,  9.7715562102146359e-04,
       1.9699067404442962e-03,  9.9213725697257197e-05,  1.6358151057562534e-09,
      -8.9770499075323632e-06,  6.3656586712399692e-04,
       1.9699067404442966e-03,  6.3656586712399627e-04, -8.9770499075323903e-06,
       1.6358151057562547e-09,  9.9213725697256980e-05,
       7.5651563805650697e-04,  9.7715562102146315e-04,  2.0905519501674579e-05,
       4.1954120810574011e-10,  3.0022125798956054e-06,
       9.0423892898165560e-05,  5.9601729446907923e-04,  6.2106961880223127e-05,
      -7.8141332726062940e-07, -1.5766645989926737e-06,
       3.7270542239466922e-06,  1.2599699723571601e-04,  7.6711882696143640e-05,
      -2.7206923690464432e-07,  6.1264779478169876e-11,
      -1.2776583204321266e-07,  1.0916599224379833e-05,  2.5199396402501990e-05,
       1.0311777925221199e-06,  1.6243951159613701e-11,
    };
    const double f1[] =
    {
      -4.9195522432173540e-08,  1.5685739030718726e-12,  1.4769458228600922e-07,
       1.0802674169515287e-05,  8.3762046049376821e-06,
       3.3888736767723010e-06,  4.5856936150936589e-11, -8.5121281577597483e-07,
       1.6013444035211209e-05,  6.8268181976822328e-05,
       7.3577602711338998e-05, -1.7788689927180866e-06,  7.8492724484237293e-07,
      -1.7712840678887841e-05,  1.7658618650808429e-04,
       3.0511179936649428e-04,  9.3609830699107460e-06, -6.7510432676475386e-11,
      -1.2765598324806747e-05,  3.2424269440705722e-05,
       3.5458992651786952e-04,  4.7597851765421919e-05,  8.8303883448615030e-10,
       1.9314968272504524e-07, -2.2962463818142971e-04,
      -3.5458992651786996e-04,  2.2962463818142974e-04, -1.9314968272508097e-07,
      -8.8303883448615154e-10, -4.7597851765421851e-05,
      -3.0511179936649438e-04, -3.2424269440705817e-05,  1.2765598324806735e-05,
       6.7510432676475619e-11, -9.3609830699107308e-06,
      -7.3577602711338781e-05, -1.7658618650808405e-04,  1.7712840678887892e-05,
      -7.8492724484237631e-07,  1.7788689927180847e-06,
      -3.3888736767722963e-06, -6.8268181976822274e-05, -1.6013444035211209e-05,
       8.5121281577597260e-07, -4.5856936150936738e-11,
       4.9195522432173143e-08, -8.3762046049376720e-06, -1.0802674169515287e-05,
      -1.4769458228601091e-07, -1.5685739030718849e-12,
    };
    const double f2[] =
    {
       1.0834595144400169e-07,  2.8074844853976726e-11,  1.4071551236371810e-06,
       8.5504722173976439e-06, -1.0066001367323685e-05,
      -5.8117617082185942e-06, -1.0238374930018830e-10,  2.8501052072984633e-07,
       6.3057248605008680e-05, -5.7530395033770669e-05,
      -9.7137704174052986e-05,  2.3069951842792511e-06, -1.1695510118502272e-06,
       6.8060999036891783e-05,  2.7939260964732569e-05,
      -5.5618349420324920e-04, -1.1574301899303677e-05,  4.5994892978466414e-10,
       3.6402688349758515e-05,  5.3135464780386854e-04,
      -4.8110806422932840e-04, -1.2202918179822416e-04, -2.4422466434550431e-09,
      -9.8051144559640598e-06,  6.1294480273016622e-04,
       4.8110806422933762e-04, -6.1294480273016254e-04,  9.8051144559646528e-06,
       2.4422466434550737e-09,  1.2202918179822499e-04,
       5.5618349420325061e-04, -5.3135464780386637e-04, -3.6402688349757959e-05,
      -4.5994892978464471e-10,  1.1574301899304035e-05,
       9.7137704174053312e-05, -2.7939260964731105e-05, -6.8060999036891905e-05,
       1.1695510118502157e-06, -2.3069951842791736e-06,
       5.8117617082186383e-06,  5.7530395033770560e-05, -6.3057248605008110e-05,
      -2.8501052072977417e-07,  1.0238374930019089e-10,
      -1.0834595144399179e-07,  1.0066001367323711e-05, -8.5504722173975829e-06,
      -1.4071551236371706e-06, -2.8074844853976218e-11,
    };
    const double f6[] =
    {
      -1.6760873656052127e-08,  6.9405324157913410e-12,  3.9141599258817577e-07,
       6.6551168467387512e-06, -7.0297789062032620e-06,
      -4.5612495563369940e-06, -1.0697208933100729e-10, -1.2131836499023975e-06,
       2.2469828589386375e-05, -1.6695288411057533e-05,
      -7.7165270659134859e-05,  2.5304169488836587e-06,  1.1365688119276509e-06,
      -1.6574594816272132e-05,  9.0072879714596268e-05,
      -1.3067050406710185e-04, -1.7045882735192163e-05, -5.6876686935148544e-10,
      -1.3656577294632045e-05,  1.6137353286379611e-04,
       2.1258539335022546e-04, -4.1000351637691106e-05, -8.2479310847183273e-10,
      -3.5225958304071601e-06, -1.6806162108901716e-04,
       2.1258539335022143e-04, -1.6806162108901589e-04, -3.5225958304067087e-06,
      -8.2479310847182528e-10, -4.1000351637691127e-05,
      -1.3067050406710394e-04,  1.6137353286379722e-04, -1.3656577294632221e-05,
      -5.6876686935149423e-10, -1.7045882735192441e-05,
      -7.7165270659135103e-05,  9.0072879714595211e-05, -1.6574594816271881e-05,
       1.1365688119277396e-06,  2.5304169488835875e-06,
      -4.5612495563370270e-06, -1.6695288411057713e-05,  2.2469828589386263e-05,
      -1.2131836499024202e-06, -1.0697208933100851e-10,
      -1.6760873656069723e-08, -7.0297789062032917e-06,  6.6551168467386817e-06,
       3.9141599258816127e-07,  6.9405324157907182e-12,
    };

    struct gkyl_range check_ghost_rng, check_other_ghost_rng;
    if (apply_in_half_x < 0) {
      // Applied the BC only on the lower half of the domain.
      int x_half_len = (ghost_rng.upper[0] - ghost_rng.lower[0] + 1)/2;
      gkyl_range_shorten_from_above(&check_ghost_rng, &ghost_rng, 0, x_half_len);
      gkyl_range_shorten_from_below(&check_other_ghost_rng, &ghost_rng, 0, x_half_len);
    }
    else if (apply_in_half_x > 0) {
      // Applied the BC only on the upper half of the domain.
      int x_half_len = (ghost_rng.upper[0] - ghost_rng.lower[0] + 1)/2;
      gkyl_range_shorten_from_below(&check_ghost_rng, &ghost_rng, 0, x_half_len);
      gkyl_range_shorten_from_above(&check_other_ghost_rng, &ghost_rng, 0, x_half_len);
    }
    else
      check_ghost_rng = ghost_rng;

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &check_ghost_rng);
    while (gkyl_range_iter_next(&iter)) {
      if (iter.idx[3]==1 && iter.idx[4]==1) {
        long linidx = gkyl_range_idx(&check_ghost_rng, iter.idx);
        double *f_c = gkyl_array_fetch(distf_ho, linidx);
        int refidx = (iter.idx[0]-1)*cells[1] + iter.idx[1]-1;
        TEST_CHECK( gkyl_compare(f0[refidx], f_c[0], 1e-13) );
        TEST_CHECK( gkyl_compare(f1[refidx], f_c[1], 1e-13) );
        TEST_CHECK( gkyl_compare(f2[refidx], f_c[2], 1e-12) );
        TEST_CHECK( gkyl_compare(f6[refidx], f_c[6], 1e-12) );
      }
    }

    if (apply_in_half_x != 0) {
      // Check that the other half is untouched. 
      int skin_idx[GKYL_MAX_DIM];
      gkyl_range_iter_init(&iter, &check_other_ghost_rng);
      while (gkyl_range_iter_next(&iter)) {
        long linidx = gkyl_range_idx(&check_other_ghost_rng, iter.idx);
        double *f_c = gkyl_array_fetch(distf_ho, linidx);
  
        for (int d=0; d<check_other_ghost_rng.ndim; d++)
          skin_idx[d] = iter.idx[d];
        if (edge == GKYL_LOWER_EDGE)
          skin_idx[bc_dir] = local.upper[bc_dir];
        else
          skin_idx[bc_dir] = local.lower[bc_dir];
  
        linidx = gkyl_range_idx(&local_ext, skin_idx);
        double *fskin_c = gkyl_array_fetch(distf_ho, linidx);
  
        for (int k=0; k<distf_ho->ncomp; k++)
          TEST_CHECK( gkyl_compare(fskin_c[k], f_c[k], 1e-15) );
      }
    }
  }

  gkyl_free(red_integ_mom_skin_ho);
  gkyl_free(red_integ_mom_ghost_ho);
  if (use_gpu) {
    gkyl_cu_free(red_integ_mom_skin);
    gkyl_cu_free(red_integ_mom_ghost);
  }
  else {
    gkyl_free(red_integ_mom_skin);
    gkyl_free(red_integ_mom_ghost);
  }
  gkyl_dg_updater_moment_gyrokinetic_release(mcalc);
  gkyl_array_release(marr);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_position_map_release(pmap);
  gkyl_velocity_map_release(gvm);
  gkyl_array_release(buff_per);
  test_bc_twistshift_array_meta_release(mt);
  gkyl_twistshift_dg_release(tsup);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(distf);

}

// CBC geometry (see rt_gk_cbc_passive_3x2v_p1.c)

static double
interp_1x_lut_cbc(double x, double *lut_grid, double *lut_val, int N)
{
  double x_min = lut_grid[0];
  double x_max = lut_grid[N-1];
  if (x <= x_min) return lut_val[0];
  if (x >= x_max) return lut_val[N-1];
  double dx = (x_max - x_min)/(N-1);
  int idx = (int)((x - x_min)/dx);
  if (idx < 0) idx = 0;
  if (idx >= N-1) idx = N-2;
  return lut_val[idx] + (lut_val[idx+1] - lut_val[idx])*(x - lut_grid[idx])/(lut_grid[idx+1] - lut_grid[idx]);
}

struct gk_cbc_app_ctx {
  double a_shift, Z_axis, R_axis, R0, a_mid, r0, B0, kappa, delta, q0, Cy, qaxis, qlcfs;
  double Lx, Ly, Lz;
  double x_min, y_min, z_min, x_max, y_max, z_max;
  int psi_lut_size;
  double *r_lut;
  double *dPsidr_int_lut;
};

struct integrand_cbc_ctx {
  struct gk_cbc_app_ctx *app_ctx;
  double r;
  double theta;
};

static double r_x_cbc(double x, double r0) { return x + r0; }

static double qprofile_cbc(double r, double a_mid, double qaxis, double qlcfs)
{
  return 1.0 + 2.78*pow(r/a_mid, 2.8);
}

static double R_rtheta_cbc(double r, double theta, void *ctx)
{
  struct gk_cbc_app_ctx *app = ctx;
  return app->R_axis - app->a_shift*r*r/(2.*app->R_axis)
    + r*cos(theta + asin(app->delta)*sin(theta));
}

static double dRdr_cbc(double r, double theta, void *ctx)
{
  struct gk_cbc_app_ctx *app = ctx;
  return -app->a_shift*r/app->R_axis + cos(theta + asin(app->delta)*sin(theta));
}

static double dRdtheta_cbc(double r, double theta, void *ctx)
{
  struct gk_cbc_app_ctx *app = ctx;
  return -r*sin(theta + asin(app->delta)*sin(theta))*(1. + asin(app->delta)*cos(theta));
}

static double dZdr_cbc(double r, double theta, void *ctx)
{
  struct gk_cbc_app_ctx *app = ctx;
  return app->kappa*sin(theta);
}

static double dZdtheta_cbc(double r, double theta, void *ctx)
{
  struct gk_cbc_app_ctx *app = ctx;
  return app->kappa*r*cos(theta);
}

static double Jr_cbc(double r, double theta, void *ctx)
{
  return R_rtheta_cbc(r, theta, ctx)
    *(dRdr_cbc(r, theta, ctx)*dZdtheta_cbc(r, theta, ctx)
      - dRdtheta_cbc(r, theta, ctx)*dZdr_cbc(r, theta, ctx));
}

static double Bphi_cbc(double R, void *ctx)
{
  struct gk_cbc_app_ctx *app = ctx;
  return app->B0*app->R0/R;
}

static double integrand_JoRsq_cbc(double t, void *int_ctx)
{
  struct integrand_cbc_ctx *inctx = int_ctx;
  return Jr_cbc(inctx->r, t, inctx->app_ctx)
    / pow(R_rtheta_cbc(inctx->r, t, inctx->app_ctx), 2);
}

static double intdPsidr_cbc(double r, void *ctx)
{
  struct gk_cbc_app_ctx *app = ctx;
  struct integrand_cbc_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral = gkyl_dbl_exp(integrand_JoRsq_cbc, &tmp_ctx, 0., 2.*M_PI, 7, 1e-10);
  return integral.res;
}

static double dPsidr_cbc(double r, double theta, void *ctx)
{
  struct gk_cbc_app_ctx *app = ctx;
  double integral_val = interp_1x_lut_cbc(r, app->r_lut, app->dPsidr_int_lut, app->psi_lut_size);
  double R = R_rtheta_cbc(r, theta, ctx);
  double Bt = Bphi_cbc(R, ctx);
  return (R*Bt/(2.*M_PI*qprofile_cbc(r, app->a_mid, app->qaxis, app->qlcfs)))*integral_val;
}

static double compute_alpha_integral_cbc(double r, double twrap, void *ctx)
{
  struct gk_cbc_app_ctx *app = ctx;
  struct integrand_cbc_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral;
  if (twrap == 0.0) return 0.0;
  if (0. < twrap) {
    integral = gkyl_dbl_exp(integrand_JoRsq_cbc, &tmp_ctx, 0., twrap, 7, 1e-10);
    return integral.res;
  } else {
    integral = gkyl_dbl_exp(integrand_JoRsq_cbc, &tmp_ctx, twrap, 0., 7, 1e-10);
    return -integral.res;
  }
}

static double alpha_cbc(double r, double theta, double phi, void *ctx)
{
  double twrap = theta;
  while (twrap < -M_PI) twrap += 2.*M_PI;
  while (twrap >  M_PI) twrap -= 2.*M_PI;
  double integral_val = compute_alpha_integral_cbc(r, twrap, ctx);
  double R  = R_rtheta_cbc(r, theta, ctx);
  double Bt = Bphi_cbc(R, ctx);
  return phi - R*Bt*integral_val/dPsidr_cbc(r, theta, ctx);
}

void bc_shift_func_lo_cbc(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_cbc_app_ctx *app = ctx;
  double r = r_x_cbc(xc[0], app->r0);
  fout[0] = app->Cy*(alpha_cbc(r, app->z_min, 0.0, ctx) - alpha_cbc(r, app->z_max, 0.0, ctx));
}

void bc_shift_func_up_cbc(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_cbc_app_ctx *app = ctx;
  double r = r_x_cbc(xc[0], app->r0);
  fout[0] = -app->Cy*(alpha_cbc(r, app->z_min, 0.0, ctx) - alpha_cbc(r, app->z_max, 0.0, ctx));
}

void init_donor_3x_cbc(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_cbc_app_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double sigx = app->Lx/7.0, sigy = app->Ly/7.0, sigz = app->Lz/10.0;
  double mux = 0.*app->Lx;
  fout[0] = exp(-pow(x-mux,2)/(2.*pow(sigx,2)) - pow(y,2)/(2.*pow(sigy,2)) - pow(z,2)/(2.*pow(sigz,2)));
}

void
test_bc_twistshift_3x_cbc_wcells(const int *cells, enum gkyl_edge_loc edge,
  bool check_distf, bool use_gpu, bool write_f)
{
  int bc_dir = 2;
  int poly_order = 1;

  // Physical parameters matching rt_gk_cbc_passive_3x2v_p1.c create_ctx().
  double eV  = GKYL_ELEMENTARY_CHARGE;
  double mp  = GKYL_PROTON_MASS;
  double qi  = eV;
  double mi  = mp; // AMU = 1 (hydrogen)
  double Te0 = 2000.*eV;

  double R_axis = 1.6714;
  double B_axis = 1.54;
  double a_mid  = 0.604;
  double R0     = R_axis + 0.5*a_mid;
  double r0     = 0.5*a_mid;
  double B0     = B_axis*(R_axis/R0);
  double qaxis  = 1.2;
  double qlcfs  = 2.0;

  double c_s      = sqrt(Te0/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s    = c_s/omega_ci;
  double q0       = qprofile_cbc(r0, a_mid, qaxis, qlcfs);
  double Cy       = r0/q0;

  double Lx = 150.*rho_s;
  double Ly = 150.*rho_s;
  Ly = 2.*M_PI*Cy/round(2.*M_PI*Cy/Ly); // adjust to integer toroidal mode number
  double Lz = 2.*M_PI - 1e-10;

  // Use a fixed LUT size sufficient for accurate geometry evaluation.
  int psi_lut_size = 200;

  struct gk_cbc_app_ctx app_ctx = {
    .a_shift = 0.0, .Z_axis = 0.0,
    .R_axis  = R_axis, .R0 = R0, .a_mid = a_mid, .r0 = r0,
    .B0      = B0, .kappa = 1.0, .delta = 0.0,
    .q0 = q0, .Cy = Cy, .qaxis = qaxis, .qlcfs = qlcfs,
    .Lx = Lx, .Ly = Ly, .Lz = Lz,
    .x_min = -Lx/2., .x_max = Lx/2.,
    .y_min = -Ly/2., .y_max = Ly/2.,
    .z_min = -Lz/2., .z_max = Lz/2.,
    .psi_lut_size = psi_lut_size,
  };

  // Populate lookup tables (avoids redundant integration in geometry evaluations).
  app_ctx.r_lut        = gkyl_malloc(psi_lut_size*sizeof(double));
  app_ctx.dPsidr_int_lut = gkyl_malloc(psi_lut_size*sizeof(double));
  double r_lut_min = 0.0, r_lut_max = 2.0*a_mid;
  for (int i=0; i<psi_lut_size; i++)
    app_ctx.r_lut[i] = i*(r_lut_max - r_lut_min)/(psi_lut_size - 1);
  for (int i=0; i<psi_lut_size; i++)
    app_ctx.dPsidr_int_lut[i] = intdPsidr_cbc(app_ctx.r_lut[i], &app_ctx);

  const double lower[] = {app_ctx.x_min, app_ctx.y_min, app_ctx.z_min};
  const double upper[] = {app_ctx.x_max, app_ctx.y_max, app_ctx.z_max};
  int ndim = 3, cdim = 3, vdim = 0;

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d=0; d<cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = cells[d];
  }

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  int ghost_conf[cdim];
  for (int d=0; d<cdim; d++) ghost_conf[d] = 1;
  struct gkyl_range local_conf, local_ext_conf;
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);

  int ghost[ndim];
  for (int d=0; d<cdim; d++) ghost[d] = ghost_conf[d];
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct skin_ghost_ranges skin_ghost;
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  struct gkyl_range skin_rng, ghost_rng;
  if (edge == GKYL_LOWER_EDGE) {
    skin_rng  = skin_ghost.upper_skin[bc_dir];
    ghost_rng = skin_ghost.lower_ghost[bc_dir];
  } else {
    skin_rng  = skin_ghost.lower_skin[bc_dir];
    ghost_rng = skin_ghost.upper_ghost[bc_dir];
  }

  struct gkyl_array *distf    = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu?
    mkarr(false, basis.num_basis, local_ext.volume) : gkyl_array_acquire(distf);

  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew(
    &(struct gkyl_proj_on_basis_inp){
      .grid = &grid, .basis = &basis,
      .num_ret_vals = 1,
      .eval = init_donor_3x_cbc,
      .ctx  = &app_ctx,
    }
  );
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  struct gkyl_msgpack_data *mt = test_bc_twistshift_array_meta_new(
    (struct test_bc_twistshift_output_meta){
      .poly_order = poly_order, .basis_type = basis.id
    }
  );
  if (write_f)
    gkyl_grid_sub_array_write(&grid, &local, mt, distf_ho,
      "ctest_bc_twistshift_3x_cbc_do.gkyl");

  // Range extended only in bc_dir.
  struct gkyl_range update_rng;
  int lower_bcdir_ext[ndim], upper_bcdir_ext[ndim];
  for (int d=0; d<ndim; d++) {
    lower_bcdir_ext[d] = local.lower[d];
    upper_bcdir_ext[d] = local.upper[d];
  }
  lower_bcdir_ext[bc_dir] = local.lower[bc_dir] - ghost[bc_dir];
  upper_bcdir_ext[bc_dir] = local.upper[bc_dir] + ghost[bc_dir];
  gkyl_sub_range_init(&update_rng, &local_ext, lower_bcdir_ext, upper_bcdir_ext);

  // Choose shift function based on edge.
  struct gkyl_twistshift_dg_inp tsinp = {
    .bc_dir   = bc_dir,
    .shift_dir = 1,
    .shear_dir = 0,
    .edge      = edge,
    .cdim      = cdim,
    .bcdir_ext_update_r = &update_rng,
    .num_ghost = ghost,
    .basis = &basis,
    .grid = &grid,
    .shift_func = (edge == GKYL_LOWER_EDGE) ? bc_shift_func_lo_cbc : bc_shift_func_up_cbc,
    .shift_func_ctx = &app_ctx,
    .use_gpu   = use_gpu,
  };

  struct gkyl_twistshift_dg *tsup = gkyl_twistshift_dg_inew(&tsinp);

  struct gkyl_array *buff_per = mkarr(use_gpu, basis.num_basis, skin_rng.volume);
  apply_periodic_bc(buff_per, distf, bc_dir, skin_ghost);

  gkyl_twistshift_dg_advance(tsup, distf, distf);
  gkyl_array_copy(distf_ho, distf);

  if (write_f) {
    double lower_ext[ndim], upper_ext[ndim];
    int cells_ext[ndim];
    for (int d=0; d<ndim; d++) {
      double dx = (upper[d] - lower[d])/cells[d];
      lower_ext[d] = lower[d] - dx*ghost[d];
      upper_ext[d] = upper[d] + dx*ghost[d];
      cells_ext[d] = cells[d] + 2*ghost[d];
    }
    struct gkyl_rect_grid grid_ext;
    gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);
    gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho,
      "ctest_bc_twistshift_3x_cbc_tar.gkyl");
  }

  if (check_distf) {
    // Check 4 x=const lines, and 4 y=const lines.
    int idx0_ref[] = {4,12,19,28};
    int idx1_ref[] = {2,7,10,13};
    int num_ref0 = sizeof(idx0_ref)/sizeof(idx0_ref[0]);
    int num_ref1 = sizeof(idx1_ref)/sizeof(idx1_ref[0]);
    double f0[num_ref0*cells[1]*basis.num_basis];
    double f1[num_ref1*cells[0]*basis.num_basis];
    if (edge == GKYL_LOWER_EDGE) {
      const double f0m[] =
      {
        // idx[0] = 4;
         1.4147068739149354e-06,  3.9678563980892827e-07, -5.3677388121514690e-07,
        -1.4021443738919494e-06, -2.1663045016556100e-07, -3.9326221053808132e-07,
         5.3200736596068834e-07,  2.1470678662411421e-07,
         2.2073750885337050e-06,  5.4194161567241438e-08,  7.7188487204183924e-07,
        -2.1877737491242868e-06,  2.2856919929920019e-08, -5.3712921129036553e-08,
        -7.6503058731221621e-07, -2.2653952048417829e-08,
         8.0151330159023854e-06,  3.5840433095651838e-07,  2.4919021860714378e-06,
        -7.9439591843807564e-06,  1.4691992411677311e-07, -3.5522172507617932e-07,
        -2.4697742655481038e-06, -1.4561528526605157e-07,
         2.4056337321823511e-05,  1.5575127200312389e-06,  6.4725989939125481e-06,
        -2.3842718696134587e-05,  5.6131534500851025e-07, -1.5436821138880369e-06,
        -6.4151227587227262e-06, -5.5633090323856011e-07,
         5.9715709406034343e-05,  5.0593586542458534e-06,  1.3398359967595644e-05,
        -5.9185438001676253e-05,  1.5246943425719446e-06, -5.0144318963554509e-06,
        -1.3279383449912431e-05, -1.5111551613699758e-06,
         1.2266163264877539e-04,  1.2843394764617348e-05,  2.1715606284302345e-05,
        -1.2157240576270994e-04,  3.0737503418026385e-06, -1.2729346299878336e-05,
        -2.1522773189704626e-05, -3.0464556495582395e-06,
         2.0858381754799334e-04,  2.6012564754977512e-05,  2.6550477975637763e-05,
        -2.0673160755237084e-04,  4.5938737655600021e-06, -2.5781574963835826e-05,
        -2.6314711551985457e-05, -4.5530804815592547e-06,
         2.9374366736149827e-04,  4.2514292127503834e-05,  2.2141578694173991e-05,
        -2.9113524373960048e-04,  4.8573878417790354e-06, -4.2136768128945044e-05,
        -2.1944962993788802e-05, -4.8142545708526606e-06,
         3.4269905278866652e-04,  5.6458355498447365e-05,  7.3950701740712856e-06,
        -3.3965590870142440e-04,  2.9844920494078483e-06, -5.5957009173406751e-05,
        -7.3294024580625872e-06, -2.9579899646787672e-06,
         3.3130814623591096e-04,  6.1194737241075030e-05, -1.1197518621951395e-05,
        -3.2836615261769482e-04, -4.4059942666451672e-07, -6.0651332170970138e-05,
         1.1098085424488793e-05,  4.3668693397105754e-07,
         2.6547493073777388e-04,  5.4303954632758055e-05, -2.4003510429959307e-05,
        -2.6311753155848911e-04, -3.5659721755455483e-06, -5.3821739239660803e-05,
         2.3790360903447407e-05,  3.5343065871729570e-06,
         1.7634549859606930e-04,  3.9539852198925795e-05, -2.6084796724879738e-05,
        -1.7477956266190929e-04, -4.8370620921015304e-06, -3.9188741022952437e-05,
         2.5853165518798813e-05,  4.7941093124383771e-06,
         9.7120871429426495e-05,  2.3660577068010910e-05, -1.9999449126191993e-05,
        -9.6258444751461542e-05, -4.1887526472941376e-06, -2.3450472766235517e-05,
         1.9821855389468108e-05,  4.1515568110413017e-06,
         4.4351391469818973e-05,  1.1650043122384966e-05, -1.1710592695891665e-05,
        -4.3957554155084504e-05, -2.6752574970654273e-06, -1.1546591538391639e-05,
         1.1606603435837950e-05,  2.6515013939545564e-06,
         1.6794598764719960e-05,  4.7243461185578233e-06, -5.4058549896015648e-06,
        -1.6645463879424616e-05, -1.3229289098014698e-06, -4.6823942490099374e-06,
         5.3578513680151396e-06,  1.3111813918058868e-06,
         5.2735603406759165e-06,  1.5789778675299268e-06, -1.9989826780664712e-06,
        -5.2267314864994830e-06, -5.1808733150377238e-07, -1.5649566523489707e-06,
         1.9812318489709252e-06,  5.1348675152926985e-07,
        // idx[0] = 12;
         8.5771895799448463e-03,  3.1197886695786712e-04,  1.2589091675382295e-04,
        -8.5010247246790835e-03,  2.3881143669295742e-04, -3.0920851601408221e-04,
        -1.2477301404636147e-04, -2.3669080751226941e-04,
         8.2661626046304885e-03,  1.0720660909669306e-03, -3.0982027320798947e-04,
        -8.1927596475759902e-03,  1.8255117221452561e-04, -1.0625462175989969e-03,
         3.0706909042867208e-04, -1.8093013032420406e-04,
         6.6575398175439011e-03,  1.4760068941306471e-03, -5.9457869281011441e-04,
        -6.5984213205230440e-03,  4.2899447104076090e-05, -1.4629000541319605e-03,
         5.8929887479281385e-04, -4.2518502955738542e-05,
         4.4799350811575711e-03,  1.3955105208072918e-03, -6.3033616055554532e-04,
        -4.4401535648606475e-03, -8.3328432799330363e-05, -1.3831184830834593e-03,
         6.2473881867683195e-04,  8.2588482030546096e-05,
         2.5178444512609479e-03,  1.0028214075435281e-03, -4.8235931958416503e-04,
        -2.4954861652017948e-03, -1.3193983394577959e-04, -9.9391642221575223e-04,
         4.7807600190536255e-04,  1.3076821726847499e-04,
         1.1813969572080783e-03,  5.6871985765870141e-04, -2.8575381884705388e-04,
        -1.1709062332455860e-03, -1.1130086733193673e-04, -5.6366966432420582e-04,
         2.8321634453203873e-04,  1.1031252326278932e-04,
         4.6254048415618595e-04,  2.5918239859045285e-04, -1.3487188316997178e-04,
        -4.5843315637685304e-04, -6.6255494613212145e-05, -2.5688087666510880e-04,
         1.3367423009662933e-04,  6.5667150364697886e-05,
         1.5102266511218463e-04,  9.5873653213552430e-05, -5.1458764419249627e-05,
        -1.4968159420277919e-04, -3.0014361827430350e-05, -9.5022301747813123e-05,
         5.1001814120131696e-05,  2.9747836352719439e-05,
         4.9288640555562537e-05,  1.9197719093399958e-05, -9.3898074536993173e-06,
        -4.8850960807531967e-05, -1.7734142076995751e-05, -1.9027244664386457e-05,
         9.3064266074420962e-06,  1.7576664111519137e-05,
         7.4312889120419612e-05, -3.9559856706301425e-05,  2.5345644730922928e-05,
        -7.3652995761239122e-05, -1.7079041586904331e-05,  3.9208567891674429e-05,
        -2.5120577143858990e-05,  1.6927380868854436e-05,
         2.4668363462348900e-04, -1.2316540125222137e-04,  8.0544516023109083e-05,
        -2.4449310086503489e-04, -3.4617348824034278e-05,  1.2207170093550459e-04,
        -7.9829286244383513e-05,  3.4309949140455846e-05,
         6.9889469196789068e-04, -2.9144378275440763e-04,  1.9257347432050185e-04,
        -6.9268855503182346e-04, -6.3595836624040557e-05,  2.8885578195009997e-04,
        -1.9086343507478647e-04,  6.3031109956076864e-05,
         1.6531830617114835e-03, -5.4719482372872844e-04,  3.7015196128372313e-04,
        -1.6385029094950184e-03, -8.0981230730739408e-05,  5.4233577121939294e-04,
        -3.6686503725170299e-04,  8.0262122955358892e-05,
         3.2662572681362470e-03, -7.8998871116982131e-04,  5.6119422549585842e-04,
        -3.2372531275876876e-03, -5.0252394101103727e-05,  7.8297366559026060e-04,
        -5.5621085925887826e-04,  4.9806156287680392e-05,
         5.3918935380395400e-03, -8.1066940565628718e-04,  6.4417804850085126e-04,
        -5.3440138931855393e-03,  4.8296611837007029e-05,  8.0347071692792157e-04,
        -6.3845779160642867e-04, -4.7867741235965884e-05,
         7.4385809600035612e-03, -4.2620301983017786e-04,  4.9868993294226453e-04,
        -7.3725268711996870e-03,  1.7454031661297211e-04,  4.2241836624211989e-04,
        -4.9426160053675874e-04, -1.7299041057102665e-04,
        // idx[0] = 19;
         5.9328259833403563e-03, -3.1592831855465119e-03,  7.6447171990714363e-04,
        -5.8801428954679998e-03, -1.4953722986172529e-04,  3.1312289675153111e-03,
        -7.5768326346009965e-04,  1.4820934951547608e-04,
         8.5438428607558271e-03, -3.1770890503348310e-03,  7.1469750720526751e-04,
        -8.4679741220696214e-03,  1.4653149000110456e-04,  3.1488767174453966e-03,
        -7.0835104235361418e-04, -1.4523030041869057e-04,
         1.0538391888353600e-02, -2.0861147480127039e-03,  4.1035235387424735e-04,
        -1.0444811691083961e-02,  4.6474412137478837e-04,  2.0675901921113540e-03,
        -4.0670845311287489e-04, -4.6061722544807500e-04,
         1.1157187397402466e-02, -1.8349902862187543e-04, -6.0616833961144544e-05,
        -1.1058112338448060e-02,  5.9620026364061495e-04,  1.8186956983165916e-04,
         6.0078560632534608e-05, -5.9090604618554376e-04,
         1.0147502273615178e-02,  1.6812950536897154e-03, -5.0484232418428394e-04,
        -1.0057393149318094e-02,  4.4927500173372542e-04, -1.6663652689123407e-03,
         5.0035935896646541e-04, -4.4528547052858725e-04,
         7.9236576295588640e-03,  2.7048082887694922e-03, -7.4744339840912392e-04,
        -7.8532960931947931e-03,  1.3680114431856099e-04, -2.6807897766548171e-03,
         7.4080615228921856e-04, -1.3558635953852334e-04,
         5.3016792692569745e-03,  2.6754726903120800e-03, -7.3909622093528999e-04,
        -5.2546007209230746e-03, -1.3556652552430490e-04, -2.6517146762999374e-03,
         7.3253309717892176e-04,  1.3436270407449278e-04,
         3.0306885085190473e-03,  1.9748354796544207e-03, -5.6082860635105409e-04,
        -3.0037762024011207e-03, -2.4639127781914962e-04, -1.9572990760248132e-03,
         5.5584848678700660e-04,  2.4420334016906899e-04,
         1.4749219675761613e-03,  1.1511708985332411e-03, -3.4089330270885629e-04,
        -1.4618247616508778e-03, -2.1626062462549767e-04, -1.1409485798989469e-03,
         3.3786619355848882e-04,  2.1434024510948494e-04,
         6.0905145791109738e-04,  5.4222206411941179e-04, -1.6846528502594774e-04,
        -6.0364312273220987e-04, -1.3353062399755835e-04, -5.3740716937438697e-04,
         1.6696932484788329e-04,  1.3234488121368636e-04,
         2.2538939529024023e-04,  1.9086017837583353e-04, -6.0622757673701499e-05,
        -2.2338795291675480e-04, -7.2883062312442410e-05, -1.8916535308061285e-04,
         6.0084431742929801e-05,  7.2235865717262487e-05,
         1.3692552150240531e-04, -2.3628334870109403e-05,  6.7351766555291856e-06,
        -1.3570963226150468e-04, -5.5072031551511001e-05,  2.3418516876840270e-05,
        -6.6753687487111750e-06,  5.4582995688048192e-05,
         2.8180957439635708e-04, -2.6245954825184062e-04,  7.9675931947129452e-05,
        -2.7930712470157682e-04, -8.6112162611900079e-05,  2.6012892546224296e-04,
        -7.8968415135433089e-05,  8.5347492513282593e-05,
         7.5858237184660477e-04, -6.8983776721580120e-04,  2.0188206214678025e-04,
        -7.5184621240646007e-04, -1.6201786760815427e-04,  6.8371205515042964e-04,
        -2.0008936327952937e-04,  1.6057916005461062e-04,
         1.7743608964103583e-03, -1.4228383602125521e-03,  3.9110700695123274e-04,
        -1.7586046933318674e-03, -2.5501449371086789e-04,  1.4102036531488827e-03,
        -3.8763400355073456e-04,  2.5274998249502926e-04,
         3.5176447536171261e-03, -2.3790827021664203e-03,  6.1388697056459030e-04,
        -3.4864083094370800e-03, -2.8116612144444690e-04,  2.3579566109232599e-03,
        -6.0843569636494840e-04,  2.7866938556775280e-04,
        // idx[0] = 28;
         3.9008280658611207e-04, -2.5603798505596852e-04,  1.8683093211897101e-05,
        -3.8661889801461870e-04,  1.1175426823790703e-05,  2.5376438530717388e-04,
        -1.8517188625419492e-05, -1.1076189799980769e-05,
         4.3680260689641130e-04, -1.8748123385802201e-04,  8.0091734113950957e-06,
        -4.3292382970210014e-04,  2.7606287551937715e-05,  1.8581641335839929e-04,
        -7.9380525007527933e-06, -2.7361145611652907e-05,
         4.4457679107663103e-04, -7.5578664486840552e-05, -3.4115616360440883e-06,
        -4.4062897965993530e-04,  3.5399985354932533e-05,  7.4907531129210281e-05,
         3.3812672026731489e-06, -3.5085635912631181e-05,
         4.1504347649095089e-04,  4.2960315696147049e-05, -1.3312487632414327e-05,
        -4.1135791888244848e-04,  3.1477256709139010e-05, -4.2578831038885959e-05,
         1.3194273655178598e-05, -3.1197740828200509e-05,
         3.5577914498131479e-04,  1.3094475217078906e-04, -2.0440326610361184e-05,
        -3.5261985057240502e-04,  1.8525522986963621e-05, -1.2978197175140562e-04,
         2.0258817911812268e-05, -1.8361017613277599e-05,
         2.7796534980573800e-04,  1.6786465218620076e-04, -2.3905953966682069e-05,
        -2.7549703656169430e-04,  2.9974039875550468e-06, -1.6637402558656323e-04,
         2.3693670734873877e-05, -2.9707872457006951e-06,
         1.9535178824326644e-04,  1.5650334125842896e-04, -2.3209649131201253e-05,
        -1.9361707776044740e-04, -8.6672787668020905e-06, -1.5511360232069719e-04,
         2.3003549038581116e-05,  8.5903139257351561e-06,
         1.2172339077131273e-04,  1.1660274450623213e-04, -1.8887207777222881e-05,
        -1.2064249540877665e-04, -1.3462333943164936e-05, -1.1556732013136769e-04,
         1.8719490667402875e-05,  1.3342789340965859e-05,
         6.7058946442705194e-05,  7.0433021024861683e-05, -1.2563048496751787e-05,
        -6.6463467597043934e-05, -1.2699023524379705e-05, -6.9807580628297237e-05,
         1.2451489487645912e-05,  1.2586257066353350e-05,
         3.4938222906600265e-05,  3.1273776945169202e-05, -6.0324138429606134e-06,
        -3.4627973883173486e-05, -1.0006175581767381e-05, -3.0996067950580566e-05,
         5.9788464217243402e-06,  9.9173214288020404e-06,
         2.4837446165027290e-05, -9.6145358765583846e-07,  1.9699547938565060e-07,
        -2.4616891346377383e-05, -9.0936978967278239e-06,  9.5291594573175467e-07,
        -1.9524617303820833e-07,  9.0129464830299899e-06,
         3.6725180241487430e-05, -3.6057780516486035e-05,  6.7954495238266452e-06,
        -3.6399063159472330e-05, -1.1694827060907848e-05,  3.5737589898260177e-05,
        -6.7351063980715334e-06,  1.1590977798611619e-05,
         7.3131479750304169e-05, -8.5125203228748116e-05,  1.4329502961988946e-05,
        -7.2482077225311949e-05, -1.6767541211874010e-05,  8.4369297261769594e-05,
        -1.4202257958371192e-05,  1.6618646595792397e-05,
         1.3604357772040060e-04, -1.5068140809657870e-04,  2.1787568445735812e-05,
        -1.3483552007980250e-04, -2.0478489031778651e-05,  1.4934336752606917e-04,
        -2.1594096332079554e-05,  2.0296641453544639e-05,
         2.2070902006211332e-04, -2.1925837439635583e-04,  2.6427797458379611e-05,
        -2.1874913909968418e-04, -1.7877051979819007e-05,  2.1731137506795665e-04,
        -2.6193120429307424e-05,  1.7718304984205783e-05,
         3.1225316368462653e-04, -2.6354343092504682e-04,  2.5533068602001434e-05,
        -3.0948037700471862e-04, -6.4354644170033844e-06,  2.6120318333164132e-04,
        -2.5306336703815414e-05,  6.3783179343098985e-06,
      };
  
      const double f1m[] =
      {
        // idx[1] = 2;
         7.5742778736170694e-07,  1.3311910271755548e-07,  2.8478156167525587e-07,
        -7.5070188055260917e-07,  5.4750238065324756e-08, -1.3193701421442880e-07,
        -2.8225272093725106e-07, -5.4264059705951382e-08,
         1.2925859849765668e-06,  1.7301030138070181e-07,  5.1100320791830013e-07,
        -1.2811079100725808e-06,  7.3239351158981957e-08, -1.7147398176909403e-07,
        -5.0646553447542282e-07, -7.2588990743278259e-08,
         1.8635824605983989e-06,  1.4918689410169481e-07,  7.1992200307670070e-07,
        -1.8470339761485296e-06,  3.9576065168750808e-08, -1.4786212471296968e-07,
        -7.1352914505999257e-07, -3.9224632424088840e-08,
         2.2073750885337050e-06,  5.4194161567241438e-08,  7.7188487204183924e-07,
        -2.1877737491242868e-06,  2.2856919929920019e-08, -5.3712921129036553e-08,
        -7.6503058731221621e-07, -2.2653952048417829e-08,
         2.3371245739281691e-06,  7.7901820168274173e-08,  7.1291202449929149e-07,
        -2.3163710679863721e-06, -2.4448859078376251e-07, -7.7210057347509384e-08,
        -7.0658141461162244e-07,  2.4231754886407394e-07,
         4.4362078536198185e-06,  1.4349007654002696e-06, -1.4918804712085053e-06,
        -4.3968146321047141e-06, -7.9055943126041467e-07, -1.4221589450057450e-06,
         1.4786326749339418e-06,  7.8353931772558917e-07,
         1.8432855790288429e-05,  7.7265468592736174e-06, -7.0218535854514018e-06,
        -1.8269173294953980e-05, -3.1418179249565271e-06, -7.6579356530322741e-06,
         6.9595000071553346e-06,  3.1139187972416376e-06,
         8.2880124644381316e-05,  3.4710110035730186e-05, -2.7056128318655973e-05,
        -8.2144154821269593e-05, -1.1646701134768506e-05, -3.4401886638953268e-05,
         2.6815871754634015e-05,  1.1543279227396033e-05,
         3.6130932397451928e-04,  1.4417479842687968e-04, -9.7273919862532121e-05,
        -3.5810092195539682e-04, -3.5891598669338311e-05, -1.4289453610402862e-04,
         9.6410134124976947e-05,  3.5572883734518748e-05,
         1.3866911075139233e-03,  4.8934908019720340e-04, -2.7681619210264059e-04,
        -1.3743773855753208e-03, -7.3086254588880316e-05, -4.8500369392349860e-04,
         2.7435808330019372e-04,  7.2437253660220146e-05,
         4.1402065573272135e-03,  1.1066127780317765e-03, -5.1267688697286180e-04,
        -4.1034418070240936e-03, -3.9937840842209967e-05, -1.0967861324518405e-03,
         5.0812435137475719e-04,  3.9583195554377857e-05,
         8.2661626046304885e-03,  1.0720660909669306e-03, -3.0982027320798947e-04,
        -8.1927596475759902e-03,  1.8255117221452561e-04, -1.0625462175989969e-03,
         3.0706909042867208e-04, -1.8093013032420406e-04,
         9.2071430199458723e-03, -6.7468017399297192e-04,  5.5519252734112100e-04,
        -9.1253842213336704e-03,  2.3350197769444222e-04,  6.6868906031594807e-04,
        -5.5026245577217254e-04, -2.3142849614555009e-04,
         4.6562742627774770e-03, -1.5981825134639602e-03,  7.7093816367743357e-04,
        -4.6149268666407072e-03, -1.2283002048105933e-04,  1.5839907623441212e-03,
        -7.6409228565316330e-04,  1.2173929832259107e-04,
         8.4265128323255962e-04, -5.4995291705314407e-04,  2.3627440955229246e-04,
        -8.3516859762456210e-04, -1.3174780227490776e-04,  5.4506937286423181e-04,
        -2.3417630899862953e-04,  1.3057789082567117e-04,
         1.1584483866390626e-04,  4.3858963858025070e-05, -1.6873108923339463e-05,
        -1.1481614444094626e-04, -4.5483891179512000e-05, -4.3469499266713753e-05,
         1.6723276873219639e-05,  4.5079997345019702e-05,
         1.9929185999929553e-03,  1.3468353897936511e-03, -4.5601060935752653e-04,
        -1.9752216194948405e-03, -2.4418116907756006e-04, -1.3348755838951901e-03,
         4.5196126641860161e-04,  2.4201285704154009e-04,
         1.0151660891080033e-02,  2.5306058282125788e-03, -7.5599768765257259e-04,
        -1.0061514838545048e-02,  2.8830814950777197e-04, -2.5081342220012403e-03,
         7.4928448003081915e-04, -2.8574799291165951e-04,
         8.5438428607558271e-03, -3.1770890503348310e-03,  7.1469750720526751e-04,
        -8.4679741220696214e-03,  1.4653149000110456e-04,  3.1488767174453966e-03,
        -7.0835104235361418e-04, -1.4523030041869057e-04,
         8.3058989180792564e-04, -7.7938148858281266e-04,  2.0964975534781683e-04,
        -8.2321431052863629e-04, -1.7368462202414781e-04,  7.7246063441240638e-04,
        -2.0778807989763133e-04,  1.7214231449146951e-04,
         4.8784698194571816e-04,  4.6603525819394176e-04, -1.2686455710601843e-04,
        -4.8351493420147626e-04, -1.1084020143927416e-04, -4.6189689706081717e-04,
         1.2573800853899969e-04,  1.0985594805166225e-04,
         4.8781825809830056e-03,  1.3781496758297112e-03, -3.2928825061709555e-04,
        -4.8348646541985394e-03,  2.0371545823871589e-04, -1.3659118012187134e-03,
         3.2636419353345472e-04, -2.0190647893989855e-04,
         2.3680721830718248e-03, -1.7400576716845859e-03,  2.7950593809693524e-04,
        -2.3470438644831496e-03, -9.5410574014376798e-05,  1.7246060788892229e-03,
        -2.7702394453451581e-04,  9.4563334660169134e-05,
         8.9363290665012558e-05,  3.9704854712887114e-05, -9.8017374005486476e-06,
        -8.8569750772238376e-05, -3.3694653872762494e-05, -3.9352278325905250e-05,
         9.7146986446126184e-06,  3.3395447657075059e-05,
         1.5558629382148325e-03,  5.9445795666028924e-04, -1.1056083435060165e-04,
        -1.5420469820210249e-03,  6.1304104069997447e-05, -5.8917921077172018e-04,
         1.0957906172357881e-04, -6.0759727829951777e-05,
         5.8935835734624667e-04, -5.5038435706920711e-04,  7.6629057291856129e-05,
        -5.8412489554986722e-04, -4.8460576140242794e-05,  5.4549698172253955e-04,
        -7.5948596518150819e-05,  4.8030249547433478e-05,
         1.2304561996882203e-04,  1.2696504265504160e-04, -2.3555535292118696e-05,
        -1.2195298330168754e-04, -2.1706603802157628e-05, -1.2583760178323854e-04,
         2.3346363754109702e-05,  2.1513850648983980e-05,
         4.3680260689641130e-04, -1.8748123385802201e-04,  8.0091734113950957e-06,
        -4.3292382970210014e-04,  2.7606287551937715e-05,  1.8581641335839929e-04,
        -7.9380525007527933e-06, -2.7361145611652907e-05,
         1.6405328706551727e-05,  2.9047681450653114e-06, -5.1336046315663870e-07,
        -1.6259650512677665e-05, -5.7121417723422041e-06, -2.8789740031395213e-06,
         5.0880185744890011e-07,  5.6614183451294898e-06,
         1.2525989018365014e-04, -1.7226749700817367e-05, -1.8059029380175158e-06,
        -1.2414759094885760e-04,  1.0448806188062650e-05,  1.7073777344845405e-05,
         1.7898666437726345e-06, -1.0356021505668008e-05,
         6.1342326303602555e-06, -1.4359128614618765e-06,  2.5224463428371958e-07,
        -6.0797610652743729e-06, -2.0408390187610221e-06,  1.4231620537237336e-06,
        -2.5000471922964324e-07,  2.0227165082305203e-06,
         2.9913061136183802e-05, -9.4004861242035489e-06, -2.0777664618398835e-07,
        -2.9647435204664162e-05,  2.4116686151255107e-06,  9.3170104520845704e-06,
         2.0593160381472413e-07, -2.3902531632099425e-06,
        // idx[1] = 7;
         3.6490020128125056e-05,  7.3675063786135985e-06,  3.2021281911722955e-06,
        -3.6165991251789855e-05,  6.9414808136517716e-07, -7.3020834272183321e-06,
        -3.1736935124293500e-06, -6.8798409400575064e-07,
         7.1180658041772527e-05,  1.2900664940831104e-05,  6.7559365584194971e-06,
        -7.0548578679769132e-05,  1.4100360605950567e-06, -1.2786107920854779e-05,
        -6.6959443050874898e-06, -1.3975150370740391e-06,
         1.2807110995553272e-04,  2.0058649867687919e-05,  1.3810778814974065e-05,
        -1.2693384728195296e-04,  2.7419450947357848e-06, -1.9880530432439409e-05,
        -1.3688140046208898e-05, -2.7175968103308825e-06,
         2.0858381754799334e-04,  2.6012564754977512e-05,  2.6550477975637763e-05,
        -2.0673160755237084e-04,  4.5938737655600021e-06, -2.5781574963835826e-05,
        -2.6314711551985457e-05, -4.5530804815592547e-06,
         2.9871110310366976e-04,  2.5205913165246253e-05,  4.5389429729346178e-05,
        -2.9605856899303685e-04,  5.9426328345683972e-06, -2.4982086384903072e-05,
        -4.4986374706053409e-05, -5.8898626625296713e-06,
         3.6469434674854203e-04,  1.1237385446411157e-05,  6.5265128683803013e-05,
        -3.6145588595932354e-04,  6.3539671552325690e-06, -1.1137598234281373e-05,
        -6.4685578816824895e-05, -6.2975443626347055e-06,
         3.5893595493224776e-04, -1.4493122690417841e-05,  7.7817247724840170e-05,
        -3.5574862826746054e-04, -1.6084567762455801e-06,  1.4364424754832286e-05,
        -7.7126235901573401e-05,  1.5941737903768229e-06,
         2.6694095601170693e-04, -3.6280004828175835e-05,  7.0464806864678920e-05,
        -2.6457053862852633e-04, -7.8696688227856493e-06,  3.5957840873301803e-05,
        -6.9839084212036697e-05,  7.7997867033233192e-06,
         1.3743518008194365e-04, -3.5616161618941028e-05,  4.3245713680199841e-05,
        -1.3621476510780793e-04, -1.0173131963280475e-05,  3.5299892546240166e-05,
        -4.2861694708412319e-05,  1.0082795249096795e-05,
         4.4281012350574053e-05, -1.7888520118666143e-05,  1.6654315330488873e-05,
        -4.3887799997592024e-05, -6.9220705758281906e-06,  1.7729671286766401e-05,
        -1.6506426152469179e-05,  6.8606030638146820e-06,
         5.1438549738468361e-05,  2.9665824357334224e-05, -1.7900189835034524e-05,
        -5.0981778944329460e-05, -1.3973723052800237e-05, -2.9402393871455844e-05,
         1.7741237376853713e-05,  1.3849637350377501e-05,
         4.6254048415618595e-04,  2.5918239859045285e-04, -1.3487188316997178e-04,
        -4.5843315637685304e-04, -6.6255494613212145e-05, -2.5688087666510880e-04,
         1.3367423009662933e-04,  6.5667150364697886e-05,
         2.8294587601195142e-03,  1.2375254640891553e-03, -5.5518378989416852e-04,
        -2.8043333603674971e-03, -1.6372300936094285e-04, -1.2265363228347245e-03,
         5.5025379591314930e-04,  1.6226916026556920e-04,
         9.1957619553724914e-03,  2.1739429633878338e-03, -7.9765225550990928e-04,
        -9.1141042198332155e-03,  1.2436740398928127e-04, -2.1546384989570115e-03,
         7.9056915818215777e-04, -1.2326302997069018e-04,
         1.2500828315738250e-02, -7.8508011247857622e-04,  3.9877045177863924e-04,
        -1.2389821817572859e-02,  4.5450998846944439e-04,  7.7810865491877395e-04,
        -3.9522939751361548e-04, -4.5047397094109647e-04,
         5.4790816862939707e-03, -2.4387569363537088e-03,  8.7033830325430211e-04,
        -5.4304278166626358e-03, -1.8988685124215858e-04,  2.4171009420032898e-03,
        -8.6260975880722460e-04,  1.8820066902513577e-04,
         5.7215777750936180e-04, -4.7524252007978134e-04,  1.6913878684932202e-04,
        -5.6707705568235720e-04, -1.2685880911775257e-04,  4.7102239909252349e-04,
        -1.6763684602125322e-04,  1.2573231159247535e-04,
         3.7427031371474126e-04,  3.3208117201649346e-04, -1.1114484687922794e-04,
        -3.7094681899556086e-04, -9.6463580458111862e-05, -3.2913231398235817e-04,
         1.1015788826100110e-04,  9.5606990478896647e-05,
         5.3016792692569745e-03,  2.6754726903120800e-03, -7.3909622093528999e-04,
        -5.2546007209230746e-03, -1.3556652552430490e-04, -2.6517146762999374e-03,
         7.3253309717892176e-04,  1.3436270407449278e-04,
         9.1961936544017395e-03, -1.7368826364828002e-03,  2.6809009052306064e-04,
        -9.1145320854046884e-03,  4.6206186478473183e-04,  1.7214592377822928e-03,
        -2.6570946890422425e-04, -4.5795878711259467e-04,
         1.4411213884148821e-03, -1.2810931792856962e-03,  2.9718607020727523e-04,
        -1.4283243293145838e-03, -2.1624336633437261e-04,  1.2697171597080945e-03,
        -2.9454707828417920e-04,  2.1432314007080162e-04,
         2.5492665841420081e-04,  2.4888323947485508e-04, -6.6278567407695297e-05,
        -2.5266292716977288e-04, -7.1894578026316508e-05, -2.4667317337616465e-04,
         6.5690018274348292e-05,  7.1256159103809388e-05,
         3.4578494082033308e-03,  1.1303046461231621e-03, -2.4223248931282091e-04,
        -3.4271439425899294e-03,  1.3900706598020271e-04, -1.1202676183793090e-03,
         2.4008148141947591e-04, -1.3777269276708681e-04,
         1.5781234722698919e-03, -1.2618335200808000e-03,  1.9187865723672629e-04,
        -1.5641098440602522e-03, -8.2967902797388205e-05,  1.2506285249561542e-03,
        -1.9017478791906594e-04,  8.2231153510297416e-05,
         9.6977901902638649e-05,  8.6316756686137649e-05, -1.9677860154673218e-05,
        -9.6116744784267082e-05, -3.0120448337275501e-05, -8.5550269806171524e-05,
         1.9503122105962595e-05,  2.9852980821632315e-05,
         1.1828731451272050e-03,  1.7396447751082442e-04, -4.5145997322417246e-05,
        -1.1723693127172242e-03,  8.1902648560112776e-05, -1.7241968488061619e-04,
         4.4745103962203978e-05, -8.1175358657596279e-05,
         1.8720681143410243e-04, -2.0805519933339648e-04,  3.3414378701706986e-05,
        -1.8554442778678481e-04, -3.2611926466637799e-05,  2.0620768343126807e-04,
        -3.3117661310317896e-05,  3.2322334796065803e-05,
         1.9535178824326644e-04,  1.5650334125842896e-04, -2.3209649131201253e-05,
        -1.9361707776044740e-04, -8.6672787668020905e-06, -1.5511360232069719e-04,
         2.3003549038581116e-05,  8.5903139257351561e-06,
         1.6951067600112064e-04, -1.5119863365618911e-04,  1.3863877085521073e-05,
        -1.6800543282288719e-04, -4.5044095319893205e-06,  1.4985600015818046e-04,
        -1.3740766807754799e-05,  4.4644106842474151e-06,
         3.6276567298503912e-05,  3.6027257115504254e-05, -5.2295393006259687e-06,
        -3.5954433868657644e-05, -3.9956008741888273e-06, -3.5707337542985901e-05,
         5.1831013502663178e-06,  3.9601202124352467e-06,
         5.1575023302046870e-05, -4.4579128212970339e-05,  2.9869584382339701e-06,
        -5.1117040632022451e-05, -3.2845211577127062e-09,  4.4183268611574342e-05,
        -2.9604344521412661e-06,  3.2553548350773803e-09,
         1.2401275268881665e-05,  1.1336295347525180e-05, -1.3473909962905102e-06,
        -1.2291152794945234e-05, -6.9063884391199019e-07, -1.1235629822256499e-05,
         1.3354262566444594e-06,  6.8450601834049220e-07,
        // idx[1] = 10;
         3.9161987831707695e-05,  8.5585033279271904e-06, -3.0940904155168843e-06,
        -3.8814232065402164e-05, -6.2790787320485707e-07, -8.4825044052978823e-06,
         3.0666151048127765e-06,  6.2233209435129094e-07,
         8.2966973231524862e-05,  1.7351396003363992e-05, -5.8734092978937982e-06,
        -8.2230232198915169e-05, -9.6285893253993122e-07, -1.7197316796773389e-05,
         5.8212538261135436e-06,  9.5430881124964793e-07,
         1.6948528988022327e-04,  3.3617119358933897e-05, -9.1922219423833112e-06,
        -1.6798027212900247e-04, -8.7142547168146892e-07, -3.3318601644412122e-05,
         9.1105956419166347e-06,  8.6368727325331598e-07,
         3.3130814623591096e-04,  6.1194737241075030e-05, -1.1197518621951395e-05,
        -3.2836615261769482e-04, -4.4059942666451672e-07, -6.0651332170970138e-05,
         1.1098085424488793e-05,  4.3668693397105754e-07,
         6.0949868976399916e-04,  1.0001588819084044e-04, -9.8223350700817221e-06,
        -6.0408638319692828e-04,  3.3693722410867784e-06, -9.9127753962567809e-05,
         9.7351134082527874e-06, -3.3394524466498374e-06,
         1.0226479509198823e-03,  1.3603218361846981e-04,  2.0902400504249218e-05,
        -1.0135669072465832e-03,  1.1408417045058132e-05, -1.3482422715671544e-04,
        -2.0716788621210686e-05, -1.1307111084061294e-05,
         1.4988806244302360e-03,  1.3106170460485647e-04,  8.5695676969146863e-05,
        -1.4855706672751183e-03,  3.0453476425095055e-05, -1.2989788565587895e-04,
        -8.4934705234475415e-05, -3.0183051642870378e-05,
         1.8030298257953226e-03,  3.5581460837069158e-05,  1.8651702121554774e-04,
        -1.7870190445899432e-03,  3.0633437103590942e-05, -3.5265499904931722e-05,
        -1.8486076285806655e-04, -3.0361414282882376e-05,
         1.6413216149696490e-03, -1.2558118114409131e-04,  2.6036761791834394e-04,
        -1.6267467915868207e-03,  5.4223079993272169e-06,  1.2446602886760324e-04,
        -2.5805557132664797e-04, -5.3741582761426463e-06,
         1.0189032686435122e-03, -2.1096043183989262e-04,  2.2864831174272861e-04,
        -1.0098554774920258e-03, -2.7657765772399670e-05,  2.0908711767234067e-04,
        -2.2661793041463409e-04,  2.7412166709047420e-05,
         3.7915828230823104e-04, -1.4271127978161287e-04,  1.1053775807032913e-04,
        -3.7579138276314626e-04, -3.4225960821842870e-05,  1.4144401340396710e-04,
        -1.0955619035909259e-04,  3.3922036636883500e-05,
         7.4312889120419612e-05, -3.9559856706301425e-05,  2.5345644730922928e-05,
        -7.3652995761239122e-05, -1.7079041586904331e-05,  3.9208567891674429e-05,
        -2.5120577143858990e-05,  1.6927380868854436e-05,
         1.6452199128976365e-04,  1.1354492403363452e-04, -5.6206288820240358e-05,
        -1.6306104728965439e-04, -3.6099553950468539e-05, -1.1253665289486835e-04,
         5.5707180830016014e-05,  3.5778992387013618e-05,
         1.7045525792129491e-03,  9.4638345870692170e-04, -4.0697719285162913e-04,
        -1.6894162692038758e-03, -1.8078271052518686e-04, -9.3797963849442858e-04,
         4.0336326328867662e-04,  1.7917737245339021e-04,
         8.0791678878996270e-03,  2.6080011965470560e-03, -9.4517642723761621e-04,
        -8.0074254311060314e-03, -2.5821417318396076e-05, -2.5848423247725079e-03,
         9.3678332538179568e-04,  2.5592124903382838e-05,
         1.3356481058149027e-02, -4.3293862331198071e-04,  1.8456831679738910e-04,
        -1.3237876422309745e-02,  5.5117088081251008e-04,  4.2909415802691434e-04,
        -1.8292936280150354e-04, -5.4627652118893928e-04,
         5.8131600430988816e-03, -2.8057188210851179e-03,  8.7511944457940491e-04,
        -5.7615395805694532e-03, -1.9564030185418527e-04,  2.7808042303636956e-03,
        -8.6734844392524053e-04,  1.9390302939028701e-04,
         4.7634326079316506e-04, -4.2920915428436507e-04,  1.4088268702260373e-04,
        -4.7211336530385869e-04, -1.1855725319570570e-04,  4.2539780642850724e-04,
        -1.3963165842325711e-04,  1.1750447291771449e-04,
         6.0905145791109738e-04,  5.4222206411941179e-04, -1.6846528502594774e-04,
        -6.0364312273220987e-04, -1.3353062399755835e-04, -5.3740716937438697e-04,
         1.6696932484788329e-04,  1.3234488121368636e-04,
         6.4242974589046806e-03,  2.4764802964899451e-03, -6.3518193397605162e-04,
        -6.3672501380333101e-03,  8.6162573411439783e-05, -2.4544893212884973e-03,
         6.2954156196167615e-04, -8.5397455668379408e-05,
         6.1192473190725665e-03, -2.6624998539887439e-03,  4.5719660344097069e-04,
        -6.0649088225232604e-03,  1.3024981958765076e-04,  2.6388570378735065e-03,
        -4.5313672895591774e-04, -1.2909321012194870e-04,
         3.6469519543160190e-04, -3.8363095857578216e-04,  9.4090664425710294e-05,
        -3.6145672710613595e-04, -9.7469050278121223e-05,  3.8022434197216746e-04,
        -9.3255145778133679e-05,  9.6603531795857238e-05,
         9.1178047998913031e-04,  7.8751323767825392e-04, -1.7362275060952033e-04,
        -9.0368393185465785e-04, -1.1501287471403601e-04, -7.8052017413353706e-04,
         1.7208099249076104e-04,  1.1399156827389669e-04,
         3.0428994134042461e-03, -6.0542157405543088e-04,  2.7297328115299833e-05,
        -3.0158786752883344e-03,  2.0376965733455402e-04,  6.0004547199625242e-04,
        -2.7054929713624570e-05, -2.0196019675148278e-04,
         2.4486683302324581e-04, -2.7365088383468937e-04,  5.4183415554651385e-05,
        -2.4269243233840965e-04, -5.7765593385112273e-05,  2.7122088275259229e-04,
        -5.3702270540301433e-05,  5.7252638877285439e-05,
         4.8684120652023891e-04,  3.9028286351598619e-04, -6.7232236617642468e-05,
        -4.8251808999280144e-04, -3.0667260128909864e-05, -3.8681717845267502e-04,
         6.6635218967111602e-05,  3.0394936962056052e-05,
         6.3412598830576885e-04, -4.3372366375867428e-04,  3.8776407472195745e-05,
        -6.2849499301653115e-04,  9.5143346917004076e-06,  4.2987222736827919e-04,
        -3.8432075633040440e-05, -9.4298480521095317e-06,
         3.4938222906600265e-05,  3.1273776945169202e-05, -6.0324138429606134e-06,
        -3.4627973883173486e-05, -1.0006175581767381e-05, -3.0996067950580566e-05,
         5.9788464217243402e-06,  9.9173214288020404e-06,
         2.4392083227262556e-04, -5.6676719861446003e-05, -1.1626687511546724e-06,
        -2.4175483200958035e-04,  1.9445449116177598e-05,  5.6173434475836851e-05,
         1.1523443323775968e-06, -1.9272774861550602e-05,
         1.0411319188242268e-05, -3.0681569078749403e-06,  5.5818318185156826e-07,
        -1.0318867387835258e-05, -3.5693133209550033e-06,  3.0409118849402331e-06,
        -5.5322655347565932e-07,  3.5376180634403095e-06,
         6.2185931257129870e-05, -1.0179624868692800e-05, -8.9508031605233602e-07,
        -6.1633724452136231e-05,  5.2188375294746095e-06,  1.0089230497954324e-05,
         8.8713206422836305e-07, -5.1724946101088791e-06,
         3.4821562995004312e-06,  1.0893657152118253e-06, -1.5562699790180428e-07,
        -3.4512349903592222e-06, -1.0573421757721680e-06, -1.0796922223669664e-06,
         1.5424504083298160e-07,  1.0479530497614316e-06,
        // idx[1] = 13;
         7.6939600011633745e-06,  1.8061169092499393e-06, -1.7409386096072833e-06,
        -7.6256381639873142e-06, -3.8744241924183407e-07, -1.7900787149551884e-06,
         1.7254791941435799e-06,  3.8400195712888607e-07,
         1.7708033426125268e-05,  4.1983811964132660e-06, -3.8745042608883230e-06,
        -1.7550787303677916e-05, -8.9714935127936357e-07, -4.1610998593045459e-06,
         3.8400989287563202e-06,  8.8918272656446436e-07,
         4.1224027239434575e-05,  9.9239799414819331e-06, -8.9605036911705433e-06,
        -4.0857960704598385e-05, -2.1501854097138079e-06, -9.8358556801655410e-06,
         8.8809350328839017e-06,  2.1310919107303748e-06,
         9.7120871429426495e-05,  2.3660577068010910e-05, -1.9999449126191993e-05,
        -9.6258444751461542e-05, -4.1887526472941376e-06, -2.3450472766235517e-05,
         1.9821855389468108e-05,  4.1515568110413017e-06,
         2.2966717127649675e-04,  5.5667830351288630e-05, -3.9586034209482710e-05,
        -2.2762774254561342e-04, -7.9335419238213832e-06, -5.5173503835342345e-05,
         3.9234512940422678e-05,  7.8630926156020864e-06,
         5.3200292659768048e-04,  1.2461168529569132e-04, -7.7661558932109947e-05,
        -5.2727877709304340e-04, -1.3162232392700854e-05, -1.2350514207585148e-04,
         7.6971929614645792e-05,  1.3045352671689596e-05,
         1.1666581574647169e-03,  2.4674867349924984e-04, -1.2191453257652898e-04,
        -1.1562983130331923e-03, -1.1089595832661965e-05, -2.4455756220003919e-04,
         1.2083193988787077e-04,  1.0991120982167166e-05,
         2.2887686678893923e-03,  3.9334003412439531e-04, -1.3148940639503572e-04,
        -2.2684445590770812e-03,  1.6659921237230303e-05, -3.8984720159573514e-04,
         1.3032178948349212e-04, -1.6511982279143148e-05,
         3.7173943219990105e-03,  3.9400146321998456e-04, -3.1956132370684377e-06,
        -3.6843841153501047e-03,  7.3337556977308329e-05, -3.9050275724631495e-04,
         3.1672364106713046e-06, -7.2686324500672754e-05,
         4.5199940871419554e-03,  2.1470803254982653e-05,  2.8930333733463666e-04,
        -4.4798568496190452e-03,  9.9813558392372950e-05, -2.1280143994497898e-05,
        -2.8673434353886968e-04, -9.8927220838835293e-05,
         3.6305853350697154e-03, -5.0640172967673838e-04,  5.0404713216872657e-04,
        -3.5983459862715781e-03,  1.2471040562830394e-05,  5.0190491704508729e-04,
        -4.9957122820147365e-04, -1.2360298577867939e-05,
         1.6531830617114835e-03, -5.4719482372872844e-04,  3.7015196128372313e-04,
        -1.6385029094950184e-03, -8.0981230730739408e-05,  5.4233577121939294e-04,
        -3.6686503725170299e-04,  8.0262122955358892e-05,
         3.5645930339666049e-04, -2.0018386894217195e-04,  1.1080594537964248e-04,
        -3.5329396922766608e-04, -5.1667624794085499e-05,  1.9840624991411976e-04,
        -1.0982199618349089e-04,  5.1208819829162221e-05,
         8.1645729682003766e-05,  3.1518177127831518e-05, -1.4204895571968811e-05,
        -8.0920721201505936e-05, -3.1610417461887078e-05, -3.1238297876381250e-05,
         1.4078757073429143e-05,  3.1329719122596009e-05,
         1.1112012333347322e-03,  7.2994071394047412e-04, -2.9209530026815726e-04,
        -1.1013338425863719e-03, -1.6061672159759519e-04, -7.2345889045835620e-04,
         2.8950151403300714e-04,  1.5919045612453774e-04,
         7.2402254655031442e-03,  2.7989480758611107e-03, -9.5793435231419509e-04,
        -7.1759327598876653e-03, -9.4784436742530046e-05, -2.7740936088929552e-03,
         9.4942796095860336e-04,  9.3942757444357306e-05,
         1.3241006243061248e-02, -4.8348417204651978e-04,  1.1607405154285738e-04,
        -1.3123427015660945e-02,  6.0299863597532274e-04,  4.7919086575499243e-04,
        -1.1504332192525146e-04, -5.9764404944013681e-04,
         5.0513947154027959e-03, -2.8225427733082581e-03,  7.8810546751073931e-04,
        -5.0065386767431759e-03, -2.3600245213967414e-04,  2.7974787870448416e-03,
        -7.8110714500572343e-04,  2.3390676654918446e-04,
         2.8180957439635708e-04, -2.6245954825184062e-04,  7.9675931947129452e-05,
        -2.7930712470157682e-04, -8.6112162611900079e-05,  2.6012892546224296e-04,
        -7.8968415135433089e-05,  8.5347492513282593e-05,
         1.1583383425450180e-03,  9.6723264576604258e-04, -2.6977156771818526e-04,
        -1.1480523774994261e-03, -1.8655819799419011e-04, -9.5864368620225443e-04,
         2.6737601469716817e-04,  1.8490157398973845e-04,
         6.9931044556163944e-03,  1.2000797468667162e-03, -3.3704778425109526e-04,
        -6.9310061676217134e-03,  3.5590508141603505e-04, -1.1894231210132782e-03,
         3.3405482304091577e-04, -3.5274466870022149e-04,
         2.8733758677926486e-03, -2.0671103848013202e-03,  3.8129908297979793e-04,
        -2.8478604871362292e-03, -1.4871125794988407e-04,  2.0487545863420117e-03,
        -3.7791317327156090e-04,  1.4739071217756880e-04,
         1.0322441124946509e-04,  1.8476743202895893e-05, -4.8168812615395016e-06,
        -1.0230778555646509e-04, -4.0962290464150215e-05, -1.8312670990347184e-05,
         4.7741076862677817e-06,  4.0598548133929244e-05,
         1.9582591910971407e-03,  9.7865979301640363e-04, -1.8723018386833398e-04,
        -1.9408699837731578e-03,  2.3128246186372299e-05, -9.6996936116868003e-04,
         1.8556759267540323e-04, -2.2922868946319787e-05,
         1.2887036948865910e-03, -9.7941769109876026e-04,  1.2664939140064258e-04,
        -1.2772600944523653e-03, -3.3426030714678910e-05,  9.7072052916803990e-04,
        -1.2552475349033199e-04,  3.3129209854213705e-05,
         6.9313334777440286e-05,  5.9715287673186067e-05, -1.2957096051837794e-05,
        -6.8697837117967416e-05, -2.1501206316806118e-05, -5.9185020013786713e-05,
         1.2842037927465642e-05,  2.1310277079216685e-05,
         7.4156655404620110e-04,  3.4961769651018976e-05, -2.0021392848953567e-05,
        -7.3498149390121789e-04,  5.6888820265120965e-05, -3.4651311534116488e-05,
         1.9843604253476427e-05, -5.6383651442472562e-05,
         7.3131479750304169e-05, -8.5125203228748116e-05,  1.4329502961988946e-05,
        -7.2482077225311949e-05, -1.6767541211874010e-05,  8.4369297261769594e-05,
        -1.4202257958371192e-05,  1.6618646595792397e-05,
         1.6327041844856348e-04,  9.2190327694195275e-05, -1.2311099334812739e-05,
        -1.6182058832945405e-04,  3.8034733974350099e-06, -9.1371683906473747e-05,
         1.2201777616986781e-05, -3.7696988144993334e-06,
         5.1551691384217362e-05, -5.7564660526755267e-05,  7.0153713767535825e-06,
        -5.1093915900023203e-05, -6.4107510814715360e-06,  5.7053490288931122e-05,
        -6.9530753600262190e-06,  6.3538240515027896e-06,
         4.0328806418053032e-05,  2.5609342610977432e-05, -2.9607290538320635e-06,
        -3.9970689382718498e-05,  7.0182568367559844e-07, -2.5381933404822401e-05,
         2.9344379828740971e-06, -6.9559351973416766e-07,
         1.1170036525578136e-05, -1.3138092309676044e-05,  1.5550589860992071e-06,
        -1.1070847367246694e-05, -1.6709719218995304e-06,  1.3021426950946600e-05,
        -1.5412501689450531e-06,  1.6561338058245150e-06,
      };
  
      int nvals0 = sizeof(f0m)/sizeof(f0m[0]);
      int nvals1 = sizeof(f1m)/sizeof(f1m[0]);
      for (int i=0; i<nvals0; i++)
        f0[i] = f0m[i];
      for (int i=0; i<nvals1; i++)
        f1[i] = f1m[i];
    }
    else {
      // Upper boundary.
      const double f0m[] =
      {
        // idx[0] = 4;
         5.2735603406759081e-06,  1.5789778675299249e-06,  1.9989826780663522e-06,
         5.2267314864994754e-06,  5.1808733150371001e-07,  1.5649566523489684e-06,
         1.9812318489708070e-06,  5.1348675152920770e-07,
         1.6794598764719933e-05,  4.7243461185578165e-06,  5.4058549896011972e-06,
         1.6645463879424589e-05,  1.3229289098012843e-06,  4.6823942490099290e-06,
         5.3578513680147745e-06,  1.3111813918057051e-06,
         4.4351391469818912e-05,  1.1650043122384951e-05,  1.1710592695890722e-05,
         4.3957554155084436e-05,  2.6752574970649754e-06,  1.1546591538391625e-05,
         1.1606603435837015e-05,  2.6515013939541096e-06,
         9.7120871429426359e-05,  2.3660577068010880e-05,  1.9999449126189980e-05,
         9.6258444751461434e-05,  4.1887526472932372e-06,  2.3450472766235483e-05,
         1.9821855389466126e-05,  4.1515568110403996e-06,
         1.7634549859606911e-04,  3.9539852198925748e-05,  2.6084796724876194e-05,
         1.7477956266190908e-04,  4.8370620921000345e-06,  3.9188741022952382e-05,
         2.5853165518795303e-05,  4.7941093124368999e-06,
         2.6547493073777355e-04,  5.4303954632757974e-05,  2.4003510429954116e-05,
         2.6311753155848884e-04,  3.5659721755435086e-06,  5.3821739239660735e-05,
         2.3790360903442251e-05,  3.5343065871709461e-06,
         3.3130814623591058e-04,  6.1194737241074962e-05,  1.1197518621945012e-05,
         3.2836615261769449e-04,  4.4059942666224667e-07,  6.0651332170970071e-05,
         1.1098085424482478e-05,  4.3668693396879427e-07,
         3.4269905278866620e-04,  5.6458355498447337e-05, -7.3950701740777095e-06,
         3.3965590870142407e-04, -2.9844920494099456e-06,  5.5957009173406704e-05,
        -7.3294024580689704e-06, -2.9579899646808611e-06,
         2.9374366736149800e-04,  4.2514292127503821e-05, -2.2141578694179419e-05,
         2.9113524373960021e-04, -4.8573878417806414e-06,  4.2136768128945003e-05,
        -2.1944962993794169e-05, -4.8142545708542361e-06,
         2.0858381754799318e-04,  2.6012564754977512e-05, -2.6550477975641544e-05,
         2.0673160755237065e-04, -4.5938737655610049e-06,  2.5781574963835798e-05,
        -2.6314711551989231e-05, -4.5530804815602338e-06,
         1.2266163264877528e-04,  1.2843394764617339e-05, -2.1715606284304565e-05,
         1.2157240576270986e-04, -3.0737503418031502e-06,  1.2729346299878326e-05,
        -2.1522773189706808e-05, -3.0464556495587460e-06,
         5.9715709406034303e-05,  5.0593586542458466e-06, -1.3398359967596715e-05,
         5.9185438001676212e-05, -1.5246943425721614e-06,  5.0144318963554475e-06,
        -1.3279383449913500e-05, -1.5111551613701918e-06,
         2.4056337321823497e-05,  1.5575127200312380e-06, -6.4725989939129826e-06,
         2.3842718696134580e-05, -5.6131534500858648e-07,  1.5436821138880390e-06,
        -6.4151227587231497e-06, -5.5633090323863634e-07,
         8.0151330159023820e-06,  3.5840433095651795e-07, -2.4919021860715822e-06,
         7.9439591843807530e-06, -1.4691992411679513e-07,  3.5522172507617869e-07,
        -2.4697742655482486e-06, -1.4561528526607401e-07,
         2.2073750885337045e-06,  5.4194161567241385e-08, -7.7188487204187968e-07,
         2.1877737491242860e-06, -2.2856919929925578e-08,  5.3712921129036765e-08,
        -7.6503058731225592e-07, -2.2653952048423441e-08,
         1.4147068739149328e-06,  3.9678563980892764e-07,  5.3677388121511292e-07,
         1.4021443738919469e-06,  2.1663045016554326e-07,  3.9326221053808074e-07,
         5.3200736596065414e-07,  2.1470678662409674e-07,
        // idx[0] = 12;
         7.4385809600035673e-03, -4.2620301983017916e-04, -4.9868993294229142e-04,
         7.3725268711996930e-03, -1.7454031661270062e-04, -4.2241836624211913e-04,
        -4.9426160053678498e-04, -1.7299041057075798e-04,
         5.3918935380395478e-03, -8.1066940565628881e-04, -6.4417804850088346e-04,
         5.3440138931855462e-03, -4.8296611836812415e-05, -8.0347071692792309e-04,
        -6.3845779160646109e-04, -4.7867741235772842e-05,
         3.2662572681362526e-03, -7.8998871116982251e-04, -5.6119422549588910e-04,
         3.2372531275876946e-03,  5.0252394101218381e-05, -7.8297366559026266e-04,
        -5.5621085925890883e-04,  4.9806156287794315e-05,
         1.6531830617114878e-03, -5.4719482372873017e-04, -3.7015196128374682e-04,
         1.6385029094950227e-03,  8.0981230730795285e-05, -5.4233577121939500e-04,
        -3.6686503725172652e-04,  8.0262122955414322e-05,
         6.9889469196789285e-04, -2.9144378275440866e-04, -1.9257347432051673e-04,
         6.9268855503182584e-04,  6.3595836624062553e-05, -2.8885578195010121e-04,
        -1.9086343507480132e-04,  6.3031109956098697e-05,
         2.4668363462349003e-04, -1.2316540125222183e-04, -8.0544516023116727e-05,
         2.4449310086503576e-04,  3.4617348824041223e-05, -1.2207170093550508e-04,
        -7.9829286244391116e-05,  3.4309949140462771e-05,
         7.4312889120419911e-05, -3.9559856706301649e-05, -2.5345644730926249e-05,
         7.3652995761239434e-05,  1.7079041586906184e-05, -3.9208567891674639e-05,
        -2.5120577143862260e-05,  1.6927380868856262e-05,
         4.9288640555562483e-05,  1.9197719093399832e-05,  9.3898074536981298e-06,
         4.8850960807531919e-05,  1.7734142076997052e-05,  1.9027244664386345e-05,
         9.3064266074409205e-06,  1.7576664111520431e-05,
         1.5102266511218419e-04,  9.5873653213552322e-05,  5.1458764419249045e-05,
         1.4968159420277876e-04,  3.0014361827435697e-05,  9.5022301747812988e-05,
         5.1001814120131087e-05,  2.9747836352724691e-05,
         4.6254048415618476e-04,  2.5918239859045258e-04,  1.3487188316997143e-04,
         4.5843315637685180e-04,  6.6255494613228503e-05,  2.5688087666510847e-04,
         1.3367423009662892e-04,  6.5667150364714081e-05,
         1.1813969572080760e-03,  5.6871985765870109e-04,  2.8575381884705361e-04,
         1.1709062332455834e-03,  1.1130086733197859e-04,  5.6366966432420550e-04,
         2.8321634453203851e-04,  1.1031252326283079e-04,
         2.5178444512609440e-03,  1.0028214075435281e-03,  4.8235931958416520e-04,
         2.4954861652017901e-03,  1.3193983394586925e-04,  9.9391642221575180e-04,
         4.7807600190536245e-04,  1.3076821726856367e-04,
         4.4799350811575659e-03,  1.3955105208072922e-03,  6.3033616055554575e-04,
         4.4401535648606405e-03,  8.3328432799491150e-05,  1.3831184830834593e-03,
         6.2473881867683195e-04,  8.2588482030706125e-05,
         6.6575398175438951e-03,  1.4760068941306469e-03,  5.9457869281011257e-04,
         6.5984213205230388e-03, -4.2899447103834747e-05,  1.4629000541319611e-03,
         5.8929887479281244e-04, -4.2518502955499150e-05,
         8.2661626046304833e-03,  1.0720660909669321e-03,  3.0982027320798264e-04,
         8.1927596475759885e-03, -1.8255117221422333e-04,  1.0625462175989991e-03,
         3.0706909042866552e-04, -1.8093013032390482e-04,
         8.5771895799448498e-03,  3.1197886695786821e-04, -1.2589091675383921e-04,
         8.5010247246790852e-03, -2.3881143669264430e-04,  3.0920851601408308e-04,
        -1.2477301404637767e-04, -2.3669080751195825e-04,
        // idx[0] = 19;
         3.5176447536171117e-03, -2.3790827021664082e-03, -6.1388697056457772e-04,
         3.4864083094370652e-03,  2.8116612144429246e-04, -2.3579566109232482e-03,
        -6.0843569636493593e-04,  2.7866938556759982e-04,
         1.7743608964103485e-03, -1.4228383602125439e-03, -3.9110700695119658e-04,
         1.7586046933318578e-03,  2.5501449371076240e-04, -1.4102036531488738e-03,
        -3.8763400355069862e-04,  2.5274998249492469e-04,
         7.5858237184659968e-04, -6.8983776721579621e-04, -2.0188206214675098e-04,
         7.5184621240645497e-04,  1.6201786760810326e-04, -6.8371205515042432e-04,
        -2.0008936327950031e-04,  1.6057916005455996e-04,
         2.8180957439635491e-04, -2.6245954825183818e-04, -7.9675931947110248e-05,
         2.7930712470157471e-04,  8.6112162611884019e-05, -2.6012892546224052e-04,
        -7.8968415135414035e-05,  8.5347492513266669e-05,
         1.3692552150240488e-04, -2.3628334870108142e-05, -6.7351766555160489e-06,
         1.3570963226150422e-04,  5.5072031551524784e-05, -2.3418516876839040e-05,
        -6.6753687486981578e-06,  5.4582995688061867e-05,
         2.2538939529024113e-04,  1.9086017837583466e-04,  6.0622757673717606e-05,
         2.2338795291675570e-04,  7.2883062312483610e-05,  1.8916535308061394e-04,
         6.0084431742945773e-05,  7.2235865717303348e-05,
         6.0905145791110009e-04,  5.4222206411941353e-04,  1.6846528502596943e-04,
         6.0364312273221280e-04,  1.3353062399761872e-04,  5.3740716937438848e-04,
         1.6696932484790475e-04,  1.3234488121374621e-04,
         1.4749219675761675e-03,  1.1511708985332432e-03,  3.4089330270886968e-04,
         1.4618247616508841e-03,  2.1626062462555974e-04,  1.1409485798989493e-03,
         3.3786619355850205e-04,  2.1434024510954647e-04,
         3.0306885085190582e-03,  1.9748354796544224e-03,  5.6082860635101636e-04,
         3.0037762024011311e-03,  2.4639127781919706e-04,  1.9572990760248154e-03,
         5.5584848678696909e-04,  2.4420334016911550e-04,
         5.3016792692569901e-03,  2.6754726903120783e-03,  7.3909622093513787e-04,
         5.2546007209230902e-03,  1.3556652552434143e-04,  2.6517146762999361e-03,
         7.3253309717877138e-04,  1.3436270407452921e-04,
         7.9236576295588831e-03,  2.7048082887694853e-03,  7.4744339840881547e-04,
         7.8532960931948122e-03, -1.3680114431851166e-04,  2.6807897766548102e-03,
         7.4080615228891281e-04, -1.3558635953847444e-04,
         1.0147502273615196e-02,  1.6812950536897039e-03,  5.0484232418384592e-04,
         1.0057393149318113e-02, -4.4927500173365929e-04,  1.6663652689123294e-03,
         5.0035935896603140e-04, -4.4528547052852176e-04,
         1.1157187397402474e-02, -1.8349902862188801e-04,  6.0616833960680126e-05,
         1.1058112338448070e-02, -5.9620026364057234e-04, -1.8186956983167152e-04,
         6.0078560632073659e-05, -5.9090604618550191e-04,
         1.0538391888353600e-02, -2.0861147480127104e-03, -4.1035235387461852e-04,
         1.0444811691083959e-02, -4.6474412137482361e-04, -2.0675901921113600e-03,
        -4.0670845311324271e-04, -4.6061722544811013e-04,
         8.5438428607558167e-03, -3.1770890503348258e-03, -7.1469750720547990e-04,
         8.4679741220696093e-03, -1.4653149000122957e-04, -3.1488767174453918e-03,
        -7.0835104235382441e-04, -1.4523030041881422e-04,
         5.9328259833403407e-03, -3.1592831855465011e-03, -7.6447171990721248e-04,
         5.8801428954679841e-03,  1.4953722986155545e-04, -3.1312289675153003e-03,
        -7.5768326346016773e-04,  1.4820934951530792e-04,
        // idx[0] = 28;
         3.1225316368462653e-04, -2.6354343092504693e-04, -2.5533068601899983e-05,
         3.0948037700471851e-04,  6.4354644170239486e-06, -2.6120318333164132e-04,
        -2.5306336703714858e-05,  6.3783179343302755e-06,
         2.2070902006211341e-04, -2.1925837439635602e-04, -2.6427797458292739e-05,
         2.1874913909968421e-04,  1.7877051979852441e-05, -2.1731137506795673e-04,
        -2.6193120429221325e-05,  1.7718304984238919e-05,
         1.3604357772040071e-04, -1.5068140809657886e-04, -2.1787568445674219e-05,
         1.3483552007980261e-04,  2.0478489031820037e-05, -1.4934336752606933e-04,
        -2.1594096332018494e-05,  2.0296641453585649e-05,
         7.3131479750304263e-05, -8.5125203228748197e-05, -1.4329502961952065e-05,
         7.2482077225312057e-05,  1.6767541211912712e-05, -8.4369297261769675e-05,
        -1.4202257958334635e-05,  1.6618646595830757e-05,
         3.6725180241487518e-05, -3.6057780516486069e-05, -6.7954495238028766e-06,
         3.6399063159472418e-05,  1.1694827060935900e-05, -3.5737589898260231e-05,
        -6.7351063980479749e-06,  1.1590977798639422e-05,
         2.4837446165027385e-05, -9.6145358765581644e-07, -1.9699547936173084e-07,
         2.4616891346377475e-05,  9.0936978967432890e-06, -9.5291594573173095e-07,
        -1.9524617301449969e-07,  9.0129464830453212e-06,
         3.4938222906600408e-05,  3.1273776945169290e-05,  6.0324138429933910e-06,
         3.4627973883173628e-05,  1.0006175581772412e-05,  3.0996067950580641e-05,
         5.9788464217568256e-06,  9.9173214288070243e-06,
         6.7058946442705397e-05,  7.0433021024861873e-05,  1.2563048496794633e-05,
         6.6463467597044124e-05,  1.2699023524377309e-05,  6.9807580628297399e-05,
         1.2451489487688373e-05,  1.2586257066350973e-05,
         1.2172339077131295e-04,  1.1660274450623232e-04,  1.8887207777269423e-05,
         1.2064249540877686e-04,  1.3462333943158600e-05,  1.1556732013136789e-04,
         1.8719490667449008e-05,  1.3342789340959595e-05,
         1.9535178824326663e-04,  1.5650334125842918e-04,  2.3209649131246132e-05,
         1.9361707776044762e-04,  8.6672787667954396e-06,  1.5511360232069738e-04,
         2.3003549038625602e-05,  8.5903139257285476e-06,
         2.7796534980573806e-04,  1.6786465218620092e-04,  2.3905953966723439e-05,
         2.7549703656169435e-04, -2.9974039875562224e-06,  1.6637402558656333e-04,
         2.3693670734914863e-05, -2.9707872457018826e-06,
         3.5577914498131474e-04,  1.3094475217078909e-04,  2.0440326610403722e-05,
         3.5261985057240491e-04, -1.8525522986956615e-05,  1.2978197175140562e-04,
         2.0258817911854427e-05, -1.8361017613270647e-05,
         4.1504347649095046e-04,  4.2960315696147001e-05,  1.3312487632465493e-05,
         4.1135791888244821e-04, -3.1477256709126155e-05,  4.2578831038885878e-05,
         1.3194273655229314e-05, -3.1197740828187762e-05,
         4.4457679107663076e-04, -7.5578664486840579e-05,  3.4115616361105143e-06,
         4.4062897965993503e-04, -3.5399985354919482e-05, -7.4907531129210335e-05,
         3.3812672027389913e-06, -3.5085635912618252e-05,
         4.3680260689641108e-04, -1.8748123385802206e-04, -8.0091734113102738e-06,
         4.3292382970209987e-04, -2.7606287551927432e-05, -1.8581641335839935e-04,
        -7.9380525006687252e-06, -2.7361145611642716e-05,
         3.9008280658611196e-04, -2.5603798505596852e-04, -1.8683093211797680e-05,
         3.8661889801461854e-04, -1.1175426823779215e-05, -2.5376438530717394e-04,
        -1.8517188625320945e-05, -1.1076189799969395e-05,
      };
  
      const double f1m[] =
      {
        // idx[1] = 2;
         1.0064242557051679e-06,  2.4638105513050012e-07,  3.0952299929100151e-07,
         9.9748727733279151e-07,  7.1551685415696771e-08,  2.4419320825719352e-07,
         3.0677445628368039e-07,  7.0916311356057116e-08,
         2.4501722734192372e-06,  6.2869377552759965e-07,  7.5054342119427020e-07,
         2.4284149116585134e-06,  2.0047632122023783e-07,  6.2311101791530323e-07,
         7.4387864708462476e-07,  1.9869610523601224e-07,
         6.2693514816208999e-06,  1.6983388600546720e-06,  2.0605121390263682e-06,
         6.2136800704019623e-06,  5.9893969448756465e-07,  1.6832577274456009e-06,
         2.0422149325370482e-06,  5.9362115107443639e-07,
         1.6794598764719933e-05,  4.7243461185578165e-06,  5.4058549896011972e-06,
         1.6645463879424589e-05,  1.3229289098012843e-06,  4.6823942490099290e-06,
         5.3578513680147745e-06,  1.3111813918057051e-06,
         4.6722524149225004e-05,  1.3757882849431861e-05,  1.2559267568821096e-05,
         4.6307631338904288e-05,  3.6001849176655706e-06,  1.3635713793213210e-05,
         1.2447742134095955e-05,  3.5682155224891871e-06,
         1.3509589455886354e-04,  4.0376357161087206e-05,  3.4643122769977236e-05,
         1.3389625227974990e-04,  8.2223573880571182e-06,  4.0017817878415452e-05,
         3.4335494215526362e-05,  8.1493434183217232e-06,
         3.8840901300543125e-04,  1.1420438032971131e-04,  8.3157991770278211e-05,
         3.8495996760614943e-04,  2.2616417530061213e-05,  1.1319025326426001e-04,
         8.2419554506144016e-05,  2.2415585293376721e-05,
         1.0598521823626484e-03,  2.8756536278117851e-04,  1.7461178578371693e-04,
         1.0504407676655194e-03,  3.4239520152090462e-05,  2.8501180208026012e-04,
         1.7306124510043641e-04,  3.3935475561208271e-05,
         2.5465331961989765e-03,  5.7389435441791746e-04,  2.7509217369229086e-04,
         2.5239201560521818e-03,  1.6496079388288131e-05,  5.6879821190705530e-04,
         2.7264937405509861e-04,  1.6349595334583635e-05,
         4.8835866960519232e-03,  7.1615900951608373e-04,  2.0939718370369375e-04,
         4.8402207810962460e-03, -6.7344077386474613e-05,  7.0979956662413596e-04,
         2.0753775107238749e-04, -6.6746066597586618e-05,
         6.5963664604376145e-03,  1.6275731066703455e-04, -2.2806276733301181e-04,
         6.5377911786327313e-03, -1.6368482665102420e-04,  1.6131203690983743e-04,
        -2.2603758559912031e-04, -1.6223131661549542e-04,
         5.3918935380395478e-03, -8.1066940565628881e-04, -6.4417804850088346e-04,
         5.3440138931855462e-03, -4.8296611836812415e-05, -8.0347071692792309e-04,
        -6.3845779160646109e-04, -4.7867741235772842e-05,
         2.2325143698531066e-03, -8.4647335003058420e-04, -4.7667889076579453e-04,
         2.2126897953495669e-03,  1.1597626131718698e-04, -8.3895672473152899e-04,
        -4.7244601490534796e-04,  1.1494640007007816e-04,
         3.7918019032058645e-04, -2.4429110051544058e-04, -1.1791962244750706e-04,
         3.7581309623385434e-04,  6.5545377999984785e-05, -2.4212181229579129e-04,
        -1.1687250428683259e-04,  6.4963339538297809e-05,
         1.3265773496585429e-04,  9.0192913254221100e-05,  3.6330610507766499e-05,
         1.3147974337672313e-04,  4.2769662240897072e-05,  8.9392006369748384e-05,
         3.6007997177928130e-05,  4.2389870573244011e-05,
         2.0740722732800190e-03,  1.3049423409384143e-03,  4.7737679678530820e-04,
         2.0556546537284495e-03,  2.3524506169704339e-04,  1.2933545424408577e-03,
         4.7313772356726925e-04,  2.3315610168174477e-04,
         1.0158496318663347e-02,  2.7739047335627711e-03,  8.7039930616657365e-04,
         1.0068289568000181e-02, -2.0057911001830609e-04,  2.7492726497568079e-03,
         8.6267022001781804e-04, -1.9879798127660263e-04,
         1.0744903537459927e-02, -2.7332594355407298e-03, -6.8024602402668875e-04,
         1.0649489530908321e-02, -3.3891226168725110e-04, -2.7089882791938406e-03,
        -6.7420548598305674e-04, -3.3590274404529956e-04,
         1.7743608964103485e-03, -1.4228383602125439e-03, -3.9110700695119658e-04,
         1.7586046933318578e-03,  2.5501449371076240e-04, -1.4102036531488738e-03,
        -3.8763400355069862e-04,  2.5274998249492469e-04,
         1.8587668557072320e-04,  1.4266611485212182e-04,  4.2787284932814561e-05,
         1.8422611334985533e-04,  6.5393262891317997e-05,  1.4139924953595248e-04,
         4.2407336776865824e-05,  6.4812575201364270e-05,
         3.7794409682195629e-03,  2.0869073873386228e-03,  4.9681034418490398e-04,
         3.7458797916072352e-03,  8.9354392597548640e-05,  2.0683757928544161e-03,
         4.9239870239866507e-04,  8.8560931719003469e-05,
         5.3589500608310107e-03, -1.9322998876550560e-03, -2.6621995559220085e-04,
         5.3113629517954333e-03, -2.1708287121334914e-04, -1.9151411971654085e-03,
        -2.6385594064330117e-04, -2.1515519020402394e-04,
         3.5447984105491342e-04, -3.8297317940784757e-04, -8.6790399124008177e-05,
         3.5133208437576662e-04,  9.1077196943933081e-05, -3.7957240383813508e-04,
        -8.6019706331675079e-05,  9.0268437680937836e-05,
         7.7591183917630939e-04,  6.6057430431279971e-04,  1.3512321260020200e-04,
         7.6902179525469223e-04,  8.3115590710258533e-05,  6.5470844979118492e-04,
         1.3392332773875466e-04,  8.2377529964629514e-05,
         1.9625514353073274e-03, -6.9400771096554106e-04, -5.0835475579179123e-05,
         1.9451241131491104e-03, -1.1538337848761807e-04, -6.8784497008563667e-04,
        -5.0384060042217249e-05, -1.1435878199937753e-04,
         8.4095026418060529e-05, -8.6935345625212553e-05, -1.7840996781592987e-05,
         8.3348268350513742e-05,  2.6621237831636578e-05, -8.6163365717888856e-05,
        -1.7682569953667403e-05,  2.6384843065315550e-05,
         5.0821419007155365e-04,  2.7840870285490151e-04,  4.2754573047174882e-05,
         5.0370128291589108e-04, -9.1435198883331253e-06,  2.7593645266618531e-04,
         4.2374915370527306e-05, -9.0623260587665063e-06,
         2.2070902006211341e-04, -2.1925837439635602e-04, -2.6427797458292739e-05,
         2.1874913909968421e-04,  1.7877051979852441e-05, -2.1731137506795673e-04,
        -2.6193120429221325e-05,  1.7718304984238919e-05,
         7.8376321342582138e-05,  7.3618252299111100e-05,  1.0968129522333047e-05,
         7.7680344983929785e-05,  7.1081816882799169e-06,  7.2964527267267229e-05,
         1.0870733284344096e-05,  7.0450615223509237e-06,
         1.0129896530548440e-04, -8.5306127837966829e-05, -6.3059363521162601e-06,
         1.0039943744042426e-04,  1.7577593965507822e-07, -8.4548615272865850e-05,
        -6.2499400697562340e-06,  1.7421506136525714e-07,
         1.9919275589572257e-05,  1.9647328326244563e-05,  2.6426654263499524e-06,
         1.9742393788355451e-05,  1.9320906579101043e-06,  1.9472861398076579e-05,
         2.6191987385919478e-06,  1.9149338253662459e-06,
         2.2776817069048982e-05, -2.1810758063425291e-05, -1.5443168402192109e-06,
         2.2574560495458119e-05,  4.7017697371218268e-07, -2.1617080027554231e-05,
        -1.5306034125837987e-06,  4.6600183442928546e-07,
        // idx[1] = 7;
         3.9161987831707655e-05,  8.5585033279271802e-06,  3.0940904155156180e-06,
         3.8814232065402117e-05,  6.2790787320463302e-07,  8.4825044052978738e-06,
         3.0666151048115204e-06,  6.2233209435106859e-07,
         8.2966973231524767e-05,  1.7351396003363978e-05,  5.8734092978907861e-06,
         8.2230232198915074e-05,  9.6285893253930273e-07,  1.7197316796773362e-05,
         5.8212538261105722e-06,  9.5430881124903129e-07,
         1.6948528988022310e-04,  3.3617119358933856e-05,  9.1922219423802619e-06,
         1.6798027212900225e-04,  8.7142547168074894e-07,  3.3318601644412088e-05,
         9.1105956419135583e-06,  8.6368727325260448e-07,
         3.3130814623591058e-04,  6.1194737241074962e-05,  1.1197518621945012e-05,
         3.2836615261769449e-04,  4.4059942666224667e-07,  6.0651332170970071e-05,
         1.1098085424482478e-05,  4.3668693396879427e-07,
         6.0949868976399829e-04,  1.0001588819084109e-04,  9.8223350700890946e-06,
         6.0408638319692741e-04, -3.3693722410864667e-06,  9.9127753962568473e-05,
         9.7351134082600786e-06, -3.3394524466495121e-06,
         1.0226479509198830e-03,  1.3603218361846954e-04, -2.0902400504305411e-05,
         1.0135669072465837e-03, -1.1408417045065111e-05,  1.3482422715671509e-04,
        -2.0716788621266325e-05, -1.1307111084068212e-05,
         1.4988806244302360e-03,  1.3106170460485647e-04, -8.5695676969186599e-05,
         1.4855706672751183e-03, -3.0453476425098579e-05,  1.2989788565587884e-04,
        -8.4934705234514663e-05, -3.0183051642873811e-05,
         1.8030298257953232e-03,  3.5581460837068996e-05, -1.8651702121551744e-04,
         1.7870190445899432e-03, -3.0633437103584870e-05,  3.5265499904931559e-05,
        -1.8486076285803674e-04, -3.0361414282876413e-05,
         1.6413216149696493e-03, -1.2558118114409342e-04, -2.6036761791839430e-04,
         1.6267467915868209e-03, -5.4223079993046918e-06, -1.2446602886760557e-04,
        -2.5805557132669790e-04, -5.3741582761204354e-06,
         1.0189032686435142e-03, -2.1096043183989230e-04, -2.2864831174279596e-04,
         1.0098554774920273e-03,  2.7657765772379680e-05, -2.0908711767234027e-04,
        -2.2661793041470085e-04,  2.7412166709027565e-05,
         3.7915828230823055e-04, -1.4271127978161298e-04, -1.1053775807032425e-04,
         3.7579138276314572e-04,  3.4225960821833308e-05, -1.4144401340396707e-04,
        -1.0955619035908770e-04,  3.3922036636873959e-05,
         7.4312889120419911e-05, -3.9559856706301649e-05, -2.5345644730926249e-05,
         7.3652995761239434e-05,  1.7079041586906184e-05, -3.9208567891674639e-05,
        -2.5120577143862260e-05,  1.6927380868856262e-05,
         1.6452199128976321e-04,  1.1354492403363433e-04,  5.6206288820234808e-05,
         1.6306104728965398e-04,  3.6099553950453624e-05,  1.1253665289486821e-04,
         5.5707180830010552e-05,  3.5778992386998920e-05,
         1.7045525792129461e-03,  9.4638345870692116e-04,  4.0697719285157964e-04,
         1.6894162692038725e-03,  1.8078271052515243e-04,  9.3797963849442793e-04,
         4.0336326328862745e-04,  1.7917737245335601e-04,
         8.0791678878996687e-03,  2.6080011965470950e-03,  9.4517642723753457e-04,
         8.0074254311060748e-03,  2.5821417318181404e-05,  2.5848423247725465e-03,
         9.3678332538171393e-04,  2.5592124903170118e-05,
         1.3356481058148951e-02, -4.3293862331202017e-04, -1.8456831679660113e-04,
         1.3237876422309669e-02, -5.5117088081241728e-04, -4.2909415802695381e-04,
        -1.8292936280072232e-04, -5.4627652118884712e-04,
         5.8131600430988686e-03, -2.8057188210851105e-03, -8.7511944457968171e-04,
         5.7615395805694384e-03,  1.9564030185384277e-04, -2.7808042303636878e-03,
        -8.6734844392551505e-04,  1.9390302938994804e-04,
         4.7634326079316609e-04, -4.2920915428436561e-04, -1.4088268702256988e-04,
         4.7211336530385978e-04,  1.1855725319557678e-04, -4.2539780642850768e-04,
        -1.3963165842322361e-04,  1.1750447291758676e-04,
         6.0905145791110009e-04,  5.4222206411941353e-04,  1.6846528502596943e-04,
         6.0364312273221280e-04,  1.3353062399761872e-04,  5.3740716937438848e-04,
         1.6696932484790475e-04,  1.3234488121374621e-04,
         6.4242974589046650e-03,  2.4764802964899277e-03,  6.3518193397627161e-04,
         6.3672501380332945e-03, -8.6162573411864465e-05,  2.4544893212884809e-03,
         6.2954156196189386e-04, -8.5397455668800133e-05,
         6.1192473190725622e-03, -2.6624998539887487e-03, -4.5719660344076035e-04,
         6.0649088225232552e-03, -1.3024981958776287e-04, -2.6388570378735108e-03,
        -4.5313672895570958e-04, -1.2909321012205950e-04,
         3.6469519543160179e-04, -3.8363095857578205e-04, -9.4090664425681766e-05,
         3.6145672710613579e-04,  9.7469050278256532e-05, -3.8022434197216740e-04,
        -9.3255145778105395e-05,  9.6603531795991395e-05,
         9.1178047998913487e-04,  7.8751323767826042e-04,  1.7362275060939670e-04,
         9.0368393185466219e-04,  1.1501287471411306e-04,  7.8052017413354356e-04,
         1.7208099249063841e-04,  1.1399156827397299e-04,
         3.0428994134042461e-03, -6.0542157405543478e-04, -2.7297328115239033e-05,
         3.0158786752883348e-03, -2.0376965733478786e-04, -6.0004547199625686e-04,
        -2.7054929713564356e-05, -2.0196019675171445e-04,
         2.4486683302324668e-04, -2.7365088383469040e-04, -5.4183415554647035e-05,
         2.4269243233841054e-04,  5.7765593385118677e-05, -2.7122088275259332e-04,
        -5.3702270540297137e-05,  5.7252638877291795e-05,
         4.8684120652024178e-04,  3.9028286351599020e-04,  6.7232236617606201e-05,
         4.8251808999280420e-04,  3.0667260129115022e-05,  3.8681717845267898e-04,
         6.6635218967075661e-05,  3.0394936962259377e-05,
         6.3412598830576864e-04, -4.3372366375867385e-04, -3.8776407472139081e-05,
         6.2849499301653093e-04, -9.5143346916887677e-06, -4.2987222736827886e-04,
        -3.8432075632984299e-05, -9.4298480520979985e-06,
         3.4938222906600408e-05,  3.1273776945169290e-05,  6.0324138429933910e-06,
         3.4627973883173628e-05,  1.0006175581772412e-05,  3.0996067950580641e-05,
         5.9788464217568256e-06,  9.9173214288070243e-06,
         2.4392083227262588e-04, -5.6676719861445929e-05,  1.1626687511291100e-06,
         2.4175483200958067e-04, -1.9445449116183467e-05, -5.6173434475836811e-05,
         1.1523443323522607e-06, -1.9272774861556413e-05,
         1.0411319188242222e-05, -3.0681569078749064e-06, -5.5818318185091001e-07,
         1.0318867387835211e-05,  3.5693133209598589e-06, -3.0409118849401984e-06,
        -5.5322655347500679e-07,  3.5376180634451215e-06,
         6.2185931257129586e-05, -1.0179624868692354e-05,  8.9508031604239662e-07,
         6.1633724452135960e-05, -5.2188375294793664e-06, -1.0089230497953875e-05,
         8.8713206421851195e-07, -5.1724946101135919e-06,
         3.4821562995004494e-06,  1.0893657152118484e-06,  1.5562699789901857e-07,
         3.4512349903592412e-06,  1.0573421757814150e-06,  1.0796922223669888e-06,
         1.5424504083022067e-07,  1.0479530497705965e-06,
        // idx[1] = 10;
         3.6490020128125022e-05,  7.3675063786135918e-06, -3.2021281911734263e-06,
         3.6165991251789828e-05, -6.9414808136536902e-07,  7.3020834272183270e-06,
        -3.1736935124304748e-06, -6.8798409400593995e-07,
         7.1180658041772459e-05,  1.2900664940831099e-05, -6.7559365584218519e-06,
         7.0548578679769091e-05, -1.4100360605954828e-06,  1.2786107920854769e-05,
        -6.6959443050898310e-06, -1.3975150370744643e-06,
         1.2807110995553264e-04,  2.0058649867687895e-05, -1.3810778814976504e-05,
         1.2693384728195280e-04, -2.7419450947362422e-06,  1.9880530432439396e-05,
        -1.3688140046211276e-05, -2.7175968103313399e-06,
         2.0858381754799318e-04,  2.6012564754977512e-05, -2.6550477975641544e-05,
         2.0673160755237065e-04, -4.5938737655610049e-06,  2.5781574963835798e-05,
        -2.6314711551989231e-05, -4.5530804815602338e-06,
         2.9871110310366873e-04,  2.5205913165246944e-05, -4.5389429729339090e-05,
         2.9605856899303593e-04, -5.9426328345682074e-06,  2.4982086384903736e-05,
        -4.4986374706046403e-05, -5.8898626625294748e-06,
         3.6469434674854192e-04,  1.1237385446410865e-05, -6.5265128683820238e-05,
         3.6145588595932343e-04, -6.3539671552322184e-06,  1.1137598234281054e-05,
        -6.4685578816841917e-05, -6.2975443626344209e-06,
         3.5893595493224776e-04, -1.4493122690417814e-05, -7.7817247724849183e-05,
         3.5574862826746054e-04,  1.6084567762459325e-06, -1.4364424754832272e-05,
        -7.7126235901582319e-05,  1.5941737903771684e-06,
         2.6694095601170731e-04, -3.6280004828176018e-05, -7.0464806864673933e-05,
         2.6457053862852676e-04,  7.8696688227861677e-06, -3.5957840873301966e-05,
        -6.9839084212031750e-05,  7.7997867033238258e-06,
         1.3743518008194373e-04, -3.5616161618941401e-05, -4.3245713680204564e-05,
         1.3621476510780801e-04,  1.0173131963281794e-05, -3.5299892546240518e-05,
        -4.2861694708416988e-05,  1.0082795249098096e-05,
         4.4281012350574168e-05, -1.7888520118666163e-05, -1.6654315330492159e-05,
         4.3887799997592119e-05,  6.9220705758272029e-06, -1.7729671286766401e-05,
        -1.6506426152472439e-05,  6.8606030638137037e-06,
         5.1438549738468456e-05,  2.9665824357334248e-05,  1.7900189835036316e-05,
         5.0981778944329555e-05,  1.3973723052796347e-05,  2.9402393871455861e-05,
         1.7741237376855478e-05,  1.3849637350373645e-05,
         4.6254048415618476e-04,  2.5918239859045258e-04,  1.3487188316997143e-04,
         4.5843315637685180e-04,  6.6255494613228503e-05,  2.5688087666510847e-04,
         1.3367423009662892e-04,  6.5667150364714081e-05,
         2.8294587601195094e-03,  1.2375254640891551e-03,  5.5518378989409295e-04,
         2.8043333603674919e-03,  1.6372300936085449e-04,  1.2265363228347241e-03,
         5.5025379591307449e-04,  1.6226916026548170e-04,
         9.1957619553724862e-03,  2.1739429633878351e-03,  7.9765225550964614e-04,
         9.1141042198332120e-03, -1.2436740398918847e-04,  2.1546384989570141e-03,
         7.9056915818189724e-04, -1.2326302997059803e-04,
         1.2500828315738304e-02, -7.8508011247854954e-04, -3.9877045177897329e-04,
         1.2389821817572911e-02, -4.5450998846959987e-04, -7.7810865491874858e-04,
        -3.9522939751394622e-04, -4.5047397094125129e-04,
         5.4790816862939342e-03, -2.4387569363537231e-03, -8.7033830325416366e-04,
         5.4304278166625993e-03,  1.8988685124218666e-04, -2.4171009420033036e-03,
        -8.6260975880708712e-04,  1.8820066902516417e-04,
         5.7215777750935844e-04, -4.7524252007977803e-04, -1.6913878684937368e-04,
         5.6707705568235394e-04,  1.2685880911773336e-04, -4.7102239909252029e-04,
        -1.6763684602130442e-04,  1.2573231159245635e-04,
         3.7427031371474034e-04,  3.3208117201649276e-04,  1.1114484687927514e-04,
         3.7094681899556005e-04,  9.6463580458143968e-05,  3.2913231398235741e-04,
         1.1015788826104790e-04,  9.5606990478928455e-05,
         5.3016792692569901e-03,  2.6754726903120783e-03,  7.3909622093513787e-04,
         5.2546007209230902e-03,  1.3556652552434143e-04,  2.6517146762999361e-03,
         7.3253309717877138e-04,  1.3436270407452921e-04,
         9.1961936544017117e-03, -1.7368826364828100e-03, -2.6809009052276823e-04,
         9.1145320854046606e-03, -4.6206186478537920e-04, -1.7214592377823002e-03,
        -2.6570946890393444e-04, -4.5795878711323652e-04,
         1.4411213884148813e-03, -1.2810931792856986e-03, -2.9718607020745780e-04,
         1.4283243293145833e-03,  2.1624336633407510e-04, -1.2697171597080966e-03,
        -2.9454707828436010e-04,  2.1432314007050667e-04,
         2.5492665841420016e-04,  2.4888323947485432e-04,  6.6278567407822636e-05,
         2.5266292716977223e-04,  7.1894578026384081e-05,  2.4667317337616389e-04,
         6.5690018274474520e-05,  7.1256159103876392e-05,
         3.4578494082033386e-03,  1.1303046461231766e-03,  2.4223248931237062e-04,
         3.4271439425899376e-03, -1.3900706598036629e-04,  1.1202676183793233e-03,
         2.4008148141902978e-04, -1.3777269276724882e-04,
         1.5781234722698934e-03, -1.2618335200808048e-03, -1.9187865723666297e-04,
         1.5641098440602533e-03,  8.2967902797066034e-05, -1.2506285249561588e-03,
        -1.9017478791900309e-04,  8.2231153509978078e-05,
         9.6977901902638622e-05,  8.6316756686137351e-05,  1.9677860154709606e-05,
         9.6116744784267068e-05,  3.0120448337405375e-05,  8.5550269806171266e-05,
         1.9503122105998658e-05,  2.9852980821761017e-05,
         1.1828731451272081e-03,  1.7396447751082816e-04,  4.5145997322302639e-05,
         1.1723693127172275e-03, -8.1902648559895069e-05,  1.7241968488062006e-04,
         4.4745103962090367e-05, -8.1175358657380482e-05,
         1.8720681143410254e-04, -2.0805519933339656e-04, -3.3414378701659159e-05,
         1.8554442778678492e-04,  3.2611926466714527e-05, -2.0620768343126817e-04,
        -3.3117661310270503e-05,  3.2322334796141819e-05,
         1.9535178824326663e-04,  1.5650334125842918e-04,  2.3209649131246132e-05,
         1.9361707776044762e-04,  8.6672787667954396e-06,  1.5511360232069738e-04,
         2.3003549038625602e-05,  8.5903139257285476e-06,
         1.6951067600112074e-04, -1.5119863365618903e-04, -1.3863877085534901e-05,
         1.6800543282288725e-04,  4.5044095319706036e-06, -1.4985600015818035e-04,
        -1.3740766807768502e-05,  4.4644106842288668e-06,
         3.6276567298503925e-05,  3.6027257115504003e-05,  5.2295393006317666e-06,
         3.5954433868657657e-05,  3.9956008741687153e-06,  3.5707337542985657e-05,
         5.1831013502720623e-06,  3.9601202124153117e-06,
         5.1575023302046816e-05, -4.4579128212969987e-05, -2.9869584382357371e-06,
         5.1117040632022390e-05,  3.2845211634800542e-09, -4.4183268611573997e-05,
        -2.9604344521430178e-06,  3.2553548407936682e-09,
         1.2401275268881753e-05,  1.1336295347525348e-05,  1.3473909962804335e-06,
         1.2291152794945326e-05,  6.9063884390899995e-07,  1.1235629822256665e-05,
         1.3354262566344723e-06,  6.8450601833752906e-07,
        // idx[1] = 13;
         6.2657435059192389e-06,  1.1641163689362830e-06, -1.6403632698369783e-06,
         6.2101041332771099e-06, -3.2689517205055183e-07,  1.1537790954125808e-06,
        -1.6257969564932065e-06, -3.2399236534047043e-07,
         1.1267655777296214e-05,  1.7224008808195605e-06, -3.0785713582669099e-06,
         1.1167599766704021e-05, -4.9786387303821501e-07,  1.7071060791162162e-06,
        -3.0512338557269196e-06, -4.9344287598800761e-07,
         1.7792497642738746e-05,  1.9998167603756296e-06, -4.8547004491311306e-06,
         1.7634501483840234e-05, -4.9857845221345014e-07,  1.9820585246863783e-06,
        -4.8115910420671554e-06, -4.9415110974917022e-07,
         2.4056337321823497e-05,  1.5575127200312380e-06, -6.4725989939129826e-06,
         2.3842718696134580e-05, -5.6131534500858648e-07,  1.5436821138880390e-06,
        -6.4151227587231497e-06, -5.5633090323863634e-07,
         2.7164067597516312e-05,  1.1098766716295953e-07, -8.8660114789262803e-06,
         2.6922852540100311e-05, -3.1627341965558876e-07,  1.1000210428983373e-07,
        -8.7872819050043051e-06, -3.1346493337835670e-07,
         2.3994825830339474e-05, -1.6841247767798238e-06, -7.7539596815715105e-06,
         2.3781753422476620e-05,  2.6317833456926379e-07, -1.6691698642554609e-06,
        -7.6851050513480680e-06,  2.6084132900645678e-07,
         1.6016882937583991e-05, -2.8230768932199529e-06, -6.0290454073971664e-06,
         1.5874654115499630e-05,  1.3614222825750313e-06, -2.7980081640083472e-06,
        -5.9755078976376938e-06,  1.3493329460690080e-06,
         9.8455354864211670e-06,  2.9617890833854716e-07, -6.9489101822206762e-07,
         9.7581078064861676e-06,  3.1357290448382672e-06,  2.9354885994377525e-07,
        -6.8872043363429996e-07,  3.1078839859612674e-06,
         3.9011083818695911e-05,  1.9902575296058611e-05,  1.4677088296293310e-05,
         3.8664668069677257e-05,  7.5064936588684329e-06,  1.9725841792302946e-05,
         1.4546756758743548e-05,  7.4398365099559416e-06,
         2.3888049043531946e-04,  1.1483546140429011e-04,  7.1312268443981844e-05,
         2.3675924806213432e-04,  2.9643345516891777e-05,  1.1381573038217466e-04,
         7.0679020390632398e-05,  2.9380114641574981e-05,
         1.2234867187375981e-03,  5.1387522103972256e-04,  2.8134955828176556e-04,
         1.2126222405791434e-03,  9.3470469345815967e-05,  5.0931204431728318e-04,
         2.7885119349853603e-04,  9.2640457988017256e-05,
         4.4799350811575659e-03,  1.3955105208072922e-03,  6.3033616055554575e-04,
         4.4401535648606405e-03,  8.3328432799491150e-05,  1.3831184830834593e-03,
         6.2473881867683195e-04,  8.2588482030706125e-05,
         9.7716297507283792e-03,  1.3456395734362109e-03,  4.2230678771116898e-04,
         9.6848583486580462e-03, -2.3479472920207153e-04,  1.3336903863049993e-03,
         4.1855673239713247e-04, -2.3270976810844599e-04,
         1.0306703041584746e-02, -1.1825534137556507e-03, -6.7086169125845014e-04,
         1.0215180225385718e-02, -2.6682678254561637e-04, -1.1720524205383232e-03,
        -6.6490448544624355e-04, -2.6445737901498439e-04,
         4.1568120046519205e-03, -1.8174644792154824e-03, -7.5296484661487636e-04,
         4.1198997991153197e-03,  1.9040708381710860e-04, -1.8013255192776273e-03,
        -7.4627857041355196e-04,  1.8871628197049720e-04,
         4.8860692316653393e-04, -3.8023488046055041e-04, -1.4753883724741325e-04,
         4.8426812719636612e-04,  1.0196923908055593e-04, -3.7685842079770358e-04,
        -1.4622870249053404e-04,  1.0106375923033906e-04,
         2.9294848655103347e-04,  2.4980931921513331e-04,  9.0157589311276736e-05,
         2.9034712408020091e-04,  8.0321792535837242e-05,  2.4759102959185763e-04,
         8.9356996100994944e-05,  7.9608540526404764e-05,
         4.4254596134883665e-03,  2.4549067862642271e-03,  7.3848760727670605e-04,
         4.3861618356085411e-03,  2.3344588158139668e-04,  2.4331073823541218e-03,
         7.3192988796788876e-04,  2.3137289816214148e-04,
         1.1157187397402474e-02, -1.8349902862188801e-04,  6.0616833960680126e-05,
         1.1058112338448070e-02, -5.9620026364057234e-04, -1.8186956983167152e-04,
         6.0078560632073659e-05, -5.9090604618550191e-04,
         3.4746815254179230e-03, -2.4230957949160658e-03, -5.5706413873520704e-04,
         3.4438265917579652e-03,  2.4825618773810390e-04, -2.4015788704276086e-03,
        -5.5211744738542956e-04,  2.4605168981583451e-04,
         1.2967070437580762e-04, -3.0372989641602543e-05, -8.4366829497824695e-06,
         1.2851923741347230e-04,  5.2883604796430712e-05, -3.0103279576495651e-05,
        -8.3617657837569207e-06,  5.2414002012475693e-05,
         2.3207575119290032e-03,  1.5289575356248154e-03,  3.4582312161757889e-04,
         2.3001493443753429e-03,  1.1872552071725908e-04,  1.5153804975608961e-03,
         3.4275223601337767e-04,  1.1767124623521591e-04,
         4.0562948476909767e-03, -1.4421644508333046e-03, -1.5972876124242848e-04,
         4.0202752276149523e-03, -1.9704447144722474e-04, -1.4293581294102611e-03,
        -1.5831038079643057e-04, -1.9529473005363505e-04,
         2.3072328684209585e-04, -2.5397968959006082e-04, -5.5574896425858917e-05,
         2.2867447987741607e-04,  6.3533756465527876e-05, -2.5172436736385032e-04,
        -5.5081395156797248e-05,  6.2969581065114454e-05,
         7.4703445644323579e-04,  5.8308309956022755e-04,  1.0745158404378290e-04,
         7.4040084169994325e-04,  4.8758475781811738e-05,  5.7790536162264500e-04,
         1.0649742134625928e-04,  4.8325503860614624e-05,
         1.1026100056389333e-03, -6.4967562903913128e-04, -5.9942438906131847e-05,
         1.0928189044033454e-03, -3.1102503922783453e-05, -6.4390655400654037e-04,
        -5.9410154159361466e-05, -3.0826315820888077e-05,
         3.6855264730397853e-05,  9.9690149435174425e-06,  2.1646015699630160e-06,
         3.6527992507042120e-05,  1.3351157106268053e-05,  9.8804907744713724e-06,
         2.1453800564651972e-06,  1.3232599746756376e-05,
         4.1504347649095046e-04,  4.2960315696147001e-05,  1.3312487632465493e-05,
         4.1135791888244821e-04, -3.1477256709126155e-05,  4.2578831038885878e-05,
         1.3194273655229314e-05, -3.1197740828187762e-05,
         4.0208731565052655e-05, -4.6888305906027664e-05, -7.6469715151629098e-06,
         3.9851680786178314e-05,  9.5815298552671653e-06, -4.6471940965075453e-05,
        -7.5790669325201704e-06,  9.4964465272319093e-06,
         9.8837846986583506e-05,  4.0861517241056017e-05,  5.4635998823661205e-06,
         9.7960173683417302e-05, -5.0402301787694723e-06,  4.0498669770145315e-05,
         5.4150834900919591e-06, -4.9954732804295155e-06,
         1.6363308286919829e-05, -1.9376034559952714e-05, -2.6746074969222150e-06,
         1.6218003231490433e-05,  3.2094122010726464e-06, -1.9203976701824692e-05,
        -2.6508571657680423e-06,  3.1809128408213089e-06,
         2.5559272241391462e-05,  6.3438828943389592e-06,  9.7748053065870137e-07,
         2.5332307656684648e-05, -1.7965596045881944e-06,  6.2875496492862865e-06,
         9.6880057057984307e-07, -1.7806062785034210e-06,
      };
  
      int nvals0 = sizeof(f0m)/sizeof(f0m[0]);
      int nvals1 = sizeof(f1m)/sizeof(f1m[0]);
      for (int i=0; i<nvals0; i++)
        f0[i] = f0m[i];
      for (int i=0; i<nvals1; i++)
        f1[i] = f1m[i];
    }
  
    struct gkyl_range check_ghost_rng = ghost_rng;
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &check_ghost_rng);
    while (gkyl_range_iter_next(&iter)) {
      for (int i=0; i<num_ref0; i++) {
        if (iter.idx[0]==idx0_ref[i]) {
          long linidx = gkyl_range_idx(&check_ghost_rng, iter.idx);
          double *f_c = gkyl_array_fetch(distf_ho, linidx);
          int refidx = iter.idx[1]-1;
          for (int k=0; k<basis.num_basis; ++k) {
            TEST_CHECK( ts_compare_coeff(f0[i*cells[1]*basis.num_basis+refidx*basis.num_basis+k], f_c[k]) );
          }
        }
      }
      for (int i=0; i<num_ref1; i++) {
        if (iter.idx[1]==idx1_ref[i]) {
          long linidx = gkyl_range_idx(&check_ghost_rng, iter.idx);
          double *f_c = gkyl_array_fetch(distf_ho, linidx);
          int refidx = iter.idx[0]-1;
          for (int k=0; k<basis.num_basis; ++k) {
            TEST_CHECK( ts_compare_coeff(f1[i*cells[0]*basis.num_basis+refidx*basis.num_basis+k], f_c[k]) );
          }
        }
      }
    }
  }


  gkyl_free(app_ctx.r_lut);
  gkyl_free(app_ctx.dPsidr_int_lut);

  gkyl_array_release(buff_per);
  test_bc_twistshift_array_meta_release(mt);
  gkyl_twistshift_dg_release(tsup);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(distf);
}

void
test_bc_twistshift_3x_cbc(bool use_gpu)
{
  const int cells0[] = {32, 16, 4};

  test_bc_twistshift_3x_cbc_wcells(cells0, GKYL_LOWER_EDGE, true, use_gpu, false);
  test_bc_twistshift_3x_cbc_wcells(cells0, GKYL_UPPER_EDGE, true, use_gpu, false);
}

void
test_bc_twistshift_3x_fig6(bool use_gpu)
{
  const int cells0[] = {1, 10, 4};

  enum gkyl_edge_loc edgelo = GKYL_LOWER_EDGE; // Lower edge.
  test_bc_twistshift_3x_fig6_wcells(cells0, edgelo, true, use_gpu, false);
}

void
test_bc_twistshift_3x2v_fig6(bool use_gpu)
{
  const int cells0[] = {1, 10, 4, 2, 1};

  enum gkyl_edge_loc edgelo = GKYL_LOWER_EDGE; // Lower edge.
  test_bc_twistshift_3x2v_fig6_wcells(cells0, edgelo, true, use_gpu, false);
}

void
test_bc_twistshift_3x_fig11(bool use_gpu)
{
  const int cells0[] = {10, 5,  4};
  const int cells1[] = {20, 10, 4};
  const int cells2[] = {40, 20, 4};
  const int cells3[] = {80, 40, 4};

  enum gkyl_edge_loc edgelo = GKYL_LOWER_EDGE; // Lower edge.
  test_bc_twistshift_3x_fig11_wcells(cells0, edgelo, 0, true, use_gpu, false);
  test_bc_twistshift_3x_fig11_wcells(cells1, edgelo, 0, false, use_gpu, false);
  test_bc_twistshift_3x_fig11_wcells(cells2, edgelo, 0, false, use_gpu, false);
  test_bc_twistshift_3x_fig11_wcells(cells3, edgelo, 0, false, use_gpu, false);

  enum gkyl_edge_loc edgeup = GKYL_UPPER_EDGE; // Upper edge.
  test_bc_twistshift_3x_fig11_wcells(cells0, edgeup, 0, true, use_gpu, false);
  test_bc_twistshift_3x_fig11_wcells(cells1, edgeup, 0, false, use_gpu, false);
  test_bc_twistshift_3x_fig11_wcells(cells2, edgeup, 0, false, use_gpu, false);
  test_bc_twistshift_3x_fig11_wcells(cells3, edgeup, 0, false, use_gpu, false);

  // Apply the TS BC on the lower half of the x domain.
  test_bc_twistshift_3x_fig11_wcells(cells0, edgelo, -1, true, use_gpu, false);
  test_bc_twistshift_3x_fig11_wcells(cells1, edgelo, -1, false, use_gpu, false);
  test_bc_twistshift_3x_fig11_wcells(cells2, edgelo, -1, false, use_gpu, false);

  // Apply the TS BC on the upper half of the x domain.
  test_bc_twistshift_3x_fig11_wcells(cells0, edgelo, 1, true, use_gpu, false);
  test_bc_twistshift_3x_fig11_wcells(cells1, edgelo, 1, false, use_gpu, false);
  test_bc_twistshift_3x_fig11_wcells(cells2, edgelo, 1, false, use_gpu, false);
}

void
test_bc_twistshift_3x_fig14(bool use_gpu)
{
  const int cells0[] = {8 , 8 , 4};
  const int cells1[] = {16, 16, 4};
  const int cells2[] = {32, 32, 4};
  const int cells3[] = {64, 32, 4};

  enum gkyl_edge_loc edgelo = GKYL_LOWER_EDGE; // Lower edge.
  test_bc_twistshift_3x_fig14_wcells(cells0, edgelo, 0, true, use_gpu, false);
  test_bc_twistshift_3x_fig14_wcells(cells1, edgelo, 0, false, use_gpu, false);
  test_bc_twistshift_3x_fig14_wcells(cells2, edgelo, 0, false, use_gpu, false);
  test_bc_twistshift_3x_fig14_wcells(cells3, edgelo, 0, false, use_gpu, false);

  enum gkyl_edge_loc edgeup = GKYL_UPPER_EDGE; // Upper edge.
  test_bc_twistshift_3x_fig14_wcells(cells0, edgeup, 0, true, use_gpu, false);
  test_bc_twistshift_3x_fig14_wcells(cells1, edgeup, 0, false, use_gpu, false);
  test_bc_twistshift_3x_fig14_wcells(cells2, edgeup, 0, false, use_gpu, false);
  test_bc_twistshift_3x_fig14_wcells(cells3, edgeup, 0, false, use_gpu, false);
}

void
test_bc_twistshift_3x2v_fig11(bool use_gpu)
{
  const int cells0[] = {10, 5,  4, 2, 1};
  const int cells1[] = {20, 10, 4, 2, 1};
  const int cells2[] = {40, 20, 4, 2, 1};
  const int cells3[] = {80, 40, 4, 2, 1};

  enum gkyl_edge_loc edgelo = GKYL_LOWER_EDGE; // Lower edge.
  test_bc_twistshift_3x2v_fig11_wcells(cells0, edgelo, 0, true, use_gpu, false);
  test_bc_twistshift_3x2v_fig11_wcells(cells1, edgelo, 0, false, use_gpu, false);
  test_bc_twistshift_3x2v_fig11_wcells(cells2, edgelo, 0, false, use_gpu, false);
  test_bc_twistshift_3x2v_fig11_wcells(cells3, edgelo, 0, false, use_gpu, false);

  enum gkyl_edge_loc edgeup = GKYL_UPPER_EDGE; // Upper edge.
  test_bc_twistshift_3x2v_fig11_wcells(cells0, edgeup, 0, true, use_gpu, false);
  test_bc_twistshift_3x2v_fig11_wcells(cells1, edgeup, 0, false, use_gpu, false);
  test_bc_twistshift_3x2v_fig11_wcells(cells2, edgeup, 0, false, use_gpu, false);
  test_bc_twistshift_3x2v_fig11_wcells(cells3, edgeup, 0, false, use_gpu, false);

  // Apply the TS BC on the lower half of the x domain.
  test_bc_twistshift_3x2v_fig11_wcells(cells0, edgelo, -1, true, use_gpu, false);
  test_bc_twistshift_3x2v_fig11_wcells(cells1, edgelo, -1, false, use_gpu, false);
  test_bc_twistshift_3x2v_fig11_wcells(cells2, edgelo, -1, false, use_gpu, false);

  // Apply the TS BC on the upper half of the x domain.
  test_bc_twistshift_3x2v_fig11_wcells(cells0, edgelo, 1, true, use_gpu, false);
  test_bc_twistshift_3x2v_fig11_wcells(cells1, edgelo, 1, false, use_gpu, false);
  test_bc_twistshift_3x2v_fig11_wcells(cells2, edgelo, 1, false, use_gpu, false);
}

void test_bc_twistshift_3x_fig6_ho(){ test_bc_twistshift_3x_fig6(false); }
void test_bc_twistshift_3x_fig11_ho(){ test_bc_twistshift_3x_fig11(false); }
void test_bc_twistshift_3x_fig14_ho(){ test_bc_twistshift_3x_fig14(false); }
void test_bc_twistshift_3x_cbc_ho(){ test_bc_twistshift_3x_cbc(false); }

void test_bc_twistshift_3x2v_fig6_ho(){ test_bc_twistshift_3x2v_fig6(false); }
void test_bc_twistshift_3x2v_fig11_ho(){ test_bc_twistshift_3x2v_fig11(false); }

#ifdef GKYL_HAVE_CUDA
void test_bc_twistshift_3x_fig6_dev(){ test_bc_twistshift_3x_fig6(true); }
void test_bc_twistshift_3x_fig11_dev(){ test_bc_twistshift_3x_fig11(true); }
void test_bc_twistshift_3x_fig14_dev(){ test_bc_twistshift_3x_fig14(true); }
void test_bc_twistshift_3x_cbc_dev(){ test_bc_twistshift_3x_cbc(true); }

void test_bc_twistshift_3x2v_fig6_dev(){ test_bc_twistshift_3x2v_fig6(true); }
void test_bc_twistshift_3x2v_fig11_dev(){ test_bc_twistshift_3x2v_fig11(true); }
#endif

TEST_LIST = {
  { "test_bc_twistshift_3x_fig6_ho", test_bc_twistshift_3x_fig6_ho },
  { "test_bc_twistshift_3x_fig11_ho", test_bc_twistshift_3x_fig11_ho },
  { "test_bc_twistshift_3x_fig14_ho", test_bc_twistshift_3x_fig14_ho },
  { "test_bc_twistshift_3x_cbc_ho", test_bc_twistshift_3x_cbc_ho },
  { "test_bc_twistshift_3x2v_fig6_ho", test_bc_twistshift_3x2v_fig6_ho },
  { "test_bc_twistshift_3x2v_fig11_ho", test_bc_twistshift_3x2v_fig11_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_bc_twistshift_3x_fig6_dev", test_bc_twistshift_3x_fig6_dev },
  { "test_bc_twistshift_3x_fig11_dev", test_bc_twistshift_3x_fig11_dev },
  { "test_bc_twistshift_3x_fig14_dev", test_bc_twistshift_3x_fig14_dev },
  { "test_bc_twistshift_3x_cbc_dev", test_bc_twistshift_3x_cbc_dev },
  { "test_bc_twistshift_3x2v_fig6_dev", test_bc_twistshift_3x2v_fig6_dev },
  { "test_bc_twistshift_3x2v_fig11_dev", test_bc_twistshift_3x2v_fig11_dev },
#endif
  { NULL, NULL },
};
