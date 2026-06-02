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
#include <gkyl_bc_twistshift.h>
#include <gkyl_velocity_map.h>
#include <gkyl_position_map.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>

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

  const int poly_order = 1;
  const double lower[] = {-2.0, -1.50, -3.0};
  const double upper[] = { 2.0,  1.50,  3.0};
  const int vdim = 0;
  const int ndim = sizeof(lower)/sizeof(lower[0]);
  const int cdim = ndim - vdim;

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
  if (write_f){}

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
  struct gkyl_bc_twistshift_inp tsinp = {
    .bc_dir = bc_dir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = edge,
    .cdim = cdim,
    .bcdir_ext_update_r = update_rng,
    .num_ghost = ghost,
    .basis = basis,
    .grid = grid,
    .shift_func = shift1_fig6,
//    .shift_func = shift_fig9,
    .shift_func_ctx = &proj_ctx,
    .use_gpu = use_gpu,
  };

  struct gkyl_bc_twistshift *tsup = gkyl_bc_twistshift_new(&tsinp);

  // First apply periodicity in z.
  struct gkyl_array *buff_per = mkarr(use_gpu, basis.num_basis, skin_rng.volume);
  apply_periodic_bc(buff_per, distf, bc_dir, skin_ghost);

  gkyl_bc_twistshift_advance(tsup, distf, distf);
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
    // gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "ctest_bc_twistshift_3x_fig6_tar.gkyl");
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
  struct gkyl_bc_twistshift *tsup_m = gkyl_bc_twistshift_new(&tsinp);
  gkyl_bc_twistshift_advance(tsup_m, distf, distf);
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
    // gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "ctest_bc_twistshift_3x_fig6_tar_shifted.gkyl");
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
  gkyl_bc_twistshift_release(tsup);
  gkyl_bc_twistshift_release(tsup_m);
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

  const int poly_order = 1;
  const double lower[] = {-2.0, -1.50, -3.0, -5.0*vt, 0.};
  const double upper[] = { 2.0,  1.50,  3.0,  5.0*vt, mass*(pow(5.0*vt,2))/(2.0*B0)};
  const int vdim = 2;
  const int ndim = sizeof(lower)/sizeof(lower[0]);
  const int cdim = ndim - vdim;

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
  if (write_f) {}
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
  struct gkyl_bc_twistshift_inp tsinp = {
    .bc_dir = bc_dir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = edge,
    .cdim = cdim,
    .bcdir_ext_update_r = update_rng,
    .num_ghost = ghost,
    .basis = basis,
    .grid = grid,
    .shift_func = shift1_fig6,
//    .shift_func = shift_fig9,
    .shift_func_ctx = &proj_ctx,
    .use_gpu = use_gpu,
  };

  struct gkyl_bc_twistshift *tsup = gkyl_bc_twistshift_new(&tsinp);

  // First apply periodicity in z.
  struct gkyl_array *buff_per = mkarr(use_gpu, basis.num_basis, skin_rng.volume);
  apply_periodic_bc(buff_per, distf, bc_dir, skin_ghost);

  gkyl_bc_twistshift_advance(tsup, distf, distf);
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
    // gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "ctest_bc_twistshift_3x2v_fig6_tar.gkyl");
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
  struct gkyl_bc_twistshift *tsup_m = gkyl_bc_twistshift_new(&tsinp);
  gkyl_bc_twistshift_advance(tsup_m, distf, distf);
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
    // gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "ctest_bc_twistshift_3x2v_fig6_tar_shifted.gkyl");
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
  gkyl_bc_twistshift_release(tsup);
  gkyl_bc_twistshift_release(tsup_m);
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

  const int poly_order = 1;
  const double lower[] = {-2.0, -1.50, -3.0};
  const double upper[] = { 2.0,  1.50,  3.0};
  const int vdim = 0;
  const int ndim = sizeof(lower)/sizeof(lower[0]);
  const int cdim = ndim - vdim;

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
  if (write_f){}
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
  struct gkyl_bc_twistshift_inp tsinp = {
    .bc_dir = bc_dir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = edge,
    .cdim = cdim,
    .bcdir_ext_update_r = update_rng,
    .num_ghost = ghost,
    .basis = basis,
    .grid = grid,
    .shift_func = shift_fig11,
    .shift_func_ctx = &proj_ctx,
    .use_gpu = use_gpu,
  };

  struct gkyl_bc_twistshift *tsup = gkyl_bc_twistshift_new(&tsinp);

  // First apply periodicity in z.
  struct gkyl_array *buff_per = mkarr(use_gpu, basis.num_basis, skin_rng.volume);
  apply_periodic_bc(buff_per, distf, bc_dir, skin_ghost);

  gkyl_bc_twistshift_advance(tsup, distf, distf);
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
    // gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "ctest_bc_twistshift_3x_fig11_tar.gkyl");
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
  gkyl_bc_twistshift_release(tsup);
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

  const int poly_order = 1;
  const double lower[] = {-2.0, -1.50, -3.0, -5.0*vt, 0.};
  const double upper[] = { 2.0,  1.50,  3.0,  5.0*vt, mass*(pow(5.0*vt,2))/(2.0*B0)};
  const int vdim = 2;
  const int ndim = sizeof(lower)/sizeof(lower[0]);
  const int cdim = ndim - vdim;

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
  if (write_f) {}

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
  struct gkyl_bc_twistshift_inp tsinp = {
    .bc_dir = bc_dir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = edge,
    .cdim = cdim,
    .bcdir_ext_update_r = update_rng,
    .num_ghost = ghost,
    .basis = basis,
    .grid = grid,
    .shift_func = shift_fig11,
    .shift_func_ctx = &proj_ctx,
    .use_gpu = use_gpu,
  };

  struct gkyl_bc_twistshift *tsup = gkyl_bc_twistshift_new(&tsinp);

  // First apply periodicity in z.
  struct gkyl_array *buff_per = mkarr(use_gpu, basis.num_basis, skin_rng.volume);
  apply_periodic_bc(buff_per, distf, bc_dir, skin_ghost);

  gkyl_bc_twistshift_advance(tsup, distf, distf);
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
    // gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "ctest_bc_twistshift_3x2v_fig11_tar.gkyl");
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
  gkyl_bc_twistshift_release(tsup);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(distf);

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
  test_bc_twistshift_3x_fig11_wcells(cells0, edgelo, 0, true, use_gpu, true);
  test_bc_twistshift_3x_fig11_wcells(cells1, edgelo, 0, false, use_gpu, false);
  test_bc_twistshift_3x_fig11_wcells(cells2, edgelo, 0, false, use_gpu, false);
  test_bc_twistshift_3x_fig11_wcells(cells3, edgelo, 0, false, use_gpu, true);

  enum gkyl_edge_loc edgeup = GKYL_UPPER_EDGE; // Upper edge.
  test_bc_twistshift_3x_fig11_wcells(cells0, edgeup, 0, true, use_gpu, false);
  test_bc_twistshift_3x_fig11_wcells(cells1, edgeup, 0, false, use_gpu, false);
  test_bc_twistshift_3x_fig11_wcells(cells2, edgeup, 0, false, use_gpu, false);
  test_bc_twistshift_3x_fig11_wcells(cells3, edgeup, 0, false, use_gpu, true);

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

void test_bc_twistshift_3x2v_fig6_ho(){ test_bc_twistshift_3x2v_fig6(false); }
void test_bc_twistshift_3x2v_fig11_ho(){ test_bc_twistshift_3x2v_fig11(false); }

#ifdef GKYL_HAVE_CUDA
void test_bc_twistshift_3x_fig6_dev(){ test_bc_twistshift_3x_fig6(true); }
void test_bc_twistshift_3x_fig11_dev(){ test_bc_twistshift_3x_fig11(true); }

void test_bc_twistshift_3x2v_fig6_dev(){ test_bc_twistshift_3x2v_fig6(true); }
void test_bc_twistshift_3x2v_fig11_dev(){ test_bc_twistshift_3x2v_fig11(true); }
#endif

TEST_LIST = {
  { "test_bc_twistshift_3x_fig6_ho", test_bc_twistshift_3x_fig6_ho },
  { "test_bc_twistshift_3x_fig11_ho", test_bc_twistshift_3x_fig11_ho },
  { "test_bc_twistshift_3x2v_fig6_ho", test_bc_twistshift_3x2v_fig6_ho },
  { "test_bc_twistshift_3x2v_fig11_ho", test_bc_twistshift_3x2v_fig11_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_bc_twistshift_3x_fig6_dev", test_bc_twistshift_3x_fig6_dev },
  { "test_bc_twistshift_3x_fig11_dev", test_bc_twistshift_3x_fig11_dev },
  { "test_bc_twistshift_3x2v_fig6_dev", test_bc_twistshift_3x2v_fig6_dev },
  { "test_bc_twistshift_3x2v_fig11_dev", test_bc_twistshift_3x2v_fig11_dev },
#endif
  { NULL, NULL },
};
