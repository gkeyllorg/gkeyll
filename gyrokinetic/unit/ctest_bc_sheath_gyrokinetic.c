// Test creation and deallocation of updater that applies the
// conducting sheath BC for gyrokinetics.
//
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_velocity_map.h>
#include <gkyl_proj_on_basis.h>
#include <float.h>

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (use_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct test_sheath_ctx {
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM]; // Grid extents.
  int cells[GKYL_MAX_DIM]; // Number of cells.
  double upar; // Mean parallel drift speed.
  double B0; // Reference magnetic field.
  double vt; // Reference thermal speed.
  double mass; // Species mass.
  double x0; // x-center
  double y0; // y-center
  double z0; // z-center
  double sigmax; // Width of distribution in x.
  double sigmay; // Width of distribution in y.
  double sigmaz; // Width of distribution in z.
};

void
eval_func_1x2v(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double vpar = xn[1], mu = xn[2];
  double z = xn[0];

  struct test_sheath_ctx *pars = ctx;

  double vt = pars->vt;
  double B0 = pars->B0;
  double m = pars->mass;
  double upar = pars->upar;
  double z0 = pars->z0;
  double sigmaz = pars->sigmaz;

  double envelope = exp( -pow(z-z0,2)/(2.0*pow(sigmaz,2)) );
  fout[0] = exp( -(pow(vpar-upar,2) + 2.0*mu*B0/m) /(2.0*pow(vt,2)) ) * envelope + 0.1;
}

void
eval_func_2x2v(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double vpar = xn[2], mu = xn[3];
  double x = xn[0], z = xn[1];

  struct test_sheath_ctx *pars = ctx;

  double vt = pars->vt;
  double B0 = pars->B0;
  double m = pars->mass;
  double upar = pars->upar;
  double x0 = pars->x0;
  double z0 = pars->z0;
  double sigmax = pars->sigmax;
  double sigmaz = pars->sigmaz;

  double envelope = exp( -pow(x-x0,2)/(2.0*pow(sigmax,2)) ) * exp( -pow(z-z0,2)/(2.0*pow(sigmaz,2)) );
  fout[0] = exp( -(pow(vpar-upar,2) + 2.0*mu*B0/m) /(2.0*pow(vt,2)) ) * envelope + 0.1;
}

void
eval_func_3x2v(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double vpar = xn[3], mu = xn[4];
  double x = xn[0], y = xn[1], z = xn[2];

  struct test_sheath_ctx *pars = ctx;

  double vt = pars->vt;
  double B0 = pars->B0;
  double m = pars->mass;
  double upar = pars->upar;
  double x0 = pars->x0;
  double y0 = pars->y0;
  double z0 = pars->z0;
  double sigmax = pars->sigmax;
  double sigmay = pars->sigmay;
  double sigmaz = pars->sigmaz;

  double envelope = exp( -pow(x-x0,2)/(2.0*pow(sigmax,2)) ) * exp( -pow(y-y0,2)/(2.0*pow(sigmay,2)) ) * exp( -pow(z-z0,2)/(2.0*pow(sigmaz,2)) );
  fout[0] = exp( -(pow(vpar-upar,2) + 2.0*mu*B0/m) /(2.0*pow(vt,2)) ) * envelope + 0.1;
}

// Checks that the distribution function values in the ghost cells are set to 
// 0 for velocities beyond the cutoff velocity and >0 otherwise.
void
check_function(double phi_mpe, double phi_wall, double charge, double mass, int cdim,
  struct gkyl_array *distf_ho, struct gkyl_rect_grid grid, struct gkyl_range ghost_r, enum gkyl_edge_loc edge){

  double delta_phi = phi_mpe - phi_wall;
  double vpar_cut;
  if (delta_phi < 0) {
    vpar_cut = charge < 0.? DBL_MAX : sqrt(charge*fabs(delta_phi)/(0.5*mass));
  }
  else {
    vpar_cut = charge < 0.? sqrt(fabs(charge)*delta_phi/(0.5*mass)) : DBL_MAX;
  }
  int num_cells = 0;
  int num_zero_cells = 0;
  int num_zero_cells_expected = 0;
  // We need the sign of charge*delta_phi to determine 
  // which part of the distribution function is cut off by the sheath BC.
  double qphi_sign = charge * delta_phi > 0 ? 1.0 : -1.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ghost_r);
  while (gkyl_range_iter_next(&iter)) {
    int *idx_g = iter.idx;
    long linidx_g = gkyl_range_idx(&ghost_r, idx_g);
    double *distf_c = gkyl_array_fetch(distf_ho, linidx_g);
    double tol = 1e-8;
    double ref_val = 0.0;

    // Get cell boundaries along vpar.
    double xc[GKYL_MAX_DIM] = {0};
    gkyl_rect_grid_cell_center(&grid, idx_g, xc);
    double cell_lower_vpar = xc[cdim] - 0.5*grid.dx[cdim];
    double cell_upper_vpar = xc[cdim] + 0.5*grid.dx[cdim];
    
    num_cells++;
    num_zero_cells += distf_c[0] > tol ? 0 : 1;
    if ( edge == GKYL_LOWER_EDGE && 
        ( cell_lower_vpar > -qphi_sign*vpar_cut || cell_upper_vpar < qphi_sign*vpar_cut ) ) {
      num_zero_cells_expected++;
      for (int k=0; k<distf_ho->ncomp; k++) {
        TEST_CHECK( gkyl_compare(distf_c[k], ref_val, tol) );
        TEST_MSG( "Expected %.9e | Got: %.9e at idx=%d,%d,%d\n", ref_val, distf_c[k], idx_g[0], idx_g[1], idx_g[2]);
      }
    }
    else if (edge == GKYL_UPPER_EDGE && 
        ( cell_lower_vpar > -qphi_sign*vpar_cut || cell_upper_vpar < qphi_sign*vpar_cut ) ) {
      num_zero_cells_expected++;
      for (int k=0; k<distf_ho->ncomp; k++) {
        TEST_CHECK( gkyl_compare(distf_c[k], ref_val, tol) );
        TEST_MSG( "Expected %.9e | Got: %.9e at idx=%d,%d,%d\n", ref_val, distf_c[k], idx_g[0], idx_g[1], idx_g[2]);
      }
    } else {
      // Check that the average cell value is bigger than 0 for cells that are not expected to be cut off by the sheath BC.
      double cell_avg = distf_c[0];
      TEST_CHECK( cell_avg > ref_val );
      TEST_MSG( "Expected > %.9e | Got: %.9e at idx=%d,%d,%d\n", ref_val, cell_avg, idx_g[0], idx_g[1], idx_g[2]);
    }
  }
  // Uncomment for debugging.
  // printf("\n");
  // printf("Edge: %s | Charge: %.1e | Delta phi: %.2e\n", 
  //   edge == GKYL_LOWER_EDGE? "lower" : "upper", charge, delta_phi);
  // printf("zero cell = %d | zero cell expected = %d | total cells = %d\n", num_zero_cells, num_zero_cells_expected, num_cells);
}

void write_out_fields(int cdim, int vdim, enum gkyl_edge_loc edge, bool use_gpu,
  struct gkyl_array *distf_ho, struct gkyl_array *phi_ho, struct gkyl_array *phiw_ho,
  struct gkyl_rect_grid grid_ext, struct gkyl_range local_ext,
  struct gkyl_rect_grid grid_conf, struct gkyl_range local_conf, struct gkyl_basis *basis, struct gkyl_basis basis_conf) {

  struct gkyl_msgpack_map_elem io_meta[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = basis->poly_order },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = basis->id }
  };
  int io_meta_len = sizeof(io_meta)/sizeof(io_meta[0]);
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create(io_meta_len, io_meta);

  struct gkyl_msgpack_map_elem io_meta_conf[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = basis_conf.poly_order },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = basis_conf.id }
  };
  int io_meta_conf_len = sizeof(io_meta_conf)/sizeof(io_meta_conf[0]);
  struct gkyl_msgpack_data *mt_conf = gkyl_msgpack_create(io_meta_conf_len, io_meta_conf);

  char fname[256];
  const char *fmt = "bc_sheath_%dx%dv_%s_%s_%s.gkyl";
  snprintf(fname, sizeof(fname), fmt, cdim, vdim, edge == GKYL_LOWER_EDGE? "lower" : "upper", use_gpu? "gpu" : "cpu", "distf_out");
  gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, fname);
  snprintf(fname, sizeof(fname), fmt, cdim, vdim, edge == GKYL_LOWER_EDGE? "lower" : "upper", use_gpu? "gpu" : "cpu", "phi_mpe");
  gkyl_grid_sub_array_write(&grid_conf, &local_conf, mt_conf, phi_ho, fname);
  snprintf(fname, sizeof(fname), fmt, cdim, vdim, edge == GKYL_LOWER_EDGE? "lower" : "upper", use_gpu? "gpu" : "cpu", "phi_wall");
  gkyl_grid_sub_array_write(&grid_conf, &local_conf, mt_conf, phiw_ho, fname);

  gkyl_msgpack_data_release(mt);
  gkyl_msgpack_data_release(mt_conf);
}

void
test_bc_sheath_gyrokinetic_1x2v(const int *cells, enum gkyl_edge_loc edge,
  double charge, double phi_mpe, bool write_fields, bool use_gpu)
{
  /*
  This test applies the sheath BC to a Maxwellian distribution function on both upper and lower edges in the z direction,
  according to a potential drop between the magnetic presheath entrance and the wall. 
  The folloing files are written out if requested:
  - "bc_sheath_1x2v_distf_in_${edge}_${cpugpu}.gkyl": Input distribution function.
  - "bc_sheath_1x2v_phi_mpe.gkyl": Potential at the magnetic presheath entrance (constant everywhere).
  - "bc_sheath_1x2v_phi_wall.gkyl": Potential at the wall (constant everywhere).
  - "bc_sheath_1x2v_distf_out_${edge}_${cpugpu}.gkyl": Distribution after the sheath BC has been applied (includes ghost cells).
  */

  double mass = 1.0; // Species mass.
  double vt = 1.0; // Reference thermal speed (for grid extents).
  double B0 = 1.0; // Magnetic field magnitude.
  double upar_distf = 1.0*vt; // Parallel flow speed in distribution function.
  double vt_distf = 1.5*vt; // Thermal speed in distribution function.
  double phi_wall = 0.0; // Potential at wall.
  double z0 = 0.2;
  double sigmaz = 1.0;
  
  int poly_order = 1;
  double lower[] = {-1.0, -5.0*vt, 0.};
  double upper[] = { 1.0,  5.0*vt, mass*(pow(5.0*vt,2))/(2.0*B0)};
  int vdim = 2;
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int cdim = ndim - vdim;
  int dir = cdim-1; // Always applied in the parallel direction.

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
  struct gkyl_basis *basis;
  basis = use_gpu? gkyl_cart_modal_gkhybrid_cu_dev_new(cdim, vdim)
                  : gkyl_cart_modal_gkhybrid_new(cdim, vdim);

  struct gkyl_basis basis_ho;
  gkyl_cart_modal_gkhybrid(&basis_ho, cdim, vdim);

  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  int ghost_conf[] = { 1, 1, 1 };
  struct gkyl_range local_conf, local_conf_ext; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_conf_ext, &local_conf);

  int ghost_vel[] = { 0, 0 };
  struct gkyl_range local_vel, local_vel_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_vel_ext, &local_vel);

  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; d++)
    ghost[d] = ghost_conf[d];
  for (int d=0; d<vdim; d++)
    ghost[cdim+d] = ghost_vel[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_vel_ext, use_gpu);

  // Extended grid for the distribution function, which includes ghost cells.
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

  // Create the skin/ghost ranges.
  struct gkyl_range skin_r, ghost_r;
  gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

  // Initialize the distribution
  struct gkyl_array *distf = mkarr(use_gpu, basis_ho.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size) : gkyl_array_acquire(distf);
  struct test_sheath_ctx proj_ctx = {
    .B0 = B0,
    .mass = mass,
    .upar = upar_distf,
    .vt = vt_distf,
    .z0 = z0,
    .sigmaz = sigmaz
  };
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis_ho,
      .num_ret_vals = 1,
      .eval = eval_func_1x2v,
      .ctx = &proj_ctx,
    }
  );
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // Initialize the electrostatic potential at MPE and WALL
  double dgnormc = pow(sqrt(2.0), cdim);
  struct gkyl_array *phi = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phi_ho = use_gpu? mkarr(false, phi->ncomp, phi->size) : gkyl_array_acquire(phi);
  gkyl_array_shiftc(phi_ho, phi_mpe * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phi, phi_ho);

  struct gkyl_array *phiw = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phiw_ho = use_gpu? mkarr(false, phiw->ncomp, phiw->size) : gkyl_array_acquire(phiw);
  gkyl_array_shiftc(phiw_ho, phi_wall * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phiw, phiw_ho);

  // Create the BC updater.
  struct gkyl_bc_sheath_gyrokinetic *bcsheath = gkyl_bc_sheath_gyrokinetic_new(dir, edge,
    basis, &skin_r, &ghost_r, gvm, cdim, 2.*charge/mass, use_gpu);

  // Advance the BC updater.
  gkyl_bc_sheath_gyrokinetic_advance(bcsheath, phi, phiw, distf, &local_conf);

  // Copy back to host.
  gkyl_array_copy(distf_ho, distf);

  // Check the result.
  check_function(phi_mpe, phi_wall, charge, mass, cdim, distf_ho, grid, ghost_r, edge);

  // Write out the distribution function after applying BC if requested.
  if (write_fields)
    write_out_fields(cdim, vdim, edge, use_gpu, distf_ho, phi_ho, phiw_ho, 
      grid_ext, local_ext, grid_conf, local_conf, &basis_ho, basis_conf);

  // Clean up.
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(phiw);
  gkyl_array_release(phiw_ho);
  if (use_gpu)
    gkyl_cart_modal_basis_release_cu(basis);
  else
    gkyl_cart_modal_basis_release(basis);

  gkyl_velocity_map_release(gvm);
  gkyl_bc_sheath_gyrokinetic_release(bcsheath);
}

void
test_bc_sheath_gyrokinetic_2x2v(const int *cells, enum gkyl_edge_loc edge,
  double charge, double phi_mpe, bool write_fields, bool use_gpu)
{
  /*
  This test applies the sheath BC to a Maxwellian distribution function on both upper and lower edges in the z direction,
  according to a potential drop between the magnetic presheath entrance and the wall.
  The following files are written out if requested:
  - "bc_sheath_2x2v_distf_in.gkyl": Input distribution function.
  - "bc_sheath_2x2v_phi_mpe.gkyl": Potential at the magnetic presheath entrance (constant everywhere).
  - "bc_sheath_2x2v_phi_wall.gkyl": Potential at the wall (constant everywhere).
  - "bc_sheath_2x2v_distf_out.gkyl": Distribution after the sheath BC has been applied (includes ghost cells).
  */

  double mass = 1.; // Species mass.
  double vt = 1.0; // Reference thermal speed (for grid extents).
  double B0 = 1.0; // Magnetic field magnitude.
  double upar_distf = 1.0*vt; // Parallel flow speed in distribution function.
  double vt_distf = 1.5*vt; // Thermal speed in distribution function.
  double phi_wall = 0.0; // Potential at wall.
  double x0 = 0.0; // x-center of distribution function
  double z0 = 0.2; // z-center of distribution function
  double sigmax = 0.5; // Width of distribution function in x.
  double sigmaz = 1.0; // Width of distribution function in z.

  int poly_order = 1;
  double lower[] = { 0.0, -1.0, -5.0*vt, 0.};
  double upper[] = { 1.0,  1.0,  5.0*vt, mass*(pow(5.0*vt,2))/(2.0*B0)};
  int vdim = 2;
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int cdim = ndim - vdim;
  int dir = cdim-1; // Always applied in the parallel direction.

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
  struct gkyl_basis *basis;
  basis = use_gpu? gkyl_cart_modal_gkhybrid_cu_dev_new(cdim, vdim)
                  : gkyl_cart_modal_gkhybrid_new(cdim, vdim);

  struct gkyl_basis basis_ho;
  gkyl_cart_modal_gkhybrid(&basis_ho, cdim, vdim);

  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  int ghost_conf[] = { 1, 1, 1 };
  struct gkyl_range local_conf, local_conf_ext; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_conf_ext, &local_conf);

  int ghost_vel[] = { 0, 0 };
  struct gkyl_range local_vel, local_vel_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_vel_ext, &local_vel);

  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; d++)
    ghost[d] = ghost_conf[d];
  for (int d=0; d<vdim; d++)
    ghost[cdim+d] = ghost_vel[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_vel_ext, use_gpu);

  // Extended grid for the distribution function, which includes ghost cells.
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

  // Create the skin/ghost ranges.
  struct gkyl_range skin_r, ghost_r;
  gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

  // Initialize the distribution.
  struct gkyl_array *distf = mkarr(use_gpu, basis_ho.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size) : gkyl_array_acquire(distf);
  struct test_sheath_ctx proj_ctx = {
    .B0 = B0,
    .mass = mass,
    .upar = upar_distf,
    .vt = vt_distf,
    .x0 = x0,
    .z0 = z0,
    .sigmax = sigmax,
    .sigmaz = sigmaz
  };
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis_ho,
      .num_ret_vals = 1,
      .eval = eval_func_2x2v,
      .ctx = &proj_ctx,
    }
  );
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // Initialize the electrostatic potential at MPE and WALL.
  double dgnormc = pow(sqrt(2.0), cdim);
  struct gkyl_array *phi = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phi_ho = use_gpu? mkarr(false, phi->ncomp, phi->size) : gkyl_array_acquire(phi);
  gkyl_array_shiftc(phi_ho, phi_mpe * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phi, phi_ho);

  struct gkyl_array *phiw = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phiw_ho = use_gpu? mkarr(false, phiw->ncomp, phiw->size) : gkyl_array_acquire(phiw);
  gkyl_array_shiftc(phiw_ho, phi_wall * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phiw, phiw_ho);

  // Create the BC updater.
  struct gkyl_bc_sheath_gyrokinetic *bcsheath = gkyl_bc_sheath_gyrokinetic_new(dir, edge,
    basis, &skin_r, &ghost_r, gvm, cdim, 2.*charge/mass, use_gpu);

  // Advance the BC updater.
  gkyl_bc_sheath_gyrokinetic_advance(bcsheath, phi, phiw, distf, &local_conf);

  // Copy back to host.
  gkyl_array_copy(distf_ho, distf);

  // Check the result.
  check_function(phi_mpe, phi_wall, charge, mass, cdim, distf_ho, grid, ghost_r, edge);

  // Write out the distribution function after applying BC if requested.
  if (write_fields)
    write_out_fields(cdim, vdim, edge, use_gpu, distf_ho, phi_ho, phiw_ho, 
      grid_ext, local_ext, grid_conf, local_conf, &basis_ho, basis_conf);

  // Clean up.
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(phiw);
  gkyl_array_release(phiw_ho);
  if (use_gpu)
    gkyl_cart_modal_basis_release_cu(basis);
  else
    gkyl_cart_modal_basis_release(basis);

  gkyl_velocity_map_release(gvm);
  gkyl_bc_sheath_gyrokinetic_release(bcsheath);
}

void
test_bc_sheath_gyrokinetic_3x2v(const int *cells, enum gkyl_edge_loc edge,
  double charge, double phi_mpe, bool write_fields, bool use_gpu)
{
  /*
  This test applies the sheath BC to a Maxwellian distribution function on both upper and lower edges in the z direction,
  according to a potential drop between the magnetic presheath entrance and the wall.
  The following files are written out if requested:
  - "bc_sheath_3x2v_distf_in.gkyl": Input distribution function.
  - "bc_sheath_3x2v_phi_mpe.gkyl": Potential at the magnetic presheath entrance (constant everywhere).
  - "bc_sheath_3x2v_phi_wall.gkyl": Potential at the wall (constant everywhere).
  - "bc_sheath_3x2v_distf_out.gkyl": Distribution after the sheath BC has been applied (includes ghost cells).
  */

  double mass = 1.; // Species mass.
  double vt = 1.0; // Reference thermal speed (for grid extents).
  double B0 = 1.0; // Magnetic field magnitude.
  double upar_distf = 1.0*vt; // Parallel flow speed in distribution function.
  double vt_distf = 1.5*vt; // Thermal speed in distribution function.
  double phi_wall = 0.0; // Potential at wall.

  int poly_order = 1;
  double lower[] = { 0.0, -2.0, -1.0, -5.0*vt, 0.};
  double upper[] = { 1.0,  2.0,  1.0,  5.0*vt, mass*(pow(5.0*vt,2))/(2.0*B0)};
  int vdim = 2;
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int cdim = ndim - vdim;
  int dir = cdim-1; // Always applied in the parallel direction.
  double x0 = 0.0; // x-center of distribution function
  double y0 = 0.0; // y-center of distribution function
  double z0 = 0.2; // z-center of distribution function
  double sigmax = 0.5; // Width of distribution function in x.
  double sigmay = 0.5; // Width of distribution function in y.
  double sigmaz = 1.0; // Width of distribution function in z.

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
  struct gkyl_basis *basis;
  basis = use_gpu? gkyl_cart_modal_gkhybrid_cu_dev_new(cdim, vdim)
                  : gkyl_cart_modal_gkhybrid_new(cdim, vdim);

  struct gkyl_basis basis_ho;
  gkyl_cart_modal_gkhybrid(&basis_ho, cdim, vdim);

  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  int ghost_conf[] = { 1, 1, 1 };
  struct gkyl_range local_conf, local_conf_ext; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_conf_ext, &local_conf);

  int ghost_vel[] = { 0, 0 };
  struct gkyl_range local_vel, local_vel_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_vel_ext, &local_vel);

  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; d++)
    ghost[d] = ghost_conf[d];
  for (int d=0; d<vdim; d++)
    ghost[cdim+d] = ghost_vel[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_vel_ext, use_gpu);

  // Extended grid for the distribution function, which includes ghost cells.
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

  // Create the skin/ghost ranges.
  struct gkyl_range skin_r, ghost_r;
  gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

  // Initialize the distribution.
  struct gkyl_array *distf = mkarr(use_gpu, basis_ho.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size) : gkyl_array_acquire(distf);
  struct test_sheath_ctx proj_ctx = {
    .B0 = B0,
    .mass = mass,
    .upar = upar_distf,
    .vt = vt_distf,
    .x0 = x0,
    .y0 = y0,
    .z0 = z0,
    .sigmax = sigmax,
    .sigmay = sigmay,
    .sigmaz = sigmaz
  };
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis_ho,
      .num_ret_vals = 1,
      .eval = eval_func_3x2v,
      .ctx = &proj_ctx,
    }
  );
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // Initialize the electrostatic potential at MPE and WALL.
  double dgnormc = pow(sqrt(2.0), cdim);
  struct gkyl_array *phi = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phi_ho = use_gpu? mkarr(false, phi->ncomp, phi->size) : gkyl_array_acquire(phi);
  gkyl_array_shiftc(phi_ho, phi_mpe * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phi, phi_ho);

  struct gkyl_array *phiw = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phiw_ho = use_gpu? mkarr(false, phiw->ncomp, phiw->size) : gkyl_array_acquire(phiw);
  gkyl_array_shiftc(phiw_ho, phi_wall * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phiw, phiw_ho);

  // Create the BC updater.
  struct gkyl_bc_sheath_gyrokinetic *bcsheath = gkyl_bc_sheath_gyrokinetic_new(dir, edge,
    basis, &skin_r, &ghost_r, gvm, cdim, 2.*charge/mass, use_gpu);

  // Advance the BC updater.
  gkyl_bc_sheath_gyrokinetic_advance(bcsheath, phi, phiw, distf, &local_conf);

  // Copy back to host.
  gkyl_array_copy(distf_ho, distf);

  // Check the result.
  check_function(phi_mpe, phi_wall, charge, mass, cdim, distf_ho, grid, ghost_r, edge);

  // Write out the distribution function after applying BC if requested.
  if (write_fields)
    write_out_fields(cdim, vdim, edge, use_gpu, distf_ho, phi_ho, phiw_ho, 
      grid_ext, local_ext, grid_conf, local_conf, &basis_ho, basis_conf);

  // Clean up.
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(phiw);
  gkyl_array_release(phiw_ho);
  if (use_gpu)
    gkyl_cart_modal_basis_release_cu(basis);
  else
    gkyl_cart_modal_basis_release(basis);

  gkyl_velocity_map_release(gvm);
  gkyl_bc_sheath_gyrokinetic_release(bcsheath);
}

void test_bc_sheath_gk_1x2v_ho()
{
  double phi_mpe; // Potential at the magnetic presheath entrance.
  double charge; // Species charge (+ or - for electrons/ions).
  bool write_fields;

  // Electrons with positive sheath entrance potential.
  phi_mpe = 1.0;
  charge = -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, false);
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, false);

  // Electrons with negative sheath entrance potential.
  phi_mpe = -1.0;
  charge = -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, false);
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, false);

  // Ions with positive sheath entrance potential.
  phi_mpe = 1.0;
  charge = 1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, false);
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, false);

  // Ions with negative sheath entrance potential.
  phi_mpe = -1.0;
  charge = 1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, false);
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, false);
}

void test_bc_sheath_gk_2x2v_ho()
{
  double phi_mpe; // Potential at the magnetic presheath entrance.
  double charge; // Species charge (+ or - for electrons/ions).
  bool write_fields;
  
  // Electrons with positive sheath entrance potential.
  phi_mpe = 1.0;
  charge = -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, false);
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, false);

  // Electrons with negative sheath entrance potential.
  phi_mpe = -1.0;
  charge = -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, false);
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, false);

  // Ions with positive sheath entrance potential.
  phi_mpe = 1.0;
  charge = 1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, false);
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, false);

  // Ions with negative sheath entrance potential.
  phi_mpe = -1.0;
  charge = 1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, false);
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, false);
}

void test_bc_sheath_gk_3x2v_ho()
{
  double phi_mpe; // Potential at the magnetic presheath entrance.
  double charge; // Species charge (+ or - for electrons/ions).
  bool write_fields;

  // Electrons with positive sheath entrance potential.
  phi_mpe = 1.0;
  charge = -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, false);
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, false);

  // Electrons with negative sheath entrance potential.
  phi_mpe = -1.0;
  charge = -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, false);
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, false);

  // Ions with positive sheath entrance potential.
  phi_mpe = 1.0;
  charge = 1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, false);
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, false);

  // Ions with negative sheath entrance potential.
  phi_mpe = -1.0;
  charge = 1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, false);
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, false);
}

#ifdef GKYL_HAVE_CUDA
void test_bc_sheath_gk_1x2v_dev()
{
  double phi_mpe;
  double charge;
  bool write_fields;

  // Electrons with positive sheath entrance potential.
  phi_mpe = 1.0;
  charge = -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, true);
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, true);

  // Electrons with negative sheath entrance potential.
  phi_mpe = -1.0;
  charge = -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, true);
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, true);

  // Ions with positive sheath entrance potential.
  phi_mpe = 1.0;
  charge = 1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, true);
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, true);

  // Ions with negative sheath entrance potential.
  phi_mpe = -1.0;
  charge = 1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, true);
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, true);
}

void test_bc_sheath_gk_2x2v_dev()
{
  double phi_mpe;
  double charge;
  bool write_fields;

  // Electrons with positive sheath entrance potential.
  phi_mpe = 1.0;
  charge = -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, true);
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, true);

  // Electrons with negative sheath entrance potential.
  phi_mpe = -1.0;
  charge = -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, true);
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, true);

  // Ions with positive sheath entrance potential.
  phi_mpe = 1.0;
  charge = 1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, true);
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, true);

  // Ions with negative sheath entrance potential.
  phi_mpe = -1.0;
  charge = 1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, true);
  test_bc_sheath_gyrokinetic_2x2v((int[]){4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, true);
}

void test_bc_sheath_gk_3x2v_dev()
{
  double phi_mpe;
  double charge;
  bool write_fields;

  // Electrons with positive sheath entrance potential.
  phi_mpe = 1.0;
  charge = -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, true);
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, true);

  // Electrons with negative sheath entrance potential.
  phi_mpe = -1.0;
  charge = -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, true);
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, true);

  // Ions with positive sheath entrance potential.
  phi_mpe = 1.0;
  charge = 1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, true);
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, true);

  // Ions with negative sheath entrance potential.
  phi_mpe = -1.0;
  charge = 1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_LOWER_EDGE, charge, phi_mpe, write_fields, true);
  test_bc_sheath_gyrokinetic_3x2v((int[]){4, 4, 4, 16, 12}, GKYL_UPPER_EDGE, charge, phi_mpe, write_fields, true);
}
#endif

TEST_LIST = {
  { "test_bc_sheath_gk_1x2v_ho", test_bc_sheath_gk_1x2v_ho },
  { "test_bc_sheath_gk_2x2v_ho", test_bc_sheath_gk_2x2v_ho },
  { "test_bc_sheath_gk_3x2v_ho", test_bc_sheath_gk_3x2v_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_bc_sheath_gk_1x2v_dev", test_bc_sheath_gk_1x2v_dev },
  { "test_bc_sheath_gk_2x2v_dev", test_bc_sheath_gk_2x2v_dev },
  { "test_bc_sheath_gk_3x2v_dev", test_bc_sheath_gk_3x2v_dev },
#endif
  { NULL, NULL },
};
