/* This is the unit test for the sheath BC updater.
  * The test does the following:
  * 1. Initializes the sheath BC updater with a given configuration space dimension, edge.
  * 2. Project a maxwellian distribution function with input parameters to the edge of the domain (skin cells).
  * 3. Apply the sheath BC updater to the distribution function.
  * 4. Check the ghost values of the distribution function (where the result of the sheath BC is stored).
  * We check that the reflected function is zero in cells that are fully outside the cutoff 
  * and non zero in cells that are fully inside. You can check these numbers by setting pars.verbose to true.
  * 5. Deallocate memory used by the updater.
  * * This is done for each configuration space dimension (1x2v, 2x2v, 3x2v), for both lower and upper edge, and
  * positive, negative potential, for electrons and ions.
  * * The test is also able to verify the interface of the surrogate model for the sheath BC by setting 
  * use_surrogate to true. In this case the check uses the surrogate directly to obtain the vcutsq.
  * * Finally, one can set write_fields to true to write the distribution function and visualize it with postgkyl.
  * */

#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_velocity_map.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_const.h>
#include <float.h>
#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>

// Path to the KANN surrogate model used in surrogate tests.
// Set to NULL (or set surr_test_enabled = false) to skip surrogate tests.
static const char *surr_model_path = "gyrokinetic/data/nn_model/nn_model_sheath_bc_conv_MPE.kann";
static const bool surr_test_enabled = true;

static struct gkyl_array *
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array *a;
  if (use_gpu) {
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  } else {
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  }
  return a;
}

struct test_sheath_ctx {
  struct gkyl_kann_net *model; // KANN surrogate model for sheath BC (NULL if not using surrogate).
  int cdim; // Configuration space dimension.
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM]; // Grid extents.
  int cells[GKYL_MAX_DIM]; // Number of cells.
  double mass; // Species mass.
  double charge; // Species charge.
  double dens; // Reference density.
  double upar; // Mean parallel drift speed.
  double temp; // Reference temperature.
  double B0; // Reference magnetic field.
  double phi_mpe; // ES potential at the sheath entrance.
  double phi_wall; // ES potential at the wall.
  double impact_angle; // Angle [rad] of incidence if B to the plate (for surrogate sheath).
  double x0; // x-center
  double y0; // y-center
  double z0; // z-center
  double sigmax; // Width of distribution in x.
  double sigmay; // Width of distribution in y.
  double sigmaz; // Width of distribution in z.
  const char *surrogate_model_path; // Path to the .kann model file (NULL if not using surrogate).
  bool verbose; // Whether to print detailed output.
};

void
surr_interpf(const float *vcut, const double *mu_new, int n, double mu_ref, double *out)
{
  const double mu_grid[20] = {0.00, 0.02, 0.08, 0.18, 0.32, 0.50, 0.72, 0.98, 1.28, 1.62,
                              2.00, 2.42, 2.88, 3.38, 3.92, 4.50, 5.12, 5.78, 6.48, 7.22};
  const int ng = 20;
  for (int i = 0; i < n; i++) {
    double mu = mu_new[i] / mu_ref;
    if (mu <= mu_grid[0]) {
      out[i] = vcut[0];
      continue;
    }
    if (mu >= mu_grid[ng - 1]) {
      out[i] = vcut[ng - 1];
      continue;
    }
    // binary search for the bracketing interval
    int lo = 0, hi = ng - 1;
    while (hi - lo > 1) {
      int mid = (lo + hi) >> 1;
      if (mu_grid[mid] <= mu) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    double t = (mu - mu_grid[lo]) / (mu_grid[hi] - mu_grid[lo]);
    out[i] = vcut[lo] + t * (vcut[hi] - vcut[lo]);
  }
}

void
eval_func_vcutsq(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct test_sheath_ctx *pars = ctx;
  double mu = xn[pars->cdim - 1]; // Mu is the last configuration/velocity dimension mapped.

  if (pars->model) {
    float alphadeg = pars->impact_angle * 180.0 / GKYL_PI;
    float gamma = sqrt(pars->mass * fmax(pars->dens, 0.0) / GKYL_EPSILON0) / pars->B0;
    float phi = -pars->charge * (pars->phi_mpe - pars->phi_wall) / pars->temp;
    struct gkyl_kn_vec *inp = gkyl_kn_vec_new(1, 3);
    inp->vals[0][0] = alphadeg;
    inp->vals[0][1] = gamma;
    inp->vals[0][2] = phi;
    struct gkyl_kn_vec *out = gkyl_kn_vec_new(1, 20);
    gkyl_kann_net_apply(pars->model, inp, out);

    // Interp the output of the surrogate to get vcutsq corresponding to input mu.
    double mu_ref = pars->temp / pars->B0; // Normalization for mu used in surrogate training.
    double vcut_norm = 0.0;
    surr_interpf(out->data, &mu, 1, mu_ref, &vcut_norm);

    double vt_sq = pars->temp / pars->mass;
    fout[0] = (vcut_norm * vcut_norm) * vt_sq;

    // Clean up temporary neural net vectors to prevent severe memory leaks
    gkyl_kn_vec_release(inp);
    gkyl_kn_vec_release(out);

  } else {
    // Physical analytical calculation of vcutsq = 2 * q * delta_phi / m
    double delta_phi = pars->phi_mpe - pars->phi_wall;
    double vcutsq_0 = fmax(0.0, -2.0 * pars->charge * delta_phi / pars->mass);
    double Lmu = pars->temp / pars->B0;
    if (Lmu == 0.0) {
      Lmu = 1.0;
    }

    // Smooth mu variation to test mu-dependent evaluation machinery
    fout[0] = vcutsq_0 * (1.0 + 0.1 * exp(-mu / Lmu));
  }
}

void
eval_func_1x2v(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  double vpar = xn[1], mu = xn[2];
  struct test_sheath_ctx *pars = ctx;

  double B0 = pars->B0;
  double m = pars->mass;
  double upar = pars->upar;
  double T = pars->temp;

  double envelope = 1;
  fout[0] = exp(-(m * pow(vpar - upar, 2) / 2.0 + mu * B0) / T) * envelope;
}

void
eval_func_2x2v(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  double vpar = xn[2], mu = xn[3];
  struct test_sheath_ctx *pars = ctx;

  double B0 = pars->B0;
  double m = pars->mass;
  double upar = pars->upar;
  double T = pars->temp;

  double envelope = 1;
  fout[0] = exp(-(m * pow(vpar - upar, 2) / 2.0 + mu * B0) / T) * envelope;
}

void
eval_func_3x2v(double t, const double *xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double vpar = xn[3], mu = xn[4];
  struct test_sheath_ctx *pars = ctx;

  double B0 = pars->B0;
  double m = pars->mass;
  double upar = pars->upar;
  double T = pars->temp;

  double envelope = 1;
  fout[0] = exp(-(m * pow(vpar - upar, 2) / 2.0 + mu * B0) / T) * envelope;
}

// Checks that the distribution function values in the ghost cells are set to
// 0 for velocities beyond the cutoff velocity squared and >0 otherwise.
void
check_function(
  struct test_sheath_ctx *pars, struct gkyl_array *distf_ho, struct gkyl_rect_grid grid,
  struct gkyl_range ghost_r, enum gkyl_edge_loc edge
)
{
  double delta_phi = pars->phi_mpe - pars->phi_wall;
  int cdim = pars->cdim;

  int num_cells = 0;
  int num_zero_cells = 0;
  int num_zero_cells_expected = 0;
  int num_uknown_cells = 0;

  // Sign of charge * delta_phi determines if there is any reflection or full absorption
  double qphi_sign = pars->charge * delta_phi > 0 ? 1.0 : -1.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ghost_r);
  while (gkyl_range_iter_next(&iter)) {
    int *idx_g = iter.idx;
    long linidx_g = gkyl_range_idx(&ghost_r, idx_g);
    const double *distf_c = gkyl_array_cfetch(distf_ho, linidx_g);
    double tol = 1e-12;
    double ref_val = 0.0;

    // Get cell boundaries along vpar.
    double xc[GKYL_MAX_DIM] = {0};
    gkyl_rect_grid_cell_center(&grid, idx_g, xc);
    double cell_lower_vpar = xc[cdim] - 0.5 * grid.dx[cdim];
    double cell_upper_vpar = xc[cdim] + 0.5 * grid.dx[cdim];

    // Compute the absolute velocity bounds of the cell (handle zero crossing safely)
    double vr_min = (cell_lower_vpar * cell_upper_vpar < 0.0) ?
                      0.0 :
                      fmin(fabs(cell_lower_vpar), fabs(cell_upper_vpar));
    double vr_max = fmax(fabs(cell_lower_vpar), fabs(cell_upper_vpar));

    // Square the cell velocities to safely match vcutsq dimensions
    double vr_min_sq = vr_min * vr_min;
    double vr_max_sq = vr_max * vr_max;

    // Compute the vpar cut radius squared at the lower and upper mu boundaries
    double cell_lower_mu = xc[cdim + 1] - 0.5 * grid.dx[cdim + 1];
    double cell_upper_mu = xc[cdim + 1] + 0.5 * grid.dx[cdim + 1];
    double vcut_coord_lower[GKYL_MAX_DIM], vcut_coord_upper[GKYL_MAX_DIM];

    // Map dimensions correctly
    for (int d = 0; d < cdim - 1; d++) {
      vcut_coord_lower[d] = xc[d];
      vcut_coord_upper[d] = xc[d];
    }
    vcut_coord_lower[cdim - 1] = cell_lower_mu;
    vcut_coord_upper[cdim - 1] = cell_upper_mu;

    double vcutsq_lower[1], vcutsq_upper[1];
    eval_func_vcutsq(0, vcut_coord_lower, vcutsq_lower, pars);
    eval_func_vcutsq(0, vcut_coord_upper, vcutsq_upper, pars);

    double vcutsq_min = fmin(vcutsq_lower[0], vcutsq_upper[0]);
    double vcutsq_max = fmax(vcutsq_lower[0], vcutsq_upper[0]);

    num_cells++;
    num_zero_cells += distf_c[0] > tol ? 0 : 1;

    if ((vr_min_sq >= vcutsq_max) || (qphi_sign > 0)) {
      // Cell is fully outside cutoff or under full absorption -> expect exact zeros
      num_zero_cells_expected++;
      for (int k = 0; k < distf_ho->ncomp; k++) {
        TEST_CHECK(gkyl_compare(distf_c[k], ref_val, tol));
        TEST_MSG(
          "Expected %.9e | Got: %.9e at idx=%d,%d,%d\n", ref_val, distf_c[k], idx_g[0], idx_g[1],
          idx_g[2]
        );
      }
    } else if (((vcutsq_min < vr_max_sq) && (vcutsq_min > vr_min_sq)) ||
               ((vcutsq_max > vr_min_sq) && (vcutsq_max < vr_max_sq))) {
      // Cell is partially cut off by the sheath boundaries
      num_uknown_cells++;
    } else {
      // Cell is fully inside the cutoff range -> average cell value must be non-zero
      double cell_avg = distf_c[0];
      TEST_CHECK(cell_avg > ref_val);
      TEST_MSG(
        "Expected > %.9e | Got: %.9e at idx=%d,%d,%d\n", ref_val, cell_avg, idx_g[0], idx_g[1],
        idx_g[2]
      );
    }
  }
}

void
write_out_fields(
  struct test_sheath_ctx *pars, int cdim, int vdim, enum gkyl_edge_loc edge, bool use_gpu,
  struct gkyl_array *distf_ho, struct gkyl_array *phi_ho, struct gkyl_array *phiw_ho,
  struct gkyl_rect_grid grid_ext, struct gkyl_range local_ext, struct gkyl_rect_grid grid_conf,
  struct gkyl_range local_conf, struct gkyl_basis *basis, struct gkyl_basis basis_conf
)
{
  struct gkyl_rect_grid grid_norm_ext = grid_ext;
  double vt = sqrt(pars->temp / pars->mass);
  grid_norm_ext.dx[cdim] /= vt;
  grid_norm_ext.lower[cdim] /= vt;
  grid_norm_ext.upper[cdim] /= vt;
  grid_norm_ext.dx[cdim + 1] *= pars->B0 / pars->temp;
  grid_norm_ext.lower[cdim + 1] *= pars->B0 / pars->temp;
  grid_norm_ext.upper[cdim + 1] *= pars->B0 / pars->temp;

  struct gkyl_msgpack_map_elem io_meta[] = {
    {.key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = basis->poly_order},
    {.key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = basis->id}
  };
  int io_meta_len = sizeof(io_meta) / sizeof(io_meta[0]);
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create(io_meta_len, io_meta);

  struct gkyl_msgpack_map_elem io_meta_conf[] = {
    {.key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = basis_conf.poly_order},
    {.key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = basis_conf.id}
  };
  int io_meta_conf_len = sizeof(io_meta_conf) / sizeof(io_meta_conf[0]);
  struct gkyl_msgpack_data *mt_conf = gkyl_msgpack_create(io_meta_conf_len, io_meta_conf);

  char fname[256];
  const char *fmt = "bc_sheath_%dx%dv_%s_%s_%s.gkyl";
  snprintf(
    fname, sizeof(fname), fmt, cdim, vdim, edge == GKYL_LOWER_EDGE ? "lower" : "upper",
    use_gpu ? "gpu" : "cpu", "distf_out"
  );
  gkyl_grid_sub_array_write(&grid_norm_ext, &local_ext, mt, distf_ho, fname);
  snprintf(
    fname, sizeof(fname), fmt, cdim, vdim, edge == GKYL_LOWER_EDGE ? "lower" : "upper",
    use_gpu ? "gpu" : "cpu", "phi_mpe"
  );
  gkyl_grid_sub_array_write(&grid_conf, &local_conf, mt_conf, phi_ho, fname);
  snprintf(
    fname, sizeof(fname), fmt, cdim, vdim, edge == GKYL_LOWER_EDGE ? "lower" : "upper",
    use_gpu ? "gpu" : "cpu", "phi_wall"
  );
  gkyl_grid_sub_array_write(&grid_conf, &local_conf, mt_conf, phiw_ho, fname);

  gkyl_msgpack_data_release(mt);
  gkyl_msgpack_data_release(mt_conf);
}

void
test_bc_sheath_gyrokinetic_1x2v(
  struct test_sheath_ctx *pars, enum gkyl_edge_loc edge, bool use_surrogate, bool write_fields,
  bool use_gpu
)
{
  double ms = pars->mass;
  double qs = pars->charge;
  double ns = pars->dens;
  double Ts = pars->temp;
  double B0 = pars->B0;
  double phi_mpe = pars->phi_mpe;
  double phi_wall = pars->phi_wall;
  double impact_angle = pars->impact_angle;
  double vt = sqrt(Ts / ms);
  double mu0 = Ts / B0;

  int poly_order = 1;
  double lower[] = {-1.0, -3.0 * vt, 0.};
  double upper[] = {1.0, 3.0 * vt, 5.0 * mu0};
  int vdim = 2;
  int ndim = sizeof(lower) / sizeof(lower[0]);
  int cdim = ndim - vdim;
  int dir = cdim - 1;

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d = 0; d < cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = pars->cells[d];
  }
  double lower_vel[vdim], upper_vel[vdim];
  int cells_vel[vdim];
  for (int d = 0; d < vdim; d++) {
    lower_vel[d] = lower[cdim + d];
    upper_vel[d] = upper[cdim + d];
    cells_vel[d] = pars->cells[cdim + d];
  }

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, pars->cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  struct gkyl_rect_grid grid_vel;
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

  struct gkyl_basis *basis;
  basis = use_gpu ? gkyl_cart_modal_gkhybrid_cu_dev_new(cdim, vdim) :
                    gkyl_cart_modal_gkhybrid_new(cdim, vdim);

  struct gkyl_basis basis_ho;
  gkyl_cart_modal_gkhybrid(&basis_ho, cdim, vdim);

  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  int ghost_conf[] = {1, 1, 1};
  struct gkyl_range local_conf, local_conf_ext;
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_conf_ext, &local_conf);

  int ghost_vel[] = {0, 0};
  struct gkyl_range local_vel, local_vel_ext;
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_vel_ext, &local_vel);

  int ghost[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; d++) {
    ghost[d] = ghost_conf[d];
  }
  for (int d = 0; d < vdim; d++) {
    ghost[cdim + d] = ghost_vel[d];
  }
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_mapc2p_inp c2p_in = {};
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(
    c2p_in, grid, grid_vel, local, local_ext, local_vel, local_vel_ext, use_gpu
  );

  double lower_ext[ndim], upper_ext[ndim];
  int cells_ext[ndim];
  for (int d = 0; d < ndim; d++) {
    double dx = (upper[d] - lower[d]) / pars->cells[d];
    lower_ext[d] = lower[d] - dx * ghost[d];
    upper_ext[d] = upper[d] + dx * ghost[d];
    cells_ext[d] = pars->cells[d] + 2 * ghost[d];
  }
  struct gkyl_rect_grid grid_ext;
  gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);

  struct gkyl_range skin_r, ghost_r;
  gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

  struct gkyl_array *distf = mkarr(use_gpu, basis_ho.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu ? mkarr(false, distf->ncomp, distf->size) :
                                          gkyl_array_acquire(distf);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &grid,
    .basis = &basis_ho,
    .num_ret_vals = 1,
    .eval = eval_func_1x2v,
    .ctx = pars,
  });
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  double dgnormc = pow(sqrt(2.0), cdim);
  struct gkyl_array *phi = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phi_ho = use_gpu ? mkarr(false, phi->ncomp, phi->size) :
                                        gkyl_array_acquire(phi);
  gkyl_array_shiftc(phi_ho, phi_mpe * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phi, phi_ho);

  struct gkyl_array *phiw = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phiw_ho = use_gpu ? mkarr(false, phiw->ncomp, phiw->size) :
                                         gkyl_array_acquire(phiw);
  gkyl_array_shiftc(phiw_ho, phi_wall * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phiw, phiw_ho);

  struct gkyl_bc_sheath_gyrokinetic *bcsheath = gkyl_bc_sheath_gyrokinetic_new(
    dir, edge, basis, &skin_r, &ghost_r, gvm, cdim, 2. * qs / ms, use_surrogate,
    pars->surrogate_model_path, NULL, NULL, use_gpu
  ); // No vcutsq diagnostic output in this test.

  if (use_surrogate) {
    pars->model = gkyl_kann_net_load(pars->surrogate_model_path, false);
  }

  int vcutsq_dim = cdim - 1 + vdim - 1;
  struct gkyl_basis vcutsq_basis;
  gkyl_cart_modal_serendip(&vcutsq_basis, vcutsq_dim, poly_order);
  struct gkyl_range vcutsq_local;
  int nc_lower[vcutsq_dim], nc_upper[vcutsq_dim];
  for (int d = 0; d < cdim - 1; d++) {
    nc_lower[d] = skin_r.lower[d];
    nc_upper[d] = skin_r.upper[d];
  }
  nc_lower[vcutsq_dim - 1] = skin_r.lower[skin_r.ndim - 1];
  nc_upper[vcutsq_dim - 1] = skin_r.upper[skin_r.ndim - 1];
  gkyl_range_init(&vcutsq_local, vcutsq_basis.ndim, nc_lower, nc_upper);
  struct gkyl_array *vcutsq = mkarr(use_gpu, vcutsq_basis.num_basis, vcutsq_local.volume);
  struct gkyl_array *vcutsq_ho = use_gpu ? mkarr(false, vcutsq->ncomp, vcutsq->size) :
                                           gkyl_array_acquire(vcutsq);

  struct gkyl_rect_grid vcut_grid;
  double vcut_lower[cdim], vcut_upper[cdim];
  for (int d = 0; d < cdim - 1; d++) {
    vcut_lower[d] = lower[d];
    vcut_upper[d] = upper[d];
  }
  vcut_lower[cdim - 1] = lower[cdim + vdim - 1];
  vcut_upper[cdim - 1] = upper[cdim + vdim - 1];
  gkyl_rect_grid_init(&vcut_grid, cdim, vcut_lower, vcut_upper, pars->cells);
  gkyl_proj_on_basis *projVcut = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &vcut_grid,
    .basis = &vcutsq_basis,
    .num_ret_vals = 1,
    .eval = eval_func_vcutsq,
    .ctx = pars,
  });
  gkyl_proj_on_basis_advance(projVcut, 0.0, &vcutsq_local, vcutsq_ho);
  gkyl_array_copy(vcutsq, vcutsq_ho);

  gkyl_bc_sheath_gyrokinetic_set_vcutsq(bcsheath, vcutsq);

  if (use_surrogate) {
    struct gkyl_array *density = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *density_ho = use_gpu ? mkarr(false, density->ncomp, density->size) :
                                              gkyl_array_acquire(density);
    gkyl_array_shiftc(density_ho, ns * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(density, density_ho);

    struct gkyl_array *temperature = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *temperature_ho = use_gpu ?
                                          mkarr(false, temperature->ncomp, temperature->size) :
                                          gkyl_array_acquire(temperature);
    gkyl_array_shiftc(temperature_ho, Ts * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(temperature, temperature_ho);

    struct gkyl_array *bmag = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *bmag_ho = use_gpu ? mkarr(false, bmag->ncomp, bmag->size) :
                                           gkyl_array_acquire(bmag);
    gkyl_array_shiftc(bmag_ho, B0 * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(bmag, bmag_ho);

    struct gkyl_array *bimpact_angle = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *bimpact_angle_ho =
      use_gpu ? mkarr(false, bimpact_angle->ncomp, bimpact_angle->size) :
                gkyl_array_acquire(bimpact_angle);
    gkyl_array_shiftc(bimpact_angle_ho, impact_angle * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(bimpact_angle, bimpact_angle_ho);

    gkyl_bc_sheath_gyrokinetic_update_vcutsq(
      bcsheath, phi, phiw, density, temperature, bmag, bimpact_angle, &local_conf
    );

    gkyl_array_release(density);
    gkyl_array_release(density_ho);
    gkyl_array_release(temperature);
    gkyl_array_release(temperature_ho);
    gkyl_array_release(bmag);
    gkyl_array_release(bmag_ho);
    gkyl_array_release(bimpact_angle);
    gkyl_array_release(bimpact_angle_ho);
  }

  gkyl_bc_sheath_gyrokinetic_advance(bcsheath, distf, &local_conf);
  gkyl_array_copy(distf_ho, distf);

  check_function(pars, distf_ho, grid, ghost_r, edge);

  if (write_fields) {
    write_out_fields(
      pars, cdim, vdim, edge, use_gpu, distf_ho, phi_ho, phiw_ho, grid_ext, local_ext, grid_conf,
      local_conf, &basis_ho, basis_conf
    );
  }

  // Clean up and prevent memory leaks.
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(phiw);
  gkyl_array_release(phiw_ho);
  if (use_gpu) {
    gkyl_cart_modal_basis_release_cu(basis);
  } else {
    gkyl_cart_modal_basis_release(basis);
  }

  gkyl_velocity_map_release(gvm);
  gkyl_bc_sheath_gyrokinetic_release(bcsheath);
  gkyl_proj_on_basis_release(projVcut);
  gkyl_array_release(vcutsq);
  gkyl_array_release(vcutsq_ho);

  if (pars->model) {
    gkyl_kann_net_release(pars->model);
    pars->model = NULL;
  }
}

void
test_bc_sheath_gyrokinetic_2x2v(
  struct test_sheath_ctx *pars, enum gkyl_edge_loc edge, bool use_surrogate, bool write_fields,
  bool use_gpu
)
{
  double ms = pars->mass;
  double qs = pars->charge;
  double ns = pars->dens;
  double Ts = pars->temp;
  double B0 = pars->B0;
  double phi_mpe = pars->phi_mpe;
  double phi_wall = pars->phi_wall;
  double impact_angle = pars->impact_angle;
  double vt = sqrt(Ts / ms);
  double mu0 = Ts / B0;

  int poly_order = 1;
  double lower[] = {0.0, -1.0, -3.0 * vt, 0.};
  double upper[] = {1.0, 1.0, 3.0 * vt, 5 * mu0};
  int vdim = 2;
  int ndim = sizeof(lower) / sizeof(lower[0]);
  int cdim = ndim - vdim;
  int dir = cdim - 1;

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d = 0; d < cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = pars->cells[d];
  }
  double lower_vel[vdim], upper_vel[vdim];
  int cells_vel[vdim];
  for (int d = 0; d < vdim; d++) {
    lower_vel[d] = lower[cdim + d];
    upper_vel[d] = upper[cdim + d];
    cells_vel[d] = pars->cells[cdim + d];
  }

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, pars->cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  struct gkyl_rect_grid grid_vel;
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

  struct gkyl_basis *basis;
  basis = use_gpu ? gkyl_cart_modal_gkhybrid_cu_dev_new(cdim, vdim) :
                    gkyl_cart_modal_gkhybrid_new(cdim, vdim);

  struct gkyl_basis basis_ho;
  gkyl_cart_modal_gkhybrid(&basis_ho, cdim, vdim);

  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  int ghost_conf[] = {1, 1, 1};
  struct gkyl_range local_conf, local_conf_ext;
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_conf_ext, &local_conf);

  int ghost_vel[] = {0, 0};
  struct gkyl_range local_vel, local_vel_ext;
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_vel_ext, &local_vel);

  int ghost[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; d++) {
    ghost[d] = ghost_conf[d];
  }
  for (int d = 0; d < vdim; d++) {
    ghost[cdim + d] = ghost_vel[d];
  }
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_mapc2p_inp c2p_in = {};
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(
    c2p_in, grid, grid_vel, local, local_ext, local_vel, local_vel_ext, use_gpu
  );

  double lower_ext[ndim], upper_ext[ndim];
  int cells_ext[ndim];
  for (int d = 0; d < ndim; d++) {
    double dx = (upper[d] - lower[d]) / pars->cells[d];
    lower_ext[d] = lower[d] - dx * ghost[d];
    upper_ext[d] = upper[d] + dx * ghost[d];
    cells_ext[d] = pars->cells[d] + 2 * ghost[d];
  }
  struct gkyl_rect_grid grid_ext;
  gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);

  struct gkyl_range skin_r, ghost_r;
  gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

  struct gkyl_array *distf = mkarr(use_gpu, basis_ho.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu ? mkarr(false, distf->ncomp, distf->size) :
                                          gkyl_array_acquire(distf);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &grid,
    .basis = &basis_ho,
    .num_ret_vals = 1,
    .eval = eval_func_2x2v,
    .ctx = pars,
  });
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  double dgnormc = pow(sqrt(2.0), cdim);
  struct gkyl_array *phi = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phi_ho = use_gpu ? mkarr(false, phi->ncomp, phi->size) :
                                        gkyl_array_acquire(phi);
  gkyl_array_shiftc(phi_ho, phi_mpe * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phi, phi_ho);

  struct gkyl_array *phiw = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phiw_ho = use_gpu ? mkarr(false, phiw->ncomp, phiw->size) :
                                         gkyl_array_acquire(phiw);
  gkyl_array_shiftc(phiw_ho, phi_wall * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phiw, phiw_ho);

  struct gkyl_bc_sheath_gyrokinetic *bcsheath = gkyl_bc_sheath_gyrokinetic_new(
    dir, edge, basis, &skin_r, &ghost_r, gvm, cdim, 2. * qs / ms, use_surrogate,
    pars->surrogate_model_path, NULL, NULL, use_gpu
  ); // No vcutsq diagnostic output in this test.

  if (use_surrogate) {
    pars->model = gkyl_kann_net_load(pars->surrogate_model_path, false);
  }

  int vcutsq_dim = cdim - 1 + vdim - 1;
  struct gkyl_basis vcutsq_basis;
  gkyl_cart_modal_serendip(&vcutsq_basis, vcutsq_dim, poly_order);
  struct gkyl_range vcutsq_local;
  int nc_lower[vcutsq_dim], nc_upper[vcutsq_dim];
  for (int d = 0; d < cdim - 1; d++) {
    nc_lower[d] = skin_r.lower[d];
    nc_upper[d] = skin_r.upper[d];
  }
  nc_lower[vcutsq_dim - 1] = skin_r.lower[skin_r.ndim - 1];
  nc_upper[vcutsq_dim - 1] = skin_r.upper[skin_r.ndim - 1];
  gkyl_range_init(&vcutsq_local, vcutsq_basis.ndim, nc_lower, nc_upper);
  struct gkyl_array *vcutsq = mkarr(use_gpu, vcutsq_basis.num_basis, vcutsq_local.volume);
  struct gkyl_array *vcutsq_ho = use_gpu ? mkarr(false, vcutsq->ncomp, vcutsq->size) :
                                           gkyl_array_acquire(vcutsq);

  struct gkyl_rect_grid vcut_grid;
  double vcut_lower[cdim], vcut_upper[cdim];
  for (int d = 0; d < cdim - 1; d++) {
    vcut_lower[d] = lower[d];
    vcut_upper[d] = upper[d];
  }
  vcut_lower[cdim - 1] = lower[cdim + vdim - 1];
  vcut_upper[cdim - 1] = upper[cdim + vdim - 1];
  gkyl_rect_grid_init(&vcut_grid, cdim, vcut_lower, vcut_upper, pars->cells);
  gkyl_proj_on_basis *projVcut = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &vcut_grid,
    .basis = &vcutsq_basis,
    .num_ret_vals = 1,
    .eval = eval_func_vcutsq,
    .ctx = pars,
  });
  gkyl_proj_on_basis_advance(projVcut, 0.0, &vcutsq_local, vcutsq_ho);
  gkyl_array_copy(vcutsq, vcutsq_ho);

  gkyl_bc_sheath_gyrokinetic_set_vcutsq(bcsheath, vcutsq);

  if (use_surrogate) {
    struct gkyl_array *density = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *density_ho = use_gpu ? mkarr(false, density->ncomp, density->size) :
                                              gkyl_array_acquire(density);
    gkyl_array_shiftc(density_ho, ns * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(density, density_ho);

    struct gkyl_array *temperature = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *temperature_ho = use_gpu ?
                                          mkarr(false, temperature->ncomp, temperature->size) :
                                          gkyl_array_acquire(temperature);
    gkyl_array_shiftc(temperature_ho, Ts * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(temperature, temperature_ho);

    struct gkyl_array *bmag = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *bmag_ho = use_gpu ? mkarr(false, bmag->ncomp, bmag->size) :
                                           gkyl_array_acquire(bmag);
    gkyl_array_shiftc(bmag_ho, B0 * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(bmag, bmag_ho);

    struct gkyl_array *bimpact_angle = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *bimpact_angle_ho =
      use_gpu ? mkarr(false, bimpact_angle->ncomp, bimpact_angle->size) :
                gkyl_array_acquire(bimpact_angle);
    gkyl_array_shiftc(bimpact_angle_ho, impact_angle * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(bimpact_angle, bimpact_angle_ho);

    gkyl_bc_sheath_gyrokinetic_update_vcutsq(
      bcsheath, phi, phiw, density, temperature, bmag, bimpact_angle, &local_conf
    );

    gkyl_array_release(density);
    gkyl_array_release(density_ho);
    gkyl_array_release(temperature);
    gkyl_array_release(temperature_ho);
    gkyl_array_release(bmag);
    gkyl_array_release(bmag_ho);
    gkyl_array_release(bimpact_angle);
    gkyl_array_release(bimpact_angle_ho);
  }

  gkyl_bc_sheath_gyrokinetic_advance(bcsheath, distf, &local_conf);
  gkyl_array_copy(distf_ho, distf);

  check_function(pars, distf_ho, grid, ghost_r, edge);

  if (write_fields) {
    write_out_fields(
      pars, cdim, vdim, edge, use_gpu, distf_ho, phi_ho, phiw_ho, grid_ext, local_ext, grid_conf,
      local_conf, &basis_ho, basis_conf
    );
  }

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(phiw);
  gkyl_array_release(phiw_ho);
  if (use_gpu) {
    gkyl_cart_modal_basis_release_cu(basis);
  } else {
    gkyl_cart_modal_basis_release(basis);
  }

  gkyl_velocity_map_release(gvm);
  gkyl_bc_sheath_gyrokinetic_release(bcsheath);
  gkyl_proj_on_basis_release(projVcut);
  gkyl_array_release(vcutsq);
  gkyl_array_release(vcutsq_ho);

  if (pars->model) {
    gkyl_kann_net_release(pars->model);
    pars->model = NULL;
  }
}

void
test_bc_sheath_gyrokinetic_3x2v(
  struct test_sheath_ctx *pars, enum gkyl_edge_loc edge, bool use_surrogate, bool write_fields,
  bool use_gpu
)
{
  double ms = pars->mass;
  double qs = pars->charge;
  double ns = pars->dens;
  double Ts = pars->temp;
  double B0 = pars->B0;
  double phi_mpe = pars->phi_mpe;
  double phi_wall = pars->phi_wall;
  double impact_angle = pars->impact_angle;
  double vt = sqrt(Ts / ms);
  double mu0 = Ts / B0;

  int poly_order = 1;
  double lower[] = {0.0, -2.0, -1.0, -3.0 * vt, 0.};
  double upper[] = {1.0, 2.0, 1.0, 3.0 * vt, 5 * mu0};
  int vdim = 2;
  int ndim = sizeof(lower) / sizeof(lower[0]);
  int cdim = ndim - vdim;
  int dir = cdim - 1;

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d = 0; d < cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = pars->cells[d];
  }
  double lower_vel[vdim], upper_vel[vdim];
  int cells_vel[vdim];
  for (int d = 0; d < vdim; d++) {
    lower_vel[d] = lower[cdim + d];
    upper_vel[d] = upper[cdim + d];
    cells_vel[d] = pars->cells[cdim + d];
  }

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, pars->cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  struct gkyl_rect_grid grid_vel;
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

  struct gkyl_basis *basis;
  basis = use_gpu ? gkyl_cart_modal_gkhybrid_cu_dev_new(cdim, vdim) :
                    gkyl_cart_modal_gkhybrid_new(cdim, vdim);

  struct gkyl_basis basis_ho;
  gkyl_cart_modal_gkhybrid(&basis_ho, cdim, vdim);

  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  int ghost_conf[] = {1, 1, 1};
  struct gkyl_range local_conf, local_conf_ext;
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_conf_ext, &local_conf);

  int ghost_vel[] = {0, 0};
  struct gkyl_range local_vel, local_vel_ext;
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_vel_ext, &local_vel);

  int ghost[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; d++) {
    ghost[d] = ghost_conf[d];
  }
  for (int d = 0; d < vdim; d++) {
    ghost[cdim + d] = ghost_vel[d];
  }
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_mapc2p_inp c2p_in = {};
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(
    c2p_in, grid, grid_vel, local, local_ext, local_vel, local_vel_ext, use_gpu
  );

  double lower_ext[ndim], upper_ext[ndim];
  int cells_ext[ndim];
  for (int d = 0; d < ndim; d++) {
    double dx = (upper[d] - lower[d]) / pars->cells[d];
    lower_ext[d] = lower[d] - dx * ghost[d];
    upper_ext[d] = upper[d] + dx * ghost[d];
    cells_ext[d] = pars->cells[d] + 2 * ghost[d];
  }
  struct gkyl_rect_grid grid_ext;
  gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);

  struct gkyl_range skin_r, ghost_r;
  gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

  struct gkyl_array *distf = mkarr(use_gpu, basis_ho.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu ? mkarr(false, distf->ncomp, distf->size) :
                                          gkyl_array_acquire(distf);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &grid,
    .basis = &basis_ho,
    .num_ret_vals = 1,
    .eval = eval_func_3x2v,
    .ctx = pars,
  });
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  double dgnormc = pow(sqrt(2.0), cdim);
  struct gkyl_array *phi = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phi_ho = use_gpu ? mkarr(false, phi->ncomp, phi->size) :
                                        gkyl_array_acquire(phi);
  gkyl_array_shiftc(phi_ho, phi_mpe * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phi, phi_ho);

  struct gkyl_array *phiw = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phiw_ho = use_gpu ? mkarr(false, phiw->ncomp, phiw->size) :
                                         gkyl_array_acquire(phiw);
  gkyl_array_shiftc(phiw_ho, phi_wall * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phiw, phiw_ho);

  struct gkyl_bc_sheath_gyrokinetic *bcsheath = gkyl_bc_sheath_gyrokinetic_new(
    dir, edge, basis, &skin_r, &ghost_r, gvm, cdim, 2. * qs / ms, use_surrogate,
    pars->surrogate_model_path, NULL, NULL, use_gpu
  ); // No vcutsq diagnostic output in this test.

  if (use_surrogate) {
    pars->model = gkyl_kann_net_load(pars->surrogate_model_path, false);
  }

  int vcutsq_dim = cdim - 1 + vdim - 1;
  struct gkyl_basis vcutsq_basis;
  gkyl_cart_modal_serendip(&vcutsq_basis, vcutsq_dim, poly_order);
  struct gkyl_range vcutsq_local;
  int nc_lower[vcutsq_dim], nc_upper[vcutsq_dim];
  for (int d = 0; d < cdim - 1; d++) {
    nc_lower[d] = skin_r.lower[d];
    nc_upper[d] = skin_r.upper[d];
  }
  nc_lower[vcutsq_dim - 1] = skin_r.lower[skin_r.ndim - 1];
  nc_upper[vcutsq_dim - 1] = skin_r.upper[skin_r.ndim - 1];
  gkyl_range_init(&vcutsq_local, vcutsq_basis.ndim, nc_lower, nc_upper);
  struct gkyl_array *vcutsq = mkarr(use_gpu, vcutsq_basis.num_basis, vcutsq_local.volume);
  struct gkyl_array *vcutsq_ho = use_gpu ? mkarr(false, vcutsq->ncomp, vcutsq->size) :
                                           gkyl_array_acquire(vcutsq);

  struct gkyl_rect_grid vcut_grid;
  double vcut_lower[cdim], vcut_upper[cdim];
  for (int d = 0; d < cdim - 1; d++) {
    vcut_lower[d] = lower[d];
    vcut_upper[d] = upper[d];
  }
  vcut_lower[cdim - 1] = lower[cdim + vdim - 1];
  vcut_upper[cdim - 1] = upper[cdim + vdim - 1];
  gkyl_rect_grid_init(&vcut_grid, cdim, vcut_lower, vcut_upper, pars->cells);
  gkyl_proj_on_basis *projVcut = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &vcut_grid,
    .basis = &vcutsq_basis,
    .num_ret_vals = 1,
    .eval = eval_func_vcutsq,
    .ctx = pars,
  });
  gkyl_proj_on_basis_advance(projVcut, 0.0, &vcutsq_local, vcutsq_ho);
  gkyl_array_copy(vcutsq, vcutsq_ho);

  gkyl_bc_sheath_gyrokinetic_set_vcutsq(bcsheath, vcutsq);

  if (use_surrogate) {
    struct gkyl_array *density = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *density_ho = use_gpu ? mkarr(false, density->ncomp, density->size) :
                                              gkyl_array_acquire(density);
    gkyl_array_shiftc(density_ho, ns * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(density, density_ho);

    struct gkyl_array *temperature = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *temperature_ho = use_gpu ?
                                          mkarr(false, temperature->ncomp, temperature->size) :
                                          gkyl_array_acquire(temperature);
    gkyl_array_shiftc(temperature_ho, Ts * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(temperature, temperature_ho);

    struct gkyl_array *bmag = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *bmag_ho = use_gpu ? mkarr(false, bmag->ncomp, bmag->size) :
                                           gkyl_array_acquire(bmag);
    gkyl_array_shiftc(bmag_ho, B0 * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(bmag, bmag_ho);

    struct gkyl_array *bimpact_angle = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *bimpact_angle_ho =
      use_gpu ? mkarr(false, bimpact_angle->ncomp, bimpact_angle->size) :
                gkyl_array_acquire(bimpact_angle);
    gkyl_array_shiftc(bimpact_angle_ho, impact_angle * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(bimpact_angle, bimpact_angle_ho);

    gkyl_bc_sheath_gyrokinetic_update_vcutsq(
      bcsheath, phi, phiw, density, temperature, bmag, bimpact_angle, &local_conf
    );

    gkyl_array_release(density);
    gkyl_array_release(density_ho);
    gkyl_array_release(temperature);
    gkyl_array_release(temperature_ho);
    gkyl_array_release(bmag);
    gkyl_array_release(bmag_ho);
    gkyl_array_release(bimpact_angle);
    gkyl_array_release(bimpact_angle_ho);
  }

  gkyl_bc_sheath_gyrokinetic_advance(bcsheath, distf, &local_conf);
  gkyl_array_copy(distf_ho, distf);

  check_function(pars, distf_ho, grid, ghost_r, edge);

  if (write_fields) {
    write_out_fields(
      pars, cdim, vdim, edge, use_gpu, distf_ho, phi_ho, phiw_ho, grid_ext, local_ext, grid_conf,
      local_conf, &basis_ho, basis_conf
    );
  }

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(phiw);
  gkyl_array_release(phiw_ho);
  if (use_gpu) {
    gkyl_cart_modal_basis_release_cu(basis);
  } else {
    gkyl_cart_modal_basis_release(basis);
  }

  gkyl_velocity_map_release(gvm);
  gkyl_bc_sheath_gyrokinetic_release(bcsheath);
  gkyl_proj_on_basis_release(projVcut);
  gkyl_array_release(vcutsq);
  gkyl_array_release(vcutsq_ho);

  if (pars->model) {
    gkyl_kann_net_release(pars->model);
    pars->model = NULL;
  }
}

void
test_bc_sheath_gk_1x2v_ho()
{
  if (!surr_test_enabled || !surr_model_path) {
    return;
  }
  bool write_fields;

  struct test_sheath_ctx pars = {
    .cdim = 1,
    .cells = {4, 16, 12},
    .mass = GKYL_ELECTRON_MASS,
    .charge = -GKYL_ELEMENTARY_CHARGE,
    .dens = 1.0e19,
    .upar = 0.0,
    .temp = 5.0 * GKYL_ELEMENTARY_CHARGE,
    .B0 = 1.0,
    .phi_mpe = 15.0,
    .phi_wall = 0.0,
    .impact_angle = 5.0 * GKYL_PI / 180.0,
    .z0 = 0.0,
    .sigmaz = 1.0,
    .surrogate_model_path = surr_model_path,
    .verbose = false,
  };

  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, true, write_fields, false);
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, false);

  write_fields = false;
  pars.phi_mpe *= -1.0;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, false);
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, false);

  write_fields = false;
  pars.charge = GKYL_ELEMENTARY_CHARGE;
  pars.mass = GKYL_PROTON_MASS;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, false);
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, false);

  write_fields = false;
  pars.phi_mpe *= -1.0;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, false);
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, false);
}

void
test_bc_sheath_gk_2x2v_ho()
{
  bool write_fields;

  struct test_sheath_ctx pars = {
    .cdim = 2,
    .cells = {4, 4, 16, 12},
    .mass = GKYL_ELECTRON_MASS,
    .charge = -GKYL_ELEMENTARY_CHARGE,
    .dens = 1.0e19,
    .upar = 0.0,
    .temp = 5.0 * GKYL_ELEMENTARY_CHARGE,
    .B0 = 1.0,
    .phi_mpe = 15.0,
    .phi_wall = 0.0,
    .impact_angle = 5.0 * GKYL_PI / 180.0,
    .x0 = 0.0,
    .sigmax = 1.0,
    .z0 = 0.0,
    .sigmaz = 1.0,
    .surrogate_model_path = surr_model_path,
    .verbose = false,
  };

  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, true, write_fields, false);
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, false);

  pars.phi_mpe *= -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, false);
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, false);

  write_fields = false;
  pars.charge = GKYL_ELEMENTARY_CHARGE;
  pars.mass = GKYL_PROTON_MASS;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, false);
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, false);

  write_fields = false;
  pars.phi_mpe *= -1.0;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, false);
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, false);
}

void
test_bc_sheath_gk_3x2v_ho()
{
  if (!surr_test_enabled || !surr_model_path) {
    return;
  }
  bool write_fields;

  struct test_sheath_ctx pars = {
    .cdim = 3,
    .cells = {4, 4, 4, 16, 12},
    .mass = GKYL_ELECTRON_MASS,
    .charge = -GKYL_ELEMENTARY_CHARGE,
    .dens = 1.0e19,
    .upar = 0.0,
    .temp = 5.0 * GKYL_ELEMENTARY_CHARGE,
    .B0 = 1.0,
    .phi_mpe = 15.0,
    .phi_wall = 0.0,
    .impact_angle = 5.0 * GKYL_PI / 180.0,
    .x0 = 0.0,
    .sigmax = 1.0,
    .y0 = 0.0,
    .sigmay = 1.0,
    .z0 = 0.0,
    .sigmaz = 1.0,
    .surrogate_model_path = surr_model_path,
    .verbose = false,
  };

  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, true, write_fields, false);
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, false);

  write_fields = false;
  pars.phi_mpe *= -1.0;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, false);
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, false);

  write_fields = false;
  pars.charge = GKYL_ELEMENTARY_CHARGE;
  pars.mass = GKYL_PROTON_MASS;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, false);
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, false);

  write_fields = false;
  pars.phi_mpe *= -1.0;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, false);
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, false);
}

#ifdef GKYL_HAVE_CUDA
void
test_bc_sheath_gk_1x2v_dev()
{
  bool write_fields;

  struct test_sheath_ctx pars = {
    .cdim = 1,
    .cells = {4, 16, 12},
    .mass = GKYL_ELECTRON_MASS,
    .charge = -GKYL_ELEMENTARY_CHARGE,
    .dens = 1.0e19,
    .upar = 0.0,
    .temp = 5.0 * GKYL_ELEMENTARY_CHARGE,
    .B0 = 1.0,
    .phi_mpe = 15.0,
    .phi_wall = 0.0,
    .impact_angle = 5.0 * GKYL_PI / 180.0,
    .z0 = 0.0,
    .sigmaz = 1.0,
    .surrogate_model_path = surr_model_path,
    .verbose = false,
  };

  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, true, write_fields, true);
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, true);

  pars.phi_mpe *= -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, true);
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, true);

  pars.charge = GKYL_ELEMENTARY_CHARGE;
  pars.mass = GKYL_PROTON_MASS;
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, true);
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, true);

  pars.phi_mpe *= -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, true);
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, true);
}

void
test_bc_sheath_gk_2x2v_dev()
{
  bool write_fields;

  struct test_sheath_ctx pars = {
    .cdim = 2,
    .cells = {4, 4, 16, 12},
    .mass = GKYL_ELECTRON_MASS,
    .charge = -GKYL_ELEMENTARY_CHARGE,
    .dens = 1.0e19,
    .upar = 0.0,
    .temp = 5.0 * GKYL_ELEMENTARY_CHARGE,
    .B0 = 1.0,
    .phi_mpe = 15.0,
    .phi_wall = 0.0,
    .impact_angle = 5.0 * GKYL_PI / 180.0,
    .x0 = 0.0,
    .sigmax = 1.0,
    .z0 = 0.0,
    .sigmaz = 1.0,
    .surrogate_model_path = surr_model_path,
    .verbose = false,
  };

  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, true, write_fields, true);
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, true);

  pars.phi_mpe *= -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, true);
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, true);

  pars.charge = GKYL_ELEMENTARY_CHARGE;
  pars.mass = GKYL_PROTON_MASS;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, true);
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, true);

  pars.phi_mpe *= -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, true);
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, true);
}

void
test_bc_sheath_gk_3x2v_dev()
{
  bool write_fields;

  struct test_sheath_ctx pars = {
    .cdim = 3,
    .cells = {4, 4, 4, 16, 12},
    .mass = GKYL_ELECTRON_MASS,
    .charge = -GKYL_ELEMENTARY_CHARGE,
    .dens = 1.0e19,
    .upar = 0.0,
    .temp = 5.0 * GKYL_ELEMENTARY_CHARGE,
    .B0 = 1.0,
    .phi_mpe = 15.0,
    .phi_wall = 0.0,
    .impact_angle = 5.0 * GKYL_PI / 180.0,
    .x0 = 0.0,
    .sigmax = 1.0,
    .y0 = 0.0,
    .sigmay = 1.0,
    .z0 = 0.0,
    .sigmaz = 1.0,
    .surrogate_model_path = surr_model_path,
    .verbose = false,
  };

  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, true, write_fields, true);
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, true);

  pars.phi_mpe *= -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, true);
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, true);

  pars.charge = GKYL_ELEMENTARY_CHARGE;
  pars.mass = GKYL_PROTON_MASS;
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, true);
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, true);

  pars.phi_mpe *= -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, false, write_fields, true);
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, false, write_fields, true);
}
#endif

TEST_LIST = {
  {"test_bc_sheath_gk_1x2v_ho", test_bc_sheath_gk_1x2v_ho},
  {"test_bc_sheath_gk_2x2v_ho", test_bc_sheath_gk_2x2v_ho},
  {"test_bc_sheath_gk_3x2v_ho", test_bc_sheath_gk_3x2v_ho},
#ifdef GKYL_HAVE_CUDA
  {"test_bc_sheath_gk_1x2v_dev", test_bc_sheath_gk_1x2v_dev},
  {"test_bc_sheath_gk_2x2v_dev", test_bc_sheath_gk_2x2v_dev},
  {"test_bc_sheath_gk_3x2v_dev", test_bc_sheath_gk_3x2v_dev},
#endif
  {NULL, NULL}
};