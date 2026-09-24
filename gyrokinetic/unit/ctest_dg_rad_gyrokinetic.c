// A test for the line radiation operator
//
#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_const.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_gk_rad_vars.h>
#include <gkyl_dg_rad_gyrokinetic_drag.h>
#include <gkyl_dg_updater_rad_gyrokinetic.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_gk_maxwellian_proj_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_position_map.h>
#include <gkyl_range.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_radiation_read.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <math.h>
#include <assert.h>

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

void
mapc2p_1x(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; 
}

void
bfield_func_1x(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0]; 
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 1.0;
}

void
mapc2p_2x(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; 
}

void
bfield_func_2x(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0], y = xc[1];
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 1.0;
}

void
mapc2p_3x(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void
bfield_func_3x(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 1.0;
}

void
eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double *arr = ctx;
  double ne = arr[1];
  fout[0] = ne*1.1;
}

void
eval_cdensity(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1e19;
}

void
eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_vthsq(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double *arr=ctx;
  double te = arr[0];
  fout[0] = te*GKYL_ELEMENTARY_CHARGE/GKYL_ELECTRON_MASS;
}

void
test_1x(int poly_order, bool use_gpu, double te, int atomic_z,
  int charge_state, int num_ne[1], int ne_interval)
{
  double mass = GKYL_ELECTRON_MASS;
  double charge = -1.0*GKYL_ELEMENTARY_CHARGE;

  double vtsq_min = 0.0;
  double vth = sqrt(te*GKYL_ELEMENTARY_CHARGE/GKYL_ELECTRON_MASS);
  // Phase space and Configuration space extents and resolution
  double lower[] = {-1.0, -4*vth, 0.0};
  double upper[] = {1.0, 4*vth, 9*vth*vth*GKYL_ELECTRON_MASS};
  int cells[] = {2, 256, 128};
  int vdim = 2;

  int ndim = sizeof(cells)/sizeof(cells[0]);
  int cdim = ndim - vdim;

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double vLower[vdim], vUpper[vdim];
  int vCells[vdim];
  for (int d=0; d<vdim; d++) {
    vLower[d] = lower[cdim+d];
    vUpper[d] = upper[cdim+d];
    vCells[d] = cells[cdim+d];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vGrid;
  gkyl_rect_grid_init(&vGrid, vdim, vLower, vUpper, vCells);

  // basis functions
  struct gkyl_basis basis, confBasis, surf_vpar_basis, surf_mu_basis;
  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
    gkyl_cart_modal_serendip(&surf_mu_basis, ndim-1, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in vpar). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
    // constant mu surface
    gkyl_cart_modal_gkhybrid(&surf_mu_basis, cdim, poly_order);
  }
  // constant vpar surface
  gkyl_cart_modal_serendip(&surf_vpar_basis, ndim-1, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], 0 , 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  int vGhost[] = {0,0};
  struct gkyl_range vLocal, vLocal_ext;
  gkyl_create_grid_ranges(&vGrid, vGhost, &vLocal_ext, &vLocal);

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = {0.0, 0.0},
    .mapc2p = mapc2p_3x, // mapping of computational to physical space
    .c2p_ctx = 0,
    .bfield_func = bfield_func_3x, // magnetic field magnitude
    .bfield_ctx = 0 ,
    .position_map = pmap,
    .grid = confGrid,
    .local = confLocal,
    .local_ext = confLocal_ext,
    .global = confLocal,
    .global_ext = confLocal_ext,
    .basis = confBasis,
  };

  int geo_ghost[3] = {1};
  geometry_input.geo_grid = gkyl_gk_geometry_augment_grid(confGrid, geometry_input);
  gkyl_cart_modal_serendip(&geometry_input.geo_basis, 3, poly_order);
  gkyl_create_grid_ranges(&geometry_input.geo_grid, geo_ghost, &geometry_input.geo_global_ext, &geometry_input.geo_global);
  memcpy(&geometry_input.geo_local, &geometry_input.geo_global, sizeof(struct gkyl_range));
  memcpy(&geometry_input.geo_local_ext, &geometry_input.geo_global_ext, sizeof(struct gkyl_range));

  struct gk_geometry* gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  // deflate geometry
  struct gk_geometry *gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_input);
  gkyl_gk_geometry_release(gk_geom_3d);

  // If we are on the gpu, copy from host
  if (use_gpu) {
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_input, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, vGrid,
    local, local_ext, vLocal, vLocal_ext, use_gpu);

  // allocate drag coefficients in vparallel and mu for each collision
  // vnu = v_par*nu(v)
  // vsqnu = 2*mu*nu(v)
  // Note that through the spatial variation of B, both these drag coefficients depend on the full phase space

  int num_collisions = 1;
  int num_ne_per_coll[num_collisions];
  num_ne_per_coll[0] = num_ne[0];
  struct all_radiation_states *rad_data = gkyl_radiation_read_rad_fit_params();
  int stat = gkyl_radiation_read_get_num_densities(*rad_data, atomic_z, charge_state,
    1, 1e30, num_ne_per_coll);

  struct gkyl_gk_rad_drag *vnu, *vsqnu, *vnu_surf, *vsqnu_surf;
  vnu = gkyl_dg_calc_gk_rad_vars_drag_new(num_collisions,
    num_ne_per_coll, basis.num_basis, local_ext.volume, use_gpu);
  vsqnu = gkyl_dg_calc_gk_rad_vars_drag_new(num_collisions,
    num_ne_per_coll, basis.num_basis, local_ext.volume, use_gpu);
  vnu_surf = gkyl_dg_calc_gk_rad_vars_drag_new(num_collisions,
    num_ne_per_coll, surf_vpar_basis.num_basis, local_ext.volume, use_gpu);
  vsqnu_surf = gkyl_dg_calc_gk_rad_vars_drag_new(num_collisions,
    num_ne_per_coll, surf_mu_basis.num_basis, local_ext.volume, use_gpu);

  double a[GKYL_MAX_RAD_DENSITIES], alpha[GKYL_MAX_RAD_DENSITIES], beta[GKYL_MAX_RAD_DENSITIES],
    gamma[GKYL_MAX_RAD_DENSITIES], v0[GKYL_MAX_RAD_DENSITIES], n_elc_d[GKYL_MAX_RAD_DENSITIES];
  int status = gkyl_radiation_read_get_fit_params(*rad_data, atomic_z, charge_state,
    a, alpha, beta, gamma, v0, num_ne, n_elc_d, 1e19, 1, 1e30);
  assert(num_ne[0] == num_ne_per_coll[0]);
  if (status == 1) {
    printf("No radiation fits exist for z=%d, charge state=%d\n",atomic_z, charge_state);
    TEST_CHECK( status==0 );
  }
  double ctx[2], Lz[1];
  ne_interval = fmin(ne_interval, num_ne[0]);
  struct gkyl_array *n_elc = mkarr(use_gpu, 1, num_ne[0]);
  struct gkyl_array *n_elc_ho = mkarr(false, 1, num_ne[0]);
  memcpy(n_elc_ho->data, n_elc_d, num_ne[0]*sizeof(double));
  gkyl_array_copy(n_elc, n_elc_ho);
  gkyl_array_release(n_elc_ho);
  
  ctx[0] = te;
  ctx[1] = n_elc_d[ne_interval-1];

  gkyl_radiation_read_get_fit_lz(*rad_data, atomic_z, charge_state, log10(n_elc_d[ne_interval-1]), ctx, Lz);
  gkyl_radiation_read_release_fit_params(rad_data);
  
  struct gkyl_dg_calc_gk_rad_vars *calc_gk_rad_vars = gkyl_dg_calc_gk_rad_vars_new(&grid,
    &confBasis, &basis, charge, mass, gk_geom, gvm, use_gpu);

  gkyl_dg_calc_gk_rad_vars_nu_advance(calc_gk_rad_vars, &confLocal, &local,
    a[ne_interval-1], alpha[ne_interval-1],beta[ne_interval-1], gamma[ne_interval-1], v0[ne_interval-1],
    vnu_surf[0].data[ne_interval-1].arr, vnu[0].data[ne_interval-1].arr,
    vsqnu_surf[0].data[ne_interval-1].arr, vsqnu[0].data[ne_interval-1].arr);

  struct gkyl_array *nvnu, *nvsqnu, *nvnu_surf, *nvsqnu_surf;
  nvnu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  nvsqnu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  nvnu_surf = mkarr(use_gpu, surf_vpar_basis.num_basis, local_ext.volume);
  nvsqnu_surf = mkarr(use_gpu, surf_mu_basis.num_basis, local_ext.volume);

  // Project n, udrift, and vt^2 based on input functions
  struct gkyl_array *m0_ho = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *udrift_ho = mkarr(false, vdim*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *vtsq_ho = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *vtsq_imp = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_density, ctx);

  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_upar, 0);
  
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_vthsq, ctx);
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_ho); 
  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &confLocal, udrift_ho);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal, vtsq_ho);
  gkyl_array_copy(vtsq_imp, vtsq_ho);
  // Projection routine expects the primitive moments as a single array.
  struct gkyl_array *prim_moms_ho = mkarr(false, 3*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(prim_moms_ho, 1.0, m0_ho, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1.0, udrift_ho, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1.0, vtsq_ho, 2*confBasis.num_basis);

  // Create distribution function array
  struct gkyl_array *f = mkarr(use_gpu, basis.num_basis, local_ext.volume);

  // Maxwellian (or bi-Maxwellian) projection updater.
  struct gkyl_gk_maxwellian_proj_on_basis_inp inp_proj = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &vLocal, 
    .gk_geom = gk_geom,
    .vel_map = gvm,
    .mass = mass,
    .bimaxwellian = false, 
    .use_gpu = use_gpu,
  };
  struct gkyl_gk_maxwellian_proj_on_basis *proj_max = gkyl_gk_maxwellian_proj_on_basis_inew( &inp_proj );

  // If on GPUs, need to copy n, udrift, and vt^2 onto device
  struct gkyl_array *prim_moms, *m0, *vtsq;  
  if (use_gpu) {
    prim_moms = mkarr(use_gpu, 3*confBasis.num_basis, confLocal_ext.volume);
    m0 = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
    vtsq = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
    
    gkyl_array_copy(prim_moms, prim_moms_ho);
    gkyl_array_copy(m0, m0_ho);
    gkyl_array_copy(vtsq, vtsq_ho);

    gkyl_gk_maxwellian_proj_on_basis_advance(proj_max,
      &local, &confLocal, prim_moms, false, f);  
  } 
  else {
    vtsq = vtsq_ho;
    m0 = m0_ho;

    gkyl_gk_maxwellian_proj_on_basis_advance(proj_max,
      &local, &confLocal, prim_moms_ho, false, f);  
  }

  struct gkyl_array *vtsq_min_normalized = mkarr(use_gpu, 1, num_ne[0]);
  gkyl_array_clear(vtsq_min_normalized, vtsq_min);

  // initialize solver 
  struct gkyl_dg_updater_collisions *slvr;
  struct gkyl_dg_rad_gyrokinetic_auxfields drag_inp = { .nvnu_surf = nvnu_surf, .nvnu = nvnu, .nvsqnu_surf = nvsqnu_surf, .nvsqnu = nvsqnu};
  slvr = gkyl_dg_updater_rad_gyrokinetic_new(&grid, &confBasis, &basis, &local, &confLocal, gvm, &drag_inp, use_gpu);

  struct gkyl_array *cflrate, *rhs, *fmax;
  cflrate = mkarr(use_gpu, 1, local_ext.volume);
  rhs = mkarr(use_gpu, basis.num_basis, local_ext.volume);

  // run hyper_dg_advance
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_array_clear(nvnu_surf, 0.0);
  gkyl_array_clear(nvnu, 0.0);
  gkyl_array_clear(nvsqnu_surf, 0.0);
  gkyl_array_clear(nvsqnu, 0.0);

  // Assumed electron and ion density are the same and uniform
  gkyl_dg_calc_gk_rad_vars_nI_nu_advance(calc_gk_rad_vars, 
    &confLocal, &local, &vnu_surf[0], &vnu[0], &vsqnu_surf[0], &vsqnu[0], 
    n_elc, m0, m0, nvnu_surf, nvnu, nvsqnu_surf, nvsqnu, vtsq_min_normalized, vtsq);

  gkyl_dg_updater_rad_gyrokinetic_advance(slvr, &local, f, cflrate, rhs);

  // Take 2nd moment of rhs to find energy loss on host
  struct gkyl_dg_updater_moment *m2_calc = gkyl_dg_updater_moment_gyrokinetic_new(&grid, &confBasis, &basis,
    &confLocal, GKYL_ELECTRON_MASS, -GKYL_ELEMENTARY_CHARGE, gvm, gk_geom, NULL, GKYL_F_MOMENT_M2, false, use_gpu);

  struct gkyl_array *m2_final = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  gkyl_dg_updater_moment_gyrokinetic_advance(m2_calc, &local, &confLocal, rhs, m2_final);

  struct gkyl_array *m2_final_host;
  m2_final_host = m2_final;
  if (use_gpu) {
    m2_final_host = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
    gkyl_array_copy(m2_final_host, m2_final);
  }

  double *m00 = gkyl_array_fetch(m0_ho, 0+ghost[0]);
  double *m20 = gkyl_array_fetch(m2_final_host, 0+ghost[0]);
  
  double cell_avg_m2 = m20[0]/pow(sqrt(2.0),cdim);
  double cell_avg_m0 = m00[0]/pow(sqrt(2.0),cdim);
  // two factors of density, one for the electrons and one for the ions
  double cell_avg0 = 1.0/2.0*GKYL_ELECTRON_MASS*cell_avg_m2/(cell_avg_m0*cell_avg_m0);

  double correct = Lz[0];

  // Fit error typically >10%, so %1 should be sufficient here
  TEST_CHECK( gkyl_compare( -correct*1e30, cell_avg0*1e30, 1e-2));
  TEST_CHECK( cell_avg0<0 );
  
  // Release memory
  gkyl_dg_calc_gk_rad_vars_drag_release(vnu, num_collisions, use_gpu);
  gkyl_dg_calc_gk_rad_vars_drag_release(vnu_surf, num_collisions, use_gpu);
  gkyl_dg_calc_gk_rad_vars_drag_release(vsqnu, num_collisions, use_gpu);
  gkyl_dg_calc_gk_rad_vars_drag_release(vsqnu_surf, num_collisions, use_gpu);
  gkyl_array_release(nvnu);
  gkyl_array_release(nvnu_surf);
  gkyl_array_release(nvsqnu);
  gkyl_array_release(nvsqnu_surf);
  gkyl_array_release(n_elc);

  gkyl_dg_calc_gk_rad_vars_release(calc_gk_rad_vars);

  gkyl_array_release(m0_ho);
  gkyl_array_release(udrift_ho); 
  gkyl_array_release(vtsq_ho);
  gkyl_array_release(vtsq_imp);
  gkyl_array_release(vtsq_min_normalized);
  gkyl_array_release(prim_moms_ho);

  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);
  gkyl_gk_maxwellian_proj_on_basis_release(proj_max);

  gkyl_array_release(f);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_dg_updater_rad_gyrokinetic_release(slvr);
  gkyl_dg_updater_moment_gyrokinetic_release(m2_calc);
  gkyl_array_release(m2_final);

  if (use_gpu) {
    gkyl_array_release(m0);
    gkyl_array_release(prim_moms);   
    gkyl_array_release(m2_final_host);   
  }

  gkyl_velocity_map_release(gvm);
  gkyl_position_map_release(pmap);
  gkyl_gk_geometry_release(gk_geom);
}

void
test_2x(int poly_order, bool use_gpu, double te)
{
  double mass = GKYL_ELECTRON_MASS;
  double charge = -1.0*GKYL_ELEMENTARY_CHARGE;

  double vth = sqrt(te*GKYL_ELEMENTARY_CHARGE/GKYL_ELECTRON_MASS);
  double vtsq_min = 0.0;
  // Phase space and Configuration space extents and resolution
  double lower[] = {-2.0, -1.0, -4*vth, 0.0};
  double upper[] = {2.0, 1.0, 4*vth, 9*vth*vth*GKYL_ELECTRON_MASS};
  int cells[] = {2, 2, 256, 128};
  int vdim = 2;

  int ndim = sizeof(cells)/sizeof(cells[0]);
  int cdim = ndim - vdim;

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double vLower[vdim], vUpper[vdim];
  int vCells[vdim];
  for (int d=0; d<vdim; d++) {
    vLower[d] = lower[cdim+d];
    vUpper[d] = upper[cdim+d];
    vCells[d] = cells[cdim+d];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vGrid;
  gkyl_rect_grid_init(&vGrid, vdim, vLower, vUpper, vCells);

  // basis functions
  struct gkyl_basis basis, confBasis, surf_vpar_basis, surf_mu_basis;
  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
    gkyl_cart_modal_serendip(&surf_mu_basis, ndim-1, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in vpar). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
    // constant mu surface
    gkyl_cart_modal_gkhybrid(&surf_mu_basis, cdim, poly_order);
  }
  // constant vpar surface
  gkyl_cart_modal_serendip(&surf_vpar_basis, ndim-1, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1, 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], confGhost[1], 0 , 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  int vGhost[] = {0, 0, 0};
  struct gkyl_range vLocal, vLocal_ext;
  gkyl_create_grid_ranges(&vGrid, vGhost, &vLocal_ext, &vLocal);

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = {0.0},
    .mapc2p = mapc2p_3x, // mapping of computational to physical space
    .c2p_ctx = 0,
    .bfield_func = bfield_func_3x, // magnetic field magnitude
    .bfield_ctx = 0 ,
    .position_map = pmap,
    .grid = confGrid,
    .local = confLocal,
    .local_ext = confLocal_ext,
    .global = confLocal,
    .global_ext = confLocal_ext,
    .basis = confBasis,
  };

  int geo_ghost[3] = {1};
  geometry_input.geo_grid = gkyl_gk_geometry_augment_grid(confGrid, geometry_input);
  gkyl_cart_modal_serendip(&geometry_input.geo_basis, 3, poly_order);
  gkyl_create_grid_ranges(&geometry_input.geo_grid, geo_ghost, &geometry_input.geo_global_ext, &geometry_input.geo_global);
  memcpy(&geometry_input.geo_local, &geometry_input.geo_global, sizeof(struct gkyl_range));
  memcpy(&geometry_input.geo_local_ext, &geometry_input.geo_global_ext, sizeof(struct gkyl_range));

  struct gk_geometry* gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  // deflate geometry
  struct gk_geometry *gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_input);
  gkyl_gk_geometry_release(gk_geom_3d);

  // If we are on the gpu, copy from host
  if (use_gpu) {
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_input, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, vGrid,
    local, local_ext, vLocal, vLocal_ext, use_gpu);

  // allocate drag coefficients in vparallel and mu for each collision
  // vnu = v_par*nu(v)
  // vsqnu = 2*mu*nu(v)
  // Note that through the spatial variation of B, both these drag coefficients depend on the full phase space

  int num_collisions = 1;
  int num_ne_per_coll[1] = {1};
  int atomic_z = 3;
  int charge_state = 0;
  struct all_radiation_states *rad_data = gkyl_radiation_read_rad_fit_params();
  int stat = gkyl_radiation_read_get_num_densities(*rad_data, atomic_z, charge_state,
    1, 1e30, num_ne_per_coll);

  struct gkyl_gk_rad_drag *vnu, *vsqnu, *vnu_surf, *vsqnu_surf;
  vnu = gkyl_dg_calc_gk_rad_vars_drag_new(num_collisions,
    num_ne_per_coll, basis.num_basis, local_ext.volume, use_gpu);
  vsqnu = gkyl_dg_calc_gk_rad_vars_drag_new(num_collisions,
    num_ne_per_coll, basis.num_basis, local_ext.volume, use_gpu);
  vnu_surf = gkyl_dg_calc_gk_rad_vars_drag_new(num_collisions,
    num_ne_per_coll, surf_vpar_basis.num_basis, local_ext.volume, use_gpu);
  vsqnu_surf = gkyl_dg_calc_gk_rad_vars_drag_new(num_collisions,
    num_ne_per_coll, surf_mu_basis.num_basis, local_ext.volume, use_gpu);

  double a[1], alpha[1], beta[1], gamma[1], v0[1], n_elc_d[1];
  int num_ne[1] = {1};
  int ne_interval = 1;
  int status = gkyl_radiation_read_get_fit_params(*rad_data, atomic_z, charge_state,
    a, alpha, beta, gamma, v0, num_ne, n_elc_d, 1e19, 1, 1e30);
  assert(num_ne[0] == num_ne_per_coll[0]);

  struct gkyl_array *n_elc = mkarr(use_gpu, 1, num_ne[0]);
  struct gkyl_array *n_elc_ho = mkarr(false, 1, num_ne[0]);
  memcpy(n_elc_ho->data, n_elc_d, num_ne[0]*sizeof(double));
  gkyl_array_copy(n_elc, n_elc_ho);
  gkyl_array_release(n_elc_ho);

  if (status == 1) {
    printf("No radiation fits exist for z=%d, charge state=%d\n",atomic_z, charge_state);
    TEST_CHECK( status==0 );
  }
  double ctx[2], Lz[1];
  ctx[0] = te;
  ctx[1] = n_elc_d[ne_interval-1];

  gkyl_radiation_read_get_fit_lz(*rad_data, atomic_z, charge_state, log10(ctx[1]), ctx, Lz);
  gkyl_radiation_read_release_fit_params(rad_data);
  
  struct gkyl_dg_calc_gk_rad_vars *calc_gk_rad_vars = gkyl_dg_calc_gk_rad_vars_new(&grid,
    &confBasis, &basis, charge, mass, gk_geom, gvm, use_gpu);

  gkyl_dg_calc_gk_rad_vars_nu_advance(calc_gk_rad_vars, &confLocal, &local,
    a[0], alpha[0], beta[0], gamma[0], v0[0],
    vnu_surf[0].data[ne_interval-1].arr, vnu[0].data[ne_interval-1].arr,
    vsqnu_surf[0].data[ne_interval-1].arr, vsqnu[0].data[ne_interval-1].arr);

  struct gkyl_array *nvnu, *nvsqnu, *nvnu_surf, *nvsqnu_surf;
  nvnu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  nvsqnu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  nvnu_surf = mkarr(use_gpu, surf_vpar_basis.num_basis, local_ext.volume);
  nvsqnu_surf = mkarr(use_gpu, surf_mu_basis.num_basis, local_ext.volume);

  // Project n, udrift, and vt^2 based on input functions
  struct gkyl_array *m0_ho = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *udrift_ho = mkarr(false, vdim*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *vtsq_ho = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *vtsq_imp = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_density, ctx);

  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_upar, 0);
  
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_vthsq, ctx);

  gkyl_array_copy(vtsq_imp, vtsq_ho);
  
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_ho); 
  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &confLocal, udrift_ho);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal, vtsq_ho);
  
  // proj_maxwellian expects the primitive moments as a single array.
  struct gkyl_array *prim_moms_ho = mkarr(false, 3*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(prim_moms_ho, 1.0, m0_ho, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1.0, udrift_ho, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1.0, vtsq_ho, 2*confBasis.num_basis);

  // Create distribution function array
  struct gkyl_array *f = mkarr(use_gpu, basis.num_basis, local_ext.volume);

  // Maxwellian (or bi-Maxwellian) projection updater.
  struct gkyl_gk_maxwellian_proj_on_basis_inp inp_proj = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &vLocal, 
    .gk_geom = gk_geom,
    .vel_map = gvm,
    .mass = mass,
    .bimaxwellian = false, 
    .use_gpu = use_gpu,
  };
  struct gkyl_gk_maxwellian_proj_on_basis *proj_max = gkyl_gk_maxwellian_proj_on_basis_inew( &inp_proj );

  // If on GPUs, need to copy n, udrift, and vt^2 onto device
  struct gkyl_array *prim_moms, *m0, *vtsq;  
  if (use_gpu) {
    prim_moms = mkarr(use_gpu, 3*confBasis.num_basis, confLocal_ext.volume);
    m0 = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
    vtsq = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
    
    gkyl_array_copy(prim_moms, prim_moms_ho);
    gkyl_array_copy(m0, m0_ho);
    gkyl_array_copy(vtsq, vtsq_ho);

    gkyl_gk_maxwellian_proj_on_basis_advance(proj_max,
      &local, &confLocal, prim_moms, false, f);  
  } 
  else {
    vtsq = vtsq_ho;
    m0 = m0_ho;

    gkyl_gk_maxwellian_proj_on_basis_advance(proj_max,
      &local, &confLocal, prim_moms_ho, false, f);  
  }

  struct gkyl_array *vtsq_min_normalized = mkarr(false, 1, 1);
  gkyl_array_clear(vtsq_min_normalized, vtsq_min);

  // initialize solver 
  struct gkyl_dg_updater_collisions *slvr;
  struct gkyl_dg_rad_gyrokinetic_auxfields drag_inp = { .nvnu_surf = nvnu_surf, .nvnu = nvnu, 
    .nvsqnu_surf = nvsqnu_surf, .nvsqnu = nvsqnu};
  slvr = gkyl_dg_updater_rad_gyrokinetic_new(&grid, &confBasis, &basis, &local, &confLocal, gvm, &drag_inp, use_gpu);

  struct gkyl_array *cflrate, *rhs, *fmax;
  cflrate = mkarr(use_gpu, 1, local_ext.volume);
  rhs = mkarr(use_gpu, basis.num_basis, local_ext.volume);

  // run hyper_dg_advance
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_array_clear(nvnu_surf, 0.0);
  gkyl_array_clear(nvnu, 0.0);
  gkyl_array_clear(nvsqnu_surf, 0.0);
  gkyl_array_clear(nvsqnu, 0.0);
 
  // Assumed electron and ion density are the same and uniform
  gkyl_dg_calc_gk_rad_vars_nI_nu_advance(calc_gk_rad_vars, 
    &confLocal, &local, &vnu_surf[0], &vnu[0], &vsqnu_surf[0], &vsqnu[0], 
    n_elc, m0, m0, nvnu_surf, nvnu, nvsqnu_surf, nvsqnu, vtsq_min_normalized, vtsq);


  gkyl_dg_updater_rad_gyrokinetic_advance(slvr, &local, f, cflrate, rhs);
  
  gkyl_grid_sub_array_write(&grid, &local, 0, rhs, "ctest_dg_rad_gyrokinetic_2x_rhs.gkyl");
  gkyl_grid_sub_array_write(&grid, &local, 0, nvnu, "ctest_dg_rad_gyrokinetic_2x_nvnu.gkyl");
  gkyl_grid_sub_array_write(&grid, &local, 0, nvsqnu, "ctest_dg_rad_gyrokinetic_2x_nvsqnu.gkyl");
  gkyl_grid_sub_array_write(&grid, &local, 0, f, "ctest_dg_rad_gyrokinetic_2x_f.gkyl");
  // Take 2nd moment of rhs to find energy loss on host
  struct gkyl_dg_updater_moment *m2_calc = gkyl_dg_updater_moment_gyrokinetic_new(&grid, &confBasis, &basis,
    &confLocal, GKYL_ELECTRON_MASS, -GKYL_ELEMENTARY_CHARGE, gvm, gk_geom, NULL, GKYL_F_MOMENT_M2, false, use_gpu);
  struct gkyl_array *m2_final = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  gkyl_dg_updater_moment_gyrokinetic_advance(m2_calc, &local, &confLocal, rhs, m2_final);

  struct gkyl_array *m2_final_host;
  m2_final_host = m2_final;
  if (use_gpu) {
    m2_final_host = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
    gkyl_array_copy(m2_final_host, m2_final);
  }

  double *m00 = gkyl_array_fetch(m0_ho, 0+ghost[0]);
  double *m20 = gkyl_array_fetch(m2_final_host, 0+ghost[0]);

  double cell_avg_m2 = m20[confLocal_ext.volume]/pow(sqrt(2.0),cdim);
  double cell_avg_m0 = m00[confLocal_ext.volume]/pow(sqrt(2.0),cdim);
  // two factors of density, one for the electrons and one for the ions
  double cell_avg0 = 1.0/2.0*GKYL_ELECTRON_MASS*cell_avg_m2/(cell_avg_m0*cell_avg_m0);

  double correct = Lz[0];

  // Fit error typically >10%, so %1 should be sufficient here
  TEST_CHECK( gkyl_compare( -correct*1e30, cell_avg0*1e30, 1e-2));
  TEST_CHECK( cell_avg0<0 );

  // Release memory
  gkyl_dg_calc_gk_rad_vars_drag_release(vnu, num_collisions, use_gpu);
  gkyl_dg_calc_gk_rad_vars_drag_release(vnu_surf, num_collisions, use_gpu);
  gkyl_dg_calc_gk_rad_vars_drag_release(vsqnu, num_collisions, use_gpu);
  gkyl_dg_calc_gk_rad_vars_drag_release(vsqnu_surf, num_collisions, use_gpu);
  gkyl_array_release(nvnu);
  gkyl_array_release(nvnu_surf);
  gkyl_array_release(nvsqnu);
  gkyl_array_release(nvsqnu_surf);
  gkyl_array_release(n_elc);

  gkyl_dg_calc_gk_rad_vars_release(calc_gk_rad_vars);

  gkyl_array_release(m0_ho);
  gkyl_array_release(udrift_ho); 
  gkyl_array_release(vtsq_ho);
  gkyl_array_release(vtsq_imp);
  gkyl_array_release(vtsq_min_normalized);
  gkyl_array_release(prim_moms_ho);

  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);
  gkyl_gk_maxwellian_proj_on_basis_release(proj_max);

  gkyl_array_release(f);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_dg_updater_rad_gyrokinetic_release(slvr);
  gkyl_dg_updater_moment_gyrokinetic_release(m2_calc);
  gkyl_array_release(m2_final);

  if (use_gpu) {
    gkyl_array_release(m0);
    gkyl_array_release(vtsq);
    gkyl_array_release(prim_moms);   
    gkyl_array_release(m2_final_host);   
  }

  gkyl_velocity_map_release(gvm);
  gkyl_position_map_release(pmap);
  gkyl_gk_geometry_release(gk_geom);
}

static int num_ne[1] = {1};
static int num_ne2[1] = {20};
void test_1x2v_p1_10eV() { test_1x(1, false, 10.0, 3, 0, num_ne, 1); }
void test_1x2v_p1_30eV() { test_1x(1, false, 30.0, 3, 0, num_ne, 1); }
void test_1x2v_p1_H() { test_1x(1, false, 30.0, 1, 0, num_ne, 1); }
void test_1x2v_p1_100eV() { test_1x(1, false, 100.0, 3, 0, num_ne, 1); }
void test_1x2v_p1_500eV() { test_1x(1, false, 500.0, 3, 0, num_ne, 1); }
void test_1x2v_p1_1000eV() { test_1x(1, false, 1000.0, 3, 0, num_ne, 1); }
void test_1x2v_p1_5000eV() { test_1x(1, false, 5000.0, 3, 0, num_ne, 1); }
void test_1x2v_p1_10000eV() { test_1x(1, false, 10000.0, 3, 0, num_ne, 1); }
void test_2x2v_p1() { test_2x(1, false, 30.0); }

void test_1x2v_p1_Li1_lowNe() { test_1x(1, false, 30.0, 3, 1, num_ne2, 1); }
void test_1x2v_p1_Li1_midNe() { test_1x(1, false, 30.0, 3, 1, num_ne2, 6); }
void test_1x2v_p1_Li1_highNe() { test_1x(1, false, 30.0, 3, 1, num_ne2, 13); }

#ifdef GKYL_HAVE_CUDA

void test_1x2v_p1_gpu() { test_1x(1, true, 30.0, 3, 0, num_ne, 1); }
void test_1x2v_p1_L1_midNe_gpu() {test_1x(1, true, 30.0, 3, 1, num_ne2, 6); }

#endif

TEST_LIST = {
  { "test_1x2v_p1_Li0_30eV", test_1x2v_p1_30eV },
  { "test_1x2v_p1_Li0_5000eV", test_1x2v_p1_5000eV },
  { "test_1x2v_p1_H", test_1x2v_p1_H },
  { "test_1x2v_p1_Li1_lowNe", test_1x2v_p1_Li1_lowNe },
  { "test_1x2v_p1_Li1_midNe", test_1x2v_p1_Li1_midNe },
  { "test_1x2v_p1_Li1_highNe", test_1x2v_p1_Li1_highNe },
  { "test_2x2v_p1", test_2x2v_p1 },

#ifdef GKYL_HAVE_CUDA
  { "test_1x2v_p1_gpu", test_1x2v_p1_gpu },
  { "test_1x2v_p1_L1_midNe_gpu", test_1x2v_p1_L1_midNe_gpu},

#endif
  { NULL, NULL },
};
