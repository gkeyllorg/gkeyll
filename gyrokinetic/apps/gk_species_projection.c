#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_util.h>

static void
load_projection_moment_from_file(struct gkyl_gyrokinetic_app *app, struct gkyl_array *arr,
  const struct gkyl_gyrokinetic_ic_import *inp)
{
  if (inp == NULL || inp->type == GKYL_IC_IMPORT_NONE) {
    return;
  }

  struct gkyl_array *arr_host = app->use_gpu ? mkarr(false, arr->ncomp, arr->size)
                                             : gkyl_array_acquire(arr);

  struct gkyl_app_restart_status rstat = { .io_status = GKYL_ARRAY_RIO_FOPEN_FAILED };
  rstat.io_status = gkyl_comm_array_read(app->comm, &app->grid, &app->local, arr_host, inp->file_name);
  assert(rstat.io_status == GKYL_ARRAY_RIO_SUCCESS);

  gkyl_array_copy(arr, arr_host);
  gkyl_array_release(arr_host);
}

static void
proj_on_basis_c2p_phase_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  int cdim = c2p_ctx->cdim; // Assumes update range is a phase range.
  gkyl_position_map_eval_mc2nu(c2p_ctx->pos_map, xcomp, xphys);
  gkyl_velocity_map_eval_c2p(c2p_ctx->vel_map, &xcomp[cdim], &xphys[cdim]);
}

static void
proj_on_basis_c2p_position_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  gkyl_position_map_eval_mc2nu(c2p_ctx->pos_map, xcomp, xphys);
}

struct func_gaussian_ctx {
  int cdim; // Configuration space dimension.
  bool is_dir_periodic[GKYL_MAX_CDIM]; // Periodicity in configuration space.
  double box_size[GKYL_MAX_CDIM]; // Size of the box in each direction
  double gaussian_mean[GKYL_MAX_CDIM]; // Center in configuration space.
  double gaussian_std_dev[GKYL_MAX_CDIM]; // Sigma in configuration space, function is constant if sigma is 0.
  double f_floor; // Floor value of the distribution.
};
static void 
func_gaussian(double t, const double* xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct func_gaussian_ctx *inp = ctx;
  double envelope = 1.0;
  for (int dir = 0; dir < inp->cdim; ++dir) {
    double dx = xn[dir] - inp->gaussian_mean[dir];
    double sigma = inp->gaussian_std_dev[dir];
    double L = inp->box_size[dir];
    if (inp->is_dir_periodic[dir]) { 
      // Periodic wrapping
      dx = fmod(dx + L/2.0, L);
      if (dx < 0) dx += L;
      dx -= L/2.0;
    }
    if (sigma > 0.0)
      envelope *= exp(-dx*dx/(2.0*sigma*sigma));
  }
  fout[0] = envelope + inp->f_floor;
}

static void
gk_species_projection_multiply_jacobians(gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gkyl_array *f)
{
  // Multiply by the gyrocenter coord jacobian (bmag).
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &s->basis, f,
      app->gk_geom->geo_corn.bmag, f, &app->local, &s->local);
  // Multiply by the velocity-space jacobian.
  gkyl_array_scale_by_cell(f, s->vel_map->jacobvel);
}

static void
gk_species_projection_calc_proj_func(gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  // The projection runs on the host or on the device (directly into f)
  // depending on how proj_func was created.
  gkyl_proj_on_basis_advance(proj->proj_func, tm, &s->local, f);
  gk_species_projection_multiply_jacobians(app, s, f);
}

static void
gk_species_projection_calc_proj_func_host(gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  // Project on the host and copy to the device.
  gkyl_proj_on_basis_advance(proj->proj_func, tm, &s->local, proj->proj_host);
  gkyl_array_copy(f, proj->proj_host);
  gk_species_projection_multiply_jacobians(app, s, f);
}

static void
project_moment_if_needed(bool from_file, struct gkyl_proj_on_basis *proj_op,
  double tm, const struct gkyl_range *conf_range, struct gkyl_array *arr, double scale_fac)
{
  if (from_file) {
    return;
  }

  gkyl_proj_on_basis_advance(proj_op, tm, conf_range, arr);
  if (scale_fac != 1.0) {
    gkyl_array_scale(arr, scale_fac);
  }
}

static void
init_moment_from_import_or_proj(struct gkyl_gyrokinetic_app *app, struct gk_proj *proj,
  bool from_file, struct gkyl_array *arr, const struct gkyl_gyrokinetic_ic_import *import_inp,
  evalf_t eval, void *ctx, double scale_fac, struct gkyl_proj_on_basis **proj_on_basis)
{
  if (from_file) {
    load_projection_moment_from_file(app, arr, import_inp);
    if (scale_fac != 1.0) {
      gkyl_array_scale(arr, scale_fac);
    }
    return;
  }

  *proj_on_basis = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp) {
      .grid = &app->grid,
      .basis = &app->basis,
      .qtype = GKYL_GAUSS_QUAD,
      .num_quad = app->basis.poly_order+1,
      .num_ret_vals = 1,
      .eval = eval,
      .ctx = ctx,
      .c2p_func = proj_on_basis_c2p_position_func,
      .c2p_func_ctx = &proj->proj_on_basis_c2p_ctx,
    }
  );
}

static void
gk_species_projection_calc_max_prim(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  if (proj->maxwellian_moms_from_file) {
    gkyl_array_copy(proj->prim_moms, proj->prim_moms_host);
    gkyl_gk_maxwellian_proj_on_basis_advance(proj->proj_max,
      &s->local, &app->local, proj->prim_moms, false, f);
    return;
  }

  project_moment_if_needed(proj->dens_from_file, proj->proj_dens,
    tm, &app->local, proj->dens, 1.0);
  project_moment_if_needed(proj->upar_from_file, proj->proj_upar,
    tm, &app->local, proj->upar, 1.0);
  project_moment_if_needed(proj->temp_from_file, proj->proj_temp,
    tm, &app->local, proj->vtsq, 1.0/s->info.mass);

  // proj_maxwellian expects the primitive moments as a single array.
  gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->dens, 0*app->basis.num_basis);
  gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->upar, 1*app->basis.num_basis);
  gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->vtsq, 2*app->basis.num_basis);

  // Copy the contents into the array we will use (potentially on GPUs).
  gkyl_array_copy(proj->prim_moms, proj->prim_moms_host);
  gkyl_gk_maxwellian_proj_on_basis_advance(proj->proj_max,
    &s->local, &app->local, proj->prim_moms, false, f);
}

static void
gk_species_projection_calc_bimax(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  if (proj->bimaxwellian_moms_from_file) {
    gkyl_array_copy(proj->prim_moms, proj->prim_moms_host);
    gkyl_gk_maxwellian_proj_on_basis_advance(proj->proj_max,
      &s->local, &app->local, proj->prim_moms, false, f);
    return;
  }

  project_moment_if_needed(proj->dens_from_file, proj->proj_dens,
    tm, &app->local, proj->dens, 1.0);
  project_moment_if_needed(proj->upar_from_file, proj->proj_upar,
    tm, &app->local, proj->upar, 1.0);
  project_moment_if_needed(proj->temppar_from_file, proj->proj_temppar,
    tm, &app->local, proj->vtsqpar, 1.0/s->info.mass);
  project_moment_if_needed(proj->tempperp_from_file, proj->proj_tempperp,
    tm, &app->local, proj->vtsqperp, 1.0/s->info.mass);

  // proj_bimaxwellian expects the primitive moments as a single array.
  gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->dens, 0*app->basis.num_basis);
  gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->upar, 1*app->basis.num_basis);
  gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->vtsqpar , 2*app->basis.num_basis);
  gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->vtsqperp , 3*app->basis.num_basis);  

  // Copy the contents into the array we will use (potentially on GPUs).
  gkyl_array_copy(proj->prim_moms, proj->prim_moms_host);
  gkyl_gk_maxwellian_proj_on_basis_advance(proj->proj_max,
    &s->local, &app->local, proj->prim_moms, false, f);
}

static void
gk_species_projection_calc_max_gauss(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  bool correct_mom_setting = s->lte.correct_all_moms;
  s->lte.correct_all_moms = false; // Turn off moment correction for the max gauss projection.
  gk_species_lte_from_moms(app, s, &s->lte, proj->prim_moms);
  gkyl_array_copy(f, s->lte.f_lte);
  s->lte.correct_all_moms = correct_mom_setting; // Reset to original setting.
}

static void
gk_species_projection_calc_none(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm) {}

static void
gk_species_projection_correct_all_moms(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  struct gkyl_gk_maxwellian_correct_status status_corr;
  status_corr = gkyl_gk_maxwellian_correct_all_moments(proj->corr_max, 
    f, proj->prim_moms, &s->local, &app->local);
}

static void
gk_species_projection_correct_all_moms_none(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm) {}

static void
init_maxwellian_bimaxwellian(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj)
  {
  bool maxwellian_moms_from_file = proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM &&
    inp.maxwellian_moms_import.type != GKYL_IC_IMPORT_NONE;
  proj->maxwellian_moms_from_file = maxwellian_moms_from_file;
  bool bimaxwellian_moms_from_file = proj->proj_id == GKYL_PROJ_BIMAXWELLIAN &&
    inp.bimaxwellian_moms_import.type != GKYL_IC_IMPORT_NONE;
  proj->bimaxwellian_moms_from_file = bimaxwellian_moms_from_file;

  // Determine whether we import the primitive moments from files or compute them from functions.
  proj->dens_from_file = inp.density_import.type != GKYL_IC_IMPORT_NONE;
  proj->upar_from_file = inp.upar_import.type != GKYL_IC_IMPORT_NONE;
  proj->temp_from_file = inp.temp_import.type != GKYL_IC_IMPORT_NONE;
  proj->temppar_from_file = inp.temppar_import.type != GKYL_IC_IMPORT_NONE;
  proj->tempperp_from_file = inp.tempperp_import.type != GKYL_IC_IMPORT_NONE;

  if (maxwellian_moms_from_file) {
    proj->dens_from_file = true;
    proj->upar_from_file = true;
    proj->temp_from_file = true;
  }
  else if (bimaxwellian_moms_from_file) {
    proj->dens_from_file = true;
    proj->upar_from_file = true;
    proj->temppar_from_file = true;
    proj->tempperp_from_file = true;
  }

  // prim_moms_host is a combined array of all the projected moments
  int prim_moms_ncomp = proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM ? 3 : 4;
  proj->prim_moms_host = mkarr(false, prim_moms_ncomp*app->basis.num_basis, app->local_ext.volume);
  proj->prim_moms = mkarr(app->use_gpu, proj->prim_moms_host->ncomp, proj->prim_moms_host->size);

  proj->dens = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  proj->upar = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) {
    proj->vtsq = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  }
  else if (proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
    proj->vtsqpar = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    proj->vtsqperp = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  }

  bool bimaxwellian = false;
  if (maxwellian_moms_from_file) {
    load_projection_moment_from_file(app, proj->prim_moms_host, &inp.maxwellian_moms_import);
  }
  else if (bimaxwellian_moms_from_file) {
    bimaxwellian = true;
    load_projection_moment_from_file(app, proj->prim_moms_host, &inp.bimaxwellian_moms_import);
  }
  else {
    init_moment_from_import_or_proj(app, proj, proj->dens_from_file, proj->dens, 
      &inp.density_import, inp.density, inp.ctx_density, 1.0, &proj->proj_dens);
    init_moment_from_import_or_proj(app, proj, proj->upar_from_file, proj->upar,
      &inp.upar_import, inp.upar, inp.ctx_upar, 1.0, &proj->proj_upar);
    if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) {
      init_moment_from_import_or_proj(app, proj, proj->temp_from_file, proj->vtsq,
        &inp.temp_import, inp.temp, inp.ctx_temp, 1.0/s->info.mass, &proj->proj_temp);
    }
    else {
      bimaxwellian = true;
      init_moment_from_import_or_proj(app, proj,
        proj->temppar_from_file, proj->vtsqpar, &inp.temppar_import,
        inp.temppar, inp.ctx_temppar, 1.0/s->info.mass, &proj->proj_temppar);
      init_moment_from_import_or_proj(app, proj,
        proj->tempperp_from_file, proj->vtsqperp, &inp.tempperp_import,
        inp.tempperp, inp.ctx_tempperp, 1.0/s->info.mass, &proj->proj_tempperp);
    }
  }

  // Maxwellian (or bi-Maxwellian) projection updater.
  struct gkyl_gk_maxwellian_proj_on_basis_inp inp_proj = {
    .phase_grid = &s->grid,
    .conf_basis = &app->basis,
    .phase_basis = &s->basis,
    .conf_range =  &app->local,
    .conf_range_ext = &app->local_ext,
    .vel_range = &s->local_vel, 
    .gk_geom = app->gk_geom,
    .vel_map = s->vel_map,
    .quad_type = inp.quad_type,
    .mass = s->info.mass,
    .bimaxwellian = bimaxwellian, 
    .divide_jacobgeo = false, // final Jacobian multiplication will be handled in advance
    .use_gpu = app->use_gpu,
  };
  proj->proj_max = gkyl_gk_maxwellian_proj_on_basis_inew( &inp_proj );

  proj->correct_all_moms = false; 
  if (inp.correct_all_moms) {
    proj->correct_all_moms = true;

    int max_iter = inp.max_iter > 0 ? inp.max_iter : 100;
    double iter_eps = inp.iter_eps > 0 ? inp.iter_eps  : 1e-12;
    bool use_last_converged = inp.use_last_converged; 

    // Maxwellian correction updater
    struct gkyl_gk_maxwellian_correct_inp inp_corr = {
      .phase_grid = &s->grid,
      .conf_basis = &app->basis,
      .phase_basis = &s->basis,
      .conf_range =  &app->local,
      .conf_range_ext = &app->local_ext,
      .vel_range = &s->local_vel, 
      .gk_geom = app->gk_geom,
      .vel_map = s->vel_map,
      .mass = s->info.mass,
      .bimaxwellian = bimaxwellian, 
      .divide_jacobgeo = false, // final Jacobian multiplication will be handled in advance
      .max_iter = max_iter,
      .eps = iter_eps,
      .use_last_converged = use_last_converged, 
      .use_gpu = app->use_gpu,
    };
    proj->corr_max = gkyl_gk_maxwellian_correct_inew( &inp_corr );
  }

}

static void
init_maxwellian_gaussian(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj)
{
  // Fill the box_size attribute of the projection (used for periodicity).
  struct func_gaussian_ctx fg_ctx;
  fg_ctx.cdim = app->cdim;
  fg_ctx.f_floor = inp.f_floor;
  for (int dir = 0; dir < app->cdim; ++dir) {
    fg_ctx.gaussian_mean[dir] = inp.gaussian_mean[dir];
    fg_ctx.gaussian_std_dev[dir] = inp.gaussian_std_dev[dir];
    fg_ctx.box_size[dir] = app->grid.upper[dir] - app->grid.lower[dir];
  }
  // By default, set all directions to non-periodic.
  for (int dir = 0; dir < GKYL_MAX_CDIM; ++dir)
    fg_ctx.is_dir_periodic[dir] = false;

  // Set periodicity for last dim if we are in IWL, and all other directions defined by the user.
  // First recover the BCs of the last config. space dimension for this species from the user input.
  struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(s->info.bcs, 2*app->cdim, app->cdim-1, GKYL_LOWER_EDGE);
  struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(s->info.bcs, 2*app->cdim, app->cdim-1, GKYL_UPPER_EDGE);
  if (bc_lo != 0 && bc_up != 0) {
    // Apply periodicity condition if both edges are IWL.
    fg_ctx.is_dir_periodic[app->cdim-1] = app->gk_geom->has_LCFS ||
      (bc_lo->type == GKYL_BC_GK_SPECIES_TWISTSHIFT && bc_up->type == GKYL_BC_GK_SPECIES_TWISTSHIFT);
  }

  // Set periodicity also according to the global app settings.
  for (int i=0; i < app->num_periodic_dir; ++i)
    fg_ctx.is_dir_periodic[app->periodic_dirs[i]] = true;

  struct gkyl_array *shape_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  struct gkyl_proj_on_basis *proj_gaussian = gkyl_proj_on_basis_new(&app->grid, &app->basis, app->poly_order + 1, 1, func_gaussian, &fg_ctx);
  proj->gaussian_profile = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  gkyl_proj_on_basis_advance(proj_gaussian, 0, &app->local, shape_ho);
  gkyl_array_copy(proj->gaussian_profile, shape_ho);

  gkyl_proj_on_basis_release(proj_gaussian);
  gkyl_array_release(shape_ho);

  // Build the integrant Jacobian * s(x), to normalize the shape function, and integrate it.
  struct gkyl_array *integrant = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  double *integral = app->use_gpu? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
  double *red_integral = app->use_gpu? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
  struct gkyl_array_integrate *int_op = gkyl_array_integrate_new(&app->grid, &app->basis, 1, 
    GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu);
  double red_integral_ho[1];

  gkyl_dg_mul_op_range(&app->basis, 0, integrant, 0, app->gk_geom->geo_int.jacobgeo, 0, proj->gaussian_profile, &app->local);
  gkyl_array_integrate_advance(int_op, integrant, 1.0, NULL, &app->local, NULL, integral);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, integral, red_integral);
  
  if (app->use_gpu) {
    gkyl_cu_memcpy(red_integral_ho, red_integral, sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_free(integral);
    gkyl_cu_free(red_integral);
  } else {
    memcpy(red_integral_ho, red_integral, sizeof(double));
    gkyl_free(integral);
    gkyl_free(red_integral);
  }
  gkyl_array_release(integrant);
  gkyl_array_integrate_release(int_op);
  // Scale the shape configuration function
  gkyl_array_scale(proj->gaussian_profile, 1.0/red_integral_ho[0]);
  // We can now build the moments of the projection.
  proj->prim_moms = mkarr(app->use_gpu, 4*app->basis.num_basis, app->local_ext.volume);      
  // Density
  gkyl_array_set_offset(proj->prim_moms, inp.total_num_particles + inp.f_floor, proj->gaussian_profile, 0*app->basis.num_basis);
  // Parallel velocity
  gkyl_array_set_offset(proj->prim_moms, 0.0, proj->gaussian_profile, 1*app->basis.num_basis);
  // Temperature
  assert(inp.temp_max > 0);
  double vdim_phys = s->info.vdim == 1? 1.0 : 3.0;
  double temp = inp.total_num_particles == 0 ? inp.temp_max/2.0 : 2./vdim_phys * inp.total_kin_energy/inp.total_num_particles;
  temp = temp > inp.temp_max ? inp.temp_max : temp; // saturate to max temperature.
  double dg_norm = pow(sqrt(2.0), app->cdim);
  gkyl_array_shiftc(proj->prim_moms, dg_norm * temp/s->info.mass, 2*app->basis.num_basis);
  // Moment correction
  proj->correct_all_moms = inp.correct_all_moms;
}

void 
gk_species_projection_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj)
{
  proj->proj_id = inp.proj_id;
  // Context for c2p function passed to proj_on_basis.
  proj->proj_on_basis_c2p_ctx.cdim = app->cdim;
  proj->proj_on_basis_c2p_ctx.vdim = s->local_vel.ndim;
  proj->proj_on_basis_c2p_ctx.vel_map = s->vel_map;
  proj->proj_on_basis_c2p_ctx.pos_map = app->position_map;
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    proj->proj_host = 0;
    proj->proj_on_basis_c2p_ctx_dev = 0;

    bool proj_on_gpu = false;
#ifdef GKYL_HAVE_CUDA
    proj_on_gpu = app->use_gpu && (inp.func_on_dev != 0);
    if (proj_on_gpu) {
      // Project on the device with the user's device function. The c2p
      // context must live on the device and hold the maps' device objects.
      struct gk_proj_on_basis_c2p_func_ctx c2p_ctx_dev_ho = {
        .cdim = app->cdim,
        .vdim = s->local_vel.ndim,
        .vel_map = s->vel_map->on_dev,
        .pos_map = app->position_map->on_dev,
      };
      proj->proj_on_basis_c2p_ctx_dev = gkyl_cu_malloc(sizeof(struct gk_proj_on_basis_c2p_func_ctx));
      gkyl_cu_memcpy(proj->proj_on_basis_c2p_ctx_dev, &c2p_ctx_dev_ho,
        sizeof(struct gk_proj_on_basis_c2p_func_ctx), GKYL_CU_MEMCPY_H2D);

      proj->proj_func = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
          .grid = &s->grid,
          .basis = &s->basis,
          .qtype = GKYL_GAUSS_QUAD,
          .num_quad = app->basis.poly_order+1,
          .num_ret_vals = 1,
          .eval = inp.func_on_dev,
          .ctx = inp.ctx_func_on_dev,
          .c2p_func = gk_species_projection_c2p_phase_func_cu_dev_ptr(),
          .c2p_func_ctx = proj->proj_on_basis_c2p_ctx_dev,
          .use_gpu = true,
        }
      );
    }
#endif
    if (!proj_on_gpu) {
      proj->proj_func = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
          .grid = &s->grid,
          .basis = &s->basis,
          .qtype = GKYL_GAUSS_QUAD,
          .num_quad = app->basis.poly_order+1,
          .num_ret_vals = 1,
          .eval = inp.func,
          .ctx = inp.ctx_func,
          .c2p_func = proj_on_basis_c2p_phase_func,
          .c2p_func_ctx = &proj->proj_on_basis_c2p_ctx,
        }
      );
      if (app->use_gpu)
        proj->proj_host = mkarr(false, s->basis.num_basis, s->local_ext.volume);
    }

    proj->projection_calc = (app->use_gpu && !proj_on_gpu) ?
      gk_species_projection_calc_proj_func_host : gk_species_projection_calc_proj_func;
    proj->moms_correct = gk_species_projection_correct_all_moms_none;
  }
  else {
    if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM || proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
      init_maxwellian_bimaxwellian(app, s, inp, proj);
      proj->projection_calc = proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM ?
        gk_species_projection_calc_max_prim : gk_species_projection_calc_bimax;
    } else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_GAUSSIAN) {
      init_maxwellian_gaussian(app, s, inp, proj);
      proj->projection_calc = gk_species_projection_calc_max_gauss;
    }
    proj->moms_correct = proj->correct_all_moms ? 
      gk_species_projection_correct_all_moms : gk_species_projection_correct_all_moms_none;
  }
}

void
gk_species_projection_calc(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  proj->projection_calc(app, s, proj, f, tm);
  proj->moms_correct(app, s, proj, f, tm);  
  // Multiply by the configuration space jacobian.
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &s->basis, f, 
    app->gk_geom->geo_int.jacobgeo, f, &app->local, &s->local);      
}

void
gk_species_projection_release(const struct gkyl_gyrokinetic_app *app, const struct gk_proj *proj)
{
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    gkyl_proj_on_basis_release(proj->proj_func);
    if (proj->proj_host) {
      gkyl_array_release(proj->proj_host);
    }
#ifdef GKYL_HAVE_CUDA
    if (proj->proj_on_basis_c2p_ctx_dev) {
      gkyl_cu_free(proj->proj_on_basis_c2p_ctx_dev);
    }
#endif
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM || proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) { 
    gkyl_array_release(proj->dens);
    gkyl_array_release(proj->upar); 
    gkyl_array_release(proj->prim_moms_host);
    gkyl_array_release(proj->prim_moms);
    if (!proj->dens_from_file) {
      gkyl_proj_on_basis_release(proj->proj_dens);
    }
    if (!proj->upar_from_file) {
      gkyl_proj_on_basis_release(proj->proj_upar);
    }
    gkyl_gk_maxwellian_proj_on_basis_release(proj->proj_max);

    if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) {
      gkyl_array_release(proj->vtsq);
      if (!proj->temp_from_file) {
        gkyl_proj_on_basis_release(proj->proj_temp);
      }
    }
    else if (proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
      gkyl_array_release(proj->vtsqpar);
      gkyl_array_release(proj->vtsqperp);
      if (!proj->temppar_from_file) {
        gkyl_proj_on_basis_release(proj->proj_temppar);
      }
      if (!proj->tempperp_from_file) {
        gkyl_proj_on_basis_release(proj->proj_tempperp);
      }
    }

    if (proj->correct_all_moms) {
        gkyl_gk_maxwellian_correct_release(proj->corr_max);    
    }
  } 
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_GAUSSIAN) {
    gkyl_array_release(proj->gaussian_profile);
    gkyl_array_release(proj->prim_moms);
  }
}
