#include <gkyl_alloc.h>
#include <gkyl_flash_priv.h>
#include <gkyl_flash.h>
#include <gkyl_null_comm.h>
#include <gkyl_util.h>

#include <assert.h>
#include <mpack.h>
#include <stdarg.h>

gkyl_flash_app*
gkyl_flash_app_new(struct gkyl_flash_inp *inp)
{
  gkyl_flash_app *app = gkyl_malloc(sizeof(struct gkyl_flash_app));
  bool use_gpu = 0;
  app->use_gpu = use_gpu;

  strcpy(app->name, inp->name);

  int cdim = app->cdim = inp->cdim;
  int vdim = app->vdim = inp->vdim;
  int pdim = cdim + vdim;
  int poly_order = app->poly_order = 1;

  // Skip initializing a fluid species? For use with an external fluid solver.
  app->skip_fluid_species = inp->skip_fluid_species;

  double cfl_frac = inp->cfl_frac == 0 ? 1.0 : inp->cfl_frac;
  app->cfl = cfl_frac;

  // Set periodic boundary condition directions
  app->num_periodic_dir = inp->num_periodic_dir;
  for (int d=0; d<cdim; ++d) 
    app->periodic_dirs[d] = inp->periodic_dirs[d];

  // Basis functions: always p=1 Tensor in configuration space, hybrid in phase space.
  gkyl_cart_modal_serendip(&app->confBasis, cdim, poly_order);
  gkyl_cart_modal_hybrid(&app->basis, cdim, vdim);

  // Initialize stats object.
  app->stat = (struct gkyl_flash_stat) { };

  // Initialize configuration space grids.
  // Configuration space grid for distribution function is half
  // the resolution of the grid for finite volume fluid quantities.
  int cells_kinetic[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    assert(inp->cells_fluid[d]%2 == 0);
    cells_kinetic[d] = (int)(inp->cells_fluid[d]/2);
  }
  gkyl_rect_grid_init(&app->grid_fluid, cdim, inp->lower_conf, inp->upper_conf, inp->cells_fluid);
  gkyl_rect_grid_init(&app->grid_kinetic, cdim, inp->lower_conf, inp->upper_conf, cells_kinetic);

  // Set number of ghost cells for kinetic and fluid quantities.
  // Note that wave_prop requires 2 ghost cells.
  for (int d=0; d<GKYL_MAX_CDIM; ++d) {
    app->nghost_fluid[d] = inp->nghost_fluid[d] != 0 ? inp->nghost_fluid[d] : 2;
    app->nghost_kinetic[d] = inp->nghost_kinetic[d] != 0 ? inp->nghost_kinetic[d] : 1;
  }
  gkyl_create_grid_ranges(&app->grid_fluid, app->nghost_fluid, &app->global_fluid_ext, &app->global_fluid);
  gkyl_create_grid_ranges(&app->grid_kinetic, app->nghost_kinetic, &app->global_kinetic_ext, &app->global_kinetic);

  // HARDCODE: No MPI
  int cuts[3] = { 1, 1, 1 };
  app->decomp_fluid = gkyl_rect_decomp_new_from_cuts(cdim, cuts, &app->global_fluid);
  app->decomp_kinetic = gkyl_rect_decomp_new_from_cuts(cdim, cuts, &app->global_kinetic);

  app->comm_fluid = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = app->decomp_fluid,
      .use_gpu = app->use_gpu
    }
  );
  app->comm_kinetic = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = app->decomp_kinetic,
      .use_gpu = app->use_gpu
    }
  );

  memcpy(&app->local_fluid, &app->global_fluid, sizeof(struct gkyl_range));
  memcpy(&app->local_fluid_ext, &app->global_fluid_ext, sizeof(struct gkyl_range));
  memcpy(&app->local_kinetic, &app->global_kinetic, sizeof(struct gkyl_range));
  memcpy(&app->local_kinetic_ext, &app->global_kinetic_ext, sizeof(struct gkyl_range));

  // Allocate memory for species
  app->num_species = inp->num_species;
  app->species = gkyl_malloc(sizeof(struct flash_species[app->num_species]));
  app->species->info = *inp->species;
  app->fluid_species = gkyl_malloc(sizeof(struct flash_fluid_species[app->num_species]));
  app->fluid_species->info = *inp->species;

  // Initialize field object
  app->field = flash_field_new(inp, app);

  // Initialize species objects
  for (int i=0; i<app->num_species; ++i) {
    flash_species_init(inp, app, &app->species[i]);

    // Initialize fluid species unless specified
    if (!app->skip_fluid_species)
      flash_fluid_species_init(inp, app, &app->fluid_species[i]);
  }

  // Updater to interpolate between finite volume and DG quantities
  app->interp_vars = gkyl_flash_interp_vars_new(&app->confBasis, &app->grid_kinetic,
    &app->grid_fluid, &app->local_kinetic_ext, &app->local_fluid_ext);

  return app;
}

double
gkyl_flash_app_compute_dt_max(gkyl_flash_app *app)
{

  flash_field_calc_bvar(app, app->field, app->field->em);

  double dt_max = DBL_MAX;
  for (int i=0; i<app->num_species; ++i) {
    // Compute PKPM moments, [rho, p_par, p_perp].
    flash_species_calc_pkpm_moms(app, &app->species[i], app->species[i].f);

    // Compute coupling moments, self-consistent pressure
    // force (div(p_par b_hat)), and primitive variables
    flash_species_calc_pkpm_vars(app, &app->species[i], app->species[i].f);

    // Compute remaining PKPM variables: acceleration variables in kinetic equation and 
    // source distribution functions for Laguerre couplings.
    flash_species_calc_pkpm_update_vars(app, &app->species[i], app->species[i].f);

    dt_max = fmin(dt_max, flash_species_compute_dt(app, &app->species[i]));
  }
  return dt_max;
}

void
gkyl_flash_app_apply_ic_species(gkyl_flash_app* app, int sidx, double t0)
{
  struct timespec wst = gkyl_wall_clock();
  flash_species_apply_ic(app, &app->species[sidx], t0);
  app->stat.init_species_tm += gkyl_time_diff_now_sec(wst);

  if (!app->skip_fluid_species){
    wst = gkyl_wall_clock();
    flash_fluid_species_apply_ic(app, &app->fluid_species[sidx], t0);
    app->stat.init_fluid_species_tm += gkyl_time_diff_now_sec(wst);
  }
}

void
gkyl_flash_app_apply_ic_field(gkyl_flash_app* app, double t0)
{
  struct timespec wst = gkyl_wall_clock();
  flash_field_apply_ic(app, app->field, t0);
  app->stat.init_field_tm += gkyl_time_diff_now_sec(wst);
}

void
gkyl_flash_app_apply_ic(gkyl_flash_app *app, double t0) 
{
  struct timespec wst = gkyl_wall_clock();

  // Initialize species distributions
  for (int i=0; i<app->num_species; ++i) {
    gkyl_flash_app_apply_ic_species(app, i, t0);
  }

  // Initialize electromagnetic field variables
  gkyl_flash_app_apply_ic_field(app, t0);

  // Apply boundary conditions
  for (int i=0; i<app->num_species; ++i) {
    flash_species_apply_bc(app, &app->species[i], app->species[i].f);

    if (!app->skip_fluid_species) {
      flash_fluid_species_apply_bc(app, &app->fluid_species[i], app->fluid_species[i].fluid);
    }
  }
  app->stat.init_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_update_status
gkyl_flash_update(gkyl_flash_app *app, double dt) 
{
  app->stat.nup += 1;

  struct timespec wst = gkyl_wall_clock();
  struct gkyl_update_status status;
  if (!app->skip_fluid_species) {
    status = gkyl_flash_update_ssp_rk3(app, dt);
  } else {
    gkyl_exit("Calling gkyl_flash_update without a fluid species! Try with skip_fluid_species = false in app initialization.");
  }
  app->tcurr += status.dt_actual;

  app->stat.total_tm += gkyl_time_diff_now_sec(wst);

  return status;
}

// Initialize msgpack metadata for I/O.
// Returned gkyl_msgpack_data must be freed with flash_io_meta_release
static struct gkyl_msgpack_data*
flash_io_meta_new(struct flash_io_meta meta) {
  struct gkyl_msgpack_data *mt = gkyl_malloc(sizeof(*mt));

  mt->meta_sz = 0;
  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mt->meta, &mt->meta_sz);

  mpack_build_map(&writer);
  
  mpack_write_cstr(&writer, "time");
  mpack_write_double(&writer, meta.time);

  mpack_write_cstr(&writer, "frame");
  mpack_write_i64(&writer, meta.frame);

  // If finite volume data, basis_type is not set.
  if (meta.basis_type) {
    mpack_write_cstr(&writer, "polyOrder");
    mpack_write_i64(&writer, meta.poly_order);

    mpack_write_cstr(&writer, "basisType");
    mpack_write_cstr(&writer, meta.basis_type);
  }

  // 
  if (meta.num_vdim > 0) {
    mpack_write_cstr(&writer, "num_cdim");
    mpack_write_i64(&writer, meta.num_cdim);

    mpack_write_cstr(&writer, "num_vdim");
    mpack_write_i64(&writer, meta.num_vdim);
  }

  mpack_complete_map(&writer);

  int status = mpack_writer_destroy(&writer);

  if (status != mpack_ok) {
    free(mt->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mt);
    mt = 0;
  }

  return mt;
}

// Free pointer for msgpack data.
static void
flash_io_meta_release(struct gkyl_msgpack_data *mt)
{
  if (!mt) return;
  MPACK_FREE(mt->meta);
  gkyl_free(mt);
}

void
gkyl_flash_app_calc_L2(gkyl_flash_app *app, double tm, int frame)
{
  for (int n=0; n<app->num_species; ++n) {
    flash_species_calc_L2(app, tm, &app->species[n]);
  }
}

// Metadata for phase-space DG data (hybrid basis).
static struct gkyl_msgpack_data*
flash_io_meta_phase_new(const gkyl_flash_app *app, double tm, int frame)
{
  return flash_io_meta_new( (struct flash_io_meta) {
      .frame = frame,
      .time = tm,
      .basis_type = app->basis.id,
      .poly_order = app->basis.poly_order,
      .num_cdim = app->cdim,
      .num_vdim = app->vdim,
    }
  );
}

// Metadata for DG data on the kinetic configuration space grid.
static struct gkyl_msgpack_data*
flash_io_meta_conf_new(const gkyl_flash_app *app, double tm, int frame)
{
  return flash_io_meta_new( (struct flash_io_meta) {
      .frame = frame,
      .time = tm,
      .basis_type = app->confBasis.id,
      .poly_order = app->confBasis.poly_order,
    }
  );
}

// Metadata for finite volume data on the fluid grid.
static struct gkyl_msgpack_data*
flash_io_meta_fv_new(const gkyl_flash_app *app, double tm, int frame)
{
  return flash_io_meta_new( (struct flash_io_meta) {
      .frame = frame,
      .time = tm,
    }
  );
}

// Write an array to <app name>-<name><suffix>_<frame>.gkyl.
static void
flash_app_write_array(const gkyl_flash_app *app, const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range, const struct gkyl_msgpack_data *mt,
  const struct gkyl_array *arr, const char *name, const char *suffix, int frame)
{
  const char *fmt = "%s-%s%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, name, suffix, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, name, suffix, frame);

  gkyl_grid_sub_array_write(grid, range, mt, arr, fileNm);
}

// Write species data to gkyl output files.
void
gkyl_flash_app_write_species(gkyl_flash_app *app, int sidx, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct flash_species *s = &app->species[sidx];

  // Write distribution function.
  struct gkyl_msgpack_data *mt = flash_io_meta_phase_new(app, tm, frame);
  flash_app_write_array(app, &s->grid, &s->local, mt, s->f, s->info.name_kinetic, "", frame);
  flash_io_meta_release(mt);

  // Write fluid momentum density if we don't have a fluid species.
  if (app->skip_fluid_species) {
    mt = flash_io_meta_conf_new(app, tm, frame);
    flash_app_write_array(app, &app->grid_kinetic, &app->local_kinetic, mt, s->fluid_rhou,
      s->info.name_kinetic, "_fluid_rhou", frame);
    flash_io_meta_release(mt);
  }
  app->stat.io_tm += gkyl_time_diff_now_sec(wst);
}

// Write fluid species data to output files.
void
gkyl_flash_app_write_fluid_species(gkyl_flash_app *app, int sidx, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct flash_fluid_species *fs = &app->fluid_species[sidx];
  struct gkyl_msgpack_data *mt = flash_io_meta_fv_new(app, tm, frame);

  // Conserved variables, [rho, rho*ux, rho*uy, rho*uz, E].
  flash_app_write_array(app, &app->grid_fluid, &app->local_fluid, mt, fs->fluid,
    fs->info.name_fluid, "", frame);

  // Primitive variables, [rho, ux, uy, uz, p].
  flash_fluid_species_prim_vars(app, fs, fs->fluid, fs->fluid_prim);
  flash_app_write_array(app, &app->grid_fluid, &app->local_fluid, mt, fs->fluid_prim,
    fs->info.name_fluid, "_prim", frame);

  flash_io_meta_release(mt);
  app->stat.io_tm += gkyl_time_diff_now_sec(wst);
}

// Write PKPM diagnostic moments to file.
void
gkyl_flash_app_write_pkpm_moms(gkyl_flash_app* app, int sidx, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct flash_species *s = &app->species[sidx];
  struct flash_fluid_species *fs = &app->fluid_species[sidx];

  // If there's a fluid species, interpolate the fluid velocity to DG
  if (!app->skip_fluid_species) {
    flash_species_interp_rhou_to_dg(app, s, fs, fs->fluid);
  }
  flash_species_calc_pkpm_diags(app, s, s->f);

  struct gkyl_msgpack_data *mt = flash_io_meta_conf_new(app, tm, frame);

  // PKPM diagnostic moments, [rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp].
  flash_app_write_array(app, &app->grid_kinetic, &app->local_kinetic, mt, s->pkpm_moms_diag,
    s->info.name_kinetic, "_pkpm_moms", frame);

  // PKPM primitive variables, [ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp, 3*T_xx/m, 3*T_yy/m, 3*T_zz/m].
  flash_app_write_array(app, &app->grid_kinetic, &app->local_kinetic, mt, s->pkpm_prim,
    s->info.name_kinetic, "_pkpm_prim", frame);

  flash_io_meta_release(mt);
  app->stat.io_tm += gkyl_time_diff_now_sec(wst);
}

// Write electromagnetic field data to file.
void
gkyl_flash_app_write_field(gkyl_flash_app *app, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = flash_io_meta_conf_new(app, tm, frame);
  flash_app_write_array(app, &app->grid_kinetic, &app->local_kinetic, mt, app->field->em,
    "field", "", frame);
  flash_io_meta_release(mt);
  app->stat.io_tm += gkyl_time_diff_now_sec(wst);
}

// Private function to handle variable argument list for printing.
static void
v_flash_app_cout(const gkyl_flash_app *app, FILE *fp, const char *fmt, va_list argp)
{
  int rank;
  gkyl_comm_get_rank(app->comm_kinetic, &rank);
  if ((rank == 0) && fp) {
    vfprintf(fp, fmt, argp);
    fflush(fp);
  }
}

void
gkyl_flash_app_cout(const gkyl_flash_app *app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_flash_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

//  Write species and field data to files. Each write method accumulates io_tm.
void
gkyl_flash_app_write(gkyl_flash_app *app, double tm, int frame)
{
  // Write species data
  for (int i=0; i<app->num_species; ++i) {
    gkyl_flash_app_write_species(app, i, tm, frame);

    if (!app->skip_fluid_species)
      gkyl_flash_app_write_fluid_species(app, i, tm, frame);

    // Compute and write PKPM diagnostic moments if requested
    if (app->species[i].info.write_pkpm_diags) {
      gkyl_flash_app_write_pkpm_moms(app, i, tm, frame);
    }
  }

  // Write field data
  gkyl_flash_app_write_field(app, tm, frame);
}

struct gkyl_flash_stat
gkyl_flash_app_stat(gkyl_flash_app *app)
{
  return app->stat;
}

// Percentage of num relative to den, or 0 if den is 0.
static double
ratio_to_percent(double num, double den)
{
  return den > 0.0 ? 100.0*num/den : 0.0;
}

void
gkyl_flash_app_print_timings(gkyl_flash_app *app, FILE *fp)
{
  const struct gkyl_flash_stat *st = &app->stat;
  double ref = st->fwd_euler_tm;

  gkyl_flash_app_cout(app, fp, "Number of update calls: %d\n", st->nup);
  gkyl_flash_app_cout(app, fp, "Number of forward Euler calls: %d\n", st->nfeuler);
  gkyl_flash_app_cout(app, fp, "Number of failed RK stages: %d\n", st->nfail);
  gkyl_flash_app_cout(app, fp, "Timing:\n");
  gkyl_flash_app_cout(app, fp, "  - Initial conditions:          %.4e sec.\n", st->init_tm);
  gkyl_flash_app_cout(app, fp, "    * Kinetic species:           %.4e sec.\n", st->init_species_tm);
  gkyl_flash_app_cout(app, fp, "    * Fluid species:             %.4e sec.\n", st->init_fluid_species_tm);
  gkyl_flash_app_cout(app, fp, "    * Field:                     %.4e sec.\n", st->init_field_tm);
  gkyl_flash_app_cout(app, fp, "  - Updates:                     %.4e sec.\n", st->total_tm);
  gkyl_flash_app_cout(app, fp, "    * Forward Euler:             %.4e sec. / %6.2f %%\n",
    st->fwd_euler_tm, ratio_to_percent(st->fwd_euler_tm, st->total_tm));
  gkyl_flash_app_cout(app, fp, "      ^ Field variables:         %.4e sec. / %6.2f %%\n",
    st->field_vars_tm, ratio_to_percent(st->field_vars_tm, ref));
  gkyl_flash_app_cout(app, fp, "      ^ Kinetic moments:         %.4e sec. / %6.2f %%\n",
    st->species_mom_tm, ratio_to_percent(st->species_mom_tm, ref));
  gkyl_flash_app_cout(app, fp, "      ^ Collision moments:       %.4e sec. / %6.2f %%\n",
    st->species_coll_mom_tm, ratio_to_percent(st->species_coll_mom_tm, ref));
  gkyl_flash_app_cout(app, fp, "      ^ Kinetic-fluid coupling:  %.4e sec. / %6.2f %%\n",
    st->coupling_tm, ratio_to_percent(st->coupling_tm, ref));
  gkyl_flash_app_cout(app, fp, "      ^ Kinetic RHS:             %.4e sec. / %6.2f %%\n",
    st->species_rhs_tm, ratio_to_percent(st->species_rhs_tm, ref));
  gkyl_flash_app_cout(app, fp, "      ^ Collision RHS:           %.4e sec. / %6.2f %%\n",
    st->species_coll_tm, ratio_to_percent(st->species_coll_tm, ref));
  gkyl_flash_app_cout(app, fp, "      ^ Fluid update:            %.4e sec. / %6.2f %%\n",
    st->fluid_rhs_tm, ratio_to_percent(st->fluid_rhs_tm, ref));
  gkyl_flash_app_cout(app, fp, "      ^ Kinetic BCs:             %.4e sec. / %6.2f %%\n",
    st->species_bc_tm, ratio_to_percent(st->species_bc_tm, ref));
  gkyl_flash_app_cout(app, fp, "      ^ Fluid BCs:               %.4e sec. / %6.2f %%\n",
    st->fluid_species_bc_tm, ratio_to_percent(st->fluid_species_bc_tm, ref));
  gkyl_flash_app_cout(app, fp, "  - I/O:                         %.4e sec.\n", st->io_tm);
}

void
gkyl_flash_app_release(gkyl_flash_app *app) 
{
  for (int i=0; i<app->num_species; ++i) {
    flash_species_release(app, &app->species[i]);

    if (!app->skip_fluid_species) {
      flash_fluid_species_release(&app->fluid_species[i]);
    }
  }
  flash_field_release(app->field);

  gkyl_rect_decomp_release(app->decomp_fluid);
  gkyl_rect_decomp_release(app->decomp_kinetic);
  gkyl_comm_release(app->comm_fluid);
  gkyl_comm_release(app->comm_kinetic);
  gkyl_flash_interp_vars_release(app->interp_vars);

  gkyl_free(app->species);
  gkyl_free(app->fluid_species);
  gkyl_free(app);
}