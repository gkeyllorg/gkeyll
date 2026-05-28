#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gyrokinetic_pol_density.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_dg_interpolate.h>
#include <gkyl_translate_dim.h>
#include <gkyl_proj_on_basis.h>

#include <assert.h>
#include <time.h>

void
gk_species_gyroaverage_disabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gkyl_array *field_in, struct gkyl_array *field_gyroavg)
{
  struct timespec wst = gkyl_wall_clock();
  // Don't perform gyroaveraging.
  app->stat.species_gyroavg_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_gyroaverage_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gkyl_array *field_in, struct gkyl_array *field_gyroavg)
{
  struct timespec wst = gkyl_wall_clock();
  // Gyroaverage input field.
  gkyl_deflated_fem_poisson_advance(species->flr_op, field_in, field_in, field_gyroavg);
  app->stat.species_gyroavg_tm += gkyl_time_diff_now_sec(wst);
}

// Begin static function definitions.
static double
gk_species_omegaH_dt(gkyl_gyrokinetic_app *app, struct gk_species *gks, const struct gkyl_array *fin)
{
  // Compute the time step under which the omega_H mode is stable, dt_omegaH =
  // omegaH_CFL/omega_H.
  // Each species computes its own omega_H as:
  //   omega_H = q_e*sqrt(jacobgeo*n_{s0}/m_s) * omegaH_gf
  // where
  //   - n_{s0} is either a reference, average or max density.
  //   - omegaH_gf = (cmag/(jacobgeo*B^_\parallel))*kpar_max / 
  //                 min(sqrt(k_x^2*eps_xx+k_x*k_y*eps_xy+k_y^2*eps_yy+)).
  // and k_x,k_y,k_par are wavenumbers in computational space, and eps_ij is
  // the polarization weight in our field equation.
 
  if (!(app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN || app->field->gkfield_id == GKYL_GK_FIELD_ADIABATIC)) {
    // Obtain the maximum density (using cell centers).
    gk_species_moment_calc(&gks->m0, gks->local, app->local, fin);
    gkyl_array_reduce_range(gks->m0_max, gks->m0.marr, GKYL_MAX, &app->local);
  
    double m0_max[1];
    if (app->use_gpu) {
      gkyl_cu_memcpy(m0_max, gks->m0_max, sizeof(double), GKYL_CU_MEMCPY_D2H);
    }
    else {
      m0_max[0] = gks->m0_max[0];
    }
    m0_max[0] *= 1.0/pow(sqrt(2.0),app->cdim);
  
    double omegaH = fabs(gks->info.charge)*sqrt(GKYL_MAX2(0.0,m0_max[0])/gks->info.mass)*app->omegaH_gf;
  
    return omegaH > 1e-20? app->cfl_omegaH/omegaH : DBL_MAX;
  }
  else {
    return DBL_MAX;
  }
}

static double
gk_species_rhs_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms)
{
  // Gyroaverage the potential if needed.
  species->gyroaverage(app, species, app->field->phi_smooth, species->gyro_phi);

  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // Collisionless terms.
  gk_species_collisionless_rhs(app, species, &species->collisionless, fin, rhs);

  // Damping term.
  gk_species_damping_advance(app, species, &species->damping, app->field->phi_smooth, fin,
    fbar_in, species->lte.f_lte, rhs, species->cflrate);

  // LBO Collisions.
  gk_species_lbo_rhs(app, species, &species->lbo, fin, rhs);

  // BGK collisions.
  gk_species_bgk_rhs(app, species, &species->bgk, fin, rhs);
  
  // Anomalous diffusion.
  gk_species_anomalous_diff_rhs(app, species, &species->anom_diff, fin, rhs);

  // Line radiation.
  gk_species_radiation_rhs(app, species, &species->rad, fin, rhs);

  // Reactions with charged species.
  gk_species_react_rhs(app, species, &species->react, fin, rhs);

  // Reactions with neutral species.
  gk_species_react_rhs(app, species, &species->react_neut, fin, rhs);

  // BGK source.
  gk_species_source_bgk_rhs(app, species, &species->bgk_src, fin, rhs);

  // Compute and store (in the ghost cell of rhs and in an array in bflux) the boundary fluxes.
  gk_species_bflux_rhs(app, &species->bflux, fin, rhs);

  // Compute moments of the boundary fluxes.
  gk_species_bflux_calc_moms(app, &species->bflux, rhs, bflux_moms);

  // Multiply CFL rate by the df/dt multiplier.
  gk_species_fdot_multiplier_advance_times_cfl(app, species, &species->fdot_mult, app->field->phi_smooth, species->cflrate);
  
  // Reduce the CFL frequency anD compute stable dt needed by this species.
  app->stat.n_species_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omega_cfl, species->cflrate, GKYL_MAX, &species->local);

  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, species->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else {
    omega_cfl_ho[0] = species->omega_cfl[0];
  }
  double dt_out = app->cfl/omega_cfl_ho[0];
  
  // Enforce the omega_H constraint on dt.
  double dt_omegaH = gk_species_omegaH_dt(app, species, fin);
  dt_out = fmin(dt_out, dt_omegaH);

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return dt_out;
}

static double
gk_species_rhs_implicit_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt)
{
  double omega_cfl = 1/DBL_MAX;
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // Compute implicit update and update rhs to new time step
  gk_species_bgk_rhs_implicit(app, species, &species->bgk, fin, dt, rhs);

  gkyl_array_accumulate(gkyl_array_scale(rhs, dt), 1.0, fin);

  app->stat.n_species_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omega_cfl, species->cflrate, GKYL_MAX, &species->local);

  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, species->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else {
    omega_cfl_ho[0] = species->omega_cfl[0];
  }
  omega_cfl = omega_cfl_ho[0];

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return app->cfl/omega_cfl;
}

static double
gk_species_rhs_static(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms)
{
  double omega_cfl = 1/DBL_MAX;
  return app->cfl/omega_cfl;
}

static double
gk_species_rhs_implicit_static(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt)
{
  double omega_cfl = 1/DBL_MAX;
  return app->cfl/omega_cfl;
}

static void
gk_species_apply_bc_dynamic(gkyl_gyrokinetic_app *app, const struct gk_species *species, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();
  
  int num_periodic_dir = species->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(species->comm, &species->local, &species->local_ext,
    num_periodic_dir, species->periodic_dirs, f); 
  
  for (int d=0; d<cdim; ++d) {
    if (species->bc_is_np[d]) {

      switch (species->lower_bc[d].type) {
        case GKYL_BC_GK_SPECIES_SHEATH:
          gkyl_bc_sheath_gyrokinetic_advance(species->bc_sheath_lo, app->field->phi_smooth, 
            app->field->phi_wall_lo, f, &app->local);
          break;
        case GKYL_BC_GK_SPECIES_IWL:
          gkyl_bc_sheath_gyrokinetic_advance(species->bc_sheath_lo, app->field->phi_smooth, 
            app->field->phi_wall_lo, f, &app->local);
          if (cdim == 3) {
            gkyl_bc_twistshift_advance(species->bc_ts_lo, f, f);
          }
          break;
        case GKYL_BC_GK_SPECIES_COPY:
        case GKYL_BC_GK_SPECIES_REFLECT:
        case GKYL_BC_GK_SPECIES_ABSORB:
          gkyl_bc_basic_gyrokinetic_advance(species->bc_lo[d], species->bc_buffer, f);
          break;
        case GKYL_BC_GK_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_gyrokinetic_advance(species->bc_lo[d], species->bc_buffer_lo_fixed, f);
          break;
        case GKYL_BC_GK_SPECIES_ZERO_FLUX:
          break; // do nothing, BCs already applied in hyper_dg loop by not updating flux
        default:
          break;
      }

      switch (species->upper_bc[d].type) {
        case GKYL_BC_GK_SPECIES_SHEATH:
          gkyl_bc_sheath_gyrokinetic_advance(species->bc_sheath_up, app->field->phi_smooth, 
            app->field->phi_wall_up, f, &app->local);
          break;
        case GKYL_BC_GK_SPECIES_IWL:
          gkyl_bc_sheath_gyrokinetic_advance(species->bc_sheath_up, app->field->phi_smooth, 
            app->field->phi_wall_up, f, &app->local);
          if (cdim == 3) {
            gkyl_bc_twistshift_advance(species->bc_ts_up, f, f);
          }
          break;
        case GKYL_BC_GK_SPECIES_COPY:
        case GKYL_BC_GK_SPECIES_REFLECT:
        case GKYL_BC_GK_SPECIES_ABSORB:
          gkyl_bc_basic_gyrokinetic_advance(species->bc_up[d], species->bc_buffer, f);
          break;
        case GKYL_BC_GK_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_gyrokinetic_advance(species->bc_up[d], species->bc_buffer_up_fixed, f);
          break;
        case GKYL_BC_GK_SPECIES_ZERO_FLUX:
          break; // do nothing, BCs already applied in hyper_dg loop by not updating flux
        default:
          break;
      }      
    }
  }

  gkyl_comm_array_sync(species->comm, &species->local, &species->local_ext, f);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_apply_bc_static(gkyl_gyrokinetic_app *app, const struct gk_species *species, struct gkyl_array *f)
{
  // do nothing
}

static void
gk_species_step_f_dynamic(struct gkyl_array* out, double dt,
  const struct gkyl_array* inp)
{
  gkyl_array_accumulate(gkyl_array_scale(out, dt), 1.0, inp);
}

static void
gk_species_step_f_static(struct gkyl_array* out, double dt,
  const struct gkyl_array* inp)
{
  // do nothing
}

static void
gk_species_combine_dynamic(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  gkyl_array_accumulate_range(gkyl_array_set_range(out, c1, arr1, rng),
    c2, arr2, rng);
}

static void
gk_species_combine_static(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  // do nothing
}

static void
gk_species_copy_range_dynamic(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  gkyl_array_copy_range(out, inp, range);
}

static void
gk_species_copy_range_static(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  // do nothing
}

static void
gk_species_write_dynamic(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  gks->write_cfl_func(app, gks, tm, frame);

  // Metadata from app and geo object.
  gkyl_msgpack_map_elem_set_double(app->io_meta_basic_len, app->io_meta_basic, "time", tm);
  gkyl_msgpack_map_elem_set_uint(app->io_meta_basic_len, app->io_meta_basic, "frame", frame);
  int io_meta_len[] = {app->io_meta_basic_len,gks->io_meta_len,app->gk_geom->io_meta_len};
  const struct gkyl_msgpack_map_elem* io_meta[] = {app->io_meta_basic,gks->io_meta,app->gk_geom->io_meta};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);
  
  // copy data from device to host before writing it out
  if (app->use_gpu) {
    gkyl_array_copy(gks->f_host, gks->f);
  }
  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->f_host, fileNm);
    
  gkyl_msgpack_data_release(mt);  

  app->stat.species_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_io += 1;
}

static void
gk_species_write_static(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  // do nothing
}

static void
gk_species_write_cfl_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  // DG metadata for cflrate.
  struct gkyl_msgpack_map_elem mpe_cfl[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = "serendipity" },
  };
  int mpe_cfl_len = sizeof(mpe_cfl)/sizeof(mpe_cfl[0]);
  // Update app basic metada with time/frame.
  gkyl_msgpack_map_elem_set_double(app->io_meta_basic_len, app->io_meta_basic, "time", tm);
  gkyl_msgpack_map_elem_set_uint(app->io_meta_basic_len, app->io_meta_basic, "frame", frame);
  // Unionize metadata.
  int io_meta_len[] = {app->io_meta_basic_len, mpe_cfl_len, app->gk_geom->io_meta_len};
  const struct gkyl_msgpack_map_elem* io_meta[] = {app->io_meta_basic, mpe_cfl, app->gk_geom->io_meta};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  const char *fmt = "%s-%s_cflrate_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);
  gkyl_array_copy(gks->cflrate_ho, gks->cflrate);
  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt,
    gks->cflrate_ho, fileNm);

  gkyl_msgpack_data_release(mt);  

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_io += 1;
}

static void
gk_species_write_cfl_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  // do nothing
}

static void
gk_species_write_mom_dynamic(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  // Update and package metadata.
  gkyl_msgpack_map_elem_set_double(app->io_meta_basic_len, app->io_meta_basic, "time", tm);
  gkyl_msgpack_map_elem_set_uint(app->io_meta_basic_len, app->io_meta_basic, "frame", frame);
  int io_meta_len[] = {app->io_meta_basic_len,app->io_meta_len,app->gk_geom->io_meta_len};
  const struct gkyl_msgpack_map_elem* io_meta[] = {app->io_meta_basic,app->io_meta,app->gk_geom->io_meta};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  for (int m=0; m<gks->info.num_diag_moments; ++m) {
    struct timespec wtm = gkyl_wall_clock();
    gk_species_moment_calc(&gks->moms[m], gks->local, app->local, gks->f);
    app->stat.n_mom += 1;
    
    // Rescale moment by inverse of Jacobian if necessary. 
    gk_species_moment_diag_jacobgeo_div(app, &gks->moms[m], gks->moms[m].marr, gks->moms[m].marr);
    app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wtm);
      
    struct timespec wst = gkyl_wall_clock();
    if (app->use_gpu) {
      gkyl_array_copy(gks->moms[m].marr_host, gks->moms[m].marr);
    }

    const char *fmt = "%s-%s_%s_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name,
      gkyl_distribution_moments_strs[gks->info.diag_moments[m]], frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name,
      gkyl_distribution_moments_strs[gks->info.diag_moments[m]], frame);
    
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
      gks->moms[m].marr_host, fileNm);
    app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
    app->stat.n_diag_io += 1;
  }
  gkyl_msgpack_data_release(mt); 

  app->stat.n_diag += 1;
}

static void
gk_species_write_mom_static(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  // do nothing
}

static void
gk_species_calc_int_mom_dt_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double dt, struct gkyl_array *fdot_int_mom)
{
  struct timespec wst = gkyl_wall_clock();
  // Compute moment of f_new to compute moment of df/dt.
  // Need to do it after the fields are updated.
  gk_species_moment_calc(&gks->integ_moms, gks->local, app->local, gks->f); 
  gkyl_array_set(fdot_int_mom, 1.0/dt, gks->integ_moms.marr);
  app->stat.fdot_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_calc_int_mom_dt_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double dt, struct gkyl_array *fdot_int_mom)
{
}

void
gk_species_calc_int_mom_dt(gkyl_gyrokinetic_app* app, struct gk_species *gks, double dt, struct gkyl_array *fdot_int_mom)
{
  gks->calc_int_mom_dt_func(app, gks, dt, fdot_int_mom);
}

static void
gk_species_calc_integrated_mom_dynamic(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  struct timespec wst = gkyl_wall_clock();

  int num_mom = gks->integ_moms.num_mom;
  double avals_global[num_mom];
  
  gk_species_moment_calc(&gks->integ_moms, gks->local, app->local, gks->f); 
  app->stat.n_mom += 1;

  // Reduce (sum) over whole domain, append to diagnostics.
  gkyl_array_reduce_range(gks->red_integ_diag, gks->integ_moms.marr, GKYL_SUM, &app->local);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, num_mom, 
    gks->red_integ_diag, gks->red_integ_diag_global);
  if (app->use_gpu) {
    gkyl_cu_memcpy(avals_global, gks->red_integ_diag_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    memcpy(avals_global, gks->red_integ_diag_global, sizeof(double[num_mom]));
  }
  gkyl_dynvec_append(gks->integ_diag, tm, avals_global);

  if (gks->info.time_rate_diagnostics) {
    // Reduce (sum) over whole domain, append to diagnostics.
    gkyl_array_accumulate(gks->fdot_mom_new, -1.0, gks->fdot_mom_old);
    gkyl_array_reduce_range(gks->red_integ_diag, gks->fdot_mom_new, GKYL_SUM, &app->local);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, num_mom,
      gks->red_integ_diag, gks->red_integ_diag_global);
    if (app->use_gpu) {
      gkyl_cu_memcpy(avals_global, gks->red_integ_diag_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      memcpy(avals_global, gks->red_integ_diag_global, sizeof(double[num_mom]));
    }
    gkyl_dynvec_append(gks->fdot_integ_diag, tm, avals_global);
  }

  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag += 1;
}

static void
gk_species_calc_integrated_mom_static(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  // do nothing
}

static void
gk_species_write_integrated_mom_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  struct timespec wst = gkyl_wall_clock();

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // Write integrated diagnostic moments.
    const char *fmt = "%s-%s_%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "integrated_moms");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "integrated_moms");
    
    if (gks->is_first_integ_write_call) {
      gkyl_dynvec_write(gks->integ_diag, fileNm);
      gks->is_first_integ_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(gks->integ_diag, fileNm);
    }
  }
  gkyl_dynvec_clear(gks->integ_diag);
  app->stat.n_diag_io += 1;
  
  if (gks->info.time_rate_diagnostics) {
    if (rank == 0) {
      // Write integrated diagnostic moments.
      const char *fmt = "%s-%s_fdot_%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "integrated_moms");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "integrated_moms");

      if (gks->is_first_fdot_integ_write_call) {
        gkyl_dynvec_write(gks->fdot_integ_diag, fileNm);
        gks->is_first_fdot_integ_write_call = false;
      }
      else {
        gkyl_dynvec_awrite(gks->fdot_integ_diag, fileNm);
      }
    }
    gkyl_dynvec_clear(gks->fdot_integ_diag);
    app->stat.n_diag_io += 1;
  }

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_write_integrated_mom_static(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  // do nothing
}

static void
gk_species_calc_L2norm_dynamic(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_array_integrate_advance(gks->integ_wfsq_op, gks->f, (2.0*M_PI)/gks->info.mass,
    app->jacobtot_inv_weak, &gks->local, &app->local, gks->L2norm_local);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, gks->L2norm_local, gks->L2norm_global);

  double L2norm_global[] = {0.0};
  if (app->use_gpu) {
    gkyl_cu_memcpy(L2norm_global, gks->L2norm_global, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else { 
    memcpy(L2norm_global, gks->L2norm_global, sizeof(double));
  }
  gkyl_dynvec_append(gks->L2norm, tm, L2norm_global);  

  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag += 1;
}

static void
gk_species_calc_L2norm_static(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  // do nothing
}

static void
gk_species_write_L2norm_dynamic(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  struct timespec wst = gkyl_wall_clock();

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // Write the L2 norm.
    const char *fmt = "%s-%s_%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "L2norm");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "L2norm");

    if (gks->is_first_L2norm_write_call) {
      gkyl_dynvec_write(gks->L2norm, fileNm);
      gks->is_first_L2norm_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(gks->L2norm, fileNm);
    }
  }
  gkyl_dynvec_clear(gks->L2norm);

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag_io += 1;
}

static void
gk_species_write_L2norm_static(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  // do nothing
}

static void
gk_species_release_dynamic(const gkyl_gyrokinetic_app* app, const struct gk_species *s)
{
  // Release various arrays and objects for a dynamic species.
  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->bc_buffer);
  gkyl_array_release(s->bc_buffer_lo_fixed);
  gkyl_array_release(s->bc_buffer_up_fixed);

  if (s->info.write_omega_cfl) {
    gkyl_array_release(s->cflrate_ho);
  }

  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    if (s->lower_bc[d].type == GKYL_BC_GK_SPECIES_SHEATH) {
      gkyl_bc_sheath_gyrokinetic_release(s->bc_sheath_lo);
    }
    else if (s->lower_bc[d].type == GKYL_BC_GK_SPECIES_IWL) { 
      gkyl_bc_sheath_gyrokinetic_release(s->bc_sheath_lo);
      if (app->cdim == 3) {
        gkyl_bc_twistshift_release(s->bc_ts_lo);
      }
    }
    else if ( (s->lower_bc[d].type == GKYL_BC_GK_SPECIES_COPY) ||
              (s->lower_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB) ||
              (s->lower_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
              (s->lower_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) ) {
      gkyl_bc_basic_gyrokinetic_release(s->bc_lo[d]);
    }
    
    if (s->upper_bc[d].type == GKYL_BC_GK_SPECIES_SHEATH) {
      gkyl_bc_sheath_gyrokinetic_release(s->bc_sheath_up);
    }
    else if (s->upper_bc[d].type == GKYL_BC_GK_SPECIES_IWL) {
      gkyl_bc_sheath_gyrokinetic_release(s->bc_sheath_up);
      if (app->cdim == 3) {
        gkyl_bc_twistshift_release(s->bc_ts_up);
      }
    }
    else if ( (s->upper_bc[d].type == GKYL_BC_GK_SPECIES_COPY) ||
              (s->upper_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB) ||
              (s->upper_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
              (s->upper_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) ) {
      gkyl_bc_basic_gyrokinetic_release(s->bc_up[d]);
    }
  }
  
  if (app->use_gpu) {
    gkyl_cu_free(s->omega_cfl);
    gkyl_cu_free(s->m0_max);
  }
  else {
    gkyl_free(s->omega_cfl);
    gkyl_free(s->m0_max);
  }

  // Release integrated moment memory.
  gk_species_moment_release(app, &s->integ_moms); 

  // Release integrated diag memory.
  gkyl_dynvec_release(s->integ_diag);
  if (app->use_gpu) {
    gkyl_cu_free(s->red_integ_diag);
    gkyl_cu_free(s->red_integ_diag_global);
  }
  else {
    gkyl_free(s->red_integ_diag);
    gkyl_free(s->red_integ_diag_global);
  }

  // Release L2 norm memory.
  gkyl_array_integrate_release(s->integ_wfsq_op);
  gkyl_dynvec_release(s->L2norm);
  if (app->use_gpu) {
    gkyl_cu_free(s->L2norm_local);
    gkyl_cu_free(s->L2norm_global);
  }
  else {
    gkyl_free(s->L2norm_local);
    gkyl_free(s->L2norm_global);
  }

  if (s->info.time_rate_diagnostics) {
    // Free df/dt diagnostics memory.
    gkyl_array_release(s->fdot_mom_old);
    gkyl_array_release(s->fdot_mom_new);
    gkyl_dynvec_release(s->fdot_integ_diag);
  }
}

static void
gk_species_release_static(const gkyl_gyrokinetic_app* app, const struct gk_species *s)
{
}

// Initialize species object.
static void
gk_species_init_dynamic(struct gkyl_gk *gk_app_inp, struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  int cdim = app->cdim, vdim = gks->info.vdim;
  int pdim = cdim+vdim;

  int ghost[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    // Full phase space grid.
    ghost[cdim+d] = 0; // No ghost-cells in velocity space.
  }

  // Allocate distribution function arrays.
  gks->f1 = mkarr(app->use_gpu, gks->basis.num_basis, gks->local_ext.volume);
  gks->fnew = mkarr(app->use_gpu, gks->basis.num_basis, gks->local_ext.volume);

  if (app->use_gpu) {
    gks->omega_cfl = gkyl_cu_malloc(sizeof(double));
    gks->m0_max = gkyl_cu_malloc(app->basis.num_basis*sizeof(double));
  }
  else {
    gks->omega_cfl = gkyl_malloc(sizeof(double));
    gks->m0_max = gkyl_malloc(app->basis.num_basis*sizeof(double));
  }

  // Allocate data for integrated moments.
  int num_diag_int_moms = gks->info.num_integrated_diag_moments;
  assert(num_diag_int_moms < 2); // 1 int moment allowed now.
  gk_species_moment_init(app, gks, &gks->integ_moms,
    num_diag_int_moms == 0? GKYL_F_MOMENT_M0M1M2PARM2PERP : gks->info.integrated_diag_moments[0], true);

  if (app->use_gpu) {
    gks->red_integ_diag = gkyl_cu_malloc(sizeof(double[gks->integ_moms.num_mom]));
    gks->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[gks->integ_moms.num_mom]));
  } 
  else {
    gks->red_integ_diag = gkyl_malloc(sizeof(double[gks->integ_moms.num_mom]));
    gks->red_integ_diag_global = gkyl_malloc(sizeof(double[gks->integ_moms.num_mom]));
  }
  // Allocate dynamic-vector to store all-reduced integrated moments.
  gks->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, gks->integ_moms.num_mom);
  gks->is_first_integ_write_call = true;

  // Allocate dynamic-vector to store Delta f integrated moments.
  if (gks->info.time_rate_diagnostics) {
    gks->fdot_mom_old = mkarr(app->use_gpu, gks->integ_moms.marr->ncomp, gks->integ_moms.marr->size);
    gks->fdot_mom_new = mkarr(app->use_gpu, gks->integ_moms.marr->ncomp, gks->integ_moms.marr->size);
    gks->fdot_integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, gks->integ_moms.num_mom);
    gks->is_first_fdot_integ_write_call = true;
  }

  // Objects for L2 norm diagnostic.
  gks->integ_wfsq_op = gkyl_array_integrate_new(&gks->grid, &gks->basis, 1, GKYL_ARRAY_INTEGRATE_OP_SQ_WEIGHTED, app->use_gpu);
  if (app->use_gpu) {
    gks->L2norm_local = gkyl_cu_malloc(sizeof(double));
    gks->L2norm_global = gkyl_cu_malloc(sizeof(double));
  } 
  else {
    gks->L2norm_local = gkyl_malloc(sizeof(double));
    gks->L2norm_global = gkyl_malloc(sizeof(double));
  }
  gks->L2norm = gkyl_dynvec_new(GKYL_DOUBLE, 1); // L2 norm.
  gks->is_first_L2norm_write_call = true;
  
  for (int dir=0; dir<app->cdim; ++dir) {
    if (gks->lower_bc[dir].type == GKYL_BC_GK_SPECIES_IWL || gks->upper_bc[dir].type == GKYL_BC_GK_SPECIES_IWL) {
      // Make the parallel direction periodic so that we sync the core before
      // applying sheath BCs in the SOL.
      gks->periodic_dirs[gks->num_periodic_dir] = app->cdim-1; // The last direction is the parallel one.
      gks->num_periodic_dir += 1;
    }
  }
  
  // Allocate buffer needed for BCs.
  bool need_bc_buffer = false;
  bool need_bc_buffer_lo_fixed = false;
  bool need_bc_buffer_up_fixed = false;
  for (int d=0; d<cdim; ++d) {
    if (gks->bc_is_np[d]) {
      if ( (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_COPY   ) ||
           (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
           (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB ) ||
           (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_COPY   ) ||
           (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
           (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB ) ) {
        need_bc_buffer = true;
      }
      if (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) {
        need_bc_buffer_lo_fixed = true;
      }
      if (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) {
        need_bc_buffer_up_fixed = true;
      }
    }
  }

  long buff_sz = 1;
  for (int dir=0; dir<cdim; ++dir) {
    long vol = GKYL_MAX2(gks->local_lower_skin[dir].volume, gks->local_upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  gks->bc_buffer = need_bc_buffer? mkarr(app->use_gpu, gks->basis.num_basis, buff_sz)
                                 : mkarr(app->use_gpu, 1, 1);
  // Buffer arrays for fixed function BCs on distribution function.
  gks->bc_buffer_lo_fixed = need_bc_buffer_lo_fixed? mkarr(app->use_gpu, gks->basis.num_basis, buff_sz)
                                                   : mkarr(app->use_gpu, 1, 1);
  gks->bc_buffer_up_fixed = need_bc_buffer_up_fixed? mkarr(app->use_gpu, gks->basis.num_basis, buff_sz)
                                                   : mkarr(app->use_gpu, 1, 1);

  for (int d=0; d<cdim; ++d) {
    // Lower BC.
    if (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_SHEATH) {
      gks->bc_sheath_lo = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_LOWER_EDGE, gks->basis_on_dev, 
        &gks->local_lower_skin[d], &gks->local_lower_ghost[d], gks->vel_map,
        cdim, 2.0*(gks->info.charge/gks->info.mass), app->use_gpu);
    }
    else if (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_IWL) {
      gks->bc_sheath_lo = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_LOWER_EDGE, gks->basis_on_dev, 
        &gks->local_lower_skin_par_sol, &gks->local_lower_ghost_par_sol, gks->vel_map,
        cdim, 2.0*(gks->info.charge/gks->info.mass), app->use_gpu);

      if (cdim == 3) {
        // For 3x2v we need a twistshift BC in the core.
        struct gkyl_bc_twistshift_inp tsinp = {
          .bc_dir = d,
          .shift_dir = 1, // y shift.
          .shear_dir = 0, // shift varies with x.
          .edge = GKYL_LOWER_EDGE,
          .cdim = cdim,
          .bcdir_ext_update_r = gks->local_par_ext_core,
          .num_ghost = ghost,
          .basis = gks->basis,
          .grid = gks->grid,
          .shift_func = gks->lower_bc[d].aux_profile,
          .shift_func_ctx = gks->lower_bc[d].aux_ctx,
          .use_gpu = app->use_gpu,
        };

        gks->bc_ts_lo = gkyl_bc_twistshift_new(&tsinp);
      }
    }
    else if ( (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_COPY) ||
              (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB) ||
              (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
              (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) ) {

      if (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT)
        assert(d == cdim-1); // MF 2025/09/25: Reflecting BCs only allowed in the parallel direction.

      gks->bc_lo[d] = gkyl_bc_basic_gyrokinetic_new(d, GKYL_LOWER_EDGE, gks->lower_bc[d].type, gks->basis_on_dev,
        &gks->local_lower_skin[d], &gks->local_lower_ghost[d], gks->f->ncomp, app->cdim, app->use_gpu);

      if (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) {
        // Project distribution desired in the ghost cell.
        struct gk_proj gk_proj_bc_lo;
        gk_species_projection_init(app, gks, gks->lower_bc[d].projection, &gk_proj_bc_lo);
        gk_species_projection_calc(app, gks, &gk_proj_bc_lo, gks->f1, 0.0); // Temporarily use f1.
        // Fill ghost cell with the boundary value.
        struct gkyl_bc_basic_gyrokinetic *gfss_bc_op = gkyl_bc_basic_gyrokinetic_new(d, GKYL_LOWER_EDGE,
          GKYL_BC_GK_SPECIES_BOUNDARY_VALUE, gks->basis_on_dev, &gks->local_lower_skin[d], &gks->local_lower_ghost[d],
	  gks->f->ncomp, app->cdim, app->use_gpu);
        gkyl_bc_basic_gyrokinetic_advance(gfss_bc_op, gks->bc_buffer_lo_fixed, gks->f1);
        gkyl_bc_basic_gyrokinetic_release(gfss_bc_op);
	// Copy ghost range into buffer.
        gkyl_bc_basic_gyrokinetic_buffer_fixed_func(gks->bc_lo[d], gks->bc_buffer_lo_fixed, gks->f1);
        gkyl_array_clear(gks->f1, 0.0);
        gk_species_projection_release(app, &gk_proj_bc_lo);
      }
    }

    // Upper BC.
    if (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_SHEATH) {
      gks->bc_sheath_up = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_UPPER_EDGE, gks->basis_on_dev, 
        &gks->local_upper_skin[d], &gks->local_upper_ghost[d], gks->vel_map,
        cdim, 2.0*(gks->info.charge/gks->info.mass), app->use_gpu);
    }
    else if (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_IWL) {
      gks->bc_sheath_up = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_UPPER_EDGE, gks->basis_on_dev, 
        &gks->local_upper_skin_par_sol, &gks->local_upper_ghost_par_sol, gks->vel_map,
        cdim, 2.0*(gks->info.charge/gks->info.mass), app->use_gpu);

      if (cdim == 3) {
        // For 3x2v we need a twistshift BC in the core.
        struct gkyl_bc_twistshift_inp tsinp = {
          .bc_dir = d,
          .shift_dir = 1, // y shift.
          .shear_dir = 0, // shift varies with x.
          .edge = GKYL_UPPER_EDGE,
          .cdim = cdim,
          .bcdir_ext_update_r = gks->local_par_ext_core,
          .num_ghost = ghost,
          .basis = gks->basis,
          .grid = gks->grid,
          .shift_func = gks->upper_bc[d].aux_profile,
          .shift_func_ctx = gks->upper_bc[d].aux_ctx,
          .use_gpu = app->use_gpu,
        };

        gks->bc_ts_up = gkyl_bc_twistshift_new(&tsinp);
      }
    }
    else if ( (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_COPY) ||
              (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB) ||
              (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
              (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) ) {

      if (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT)
        assert(d == cdim-1); // MF 2025/09/25: Reflecting BCs only allowed in the parallel direction.

      gks->bc_up[d] = gkyl_bc_basic_gyrokinetic_new(d, GKYL_UPPER_EDGE, gks->upper_bc[d].type, gks->basis_on_dev,
        &gks->local_upper_skin[d], &gks->local_upper_ghost[d], gks->f->ncomp, app->cdim, app->use_gpu);

      if (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) {
        // Project distribution desired in the ghost cell.
        struct gk_proj gk_proj_bc_up;
        gk_species_projection_init(app, gks, gks->upper_bc[d].projection, &gk_proj_bc_up);
        gk_species_projection_calc(app, gks, &gk_proj_bc_up, gks->f1, 0.0); // Temporarily use f1.
        // Fill ghost cell with the boundary value.
        struct gkyl_bc_basic_gyrokinetic *gfss_bc_op = gkyl_bc_basic_gyrokinetic_new(d, GKYL_UPPER_EDGE,
          GKYL_BC_GK_SPECIES_BOUNDARY_VALUE, gks->basis_on_dev, &gks->local_upper_skin[d], &gks->local_upper_ghost[d],
	  gks->f->ncomp, app->cdim, app->use_gpu);
        gkyl_bc_basic_gyrokinetic_advance(gfss_bc_op, gks->bc_buffer_up_fixed, gks->f1);
        gkyl_bc_basic_gyrokinetic_release(gfss_bc_op);
	// Copy ghost range into buffer.
        gkyl_bc_basic_gyrokinetic_buffer_fixed_func(gks->bc_up[d], gks->bc_buffer_up_fixed, gks->f1);
        gkyl_array_clear(gks->f1, 0.0);
        gk_species_projection_release(app, &gk_proj_bc_up);
      }
    }
  }

  // Set function pointers.
  gks->rhs_func = gk_species_rhs_dynamic;
  gks->rhs_implicit_func = gk_species_rhs_implicit_dynamic;
  gks->bc_func = gk_species_apply_bc_dynamic;
  gks->release_func = gk_species_release_dynamic;
  gks->step_f_func = gk_species_step_f_dynamic;
  gks->combine_func = gk_species_combine_dynamic;
  gks->copy_func = gk_species_copy_range_dynamic;
  gks->write_func = gk_species_write_dynamic;
  if (gks->info.write_omega_cfl) {
    gks->cflrate_ho = mkarr(false, gks->cflrate->ncomp, gks->cflrate->size);
    gks->write_cfl_func = gk_species_write_cfl_enabled;
  }
  else 
    gks->write_cfl_func = gk_species_write_cfl_disabled;
  gks->write_mom_func = gk_species_write_mom_dynamic;
  gks->calc_integrated_mom_func = gk_species_calc_integrated_mom_dynamic;
  gks->write_integrated_mom_func = gk_species_write_integrated_mom_dynamic;
  gks->calc_L2norm_func = gk_species_calc_L2norm_dynamic;
  gks->write_L2norm_func = gk_species_write_L2norm_dynamic;
  if (gks->info.time_rate_diagnostics)
    gks->calc_int_mom_dt_func = gk_species_calc_int_mom_dt_enabled;
  else
    gks->calc_int_mom_dt_func = gk_species_calc_int_mom_dt_disabled;
}

// Initialize static species object.
static void
gk_species_init_static(struct gkyl_gk *gk_app_inp, struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // Allocate distribution function arrays.
  gks->f1 = gks->f;
  gks->fnew = gks->f;
  
  // Set function pointers.
  gks->rhs_func = gk_species_rhs_static;
  gks->rhs_implicit_func = gk_species_rhs_implicit_static;
  gks->bc_func = gk_species_apply_bc_static;
  gks->release_func = gk_species_release_static;
  gks->step_f_func = gk_species_step_f_static;
  gks->combine_func = gk_species_combine_static;
  gks->copy_func = gk_species_copy_range_static;
  gks->write_func = gk_species_write_static;
  gks->write_cfl_func = gk_species_write_cfl_disabled;
  gks->write_mom_func = gk_species_write_mom_static;
  gks->calc_integrated_mom_func = gk_species_calc_integrated_mom_static;
  gks->write_integrated_mom_func = gk_species_write_integrated_mom_static;
  gks->calc_L2norm_func = gk_species_calc_L2norm_static;
  gks->write_L2norm_func = gk_species_write_L2norm_static;
  gks->calc_int_mom_dt_func = gk_species_calc_int_mom_dt_disabled;
}

void
gk_species_file_import_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, 
  struct gkyl_gyrokinetic_ic_import inp)
{
  // Import initial condition from a file. Intended options include importing:
  //   1) ICs with same grid.
  //   2) ICs one dimensionality lower (e.g. 2x2v for 3x2v sim).
  //   3) ICs with same grid extents but different resolution (NYI).

  struct gkyl_rect_grid grid = gks->grid;
  int pdim = grid.ndim;
  int cdim = app->cdim;
  int vdim = pdim - cdim;
  int poly_order = app->poly_order;

  struct gkyl_rect_grid grid_do; // Donor grid.
  struct gkyl_array_header_info hdr;
  int pdim_do, vdim_do, cdim_do;
  bool same_res = true;

  // Read the header of the input file, extract needed info an create a grid
  // and other things needed.
  FILE *fp;
  with_file(fp, inp.file_name, "r") {

    int status = gkyl_grid_sub_array_header_read_fp(&grid_do, &hdr, fp);

    pdim_do = grid_do.ndim;
    vdim_do = vdim; // Assume velocity space dimensionality is the same.
    cdim_do = pdim_do - vdim_do;

    // Perform some basic checks.
    if (pdim_do == pdim) {
      // Check if the grid resolution is the same.
      for (int d=0; d<pdim; d++)
        same_res = same_res && (grid_do.cells[d] == grid.cells[d]);
    }
    else {
      // Assume the loaded file has one lower conf-space dimension.
      // Primarily meant for loading:
      //   - 1x2v for a 2x2v sim.
      //   - 2x2v for a 3x2v sim.
      assert(pdim_do == pdim-1);
      double tol_eq = 1e-10;
      for (int d=0; d<cdim_do-1; d++) {
        assert(gkyl_compare_double(grid_do.lower[d], grid.lower[d], tol_eq));
        assert(gkyl_compare_double(grid_do.upper[d], grid.upper[d], tol_eq));
        assert(gkyl_compare_double(grid_do.dx[d]   , grid.dx[d]   , tol_eq));
        assert(grid_do.cells[d] == grid.cells[d]);
      }
      assert(gkyl_compare_double(grid_do.lower[cdim_do-1], grid.lower[cdim-1], tol_eq));
      assert(gkyl_compare_double(grid_do.upper[cdim_do-1], grid.upper[cdim-1], tol_eq));
      assert(gkyl_compare_double(grid_do.dx[cdim_do-1]   , grid.dx[cdim-1]   , tol_eq));
      assert(grid_do.cells[cdim_do-1] == grid.cells[cdim-1]);
      for (int d=0; d<vdim; d++) {
        assert(gkyl_compare_double(grid_do.lower[cdim_do+d], grid.lower[cdim+d], tol_eq));
        assert(gkyl_compare_double(grid_do.upper[cdim_do+d], grid.upper[cdim+d], tol_eq));
        assert(gkyl_compare_double(grid_do.dx[cdim_do+d]   , grid.dx[cdim+d]   , tol_eq));
        assert(grid_do.cells[cdim_do+d] == grid.cells[cdim+d]);
      }
    }

    // Read basis info from header, check its consistency.
    struct gkyl_msgpack_map_elem elem_list[] = {
      { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
      { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = 0 },
    };
    int elem_list_len = sizeof(elem_list)/sizeof(elem_list[0]);
    gkyl_msgpack_to_map_elem_list(&(struct gkyl_msgpack_data) {
        .meta = hdr.meta,
        .meta_sz = hdr.meta_size
      }, elem_list_len, elem_list);
    assert(strcmp(gks->basis.id, gkyl_msgpack_map_elem_get_string(elem_list_len, elem_list, "basis_type")) == 0);
    assert(poly_order == gkyl_msgpack_map_elem_get_uint(elem_list_len, elem_list, "poly_order"));
    gkyl_msgpack_map_elem_release_string(elem_list_len, elem_list, "basis_type");

    gkyl_grid_sub_array_header_release(&hdr);
  }

  // Donor basis.
  struct gkyl_basis basis_do;
  switch (gks->basis.b_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (poly_order > 1) {
        gkyl_cart_modal_serendip(&basis_do, pdim_do, poly_order);
      }
      else if (poly_order == 1) {
        // p=2 in vparallel
        gkyl_cart_modal_gkhybrid(&basis_do, cdim_do, vdim_do); 
      }
      break;
    default:
      assert(false);
      break;
  }

  // Donor global range.
  int ghost_do[pdim_do];
  for (int d=0; d<cdim_do; d++) ghost_do[d] = 1;
  for (int d=0; d<vdim_do; d++) ghost_do[cdim_do+d] = 0;
  struct gkyl_range global_ext_do, global_do;
  gkyl_create_grid_ranges(&grid_do, ghost_do, &global_ext_do, &global_do);

  // Create a donor communicator.
  int cuts_tar[GKYL_MAX_DIM] = {-1}, cuts_do[GKYL_MAX_CDIM] = {-1};
  gkyl_rect_decomp_get_cuts(app->decomp, cuts_tar);
  if (cdim_do == cdim-1) {
    for (int d=0; d<cdim_do-1; d++) {
      cuts_do[d] = cuts_tar[d];
    }
    cuts_do[cdim_do-1] = cuts_tar[cdim-1];
  }
  else {
    for (int d=0; d<cdim; d++) {
      cuts_do[d] = cuts_tar[d];
    }
  }
  // Set velocity space cuts to 1 as we do not use MPI in vel-space.
  for (int d=0; d<vdim; d++) {
    cuts_do[cdim_do+d] = 1;
  }

  struct gkyl_rect_decomp *decomp_do = gkyl_rect_decomp_new_from_cuts(pdim_do, cuts_do, &global_do);
  struct gkyl_comm* comm_do = gkyl_comm_split_comm(gks->comm, 0, decomp_do);

  // Donor local range.
  int my_rank = 0;
  gkyl_comm_get_rank(comm_do, &my_rank);

  struct gkyl_range local_ext_do, local_do;
  gkyl_create_ranges(&decomp_do->ranges[my_rank], ghost_do, &local_ext_do, &local_do);

  // Donor array.
  struct gkyl_array *fdo = mkarr(app->use_gpu, basis_do.num_basis, local_ext_do.volume);
  struct gkyl_array *fdo_host = app->use_gpu? mkarr(false, basis_do.num_basis, local_ext_do.volume)
                                            : gkyl_array_acquire(fdo);

  // Read donor distribution function and Jacobian inverse.
  struct gkyl_app_restart_status rstat;
  rstat.io_status = gkyl_comm_array_read(comm_do, &grid_do, &local_do, fdo_host, inp.file_name);

  bool scale_by_jacobtot = false;
  with_file(fp, inp.jacobtot_inv_file_name, "r") {
    // Configuration space donor grid and basis
    struct gkyl_rect_grid conf_grid_do;
    gkyl_rect_grid_init(&conf_grid_do, cdim_do, grid_do.lower, grid_do.upper, grid_do.cells);
    struct gkyl_basis conf_basis_do;
    gkyl_cart_modal_serendip(&conf_basis_do, cdim_do, poly_order);
    // Configuration space donor ranges and decomposition.
    struct gkyl_range conf_local_ext_do, conf_local_do, conf_global_ext_do, conf_global_do;
    gkyl_create_grid_ranges(&conf_grid_do, ghost_do, &conf_global_ext_do, &conf_global_do);
    struct gkyl_rect_decomp *conf_decomp_do = gkyl_rect_decomp_new_from_cuts(cdim_do, cuts_do, &conf_global_do);
    gkyl_create_ranges(&conf_decomp_do->ranges[my_rank], ghost_do, &conf_local_ext_do, &conf_local_do);
    // Read and multiply by the reciprocal of jacobtot.
    struct gkyl_array *jacobtot_inv_do_host = mkarr(false, conf_basis_do.num_basis, conf_local_ext_do.volume);
    rstat.io_status = gkyl_comm_array_read(comm_do, &conf_grid_do, &conf_local_do, jacobtot_inv_do_host, inp.jacobtot_inv_file_name);
    gkyl_dg_mul_conf_phase_op_range(&conf_basis_do, &basis_do, fdo_host, jacobtot_inv_do_host, fdo_host, &conf_local_ext_do, &local_ext_do);
    gkyl_array_release(jacobtot_inv_do_host);
    gkyl_rect_decomp_release(conf_decomp_do);
    scale_by_jacobtot = true;
  }

  bool scale_by_jacobvel = false;
  with_file(fp, inp.jacobvel_file_name, "r") {
    // Read and multiply by the reciprocal of jacobvel.
    struct gkyl_array *jacobvel_do_host = mkarr(false, 1, local_ext_do.volume);
    rstat.io_status = gkyl_comm_array_read(comm_do, &grid_do, &local_do, jacobvel_do_host, inp.jacobvel_file_name);
    gkyl_array_divide_by_cell(fdo_host, jacobvel_do_host);
    gkyl_array_release(jacobvel_do_host);
    scale_by_jacobvel = true;
  }

  if (app->use_gpu) {
    gkyl_array_copy(fdo, fdo_host);
  }

  if (pdim_do == pdim-1) {
    struct gkyl_translate_dim* transdim = gkyl_translate_dim_new(cdim_do,
      basis_do, cdim, gks->basis, -1, GKYL_NO_EDGE, app->use_gpu);
    gkyl_translate_dim_advance(transdim, &local_do, &gks->local, fdo, 1, gks->f);
    gkyl_translate_dim_release(transdim);
  }
  else {
    if (same_res) {
      gkyl_array_copy(gks->f, fdo);
    }
    else {
      // Interpolate the donor distribution to the target grid.
      struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(app->cdim, &gks->basis,
        &grid_do, &grid, &local_do, &gks->local, ghost_do, app->use_gpu);
      gkyl_dg_interpolate_advance(interp, fdo, gks->f);
      gkyl_dg_interpolate_release(interp);
    }
  }

  if (inp.type == GKYL_IC_IMPORT_AF) {
    // Scale f by a conf-space factor.
    gkyl_proj_on_basis *proj_conf_scale = gkyl_proj_on_basis_new(&app->grid, &app->basis,
      poly_order+1, 1, inp.conf_scale, inp.conf_scale_ctx);
    struct gkyl_array *xfac = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    struct gkyl_array *xfac_ho = app->use_gpu? mkarr(false, app->basis.num_basis, app->local_ext.volume)
                                             : gkyl_array_acquire(xfac);
    gkyl_proj_on_basis_advance(proj_conf_scale, 0.0, &app->local, xfac_ho);
    gkyl_array_copy(xfac, xfac_ho);
    gkyl_dg_mul_conf_phase_op_range(&app->basis, &gks->basis, gks->f, xfac, gks->f, &app->local, &gks->local);
    gkyl_proj_on_basis_release(proj_conf_scale);
    gkyl_array_release(xfac_ho);
    gkyl_array_release(xfac);
  }

  // Multiply f by the Jacobian.
  if (scale_by_jacobtot)
    gkyl_dg_mul_conf_phase_op_range(&app->basis, &gks->basis, gks->f, app->gk_geom->geo_int.jacobtot, gks->f, &app->local, &gks->local);

  // Multiply f by the velocity space Jacobian.
  if (scale_by_jacobvel)
    gkyl_array_scale_by_cell(gks->f, gks->vel_map->jacobvel);

  if (inp.enforce_positivity) {
    // Positivity enforcing by shifting f (ps=positivity shift).
    struct gkyl_positivity_shift_gyrokinetic *pos_shift_op = gkyl_positivity_shift_gyrokinetic_new(app->basis,
      gks->basis, gks->grid, gks->info.mass, app->gk_geom, gks->vel_map, &app->local_ext, app->use_gpu);

    gkyl_positivity_shift_gyrokinetic_advance(pos_shift_op, &app->local, &gks->local,
      gks->f, gks->m0.marr, gks->m0.marr);

    gkyl_positivity_shift_gyrokinetic_release(pos_shift_op);
  }

  gkyl_rect_decomp_release(decomp_do);
  gkyl_comm_release(comm_do);
  gkyl_array_release(fdo);
  gkyl_array_release(fdo_host);
}

static bool
gk_species_do_I_recycle(struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // Check whether one of the neutral species has recycling BCs that depend on
  // this gyrokinetic species.
  bool recycling_bcs = false;
  int neuts = app->num_neut_species;
  for (int i=0; i<neuts; ++i) {
    for (int k=0; k<2*app->cdim; k++) {
      const struct gkyl_gyrokinetic_bc *bc = &app->neut_species[i].info.bcs[k];
      if (bc->type == GKYL_BC_GK_SPECIES_RECYCLE) {
        for (int j=0; j<bc->emission.num_species; j++)
          recycling_bcs = recycling_bcs || 0 == strcmp(gks->info.name, bc->emission.in_species[j]);
      }
    }
  }
  return recycling_bcs;
}

static bool
gk_species_do_I_recycle_react_scale(struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // Check whether one of the neutral species has a recycle_react_scale
  // operation thats depend on this gyrokinetic species.
  bool has_sca = false;
  int neuts = app->num_neut_species;
  for (int i=0; i<neuts; ++i) {
    struct gk_neut_species *ns = &app->neut_species[i];
    struct gkyl_gyrokinetic_scaling_inp *sca_inp = &ns->info.scaling;
    if ((sca_inp->num_boundaries > 0) && (0 == strcmp(gks->info.name, sca_inp->impacting_ion_name))) {
      has_sca = true;
      break;
    }
  }
  return has_sca;
}

static bool
gk_species_do_I_adapt_src(struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // Check whether one of the species adapts its source depending on
  // this gyrokinetic species.
  bool adapt_src = false;
  for (int i=0; i<app->num_species; ++i) {
    struct gkyl_gyrokinetic_source *source = &app->species[i].info.source;
    if (source->num_adapt_sources > 0) {
      for (int j=0; j<source->num_adapt_sources; j++)
        adapt_src = adapt_src || 0 == strcmp(gks->info.name, source->adapt[j].adapt_to_species);
    }
  }
  return adapt_src;
}

void
gk_species_init(struct gkyl_gk *gk_app_inp, struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  int cdim = app->cdim, vdim = gks->info.vdim;
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  assert(vdim > 0); // Ensure user provided vdim in input file.

  for (int d=0; d<cdim; ++d) {
    cells[d] = gk_app_inp->cells[d];
    lower[d] = gk_app_inp->lower[d];
    upper[d] = gk_app_inp->upper[d];
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    cells[cdim+d] = gks->info.cells[d];
    lower[cdim+d] = gks->info.lower[d];
    upper[cdim+d] = gks->info.upper[d];
    ghost[cdim+d] = 0; // No ghost-cells in velocity space.

    cells_vel[d] = gks->info.cells[d];
    lower_vel[d] = gks->info.lower[d];
    upper_vel[d] = gks->info.upper[d];
    ghost_vel[d] = 0; // No ghost-cells in velocity space.
  }

  // Define phase basis
  if (app->use_gpu) {
    // allocate device basis if we are using GPUs
    gks->basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  }
  else {
    gks->basis_on_dev = &gks->basis;
  }

  enum gkyl_basis_type b_type = app->basis.b_type;

  // basis functions
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (app->poly_order > 1) {
        gkyl_cart_modal_serendip(&gks->basis, pdim, app->poly_order);
      }
      else if (app->poly_order == 1) {
        gkyl_cart_modal_gkhybrid(&gks->basis, cdim, vdim); // p=2 in vparallel
      }

      if (app->use_gpu) {
        if (app->poly_order > 1) {
          gkyl_cart_modal_serendip_cu_dev(gks->basis_on_dev, pdim, app->poly_order);
	}
        else if (app->poly_order == 1) {
          gkyl_cart_modal_gkhybrid_cu_dev(gks->basis_on_dev, cdim, vdim); // p=2 in vparallel
        }
      }
      break;
    default:
      assert(false);
      break;
  }
  
  // Full phase space grid.
  gkyl_rect_grid_init(&gks->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&gks->grid, ghost, &gks->global_ext, &gks->global);
  
  // Velocity space grid.
  gkyl_rect_grid_init(&gks->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&gks->grid_vel, ghost_vel, &gks->local_ext_vel, &gks->local_vel);

  // Phase-space communicator.
  gks->comm = gkyl_comm_extend_comm(app->comm, &gks->local_vel);

  // Create local and local_ext from app local range.
  struct gkyl_range local;
  // Local = conf-local X local_vel.
  gkyl_range_ten_prod(&local, &app->local, &gks->local_vel);
  gkyl_create_ranges(&local, ghost, &gks->local_ext, &gks->local);

  // Velocity space mapping.
  gks->vel_map = gkyl_velocity_map_new(gks->info.mapc2p, gks->grid, gks->grid_vel,
    gks->local, gks->local_ext, gks->local_vel, gks->local_ext_vel, app->use_gpu);

  // Write out the velocity space mapping and its Jacobian.
  gkyl_velocity_map_write(gks->vel_map, gks->comm, app->name, gks->info.name);
  
  // Keep a copy of num_periodic_dir and periodic_dirs in species so we can
  // modify it in GK_IWL BCs without modifying the app's.
  gks->num_periodic_dir = app->num_periodic_dir;
  for (int d=0; d<gks->num_periodic_dir; ++d)
    gks->periodic_dirs[d] = app->periodic_dirs[d];

  for (int d=0; d<app->cdim; ++d) gks->bc_is_np[d] = true;
  for (int d=0; d<gks->num_periodic_dir; ++d)
    gks->bc_is_np[gks->periodic_dirs[d]] = false;

  // Store the BCs from the input file.
  for (int d=0; d<app->cdim; ++d) {
    if (gks->bc_is_np[d]) {
      struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(gks->info.bcs, 2*app->cdim, d, GKYL_LOWER_EDGE);
      if (bc_lo != 0)
        gks->lower_bc[d] = *bc_lo;
      else
        gks->lower_bc[d].type = GKYL_BC_GK_SKIP;
  
      struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(gks->info.bcs, 2*app->cdim, d, GKYL_UPPER_EDGE);
      if (bc_up != 0)
        gks->upper_bc[d] = *bc_up;
      else
        gks->upper_bc[d].type = GKYL_BC_GK_SKIP;
    }
    else {
      gks->lower_bc[d].type = GKYL_BC_GK_SPECIES_PERIODIC;
      gks->upper_bc[d].type = GKYL_BC_GK_SPECIES_PERIODIC;
    }
  }
 
  // Metadata for gk_species app.
  struct gkyl_msgpack_map_elem io_meta[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = gks->basis.poly_order },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = gks->basis.id }
  };
  gks->io_meta_len = sizeof(io_meta)/sizeof(io_meta[0]);
  gks->io_meta = gkyl_msgpack_map_elem_clone(gks->io_meta_len, io_meta);

  // Allocate distribution function arrays.
  gks->f = mkarr(app->use_gpu, gks->basis.num_basis, gks->local_ext.volume);

  gks->f_host = gks->f;
  if (app->use_gpu)
    gks->f_host = mkarr(false, gks->basis.num_basis, gks->local_ext.volume);

  // Allocate cflrate (scalar array).
  gks->cflrate = mkarr(app->use_gpu, 1, gks->local_ext.volume);

  // Allocate data for density (for charge density or upar calculation).
  gk_species_moment_init(app, gks, &gks->m0, GKYL_F_MOMENT_M0, false);

  if (gks->info.flr.type) {
    // Create operator needed for FLR effects.
    assert(app->cdim > 1);
    // Pointer to function performing the gyroaverage.
    gks->gyroaverage = gk_species_gyroaverage_enabled;
    // Gyroaveraged M0 and phi.
    gks->m0_gyroavg = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gks->gyro_phi = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    double gyroradius_bmag = gks->info.flr.bmag ? gks->info.flr.bmag : app->bmag_ref;

    double flr_weight = gks->info.flr.Tperp*gks->info.mass/(pow(gks->info.charge*gyroradius_bmag,2.0));
    // Initialize the weight in the Laplacian operator.
    gks->flr_rhoSqD2 = mkarr(app->use_gpu, (2*(app->cdim-1)-1)*app->basis.num_basis, app->local_ext.volume);
    gkyl_array_set_offset(gks->flr_rhoSqD2, flr_weight, app->gk_geom->geo_int.gxxj, 0*app->basis.num_basis);
    if (app->cdim > 2) {
      gkyl_array_set_offset(gks->flr_rhoSqD2, flr_weight, app->gk_geom->geo_int.gxyj, 1*app->basis.num_basis);
      gkyl_array_set_offset(gks->flr_rhoSqD2, flr_weight, app->gk_geom->geo_int.gyyj, 2*app->basis.num_basis);
    }
    // Initialize the factor multiplying the field in the FLR operator.
    gks->flr_kSq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkyl_array_shiftc(gks->flr_kSq, -pow(sqrt(2.0),app->cdim), 0); // Sets kSq=-1.

    // If domain is not periodic use Dirichlet BCs.
    struct gkyl_poisson_bc flr_bc;
    for (int d=0; d<app->cdim-1; d++) {
      struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(app->field->info.poisson_bcs, 2*app->cdim, d, GKYL_LOWER_EDGE);
      if (bc_lo != 0) {
        if (bc_lo->type == GKYL_BC_GK_FIELD_PERIODIC)
          flr_bc.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
        else
          flr_bc.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_DIRICHLET_VARYING);
      }

      struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(app->field->info.poisson_bcs, 2*app->cdim, d, GKYL_UPPER_EDGE);
      if (bc_up != 0) {
        if (bc_up->type == GKYL_BC_GK_FIELD_PERIODIC)
          flr_bc.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
        else
          flr_bc.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_DIRICHLET_VARYING);
      }
    }
    // Deflated Poisson solve is performed on range assuming decomposition is *only* in z.
    gks->flr_op = gkyl_deflated_fem_poisson_new(app->grid, app->basis_on_dev, app->basis,
      app->local, app->local, gks->flr_rhoSqD2, gks->flr_kSq, flr_bc, NULL, app->use_gpu);
  }
  else {
    gks->gyroaverage = gk_species_gyroaverage_disabled;
    gks->m0_gyroavg = gkyl_array_acquire(gks->m0.marr);
    gks->gyro_phi = gkyl_array_acquire(app->field->phi_smooth);
  }

  // Initialize the collisionless solver.
  gks->collisionless = (struct gk_collisionless) { };
  gk_species_collisionless_init(app, gks, &gks->collisionless);

  // Initialize an anomalous diffusion term.
  gks->anom_diff = (struct gk_anomalous_diff) { };
  gk_species_anomalous_diff_init(app, gks, &gks->anom_diff);

  // Damping term -nu*f on RHS.
  gks->damping = (struct gk_damping) { };
  gk_species_damping_init(app, gks, &gks->damping);

  // Function multiplying df/dt.
  gk_species_fdot_multiplier_init(app, gks, &gks->fdot_mult);

  // Allocate data for diagnostic moments.
  int ndm = gks->info.num_diag_moments;
  gks->moms = gkyl_malloc(sizeof(struct gk_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    gk_species_moment_init(app, gks, &gks->moms[m], gks->info.diag_moments[m], false);

  // Initialize projection routine for initial conditions.
  if (gks->info.init_from_file.type == 0) {
    gk_species_projection_init(app, gks, gks->info.projection, &gks->proj_init);
  }
  else {
    gk_species_file_import_init(app, gks, gks->info.init_from_file);
  }

  // Create skin/ghost ranges.
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&gks->local_lower_skin[dir], &gks->local_lower_ghost[dir],
      dir, GKYL_LOWER_EDGE, &gks->local_ext, ghost);
    gkyl_skin_ghost_ranges(&gks->local_upper_skin[dir], &gks->local_upper_ghost[dir],
      dir, GKYL_UPPER_EDGE, &gks->local_ext, ghost);
    gkyl_skin_ghost_ranges(&gks->global_lower_skin[dir], &gks->global_lower_ghost[dir],
      dir, GKYL_LOWER_EDGE, &gks->global_ext, ghost); 
    gkyl_skin_ghost_ranges(&gks->global_upper_skin[dir], &gks->global_upper_ghost[dir],
      dir, GKYL_UPPER_EDGE, &gks->global_ext, ghost);
  }

  if (gk_app_inp->geometry.has_LCFS) {
    // Simulation spans the last-closed flux surface (LCFS). Create core and SOL global ranges.
    int idx_LCFS_lo = app->gk_geom->idx_LCFS_lo;
    // Length of lower and upper x ranges (one is core, the other SOL).
    int len_lo = idx_LCFS_lo;
    int len_up = gks->global.upper[0]-len_lo;
    // Lower and upper x ranges.
    struct gkyl_range *global_lo_r, *global_ext_lo_r, *global_up_r, *global_ext_up_r;
    struct gkyl_range *local_lo_r, *local_ext_lo_r, *local_up_r, *local_ext_up_r;
    struct gkyl_range *local_lower_skin_par_lo_r, *local_upper_skin_par_lo_r,
                      *local_lower_ghost_par_lo_r, *local_upper_ghost_par_lo_r;
    struct gkyl_range *local_lower_skin_par_up_r, *local_upper_skin_par_up_r,
                      *local_lower_ghost_par_up_r, *local_upper_ghost_par_up_r;
    if (app->gk_geom->geqdsk_sign_convention == 0) {
      // x increases towards SOL.
      global_lo_r                = &gks->global_core;
      global_up_r                = &gks->global_sol;
      global_ext_lo_r            = &gks->global_ext_core;
      global_ext_up_r            = &gks->global_ext_sol;
      local_lo_r                 = &gks->local_core;
      local_up_r                 = &gks->local_sol;
      local_ext_lo_r             = &gks->local_ext_core;
      local_ext_up_r             = &gks->local_ext_sol;
      local_lower_skin_par_lo_r  = &gks->local_lower_skin_par_core;
      local_lower_ghost_par_lo_r = &gks->local_lower_ghost_par_core;
      local_lower_skin_par_up_r  = &gks->local_lower_skin_par_sol;
      local_lower_ghost_par_up_r = &gks->local_lower_ghost_par_sol;
      local_upper_skin_par_lo_r  = &gks->local_upper_skin_par_core;
      local_upper_ghost_par_lo_r = &gks->local_upper_ghost_par_core;
      local_upper_skin_par_up_r  = &gks->local_upper_skin_par_sol;
      local_upper_ghost_par_up_r = &gks->local_upper_ghost_par_sol;
    }
    else {
      // x increases towards core.
      global_lo_r                = &gks->global_sol;
      global_ext_lo_r            = &gks->global_ext_sol;
      global_up_r                = &gks->global_core;
      global_ext_up_r            = &gks->global_ext_core;
      local_lo_r                 = &gks->local_sol;
      local_up_r                 = &gks->local_core;
      local_ext_lo_r             = &gks->local_ext_sol;
      local_ext_up_r             = &gks->local_ext_core;
      local_lower_skin_par_lo_r  = &gks->local_lower_skin_par_sol;
      local_upper_skin_par_lo_r  = &gks->local_upper_skin_par_sol;
      local_lower_ghost_par_lo_r = &gks->local_lower_ghost_par_sol;
      local_upper_ghost_par_lo_r = &gks->local_upper_ghost_par_sol;
      local_lower_skin_par_up_r  = &gks->local_lower_skin_par_core;
      local_upper_skin_par_up_r  = &gks->local_upper_skin_par_core;
      local_lower_ghost_par_up_r = &gks->local_lower_ghost_par_core;
      local_upper_ghost_par_up_r = &gks->local_upper_ghost_par_core;
    }

    // Lower and upper x ranges.
    gkyl_range_shorten_from_above(global_lo_r, &gks->global, 0, len_lo);
    gkyl_range_shorten_from_below(global_up_r, &gks->global, 0, len_up);
    gkyl_range_shorten_from_above(local_lo_r, &gks->local, 0, len_lo);
    gkyl_range_shorten_from_below(local_up_r, &gks->local, 0, len_up);

    // Extended lower and upper x ranges.
    int len_lo_ext = idx_LCFS_lo+1;
    int len_up_ext = gks->global_ext.upper[0]-len_lo;
    gkyl_range_shorten_from_above(global_ext_lo_r, &gks->global_ext, 0, len_lo_ext);
    gkyl_range_shorten_from_below(global_ext_up_r, &gks->global_ext, 0, len_up_ext);
    gkyl_range_shorten_from_above(local_ext_lo_r, &gks->local_ext, 0, len_lo_ext);
    gkyl_range_shorten_from_below(local_ext_up_r, &gks->local_ext, 0, len_up_ext);

    // Parallel skin and ghost ranges, limited to the lower and upper x range.
    int par_dir = app->cdim-1;
    for (int e=0; e<2; e++) {
      gkyl_range_shorten_from_above(e==0? local_lower_skin_par_lo_r        : local_upper_skin_par_lo_r,
                                    e==0? &gks->local_lower_skin[par_dir]  : &gks->local_upper_skin[par_dir], 0, len_lo);
      gkyl_range_shorten_from_above(e==0? local_lower_ghost_par_lo_r       : local_upper_ghost_par_lo_r,
                                    e==0? &gks->local_lower_ghost[par_dir] : &gks->local_upper_ghost[par_dir], 0, len_lo);
      gkyl_range_shorten_from_below(e==0? local_lower_skin_par_up_r        : local_upper_skin_par_up_r,
                                    e==0? &gks->local_lower_skin[par_dir]  : &gks->local_upper_skin[par_dir], 0, len_up);
      gkyl_range_shorten_from_below(e==0? local_lower_ghost_par_up_r       : local_upper_ghost_par_up_r,
                                    e==0? &gks->local_lower_ghost[par_dir] : &gks->local_upper_ghost[par_dir], 0, len_up);
    }

    // Core range extended in the parallel direction.
    int ndim = gks->local.ndim;
    int lower_bcdir_ext[ndim], upper_bcdir_ext[ndim];
    for (int i=0; i<ndim; i++) {
      lower_bcdir_ext[i] = gks->local_core.lower[i];
      upper_bcdir_ext[i] = gks->local_core.upper[i];
    }
    lower_bcdir_ext[par_dir] = gks->local_ext_core.lower[par_dir];
    upper_bcdir_ext[par_dir] = gks->local_ext_core.upper[par_dir];
    gkyl_sub_range_init(&gks->local_par_ext_core, &gks->local_ext_core, lower_bcdir_ext, upper_bcdir_ext);
  }

  // Initialize boundary fluxes.
  gks->bflux = (struct gk_boundary_fluxes) { };
  // Additional bflux moments to step in time.
  struct gkyl_phase_diagnostics_inp add_bflux_moms_inp = (struct gkyl_phase_diagnostics_inp) { };
  enum gkyl_species_bflux_type bflux_type = GK_SPECIES_BFLUX_NONE;
  // Check if using Boltzmann elc.
  bool ion_in_boltz_elc_field = (gks->info.charge > 0.0) && 
    app->field->update_field && (app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN);
  // Check if sources are adaptive.
  bool adaptive_sources = gk_species_do_I_adapt_src(app, gks);
  // Check if other species use the recycle_react_scale operation.
  bool recycle_react_scale = gk_species_do_I_recycle_react_scale(app, gks);
  // Check if other species have recycling BCs.
  bool recycling_bcs = gk_species_do_I_recycle(app, gks);
  if (gks->info.boundary_flux_diagnostics.num_diag_moments > 0 ||
      gks->info.boundary_flux_diagnostics.num_integrated_diag_moments > 0) {
    bflux_type = GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS_DIAGS;
  }
  else {
    // Set bflux_type to 
    //   - GK_SPECIES_BFLUX_CALC_FLUX to only put bfluxes in ghost cells of rhs.
    //   - GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS to calc bfluxes and step its moments.
    // The latter also requires that you place the moment you desire in add_bflux_moms_inp below.
    
    if (recycling_bcs) {
      bflux_type = GK_SPECIES_BFLUX_CALC_FLUX;
    }
    if (ion_in_boltz_elc_field || adaptive_sources || recycle_react_scale) {
      // This is within an if-statement instead of an else if because it
      // superseeds (and is a superset) of GK_SPECIES_BFLUX_CALC_FLUX.
      bflux_type = GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS;
    }
  }
  int *nmom_extra = &add_bflux_moms_inp.num_diag_moments;
  if (ion_in_boltz_elc_field || adaptive_sources || recycle_react_scale) {
    add_bflux_moms_inp.diag_moments[nmom_extra[0]++] = GKYL_F_MOMENT_M0;
  }
  if (adaptive_sources) {
    add_bflux_moms_inp.diag_moments[nmom_extra[0]++] = GKYL_F_MOMENT_M2;
  }
  // Introduce new moments into moms_inp if needed.
  gk_species_bflux_init(app, gks, &gks->bflux, bflux_type, add_bflux_moms_inp);
  
  // Initialize empty structs. New methods will fill them if specified.
  gks->src = (struct gk_source) { };
  gks->rad = (struct gk_rad_drag) { };

  // Initialize a Maxwellian/LTE (local thermodynamic equilibrium) projection routine
  // Projection routine optionally corrects all the Maxwellian/LTE moments
  // This routine is utilized by both reactions and BGK collisions
  gks->lte = (struct gk_lte) { };
  bool correct_all_moms = gks->info.correct.correct_all_moms;
  int max_iter = gks->info.correct.max_iter > 0 ? gks->info.correct.max_iter : 50;
  double iter_eps = gks->info.correct.iter_eps > 0 ? gks->info.correct.iter_eps : 1e-10;
  bool use_last_converged = gks->info.correct.use_last_converged;
  struct correct_all_moms_inp corr_inp = {
    .correct_all_moms = correct_all_moms, 
    .max_iter = max_iter,
    .iter_eps = iter_eps, 
    .use_last_converged = use_last_converged
  };
  gk_species_lte_init(app, gks, &gks->lte, corr_inp);

  // Initialize the object that scales the species (meant for fluid neutrals for now).
  gks->sca = (struct gk_scaling) { };
  gk_species_scaling_init(app, gks, &gks->sca);

  // Initialize reactions with charged species.
  gks->react = (struct gk_react) { };
  gk_species_react_init(app, gks, gks->info.react, &gks->react, true);

  // Initialize reactions with neutral species.
  gks->react_neut = (struct gk_react) { };
  gk_species_react_init(app, gks, gks->info.react_neut, &gks->react_neut, false);

  // Initialize LBO collisions.
  gks->lbo = (struct gk_lbo_collisions) { };
  gk_species_lbo_init(app, gks, &gks->lbo);

  // Initialize BGK collisions.
  gks->bgk = (struct gk_bgk_collisions) { };
  gk_species_bgk_init(app, gks, &gks->bgk);

  // Initialize a BGK source.
  gks->bgk_src = (struct gk_source_bgk) { };
  gk_species_source_bgk_init(app, gks, &gks->bgk_src);

  // Initialize positivity enforcing operator.
  gks->positivity = (struct gk_positivity) { };
  gk_species_positivity_init(app, gks, &gks->positivity);

  if (gks->info.is_static) {
    gk_species_init_static(gk_app_inp, app, gks);
  }
  else {
    gk_species_init_dynamic(gk_app_inp, app, gks);
  }
}

void
gk_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_species *gks, double t0)
{
  if (gks->info.init_from_file.type == 0)
    gk_species_projection_calc(app, gks, &gks->proj_init, gks->f, t0);

  // Initialize fbar for low-pass filter damping
  gk_species_damping_set_fbar_to_f(gks, &gks->damping, gks->f);

  // We are pre-computing source for now as it is time-independent.
  gk_species_source_calc(app, gks, &gks->src, gks->lte.f_lte, t0);
}

void
gk_species_apply_ic_cross(gkyl_gyrokinetic_app *app, struct gk_species *gks_self, double t0)
{
  // IC setup step that depends on the IC of other species.
  if (app->field->init_phi_pol && gks_self->info.scale_with_polarization) {
    // Scale the distribution function so its guiding center density is computed
    // given the polarization density and the guiding center density of other species.

    // Compute the polarization density.
    struct gkyl_array *npol = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    struct gkyl_gyrokinetic_pol_density* npol_op = gkyl_gyrokinetic_pol_density_new(app->basis, app->grid, app->use_gpu);
    gkyl_gyrokinetic_pol_density_advance(npol_op, &app->local, app->field->epsilon, app->field->phi_pol, npol);
    gkyl_gyrokinetic_pol_density_release(npol_op);

    // Calculate the guiding center density of this species: (npol - q_other*n^G_other)/q_self.
    for (int i=0; i<app->num_species; ++i) {
      struct gk_species *gks = &app->species[i];
      gk_species_moment_calc(&gks->m0, gks->local, app->local, gks->f);
      if (strcmp(gks->info.name, gks_self->info.name)) {
        gkyl_array_accumulate(npol, -gks->info.charge, gks->m0.marr);
      }
    }
    gkyl_array_scale(npol, 1./gks_self->info.charge);

    // Scale the distribution function so it has this guiding center density.
    struct gkyl_dg_bin_op_mem *div_mem = app->use_gpu? gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis)
                                                     : gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);
    struct gkyl_array *den_mod = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    gkyl_dg_div_op_range(div_mem, app->basis,
      0, den_mod, 0, npol, 0, gks_self->m0.marr, &app->local);
    gkyl_dg_mul_conf_phase_op_range(&app->basis, &gks_self->basis, gks_self->f,
      den_mod, gks_self->f, &app->local_ext, &gks_self->local_ext);

    gkyl_array_release(den_mod);
    gkyl_dg_bin_op_mem_release(div_mem);
    gkyl_array_release(npol);
  }

  // Store initial density in scaling operator.
  gk_species_scaling_apply_ic_cross(app, gks_self, &gks_self->sca);
}

double
gk_species_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms)
{
  return species->rhs_func(app, species, fin, fbar_in, rhs, bflux_moms);
}

double
gk_species_rhs_implicit(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt)
{
  return species->rhs_implicit_func(app, species, fin, rhs, bflux_moms, dt);
}

void
gk_species_step_f(struct gk_species *species, struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  species->step_f_func(out, a, inp);
}

void
gk_species_combine(struct gk_species *species, struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  species->combine_func(out, c1, arr1, c2, arr2, rng);
}

void
gk_species_copy_range(struct gk_species *species, struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  species->copy_func(out, inp, range);
}

void
gk_species_apply_bc(gkyl_gyrokinetic_app *app, const struct gk_species *species, struct gkyl_array *f)
{
  species->bc_func(app, species, f);
}

void
gk_species_n_iter_corr(gkyl_gyrokinetic_app *app, const struct gk_species *s, int sidx)
{
  app->stat.num_corr[sidx] = s->lte.num_corr;
  app->stat.n_iter_corr[sidx] = s->lte.n_iter;
}

// write functions
void
gk_species_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  if (frame == 0) {
    gk_species_write_dynamic(app, gks, tm, frame);
  }
  else
    gks->write_func(app, gks, tm, frame);
}

void
gk_species_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  if (frame == 0) {
    gk_species_write_mom_dynamic(app, gks, tm, frame);
  }
  else
    gks->write_mom_func(app, gks, tm, frame);
}

void
gk_species_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  gks->calc_integrated_mom_func(app, gks, tm);
}

void
gk_species_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  gks->write_integrated_mom_func(app, gks);
}

void
gk_species_calc_L2norm(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  gks->calc_L2norm_func(app, gks, tm);
}

void
gk_species_write_L2norm(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  gks->write_L2norm_func(app, gks);
}

void
gk_species_release(const gkyl_gyrokinetic_app* app, const struct gk_species *gks)
{
  // Release resources for charged species.
  gkyl_array_release(gks->f);

  gkyl_array_release(gks->cflrate);

  if (gks->info.init_from_file.type == 0) {
    gk_species_projection_release(app, &gks->proj_init);
  }

  gkyl_comm_release(gks->comm);

  gkyl_msgpack_map_elem_release(gks->io_meta_len, gks->io_meta);

  if (app->use_gpu) {
    gkyl_array_release(gks->f_host);
    gkyl_cu_free(gks->basis_on_dev);
  }

  gkyl_velocity_map_release(gks->vel_map);

  gk_species_collisionless_release(app, &gks->collisionless);

  gk_species_scaling_release(app, &gks->sca);

  gk_species_anomalous_diff_release(app, &gks->anom_diff);

  // Release moment data.
  gk_species_moment_release(app, &gks->m0);
  for (int i=0; i<gks->info.num_diag_moments; ++i) {
    gk_species_moment_release(app, &gks->moms[i]);
  }
  gkyl_free(gks->moms);

  gk_species_source_release(app, &gks->src);

  gk_species_damping_release(app, &gks->damping);

  gk_species_fdot_multiplier_release(app, &gks->fdot_mult);

  gk_species_lbo_release(app, &gks->lbo);

  gk_species_bgk_release(app, &gks->bgk);

  gk_species_source_bgk_release(app, &gks->bgk_src);

  gk_species_positivity_release(app, &gks->positivity);

  gk_species_radiation_release(app, &gks->rad);

  gk_species_react_release(app, &gks->react);

  gk_species_react_release(app, &gks->react_neut);  

  // Free boundary flux memory.
  gk_species_bflux_release(app, gks, &gks->bflux);
  
  gk_species_lte_release(app, &gks->lte);

  gkyl_array_release(gks->m0_gyroavg);
  gkyl_array_release(gks->gyro_phi);
  if (gks->info.flr.type) {
    gkyl_array_release(gks->flr_rhoSqD2);
    gkyl_array_release(gks->flr_kSq);
    gkyl_deflated_fem_poisson_release(gks->flr_op);
  }

  gks->release_func(app, gks);
}
