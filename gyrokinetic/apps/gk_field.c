#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_position_map.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_comm_io.h>

#include <assert.h>
#include <float.h>
#include <time.h>

// Functions related to setting the potential by adjusting the polarization density
static void
eval_on_nodes_c2p_position_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gkyl_position_map *gpm = ctx;
  gkyl_position_map_eval_mc2nu(gpm, xcomp, xphys);
}

static struct gkyl_app_restart_status
header_from_file(gkyl_gyrokinetic_app *app, const char *fname)
{
  struct gkyl_app_restart_status rstat = { .io_status = GKYL_ARRAY_RIO_FOPEN_FAILED };
  
  FILE *fp = 0;
  with_file(fp, fname, "r") {
    struct gkyl_rect_grid grid;
    struct gkyl_array_header_info hdr;
    rstat.io_status = gkyl_grid_sub_array_header_read_fp(&grid, &hdr, fp);

    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      if (hdr.etype != GKYL_DOUBLE)
        rstat.io_status = GKYL_ARRAY_RIO_DATA_MISMATCH;
    }

    struct gkyl_msgpack_map_elem elem_list[] = {
      { .key = "frame", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
      { .key = "time", .elem_type = GKYL_MP_DOUBLE, .cval = 0 },
    };
    int elem_list_len = sizeof(elem_list)/sizeof(elem_list[0]);
    gkyl_msgpack_to_map_elem_list(&(struct gkyl_msgpack_data) {
        .meta = hdr.meta,
        .meta_sz = hdr.meta_size
      }, elem_list_len, elem_list);

    rstat.frame = gkyl_msgpack_map_elem_get_uint(elem_list_len, elem_list, "frame");
    rstat.stime = gkyl_msgpack_map_elem_get_double(elem_list_len, elem_list, "time");

    gkyl_grid_sub_array_header_release(&hdr);
  }
  
  return rstat;
}

static void
gk_field_polarization_potential_new(struct gk_field *f, struct gkyl_gyrokinetic_app *app)
{
  // Project the initial potential onto a p+1 tensor basis and compute the polarization
  // density to use use by species in calculating the initial ion density.
  f->init_phi_pol = true;
  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, app->cdim, app->poly_order + 1);

  f->phi_pol = mkarr(app->use_gpu, phi_pol_basis.num_basis, app->local_ext.volume);
  struct gkyl_array *phi_pol_ho = app->use_gpu ? mkarr(false, f->phi_pol->ncomp, f->phi_pol->size)
                                               : gkyl_array_acquire(f->phi_pol);

  struct gkyl_eval_on_nodes *phi_pol_proj = gkyl_eval_on_nodes_inew(&(struct gkyl_eval_on_nodes_inp){
      .grid = &app->grid,
      .basis = &phi_pol_basis,
      .num_ret_vals = 1,
      .eval = f->info.polarization_potential,
      .ctx = f->info.polarization_potential_ctx,
      .c2p_func = eval_on_nodes_c2p_position_func,
      .c2p_func_ctx = app->position_map,
  });

  gkyl_eval_on_nodes_advance(phi_pol_proj, 0.0, &app->local, phi_pol_ho);
  gkyl_array_copy(f->phi_pol, phi_pol_ho);

  gkyl_eval_on_nodes_release(phi_pol_proj);
  gkyl_array_release(phi_pol_ho);
}

static void
gk_field_polarization_potential_from_file_new(struct gk_field *f, struct gkyl_gyrokinetic_app *app,
  struct gkyl_gyrokinetic_ic_import inp)
{
  f->init_phi_pol = true;
  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, app->cdim, app->poly_order + 1);

  f->phi_pol = mkarr(app->use_gpu, phi_pol_basis.num_basis, app->local_ext.volume);
  struct gkyl_array *phi_pol_ho = app->use_gpu ? mkarr(false, f->phi_pol->ncomp, f->phi_pol->size)
                                               : gkyl_array_acquire(f->phi_pol);

  struct gkyl_app_restart_status rstat = header_from_file(app, inp.file_name);
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status = gkyl_comm_array_read(app->comm, &app->grid, &app->local, phi_pol_ho, inp.file_name);
    assert(rstat.io_status == GKYL_ARRAY_RIO_SUCCESS);
    gkyl_array_copy(f->phi_pol, phi_pol_ho);
  }
  else {
    assert(false);
  }

  gkyl_array_release(phi_pol_ho);
}

static void
gk_field_polarization_potential_release(struct gk_field *f)
{
  gkyl_array_release(f->phi_pol);
}

// Functions related to the field energy allocations, diagnostics, and release

static void
gk_field_calc_energy_dt_active(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_integrate_advance(field->calc_em_energy, field->phi_smooth, 
    1.0/dt, field->es_energy_fac, &app->local, &app->local, energy_reduced);
  app->stat.phidot_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_field_calc_apar_energy_dt_active(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_integrate_advance(field->calc_em_energy, field->apar, 
    1.0/dt, field->apar_energy_fac, &app->local, &app->local, energy_reduced);
  app->stat.phidot_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_field_calc_energy_dt_none(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced)
{
  // Do nothing.
}

static void
gk_field_calc_energy_enabled(struct gkyl_gyrokinetic_app *app, const struct gk_field *field, double tm)
{
  gkyl_array_integrate_advance(field->calc_em_energy, field->phi_smooth,
                               1.0, field->es_energy_fac, &app->local, &app->local, field->em_energy_red);

  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red, field->em_energy_red_global);

  double energy_global[1] = {0.0};
  if (app->use_gpu) {
    gkyl_cu_memcpy(energy_global, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  } else {
    energy_global[0] = field->em_energy_red_global[0];
  }

  if (app->cdim == 1) {
    energy_global[0] *= field->es_energy_fac_1d;
  }

  gkyl_dynvec_append(field->integ_energy, tm, energy_global);

  if (field->info.time_rate_diagnostics) {
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red_old, field->em_energy_red_global);
    double energy_dot_global_old[1] = {0.0};
    if (app->use_gpu) {
      gkyl_cu_memcpy(energy_dot_global_old, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    } else {
      energy_dot_global_old[0] = field->em_energy_red_global[0];
    }
    if (app->cdim == 1) {
      energy_dot_global_old[0] *= field->es_energy_fac_1d;
    }

    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red_new, field->em_energy_red_global);
    double energy_dot_global_new[1] = {0.0};
    if (app->use_gpu) {
      gkyl_cu_memcpy(energy_dot_global_new, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    } else {
      energy_dot_global_new[0] = field->em_energy_red_global[0];
    }
    if (app->cdim == 1) {
      energy_dot_global_new[0] *= field->es_energy_fac_1d;
    }

    double energy_dot_global[1] = {0.0};
    energy_dot_global[0] = energy_dot_global_new[0] - energy_dot_global_old[0];

    gkyl_dynvec_append(field->integ_energy_dot, tm, energy_dot_global);
  }

  if (field->is_em) {
    // Calculate Aparallel energy.
    gkyl_array_integrate_advance(field->calc_em_energy, field->apar, 
      1.0, field->apar_energy_fac, &app->local, &app->local, field->em_energy_red);

    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red, field->em_energy_red_global);

    double energy_global[1] = { 0.0 };
    if (app->use_gpu)
      gkyl_cu_memcpy(energy_global, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    else
      energy_global[0] = field->em_energy_red_global[0];

    if (app->cdim == 1)
      energy_global[0] *= field->apar_energy_fac_1d; 

    gkyl_dynvec_append(field->integ_apar_energy, tm, energy_global);
    
    // Calculate d(Aparallel)/dt energy.
    gkyl_array_integrate_advance(field->calc_em_energy, field->apardot, 
      1.0, field->apar_energy_fac, &app->local, &app->local, field->em_energy_red);

    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red, field->em_energy_red_global);

    double energy_dot_global[1] = { 0.0 };
    if (app->use_gpu)
      gkyl_cu_memcpy(energy_dot_global, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    else
      energy_dot_global[0] = field->em_energy_red_global[0];
    
    // if (app->cdim == 1)
    //   energy_dot_global[0] *= field->apar_energy_fac_1d;

    gkyl_dynvec_append(field->integ_apardot_energy, tm, energy_dot_global);

    if (field->info.time_rate_diagnostics) {
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->apar_energy_red_old, field->em_energy_red_global);
      double energy_dot_global_old[1] = { 0.0 };
      if (app->use_gpu)
        gkyl_cu_memcpy(energy_dot_global_old, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
      else
        energy_dot_global_old[0] = field->em_energy_red_global[0];
      if (app->cdim == 1)
        energy_dot_global_old[0] *= field->apar_energy_fac_1d; 
      
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->apar_energy_red_new, field->em_energy_red_global);
      double energy_dot_global_new[1] = { 0.0 };
      if (app->use_gpu)
        gkyl_cu_memcpy(energy_dot_global_new, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
      else
        energy_dot_global_new[0] = field->em_energy_red_global[0];
      if (app->cdim == 1)
        energy_dot_global_new[0] *= field->apar_energy_fac_1d; 

      double energy_dot_global[1] = { 0.0 };
      energy_dot_global[0] = energy_dot_global_new[0] - energy_dot_global_old[0];

      gkyl_dynvec_append(field->integ_apar_energy_dot, tm, energy_dot_global);
    }
  }
}

static void
gk_field_calc_energy_disabled(struct gkyl_gyrokinetic_app *app, const struct gk_field *field, double tm)
{
  // Do nothing.
}

static void
gk_field_combine_enabled(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  gkyl_array_accumulate_range(gkyl_array_set_range(out, c1, arr1, rng),
    c2, arr2, rng);
}

static void
gk_field_combine_disabled(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  // Do nothing.
}

static void
gk_field_step_apar_enabled(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array* out, double dt,
  const struct gkyl_array* inp)
{
  // Apar^{n+1} = Apar^{n} + dt*dApar/dt
  gkyl_array_accumulate(gkyl_array_scale(out, dt), 1.0, inp);

  // Smooth Apar along z.
  field->fem_projection_par_apar_func(app, field, out, out);
}

static void
gk_field_step_apar_disabled(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array* out, double dt,
  const struct gkyl_array* inp)
{
  // Do nothing.
}

static void
gk_field_copy_range_enabled(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  gkyl_array_copy_range(out, inp, range);
}

static void 
gk_field_copy_range_disabled(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  // Do nothing.
}

static void
gk_field_time_rate_diags_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  f->calc_energy_dt_func = gk_field_calc_energy_dt_active;
  if (app->use_gpu) {
    f->em_energy_red_new = gkyl_cu_malloc(sizeof(double[1]));
    f->em_energy_red_old = gkyl_cu_malloc(sizeof(double[1]));
    gkyl_cu_memset(f->em_energy_red_new, 0, sizeof(double[1]));
    gkyl_cu_memset(f->em_energy_red_old, 0, sizeof(double[1]));
  } else {
    f->em_energy_red_new = gkyl_malloc(sizeof(double[1]));
    f->em_energy_red_old = gkyl_malloc(sizeof(double[1]));
    memset(f->em_energy_red_new, 0, sizeof(double[1]));
    memset(f->em_energy_red_old, 0, sizeof(double[1]));
  }
  f->integ_energy_dot = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  f->is_first_energy_dot_write_call = true;
  if (f->is_em){
    f->calc_apar_energy_dt_func = gk_field_calc_apar_energy_dt_active;
    if (app->use_gpu) {
      f->apar_energy_red_new = gkyl_cu_malloc(sizeof(double[1]));
      f->apar_energy_red_old = gkyl_cu_malloc(sizeof(double[1]));
      gkyl_cu_memset(f->apar_energy_red_new, 0, sizeof(double[1]));
      gkyl_cu_memset(f->apar_energy_red_old, 0, sizeof(double[1]));
    } else {
      f->apar_energy_red_new = gkyl_malloc(sizeof(double[1]));
      f->apar_energy_red_old = gkyl_malloc(sizeof(double[1]));
      memset(f->apar_energy_red_new, 0, sizeof(double[1]));
      memset(f->apar_energy_red_old, 0, sizeof(double[1]));
    }
    f->integ_apar_energy_dot = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  }
  f->is_first_apar_energy_dot_write_call = true;
}

static void
gk_field_energy_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  // Allocate energy reduction arrays.
  if (app->use_gpu) {
    f->em_energy_red = gkyl_cu_malloc(sizeof(double[1]));
    f->em_energy_red_global = gkyl_cu_malloc(sizeof(double[1]));
  } else {
    f->em_energy_red = gkyl_malloc(sizeof(double[1]));
    f->em_energy_red_global = gkyl_malloc(sizeof(double[1]));
  }

  f->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  f->is_first_energy_write_call = true;

  f->calc_energy_func = gk_field_calc_energy_enabled;
  f->calc_energy_dt_func = gk_field_calc_energy_dt_none;
  f->calc_apar_energy_dt_func = gk_field_calc_energy_dt_none;
  
  if (f->info.time_rate_diagnostics) {
    gk_field_time_rate_diags_new(app, f);
  }

  // Factors for ES energy.
  f->es_energy_fac = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);
  f->es_energy_fac_1d = 0.0;

  if (f->is_em) {
    f->integ_apar_energy = gkyl_dynvec_new(GKYL_DOUBLE, 1);
    f->integ_apardot_energy = gkyl_dynvec_new(GKYL_DOUBLE, 1);
    f->apar_energy_fac = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);
    f->apar_energy_fac_1d = 0;
  }
}

static void
gk_field_time_rate_diags_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  f->calc_energy_dt_func = gk_field_calc_energy_dt_none;
  if (app->use_gpu) {
    gkyl_cu_free(f->em_energy_red_new);
    gkyl_cu_free(f->em_energy_red_old);
  } else {
    gkyl_free(f->em_energy_red_new);
    gkyl_free(f->em_energy_red_old);
  }
  gkyl_dynvec_release(f->integ_energy_dot);
  if (f->is_em) {
    if (app->use_gpu) {
      gkyl_cu_free(f->apar_energy_red_new);
      gkyl_cu_free(f->apar_energy_red_old);
    } else {
      gkyl_free(f->apar_energy_red_new);
      gkyl_free(f->apar_energy_red_old);
    }
    gkyl_dynvec_release(f->integ_apar_energy_dot);
  }
}

static void
gk_field_energy_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_dynvec_release(f->integ_energy);

  if (app->use_gpu) {
    gkyl_cu_free(f->em_energy_red);
    gkyl_cu_free(f->em_energy_red_global);
  } else {
    gkyl_free(f->em_energy_red);
    gkyl_free(f->em_energy_red_global);
  }

  if (f->info.time_rate_diagnostics) {
    gk_field_time_rate_diags_release(app, f);
  }

  gkyl_array_release(f->es_energy_fac);

  if (f->is_em) {
    gkyl_dynvec_release(f->integ_apar_energy);
    gkyl_dynvec_release(f->integ_apardot_energy);
    gkyl_array_release(f->apar_energy_fac);
  }
}

// Initialize field object.
struct gk_field* 
gk_field_new(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app)
{
  struct gk_field *f = gkyl_malloc(sizeof(struct gk_field));

  f->info = gk->field;

  f->gkfield_id = f->info.gkfield_id ? f->info.gkfield_id : GKYL_GK_FIELD;
  
  f->is_em = f->info.gkfield_id == GKYL_GK_FIELD_EM;
  // Ensure that if any species is electromagnetic, the field is electromagnetic.
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    f->is_em = f->is_em || (s->info.collisionless.type == GKYL_GK_COLLISIONLESS_EM);
  }
  if (f->is_em) {
    assert(f->info.mu0 > 0.0);
    f->step_apar = gk_field_step_apar_enabled;
    f->em_combine_func = gk_field_combine_enabled;
    f->em_copy_func = gk_field_copy_range_enabled;
  } else {
    f->step_apar = gk_field_step_apar_disabled;
    f->em_combine_func = gk_field_combine_disabled;
    f->em_copy_func = gk_field_copy_range_disabled;
  }
  f->calc_init_apar = f->info.calc_init_apar;

  f->calc_init_field = !f->info.zero_init_field;
  f->update_field = !f->info.is_static;
  // The combination update_field=true, calc_init_field=false is not allowed.
  assert(!(f->update_field && (!f->calc_init_field)));

  // Initialize polarization potential if needed.
  f->init_phi_pol = false;
  if (f->info.polarization_potential) {
    gk_field_polarization_potential_new(f, app);
  } else if (f->info.polarization_potential_import.type != GKYL_IC_IMPORT_NONE) {
    gk_field_polarization_potential_from_file_new(f, app, f->info.polarization_potential_import);
  }
  
  // Initialize energy diagnostics.
  gk_field_energy_new(app, f);

  // Initialize the field solver
  if (f->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
    gk_field_fem_new_boltzmann(app, f);
  } else {
    if (app->cdim == 1) {
      gk_field_fem_new_1x(app, f);
    } else if (app->cdim > 1) {
      gk_field_fem_new_2x3x(app, f);
    }
  }

  // Initialize biased walls.
  gk_field_biased_wall_new(app, f);

  return f;
}

// RHS functions for calculations

void
gk_field_rhs(gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Compute the electrostatic potential.
  struct timespec wst = gkyl_wall_clock();
  field->rhs_phi_func(app, field); // set inside the initialization function
  app->stat.field_phi_solve_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_em_rhs(gkyl_gyrokinetic_app *app, struct gk_field *field, const struct gkyl_array *f_in[],  struct gkyl_array *rhs_in[])
{
  // Compute Apardot through Ohm's law.
  struct timespec wst = gkyl_wall_clock();
  field->em_rhs_func(app, field, f_in, rhs_in);
  app->stat.field_tm += gkyl_time_diff_now_sec(wst);
}

void gk_field_calc_apar_ic(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *out)
{
  struct timespec wst = gkyl_wall_clock();
  field->ampere_solve(app, field, out);
  // Smooth Apar after solving Ampere's law.
  gk_field_fem_projection_par(app, field, out, out);
  app->stat.field_apar_solve_tm += gkyl_time_diff_now_sec(wst);
}

void gk_field_step_apar(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array* out, double a, const struct gkyl_array* inp)
{
  struct timespec wst = gkyl_wall_clock();
  field->step_apar(app, field, out, a, inp);
  app->stat.field_apar_solve_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_em_combine(struct gk_field *field, struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  field->em_combine_func(out, c1, arr1, c2, arr2, rng);
}

void
gk_field_em_copy_range(struct gk_field *field, struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  field->em_copy_func(out, inp, range);
}

void
gk_field_calc_energy(gkyl_gyrokinetic_app *app, double tm, const struct gk_field *field)
{
  struct timespec wst = gkyl_wall_clock();
  field->calc_energy_func(app, field, tm);
  app->stat.field_diag_calc_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_calc_energy_dt(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced)
{
  field->calc_energy_dt_func(app, field, dt, energy_reduced);
}

void
gk_field_calc_apar_energy_dt(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced)
{
  field->calc_apar_energy_dt_func(app, field, dt, energy_reduced);
}

void gk_field_accumulate_rho_c_adiabatic(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gk_species *s, struct gkyl_array **bflux)
{
  // Gyroaverage the density if needed.
  s->gyroaverage(app, s, s->m0.marr, s->m0_gyroavg);
  gkyl_array_accumulate_range(field->rho_c, s->info.charge, s->m0_gyroavg, &app->local);
  // Add the background (electron) charge density.
  double n_s0 = field->info.electron_density;
  double q_s = field->info.electron_charge;
  double dg_norm = pow(sqrt(2.0), app->basis.ndim);
  gkyl_array_shiftc_range(field->rho_c, q_s * n_s0 * dg_norm, 0, &app->local);
}

void gk_field_accumulate_rho_c_poisson(gkyl_gyrokinetic_app *app,
  struct gk_field *field, struct gk_species *s, struct gkyl_array **bflux)
{
  // Gyroaverage the density if needed.
  s->gyroaverage(app, s, s->m0.marr, s->m0_gyroavg);
  gkyl_array_accumulate_range(field->rho_c, s->info.charge, s->m0_gyroavg, &app->local);
}

void
gk_field_accumulate_rho_c(gkyl_gyrokinetic_app *app, struct gk_field *field, 
  const struct gkyl_array *fin[], struct gkyl_array **bflux[])
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_clear(field->rho_c, 0.0);
  for (int i = 0; i < app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    gk_species_moment_calc(&s->m0, s->local, app->local, fin[i]);
    field->accumulate_rhoc_func(app, field, s, bflux[i]);
  } 
  app->stat.field_phi_rhs_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_accumulate_current_dens(gkyl_gyrokinetic_app *app, struct gk_field *field, 
  const struct gkyl_array *fin[])
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_clear(field->currentDens, 0.0);
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    gk_species_moment_calc(&s->m1, s->local, app->local, fin[i]);
    s->gyroaverage(app, s, s->m1.marr, s->m1_gyroavg);
    gkyl_array_accumulate_range(field->currentDens, s->info.charge, s->m1_gyroavg, &app->local);
  } 
  app->stat.field_apar_rhs_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_accumulate_current_dens_dot(gkyl_gyrokinetic_app *app, struct gk_field *field, 
  struct gkyl_array *rhs_in[])
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_clear(field->currentDensdot, 0.0);
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    gk_species_moment_calc(&s->m1, s->local, app->local, rhs_in[i]);
    s->gyroaverage(app, s, s->m1.marr, s->m1_gyroavg);
    gkyl_array_accumulate_range(field->currentDensdot, s->info.charge, s->m1_gyroavg, &app->local);
  } 
  app->stat.field_apar_rhs_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_accumulate_ohms_kSq(gkyl_gyrokinetic_app *app, struct gk_field *field, 
  const struct gkyl_array *fin[])
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_clear(field->dApartdtSlvr_kSq, 0.0);
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    gk_species_moment_calc(&s->m0, s->local, app->local, fin[i]);
    // Gyroaverage the density if needed.
    s->gyroaverage(app, s, s->m0.marr, s->m0_gyroavg);
    // Use m0 to update also update the kSq matrix for the Ohm's law solver.
    double fac = -s->info.charge*s->info.charge/s->info.mass;
    gkyl_array_accumulate_range(field->dApartdtSlvr_kSq, fac, s->m0_gyroavg, &app->local);
  }
  app->stat.field_apar_rhs_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_fem_projection_par(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Project a DG field onto the parallel FEM basis to make it
  // continuous along z (or to solve a Poisson equation in 1x).

  // Gather the DG array into a global (in z) array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_dg, field->rho_c_global_dg);

  // Smooth the the DG array.
  gkyl_fem_parproj_set_rhs(field->fem_parproj, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj, field->phi_fem);

  // Copy global, continuous FEM array to a local array.
  gkyl_array_copy_range_to_range(arr_fem, field->phi_fem, &app->local, &field->global_sub_range);
}

void
gk_field_file_import_init(struct gkyl_gyrokinetic_app *app, struct gkyl_gyrokinetic_ic_import inp)
{
  // Import initial condition from a file.
  struct gkyl_app_restart_status rstat = header_from_file(app, inp.file_name);

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    struct gkyl_app_restart_status rstat;
    rstat.io_status = gkyl_comm_array_read(app->comm, &app->grid, &app->local, app->field->phi_host, inp.file_name);
    gkyl_array_copy(app->field->phi_smooth, app->field->phi_host);
  }
  else {
    assert(false);
  }
}

void
gk_field_project_init(struct gkyl_gyrokinetic_app *app)
{
  // Project the initial field.
  struct gkyl_eval_on_nodes *phi_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis,
    1, app->field->info.init_field_profile, app->field->info.init_field_profile_ctx);
  gkyl_eval_on_nodes_advance(phi_proj, 0.0, &app->local, app->field->phi_host);
  gkyl_eval_on_nodes_release(phi_proj);
  gkyl_array_copy(app->field->phi_smooth, app->field->phi_host);
}


// Release resources for field.
void
gk_field_release(const gkyl_gyrokinetic_app* app, struct gk_field *f)
{
  // Release polarization potential if used.
  if (f->init_phi_pol) {
    gk_field_polarization_potential_release(f);
  }

  // Release solver-specific resources.
  f->release_func(app, f);

  // Release energy diagnostics.
  gk_field_energy_release(app, f);

  // Release biased walls.
  gk_field_biased_wall_release(app, f);

  gkyl_free(f);
}
