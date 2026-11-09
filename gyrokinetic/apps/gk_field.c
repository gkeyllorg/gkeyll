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
  struct gkyl_app_restart_status rstat = {.io_status = GKYL_ARRAY_RIO_FOPEN_FAILED};

  FILE *fp = 0;
  with_file(fp, fname, "r")
  {
    struct gkyl_rect_grid grid;
    struct gkyl_array_header_info hdr;
    rstat.io_status = gkyl_grid_sub_array_header_read_fp(&grid, &hdr, fp);

    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      if (hdr.etype != GKYL_DOUBLE) {
        rstat.io_status = GKYL_ARRAY_RIO_DATA_MISMATCH;
      }
    }

    struct gkyl_msgpack_map_elem elem_list[] = {
      {.key = "frame", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0},
      {.key = "time", .elem_type = GKYL_MP_DOUBLE, .cval = 0}
    };
    int elem_list_len = sizeof(elem_list) / sizeof(elem_list[0]);
    gkyl_msgpack_to_map_elem_list(
      &(struct gkyl_msgpack_data){.meta = hdr.meta, .meta_sz = hdr.meta_size}, elem_list_len,
      elem_list
    );

    rstat.frame = gkyl_msgpack_map_elem_get_uint(elem_list_len, elem_list, "frame");
    rstat.stime = gkyl_msgpack_map_elem_get_double(elem_list_len, elem_list, "time");

    gkyl_grid_sub_array_header_release(&hdr);
  }

  return rstat;
}

static void
gk_field_polarization_potential_new(struct gk_field *gkf, struct gkyl_gyrokinetic_app *app)
{
  // Project the initial potential onto a p+1 tensor basis and compute the polarization
  // density to use use by species in calculating the initial ion density.
  gkf->init_phi_pol = true;
  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, app->cdim, app->poly_order + 1);

  gkf->phi_pol = mkarr(app->use_gpu, phi_pol_basis.num_basis, app->local_ext.volume);
  struct gkyl_array *phi_pol_ho = app->use_gpu ?
                                    mkarr(false, gkf->phi_pol->ncomp, gkf->phi_pol->size) :
                                    gkyl_array_acquire(gkf->phi_pol);

  struct gkyl_eval_on_nodes *phi_pol_proj =
    gkyl_eval_on_nodes_inew(&(struct gkyl_eval_on_nodes_inp){
      .grid = &app->grid,
      .basis = &phi_pol_basis,
      .num_ret_vals = 1,
      .eval = gkf->info.polarization_potential,
      .ctx = gkf->info.polarization_potential_ctx,
      .c2p_func = eval_on_nodes_c2p_position_func,
      .c2p_func_ctx = app->position_map,
    });

  gkyl_eval_on_nodes_advance(phi_pol_proj, 0.0, &app->local, phi_pol_ho);
  gkyl_array_copy(gkf->phi_pol, phi_pol_ho);

  gkyl_eval_on_nodes_release(phi_pol_proj);
  gkyl_array_release(phi_pol_ho);
}

static void
gk_field_polarization_potential_from_file_new(
  struct gk_field *gkf, struct gkyl_gyrokinetic_app *app, struct gkyl_gyrokinetic_ic_import inp
)
{
  gkf->init_phi_pol = true;
  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, app->cdim, app->poly_order + 1);

  gkf->phi_pol = mkarr(app->use_gpu, phi_pol_basis.num_basis, app->local_ext.volume);
  struct gkyl_array *phi_pol_ho = app->use_gpu ?
                                    mkarr(false, gkf->phi_pol->ncomp, gkf->phi_pol->size) :
                                    gkyl_array_acquire(gkf->phi_pol);

  struct gkyl_app_restart_status rstat = header_from_file(app, inp.file_name);
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, phi_pol_ho, inp.file_name);
    assert(rstat.io_status == GKYL_ARRAY_RIO_SUCCESS);
    gkyl_array_copy(gkf->phi_pol, phi_pol_ho);
  } else {
    assert(false);
  }

  gkyl_array_release(phi_pol_ho);
}

static void
gk_field_polarization_potential_release(struct gk_field *gkf)
{
  gkyl_array_release(gkf->phi_pol);
}

// Functions related to the field energy allocations, diagnostics, and release

// Allocate a zeroed double for a reduction (on device if using GPUs).
static double *
gk_field_energy_red_new(bool use_gpu)
{
  double *out = use_gpu ? gkyl_cu_malloc(sizeof(double[1])) : gkyl_malloc(sizeof(double[1]));
  if (use_gpu) {
    gkyl_cu_memset(out, 0, sizeof(double[1]));
  } else {
    memset(out, 0, sizeof(double[1]));
  }
  return out;
}

static void
gk_field_energy_red_release(bool use_gpu, double *buf)
{
  if (use_gpu) {
    gkyl_cu_free(buf);
  } else {
    gkyl_free(buf);
  }
}

// factor times the field energy in the local range, without reducing it.
static void
gk_field_calc_energy_local(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, double factor, double *out
)
{
  gkyl_array_integrate_advance(
    field->calc_em_energy, field->phi_smooth, factor, field->es_energy_fac, &app->local,
    &app->local, out
  );
  field->energy_adiab_func(app, field, factor, out);
}

// Sum a local energy over the communicator.
static double
gk_field_energy_reduce(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, const double *energy_local
)
{
  gkyl_comm_allreduce(
    app->comm, GKYL_DOUBLE, GKYL_SUM, 1, energy_local, field->em_energy_red_global
  );
  double energy;
  if (app->use_gpu) {
    gkyl_cu_memcpy(&energy, field->em_energy_red_global, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    energy = field->em_energy_red_global[0];
  }
  return energy;
}

// The 1x field energy uses a scalar factor.
static double
gk_field_em_energy_reduce(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, const double *energy_local
)
{
  double energy = gk_field_energy_reduce(app, field, energy_local);
  return app->cdim == 1 ? energy * field->es_energy_fac_1d : energy;
}

static void
gk_field_energy_bg_none(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, double factor, double *out
)
{
}

static void
gk_field_calc_energy_dt_active(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, bool is_old
)
{
  struct timespec wst = gkyl_wall_clock();
  gk_field_calc_energy_local(
    app, field, 1.0 / dt, is_old ? field->em_energy_red_old : field->em_energy_red_new
  );
  field->energy_bg_func(
    app, field, 1.0 / dt, is_old ? field->bg_energy_red_old : field->bg_energy_red_new
  );
  app->stat.phidot_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_field_calc_energy_dt_none(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, bool is_old
)
{
  // Do nothing.
}

static void
gk_field_calc_energy_enabled(
  struct gkyl_gyrokinetic_app *app, const struct gk_field *field, double tm
)
{
  gk_field_calc_energy_local(app, field, 1.0, field->em_energy_red);
  double energy_global[1] = {gk_field_em_energy_reduce(app, field, field->em_energy_red)};
  gkyl_dynvec_append(field->integ_energy, tm, energy_global);

  if (field->info.time_rate_diagnostics) {
    double energy_dot_global[1] = {
      gk_field_em_energy_reduce(app, field, field->em_energy_red_new) -
      gk_field_em_energy_reduce(app, field, field->em_energy_red_old)
    };
    gkyl_dynvec_append(field->integ_energy_dot, tm, energy_dot_global);
  }
}

static void
gk_field_calc_energy_disabled(
  struct gkyl_gyrokinetic_app *app, const struct gk_field *field, double tm
)
{
  // Do nothing.
}

static void
gk_field_energy_dt_new(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  // Energies/dt at the beginning and end of a time step, for the energy rate
  // diagnostic and the energy balance of the app.
  gkf->calc_energy_dt_func = gk_field_calc_energy_dt_active;
  gkf->em_energy_red_old = gk_field_energy_red_new(app->use_gpu);
  gkf->em_energy_red_new = gk_field_energy_red_new(app->use_gpu);
  gkf->bg_energy_red = gk_field_energy_red_new(app->use_gpu);
  gkf->bg_energy_red_old = gk_field_energy_red_new(app->use_gpu);
  gkf->bg_energy_red_new = gk_field_energy_red_new(app->use_gpu);
  if (gkf->info.time_rate_diagnostics) {
    gkf->integ_energy_dot = gkyl_dynvec_new(GKYL_DOUBLE, 1);
    gkf->is_first_energy_dot_write_call = true;
  }
}

static void
gk_field_energy_new(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  // Allocate energy reduction arrays.
  gkf->em_energy_red = gk_field_energy_red_new(app->use_gpu);
  gkf->em_energy_red_global = gk_field_energy_red_new(app->use_gpu);

  gkf->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  gkf->is_first_energy_write_call = true;

  gkf->calc_energy_func = gk_field_calc_energy_enabled;
  gkf->calc_energy_dt_func = gk_field_calc_energy_dt_none;
  gkf->energy_bg_func = gk_field_energy_bg_none;

  // The energy balance of the app needs the energy rates and the background-charge energy.
  gkf->energy_balance = false;
  for (int d = 0; d < app->bal.num_diags; d++) {
    gkf->energy_balance = gkf->energy_balance || app->bal.diag_types[d] == GKYL_BALANCE_GK_ENERGY;
  }
  if (gkf->info.time_rate_diagnostics || gkf->energy_balance) {
    gk_field_energy_dt_new(app, gkf);
  }

  // Factors for ES energy.
  gkf->es_energy_fac =
    mkarr(app->use_gpu, (2 * (app->cdim / 3) + 1) * app->basis.num_basis, app->local_ext.volume);
  gkf->es_energy_fac_1d = 0.0;
}

static void
gk_field_energy_dt_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gkf->calc_energy_dt_func = gk_field_calc_energy_dt_none;
  gk_field_energy_red_release(app->use_gpu, gkf->em_energy_red_old);
  gk_field_energy_red_release(app->use_gpu, gkf->em_energy_red_new);
  gk_field_energy_red_release(app->use_gpu, gkf->bg_energy_red);
  gk_field_energy_red_release(app->use_gpu, gkf->bg_energy_red_old);
  gk_field_energy_red_release(app->use_gpu, gkf->bg_energy_red_new);
  if (gkf->info.time_rate_diagnostics) {
    gkyl_dynvec_release(gkf->integ_energy_dot);
  }
}

static void
gk_field_energy_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gkyl_dynvec_release(gkf->integ_energy);
  gk_field_energy_red_release(app->use_gpu, gkf->em_energy_red);
  gk_field_energy_red_release(app->use_gpu, gkf->em_energy_red_global);

  if (gkf->info.time_rate_diagnostics || gkf->energy_balance) {
    gk_field_energy_dt_release(app, gkf);
  }

  gkyl_array_release(gkf->es_energy_fac);
}

static void
gk_field_accumulate_rho_c_bg_none(gkyl_gyrokinetic_app *app, struct gk_field *field)
{
}

static void
gk_field_adiabatic_energy_none(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, double factor, double *out
)
{
}

// Initialize field object.
struct gk_field *
gk_field_new(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app)
{
  struct gk_field *gkf = gkyl_malloc(sizeof(struct gk_field));

  gkf->info = gk->field;

  gkf->gkfield_id = gkf->info.gkfield_id ? gkf->info.gkfield_id : GKYL_GK_FIELD_ES;

  gkf->calc_init_field = !gkf->info.zero_init_field;
  gkf->update_field = !gkf->info.is_static;
  // The combination update_field=true, calc_init_field=false is not allowed.
  assert(!(gkf->update_field && (!gkf->calc_init_field)));

  // Initialize polarization potential if needed.
  gkf->init_phi_pol = false;
  if (gkf->info.polarization_potential) {
    gk_field_polarization_potential_new(gkf, app);
  } else if (gkf->info.polarization_potential_import.type != GKYL_IC_IMPORT_NONE) {
    gk_field_polarization_potential_from_file_new(gkf, app, gkf->info.polarization_potential_import);
  }

  // Initialize energy diagnostics.
  gk_field_energy_new(app, gkf);

  // Species with an adiabatic response enter the field equation through their background
  // density and K = sum_s (q_s^2 n_s0/T_s)*J instead of their charge density (see
  // gk_field_adiabatic.c).
  memset(&gkf->gk_adiab, 0, sizeof(gkf->gk_adiab));
  gkf->num_rhoc_species = 0;
  for (int i = 0; i < app->num_species; ++i) {
    if (app->species[i].info.scaling.type == GKYL_GK_SPECIES_SCALING_ADIABATIC) {
      gkf->gk_adiab.species_idx[gkf->gk_adiab.num_species++] = i;
    } else {
      gkf->rhoc_species_idx[gkf->num_rhoc_species++] = i;
    }
  }
  gkf->has_adiabatic_species = gkf->gk_adiab.num_species > 0;
  gkf->accumulate_rhoc_bg_func = gk_field_accumulate_rho_c_bg_none;
  gkf->energy_adiab_func = gk_field_adiabatic_energy_none;

  // Initialize the field solver
  if (gkf->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
    gk_field_fem_new_boltzmann(app, gkf);
  } else {
    if (app->cdim == 1) {
      gk_field_fem_new_1x(app, gkf);
    } else if (app->cdim > 1) {
      gk_field_fem_new_2x3x(app, gkf);
    }
  }

  // Set up the adiabatic response, which may replace the field solve.
  if (gkf->has_adiabatic_species) {
    gk_field_adiabatic_new(app, gkf);
  }

  return gkf;
}

void
gk_field_fem_parproj_new(
  struct gkyl_gyrokinetic_app *app, struct gk_field *gkf, const struct gkyl_array *weight
)
{
  enum gkyl_fem_parproj_bc_type fem_parproj_bc = GKYL_FEM_PARPROJ_NONE;
  for (int d = 0; d < app->num_periodic_dir; ++d) {
    if (app->periodic_dirs[d] == app->cdim - 1) {
      fem_parproj_bc = GKYL_FEM_PARPROJ_PERIODIC;
    }
  }

  // Gather the weight for the (global) projection in z.
  struct gkyl_array *weight_global = 0;
  if (weight) {
    weight_global = mkarr(app->use_gpu, weight->ncomp, app->global_ext.volume);
    gkyl_comm_array_allgather(app->comm, &app->local, &app->global, weight, weight_global);
  }

  gkf->fem_parproj = gkyl_fem_parproj_new(
    &app->global, &app->grid, &app->basis, fem_parproj_bc, 0, weight_global, 0, app->use_gpu
  );

  if (weight) {
    gkyl_array_release(weight_global);
  }
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
gk_field_calc_energy(gkyl_gyrokinetic_app *app, double tm, const struct gk_field *field)
{
  struct timespec wst = gkyl_wall_clock();
  field->calc_energy_func(app, field, tm);
  app->stat.field_diag_calc_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_calc_energy_dt(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, bool is_old
)
{
  field->calc_energy_dt_func(app, field, dt, is_old);
}

void
gk_field_calc_energy_balance(
  gkyl_gyrokinetic_app *app, const struct gk_field *field, double *energy, double *energy_dot,
  double *bg_energy, double *bg_energy_dot
)
{
  // Energies now, and their rates of change from the energies/dt stored by
  // gk_field_calc_energy_dt at the beginning and end of the last time step.
  gk_field_calc_energy_local(app, field, 1.0, field->em_energy_red);
  *energy = gk_field_em_energy_reduce(app, field, field->em_energy_red);
  *energy_dot = gk_field_em_energy_reduce(app, field, field->em_energy_red_new) -
                gk_field_em_energy_reduce(app, field, field->em_energy_red_old);
  field->energy_bg_func(app, field, 1.0, field->bg_energy_red);
  *bg_energy = gk_field_energy_reduce(app, field, field->bg_energy_red);
  *bg_energy_dot = gk_field_energy_reduce(app, field, field->bg_energy_red_new) -
                   gk_field_energy_reduce(app, field, field->bg_energy_red_old);
}

void
gk_field_accumulate_rho_c_poisson(
  gkyl_gyrokinetic_app *app, struct gk_field *field, struct gk_species *gks,
  struct gkyl_array **bflux
)
{
  // Gyroaverage the density if needed.
  gks->gyroaverage(app, gks, gks->m0.marr, gks->m0_gyroavg);
  gkyl_array_accumulate_range(field->rho_c, gks->info.charge, gks->m0_gyroavg, &app->local);
}

void
gk_field_accumulate_rho_c(
  gkyl_gyrokinetic_app *app, struct gk_field *field, const struct gkyl_array *fin[],
  struct gkyl_array **bflux[]
)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_clear(field->rho_c, 0.0);
  for (int k = 0; k < field->num_rhoc_species; ++k) {
    int i = field->rhoc_species_idx[k];
    struct gk_species *gks = &app->species[i];
    gk_species_moment_calc(app, &gks->m0, gks->local, app->local, fin[i]);
    field->accumulate_rhoc_func(app, field, gks, bflux[i]);
  }
  field->accumulate_rhoc_bg_func(app, field);
  app->stat.field_phi_rhs_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_fem_projection_par(
  gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *arr_dg,
  struct gkyl_array *arr_fem
)
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
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, app->field->phi_host, inp.file_name);
    gkyl_array_copy(app->field->phi_smooth, app->field->phi_host);
  } else {
    assert(false);
  }
}

void
gk_field_project_init(struct gkyl_gyrokinetic_app *app)
{
  // Project the initial field.
  struct gkyl_eval_on_nodes *phi_proj = gkyl_eval_on_nodes_new(
    &app->grid, &app->basis, 1, app->field->info.init_field_profile,
    app->field->info.init_field_profile_ctx
  );
  gkyl_eval_on_nodes_advance(phi_proj, 0.0, &app->local, app->field->phi_host);
  gkyl_eval_on_nodes_release(phi_proj);
  gkyl_array_copy(app->field->phi_smooth, app->field->phi_host);
}

// Release resources for field.
void
gk_field_release(const gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  // Release polarization potential if used.
  if (gkf->init_phi_pol) {
    gk_field_polarization_potential_release(gkf);
  }

  // Release solver-specific resources.
  gkf->release_func(app, gkf);

  if (gkf->has_adiabatic_species) {
    gk_field_adiabatic_release(app, gkf);
  }

  // Release energy diagnostics.
  gk_field_energy_release(app, gkf);

  gkyl_free(gkf);
}
