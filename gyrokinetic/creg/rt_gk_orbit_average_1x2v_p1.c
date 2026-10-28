#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gkyl_alloc.h>
#include <gkyl_dynvec.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_util.h>

#include <rt_arg_parse.h>

// Normalized electrostatic oscillator with reflecting endpoints. This input
// exercises the app's frozen-field OAP path, including MPI and CUDA dispatch.
// The zero-layer unit tests independently verify the projector's moments.
static void
orbit_mapc2p(double time, const double *xc, double *xp, void *ctx)
{
  xp[0] = xc[0];
  xp[1] = xc[1];
  xp[2] = xc[2];
}

static void
orbit_bfield(double time, const double *xc, double *out, void *ctx)
{
  out[0] = 0.0;
  out[1] = 0.0;
  out[2] = 1.0;
}

static void
orbit_phi(double time, const double *xn, double *out, void *ctx)
{
  out[0] = 0.5 * xn[0] * xn[0];
}

static void
orbit_initial(double time, const double *xn, double *out, void *ctx)
{
  // The positive polynomial has a resolved orbit-varying component in z.
  out[0] = 1.0 + 0.15 * xn[0] + 0.1 * xn[1] * xn[1] + 0.05 * xn[2];
}

static bool
orbit_check_moments(gkyl_gyrokinetic_app *app, const char *app_name, int num_steps)
{
  char file_name[256];
  snprintf(file_name, sizeof(file_name), "%s-ion_integrated_moms.gkyl", app_name);
  gkyl_dynvec moments = gkyl_dynvec_new(GKYL_DOUBLE, 3);
  bool valid = gkyl_dynvec_read(moments, file_name) &&
               gkyl_dynvec_size(moments) == (size_t)(num_steps + 4);
  if (valid) {
    double initial[3], current[3], before_disabled[3], after_disabled[3], before_rejection[3],
      after_rejection[3];
    gkyl_dynvec_get(moments, 0, initial);
    // The diagnostic contains particle number, momentum and H. With initial
    // collisionless NONE, orbit relaxation alone must conserve number and H.
    for (int step = 1; step <= num_steps + 3; ++step) {
      gkyl_dynvec_get(moments, step, current);
      for (int moment = 0; moment < 3; moment += 2) {
        valid = valid && fabs(current[moment] - initial[moment]) <=
                           2.0e-11 * fmax(1.0, fabs(initial[moment]));
      }
    }
    gkyl_dynvec_get(moments, num_steps, before_disabled);
    gkyl_dynvec_get(moments, num_steps + 1, after_disabled);
    gkyl_dynvec_get(moments, num_steps + 2, before_rejection);
    gkyl_dynvec_get(moments, num_steps + 3, after_rejection);
    for (int moment = 0; moment < 3; ++moment) {
      valid = valid && fabs(after_disabled[moment] - before_disabled[moment]) <=
                         2.0e-11 * fmax(1.0, fabs(before_disabled[moment]));
      valid = valid && fabs(after_rejection[moment] - before_rejection[moment]) <=
                         2.0e-11 * fmax(1.0, fabs(before_rejection[moment]));
    }
    gkyl_gyrokinetic_app_cout(
      app, stdout, "Orbit-only moments: initial N=%.16g H=%.16g; final N=%.16g H=%.16g.\n",
      initial[0], initial[2], current[0], current[2]
    );
  }
  gkyl_dynvec_release(moments);
  return valid;
}

static bool
orbit_check_collisionless_independence(const struct gkyl_gk *reference)
{
  // Orbit reset must neither suppress an existing ES solver nor alter its CFL
  // condition. Collisionless reset independently disables and restores ES.
  struct gkyl_gk *inp = gkyl_malloc(sizeof(*inp));
  *inp = *reference;
  snprintf(inp->name, sizeof(inp->name), "%.100s_es", reference->name);
  inp->species[0].collisionless.type = GKYL_GK_COLLISIONLESS_ES;
  inp->species[0].orbit_average.frequency = 0.0;
  inp->species[0].boundary_flux_diagnostics = (struct gkyl_phase_diagnostics_inp){
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = {GKYL_F_MOMENT_M0M1M2PARM2PERP},
  };
  struct gkyl_gyrokinetic_orbit_average options = inp->species[0].orbit_average;
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(inp);
  gkyl_gyrokinetic_app_apply_ic(app, 0.0);
  gkyl_gyrokinetic_app_calc_integrated_mom(app, 0.0);
  const double dt = 0.1;
  struct gkyl_update_status update = gkyl_gyrokinetic_update(app, dt);
  struct gkyl_gyrokinetic_orbit_average_stat stat =
    gkyl_gyrokinetic_app_species_orbit_average_stat(app, "ion");
  bool valid = update.success && update.dt_actual > 0.0 && update.dt_actual < dt && stat.prepared &&
               stat.num_selected_cells > 0;
  double time = update.dt_actual;
  gkyl_gyrokinetic_app_calc_integrated_mom(app, time);

  options.enabled = false;
  gkyl_gyrokinetic_app_reset_species_orbit_average(app, "ion", options);
  update = gkyl_gyrokinetic_update(app, dt);
  stat = gkyl_gyrokinetic_app_species_orbit_average_stat(app, "ion");
  valid = valid && update.success && update.dt_actual > 0.0 && update.dt_actual < dt &&
          !stat.enabled && !stat.prepared;
  time += update.dt_actual;
  gkyl_gyrokinetic_app_calc_integrated_mom(app, time);

  options.enabled = true;
  options.frequency = 1.0e6;
  gkyl_gyrokinetic_app_reset_species_orbit_average(app, "ion", options);
  update = gkyl_gyrokinetic_update(app, dt);
  stat = gkyl_gyrokinetic_app_species_orbit_average_stat(app, "ion");
  valid = valid && update.success && update.dt_actual > 0.0 && update.dt_actual < dt &&
          stat.prepared && stat.num_selected_cells > 0 && stat.min_blend > 0.0;
  time += update.dt_actual;
  gkyl_gyrokinetic_app_calc_integrated_mom(app, time);

  // Turning off both terms must remove the advection CFL and leave f unchanged.
  options.enabled = false;
  gkyl_gyrokinetic_app_reset_species_orbit_average(app, "ion", options);
  gkyl_gyrokinetic_app_reset_species_collisionless(
    app, time, "ion", (struct gkyl_gyrokinetic_collisionless){.type = GKYL_GK_COLLISIONLESS_NONE}
  );
  update = gkyl_gyrokinetic_update(app, dt);
  valid = valid && update.success && fabs(update.dt_actual - dt) <= 1.0e-12 * dt;
  time += update.dt_actual;
  gkyl_gyrokinetic_app_calc_integrated_mom(app, time);

  // Restore only the projection limit: it must ignore a zero rate and not restrict dt.
  options.enabled = true;
  options.frequency = 0.0;
  options.frequency_is_infinite = true;
  gkyl_gyrokinetic_app_reset_species_orbit_average(app, "ion", options);
  update = gkyl_gyrokinetic_update(app, dt);
  stat = gkyl_gyrokinetic_app_species_orbit_average_stat(app, "ion");
  valid = valid && update.success && fabs(update.dt_actual - dt) <= 1.0e-12 * dt && stat.prepared &&
          stat.num_selected_cells > 0 && stat.min_blend > 0.0;
  time += update.dt_actual;
  gkyl_gyrokinetic_app_calc_integrated_mom(app, time);

  // Restore ES with zero orbit rate so its RHS and CFL can be checked alone.
  options.frequency = 0.0;
  options.frequency_is_infinite = false;
  gkyl_gyrokinetic_app_reset_species_orbit_average(app, "ion", options);
  gkyl_gyrokinetic_app_reset_species_collisionless(
    app, time, "ion",
    (struct gkyl_gyrokinetic_collisionless){.type = GKYL_GK_COLLISIONLESS_ES, .scale_factor = 1.0}
  );
  update = gkyl_gyrokinetic_update(app, dt);
  valid = valid && update.success && update.dt_actual > 0.0 && update.dt_actual < dt;
  time += update.dt_actual;
  gkyl_gyrokinetic_app_calc_integrated_mom(app, time);
  gkyl_gyrokinetic_app_write_integrated_mom(app);
  gkyl_comm_barrier(reference->parallelism.comm);
  int rank = 0;
  gkyl_comm_get_rank(reference->parallelism.comm, &rank);
  if (rank == 0) {
    char file_name[256];
    snprintf(file_name, sizeof(file_name), "%s-ion_integrated_moms.gkyl", inp->name);
    gkyl_dynvec moments = gkyl_dynvec_new(GKYL_DOUBLE, 3);
    bool readable = gkyl_dynvec_read(moments, file_name) && gkyl_dynvec_size(moments) == 7;
    valid = valid && readable;
    if (readable) {
      double initial[3], after_both[3], after_es[3];
      gkyl_dynvec_get(moments, 0, initial);
      gkyl_dynvec_get(moments, 1, after_both);
      gkyl_dynvec_get(moments, 2, after_es);
      // Zero orbit frequency isolates the ES RHS while both flags are enabled.
      valid = valid && fabs(after_both[1] - initial[1]) > 1.0e-10 &&
              fabs(after_es[1] - after_both[1]) > 1.0e-10;

      double before_none[3], after_none[3], after_orbit[3], after_restored[3];
      gkyl_dynvec_get(moments, 3, before_none);
      gkyl_dynvec_get(moments, 4, after_none);
      gkyl_dynvec_get(moments, 5, after_orbit);
      gkyl_dynvec_get(moments, 6, after_restored);
      for (int moment = 0; moment < 3; ++moment) {
        valid = valid && fabs(after_none[moment] - before_none[moment]) <=
                           2.0e-11 * fmax(1.0, fabs(before_none[moment]));
      }
      for (int moment = 0; moment < 3; moment += 2) {
        valid = valid && fabs(after_orbit[moment] - after_none[moment]) <=
                           2.0e-11 * fmax(1.0, fabs(after_none[moment]));
      }
      valid = valid && fabs(after_restored[1] - after_orbit[1]) > 1.0e-10;
    }
    gkyl_dynvec_release(moments);
  }
  int local_result = valid ? 0 : 1, global_result = 0;
  gkyl_comm_allreduce_host(
    reference->parallelism.comm, GKYL_INT, GKYL_MAX, 1, &local_result, &global_result
  );
  gkyl_gyrokinetic_app_cout(
    app, global_result ? stderr : stdout, "Independent ES/orbit controls: %s.\n",
    global_result ? "FAILED" : "passed"
  );
  gkyl_gyrokinetic_app_release(app);
  gkyl_free(inp);
  return global_result == 0;
}

static bool
orbit_check_boltzmann_initialization(const struct gkyl_gk *reference)
{
  // Keep this second large application input on the heap. A Boltzmann field
  // needs physical ion sheath fluxes during initialization, even when orbit
  // relaxation is enabled alongside ES in its frozen-field evolution.
  struct gkyl_gk *inp = gkyl_malloc(sizeof(*inp));
  *inp = *reference;
  snprintf(inp->name, sizeof(inp->name), "%.100s_boltzmann_init", reference->name);
  inp->species[0].bcs[0].type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING;
  inp->species[0].bcs[1].type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING;
  inp->species[0].collisionless.type = GKYL_GK_COLLISIONLESS_ES;
  inp->species[0].boundary_flux_diagnostics = (struct gkyl_phase_diagnostics_inp){
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = {GKYL_F_MOMENT_M0},
  };
  inp->species[0].time_rate_diagnostics = false;
  inp->species[0].num_integrated_diag_moments = 0;
  inp->field = (struct gkyl_gyrokinetic_field){
    .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,
    .is_static = true,
    .electron_mass = 1.0 / 1836.0,
    .electron_charge = -1.0,
    .electron_temp = 1.0,
  };
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(inp);
  gkyl_gyrokinetic_app_apply_ic(app, 0.0);
  struct gkyl_update_status update = gkyl_gyrokinetic_update(app, 0.1);
  struct gkyl_gyrokinetic_orbit_average_stat stat =
    gkyl_gyrokinetic_app_species_orbit_average_stat(app, "ion");
  // An empty trapped set is legitimate here. Preparation still validates
  // finite potential and positive geometry throughout the global domain.
  bool valid = update.success && stat.enabled && stat.prepared && update.dt_actual > 0.0 &&
               update.dt_actual < 0.1;
  double time = update.dt_actual;
  gkyl_gyrokinetic_app_calc_integrated_mom(app, time);

  // The sheath flux computed above must not survive an ES -> NONE transition.
  // The field remains frozen, as required by the orbit-averaged phase.
  gkyl_gyrokinetic_app_reset_species_collisionless(
    app, time, "ion", (struct gkyl_gyrokinetic_collisionless){.type = GKYL_GK_COLLISIONLESS_NONE}
  );
  update = gkyl_gyrokinetic_update(app, 0.1);
  valid = valid && update.success && fabs(update.dt_actual - 0.1) < 1.0e-13;
  time += update.dt_actual;
  gkyl_gyrokinetic_app_calc_integrated_mom(app, time);

  gkyl_gyrokinetic_app_reset_species_collisionless(
    app, time, "ion",
    (struct gkyl_gyrokinetic_collisionless){.type = GKYL_GK_COLLISIONLESS_ES, .scale_factor = 1.0}
  );
  update = gkyl_gyrokinetic_update(app, 0.1);
  valid = valid && update.success && update.dt_actual > 0.0 && update.dt_actual < 0.1;
  time += update.dt_actual;
  gkyl_gyrokinetic_app_calc_integrated_mom(app, time);
  gkyl_gyrokinetic_app_write_integrated_mom(app);
  gkyl_comm_barrier(reference->parallelism.comm);
  int rank = 0;
  gkyl_comm_get_rank(reference->parallelism.comm, &rank);
  if (rank == 0) {
    const char *edges[] = {"lower", "upper"};
    for (int edge = 0; edge < 2; ++edge) {
      char file_name[256];
      snprintf(
        file_name, sizeof(file_name), "%s-ion_bflux_x%s_integrated_M0.gkyl", inp->name, edges[edge]
      );
      gkyl_dynvec flux = gkyl_dynvec_new(GKYL_DOUBLE, 1);
      bool readable = gkyl_dynvec_read(flux, file_name) && gkyl_dynvec_size(flux) == 3;
      valid = valid && readable;
      if (readable) {
        double before[1], disabled[1], restored[1];
        gkyl_dynvec_get(flux, 0, before);
        gkyl_dynvec_get(flux, 1, disabled);
        gkyl_dynvec_get(flux, 2, restored);
        valid = valid && fabs(before[0]) > 1.0e-10 && fabs(disabled[0]) < 1.0e-13 &&
                fabs(restored[0]) > 1.0e-10;
      }
      gkyl_dynvec_release(flux);
    }
  }
  int local_result = valid ? 0 : 1, global_result = 0;
  gkyl_comm_allreduce_host(
    reference->parallelism.comm, GKYL_INT, GKYL_MAX, 1, &local_result, &global_result
  );
  gkyl_gyrokinetic_app_cout(
    app, global_result ? stderr : stdout,
    "Boltzmann/sheath initialization and ES/NONE/ES boundary fluxes: %s, selected=%ld.\n",
    global_result ? "FAILED" : "passed", stat.num_selected_cells
  );
  gkyl_gyrokinetic_app_release(app);
  gkyl_free(inp);
  return global_result == 0;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args args = parse_app_args(argc, argv);
#ifdef GKYL_HAVE_MPI
  if (args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif
  if (args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(args.use_mpi, args.use_gpu, stderr);
  const double frequency = 1.0e6, dt = 0.1;
  struct gkyl_gyrokinetic_orbit_average options = {
    .enabled = true,
    .frequency = frequency,
    .num_energy_nodes = 8,
    .positivity = true,
  };
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = 1.0,
    .mass = 1.0,
    .vdim = 2,
    .lower = {-2.5, 0.0},
    .upper = {2.5, 1.0},
    .cells = {APP_ARGS_CHOOSE(args.vcells[0], 12), APP_ARGS_CHOOSE(args.vcells[1], 2)},
    .polarization_density = 1.0,
    .projection = {.proj_id = GKYL_PROJ_FUNC, .func = orbit_initial},
    .collisionless = {.type = GKYL_GK_COLLISIONLESS_NONE},
    .orbit_average = options,
    .bcs =
      {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_REFLECT},
       {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_REFLECT}},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
    .time_rate_diagnostics = true,
  };
  struct gkyl_gk inp = {
    .cdim = 1,
    .lower = {-1.0},
    .upper = {1.0},
    .cells = {APP_ARGS_CHOOSE(args.xcells[0], 12)},
    .poly_order = 1,
    .basis_type = args.basis_type,
    .cfl_frac = 1.0,
    .geometry =
      {
        .geometry_id = GKYL_GEOMETRY_MAPC2P,
        .world = {0.0, 0.0},
        .mapc2p = orbit_mapc2p,
        .bfield_func = orbit_bfield,
      },
    .num_species = 1,
    .species = {ion},
    .field =
      {
        .gkfield_id = GKYL_GK_FIELD_ES,
        .is_static = true,
        .polarization_bmag = 1.0,
        .kperpSq = 1.0,
        .init_field_profile = orbit_phi,
      },
    .parallelism = {.use_gpu = args.use_gpu, .cuts = {args.cuts[0]}, .comm = comm},
  };
  snprintf(inp.name, sizeof(inp.name), "%s", args.app_name);
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&inp);
  gkyl_gyrokinetic_app_apply_ic(app, 0.0);
  gkyl_gyrokinetic_app_calc_integrated_mom(app, 0.0);

  int result = 0;
  const int num_steps = args.num_steps < 3 ? args.num_steps : 3;
  double time = 0.0;
  for (int step = 0; step < num_steps; ++step) {
    struct gkyl_update_status update = gkyl_gyrokinetic_update(app, dt);
    struct gkyl_gyrokinetic_orbit_average_stat stat =
      gkyl_gyrokinetic_app_species_orbit_average_stat(app, "ion");
    if (!update.success || fabs(update.dt_actual - dt) > 1.0e-12 * dt || !stat.prepared ||
        stat.num_selected_cells == 0 || stat.num_modes == 0 ||
        !(stat.min_blend > 0.0 && stat.min_blend <= 1.0)) {
      gkyl_gyrokinetic_app_cout(
        app, stderr, "Orbit-average app regression failed at step %d.\n", step
      );
      result = 1;
      break;
    }
    time += update.dt_actual;
    gkyl_gyrokinetic_app_calc_integrated_mom(app, time);
    gkyl_gyrokinetic_app_cout(
      app, stdout, "OAP step %d: dt=%.12g, nu*dt=%.12g, selected=%ld, modes=%ld, blend=%.12g.\n",
      step + 1, update.dt_actual, frequency * update.dt_actual, stat.num_selected_cells,
      stat.num_modes, stat.min_blend
    );
  }
  if (!result && num_steps > 0) {
    // The two selectors are independent: disabling only orbit relaxation
    // leaves collisionless NONE and hence an unchanged distribution.
    options.enabled = false;
    gkyl_gyrokinetic_app_reset_species_orbit_average(app, "ion", options);
    struct gkyl_gyrokinetic_orbit_average_stat stat =
      gkyl_gyrokinetic_app_species_orbit_average_stat(app, "ion");
    struct gkyl_update_status update = gkyl_gyrokinetic_update(app, dt);
    if (stat.enabled || stat.prepared || !update.success ||
        fabs(update.dt_actual - dt) > 1.0e-12 * dt) {
      gkyl_gyrokinetic_app_cout(app, stderr, "Both-disabled regression failed.\n");
      result = 1;
    }
    time += update.dt_actual;
    gkyl_gyrokinetic_app_calc_integrated_mom(app, time);

    // Re-enabling orbit relaxation explicitly rebuilds its frozen cache.
    options.enabled = true;
    gkyl_gyrokinetic_app_reset_species_orbit_average(app, "ion", options);
    stat = gkyl_gyrokinetic_app_species_orbit_average_stat(app, "ion");
    if (!stat.enabled || stat.prepared) {
      result = 1;
    }
    update = gkyl_gyrokinetic_update(app, dt);
    stat = gkyl_gyrokinetic_app_species_orbit_average_stat(app, "ion");
    if (!update.success || fabs(update.dt_actual - dt) > 1.0e-12 * dt || !stat.prepared ||
        stat.num_selected_cells == 0) {
      gkyl_gyrokinetic_app_cout(app, stderr, "Orbit reset/rebuild regression failed.\n");
      result = 1;
    }
    time += update.dt_actual;
    gkyl_gyrokinetic_app_calc_integrated_mom(app, time);
    gkyl_gyrokinetic_app_cout(app, stdout, "Orbit reset rebuilt with dt=%.12g.\n", update.dt_actual);

    // Preparation must reject invalid parameters before starting a time step.
    options.frequency = -1.0;
    gkyl_gyrokinetic_app_reset_species_orbit_average(app, "ion", options);
    update = gkyl_gyrokinetic_update(app, dt);
    stat = gkyl_gyrokinetic_app_species_orbit_average_stat(app, "ion");
    if (update.success || update.dt_actual != 0.0 || stat.prepared) {
      gkyl_gyrokinetic_app_cout(app, stderr, "Invalid orbit preparation was accepted.\n");
      result = 1;
    }
    gkyl_gyrokinetic_app_calc_integrated_mom(app, time);
    options.enabled = false;
    gkyl_gyrokinetic_app_reset_species_orbit_average(app, "ion", options);
  }
  gkyl_gyrokinetic_app_write_integrated_mom(app);
  gkyl_comm_barrier(comm);
  int rank = 0;
  gkyl_comm_get_rank(comm, &rank);
  if (!result && num_steps > 0 && rank == 0 &&
      !orbit_check_moments(app, args.app_name, num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stderr, "OAP integrated particle/energy conservation failed.\n");
    result = 1;
  }
  int global_result = result;
  gkyl_comm_allreduce_host(comm, GKYL_INT, GKYL_MAX, 1, &result, &global_result);
  gkyl_gyrokinetic_app_release(app);
  if (!global_result && num_steps > 0 && !orbit_check_collisionless_independence(&inp)) {
    result = 1;
  }
  gkyl_comm_allreduce_host(comm, GKYL_INT, GKYL_MAX, 1, &result, &global_result);
  if (!global_result && num_steps > 0 && !orbit_check_boltzmann_initialization(&inp)) {
    result = 1;
  }
  gkyl_comm_allreduce_host(comm, GKYL_INT, GKYL_MAX, 1, &result, &global_result);
  gkyl_gyrokinetic_comms_release(comm);
#ifdef GKYL_HAVE_MPI
  if (args.use_mpi) {
    MPI_Finalize();
  }
#endif
  return global_result;
}
