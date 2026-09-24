#include <gkyl_amr_core.h>
#include <gkyl_amr_block_priv.h>
#include <gkyl_amr_block_coupled_priv.h>
#include <gkyl_amr_patch_priv.h>
#include <gkyl_amr_patch_coupled_priv.h>

void
euler_mixture1d_run_single(int argc, char **argv, struct euler_mixture1d_single_init* init)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  int base_Nx = init->base_Nx;
  int ref_factor = init->ref_factor;

  double coarse_x1 = init->coarse_x1;
  double coarse_x2 = init->coarse_x2;

  double refined_x1 = init->refined_x1;
  double refined_x2 = init->refined_x2;

  evalf_t eval = init->eval;
  int num_species = init->num_species;
  double *gas_gamma_s = init->gas_gamma_s;

  char euler_mixture_output[64];
  strcpy(euler_mixture_output, init->euler_mixture_output);

  bool low_order_flux = init->low_order_flux;
  int num_frames = init->num_frames;

  double cfl_frac = init->cfl_frac;
  double t_end = init->t_end;
  double dt_failure_tol = init->dt_failure_tol;
  int num_failures_max = init->num_failures_max;

  int ndim = 1;
  int num_patches = 3;
  int Nx = base_Nx;

  struct euler_patch_data mesh_pdata[num_patches];
  struct gkyl_job_pool *mesh_job_pool = gkyl_thread_pool_new(app_args.num_threads);

  gkyl_rect_grid_init(&mesh_pdata[0].grid, 1, (double []) { refined_x1 }, (double []) { refined_x2 }, (int []) { Nx * ref_factor } );
  
  gkyl_rect_grid_init(&mesh_pdata[1].grid, 1, (double []) { coarse_x1 }, (double []) { refined_x1 }, (int []) { Nx } );
  gkyl_rect_grid_init(&mesh_pdata[2].grid, 1, (double []) { refined_x2 }, (double []) { coarse_x2 }, (int []) { Nx } );
  
  for (int i = 0; i < num_patches; i++) {
    mesh_pdata[i].fv_proj = gkyl_fv_proj_new(&mesh_pdata[i].grid, 1, 4 + (2 * num_species), eval, 0);
  }

  for (int i = 0; i < num_patches; i++) {
    gkyl_create_grid_ranges(&mesh_pdata[i].grid, (int []) { 2 }, &mesh_pdata[i].ext_range, &mesh_pdata[i].range );
    mesh_pdata[i].geom = gkyl_wave_geom_new(&mesh_pdata[i].grid, &mesh_pdata[i].ext_range, 0, 0, false);
  }

  for (int i = 0; i < num_patches; i++) {
    mesh_pdata[i].euler = gkyl_wv_euler_mixture_new(num_species, gas_gamma_s, app_args.use_gpu);

    mesh_pdata[i].slvr[0] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
        .grid = &mesh_pdata[i].grid,
        .equation = mesh_pdata[i].euler,
        .limiter = GKYL_MONOTONIZED_CENTERED,
        .num_up_dirs = 1,
        .update_dirs = { 0 },
        .cfl = cfl_frac,
        .geom = mesh_pdata[i].geom,
      }
    );
  }

  struct gkyl_block_topo *ptopo = create_patch_topo();

  for (int i = 0; i < num_patches; i++) {
    euler_mixture_patch_bc_updaters_init(mesh_pdata[i].euler, &mesh_pdata[i], &ptopo->conn[i]);
  }

  for (int i = 0; i < num_patches; i++) {
    mesh_pdata[i].fdup = gkyl_array_new(GKYL_DOUBLE, 4 + (2 * num_species), mesh_pdata[i].ext_range.volume);

    for (int d = 0; d < ndim + 1; d++) {
      mesh_pdata[i].f[d] = gkyl_array_new(GKYL_DOUBLE, 4 + (2 * num_species), mesh_pdata[i].ext_range.volume);
    }
  }

#ifdef AMR_USETHREADS
  for (int i = 0; i < num_patches; i++) {
    gkyl_job_pool_add_work(mesh_job_pool, euler_init_job_func_patch, &mesh_pdata[i]);
  }
  gkyl_job_pool_wait(mesh_job_pool);
#else
  for (int i = 0; i < num_patches; i++) {
    euler_init_job_func_patch(&mesh_pdata[i]);
  }
#endif

  char amr0[64];
  snprintf(amr0, 64, "%s_0", euler_mixture_output);
  euler_write_sol_patch(amr0, num_patches, mesh_pdata);

  double coarse_t_curr = 0.0;
  double fine_t_curr = 0.0;
  double coarse_dt = euler_max_dt_patch(num_patches, mesh_pdata);

  double fine_dt = (1.0 / ref_factor) * coarse_dt;

  struct sim_stats stats = { };

  struct timespec tm_start = gkyl_wall_clock();

  long coarse_step = 1;
  long num_steps = app_args.num_steps;

  double io_trigger = t_end / num_frames;

  double dt_init = -1.0;
  int num_failures = 0;

  while ((coarse_t_curr < t_end) && (coarse_step <= num_steps)) {
    printf("Taking coarse (level 0) time-step %ld at t = %g; ", coarse_step, coarse_t_curr);
    struct gkyl_update_status coarse_status = euler_update_patch(mesh_job_pool, ptopo, mesh_pdata, coarse_t_curr, coarse_dt, &stats);
    printf(" dt = %g\n", coarse_status.dt_actual);

    if (!coarse_status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }

    for (long fine_step = 1; fine_step < ref_factor + 1; fine_step++) {
      printf("   Taking fine (level 1) time-step %ld at t = %g", fine_step, fine_t_curr);
      printf(" dt = %g\n", (1.0 / ref_factor) * coarse_status.dt_actual);

      fine_t_curr += (1.0 / ref_factor) * coarse_status.dt_actual;
      fine_dt = (1.0 / ref_factor) * coarse_status.dt_suggested;
    }

    for (int i = 1; i < num_frames; i++) {
      if (coarse_t_curr < (i * io_trigger) && (coarse_t_curr + coarse_status.dt_actual) > (i * io_trigger)) {
        char buf[64];
        snprintf(buf, 64, "%s_%d", euler_mixture_output, i);

        euler_write_sol_patch(buf, num_patches, mesh_pdata);
      }
    }

    coarse_t_curr += coarse_status.dt_actual;
    coarse_dt = coarse_status.dt_suggested;

    if (dt_init < 0.0) {
      dt_init = coarse_status.dt_actual;
    }
    else if (coarse_status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      printf("WARNING: Time-step dt = %g", coarse_status.dt_actual);
      printf(" is below %g*dt_init ...", dt_failure_tol);
      printf(" num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        printf("ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        printf("%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    coarse_step += 1;
  }

  double tm_total_sec = gkyl_time_diff_now_sec(tm_start);

  char buf[64];
  snprintf(buf, 64, "%s_%d", euler_mixture_output, num_frames);

  euler_write_sol_patch(buf, num_patches, mesh_pdata);

  printf("\n");
  printf("Number of update calls %ld\n", (coarse_step - 1));
  printf("Number of failed time-steps %d\n", stats.nfail);
  printf("Total updates took %g secs\n", tm_total_sec);

  for (int i = 0; i < num_patches; i++) {
    gkyl_fv_proj_release(mesh_pdata[i].fv_proj);
    gkyl_wv_eqn_release(mesh_pdata[i].euler);
    euler_patch_bc_updaters_release(&mesh_pdata[i]);
    gkyl_wave_geom_release(mesh_pdata[i].geom);

    gkyl_wave_prop_release(mesh_pdata[i].slvr[0]);
    gkyl_array_release(mesh_pdata[i].fdup);
    gkyl_array_release(mesh_pdata[i].f[0]);
  }

  gkyl_block_topo_release(ptopo);
  gkyl_job_pool_release(mesh_job_pool);
}

void
euler_mixture1d_run_double(int argc, char **argv, struct euler_mixture1d_double_init* init)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  int base_Nx = init->base_Nx;
  int ref_factor1 = init->ref_factor1;
  int ref_factor2 = init->ref_factor2;

  double coarse_x1 = init->coarse_x1;
  double coarse_x2 = init->coarse_x2;

  double intermediate_x1 = init->intermediate_x1;
  double intermediate_x2 = init->intermediate_x2;

  double refined_x1 = init->refined_x1;
  double refined_x2 = init->refined_x2;

  evalf_t eval = init->eval;
  int num_species = init->num_species;
  double *gas_gamma_s = init->gas_gamma_s;

  char euler_mixture_output[64];
  strcpy(euler_mixture_output, init->euler_mixture_output);

  bool low_order_flux = init->low_order_flux;
  int num_frames = init->num_frames;

  double cfl_frac = init->cfl_frac;
  double t_end = init->t_end;
  double dt_failure_tol = init->dt_failure_tol;
  int num_failures_max = init->num_failures_max;

  int ndim = 1;
  int num_patches = 5;
  int Nx = base_Nx;

  struct euler_patch_data mesh_pdata[num_patches];
  struct gkyl_job_pool *mesh_job_pool = gkyl_thread_pool_new(app_args.num_threads);

  gkyl_rect_grid_init(&mesh_pdata[0].grid, 1, (double []) { refined_x1 }, (double []) { refined_x2 }, (int []) { Nx * (ref_factor1 * ref_factor2) } );
  
  gkyl_rect_grid_init(&mesh_pdata[1].grid, 1, (double []) { intermediate_x1 }, (double []) { refined_x1 }, (int []) { Nx * ref_factor1 } );
  gkyl_rect_grid_init(&mesh_pdata[2].grid, 1, (double []) { refined_x2 }, (double []) { intermediate_x2 }, (int []) { Nx * ref_factor1 } );

  gkyl_rect_grid_init(&mesh_pdata[3].grid, 1, (double []) { coarse_x1 }, (double []) { intermediate_x1 }, (int []) { Nx } );
  gkyl_rect_grid_init(&mesh_pdata[4].grid, 1, (double []) { intermediate_x2 }, (double []) { coarse_x2 }, (int []) { Nx } );
  
  for (int i = 0; i < num_patches; i++) {
    mesh_pdata[i].fv_proj = gkyl_fv_proj_new(&mesh_pdata[i].grid, 1, 4 + (2 * num_species), eval, 0);
  }

  for (int i = 0; i < num_patches; i++) {
    gkyl_create_grid_ranges(&mesh_pdata[i].grid, (int []) { 2 }, &mesh_pdata[i].ext_range, &mesh_pdata[i].range );
    mesh_pdata[i].geom = gkyl_wave_geom_new(&mesh_pdata[i].grid, &mesh_pdata[i].ext_range, 0, 0, false);
  }

  for (int i = 0; i < num_patches; i++) {
    mesh_pdata[i].euler = gkyl_wv_euler_mixture_new(num_species, gas_gamma_s, app_args.use_gpu);

    mesh_pdata[i].slvr[0] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
        .grid = &mesh_pdata[i].grid,
        .equation = mesh_pdata[i].euler,
        .limiter = GKYL_MONOTONIZED_CENTERED,
        .num_up_dirs = 1,
        .update_dirs = { 0 },
        .cfl = cfl_frac,
        .geom = mesh_pdata[i].geom,
      }
    );
  }

  struct gkyl_block_topo *ptopo = create_nested_patch_topo();

  for (int i = 0; i < num_patches; i++) {
    euler_mixture_nested_patch_bc_updaters_init(mesh_pdata[i].euler, &mesh_pdata[i], &ptopo->conn[i]);
  }

  for (int i = 0; i < num_patches; i++) {
    mesh_pdata[i].fdup = gkyl_array_new(GKYL_DOUBLE, 4 + (2 * num_species), mesh_pdata[i].ext_range.volume);

    for (int d = 0; d < ndim + 1; d++) {
      mesh_pdata[i].f[d] = gkyl_array_new(GKYL_DOUBLE, 4 + (2 * num_species), mesh_pdata[i].ext_range.volume);
    }
  }

#ifdef AMR_USETHREADS
  for (int i = 0; i < num_patches; i++) {
    gkyl_job_pool_add_work(mesh_job_pool, euler_init_job_func_patch, &mesh_pdata[i]);
  }
  gkyl_job_pool_wait(mesh_job_pool);
#else
  for (int i = 0; i < num_patches; i++) {
    euler_init_job_func_patch(&mesh_pdata[i]);
  }
#endif

  char amr0[64];
  snprintf(amr0, 64, "%s_0", euler_mixture_output);
  euler_write_sol_patch(amr0, num_patches, mesh_pdata);

  double coarse_t_curr = 0.0;
  double intermediate_t_curr = 0.0;
  double fine_t_curr = 0.0;
  double coarse_dt = euler_max_dt_patch(num_patches, mesh_pdata);

  double intermediate_dt = (1.0 / ref_factor1) * coarse_dt;
  double fine_dt = (1.0 / (ref_factor1 * ref_factor2)) * coarse_dt;

  struct sim_stats stats = { };

  struct timespec tm_start = gkyl_wall_clock();

  long coarse_step = 1;
  long num_steps = app_args.num_steps;

  double io_trigger = t_end / num_frames;

  double dt_init = -1.0;
  int num_failures = 0;

  while ((coarse_t_curr < t_end) && (coarse_step <= num_steps)) {
    printf("Taking coarse (level 0) time-step %ld at t = %g; ", coarse_step, coarse_t_curr);
    struct gkyl_update_status coarse_status = euler_update_patch(mesh_job_pool, ptopo, mesh_pdata, coarse_t_curr, coarse_dt, &stats);
    printf(" dt = %g\n", coarse_status.dt_actual);

    if (!coarse_status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }

    for (long intermediate_step = 1; intermediate_step < ref_factor1 + 1; intermediate_step++) {
      printf("   Taking intermediate (level 1) time-step %ld at t = %g", intermediate_step, intermediate_t_curr);
      printf(" dt = %g\n", (1.0 / ref_factor1) * coarse_status.dt_actual);

      for (long fine_step = 1; fine_step < ref_factor2 + 1; fine_step++) {
        printf("      Taking fine (level 2) time-step %ld at t = %g", fine_step, fine_t_curr);
        printf(" dt = %g\n", (1.0 / (ref_factor1 * ref_factor2)) * coarse_status.dt_actual);

        fine_t_curr += (1.0 / (ref_factor1 * ref_factor2)) * coarse_status.dt_actual;
        fine_dt = (1.0 / (ref_factor1 * ref_factor2)) * coarse_status.dt_suggested;
      }

      intermediate_t_curr += (1.0 / ref_factor1) * coarse_status.dt_actual;
      intermediate_dt = (1.0 / ref_factor1) * coarse_status.dt_suggested;
    }

    for (int i = 1; i < num_frames; i++) {
      if (coarse_t_curr < (i * io_trigger) && (coarse_t_curr + coarse_status.dt_actual) > (i * io_trigger)) {
        char buf[64];
        snprintf(buf, 64, "%s_%d", euler_mixture_output, i);

        euler_write_sol_patch(buf, num_patches, mesh_pdata);
      }
    }

    coarse_t_curr += coarse_status.dt_actual;
    coarse_dt = coarse_status.dt_suggested;

    if (dt_init < 0.0) {
      dt_init = coarse_status.dt_actual;
    }
    else if (coarse_status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      printf("WARNING: Time-step dt = %g", coarse_status.dt_actual);
      printf(" is below %g*dt_init ...", dt_failure_tol);
      printf(" num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        printf("ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        printf("%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    coarse_step += 1;
  }

  double tm_total_sec = gkyl_time_diff_now_sec(tm_start);

  char buf[64];
  snprintf(buf, 64, "%s_%d", euler_mixture_output, num_frames);

  euler_write_sol_patch(buf, num_patches, mesh_pdata);

  printf("\n");
  printf("Number of update calls %ld\n", (coarse_step - 1));
  printf("Number of failed time-steps %d\n", stats.nfail);
  printf("Total updates took %g secs\n", tm_total_sec);

  for (int i = 0; i < num_patches; i++) {
    gkyl_fv_proj_release(mesh_pdata[i].fv_proj);
    gkyl_wv_eqn_release(mesh_pdata[i].euler);
    euler_patch_bc_updaters_release(&mesh_pdata[i]);
    gkyl_wave_geom_release(mesh_pdata[i].geom);

    gkyl_wave_prop_release(mesh_pdata[i].slvr[0]);
    gkyl_array_release(mesh_pdata[i].fdup);
    gkyl_array_release(mesh_pdata[i].f[0]);
  }

  gkyl_block_topo_release(ptopo);
  gkyl_job_pool_release(mesh_job_pool);
}

void
euler_mixture2d_run_single(int argc, char **argv, struct euler_mixture2d_single_init* init)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  int base_Nx = init->base_Nx;
  int base_Ny = init->base_Ny;
  int ref_factor = init->ref_factor;

  double coarse_x1 = init->coarse_x1;
  double coarse_y1 = init->coarse_y1;
  double coarse_x2 = init->coarse_x2;
  double coarse_y2 = init->coarse_y2;

  double refined_x1 = init->refined_x1;
  double refined_y1 = init->refined_y1;
  double refined_x2 = init->refined_x2;
  double refined_y2 = init->refined_y2;

  evalf_t eval = init->eval;
  int num_species = init->num_species;
  double *gas_gamma_s = init->gas_gamma_s;

  bool copy_x = init->copy_x;
  bool copy_y = init->copy_y;
  
  bool wall_x = init->wall_x;
  bool wall_y = init->wall_y;

  char euler_mixture_output[64];
  strcpy(euler_mixture_output, init->euler_mixture_output);
  
  bool low_order_flux = init->low_order_flux;
  int num_frames = init->num_frames;

  double cfl_frac = init->cfl_frac;
  double t_end = init->t_end;
  double dt_failure_tol = init->dt_failure_tol;
  int num_failures_max = init->num_failures_max;

  int ndim = 2;
  int num_blocks = 9;
  int Nx = base_Nx;
  int Ny = base_Ny;

  struct euler_block_data mesh_bdata[num_blocks];
  struct gkyl_job_pool *mesh_job_pool = gkyl_thread_pool_new(app_args.num_threads);

  gkyl_rect_grid_init(&mesh_bdata[0].grid, 2, (double []) { refined_x1, refined_y1 }, (double []) { refined_x2, refined_y2 },
    (int []) { Nx * ref_factor, Ny * ref_factor });

  gkyl_rect_grid_init(&mesh_bdata[1].grid, 2, (double []) { coarse_x1, refined_y2 }, (double []) { refined_x1, coarse_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[2].grid, 2, (double []) { refined_x1, refined_y2 }, (double []) { refined_x2, coarse_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[3].grid, 2, (double []) { refined_x2, refined_y2 }, (double []) { coarse_x2, coarse_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[4].grid, 2, (double []) { coarse_x1, refined_y1 }, (double []) { refined_x1, refined_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[5].grid, 2, (double []) { refined_x2, refined_y1 }, (double []) { coarse_x2, refined_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[6].grid, 2, (double []) { coarse_x1, coarse_y1 }, (double []) { refined_x1, refined_y1 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[7].grid, 2, (double []) { refined_x1, coarse_y1 }, (double []) { refined_x2, refined_y1 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[8].grid, 2, (double []) { refined_x2, coarse_y1 }, (double []) { coarse_x2, refined_y1 },
    (int []) { Nx, Ny });

  for (int i = 0; i < num_blocks; i++) {
    mesh_bdata[i].fv_proj = gkyl_fv_proj_new(&mesh_bdata[i].grid, 2, 4 + (2 * num_species), eval, 0);
  }
  
  for (int i = 0; i < num_blocks; i++) {
    gkyl_create_grid_ranges(&mesh_bdata[i].grid, (int []) { 2, 2 }, &mesh_bdata[i].ext_range, &mesh_bdata[i].range);
    mesh_bdata[i].geom = gkyl_wave_geom_new(&mesh_bdata[i].grid, &mesh_bdata[i].ext_range, 0, 0, false);

    mesh_bdata[i].copy_x = copy_x;
    mesh_bdata[i].copy_y = copy_y;

    mesh_bdata[i].wall_x = wall_x;
    mesh_bdata[i].wall_y = wall_y;
  }

  for (int i = 0; i < num_blocks; i++) {
    mesh_bdata[i].euler = gkyl_wv_euler_mixture_new(num_species, gas_gamma_s, app_args.use_gpu);

    for (int d = 0; d < ndim; d++) {
      mesh_bdata[i].slvr[d] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
          .grid = &mesh_bdata[i].grid,
          .equation = mesh_bdata[i].euler,
          .limiter = GKYL_MONOTONIZED_CENTERED,
          .num_up_dirs = 1,
          .update_dirs = { d },
          .cfl = cfl_frac,
          .geom = mesh_bdata[i].geom,
        }
      );
    }
  }

  struct gkyl_block_topo *btopo = create_block_topo();

  for (int i = 0; i < num_blocks; i++) {
    euler_mixture_block_bc_updaters_init(mesh_bdata[i].euler, &mesh_bdata[i], &btopo->conn[i]);
  }

  for (int i = 0; i < num_blocks; i++) {
    mesh_bdata[i].fdup = gkyl_array_new(GKYL_DOUBLE, 4 + (2 * num_species), mesh_bdata[i].ext_range.volume);

    for (int d = 0; d < ndim + 1; d++) {
      mesh_bdata[i].f[d] = gkyl_array_new(GKYL_DOUBLE, 4 + (2 * num_species), mesh_bdata[i].ext_range.volume);
    }
  }

#ifdef AMR_USETHREADS
  for (int i = 0; i < num_blocks; i++) {
    gkyl_job_pool_add_work(mesh_job_pool, euler_init_job_func_block, &mesh_bdata[i]);
  }
  gkyl_job_pool_wait(mesh_job_pool);
#else
  for (int i = 0; i < num_blocks; i++) {
    euler_init_job_func_block(&mesh_bdata[i]);
  }
#endif

  char amr0[64];
  snprintf(amr0, 64, "%s_0", euler_mixture_output);
  euler_write_sol_block(amr0, num_blocks, mesh_bdata);

  double coarse_t_curr = 0.0;
  double fine_t_curr = 0.0;
  double coarse_dt = euler_max_dt_block(num_blocks, mesh_bdata);

  double fine_dt = (1.0 / ref_factor) * coarse_dt;

  struct sim_stats stats = { };

  struct timespec tm_start = gkyl_wall_clock();

  long coarse_step = 1;
  long num_steps = app_args.num_steps;

  double io_trigger = t_end / num_frames;

  double dt_init = -1.0;
  int num_failures = 0;

  while ((coarse_t_curr < t_end) && (coarse_step <= num_steps)) {
    printf("Taking coarse (level 0) time-step %ld at t = %g; ", coarse_step, coarse_t_curr);
    struct gkyl_update_status coarse_status = euler_update_block(mesh_job_pool, btopo, mesh_bdata, coarse_t_curr, coarse_dt, &stats);
    printf(" dt = %g\n", coarse_status.dt_actual);

    if (!coarse_status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }

    for (long fine_step = 1; fine_step < ref_factor + 1; fine_step++) {
      printf("   Taking fine (level 1) time-step %ld at t = %g", fine_step, fine_t_curr);
      printf(" dt = %g\n", (1.0 / ref_factor) * coarse_status.dt_actual);

      fine_t_curr += (1.0 / ref_factor) * coarse_status.dt_actual;
      fine_dt = (1.0 / ref_factor) * coarse_status.dt_suggested;
    }

    for (int i = 1; i < num_frames; i++) {
      if (coarse_t_curr < (i * io_trigger) && (coarse_t_curr + coarse_status.dt_actual) > (i * io_trigger)) {
        char buf[64];
        snprintf(buf, 64, "%s_%d", euler_mixture_output, i);

        euler_write_sol_block(buf, num_blocks, mesh_bdata);
      }
    }

    coarse_t_curr += coarse_status.dt_actual;
    coarse_dt = coarse_status.dt_suggested;

    if (dt_init < 0.0) {
      dt_init = coarse_status.dt_actual;
    }
    else if (coarse_status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      printf("WARNING: Time-step dt = %g", coarse_status.dt_actual);
      printf(" is below %g*dt_init ...", dt_failure_tol);
      printf(" num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        printf("ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        printf("%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    coarse_step += 1;
  }

  double tm_total_sec = gkyl_time_diff_now_sec(tm_start);

  char buf[64];
  snprintf(buf, 64, "%s_%d", euler_mixture_output, num_frames);

  euler_write_sol_block(buf, num_blocks, mesh_bdata);

  printf("\n");
  printf("Number of update calls %ld\n", (coarse_step - 1));
  printf("Number of failed time-steps %d\n", stats.nfail);
  printf("Total updates took %g secs\n", tm_total_sec);

  for (int i = 0; i < num_blocks; i++) {
    gkyl_fv_proj_release(mesh_bdata[i].fv_proj);
    gkyl_wv_eqn_release(mesh_bdata[i].euler);
    euler_block_bc_updaters_release(&mesh_bdata[i]);
    gkyl_wave_geom_release(mesh_bdata[i].geom);

    for (int d = 0; d < ndim; d++) {
      gkyl_wave_prop_release(mesh_bdata[i].slvr[d]);
    }
    
    gkyl_array_release(mesh_bdata[i].fdup);
    
    for (int d = 0; d < ndim; d++) {
      gkyl_array_release(mesh_bdata[i].f[d]);
    }
  }

  gkyl_block_topo_release(btopo);
  gkyl_job_pool_release(mesh_job_pool);
}

void
euler_mixture2d_run_double(int argc, char **argv, struct euler_mixture2d_double_init* init)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  int base_Nx = init->base_Nx;
  int base_Ny = init->base_Ny;
  int ref_factor1 = init->ref_factor1;
  int ref_factor2 = init->ref_factor2;

  double coarse_x1 = init->coarse_x1;
  double coarse_y1 = init->coarse_y1;
  double coarse_x2 = init->coarse_x2;
  double coarse_y2 = init->coarse_y2;

  double intermediate_x1 = init->intermediate_x1;
  double intermediate_y1 = init->intermediate_y1;
  double intermediate_x2 = init->intermediate_x2;
  double intermediate_y2 = init->intermediate_y2;

  double refined_x1 = init->refined_x1;
  double refined_y1 = init->refined_y1;
  double refined_x2 = init->refined_x2;
  double refined_y2 = init->refined_y2;

  evalf_t eval = init->eval;
  int num_species = init->num_species;
  double *gas_gamma_s = init->gas_gamma_s;

  bool copy_x = init->copy_x;
  bool copy_y = init->copy_y;

  bool wall_x = init->wall_x;
  bool wall_y = init->wall_y;

  char euler_mixture_output[64];
  strcpy(euler_mixture_output, init->euler_mixture_output);
  
  bool low_order_flux = init->low_order_flux;
  int num_frames = init->num_frames;

  double cfl_frac = init->cfl_frac;
  double t_end = init->t_end;
  double dt_failure_tol = init->dt_failure_tol;
  int num_failures_max = init->num_failures_max;

  int ndim = 2;
  int num_blocks = 25;
  int Nx = base_Nx;
  int Ny = base_Ny;

  struct euler_block_data mesh_bdata[num_blocks];
  struct gkyl_job_pool *mesh_job_pool = gkyl_thread_pool_new(app_args.num_threads);

  gkyl_rect_grid_init(&mesh_bdata[0].grid, 2, (double []) { refined_x1, refined_y1 }, (double []) { refined_x2, refined_y2 },
    (int []) { Nx * (ref_factor1 * ref_factor2), Ny * (ref_factor1 * ref_factor2) });

  gkyl_rect_grid_init(&mesh_bdata[1].grid, 2, (double []) { intermediate_x1, refined_y2 }, (double []) { refined_x1, intermediate_y2 },
    (int []) { Nx * ref_factor1, Ny * ref_factor1 });
  gkyl_rect_grid_init(&mesh_bdata[2].grid, 2, (double []) { refined_x1, refined_y2 }, (double []) { refined_x2, intermediate_y2 },
    (int []) { Nx * ref_factor1, Ny * ref_factor1 });
  gkyl_rect_grid_init(&mesh_bdata[3].grid, 2, (double []) { refined_x2, refined_y2 }, (double []) { intermediate_x2, intermediate_y2 },
    (int []) { Nx * ref_factor1, Ny * ref_factor1 });
  gkyl_rect_grid_init(&mesh_bdata[4].grid, 2, (double []) { intermediate_x1, refined_y1 }, (double []) { refined_x1, refined_y2 },
    (int []) { Nx * ref_factor1, Ny * ref_factor1 });
  gkyl_rect_grid_init(&mesh_bdata[5].grid, 2, (double []) { refined_x2, refined_y1 }, (double []) { intermediate_x2, refined_y2 },
    (int []) { Nx * ref_factor1, Ny * ref_factor1 });
  gkyl_rect_grid_init(&mesh_bdata[6].grid, 2, (double []) { intermediate_x1, intermediate_y1 }, (double []) { refined_x1, refined_y1 },
    (int []) { Nx * ref_factor1, Ny * ref_factor1 });
  gkyl_rect_grid_init(&mesh_bdata[7].grid, 2, (double []) { refined_x1, intermediate_y1 }, (double []) { refined_x2, refined_y1 },
    (int []) { Nx * ref_factor1, Ny * ref_factor1 });
  gkyl_rect_grid_init(&mesh_bdata[8].grid, 2, (double []) { refined_x2, intermediate_y1 }, (double []) { intermediate_x2, refined_y1 },
    (int []) { Nx * ref_factor1, Ny * ref_factor1 });

  gkyl_rect_grid_init(&mesh_bdata[9].grid, 2, (double []) { coarse_x1, intermediate_y2 }, (double []) { intermediate_x1, coarse_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[10].grid, 2, (double []) { intermediate_x1, intermediate_y2 }, (double []) { refined_x1, coarse_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[11].grid, 2, (double []) { refined_x1, intermediate_y2 }, (double []) { refined_x2, coarse_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[12].grid, 2, (double []) { refined_x2, intermediate_y2 }, (double []) { intermediate_x2, coarse_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[13].grid, 2, (double []) { intermediate_x2, intermediate_y2 }, (double []) { coarse_x2, coarse_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[14].grid, 2, (double []) { coarse_x1, refined_y2 }, (double []) { intermediate_x1, intermediate_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[15].grid, 2, (double []) { intermediate_x2, refined_y2 }, (double []) { coarse_x2, intermediate_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[16].grid, 2, (double []) { coarse_x1, refined_y1 }, (double []) { intermediate_x1, refined_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[17].grid, 2, (double []) { intermediate_x2, refined_y1 }, (double []) { coarse_x2, refined_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[18].grid, 2, (double []) { coarse_x1, intermediate_y1 }, (double []) { intermediate_x1, refined_y1 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[19].grid, 2, (double []) { intermediate_x2, intermediate_y1 }, (double []) { coarse_x2, refined_y1 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[20].grid, 2, (double []) { coarse_x1, coarse_y1 }, (double []) { intermediate_x1, intermediate_y1 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[21].grid, 2, (double []) { intermediate_x1, coarse_y1 }, (double []) { refined_x1, intermediate_y1 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[22].grid, 2, (double []) { refined_x1, coarse_y1 }, (double []) { refined_x2, intermediate_y1 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[23].grid, 2, (double []) { refined_x2, coarse_y1 }, (double []) { intermediate_x2, intermediate_y1 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&mesh_bdata[24].grid, 2, (double []) { intermediate_x2, coarse_y1 }, (double []) { coarse_x2, intermediate_y1 },
    (int []) { Nx, Ny });

  for (int i = 0; i < num_blocks; i++) {
    mesh_bdata[i].fv_proj = gkyl_fv_proj_new(&mesh_bdata[i].grid, 2, 4 + (2 * num_species), eval, 0);
  }
  
  for (int i = 0; i < num_blocks; i++) {
    gkyl_create_grid_ranges(&mesh_bdata[i].grid, (int []) { 2, 2 }, &mesh_bdata[i].ext_range, &mesh_bdata[i].range);
    mesh_bdata[i].geom = gkyl_wave_geom_new(&mesh_bdata[i].grid, &mesh_bdata[i].ext_range, 0, 0, false);

    mesh_bdata[i].copy_x = copy_x;
    mesh_bdata[i].copy_y = copy_y;

    mesh_bdata[i].wall_x = wall_x;
    mesh_bdata[i].wall_y = wall_y;
  }

  for (int i = 0; i < num_blocks; i++) {
    mesh_bdata[i].euler = gkyl_wv_euler_mixture_new(num_species, gas_gamma_s, app_args.use_gpu);

    for (int d = 0; d < ndim; d++) {
      mesh_bdata[i].slvr[d] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
          .grid = &mesh_bdata[i].grid,
          .equation = mesh_bdata[i].euler,
          .limiter = GKYL_MONOTONIZED_CENTERED,
          .num_up_dirs = 1,
          .update_dirs = { d },
          .cfl = cfl_frac,
          .geom = mesh_bdata[i].geom,
        }
      );
    }
  }

  struct gkyl_block_topo *btopo = create_nested_block_topo();

  for (int i = 0; i < num_blocks; i++) {
    euler_mixture_nested_block_bc_updaters_init(mesh_bdata[i].euler, &mesh_bdata[i], &btopo->conn[i]);
  }

  for (int i = 0; i < num_blocks; i++) {
    mesh_bdata[i].fdup = gkyl_array_new(GKYL_DOUBLE, 4 + (2 * num_species), mesh_bdata[i].ext_range.volume);

    for (int d = 0; d < ndim + 1; d++) {
      mesh_bdata[i].f[d] = gkyl_array_new(GKYL_DOUBLE, 4 + (2 * num_species), mesh_bdata[i].ext_range.volume);
    }
  }

#ifdef AMR_USETHREADS
  for (int i = 0; i < num_blocks; i++) {
    gkyl_job_pool_add_work(mesh_job_pool, euler_init_job_func_block, &mesh_bdata[i]);
  }
  gkyl_job_pool_wait(mesh_job_pool);
#else
  for (int i = 0; i < num_blocks; i++) {
    euler_init_job_func_block(&mesh_bdata[i]);
  }
#endif

  char amr0[64];
  snprintf(amr0, 64, "%s_0", euler_mixture_output);
  euler_write_sol_block(amr0, num_blocks, mesh_bdata);

  double coarse_t_curr = 0.0;
  double intermediate_t_curr = 0.0;
  double fine_t_curr = 0.0;
  double coarse_dt = euler_max_dt_block(num_blocks, mesh_bdata);

  double intermediate_dt = (1.0 / ref_factor1) * coarse_dt;
  double fine_dt = (1.0 / (ref_factor1 * ref_factor2)) * coarse_dt;

  struct sim_stats stats = { };

  struct timespec tm_start = gkyl_wall_clock();

  long coarse_step = 1;
  long num_steps = app_args.num_steps;

  double io_trigger = t_end / num_frames;

  double dt_init = -1.0;
  int num_failures = 0;

  while ((coarse_t_curr < t_end) && (coarse_step <= num_steps)) {
    printf("Taking coarse (level 0) time-step %ld at t = %g; ", coarse_step, coarse_t_curr);
    struct gkyl_update_status coarse_status = euler_update_block(mesh_job_pool, btopo, mesh_bdata, coarse_t_curr, coarse_dt, &stats);
    printf(" dt = %g\n", coarse_status.dt_actual);

    if (!coarse_status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }

    for (long intermediate_step = 1; intermediate_step < ref_factor1 + 1; intermediate_step++) {
      printf("   Taking intermediate (level 1) time-step %ld at t = %g", intermediate_step, intermediate_t_curr);
      printf(" dt = %g\n", (1.0 / ref_factor1) * coarse_status.dt_actual);

      for (long fine_step = 1; fine_step < ref_factor2 + 1; fine_step++) {
        printf("      Taking fine (level 2) time-step %ld at t = %g", fine_step, fine_t_curr);
        printf(" dt = %g\n", (1.0 / (ref_factor1 * ref_factor2)) * coarse_status.dt_actual);

        fine_t_curr += (1.0 / (ref_factor1 * ref_factor2)) * coarse_status.dt_actual;
        fine_dt = (1.0 / (ref_factor1 * ref_factor2)) * coarse_status.dt_suggested;
      }

      intermediate_t_curr += (1.0 / ref_factor1) * coarse_status.dt_actual;
      intermediate_dt = (1.0 / ref_factor1) * coarse_status.dt_suggested;
    }

    for (int i = 1; i < num_frames; i++) {
      if (coarse_t_curr < (i * io_trigger) && (coarse_t_curr + coarse_status.dt_actual) > (i * io_trigger)) {
        char buf[64];
        snprintf(buf, 64, "%s_%d", euler_mixture_output, i);

        euler_write_sol_block(buf, num_blocks, mesh_bdata);
      }
    }

    coarse_t_curr += coarse_status.dt_actual;
    coarse_dt = coarse_status.dt_suggested;

    if (dt_init < 0.0) {
      dt_init = coarse_status.dt_actual;
    }
    else if (coarse_status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      printf("WARNING: Time-step dt = %g", coarse_status.dt_actual);
      printf(" is below %g*dt_init ...", dt_failure_tol);
      printf(" num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        printf("ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        printf("%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    coarse_step += 1;
  }

  double tm_total_sec = gkyl_time_diff_now_sec(tm_start);

  char buf[64];
  snprintf(buf, 64, "%s_%d", euler_mixture_output, num_frames);

  euler_write_sol_block(buf, num_blocks, mesh_bdata);

  printf("\n");
  printf("Number of update calls %ld\n", (coarse_step - 1));
  printf("Number of failed time-steps %d\n", stats.nfail);
  printf("Total updates took %g secs\n", tm_total_sec);

  for (int i = 0; i < num_blocks; i++) {
    gkyl_fv_proj_release(mesh_bdata[i].fv_proj);
    gkyl_wv_eqn_release(mesh_bdata[i].euler);
    euler_block_bc_updaters_release(&mesh_bdata[i]);
    gkyl_wave_geom_release(mesh_bdata[i].geom);

    for (int d = 0; d < ndim; d++) {
      gkyl_wave_prop_release(mesh_bdata[i].slvr[d]);
    }
    
    gkyl_array_release(mesh_bdata[i].fdup);
    
    for (int d = 0; d < ndim; d++) {
      gkyl_array_release(mesh_bdata[i].f[d]);
    }
  }

  gkyl_block_topo_release(btopo);
  gkyl_job_pool_release(mesh_job_pool);
}