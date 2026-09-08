#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_position_map.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_comm_io.h>
#include <gkyl_array_average.h>

#include <assert.h>
#include <float.h>
#include <time.h>

static void
gk_field_fem_projection_par_none(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Do nothing.
}

static void
gk_field_fem_projection_par_rho_ts_2x(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Project a DG charge density onto the parallel FEM basis to make it
  // continuous along z using different BCs in the core and SOL.

  // Gather the DG array into a global (in z) array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_dg, field->rho_c_global_dg);

  // Smooth the DG array.
  gkyl_fem_parproj_set_rhs(field->fem_parproj_rho_core, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_rho_core, field->phi_fem);

  // Copy global, continuous FEM array to a local array.
  gkyl_array_copy_range_to_range(arr_fem, field->phi_fem, &app->local, &field->global_sub_range);
}

static void
gk_field_fem_projection_par_phi_ts_2x(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Project a DG potential onto the parallel FEM basis to make it
  // continuous along z using different BCs in the core and SOL.

  // Gather the DG array into a global (in z) array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_dg, field->rho_c_global_dg);

  // Smooth the the DG array.
  gkyl_fem_parproj_set_rhs(field->fem_parproj_phi_core, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_phi_core, field->phi_fem);

  // Copy global, continuous FEM array to a local array.
  gkyl_array_copy_range_to_range(arr_fem, field->phi_fem, &app->local, &field->global_sub_range);
}

static void
gk_field_fem_projection_par_phi_ts_3x(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Project a DG field onto the parallel FEM basis to make it
  // continuous along z (or to solve a Poisson equation in 1x),
  // using twistshift BCs in the parallel direction.

  // Gather the DG array into a global (in z) array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_dg, field->rho_c_global_dg);

  // Apply TS BC in the core lower parallel boundary, and
  // fill core upper parallel boundary ghost with skin boundary value.
  int par_dir = app->cdim-1; // Parallel direction index.
  gkyl_array_copy_range_to_range(field->rho_c_global_dg, field->rho_c_global_dg,
    &app->global_lower_ghost[par_dir], &app->global_upper_skin[par_dir]);
  gkyl_bc_twistshift_advance(field->bc_ts_lo, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_bc_basic_gyrokinetic_advance(field->gfss_bc_op_core_up, field->bc_buffer, field->rho_c_global_dg);

  // Smooth the the DG array.
  gkyl_fem_parproj_set_rhs(field->fem_parproj_phi_core, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_phi_core, field->phi_fem);

  // Copy global, continuous FEM array to a local array.
  gkyl_array_copy_range_to_range(arr_fem, field->phi_fem, &app->local, &field->global_sub_range);
}

static void
gk_field_fem_projection_par_rho_iwl_2x(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Project a DG charge density onto the parallel FEM basis to make it
  // continuous along z using different BCs in the core and SOL.

  // Gather the DG array into a global (in z) array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_dg, field->rho_c_global_dg);

  // Smooth the DG array.
  gkyl_fem_parproj_set_rhs(field->fem_parproj_rho_core, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_rho_core, field->phi_fem);
  gkyl_fem_parproj_set_rhs(field->fem_parproj_rho_sol, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_rho_sol, field->phi_fem);

  // Copy global, continuous FEM array to a local array.
  gkyl_array_copy_range_to_range(arr_fem, field->phi_fem, &app->local, &field->global_sub_range);
}

static void
gk_field_fem_projection_par_phi_iwl_2x(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Project a DG potential onto the parallel FEM basis to make it
  // continuous along z using different BCs in the core and SOL.

  // Gather the DG array into a global (in z) array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_dg, field->rho_c_global_dg);

  // Smooth the the DG array.
  gkyl_fem_parproj_set_rhs(field->fem_parproj_phi_core, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_phi_core, field->phi_fem);
  gkyl_fem_parproj_set_rhs(field->fem_parproj_phi_sol, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_phi_sol, field->phi_fem);

  // Copy global, continuous FEM array to a local array.
  gkyl_array_copy_range_to_range(arr_fem, field->phi_fem, &app->local, &field->global_sub_range);
}

static void
gk_field_fem_projection_par_phi_iwl_3x(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Project a DG field onto the parallel FEM basis to make it
  // continuous along z (or to solve a Poisson equation in 1x),
  // using different BCs in the core and SOL.

  // Gather the DG array into a global (in z) array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_dg, field->rho_c_global_dg);

  // Apply TS BC in the core lower parallel boundary, and
  // fill core upper parallel boundary ghost with skin boundary value.
  gkyl_array_copy_range_to_range(field->rho_c_global_dg, field->rho_c_global_dg,
    &app->global_lower_ghost_par_core, &app->global_upper_skin_par_core);
  gkyl_bc_twistshift_advance(field->bc_ts_lo, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_bc_basic_gyrokinetic_advance(field->gfss_bc_op_core_up, field->bc_buffer, field->rho_c_global_dg);

  // Smooth the the DG array.
  gkyl_fem_parproj_set_rhs(field->fem_parproj_phi_core, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_phi_core, field->phi_fem);
  gkyl_fem_parproj_set_rhs(field->fem_parproj_phi_sol, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_phi_sol, field->phi_fem);

  // Copy global, continuous FEM array to a local array.
  gkyl_array_copy_range_to_range(arr_fem, field->phi_fem, &app->local, &field->global_sub_range);
}

static void
gk_field_2x3x_fill_fem_parproj_bias_lines(struct gkyl_gyrokinetic_app *app, struct gk_field *f, struct gkyl_poisson_bc *poisson_bcs)
{
  // Create a bias line list that includes the perpendicular BCs if
  // they are Dirichlet, and the bias lines from the input file.
  
  // Create a temporary biased line list.
  int par_dir = app->cdim-1; // Parallel direction index.
  int num_bias_line = 0;
  int bl_idx = 0;
  int num_bias_line_in = 0;
  if (f->info.bias_line_list)
    num_bias_line_in = f->info.bias_line_list->num_bias_line;

  int num_bc_bias = 4; // Number of biases for Dirichlet BCs below.
  size_t bl_sz = (num_bc_bias+num_bias_line_in) * sizeof(struct gkyl_poisson_bias_line);
  struct gkyl_poisson_bias_line *bias_lines_buff = gkyl_malloc(bl_sz);
  if (poisson_bcs->lo_type[0] == GKYL_POISSON_DIRICHLET) {
    // psi increases towards SOL.
    struct gkyl_poisson_bias_line *bl;
    // Core x boundary at lower z boundary.
    bl = &(bias_lines_buff[bl_idx]);
    bl->perp_dirs[0] = 0;
    bl->perp_dirs[1] = par_dir;
    bl->perp_coords[0] = app->grid.lower[0];
    bl->perp_coords[1] = app->grid.lower[par_dir];
    bl->val = poisson_bcs->lo_value[0].v[0];
    num_bias_line++;
    bl_idx++;
    // Core x boundary at upper z boundary.
    bl = &(bias_lines_buff[bl_idx]);
    bl->perp_dirs[0] = 0;
    bl->perp_dirs[1] = par_dir;
    bl->perp_coords[0] = app->grid.lower[0];
    bl->perp_coords[1] = app->grid.upper[par_dir];
    bl->val = poisson_bcs->lo_value[0].v[0];
    num_bias_line++;
    bl_idx++;
  } 
  if (poisson_bcs->up_type[0] == GKYL_POISSON_DIRICHLET) {
    // psi increases towards core.
    struct gkyl_poisson_bias_line *bl;
    // Core x boundary at lower z boundary.
    bl = &(bias_lines_buff[bl_idx]);
    bl->perp_dirs[0] = 0;
    bl->perp_dirs[1] = par_dir;
    bl->perp_coords[0] = app->grid.upper[0];
    bl->perp_coords[1] = app->grid.lower[par_dir];
    bl->val = poisson_bcs->up_value[0].v[0];
    num_bias_line++;
    bl_idx++;
    // Core x boundary at upper z boundary.
    bl = &(bias_lines_buff[bl_idx]);
    bl->perp_dirs[0] = 0;
    bl->perp_dirs[1] = par_dir;
    bl->perp_coords[0] = app->grid.upper[0];
    bl->perp_coords[1] = app->grid.upper[par_dir];
    bl->val = poisson_bcs->up_value[0].v[0];
    num_bias_line++;
    bl_idx++;
  } 
  for (int i=0; i<num_bias_line_in; i++) {
    struct gkyl_poisson_bias_line *bl = &(bias_lines_buff[bl_idx]);
    struct gkyl_poisson_bias_line *bl_inp = &(f->info.bias_line_list->bl[i]);
    bl->perp_dirs[0]   = bl_inp->perp_dirs[0]  ;
    bl->perp_dirs[1]   = bl_inp->perp_dirs[1]  ;
    bl->perp_coords[0] = bl_inp->perp_coords[0];
    bl->perp_coords[1] = bl_inp->perp_coords[1];
    bl->val = bl_inp->val;
    num_bias_line++;
    bl_idx++;
  }
  
  // Copy temporary bias line list into app.
  bl_sz = GKYL_MAX2(1,num_bias_line) * sizeof(struct gkyl_poisson_bias_line); // max avoids allocating 0 memory.
  f->fem_parproj_bias_line_list.num_bias_line = num_bias_line;
  f->fem_parproj_bias_line_list.bl = gkyl_malloc(bl_sz);
  if (num_bias_line)
    memcpy(f->fem_parproj_bias_line_list.bl, bias_lines_buff, bl_sz);

  gkyl_free(bias_lines_buff);
}

static void
gk_field_2x3x_add_TS_updaters(struct gkyl_gyrokinetic_app *app, struct gk_field *f, struct gkyl_poisson_bc *poisson_bcs)
{
  // Allocation ranges and updaters for TS field solve.

  // Parallel smoother for the charge density.
  enum gkyl_fem_parproj_bc_type fem_parproj_bc_rho_core;
  enum gkyl_fem_parproj_bc_type fem_parproj_bc_phi_core;

  if (app->cdim == 2) {
    fem_parproj_bc_rho_core = GKYL_FEM_PARPROJ_PERIODIC;
    fem_parproj_bc_phi_core = GKYL_FEM_PARPROJ_PERIODIC;

    f->fem_projection_par_rho_func = gk_field_fem_projection_par_rho_ts_2x;
    f->fem_projection_par_phi_func = gk_field_fem_projection_par_phi_ts_2x;
  }
  else if (app->cdim == 3) {
    // Here fem_parproj_bc_rho is not actually relevant because we don't use f->fem_parproj_rho.
    fem_parproj_bc_rho_core = 0;
    fem_parproj_bc_phi_core = GKYL_FEM_PARPROJ_DIRICHLET_GHOST;

    f->fem_projection_par_rho_func = gk_field_fem_projection_par;
    f->fem_projection_par_phi_func = gk_field_fem_projection_par_phi_ts_3x;

    int par_dir = app->cdim-1; // Parallel direction index.
    int ghost[] = {1, 1, 1};
    // TS BC updater for lower edge.
    struct gkyl_bc_twistshift_inp T_LU_lo = {
      .bc_dir = par_dir,
      .shift_dir = 1, // y shift.
      .shear_dir = 0, // shift varies with x.
      .edge = GKYL_LOWER_EDGE,
      .cdim = app->cdim,
      .bcdir_ext_update_r = &app->global_par_ext,
      .num_ghost = ghost, // one ghost per config direction
      .basis = &app->basis,
      .grid = &app->grid,
      .use_gpu = app->use_gpu,
    };
    if (app->gk_geom->geometry_id == GKYL_GEOMETRY_TOKAMAK)
      T_LU_lo.shift_dg = app->delta_ts_x_lo;
    else {
      T_LU_lo.shift_func     = app->gk_geom->parallel_lower_bc_shift_func;
      T_LU_lo.shift_func_ctx = app->gk_geom->parallel_lower_bc_shift_ctx;
    }
    f->bc_ts_lo = gkyl_bc_twistshift_inew(&T_LU_lo);

    // TS BC updater for upper edge.
    struct gkyl_bc_twistshift_inp T_UL_up = {
      .bc_dir = par_dir,
      .shift_dir = 1, // y shift.
      .shear_dir = 0, // shift varies with x.
      .edge = GKYL_UPPER_EDGE,
      .cdim = app->cdim,
      .bcdir_ext_update_r = &app->global_par_ext,
      .num_ghost = ghost, // one ghost per config direction
      .basis = &app->basis,
      .grid = &app->grid,
      .use_gpu = app->use_gpu,
    };
    if (app->gk_geom->geometry_id == GKYL_GEOMETRY_TOKAMAK)
      T_UL_up.shift_dg = app->delta_ts_x_up;
    else {
      T_UL_up.shift_func     = app->gk_geom->parallel_upper_bc_shift_func;
      T_UL_up.shift_func_ctx = app->gk_geom->parallel_upper_bc_shift_ctx;
    }
    f->bc_ts_up = gkyl_bc_twistshift_inew(&T_UL_up);

    long buff_sz = app->global_lower_ghost[par_dir].volume;
    f->bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);

    f->gfss_bc_op_core_up = gkyl_bc_basic_gyrokinetic_new(par_dir, GKYL_UPPER_EDGE, GKYL_BC_GK_FIELD_BOUNDARY_VALUE,
      app->basis_on_dev, &app->global_upper_skin[par_dir], &app->global_upper_ghost[par_dir],
      app->basis.num_basis, app->cdim, app->use_gpu);

    f->gfss_bc_op_core_lo = gkyl_bc_basic_gyrokinetic_new(par_dir, GKYL_LOWER_EDGE, GKYL_BC_GK_FIELD_BOUNDARY_VALUE,
      app->basis_on_dev, &app->global_lower_skin[par_dir], &app->global_lower_ghost[par_dir],
      app->basis.num_basis, app->cdim, app->use_gpu);

  }

  // Parallel smoother for the charge density.
  f->fem_parproj_rho_core = gkyl_fem_parproj_new(&app->global, &app->grid, &app->basis,
    fem_parproj_bc_rho_core, 0, 0, 0, app->use_gpu);

  // Fill bias line list for fem_parproj_phi.
  gk_field_2x3x_fill_fem_parproj_bias_lines(app, f, poisson_bcs);
    
  // Parallel smoother for the potential.
  f->fem_parproj_phi_core = gkyl_fem_parproj_new(&app->global, &app->grid, &app->basis,
    fem_parproj_bc_phi_core, &f->fem_parproj_bias_line_list, 0, 0, app->use_gpu);
}

static void
gk_field_2x3x_add_IWL_updaters(struct gkyl_gyrokinetic_app *app, struct gk_field *f, struct gkyl_poisson_bc *poisson_bcs)
{
  // Allocation ranges and updaters for IWL field solve.

  // Parallel smoother for the charge density.
  enum gkyl_fem_parproj_bc_type fem_parproj_bc_rho_core, fem_parproj_bc_rho_sol;
  enum gkyl_fem_parproj_bc_type fem_parproj_bc_phi_core, fem_parproj_bc_phi_sol;

  if (app->cdim == 2) {
    fem_parproj_bc_rho_core = GKYL_FEM_PARPROJ_PERIODIC;
    fem_parproj_bc_rho_sol  = GKYL_FEM_PARPROJ_NONE;
    fem_parproj_bc_phi_core = GKYL_FEM_PARPROJ_PERIODIC;
    fem_parproj_bc_phi_sol  = GKYL_FEM_PARPROJ_NONE;

    f->fem_projection_par_rho_func = gk_field_fem_projection_par_rho_iwl_2x;
    f->fem_projection_par_phi_func = gk_field_fem_projection_par_phi_iwl_2x;
  }
  else if (app->cdim == 3) {
    // Here fem_parproj_bc_rho is not actually relevant because we don't use f->fem_parproj_rho.
    fem_parproj_bc_rho_core = 0;
    fem_parproj_bc_rho_sol  = 0;
    fem_parproj_bc_phi_core = GKYL_FEM_PARPROJ_DIRICHLET_GHOST;
    fem_parproj_bc_phi_sol  = GKYL_FEM_PARPROJ_DIRICHLET_SKIN;

    f->fem_projection_par_rho_func = gk_field_fem_projection_par;
    f->fem_projection_par_phi_func = gk_field_fem_projection_par_phi_iwl_3x;

    int par_dir = app->cdim-1; // Parallel direction index.
    // TS BC updater for up to low TS for the lower edge. This sets ghost_L = T_LU(ghost_L).
    int ghost[] = {1, 1, 1};
    struct gkyl_bc_twistshift_inp T_LU_lo = {
      .bc_dir = par_dir,
      .shift_dir = 1, // y shift.
      .shear_dir = 0, // shift varies with x.
      .edge = GKYL_LOWER_EDGE,
      .cdim = app->cdim,
      .bcdir_ext_update_r = &app->global_par_ext_core,
      .num_ghost = ghost, // one ghost per config direction
      .basis = &app->basis,
      .grid = &app->grid,
      .use_gpu = app->use_gpu,
    };
    if (app->gk_geom->geometry_id == GKYL_GEOMETRY_TOKAMAK)
      T_LU_lo.shift_dg = app->delta_ts_x_lo;
    else {
      T_LU_lo.shift_func     = app->gk_geom->parallel_lower_bc_shift_func;
      T_LU_lo.shift_func_ctx = app->gk_geom->parallel_lower_bc_shift_ctx;
    }
    f->bc_ts_lo = gkyl_bc_twistshift_inew(&T_LU_lo);

    long buff_sz = GKYL_MAX2(app->global_lower_ghost_par_sol.volume, app->global_lower_ghost_par_core.volume);
    f->bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);

    f->gfss_bc_op_core_up = gkyl_bc_basic_gyrokinetic_new(par_dir, GKYL_UPPER_EDGE, GKYL_BC_GK_FIELD_BOUNDARY_VALUE,
      app->basis_on_dev, &app->global_upper_skin_par_core, &app->global_upper_ghost_par_core,
      app->basis.num_basis, app->cdim, app->use_gpu);
  }

  // Parallel smoother for the charge density.
  f->fem_parproj_rho_core = gkyl_fem_parproj_new(&app->global_core, &app->grid, &app->basis,
    fem_parproj_bc_rho_core, 0, 0, 0, app->use_gpu);
  f->fem_parproj_rho_sol = gkyl_fem_parproj_new(&app->global_sol, &app->grid, &app->basis,
    fem_parproj_bc_rho_sol, 0, 0, 0, app->use_gpu);

  // Fill bias line list for fem_parproj_phi.
  gk_field_2x3x_fill_fem_parproj_bias_lines(app, f, poisson_bcs);
    
  // Parallel smoother for the potential.
  f->fem_parproj_phi_core = gkyl_fem_parproj_new(&app->global_core, &app->grid, &app->basis,
    fem_parproj_bc_phi_core, &f->fem_parproj_bias_line_list, 0, 0, app->use_gpu);
  f->fem_parproj_phi_sol = gkyl_fem_parproj_new(&app->global_sol, &app->grid, &app->basis,
    fem_parproj_bc_phi_sol, &f->fem_parproj_bias_line_list, 0, 0, app->use_gpu);
    
}

static void
gk_field_rhs_poisson_perp_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Smooth the charge density along z.
  field->fem_projection_par_rho_func(app, field, field->rho_c, field->rho_c);

  // Solve the Poisson equation.
  gkyl_fem_poisson_perp_set_rhs(field->fem_poisson_perp, field->rho_c);
  gkyl_fem_poisson_perp_solve(field->fem_poisson_perp, field->phi_smooth);

  // Smooth the potential along z.
  field->fem_projection_par_phi_func(app, field, field->phi_smooth, field->phi_smooth);

  // Finish the Poisson solve with FLR effects.
  field->invert_flr(app, field, field->phi_smooth);
}

static void
gk_field_ohm_solve(struct gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *out){
  struct timespec wst = gkyl_wall_clock();
  
  gkyl_fem_poisson_perp_update_lhs(field->fem_apardot_solver, field->lapWeightAmpere, field->dApartdtSlvr_kSq);
  gkyl_fem_poisson_perp_set_rhs(field->fem_apardot_solver, field->currentDensdot);
  gkyl_fem_poisson_perp_solve(field->fem_apardot_solver, out);

  field->invert_flr(app, field, field->apardot);

  field->remove_em_zonal_func(app, field, field->apardot);

  app->stat.field_apar_solve_tm += gkyl_time_diff_now_sec(wst);
}

static void 
gk_field_ampere_solve_enabled(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *out){
  struct timespec wst = gkyl_wall_clock();

  gkyl_fem_poisson_perp_set_rhs(field->fem_apar_solver, field->currentDens);
  gkyl_fem_poisson_perp_solve(field->fem_apar_solver, out);

  // Smooth along z.
  // gk_field_fem_projection_par(app, field, field->apar, field->apar_smooth_aux);
  // gkyl_array_copy_range(field->apar, field->apar_smooth_aux, &app->local_ext);

  field->invert_flr(app, field, out);

  app->stat.field_apar_solve_tm += gkyl_time_diff_now_sec(wst);
}

static void 
gk_field_ampere_solve_none(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *out){
  // Do nothing.
}

static void
gk_field_em_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_field *field, const struct gkyl_array *f_in[],  struct gkyl_array *rhs_in[])
{
  gk_field_accumulate_current_dens_dot(app, field, rhs_in);
  gk_field_accumulate_ohms_kSq(app, field, f_in);
  gk_field_ohm_solve(app, field, field->apardot);
}

static void
gk_field_em_rhs_none(gkyl_gyrokinetic_app *app, struct gk_field *field, const struct gkyl_array *f_in[],  struct gkyl_array *rhs_in[])
{
  // Do nothing.
}

static void
gk_field_remove_em_zonal_enabled(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *field_array)
{
  gkyl_array_average_advance(field->fs_avg_op, field_array, field->fs_avg);
  // Hack the conf to phase space dg op to do the 1x substraction over the entire conf space.
  gkyl_dg_mul_conf_phase_op_accumulate_range( &field->fs_avg_basis, &app->basis, 
    field_array, -1.0, field->fs_avg, field->fs_avg_conf_one, &field->fs_avg_range, &app->local);
}

static void
gk_field_em_zonal_component_none(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *field_array)
{
  // Do nothing.
}


static void
gk_field_fem_release_2x3x(const gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_array_release(f->rho_c);
  gkyl_array_release(f->rho_c_global_dg);
  gkyl_array_release(f->rho_c_global_smooth);
  gkyl_array_release(f->phi_fem);
  gkyl_array_release(f->phi_smooth);
  gkyl_array_release(f->apar);
  gkyl_array_release(f->apardot);
  gkyl_array_release(f->amperesol);
  gkyl_array_release(f->apar_curr);
  gkyl_array_release(f->apar1);
  gkyl_array_release(f->aparnew);

  if (f->is_em) {
    gkyl_array_release(f->apar_smooth_aux);
    gkyl_array_release(f->currentDens);
    gkyl_array_release(f->currentDensdot);
    gkyl_array_release(f->lapWeightAmpere);
    gkyl_array_release(f->dApartdtSlvr_kSq);
    gkyl_fem_parproj_release(f->fem_apar_parproj);
    gkyl_fem_poisson_perp_release(f->fem_apar_solver);
    gkyl_fem_poisson_perp_release(f->fem_apardot_solver);
    if (app->use_gpu) {
      gkyl_array_release(f->apar_host);
      gkyl_array_release(f->apardot_host);
      gkyl_array_release(f->amperesol_host);
    }
    if (f->info.remove_em_zonal) {
     gkyl_array_average_release(f->fs_avg_op);
     gkyl_array_release(f->fs_avg);
     gkyl_array_release(f->fs_avg_conf_one);
     if (app->use_gpu) {
      gkyl_cu_free(f->fs_avg_subdir);
     } else {
      gkyl_free(f->fs_avg_subdir);
     }
    }
  }

  if (app->use_gpu) {
    gkyl_array_release(f->phi_host);
  }

  gkyl_array_release(f->epsilon);

  gkyl_fem_poisson_perp_release(f->fem_poisson_perp);
  if (f->is_dirichletvar) {
    gkyl_array_release(f->phi_bc);
  }
  
  gkyl_fem_parproj_release(f->fem_parproj);

  gkyl_array_integrate_release(f->calc_em_energy);

  if (app->gk_geom->has_LCFS) {
    // Release updaters for solve with a LCFS/separatrix.
    gkyl_free(f->fem_parproj_bias_line_list.bl);
    gkyl_fem_parproj_release(f->fem_parproj_rho_core);
    gkyl_fem_parproj_release(f->fem_parproj_phi_core);
    gkyl_fem_parproj_release(f->fem_parproj_rho_sol);
    gkyl_fem_parproj_release(f->fem_parproj_phi_sol);

    if (app->cdim == 3) {
      gkyl_bc_twistshift_release(f->bc_ts_lo);
      gkyl_bc_basic_gyrokinetic_release(f->gfss_bc_op_core_up);
      gkyl_array_release(f->bc_buffer);
    }
  }
  else if (f->bc_par_phi == GKYL_BC_GK_FIELD_TWISTSHIFT) {
    // Release TS updaters.
    gkyl_free(f->fem_parproj_bias_line_list.bl);
    gkyl_fem_parproj_release(f->fem_parproj_rho_core);
    gkyl_fem_parproj_release(f->fem_parproj_phi_core);

    if (app->cdim == 3) {
      gkyl_bc_twistshift_release(f->bc_ts_lo);
      gkyl_bc_twistshift_release(f->bc_ts_up);
      gkyl_bc_basic_gyrokinetic_release(f->gfss_bc_op_core_up);
      gkyl_bc_basic_gyrokinetic_release(f->gfss_bc_op_core_lo);
      gkyl_array_release(f->bc_buffer);
    }
  }
  
  
  if (f->use_flr) {
    gk_field_flr_release(app, f);
  }
}

void
gk_field_fem_new_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  // Create global subrange we'll copy the field solver solution from (into local).
  gkyl_sub_range_intersect(&f->global_sub_range, &app->global, &app->local);

  // Allocate arrays for charge density.
  f->rho_c = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->rho_c_global_dg = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  f->rho_c_global_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);

  // Allocate arrays for electrostatic potential.
  f->phi_fem = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  f->phi_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->apar = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->apar_curr = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->apar1 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->aparnew = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->apardot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->amperesol = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  // Setup electromagnetic variables if needed.
  if (f->is_em) {
    f->apar_host = f->apar;
    f->apardot_host = f->apardot;
    f->amperesol_host = f->amperesol;
    if (app->use_gpu) {
      f->apar_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      f->apardot_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      f->amperesol_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    }

    f->apar_smooth_aux = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->currentDens = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->currentDensdot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->lapWeightAmpere = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);
    f->dApartdtSlvr_kSq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  }

  // Allocate phi_host for I/O.
  f->phi_host = f->phi_smooth;
  if (app->use_gpu) {
    f->phi_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  }

  if (f->gkfield_id == GKYL_GK_FIELD_ADIABATIC) {
    f->accumulate_rhoc_func = gk_field_accumulate_rho_c_adiabatic;
  } else {
    f->accumulate_rhoc_func = gk_field_accumulate_rho_c_poisson;
  }

  double polarization_weight = 0.0;
  double polarization_bmag = f->info.polarization_bmag ? f->info.polarization_bmag : app->bmag_ref;
  // Linearized polarization density
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    polarization_weight += s->info.polarization_density*s->info.mass/pow(polarization_bmag,2);
  }
  // Allocate array for the polarization weight times geometric coefficients.
  f->epsilon = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);
  
  // Initialize the polarization weight.
  struct gkyl_array *Jgij[3] = {app->gk_geom->geo_int.gxxj, app->gk_geom->geo_int.gxyj, app->gk_geom->geo_int.gyyj};
  for (int i=0; i<app->cdim-2/app->cdim; i++) {
    gkyl_array_set_offset(f->epsilon, polarization_weight, Jgij[i], i*app->basis.num_basis);
  }

  if (f->is_em) {
    // Compute the weights in front of the laplacian operator (metric x 1/mu0)
    for (int i=0; i<app->cdim-2/app->cdim; i++) {
      gkyl_array_set_offset(f->lapWeightAmpere, 1.0/f->info.mu0, Jgij[i], i*app->basis.num_basis);
    }
    gkyl_array_set(f->apar_energy_fac, -0.5, f->lapWeightAmpere);
  }

  bool bc_is_np[GKYL_MAX_CDIM]; // Is the BC in this direction non-periodic?
  for (int d=0; d<app->cdim; ++d) bc_is_np[d] = true;
  for (int d=0; d<app->num_periodic_dir; ++d) {
    bc_is_np[app->periodic_dirs[d]] = false;
  }

  // Translate input file BCs into Poisson BCs.
  struct gkyl_poisson_bc poisson_bcs = { };
  for (int d=0; d<app->cdim-1; d++) {
    if (bc_is_np[d]) {
      struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2*app->cdim, d, GKYL_LOWER_EDGE);
      if (bc_lo != 0) {
        poisson_bcs.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(bc_lo->type);
        for (int i=0; i<3; i++) {
          poisson_bcs.lo_value[d].v[i] = bc_lo->value[i];
        }
      }

      struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2*app->cdim, d, GKYL_UPPER_EDGE);
      if (bc_up != 0) {
        poisson_bcs.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(bc_up->type);
        for (int i=0; i<3; i++) {
          poisson_bcs.up_value[d].v[i] = bc_up->value[i];
        }
      }
    } else {
      poisson_bcs.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
      poisson_bcs.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
    }
  }

  // Initialize the Poisson solver.
  f->fem_poisson_perp = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->basis,
    &poisson_bcs, f->info.bias_line_list, f->epsilon, NULL, app->use_gpu);

  f->phi_bc = 0;
  f->is_dirichletvar = false;
  for (int i=0; i<2*app->cdim; i++) {
    f->is_dirichletvar = f->is_dirichletvar ||
                          (f->info.poisson_bcs[i].type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING ||
                           f->info.poisson_bcs[i].type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING);
  }

  if (f->is_dirichletvar) {
    // Project the spatially varying BC if the user specifies it.
    f->phi_bc = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
    struct gkyl_array *phi_bc_ho = mkarr(false, f->phi_bc->ncomp, f->phi_bc->size);

    for (int d=0; d<app->cdim; d++) {
      struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2*app->cdim, d, GKYL_LOWER_EDGE);
      if (bc_lo != 0) {
        if (bc_lo->type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING) {
          gkyl_eval_on_nodes *phibc_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis, 
            1, bc_lo->aux_profile, bc_lo->aux_ctx);
          gkyl_eval_on_nodes_advance(phibc_proj, 0.0, &app->local_lower_skin[d], phi_bc_ho);
          gkyl_eval_on_nodes_release(phibc_proj);
        }
      }
      struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2*app->cdim, d, GKYL_UPPER_EDGE);
      if (bc_up != 0) {
        if (bc_up->type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING) {
          gkyl_eval_on_nodes *phibc_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis, 
            1, bc_up->aux_profile, bc_up->aux_ctx);
          gkyl_eval_on_nodes_advance(phibc_proj, 0.0, &app->local_lower_skin[d], phi_bc_ho);
          gkyl_eval_on_nodes_release(phibc_proj);
        }
      }
    }
    gkyl_array_copy(f->phi_bc, phi_bc_ho);
    gkyl_array_release(phi_bc_ho);
  }

  // Potential smoothing (in z) updater
  enum gkyl_fem_parproj_bc_type fem_parproj_bc = GKYL_FEM_PARPROJ_NONE;
  for (int d=0; d<app->num_periodic_dir; ++d)
    if (app->periodic_dirs[d] == app->cdim-1) fem_parproj_bc = GKYL_FEM_PARPROJ_PERIODIC;

  f->fem_parproj = gkyl_fem_parproj_new(&app->global, &app->grid, &app->basis,
    fem_parproj_bc, 0, 0, 0, app->use_gpu);

  f->fem_projection_par_rho_func = gk_field_fem_projection_par;
  f->fem_projection_par_phi_func = gk_field_fem_projection_par;

  // Setup EM solvers.
  f->ampere_solve = gk_field_ampere_solve_none;
  f->em_rhs_func = gk_field_em_rhs_none;
  f->remove_em_zonal_func = gk_field_em_zonal_component_none;
  f->fem_projection_par_apar_func = gk_field_fem_projection_par_none;
  if (f->is_em) {
    // Translate input file BCs into Ampere BCs.
    for (int d=0; d<app->cdim-1; d++) {
      if (bc_is_np[d]) {
        struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(f->info.ampere_bcs, 2*app->cdim, d, GKYL_LOWER_EDGE);
        if (bc_lo != 0) {
          f->ampere_bcs.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(bc_lo->type);
          for (int i=0; i<3; i++)
            f->ampere_bcs.lo_value[d].v[i] = bc_lo->value[i];
        }

        struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(f->info.ampere_bcs, 2*app->cdim, d, GKYL_UPPER_EDGE);
        if (bc_up != 0) {
          f->ampere_bcs.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(bc_up->type);
          for (int i=0; i<3; i++)
            f->ampere_bcs.up_value[d].v[i] = bc_up->value[i];
        }
      } else {
        f->ampere_bcs.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
        f->ampere_bcs.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
      }
    }
    // FEM solver for Ampere equation (solved at IC only).
    f->fem_apar_solver = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->basis,
        &f->ampere_bcs, 0, f->lapWeightAmpere, NULL, app->use_gpu);
    // FEM solver for Ohm's law (evolves d(Apart)/dt).
    f->fem_apardot_solver = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->basis,
      &f->ampere_bcs, 0, f->lapWeightAmpere, f->dApartdtSlvr_kSq, app->use_gpu);

    // We have to advance the updater once to be able to update the kSq later.
    gkyl_fem_poisson_perp_set_rhs(f->fem_apardot_solver, f->currentDensdot);
    gkyl_fem_poisson_perp_solve(f->fem_apardot_solver, f->apardot);
    
    // FEM smoother for Aparallel.
    f->fem_parproj_ampere_bc = GKYL_FEM_PARPROJ_NONE;
    f->fem_apar_parproj = gkyl_fem_parproj_new(&app->global, &app->grid, 
       &app->basis, f->fem_parproj_ampere_bc, 0, 0, 0, app->use_gpu);

    f->ampere_solve = gk_field_ampere_solve_enabled;
    f->em_rhs_func = f->info.is_apar_static ? gk_field_em_rhs_none : gk_field_em_rhs_enabled;

    if (f->info.remove_em_zonal) {
      f->remove_em_zonal_func = gk_field_remove_em_zonal_enabled;
      // define the reduced range and basis for averaging along x only
      struct gkyl_rect_grid grid_x;
      gkyl_rect_grid_init(&grid_x, 1, &app->grid.lower[0], &app->grid.upper[0], &app->grid.cells[0]);
      int ghost_x[] = {1};
      gkyl_create_grid_ranges(&grid_x, ghost_x, &f->fs_avg_range_ext, &f->fs_avg_range);
      gkyl_cart_modal_serendip(&f->fs_avg_basis, 1, app->basis.poly_order);
      // create and run the array average updater to average y and z
      // Not that in 2x, it will only average over z.
      int avg_dim_yz[] = {0,1,1};
      struct gkyl_array_average_inp inp_avg_xyz_to_x = {
        .grid = &app->grid,
        .basis = app->basis,
        .basis_avg = f->fs_avg_basis,
        .local = &app->local,
        .local_avg = &f->fs_avg_range,
        .local_avg_ext = &f->fs_avg_range_ext,
        .weight = app->gk_geom->geo_int.jacobgeo,
        .avg_dim = avg_dim_yz,
        .use_gpu = app->use_gpu
      };
      int fs_avg_subdir[GKYL_MAX_CDIM];
      fs_avg_subdir[0] = 0; // Average along x only, so subdir is just 0.
      if (app->use_gpu) {
        f->fs_avg_subdir = gkyl_cu_malloc(GKYL_MAX_CDIM*sizeof(int));
        gkyl_cu_memcpy(f->fs_avg_subdir, fs_avg_subdir, GKYL_MAX_CDIM*sizeof(int), GKYL_CU_MEMCPY_H2D);
      } else {
        f->fs_avg_subdir = gkyl_malloc(GKYL_MAX_CDIM*sizeof(int));
        for (int i=0; i<GKYL_MAX_CDIM; i++)
          f->fs_avg_subdir[i] = fs_avg_subdir[i];
      }
      f->fs_avg_op = gkyl_array_average_inew(&inp_avg_xyz_to_x);
      f->fs_avg = mkarr(app->use_gpu, f->fs_avg_basis.num_basis, f->fs_avg_range_ext.volume);
      f->fs_avg_conf_one = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      double dg_norm = pow(sqrt(2.0), app->basis.ndim);
      gkyl_array_shiftc_range(f->fs_avg_conf_one, dg_norm, 0, &app->local);
    }
    f->fem_projection_par_apar_func = f->fem_projection_par_phi_func;
  }

  // Updater for field energy calculation.
  gkyl_array_set(f->es_energy_fac, 0.5, f->epsilon);
  f->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->basis, 
    1, GKYL_ARRAY_INTEGRATE_OP_EPS_GRADPERP_SQ, app->use_gpu);

  // Create operator needed for FLR effects.
  f->use_flr = false;
  f->invert_flr = gk_field_invert_flr_none;
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    if (s->info.flr.type) {
      f->use_flr = f->use_flr || s->info.flr.type;
    }
  }
  if (f->use_flr) {
    gk_field_flr_new(app, f);
  }

  f->bc_par_phi = 0;
  // Deterime if we need IWL or TWISTSHIFT BCs on phi fro the species BCs.
  for (int s=0; s<app->num_species; s++) {
    struct gk_species *gks = &app->species[s];
    for (int i = 0; i < 2*app->cdim; i++) {
      if ( gks->info.bcs[i].dir == app->cdim-1 && gks->info.bcs[i].type == GKYL_BC_GK_SPECIES_TWISTSHIFT ) {
        f->bc_par_phi = GKYL_BC_GK_FIELD_TWISTSHIFT;
        break;
      }
    }
    if (f->bc_par_phi)
      break;
  }

  if (app->gk_geom->has_LCFS) {
    // Updaters to enforce twist-and-shift and sheath BCs.
    gk_field_2x3x_add_IWL_updaters(app, f, &poisson_bcs);
  }
  else if (f->bc_par_phi == GKYL_BC_GK_FIELD_TWISTSHIFT) {
    // Updaters to enforce twist-and-shift BCs.
    gk_field_2x3x_add_TS_updaters(app, f, &poisson_bcs);
  }

  // Set the pointer to the function that computes phi.
  f->rhs_phi_func = gk_field_rhs_poisson_perp_2x3x;

  // Set pointer to function that releases memory.
  f->release_func = gk_field_fem_release_2x3x;
}
