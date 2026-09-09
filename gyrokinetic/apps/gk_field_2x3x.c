#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_position_map.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_comm_io.h>
#include <gkyl_dg_interpolate.h>

#include <assert.h>
#include <float.h>
#include <time.h>

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
  // Fill upper parallel boundary ghost with skin boundary value.
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
gk_field_2x3x_add_TS_updaters(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf, struct gkyl_poisson_bc *poisson_bcs)
{
  // Allocation ranges and updaters for TS field solve.

  // Parallel smoother for the charge density.
  enum gkyl_fem_parproj_bc_type fem_parproj_bc_rho_core;
  enum gkyl_fem_parproj_bc_type fem_parproj_bc_phi_core;

  if (app->cdim == 2) {
    fem_parproj_bc_rho_core = GKYL_FEM_PARPROJ_PERIODIC;
    fem_parproj_bc_phi_core = GKYL_FEM_PARPROJ_PERIODIC;

    gkf->fem_projection_par_rho_func = gk_field_fem_projection_par_rho_ts_2x;
    gkf->fem_projection_par_phi_func = gk_field_fem_projection_par_phi_ts_2x;
  }
  else if (app->cdim == 3) {
    // Here fem_parproj_bc_rho is not actually relevant because we don't use gkf->fem_parproj_rho.
    fem_parproj_bc_rho_core = 0;
    fem_parproj_bc_phi_core = GKYL_FEM_PARPROJ_DIRICHLET_GHOST;

    gkf->fem_projection_par_rho_func = gk_field_fem_projection_par;
    gkf->fem_projection_par_phi_func = gk_field_fem_projection_par_phi_ts_3x;

    int num_ghost[] = {1, 1, 1};
    int par_dir = app->cdim-1; // Parallel direction index.

    // TS BC updater for lower edge.
    struct gkyl_bc_twistshift_inp T_LU_lo = {
      .bc_dir = par_dir,
      .shift_dir = 1, // y shift.
      .shear_dir = 0, // shift varies with x.
      .edge = GKYL_LOWER_EDGE,
      .cdim = app->cdim,
      .bcdir_ext_update_r = &app->global_par_ext,
      .num_ghost = num_ghost, // one ghost per config direction
      .basis = &app->basis,
      .grid = &app->grid,
      .use_gpu = app->use_gpu,
      .upsample_factor = app->ts_upsample_factor,
      .filter_half_width = app->ts_filter_half_width,
      .filter_cutoff_wavelength = app->ts_filter_cutoff_wavelength,
    };
    if (app->gk_geom->geometry_id == GKYL_GEOMETRY_TOKAMAK)
      T_LU_lo.shift_dg = app->delta_ts_x_lo;
    else {
      T_LU_lo.shift_func     = app->gk_geom->parallel_lower_bc_shift_func;
      T_LU_lo.shift_func_ctx = app->gk_geom->parallel_lower_bc_shift_ctx;
    }
    gkf->bc_ts_lo = gkyl_bc_twistshift_inew(&T_LU_lo);

    // TS BC updater for upper edge.
    struct gkyl_bc_twistshift_inp T_UL_up = {
      .bc_dir = par_dir,
      .shift_dir = 1, // y shift.
      .shear_dir = 0, // shift varies with x.
      .edge = GKYL_UPPER_EDGE,
      .cdim = app->cdim,
      .bcdir_ext_update_r = &app->global_par_ext,
      .num_ghost = num_ghost, // one ghost per config direction
      .basis = &app->basis,
      .grid = &app->grid,
      .use_gpu = app->use_gpu,
      .upsample_factor = app->ts_upsample_factor,
      .filter_half_width = app->ts_filter_half_width,
      .filter_cutoff_wavelength = app->ts_filter_cutoff_wavelength,
    };
    if (app->gk_geom->geometry_id == GKYL_GEOMETRY_TOKAMAK)
      T_UL_up.shift_dg = app->delta_ts_x_up;
    else {
      T_UL_up.shift_func     = app->gk_geom->parallel_upper_bc_shift_func;
      T_UL_up.shift_func_ctx = app->gk_geom->parallel_upper_bc_shift_ctx;
    }
    gkf->bc_ts_up = gkyl_bc_twistshift_inew(&T_UL_up);

    long buff_sz = app->global_lower_ghost[par_dir].volume;
    gkf->bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);

    gkf->gfss_bc_op_core_up = gkyl_bc_basic_gyrokinetic_new(par_dir, GKYL_UPPER_EDGE, GKYL_BC_GK_FIELD_BOUNDARY_VALUE,
      app->basis_on_dev, &app->global_upper_skin[par_dir], &app->global_upper_ghost[par_dir],
      app->basis.num_basis, app->cdim, app->use_gpu);

    gkf->gfss_bc_op_core_lo = gkyl_bc_basic_gyrokinetic_new(par_dir, GKYL_LOWER_EDGE, GKYL_BC_GK_FIELD_BOUNDARY_VALUE,
      app->basis_on_dev, &app->global_lower_skin[par_dir], &app->global_lower_ghost[par_dir],
      app->basis.num_basis, app->cdim, app->use_gpu);

  }

  // Parallel smoother for the charge density.
  gkf->fem_parproj_rho_core = gkyl_fem_parproj_new(&app->global, &app->grid, &app->basis,
    fem_parproj_bc_rho_core, 0, 0, 0, app->use_gpu);

  // Fill bias line list for fem_parproj_phi.
  gk_field_2x3x_fill_fem_parproj_bias_lines(app, gkf, poisson_bcs);
    
  // Parallel smoother for the potential.
  gkf->fem_parproj_phi_core = gkyl_fem_parproj_new(&app->global, &app->grid, &app->basis,
    fem_parproj_bc_phi_core, &gkf->fem_parproj_bias_line_list, 0, 0, app->use_gpu);
}

static void
gk_field_2x3x_add_IWL_updaters(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf, struct gkyl_poisson_bc *poisson_bcs)
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

    gkf->fem_projection_par_rho_func = gk_field_fem_projection_par_rho_iwl_2x;
    gkf->fem_projection_par_phi_func = gk_field_fem_projection_par_phi_iwl_2x;
  }
  else if (app->cdim == 3) {
    // Here fem_parproj_bc_rho is not actually relevant because we don't use gkf->fem_parproj_rho.
    fem_parproj_bc_rho_core = 0;
    fem_parproj_bc_rho_sol  = 0;
    fem_parproj_bc_phi_core = GKYL_FEM_PARPROJ_DIRICHLET_GHOST;
    fem_parproj_bc_phi_sol  = GKYL_FEM_PARPROJ_DIRICHLET_SKIN;

    gkf->fem_projection_par_rho_func = gk_field_fem_projection_par;
    gkf->fem_projection_par_phi_func = gk_field_fem_projection_par_phi_iwl_3x;

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
      .upsample_factor = app->ts_upsample_factor,
      .filter_half_width = app->ts_filter_half_width,
      .filter_cutoff_wavelength = app->ts_filter_cutoff_wavelength,
    };
    if (app->gk_geom->geometry_id == GKYL_GEOMETRY_TOKAMAK)
      T_LU_lo.shift_dg = app->delta_ts_x_lo;
    else {
      T_LU_lo.shift_func     = app->gk_geom->parallel_lower_bc_shift_func;
      T_LU_lo.shift_func_ctx = app->gk_geom->parallel_lower_bc_shift_ctx;
    }
    gkf->bc_ts_lo = gkyl_bc_twistshift_inew(&T_LU_lo);

    long buff_sz = GKYL_MAX2(app->global_lower_ghost_par_sol.volume, app->global_lower_ghost_par_core.volume);
    gkf->bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);

    gkf->gfss_bc_op_core_up = gkyl_bc_basic_gyrokinetic_new(par_dir, GKYL_UPPER_EDGE, GKYL_BC_GK_FIELD_BOUNDARY_VALUE,
      app->basis_on_dev, &app->global_upper_skin_par_core, &app->global_upper_ghost_par_core,
      app->basis.num_basis, app->cdim, app->use_gpu);
  }

  // Parallel smoother for the charge density.
  gkf->fem_parproj_rho_core = gkyl_fem_parproj_new(&app->global_core, &app->grid, &app->basis,
    fem_parproj_bc_rho_core, 0, 0, 0, app->use_gpu);
  gkf->fem_parproj_rho_sol = gkyl_fem_parproj_new(&app->global_sol, &app->grid, &app->basis,
    fem_parproj_bc_rho_sol, 0, 0, 0, app->use_gpu);

  // Fill bias line list for fem_parproj_phi.
  gk_field_2x3x_fill_fem_parproj_bias_lines(app, gkf, poisson_bcs);
    
  // Parallel smoother for the potential.
  gkf->fem_parproj_phi_core = gkyl_fem_parproj_new(&app->global_core, &app->grid, &app->basis,
    fem_parproj_bc_phi_core, &gkf->fem_parproj_bias_line_list, 0, 0, app->use_gpu);
  gkf->fem_parproj_phi_sol = gkyl_fem_parproj_new(&app->global_sol, &app->grid, &app->basis,
    fem_parproj_bc_phi_sol, &gkf->fem_parproj_bias_line_list, 0, 0, app->use_gpu);
    
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
gk_field_fem_release_2x3x(const gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_array_release(f->rho_c);
  gkyl_array_release(f->rho_c_global_dg);
  gkyl_array_release(f->rho_c_global_smooth);
  gkyl_array_release(f->phi_fem);
  gkyl_array_release(f->phi_smooth);

  if (f->gkfield_id == GKYL_GK_FIELD_EM) {
    gkyl_array_release(f->apar_fem);
    gkyl_array_release(f->apardot_fem);
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

  // Allocate electromagnetic arrays if needed.
  if (f->gkfield_id == GKYL_GK_FIELD_EM) {
    f->apar_fem = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->apardot_fem = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
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
