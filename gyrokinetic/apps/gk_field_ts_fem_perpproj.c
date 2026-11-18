#include <gkyl_deflate_zsurf.h>
#include <gkyl_fem_proj.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>
#include <gkyl_translate_dim.h>

struct gk_field_ts_fem_perpproj {
  struct gkyl_rect_grid grid;
  struct gkyl_basis basis;
  struct gkyl_basis *boundary_basis;
  struct gkyl_range local, local_ext, ghost;
  struct gkyl_range sol_skin, core_skin;
  struct gkyl_array *boundary_values, *boundary_buffer;
  struct gkyl_deflate_zsurf *deflate_sol;
  struct gkyl_bc_basic_gyrokinetic *copy_lcfs;
  void (*prepare_boundary)(const gkyl_gyrokinetic_app *, struct gk_field_ts_fem_perpproj *, const struct gkyl_array *);
  struct gkyl_array *buff_dg, *buff_fem;
  struct gkyl_deflate_zsurf *deflate;
  struct gkyl_translate_dim *inflate;
  struct gkyl_fem_proj *projector;
};

// Apply the default projection only when the field uses twist-shift boundaries.
static bool
gk_field_ts_fem_perpproj_requested(const gkyl_gyrokinetic_app *app, const struct gk_field *field)
{
  return !app->core_parallel_bcs.ts_field_dg && (field->bc_par_phi == GKYL_BC_GK_FIELD_TWISTSHIFT ||
                                                 (app->cdim == 3 && app->gk_geom->has_LCFS));
}

static void
gk_field_ts_fem_perpproj_disabled(
  const gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *potential
)
{
}

void
gk_field_ts_fem_perpproj_check_input(const gkyl_gyrokinetic_app *app, const struct gk_field *field)
{
  if (!gk_field_ts_fem_perpproj_requested(app, field)) {
    return;
  }
  if (app->cdim != 3 || app->poly_order != 1 || app->basis.b_type != GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_exit("ts_field_dg=false requires a 3D p1 serendipity configuration basis.");
  }
  if (field->gkfield_id != GKYL_GK_FIELD_ES) {
    gkyl_exit("ts_field_dg=false requires an electrostatic field.");
  }
}

static void
gk_field_ts_fem_perpproj_validate(const gkyl_gyrokinetic_app *app, const struct gk_field *field)
{
  if (!app->gk_geom->has_LCFS && field->bc_par_phi != GKYL_BC_GK_FIELD_TWISTSHIFT) {
    gkyl_exit("ts_field_dg=false requires twist-shift field boundaries.");
  }
  if (field->use_flr) {
    gkyl_exit("ts_field_dg=false does not support FLR.");
  }
  if (app->gk_geom->has_LCFS && (!app->global_core.volume || !app->global_sol.volume)) {
    gkyl_exit("ts_field_dg=false requires nonempty core and SOL ranges in an IWL simulation.");
  }
  const struct gkyl_poisson_bc *bcs = &field->poisson_bcs;
  if (bcs->lo_type[1] != GKYL_POISSON_PERIODIC || bcs->up_type[1] != GKYL_POISSON_PERIODIC) {
    gkyl_exit("ts_field_dg=false requires periodic y boundary conditions.");
  }
  if ((bcs->lo_type[0] == GKYL_POISSON_PERIODIC) != (bcs->up_type[0] == GKYL_POISSON_PERIODIC)) {
    gkyl_exit("ts_field_dg=false requires paired periodic x boundary conditions.");
  }
  if (app->gk_geom->has_LCFS && bcs->lo_type[0] == GKYL_POISSON_PERIODIC) {
    gkyl_exit("ts_field_dg=false cannot use periodic x with an LCFS.");
  }
  for (int d = 0; d < 2; d++) {
    if (app->local.lower[d] != app->global.lower[d] ||
        app->local.upper[d] != app->global.upper[d]) {
      gkyl_exit("ts_field_dg=false supports decomposition only along z.");
    }
  }
}

// The deflated ghost potential already supplies the physical boundary data.
static void
gk_field_ts_fem_perpproj_boundary_disabled(
  const gkyl_gyrokinetic_app *app, struct gk_field_ts_fem_perpproj *proj,
  const struct gkyl_array *potential
)
{
}

// Match the lower physical SOL boundary at the LCFS. Evaluate the SOL skin's
// lower z face, then extend its LCFS value into the adjacent core boundary cell.
// This operates only on a separate boundary-data buffer, not on the potential.
static void
gk_field_ts_fem_perpproj_boundary_sol(
  const gkyl_gyrokinetic_app *app, struct gk_field_ts_fem_perpproj *proj,
  const struct gkyl_array *potential
)
{
  gkyl_deflate_zsurf_advance(
    proj->deflate_sol, app->global.lower[2], &app->global_ext, &proj->sol_skin, potential,
    proj->boundary_values, 1
  );
  gkyl_bc_basic_gyrokinetic_advance(proj->copy_lcfs, proj->boundary_buffer, proj->boundary_values);
}

// Allocate the optional varying-Dirichlet data at the core/SOL interface.
static void
gk_field_ts_fem_perpproj_init_sol(
  const gkyl_gyrokinetic_app *app, struct gk_field_ts_fem_perpproj *proj, int side
)
{
  int lower[] = {side ? proj->local.upper[0] : proj->local.lower[0], proj->local.lower[1]};
  int upper[] = {lower[0], proj->local.upper[1]};
  gkyl_sub_range_init(&proj->core_skin, &proj->local_ext, lower, upper);

  lower[0] += side ? 1 : -1;
  upper[0] = lower[0];
  gkyl_sub_range_init(&proj->sol_skin, &proj->local_ext, lower, upper);

  proj->boundary_values = mkarr(app->use_gpu, proj->basis.num_basis, proj->local_ext.volume);
  proj->boundary_buffer = mkarr(app->use_gpu, proj->basis.num_basis, proj->sol_skin.volume);

  proj->deflate_sol = gkyl_deflate_zsurf_new(&app->basis, &proj->basis, 0, app->use_gpu);

  proj->boundary_basis = app->use_gpu ? gkyl_cart_modal_serendip_cu_dev_new(2, 1) : &proj->basis;

  proj->copy_lcfs = gkyl_bc_basic_gyrokinetic_new(
    0, side ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE, GKYL_BC_GK_FIELD_BOUNDARY_VALUE,
    proj->boundary_basis, &proj->sol_skin, &proj->core_skin, proj->basis.num_basis, 2, app->use_gpu
  );

  proj->prepare_boundary = gk_field_ts_fem_perpproj_boundary_sol;
}

// Select y-directed input lines on the lower z boundary. Boundary lines become
// ordinary Dirichlet BCs; only strictly internal lines need bias-plane kernels.
static struct gkyl_fem_proj_bias_plane_list
gk_field_ts_fem_perpproj_init_bias(
  const gkyl_gyrokinetic_app *app, const struct gk_field *field,
  struct gk_field_ts_fem_perpproj *proj, struct gkyl_fem_proj_bc *bcs
)
{
  int lcfs_side = -1;
  if (app->gk_geom->has_LCFS) {
    lcfs_side = proj->local.lower[0] == app->global.lower[0] ? 1 : 0;
    if (lcfs_side) {
      bcs->up_type[0] = GKYL_FEM_PROJ_DIRICHLET_VARYING;
    } else {
      bcs->lo_type[0] = GKYL_FEM_PROJ_DIRICHLET_VARYING;
    }
  }
  struct gkyl_fem_proj_bias_plane_list selected = {0};
  const struct gkyl_poisson_bias_line_list *list = field->info.bias_line_list;
  if (list && list->num_bias_line) {
    struct gkyl_fem_proj_bias_plane *planes = gkyl_malloc(list->num_bias_line * sizeof(*planes));
    selected.bp = planes;
    for (int b = 0; b < list->num_bias_line; b++) {
      const struct gkyl_poisson_bias_line *line = &list->bl[b];
      if (line->perp_dirs[0] != 0 || line->perp_dirs[1] != 2) {
        gkyl_exit("ts_field_dg=false supports bias lines parallel to y only.");
      }
      if (fabs((line->perp_coords[1] - app->grid.lower[2]) / app->grid.dx[2]) > 1e-9) {
        continue;
      }
      double position = (line->perp_coords[0] - proj->grid.lower[0]) / proj->grid.dx[0];
      if (position < proj->local.lower[0] - 1 - 1e-9 || position > proj->local.upper[0] + 1e-9) {
        continue;
      }
      if (!isfinite(position) || !isfinite(line->val) || fabs(position - round(position)) > 1e-9) {
        gkyl_exit("ts_field_dg=false requires finite, grid-aligned bias lines.");
      }
      int face_idx = (int)round(position);
      int side = face_idx == proj->local.lower[0] - 1 ? 0 :
                 face_idx == proj->local.upper[0]     ? 1 :
                                                        -1;
      if (side >= 0) {
        enum gkyl_fem_proj_bc_type *type = side ? &bcs->up_type[0] : &bcs->lo_type[0];
        double *value = side ? &bcs->up_value[0] : &bcs->lo_value[0];
        if (*type == GKYL_FEM_PROJ_PERIODIC) {
          gkyl_exit("ts_field_dg=false: a boundary bias conflicts with periodic x.");
        }
        if (*type == GKYL_FEM_PROJ_DIRICHLET && *value != line->val) {
          gkyl_exit("ts_field_dg=false: conflicting radial boundary biases.");
        }
        *type = GKYL_FEM_PROJ_DIRICHLET;
        *value = line->val;
      } else {
        planes[selected.num_bias_plane++] = (struct gkyl_fem_proj_bias_plane){
          .dir = 0,
          .loc = line->perp_coords[0],
          .val = line->val,
        };
      }
    }
  }
  if (lcfs_side >= 0 &&
      (lcfs_side ? bcs->up_type[0] : bcs->lo_type[0]) == GKYL_FEM_PROJ_DIRICHLET_VARYING) {
    gk_field_ts_fem_perpproj_init_sol(app, proj, lcfs_side);
  }
  return selected;
}

static void
gk_field_ts_fem_perpproj_enabled(
  const gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *potential
)
{
  struct gk_field_ts_fem_perpproj *proj = field->ts_fem_perpproj;
  // The physical lower boundary is the UPPER face of the lower ghost cell.
  gkyl_deflate_zsurf_advance(
    proj->deflate, app->global.lower[2] - 1, &app->global_ext, &proj->local, potential,
    proj->buff_dg, 1
  );
  if (proj->boundary_values) {
    gkyl_array_copy(proj->boundary_values, proj->buff_dg);
  }
  proj->prepare_boundary(app, proj, potential);
  gkyl_fem_proj_set_rhs(
    proj->projector, proj->buff_dg, proj->boundary_values ? proj->boundary_values : proj->buff_dg
  );
  gkyl_fem_proj_solve(proj->projector, proj->buff_fem);
  // Constant-in-z inflation supplies the corrected boundary value without touching the RHS.
  gkyl_translate_dim_advance(
    proj->inflate, &proj->local, &proj->ghost, proj->buff_fem, 1, potential
  );
}

void
gk_field_ts_fem_perpproj_new(const gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  field->ts_fem_perpproj = NULL;
  field->ts_fem_perpproj_func = gk_field_ts_fem_perpproj_disabled;
  if (!gk_field_ts_fem_perpproj_requested(app, field)) {
    return;
  }

  // Abort if it's for a simulation not yet supported.
  gk_field_ts_fem_perpproj_validate(app, field);

  struct gk_field_ts_fem_perpproj *proj = gkyl_malloc(sizeof(*proj));
  *proj = (struct gk_field_ts_fem_perpproj){0};
  proj->prepare_boundary = gk_field_ts_fem_perpproj_boundary_disabled;

  // Perpendicular grid, basis and range.
  gkyl_rect_grid_init(&proj->grid, 2, app->grid.lower, app->grid.upper, app->grid.cells);
  gkyl_cart_modal_serendip(&proj->basis, 2, 1);
  gkyl_create_grid_ranges(&proj->grid, (int[]){1, 1}, &proj->local_ext, &proj->local);

  proj->ghost = app->global_lower_ghost[2];
  if (app->gk_geom->has_LCFS) {
    gkyl_sub_range_init(
      &proj->local, &proj->local_ext, app->global_core.lower, app->global_core.upper
    );
    proj->ghost = app->global_lower_ghost_par_core;
  }

  // DG and FEM buffers.
  proj->buff_dg = mkarr(app->use_gpu, proj->basis.num_basis, proj->local_ext.volume);
  proj->buff_fem = mkarr(app->use_gpu, proj->basis.num_basis, proj->local_ext.volume);

  // Consistent mass projection; NULL weights give the unweighted dx dy inner product.
  bool periodic_x = field->poisson_bcs.lo_type[0] == GKYL_POISSON_PERIODIC;
  bool dirichlet_lo = field->poisson_bcs.lo_type[0] == GKYL_POISSON_DIRICHLET;
  bool dirichlet_up = field->poisson_bcs.up_type[0] == GKYL_POISSON_DIRICHLET;
  struct gkyl_fem_proj_inp inp = {
    .grid = &proj->grid,
    .solve_range = &proj->local,
    .basis = &proj->basis,
    .bcs =
      {
        .lo_type =
          {periodic_x   ? GKYL_FEM_PROJ_PERIODIC :
           dirichlet_lo ? GKYL_FEM_PROJ_DIRICHLET :
                          GKYL_FEM_PROJ_DIRICHLET_VARYING,
           GKYL_FEM_PROJ_PERIODIC},
        .up_type =
          {periodic_x   ? GKYL_FEM_PROJ_PERIODIC :
           dirichlet_up ? GKYL_FEM_PROJ_DIRICHLET :
                          GKYL_FEM_PROJ_DIRICHLET_VARYING,
           GKYL_FEM_PROJ_PERIODIC},
        .lo_value = {field->poisson_bcs.lo_value[0].v[0]},
        .up_value = {field->poisson_bcs.up_value[0].v[0]},
      },
    .use_gpu = app->use_gpu,
  };
  struct gkyl_fem_proj_bias_plane_list biases =
    gk_field_ts_fem_perpproj_init_bias(app, field, proj, &inp.bcs);
  inp.bias_plane_list = &biases;
  proj->projector = gkyl_fem_proj_inew(&inp);
  gkyl_free((void *)biases.bp);

  // Inflation and deflation operators.
  proj->deflate = gkyl_deflate_zsurf_new(&app->basis, &proj->basis, 1, app->use_gpu);
  proj->inflate =
    gkyl_translate_dim_new(2, proj->basis, 3, app->basis, 2, GKYL_UPPER_EDGE, app->use_gpu);

  // Set function pointers.
  field->ts_fem_perpproj = proj;
  field->ts_fem_perpproj_func = gk_field_ts_fem_perpproj_enabled;
}

void
gk_field_ts_fem_perpproj_advance(
  const gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *potential
)
{
  field->ts_fem_perpproj_func(app, field, potential);
}

void
gk_field_ts_fem_perpproj_release(struct gk_field *field)
{
  struct gk_field_ts_fem_perpproj *proj = field->ts_fem_perpproj;
  if (!proj) {
    return;
  }
  gkyl_fem_proj_release(proj->projector);
  gkyl_deflate_zsurf_release(proj->deflate);
  gkyl_translate_dim_release(proj->inflate);
  gkyl_array_release(proj->buff_dg);
  gkyl_array_release(proj->buff_fem);
  if (proj->boundary_values) {
    gkyl_bc_basic_gyrokinetic_release(proj->copy_lcfs);
    gkyl_deflate_zsurf_release(proj->deflate_sol);
    gkyl_array_release(proj->boundary_values);
    gkyl_array_release(proj->boundary_buffer);
    if (proj->boundary_basis != &proj->basis) {
      gkyl_cart_modal_basis_release_cu(proj->boundary_basis);
    }
  }
  gkyl_free(proj);
}
