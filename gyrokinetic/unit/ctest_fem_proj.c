#include <acutest.h>
#include <gkyl_fem_proj.h>
#include <gkyl_fem_proj_priv.h>
#include <gkyl_fem_poisson.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_array_ops.h>
#include <math.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

static struct gkyl_array *
mkarr(bool gpu, int nb, long size)
{
  return gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nb, size) : gkyl_array_new(GKYL_DOUBLE, nb, size);
}

// Shared storage for one projection test; reference indices use physical node matching.
struct fproj_case {
  struct gkyl_rect_grid grid;
  struct gkyl_range ext, local;
  struct gkyl_basis basis;
  struct gkyl_fem_proj_bc bcs;
  bool grid_indexed; // Physical coordinates use the full grid and 1-based cell indices.
  struct gkyl_fem_proj_bias_plane_list biases; // Optional constraints for subrange tests.
  bool use_gpu; // Execution backend for the projector arrays.
  int direction; // Direction of the distinct-value affine boundary field.
  bool has_left, has_right, has_dirichlet;
  int num_ref_nodes, ids[9][8]; // Independent global numbering for the test checks.
  double nodes[16], shape_coeff[8][8];
  double prescribed[64]; // Expected values at constrained reference nodes.
  bool constrained[64];
  struct gkyl_array *input_ho, *weight_left_ho, *weight_right_ho, *boundary_ho;
  struct gkyl_array *result_ho, *previous_ho; // Host copies used by the checks.
  struct gkyl_array *rhs, *out, *weight_left, *weight_right, *boundary;
  struct gkyl_fem_proj *up;
};

enum fproj_boundary_pair {
  FPROJ_CONSTANT_VARYING,
  FPROJ_VARYING_CONSTANT,
  FPROJ_CONSTANT_CONSTANT
};

// A continuous nodal field compatible with every constant Dirichlet face.
static double
fproj_boundary_function(const struct fproj_case *ctx, const double *point)
{
  double variation = 0.2;
  for (int d = 0; d < ctx->grid.ndim; d++) {
    if (ctx->bcs.lo_type[d] == GKYL_FEM_PROJ_DIRICHLET) {
      variation *= point[d] - ctx->grid.lower[d] -
                   (ctx->grid_indexed ? ctx->local.lower[d] - 1 : 0) * ctx->grid.dx[d];
    }
    if (ctx->bcs.up_type[d] == GKYL_FEM_PROJ_DIRICHLET) {
      variation *= point[d] - (ctx->grid_indexed ?
                                 ctx->grid.lower[d] + ctx->local.upper[d] * ctx->grid.dx[d] :
                                 ctx->grid.upper[d]);
    }
    variation *= 0.7 + 0.1 * point[d];
  }
  return 1.25 + variation;
}

// An affine field with different end values and no transverse variation.
static double
fproj_distinct_boundary_function(const struct fproj_case *ctx, const double *point)
{
  return 0.5 + 0.3 * point[ctx->direction];
}

// Enumerate the nine ordered nonperiodic pairs and the paired-periodic case.
static struct gkyl_fem_proj_bc
fproj_boundary_conditions(int ndim, int code)
{
  struct gkyl_fem_proj_bc bcs = {0};
  const enum gkyl_fem_proj_bc_type types[] = {
    GKYL_FEM_PROJ_NONE, GKYL_FEM_PROJ_DIRICHLET, GKYL_FEM_PROJ_DIRICHLET_VARYING
  };

  for (int d = 0; d < ndim; d++) {
    int pair = code % 10;
    code /= 10;
    bcs.lo_type[d] = pair == 9 ? GKYL_FEM_PROJ_PERIODIC : types[pair % 3];
    bcs.up_type[d] = pair == 9 ? GKYL_FEM_PROJ_PERIODIC : types[pair / 3];
    bcs.lo_value[d] = bcs.up_value[d] = 1.25;
  }

  return bcs;
}

// Allocate one case's grid, basis, host reference arrays, and CPU/GPU solver arrays.
// Each test family supplies its own boundary conditions before creating the projector.
static void
fproj_case_init(
  struct fproj_case *ctx, int ndim, int poly_order, int weights, int single_mask, bool use_gpu
)
{
  *ctx = (struct fproj_case){
    .use_gpu = use_gpu,
    .has_left = weights && weights != 4,
    .has_right = weights && weights != 5,
  };
  int cells[] = {(single_mask & 1) ? 1 : 3, (single_mask & 2) ? 1 : 3};
  gkyl_rect_grid_init(&ctx->grid, ndim, (double[]){-0.2, -0.5}, (double[]){1.3, 0.5}, cells);

  // Non-unit lower indices and an embedded range catch incorrect flat indexing.
  int lower[] = {-1, 5}, upper[] = {cells[0] - 2, cells[1] + 4};
  int extlo[] = {-2, 4}, extup[] = {cells[0] - 1, cells[1] + 5};
  gkyl_range_init(&ctx->ext, ndim, extlo, extup);
  gkyl_sub_range_init(&ctx->local, &ctx->ext, lower, upper);

  gkyl_cart_modal_serendip(&ctx->basis, ndim, poly_order);

  int nb = ctx->basis.num_basis;
  ctx->basis.node_list(ctx->nodes);
  for (int n = 0; n < nb; n++) {
    double values[8] = {0};
    values[n] = 1;
    ctx->basis.nodal_to_modal(values, ctx->shape_coeff[n]);
  }

  ctx->input_ho = mkarr(false, nb, ctx->ext.volume);
  ctx->weight_left_ho = mkarr(false, nb, ctx->ext.volume);
  ctx->weight_right_ho = mkarr(false, nb, ctx->ext.volume);
  ctx->boundary_ho = mkarr(false, nb, ctx->ext.volume);
  ctx->result_ho = mkarr(false, nb, ctx->ext.volume);
  ctx->previous_ho = mkarr(false, nb, ctx->ext.volume);
  ctx->rhs = mkarr(use_gpu, nb, ctx->ext.volume);
  ctx->out = mkarr(use_gpu, nb, ctx->ext.volume);
  ctx->weight_left = mkarr(use_gpu, nb, ctx->ext.volume);
  ctx->weight_right = mkarr(use_gpu, nb, ctx->ext.volume);
  ctx->boundary = mkarr(use_gpu, nb, ctx->ext.volume);

  struct gkyl_range_iter iter;
  int cell = 0;
  gkyl_range_iter_init(&iter, &ctx->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&ctx->local, iter.idx);
    double *input_coeff = gkyl_array_fetch(ctx->input_ho, loc),
           *left_coeff = gkyl_array_fetch(ctx->weight_left_ho, loc),
           *right_coeff = gkyl_array_fetch(ctx->weight_right_ho, loc);

    left_coeff[0] = (weights == 1 ? 1.0 : 2.0) * pow(sqrt(2.0), ndim);
    if (weights >= 2) {
      for (int k = 1; k < nb; k++) {
        left_coeff[k] = 0.015 * cos(k + cell);
      }
    }
    for (int k = 0; k < nb; k++) {
      input_coeff[k] = sin(0.7 * (cell + 1) * (k + 1));
      right_coeff[k] = weights == 3 ? 0.6 * left_coeff[k] : left_coeff[k];
    }

    cell++;
  }
}

// Match physical nodes independently of the production indexing kernels, evaluate
// this family's boundary field, and create the projector with its fixed weights.
static void
fproj_case_init_projector(
  struct fproj_case *ctx, double (*boundary_function)(const struct fproj_case *, const double *)
)
{
  int ndim = ctx->grid.ndim, nb = ctx->basis.num_basis;
  double coords[64][2] = {{0}};
  struct gkyl_range_iter iter;
  int cell = 0;

  gkyl_range_iter_init(&iter, &ctx->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&ctx->local, iter.idx);
    double values[8];
    for (int n = 0; n < nb; n++) {
      double xy[2] = {0, 0};
      for (int d = 0; d < ndim; d++) {
        xy[d] = ctx->grid.lower[d] + (iter.idx[d] - (ctx->grid_indexed ? 1 : ctx->local.lower[d]) +
                                      0.5 * (ctx->nodes[n * ndim + d] + 1)) *
                                       ctx->grid.dx[d];
        double lower =
          ctx->grid.lower[d] + (ctx->grid_indexed ? ctx->local.lower[d] - 1 : 0) * ctx->grid.dx[d];
        double upper = lower + (ctx->local.upper[d] - ctx->local.lower[d] + 1) * ctx->grid.dx[d];
        if (ctx->bcs.lo_type[d] == GKYL_FEM_PROJ_PERIODIC && fabs(xy[d] - upper) < 1e-13) {
          xy[d] = lower;
        }
      }

      int global = 0;
      while (global < ctx->num_ref_nodes &&
             (fabs(coords[global][0] - xy[0]) > 1e-12 || fabs(coords[global][1] - xy[1]) > 1e-12)) {
        global++;
      }

      if (global == ctx->num_ref_nodes) {
        TEST_ASSERT(ctx->num_ref_nodes < 64);
        coords[ctx->num_ref_nodes][0] = xy[0];
        coords[ctx->num_ref_nodes++][1] = xy[1];
      }

      ctx->ids[cell][n] = global;
      values[n] = boundary_function(ctx, xy);
      for (int d = 0; d < ndim; d++) {
        for (int side = 0; side < 2; side++) {
          int type = side ? ctx->bcs.up_type[d] : ctx->bcs.lo_type[d];
          bool at = iter.idx[d] == (side ? ctx->local.upper[d] : ctx->local.lower[d]) &&
                    ctx->nodes[n * ndim + d] == (side ? 1 : -1);
          if (at && (type == GKYL_FEM_PROJ_DIRICHLET || type == GKYL_FEM_PROJ_DIRICHLET_VARYING)) {
            ctx->has_dirichlet = true;
            ctx->constrained[global] = true;
            ctx->prescribed[global] = values[n];
          }
        }
      }
      for (int b = 0; b < ctx->biases.num_bias_plane; b++) {
        const struct gkyl_fem_proj_bias_plane *plane = &ctx->biases.bp[b];
        int dir = plane->dir;
        double location = plane->loc;
        double lower = ctx->grid.lower[dir] + (ctx->local.lower[dir] - 1) * ctx->grid.dx[dir];
        double upper = ctx->grid.lower[dir] + ctx->local.upper[dir] * ctx->grid.dx[dir];
        if (ctx->bcs.lo_type[dir] == GKYL_FEM_PROJ_PERIODIC && fabs(location - upper) < 1e-12) {
          location = lower;
        }
        if (fabs(xy[dir] - location) < 1e-12) {
          ctx->has_dirichlet = true;
          ctx->constrained[global] = true;
          ctx->prescribed[global] = plane->val;
        }
      }
    }

    ctx->basis.nodal_to_modal(values, gkyl_array_fetch(ctx->boundary_ho, loc));
    cell++;
  }

  gkyl_array_copy(ctx->weight_left, ctx->weight_left_ho);
  gkyl_array_copy(ctx->weight_right, ctx->weight_right_ho);
  gkyl_array_copy(ctx->boundary, ctx->boundary_ho);

  struct gkyl_fem_proj_inp inp = {
    .grid = &ctx->grid,
    .basis = &ctx->basis,
    .solve_range = &ctx->local,
    .use_gpu = ctx->use_gpu,
    .bcs = ctx->bcs,
    .bias_plane_list = &ctx->biases,
    .weight_left = ctx->has_left ? ctx->weight_left : NULL,
    .weight_right = ctx->has_right ? ctx->weight_right : NULL,
  };
  ctx->up = gkyl_fem_proj_inew(&inp);
}

// Check weak residuals, boundary values, continuity, conservation, and ghost cells.
// Also check field preservation when the input is continuous and the weights match.
static void
fproj_check_solution(struct fproj_case *ctx, bool check_preservation)
{
  int ndim = ctx->grid.ndim, nb = ctx->basis.num_basis;
  // Four-point quadrature integrates all weighted p2 products exactly.
  const double *quad_x = gkyl_gauss_ordinates[4], *quad_w = gkyl_gauss_weights[4];
  gkyl_array_copy(ctx->result_ho, ctx->out);
  int cell = 0;
  struct gkyl_range_iter iter;
  double seen[64], residual[64] = {0}, residual_scale[64] = {0};
  bool visited[64] = {0};
  double source_integral = 0, output_integral = 0;

  gkyl_range_iter_init(&iter, &ctx->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&ctx->local, iter.idx);

    const double *got = gkyl_array_cfetch(ctx->result_ho, loc);
    for (int n = 0; n < nb; n++) {
      double at_node = ctx->basis.eval_expand(&ctx->nodes[n * ndim], got);
      int global = ctx->ids[cell][n];
      if (ctx->constrained[global]) {
        TEST_CHECK(
          fabs(at_node - ctx->prescribed[global]) < 1e-11 * (1 + fabs(ctx->prescribed[global]))
        );
      }
      if (visited[global]) {
        TEST_CHECK(fabs(at_node - seen[global]) < 1e-11);
      }
      seen[global] = at_node;
      visited[global] = true;
    }

    if (check_preservation) {
      const double *input_coeff = gkyl_array_cfetch(ctx->input_ho, loc);
      for (int n = 0; n < nb; n++) {
        TEST_CHECK(fabs(got[n] - input_coeff[n]) < 1e-11);
      }
    }

    for (int qx = 0; qx < 4; qx++) {
      for (int qy = 0; qy < (ndim == 2 ? 4 : 1); qy++) {
        double point[] = {quad_x[qx], quad_x[qy]};

        double left_weight =
          ctx->has_left ?
            ctx->basis.eval_expand(point, gkyl_array_cfetch(ctx->weight_left_ho, loc)) :
            1;
        double right_weight =
          ctx->has_right ?
            ctx->basis.eval_expand(point, gkyl_array_cfetch(ctx->weight_right_ho, loc)) :
            1;

        double input_value = ctx->basis.eval_expand(point, gkyl_array_cfetch(ctx->input_ho, loc));
        double projected = ctx->basis.eval_expand(point, got);

        double measure = quad_w[qx] * (ndim == 2 ? quad_w[qy] : 1);
        source_integral += measure * right_weight * input_value;
        output_integral += measure * left_weight * projected;

        // Assemble the weak residual using coordinate-based node matching.
        // The uniform cell Jacobian is common to every contribution and cancels.
        for (int n = 0; n < nb; n++) {
          double shape = ctx->basis.eval_expand(point, ctx->shape_coeff[n]);
          int global = ctx->ids[cell][n];
          residual[global] +=
            measure * shape * (left_weight * projected - right_weight * input_value);
          residual_scale[global] +=
            measure * fabs(shape) *
            (fabs(left_weight * projected) + fabs(right_weight * input_value));
        }
      }
    }
    cell++;
  }

  // Dirichlet nodes satisfy prescribed values; all other nodes satisfy the projection equation.
  for (int n = 0; n < ctx->num_ref_nodes; n++) {
    if (!ctx->constrained[n]) {
      TEST_CHECK(fabs(residual[n]) < 1e-11 * (1 + residual_scale[n]));
      TEST_MSG("node=%d residual=%.17g scale=%.17g", n, residual[n], residual_scale[n]);
    }
  }
  if (!ctx->has_dirichlet) {
    TEST_CHECK(fabs(output_integral - source_integral) < 1e-10);
  }

  gkyl_range_iter_init(&iter, &ctx->ext);
  while (gkyl_range_iter_next(&iter)) {
    if (!gkyl_range_contains_idx(&ctx->local, iter.idx)) {
      const double *v = gkyl_array_cfetch(ctx->result_ho, gkyl_range_idx(&ctx->ext, iter.idx));
      for (int n = 0; n < nb; n++) {
        TEST_CHECK(v[n] == 12345.0);
      }
    }
  }
}

// Project the previous solution in place and check that it is unchanged.
// Call only when the left and right weights match.
static void
fproj_check_idempotence(struct fproj_case *ctx)
{
  int nb = ctx->basis.num_basis;
  gkyl_array_copy(ctx->previous_ho, ctx->result_ho);

  // Here we check that in-place projection is supported. Also checks that
  // this is indeed a projection operator (idempotence): that projecting a
  // second time leaves the answer unchanged.
  gkyl_fem_proj_set_rhs(ctx->up, ctx->out, ctx->boundary);
  gkyl_fem_proj_solve(ctx->up, ctx->out);

  gkyl_array_copy(ctx->result_ho, ctx->out);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ctx->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&ctx->local, iter.idx);
    const double *a = gkyl_array_cfetch(ctx->previous_ho, loc),
                 *b = gkyl_array_cfetch(ctx->result_ho, loc);
    for (int n = 0; n < nb; n++) {
      TEST_CHECK(fabs(a[n] - b[n]) < 1e-11);
    }
  }
}

// Release the projector before releasing its weight arrays and the test buffers.
static void
fproj_case_release(struct fproj_case *ctx)
{
  gkyl_fem_proj_release(ctx->up);
  gkyl_array_release(ctx->input_ho);
  gkyl_array_release(ctx->weight_left_ho);
  gkyl_array_release(ctx->weight_right_ho);
  gkyl_array_release(ctx->boundary_ho);
  gkyl_array_release(ctx->result_ho);
  gkyl_array_release(ctx->previous_ho);
  gkyl_array_release(ctx->rhs);
  gkyl_array_release(ctx->out);
  gkyl_array_release(ctx->weight_left);
  gkyl_array_release(ctx->weight_right);
  gkyl_array_release(ctx->boundary);
}

// General solver case: irregular DG, compatible continuous, and constant inputs.
static void
fproj_solver_case(int ndim, int poly_order, int scenario, int weights, int single_mask, bool use_gpu)
{
  TEST_CASE_(
    "dim=%d p=%d bc=%d weights=%d single_mask=%d use_gpu=%d", ndim, poly_order, scenario, weights,
    single_mask, use_gpu
  );
  struct fproj_case ctx;
  fproj_case_init(&ctx, ndim, poly_order, weights, single_mask, use_gpu);

  ctx.bcs = fproj_boundary_conditions(ndim, scenario);
  fproj_case_init_projector(&ctx, fproj_boundary_function);

  for (int trial = 0; trial < 3; trial++) {
    if (trial == 1) {
      gkyl_array_copy(ctx.input_ho, ctx.boundary_ho);
    }
    if (trial == 2) {
      gkyl_array_clear(ctx.input_ho, 0);
      gkyl_array_shiftc(ctx.input_ho, 1.25 * pow(sqrt(2.0), ndim), 0);
      gkyl_array_clear(ctx.boundary_ho, 0);
      gkyl_array_shiftc(ctx.boundary_ho, 1.25 * pow(sqrt(2.0), ndim), 0);
      gkyl_array_copy(ctx.boundary, ctx.boundary_ho);
      for (int n = 0; n < ctx.num_ref_nodes; n++) {
        if (ctx.constrained[n]) {
          ctx.prescribed[n] = 1.25;
        }
      }
    }
    gkyl_array_copy(ctx.rhs, ctx.input_ho);
    gkyl_array_clear(ctx.out, 12345.0);

    // Perform projection.
    gkyl_fem_proj_set_rhs(ctx.up, ctx.rhs, ctx.boundary);
    gkyl_fem_proj_solve(ctx.up, ctx.out);

    // Check result.
    fproj_check_solution(&ctx, trial > 0 && weights <= 2);
    if (weights <= 2) {
      fproj_check_idempotence(&ctx);
    }
  }

  fproj_case_release(&ctx);
}

// Distinct-value case: irregular DG and compatible affine inputs. A constant
// input cannot satisfy the different prescribed values on the two end faces.
static void
fproj_distinct_boundary_case(
  int ndim, int poly_order, int direction, enum fproj_boundary_pair boundary_pair, int weights,
  int single_mask, bool use_gpu
)
{
  TEST_CASE_(
    "dim=%d p=%d direction=%d boundary_pair=%d weights=%d single_mask=%d use_gpu=%d", ndim,
    poly_order, direction, boundary_pair, weights, single_mask, use_gpu
  );
  struct fproj_case ctx;
  fproj_case_init(&ctx, ndim, poly_order, weights, single_mask, use_gpu);

  ctx.direction = direction;
  for (int d = 0; d < ndim; d++) {
    ctx.bcs.lo_type[d] = ctx.bcs.up_type[d] = GKYL_FEM_PROJ_PERIODIC;
  }
  ctx.bcs.lo_type[direction] = boundary_pair == FPROJ_VARYING_CONSTANT ?
                                 GKYL_FEM_PROJ_DIRICHLET_VARYING :
                                 GKYL_FEM_PROJ_DIRICHLET;
  ctx.bcs.up_type[direction] = boundary_pair == FPROJ_CONSTANT_VARYING ?
                                 GKYL_FEM_PROJ_DIRICHLET_VARYING :
                                 GKYL_FEM_PROJ_DIRICHLET;
  ctx.bcs.lo_value[direction] = 0.5 + 0.3 * ctx.grid.lower[direction];
  ctx.bcs.up_value[direction] = 0.5 + 0.3 * ctx.grid.upper[direction];

  fproj_case_init_projector(&ctx, fproj_distinct_boundary_function);

  for (int trial = 0; trial < 2; trial++) {
    if (trial == 1) {
      gkyl_array_copy(ctx.input_ho, ctx.boundary_ho);
    }
    gkyl_array_copy(ctx.rhs, ctx.input_ho);
    gkyl_array_clear(ctx.out, 12345.0);

    // Perform projection.
    gkyl_fem_proj_set_rhs(ctx.up, ctx.rhs, ctx.boundary);
    gkyl_fem_proj_solve(ctx.up, ctx.out);

    // Check result.
    fproj_check_solution(&ctx, trial > 0 && weights <= 2);
    if (weights <= 2) {
      fproj_check_idempotence(&ctx);
    }
  }

  fproj_case_release(&ctx);
}

// Test projection for a fixed dimension, polynomial order, and CPU/GPU backend.
// In each direction, test all nine lower/upper pairs of NONE, constant Dirichlet,
// and varying Dirichlet conditions, plus paired periodic boundaries. Combine
// these with all six weight settings (implicit/explicit unity, equal/unequal
// varying weights, and either side weighted alone) and one or three cells per
// direction. Each configuration checks an irregular DG input, a compatible
// continuous field, and a constant field using fproj_solver_case.
static void
fproj_solver_combos(int ndim, int poly_order, bool use_gpu)
{
  for (int bc = 0; bc < (ndim == 1 ? 10 : 100); bc++) {
    for (int weights = 0; weights < 6; weights++) {
      for (int single_mask = 0; single_mask < (1 << ndim); single_mask++) {
        fproj_solver_case(ndim, poly_order, bc, weights, single_mask, use_gpu);
      }
    }
  }
}

// Check that lower and upper Dirichlet values are applied to the correct faces
// for a fixed dimension, polynomial order, and CPU/GPU backend. In turn, impose
// different values in each direction using constant/varying, varying/constant,
// and constant/constant boundary pairs; the other direction, if present, is
// periodic. Test one or three cells per direction, including cells touching both
// boundaries, with implicit unity weights and unequal varying weights. Each
// configuration checks an irregular DG input and a compatible continuous field;
// a constant field cannot satisfy the distinct prescribed values.
static void
fproj_distinct_boundary_combos(int ndim, int poly_order, bool use_gpu)
{
  for (int single_mask = 0; single_mask < (1 << ndim); single_mask++) {
    for (int direction = 0; direction < ndim; direction++) {
      for (int pair = FPROJ_CONSTANT_VARYING; pair <= FPROJ_CONSTANT_CONSTANT; pair++) {
        fproj_distinct_boundary_case(
          ndim, poly_order, direction, (enum fproj_boundary_pair)pair, 0, single_mask, use_gpu
        );
        fproj_distinct_boundary_case(
          ndim, poly_order, direction, (enum fproj_boundary_pair)pair, 3, single_mask, use_gpu
        );
      }
    }
  }
}

static void
test_fem_proj_solver_1x_p1_ho(void)
{
  fproj_solver_combos(1, 1, false);
}

static void
test_fem_proj_solver_1x_p2_ho(void)
{
  fproj_solver_combos(1, 2, false);
}

static void
test_fem_proj_solver_2x_p1_ho(void)
{
  fproj_solver_combos(2, 1, false);
}

static void
test_fem_proj_solver_2x_p2_ho(void)
{
  fproj_solver_combos(2, 2, false);
}

static void
test_fem_proj_distinct_boundary_1x_p1_ho(void)
{
  fproj_distinct_boundary_combos(1, 1, false);
}

static void
test_fem_proj_distinct_boundary_1x_p2_ho(void)
{
  fproj_distinct_boundary_combos(1, 2, false);
}

static void
test_fem_proj_distinct_boundary_2x_p1_ho(void)
{
  fproj_distinct_boundary_combos(2, 1, false);
}

static void
test_fem_proj_distinct_boundary_2x_p2_ho(void)
{
  fproj_distinct_boundary_combos(2, 2, false);
}

#ifdef GKYL_HAVE_CUDA
static void
test_fem_proj_solver_1x_p1_dev(void)
{
  fproj_solver_combos(1, 1, true);
}

static void
test_fem_proj_solver_1x_p2_dev(void)
{
  fproj_solver_combos(1, 2, true);
}

static void
test_fem_proj_solver_2x_p1_dev(void)
{
  fproj_solver_combos(2, 1, true);
}

static void
test_fem_proj_solver_2x_p2_dev(void)
{
  fproj_solver_combos(2, 2, true);
}

static void
test_fem_proj_distinct_boundary_1x_p1_dev(void)
{
  fproj_distinct_boundary_combos(1, 1, true);
}

static void
test_fem_proj_distinct_boundary_1x_p2_dev(void)
{
  fproj_distinct_boundary_combos(1, 2, true);
}

static void
test_fem_proj_distinct_boundary_2x_p1_dev(void)
{
  fproj_distinct_boundary_combos(2, 1, true);
}

static void
test_fem_proj_distinct_boundary_2x_p2_dev(void)
{
  fproj_distinct_boundary_combos(2, 2, true);
}

#endif

// Verify compatibility with the previous zero-epsilon Poisson projection.
static void
fproj_poisson_comparison(bool gpu)
{
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, (double[]){0, -1}, (double[]){1, 1}, (int[]){3, 4});

  struct gkyl_range ext, local;
  gkyl_create_grid_ranges(&grid, (int[]){1, 1}, &ext, &local);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, 1);

  struct gkyl_array *host = mkarr(false, 4, ext.volume), *a = mkarr(false, 4, ext.volume),
                    *b = mkarr(false, 4, ext.volume);

  for (long i = 0; i < ext.volume; i++) {
    double *v = gkyl_array_fetch(host, i);
    for (int k = 0; k < 4; k++) {
      v[k] = sin(0.3 * (i + 1) * (k + 1));
    }
  }

  struct gkyl_array *rhs = mkarr(gpu, 4, ext.volume), *out = mkarr(gpu, 4, ext.volume);
  struct gkyl_array *epsilon = mkarr(gpu, 4, ext.volume), *ksq = mkarr(gpu, 4, ext.volume);

  gkyl_array_copy(rhs, host);
  gkyl_array_clear(epsilon, 0);
  gkyl_array_clear(ksq, 0);
  gkyl_array_shiftc(ksq, -2, 0);

  struct gkyl_poisson_bc oldbc = {
    .lo_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_PERIODIC},
    .up_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_PERIODIC},
  };
  oldbc.lo_value[0].v[0] = 1.25;
  oldbc.up_value[0].v[0] = -0.75;

  struct gkyl_fem_poisson *old =
    gkyl_fem_poisson_new(&local, &grid, basis, &oldbc, NULL, epsilon, ksq, true, gpu);

  // Solve Poisson projection.
  gkyl_fem_poisson_set_rhs(old, rhs, NULL);
  gkyl_fem_poisson_solve(old, out);
  gkyl_array_copy(a, out);

  // Project with fem_proj.
  struct gkyl_fem_proj_inp inp = {
    .grid = &grid,
    .basis = &basis,
    .solve_range = &local,
    .use_gpu = gpu,
    .bcs =
      {
        .lo_type = {GKYL_FEM_PROJ_DIRICHLET, GKYL_FEM_PROJ_PERIODIC},
        .up_type = {GKYL_FEM_PROJ_DIRICHLET, GKYL_FEM_PROJ_PERIODIC},
        .lo_value = {1.25},
        .up_value = {-0.75},
      },
  };

  struct gkyl_fem_proj *up = gkyl_fem_proj_inew(&inp);
  gkyl_fem_proj_set_rhs(up, rhs, NULL);
  gkyl_fem_proj_solve(up, out);
  gkyl_array_copy(b, out);

  // Check for equivalence.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *av = gkyl_array_cfetch(a, loc), *bv = gkyl_array_cfetch(b, loc);
    for (int k = 0; k < 4; k++) {
      TEST_CHECK(fabs(av[k] - bv[k]) < 1e-11);
    }
  }

  gkyl_fem_proj_release(up);
  gkyl_fem_poisson_release(old);
  gkyl_array_release(host);
  gkyl_array_release(a);
  gkyl_array_release(b);
  gkyl_array_release(rhs);
  gkyl_array_release(out);
  gkyl_array_release(epsilon);
  gkyl_array_release(ksq);
}
static void
test_fem_proj_poisson_ho(void)
{
  fproj_poisson_comparison(false);
}
#ifdef GKYL_HAVE_CUDA
static void
test_fem_proj_poisson_dev(void)
{
  fproj_poisson_comparison(true);
}
#endif

// Exercise constraints on a true portion of the full grid, with untouched cells
// on both sides. Reference residuals use physical node matching and quadrature.
static void
fproj_bias_case(
  int ndim, int poly_order, int direction, int weights, int single_mask, int scenario, bool use_gpu
)
{
  TEST_CASE_(
    "bias dim=%d p=%d dir=%d weights=%d single=%d scenario=%d gpu=%d", ndim, poly_order, direction,
    weights, single_mask, scenario, use_gpu
  );

  struct fproj_case ctx;
  fproj_case_init(&ctx, ndim, poly_order, weights, single_mask, use_gpu);

  ctx.grid_indexed = true;
  int cells[2], upper[2];
  double lower[2], high[2];
  for (int d = 0; d < ndim; d++) {
    cells[d] = ctx.grid.cells[d] + 2;
    upper[d] = cells[d] - 1;
    lower[d] = ctx.grid.lower[d] - ctx.grid.dx[d];
    high[d] = ctx.grid.upper[d] + ctx.grid.dx[d];
    if (d != direction && scenario != 4) {
      ctx.bcs.lo_type[d] = ctx.bcs.up_type[d] = GKYL_FEM_PROJ_PERIODIC;
    }
  }

  // Array sizes and cell spacing are unchanged from fproj_case_init.
  gkyl_rect_grid_init(&ctx.grid, ndim, lower, high, cells);
  gkyl_range_init(&ctx.ext, ndim, (int[]){1, 1}, cells);
  gkyl_sub_range_init(&ctx.local, &ctx.ext, (int[]){2, 2}, upper);

  double lo = lower[direction] + ctx.grid.dx[direction];
  double hi = high[direction] - ctx.grid.dx[direction];
  double inside = lo + ctx.grid.dx[direction];
  struct gkyl_fem_proj_bias_plane planes[4] = {{.dir = direction, .loc = inside, .val = 0.7}};
  ctx.biases = (struct gkyl_fem_proj_bias_plane_list){.num_bias_plane = 1, .bp = planes};
  if (scenario == 0) { // Unbiased subrange must retain the original projection equations.
    ctx.biases.num_bias_plane = 0;
  } else if (scenario == 2) { // Distinct, parallel biases, including an end plane.
    planes[1] = (struct gkyl_fem_proj_bias_plane){.dir = direction, .loc = lo, .val = -0.4};
    ctx.biases.num_bias_plane = 2;
  } else if (scenario == 3) { // Both physical copies of a periodic end plane.
    ctx.bcs.lo_type[direction] = ctx.bcs.up_type[direction] = GKYL_FEM_PROJ_PERIODIC;
    planes[0].loc = lo;
    planes[1] = planes[0];
    planes[1].loc = hi;
    ctx.biases.num_bias_plane = 2;
  } else if (scenario == 4 && ndim == 2) { // Compatible crossing planes.
    int other = 1 - direction;
    planes[1] = (struct gkyl_fem_proj_bias_plane){
      .dir = other,
      .loc = lower[other] + 2 * ctx.grid.dx[other],
      .val = 0.7,
    };
    ctx.biases.num_bias_plane = 2;
  } else if (scenario == 5) { // A bias agrees with a constant outer boundary.
    planes[0].loc = lo;
    planes[0].val = 1.25;
    ctx.bcs.lo_type[direction] = GKYL_FEM_PROJ_DIRICHLET;
    ctx.bcs.lo_value[direction] = 1.25;
  }
  // Outside planes must have no effect, even if their values differ.
  planes[ctx.biases.num_bias_plane++] =
    (struct gkyl_fem_proj_bias_plane){.dir = direction, .loc = lower[direction], .val = 5.0};

  fproj_case_init_projector(&ctx, fproj_boundary_function);

  gkyl_array_copy(ctx.rhs, ctx.input_ho);
  gkyl_array_clear(ctx.out, 12345.0);

  // Perform projection.
  gkyl_fem_proj_set_rhs(ctx.up, ctx.rhs, ctx.boundary);
  gkyl_fem_proj_solve(ctx.up, ctx.out);

  // Check result.
  fproj_check_solution(&ctx, false);
  if (weights <= 2) {
    fproj_check_idempotence(&ctx);
  }

  fproj_case_release(&ctx);
}

// Cover both directions, orders and backends, including periodic aliases in
// one-cell directions and independent left/right weights on a partial range.
static void
fproj_bias_combos(int ndim, bool use_gpu)
{
  for (int poly_order = 1; poly_order <= 2; poly_order++) {
    for (int direction = 0; direction < ndim; direction++) {
      for (int weights = 0; weights < 6; weights++) {
        for (int single = 0; single < (1 << ndim); single++) {
          for (int scenario = 0; scenario < 6; scenario++) {
            fproj_bias_case(ndim, poly_order, direction, weights, single, scenario, use_gpu);
          }
        }
      }
    }
  }
}

static void
test_fem_proj_bias_1x_ho(void)
{
  fproj_bias_combos(1, false);
}
static void
test_fem_proj_bias_2x_ho(void)
{
  fproj_bias_combos(2, false);
}
#ifdef GKYL_HAVE_CUDA
static void
test_fem_proj_bias_1x_dev(void)
{
  fproj_bias_combos(1, true);
}
static void
test_fem_proj_bias_2x_dev(void)
{
  fproj_bias_combos(2, true);
}
#endif

static void
test_fem_proj_invalid_input(void)
{
  for (int scenario = 0; scenario < 7; scenario++) {
    pid_t child = fork();
    TEST_ASSERT(child >= 0);
    if (child == 0) {
      freopen("/dev/null", "w", stderr);
      int ndim = scenario == 0 ? 3 : 2;
      struct gkyl_rect_grid grid;
      gkyl_rect_grid_init(&grid, ndim, (double[]){0, 0, 0}, (double[]){1, 1, 1}, (int[]){2, 2, 2});

      struct gkyl_range range;
      gkyl_range_init(&range, ndim, (int[]){1, 1, 1}, (int[]){2, 2, 2});

      struct gkyl_basis basis;
      gkyl_cart_modal_serendip(&basis, ndim, 1);

      struct gkyl_fem_proj_inp inp = {.grid = &grid, .basis = &basis, .solve_range = &range};
      if (scenario == 1) {
        inp.bcs.lo_type[0] = GKYL_FEM_PROJ_PERIODIC;
      }
      struct gkyl_fem_proj_bias_plane planes[] = {
        {.dir = 0, .loc = 0.5, .val = 0.7}, {.dir = 0, .loc = 0.5, .val = -0.4}
      };
      struct gkyl_fem_proj_bias_plane_list list = {.num_bias_plane = 1, .bp = planes};
      if (scenario >= 2) {
        inp.bias_plane_list = &list;
      }
      if (scenario == 2) {
        planes[0].loc = 0.3; // Not a grid face.
      }
      if (scenario == 3) {
        list.num_bias_plane = 2; // Conflicting duplicate.
      }
      if (scenario == 4) { // Conflicting intersection.
        list.num_bias_plane = 2;
        planes[1].dir = 1;
      }
      if (scenario == 5) { // Inconsistent outer Dirichlet data.
        inp.bcs.lo_type[1] = GKYL_FEM_PROJ_DIRICHLET;
      }
      if (scenario == 6) {
        planes[0].dir = 2; // Invalid plane direction.
      }
      struct gkyl_fem_proj *up = gkyl_fem_proj_inew(&inp);
      gkyl_fem_proj_release(up);
      _exit(0);
    }
    int status;
    waitpid(child, &status, 0);
    TEST_CHECK(WIFEXITED(status) && WEXITSTATUS(status) != 0);
  }
}

// Check installed kernels and dispatch against independently enumerated cell owners.
// The reference traverses cells with the last coordinate fastest, as in fem_poisson.
static void
test_fem_proj_node_indexing(void)
{
  for (int ndim = 1; ndim <= 2; ndim++) {
    for (int order = 1; order <= 2; order++) {
      struct gkyl_basis basis;
      gkyl_cart_modal_serendip(&basis, ndim, order);

      double nodes[16];
      basis.node_list(nodes);

      for (int mask = 0; mask < (1 << ndim); mask++) {
        struct gkyl_fem_proj_bc bcs = {0};
        for (int d = 0; d < ndim; d++) {
          if (mask & (1 << d)) {
            bcs.lo_type[d] = bcs.up_type[d] = GKYL_FEM_PROJ_PERIODIC;
          }
        }
        struct fem_proj_kernels kernels =
          fem_proj_choose_kernels(ndim, order, false, false, mask, &bcs);
        bool periodic[2] = {(mask & 1) != 0, (mask & 2) != 0};

        for (int nx = 1; nx <= 3; nx++) {
          for (int ny = 1; ny <= (ndim == 2 ? 3 : 1); ny++) {
            TEST_CASE_("dim=%d p=%d periodic=%d cells=%dx%d", ndim, order, mask, nx, ny);
            int cells[2] = {nx, ny}, coordinates[64][2] = {{0}}, count = 0;

            for (int ix = 0; ix < nx; ix++) {
              for (int iy = 0; iy < ny; iy++) {
                int cell[2] = {ix, iy};
                for (int n = 0; n < basis.num_basis; n++) {
                  bool owns = true;
                  for (int d = 0; d < ndim; d++) {
                    if (nodes[n * ndim + d] == 1 && (cell[d] != cells[d] - 1 || periodic[d])) {
                      owns = false;
                    }
                  }
                  if (owns) {
                    TEST_ASSERT(count < 64);
                    for (int d = 0; d < ndim; d++) {
                      coordinates[count][d] = 2 * cell[d] + (int)nodes[n * ndim + d] + 1;
                    }
                    count++;
                  }
                }
              }
            }

            TEST_CHECK(kernels.num_nodes(cells, periodic) == count);

            int lower[] = {-2, 5}, upper[] = {nx - 3, ny + 4};
            int extlo[] = {-3, 4}, extup[] = {nx - 2, ny + 5};

            struct gkyl_range ext, range;
            gkyl_range_init(&ext, ndim, extlo, extup);
            gkyl_sub_range_init(&range, &ext, lower, upper);

            bool visited[64] = {false};
            struct gkyl_range_iter iter;
            gkyl_range_iter_init(&iter, &range);
            while (gkyl_range_iter_next(&iter)) {
              long global_idx[FEM_PROJ_MAX_NODES];
              fem_proj_global_indices(&kernels, &range, cells, iter.idx, global_idx);

              for (int n = 0; n < basis.num_basis; n++) {
                int point[2] = {0, 0};
                for (int d = 0; d < ndim; d++) {
                  point[d] = 2 * (iter.idx[d] - lower[d]) + (int)nodes[n * ndim + d] + 1;
                  if (periodic[d]) {
                    point[d] %= 2 * cells[d];
                  }
                }
                int expected = 0;
                while (expected < count && (coordinates[expected][0] != point[0] ||
                                            coordinates[expected][1] != point[1])) {
                  expected++;
                }
                TEST_ASSERT(expected < count);
                TEST_CHECK(global_idx[n] == expected);
                TEST_ASSERT(global_idx[n] >= 0 && global_idx[n] < count);
                visited[global_idx[n]] = true;
              }
            }

            for (int n = 0; n < count; n++) {
              TEST_CHECK(visited[n]);
            }
          }
        }
      }
    }
  }
}

TEST_LIST = {
  {"fem_proj_node_indexing", test_fem_proj_node_indexing},
  {"fem_proj_solver_1x_p1_ho", test_fem_proj_solver_1x_p1_ho},
  {"fem_proj_solver_1x_p2_ho", test_fem_proj_solver_1x_p2_ho},
  {"fem_proj_solver_2x_p1_ho", test_fem_proj_solver_2x_p1_ho},
  {"fem_proj_solver_2x_p2_ho", test_fem_proj_solver_2x_p2_ho},
  {"fem_proj_distinct_boundary_1x_p1_ho", test_fem_proj_distinct_boundary_1x_p1_ho},
  {"fem_proj_distinct_boundary_1x_p2_ho", test_fem_proj_distinct_boundary_1x_p2_ho},
  {"fem_proj_distinct_boundary_2x_p1_ho", test_fem_proj_distinct_boundary_2x_p1_ho},
  {"fem_proj_distinct_boundary_2x_p2_ho", test_fem_proj_distinct_boundary_2x_p2_ho},
  {"fem_proj_poisson_comparison_ho", test_fem_proj_poisson_ho},
  {"fem_proj_bias_1x_ho", test_fem_proj_bias_1x_ho},
  {"fem_proj_bias_2x_ho", test_fem_proj_bias_2x_ho},
  {"fem_proj_invalid_input", test_fem_proj_invalid_input},
#ifdef GKYL_HAVE_CUDA
  {"fem_proj_solver_1x_p1_dev", test_fem_proj_solver_1x_p1_dev},
  {"fem_proj_solver_1x_p2_dev", test_fem_proj_solver_1x_p2_dev},
  {"fem_proj_solver_2x_p1_dev", test_fem_proj_solver_2x_p1_dev},
  {"fem_proj_solver_2x_p2_dev", test_fem_proj_solver_2x_p2_dev},
  {"fem_proj_distinct_boundary_1x_p1_dev", test_fem_proj_distinct_boundary_1x_p1_dev},
  {"fem_proj_distinct_boundary_1x_p2_dev", test_fem_proj_distinct_boundary_1x_p2_dev},
  {"fem_proj_distinct_boundary_2x_p1_dev", test_fem_proj_distinct_boundary_2x_p1_dev},
  {"fem_proj_distinct_boundary_2x_p2_dev", test_fem_proj_distinct_boundary_2x_p2_dev},
  {"fem_proj_bias_1x_dev", test_fem_proj_bias_1x_dev},
  {"fem_proj_bias_2x_dev", test_fem_proj_bias_2x_dev},
  {"fem_proj_poisson_comparison_dev", test_fem_proj_poisson_dev},
#endif
  {NULL, NULL}
};
