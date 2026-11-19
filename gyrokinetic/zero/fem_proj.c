#include <gkyl_fem_proj_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <limits.h>
#include <math.h>
#include <string.h>

// Validate parent-array addressing as well as host/device placement.
static void
fem_proj_check_array(const struct gkyl_fem_proj *up, const struct gkyl_array *arr)
{
  if (!arr || arr->type != GKYL_DOUBLE || arr->ncomp != up->basis.num_basis ||
      arr->size < up->expect_arr_sz || gkyl_array_is_cu_dev(arr) != up->use_gpu) {
    gkyl_exit("fem_proj: incompatible DG array layout or memory location.");
  }
}

// Count global FEM nodes and check size limits before constructing the matrix.
static void
fem_proj_init_num_nodes(struct gkyl_fem_proj *up)
{
  int ndim = up->grid.ndim;
  bool periodic[GKYL_MAX_CDIM];
  // Bound the node count before evaluating the generated count formula.
  long max_nodes = 1;
  for (int d = 0; d < ndim; d++) {
    periodic[d] = (up->periodic_mask & (1 << d)) != 0;
    long extent = (long)up->basis.poly_order * up->num_cells[d] + 1;
    if (max_nodes > LONG_MAX / extent) {
      gkyl_exit("fem_proj: node grid is too large.");
    }
    max_nodes *= extent;
  }

  up->num_nodes = up->kernels.num_nodes(up->num_cells, periodic);

  if (up->num_nodes <= 0 || up->num_nodes > INT_MAX) {
    gkyl_exit("fem_proj: unsupported matrix size.");
  }
}

// No constraints beyond the ordinary boundary conditions.
static void
fem_proj_bias_src_disabled(struct gkyl_fem_proj *up, double *rhs)
{
}

// Assign fixed biases after all volume contributions have been accumulated.
static void
fem_proj_bias_src_enabled(struct gkyl_fem_proj *up, double *rhs)
{
  for (int b = 0; b < up->num_bias_face; b++) {
    const struct fem_proj_bias_face *bias = &up->bias_faces[b];
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &bias->range);
    while (gkyl_range_iter_next(&iter)) {
      long global_idx[FEM_PROJ_MAX_NODES];
      fem_proj_global_indices(&up->kernels, &up->range, up->num_cells, iter.idx, global_idx);
      up->kernels.bias_src[bias->face](bias->value, global_idx, rhs);
    }
  }
}

// Record the cells on one side of a plane, preserving parent-array strides.
static void
fem_proj_add_bias_face(struct gkyl_fem_proj *up, int dir, int cell, int side, double value)
{
  if (cell < up->range.lower[dir] || cell > up->range.upper[dir]) {
    return;
  }
  int lower[GKYL_MAX_CDIM], upper[GKYL_MAX_CDIM];
  gkyl_copy_int_arr(up->range.ndim, up->range.lower, lower);
  gkyl_copy_int_arr(up->range.ndim, up->range.upper, upper);
  lower[dir] = upper[dir] = cell;
  struct fem_proj_bias_face *bias = &up->bias_faces[up->num_bias_face++];
  gkyl_sub_range_init(&bias->range, &up->range, lower, upper);
  bias->face = 2 * dir + side;
  bias->value = value;
}

// Validate constant constraints sharing nodes with a new plane. Varying boundary
// values are supplied later and must satisfy the same consistency contract.
static void
fem_proj_check_bias_value(const struct gkyl_fem_proj *up, int dir, int face_idx, double value)
{
  for (int d = 0; d < up->range.ndim; d++) {
    for (int side = 0; side < 2; side++) {
      int type = side ? up->bcs.up_type[d] : up->bcs.lo_type[d];
      int boundary_idx = side ? up->range.upper[d] : up->range.lower[d] - 1;
      if (type == GKYL_FEM_PROJ_DIRICHLET && (d != dir || face_idx == boundary_idx) &&
          value != up->bc_values[2 * d + side]) {
        gkyl_exit("fem_proj: bias conflicts with a constant Dirichlet boundary.");
      }
    }
  }
  for (int b = 0; b < up->num_bias_face; b++) {
    const struct fem_proj_bias_face *bias = &up->bias_faces[b];
    int other_dir = bias->face / 2;
    int other_idx = bias->range.lower[other_dir] - !(bias->face % 2);
    if ((up->periodic_mask & (1 << other_dir)) && other_idx == up->range.upper[other_dir]) {
      other_idx = up->range.lower[other_dir] - 1;
    }
    if ((other_dir != dir || other_idx == face_idx) && bias->value != value) {
      gkyl_exit("fem_proj: intersecting bias planes prescribe different values.");
    }
  }
}

// Convert physical, grid-aligned planes into their incident cell faces once.
// Periodic end planes constrain both copies of the identified boundary nodes.
static void
fem_proj_init_bias(struct gkyl_fem_proj *up, const struct gkyl_fem_proj_bias_plane_list *list)
{
  up->bias_src = fem_proj_bias_src_disabled;
  if (!list || list->num_bias_plane == 0) {
    return;
  }
  if (list->num_bias_plane < 0 || list->num_bias_plane > INT_MAX / 2 || !list->bp) {
    gkyl_exit("fem_proj: invalid bias plane list.");
  }
  for (int d = 0; d < up->range.ndim; d++) {
    if (up->range.lower[d] < 1 || up->range.upper[d] > up->grid.cells[d]) {
      gkyl_exit("fem_proj: physical bias coordinates require solve indices within the grid.");
    }
  }
  up->bias_faces = gkyl_malloc(2 * (size_t)list->num_bias_plane * sizeof(*up->bias_faces));
  for (int b = 0; b < list->num_bias_plane; b++) {
    const struct gkyl_fem_proj_bias_plane *plane = &list->bp[b];
    int dir = plane->dir;
    if (dir < 0 || dir >= up->range.ndim || !isfinite(plane->loc) || !isfinite(plane->val)) {
      gkyl_exit("fem_proj: invalid bias direction, location or value.");
    }
    double position = (plane->loc - up->grid.lower[dir]) / up->grid.dx[dir];
    double tolerance = 1e-9;
    if (position < up->range.lower[dir] - 1 - tolerance ||
        position > up->range.upper[dir] + tolerance) {
      continue;
    }
    if (fabs(position - round(position)) > tolerance) {
      gkyl_exit("fem_proj: a bias plane must coincide with a grid face.");
    }
    int face_idx = (int)round(position);
    bool periodic = up->periodic_mask & (1 << dir);
    if (periodic && face_idx == up->range.upper[dir]) {
      face_idx = up->range.lower[dir] - 1;
    }
    fem_proj_check_bias_value(up, dir, face_idx, plane->val);
    fem_proj_add_bias_face(up, dir, face_idx, 1, plane->val);
    fem_proj_add_bias_face(up, dir, face_idx + 1, 0, plane->val);
    if (periodic && face_idx == up->range.lower[dir] - 1) {
      fem_proj_add_bias_face(up, dir, up->range.upper[dir], 1, plane->val);
    }
  }
  if (up->num_bias_face) {
    up->bias_src = fem_proj_bias_src_enabled;
#ifdef GKYL_HAVE_CUDA
    if (up->use_gpu) {
      up->bias_src = fem_proj_bias_src_cu;
    }
#endif
  }
}

// Clear every incident cell's contributions to constrained matrix rows. The
// generated kernels preserve unit diagonals even when local columns alias.
static void
fem_proj_bias_matrix(struct gkyl_fem_proj *up, struct gkyl_mat_triples *matrix)
{
  for (int b = 0; b < up->num_bias_face; b++) {
    const struct fem_proj_bias_face *bias = &up->bias_faces[b];
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &bias->range);
    while (gkyl_range_iter_next(&iter)) {
      long global_idx[FEM_PROJ_MAX_NODES];
      fem_proj_global_indices(&up->kernels, &up->range, up->num_cells, iter.idx, global_idx);
      up->kernels.bias_lhs[bias->face](global_idx, matrix);
    }
  }
}

// Matrix construction is on the host, including for a GPU solve.
static void
fem_proj_init_matrix(struct gkyl_fem_proj *up)
{
  struct gkyl_array *weight = NULL;
  if (up->weight_left) {
    weight = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->weight_left->size);
    gkyl_array_copy(weight, up->weight_left);
  }

  struct gkyl_mat_triples *matrix = gkyl_mat_triples_new(up->num_nodes, up->num_nodes);
  if (up->use_gpu) {
    gkyl_mat_triples_set_rowmaj_order(matrix);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->range, iter.idx);

    long global_idx[FEM_PROJ_MAX_NODES];
    fem_proj_global_indices(&up->kernels, &up->range, up->num_cells, iter.idx, global_idx);

    // Uniform cell volume factor cancels from the unconstrained equations.
    int location = fem_proj_cell_location(&up->range, iter.idx);
    up->kernels.mass[location](weight ? gkyl_array_cfetch(weight, loc) : NULL, global_idx, matrix);
  }

  fem_proj_bias_matrix(up, matrix);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->prob_cu = gkyl_culinsolver_prob_new(1, up->num_nodes, up->num_nodes, 1);
    gkyl_culinsolver_amat_from_triples(up->prob_cu, &matrix);
  } else
#endif
  {
    up->prob = gkyl_superlu_prob_new(1, up->num_nodes, up->num_nodes, 1);
    gkyl_superlu_amat_from_triples(up->prob, &matrix);
    gkyl_superlu_ludecomp(up->prob);
    up->rhs = gkyl_malloc(up->num_nodes * sizeof(double));
  }

  gkyl_mat_triples_release(matrix);
  if (weight) {
    gkyl_array_release(weight);
  }
}

gkyl_fem_proj *
gkyl_fem_proj_inew(const struct gkyl_fem_proj_inp *inp)
{
  // Consistency checks.
  if (!inp || !inp->grid || !inp->basis || !inp->solve_range) {
    gkyl_exit("fem_proj: grid, basis and solve range are required.");
  }
  int ndim = inp->grid->ndim;
  if (ndim == 3) {
    gkyl_exit("fem_proj: 3D kernels have not been generated yet.");
  }
  if (ndim < 1 || ndim > 2 || inp->basis->ndim != ndim || inp->solve_range->ndim != ndim ||
      inp->basis->b_type != GKYL_BASIS_MODAL_SERENDIPITY ||
      (inp->basis->poly_order != 1 && inp->basis->poly_order != 2)) {
    gkyl_exit("fem_proj: requires 1D/2D serendipity p1/p2.");
  }
  for (int d = 0; d < ndim; d++) {
    if (inp->solve_range->upper[d] < inp->solve_range->lower[d]) {
      gkyl_exit("fem_proj: empty solve range.");
    }
    int lo = inp->bcs.lo_type[d], hi = inp->bcs.up_type[d];
    if (lo < GKYL_FEM_PROJ_NONE || lo > GKYL_FEM_PROJ_DIRICHLET_VARYING ||
        hi < GKYL_FEM_PROJ_NONE || hi > GKYL_FEM_PROJ_DIRICHLET_VARYING ||
        ((lo == GKYL_FEM_PROJ_PERIODIC) != (hi == GKYL_FEM_PROJ_PERIODIC))) {
      gkyl_exit("fem_proj: invalid or unpaired periodic boundary condition.");
    }
    long count = (long)inp->solve_range->upper[d] - inp->solve_range->lower[d] + 1;
    if (count > INT_MAX) {
      gkyl_exit("fem_proj: too many cells in a direction.");
    }
  }
#ifndef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    gkyl_exit("fem_proj: CUDA support is not compiled.");
  }
#endif

  // Allocate updater and copy some inputs.
  struct gkyl_fem_proj *up = gkyl_malloc(sizeof(*up));

  *up = (struct gkyl_fem_proj){
    .grid = *inp->grid,
    .range = *inp->solve_range,
    .basis = *inp->basis,
    .bcs = inp->bcs,
    .use_gpu = inp->use_gpu,
  };

  for (int d = 0; d < ndim; d++) {
    int lo = inp->bcs.lo_type[d], hi = inp->bcs.up_type[d];
    up->varying |= lo == GKYL_FEM_PROJ_DIRICHLET_VARYING || hi == GKYL_FEM_PROJ_DIRICHLET_VARYING;

    long count = (long)up->range.upper[d] - up->range.lower[d] + 1;
    up->num_cells[d] = count;
    up->bc_values[2 * d] = inp->bcs.lo_value[d];
    up->bc_values[2 * d + 1] = inp->bcs.up_value[d];
    up->periodic_mask |= (lo == GKYL_FEM_PROJ_PERIODIC) << d;
  }

  // Expected size of arrays (e.g. weights, RHS source array).
  up->expect_arr_sz = gkyl_range_idx(&up->range, up->range.upper) + 1;

  if (inp->weight_left) {
    fem_proj_check_array(up, inp->weight_left);
    up->weight_left = gkyl_array_acquire(inp->weight_left);
  }
  if (inp->weight_right) {
    fem_proj_check_array(up, inp->weight_right);
    up->weight_right = gkyl_array_acquire(inp->weight_right);
  }

  up->kernels = fem_proj_choose_kernels(
    ndim, up->basis.poly_order, up->weight_left != NULL, up->weight_right != NULL,
    up->periodic_mask, &up->bcs
  );

  fem_proj_init_num_nodes(up);

  fem_proj_init_bias(up, inp->bias_plane_list);
  fem_proj_init_matrix(up);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    fem_proj_init_cu(up);
  }
#endif

  return up;
}

void
gkyl_fem_proj_set_rhs(
  gkyl_fem_proj *up, const struct gkyl_array *rhsin, const struct gkyl_array *boundary_values
)
{
  // Consistency check.
  fem_proj_check_array(up, rhsin);
  if (up->varying) {
    fem_proj_check_array(up, boundary_values);
  }

  up->rhs_ready = true;
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    fem_proj_set_rhs_cu(up, rhsin, boundary_values);
    return;
  }
#endif

  // Zero out rhs vector.
  memset(up->rhs, 0, up->num_nodes * sizeof(double));

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->range, iter.idx);

    long global_idx[FEM_PROJ_MAX_NODES];
    fem_proj_global_indices(&up->kernels, &up->range, up->num_cells, iter.idx, global_idx);

    int location = fem_proj_cell_location(&up->range, iter.idx);
    up->kernels.src[location](
      up->weight_right ? gkyl_array_cfetch(up->weight_right, loc) : NULL,
      gkyl_array_cfetch(rhsin, loc), up->bc_values,
      up->varying ? gkyl_array_cfetch(boundary_values, loc) : NULL, global_idx, up->rhs
    );
  }

  up->bias_src(up, up->rhs);
  gkyl_superlu_brhs_from_array(up->prob, up->rhs);
}

void
gkyl_fem_proj_solve(gkyl_fem_proj *up, struct gkyl_array *out)
{
  // Consistency check.
  fem_proj_check_array(up, out);
  if (!up->rhs_ready) {
    gkyl_exit("fem_proj: set_rhs must precede solve.");
  }

  up->rhs_ready = false;

  // Solve linear problem.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    fem_proj_solve_cu(up, out);
    return;
  }
#endif

  gkyl_superlu_solve(up->prob);
  const double *solution = gkyl_superlu_get_rhs_ptr(up->prob, 0);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->range, iter.idx);

    long global_idx[FEM_PROJ_MAX_NODES];
    fem_proj_global_indices(&up->kernels, &up->range, up->num_cells, iter.idx, global_idx);

    double nodal[FEM_PROJ_MAX_NODES];
    for (int n = 0; n < up->basis.num_basis; n++) {
      nodal[n] = solution[global_idx[n]];
    }
    up->kernels.sol(nodal, gkyl_array_fetch(out, loc));
  }
}

void
gkyl_fem_proj_release(gkyl_fem_proj *up)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_culinsolver_prob_release(up->prob_cu);
    gkyl_cu_free(up->kernels_cu);
    gkyl_cu_free(up->bc_values_cu);
    if (up->bias_faces_cu) {
      gkyl_cu_free(up->bias_faces_cu);
    }
  } else
#endif
  {
    gkyl_superlu_prob_release(up->prob);
    gkyl_free(up->rhs);
  }
  if (up->weight_left) {
    gkyl_array_release(up->weight_left);
  }
  if (up->weight_right) {
    gkyl_array_release(up->weight_right);
  }
  gkyl_free(up->bias_faces);
  gkyl_free(up);
}
