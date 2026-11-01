extern "C" {
#include <gkyl_fem_proj_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
}

__global__ static void
fem_proj_set_kernels(
  struct fem_proj_kernels *kernels, int ndim, int order, bool left, bool right, int periodic_mask,
  struct gkyl_fem_proj_bc bcs
)
{
  *kernels = fem_proj_choose_kernels(ndim, order, left, right, periodic_mask, &bcs);
}

void
fem_proj_init_cu(struct gkyl_fem_proj *up)
{
  if (up->num_bias_face) {
    for (int b = 0; b < up->num_bias_face; b++) {
      up->num_bias_blocks = GKYL_MAX2(up->num_bias_blocks, up->bias_faces[b].range.nblocks);
    }
    size_t size = up->num_bias_face * sizeof(*up->bias_faces);
    up->bias_faces_cu = (struct fem_proj_bias_face *)gkyl_cu_malloc(size);
    gkyl_cu_memcpy(up->bias_faces_cu, up->bias_faces, size, GKYL_CU_MEMCPY_H2D);
  }
  up->bc_values_cu = (double *)gkyl_cu_malloc(sizeof(up->bc_values));
  gkyl_cu_memcpy(up->bc_values_cu, up->bc_values, sizeof(up->bc_values), GKYL_CU_MEMCPY_H2D);

  up->kernels_cu = (struct fem_proj_kernels *)gkyl_cu_malloc(sizeof(*up->kernels_cu));
  fem_proj_set_kernels<<<1, 1>>>(
    up->kernels_cu, up->grid.ndim, up->basis.poly_order, up->weight_left != NULL,
    up->weight_right != NULL, up->periodic_mask, up->bcs
  );
}

__global__ static void
fem_proj_load_kernel(
  struct gkyl_range range, const struct fem_proj_kernels *kernels, const struct gkyl_array *weight,
  const struct gkyl_array *input, const double *bc_values, const struct gkyl_array *boundary_values,
  double *rhs
)
{
  int num_cells[GKYL_MAX_CDIM];
  for (int d = 0; d < range.ndim; d++) {
    num_cells[d] = range.upper[d] - range.lower[d] + 1;
  }

  for (long cell = threadIdx.x + blockIdx.x * blockDim.x; cell < range.volume;
       cell += blockDim.x * gridDim.x) {
    int idx[GKYL_MAX_CDIM];
    gkyl_sub_range_inv_idx(&range, cell, idx);

    long loc = gkyl_range_idx(&range, idx);

    long global_idx[FEM_PROJ_MAX_NODES];
    fem_proj_global_indices(kernels, &range, num_cells, idx, global_idx);

    int location = fem_proj_cell_location(&range, idx);
    kernels->src[location](
      weight ? (const double *)gkyl_array_cfetch(weight, loc) : NULL,
      (const double *)gkyl_array_cfetch(input, loc), bc_values,
      boundary_values ? (const double *)gkyl_array_cfetch(boundary_values, loc) : NULL, global_idx,
      rhs
    );
  }
}

// Separate launch ensures that no volume additions follow bias assignments.
__global__ static void
fem_proj_bias_kernel(
  struct gkyl_range range, const struct fem_proj_kernels *kernels, int num_bias_face,
  const struct fem_proj_bias_face *bias_faces, double *rhs
)
{
  int num_cells[GKYL_MAX_CDIM];
  for (int d = 0; d < range.ndim; d++) {
    num_cells[d] = range.upper[d] - range.lower[d] + 1;
  }
  for (int b = 0; b < num_bias_face; b++) {
    const struct fem_proj_bias_face *bias = &bias_faces[b];
    for (long cell = threadIdx.x + blockIdx.x * blockDim.x; cell < bias->range.volume;
         cell += blockDim.x * gridDim.x) {
      int idx[GKYL_MAX_CDIM];
      gkyl_sub_range_inv_idx(&bias->range, cell, idx);
      long global_idx[FEM_PROJ_MAX_NODES];
      fem_proj_global_indices(kernels, &range, num_cells, idx, global_idx);
      kernels->bias_src[bias->face](bias->value, global_idx, rhs);
    }
  }
}

void
fem_proj_bias_src_cu(struct gkyl_fem_proj *up, double *rhs)
{
  fem_proj_bias_kernel<<<up->num_bias_blocks, up->range.nthreads>>>(
    up->range, up->kernels_cu, up->num_bias_face, up->bias_faces_cu, rhs
  );
}

__global__ static void
fem_proj_solution_kernel(
  struct gkyl_range range, const struct fem_proj_kernels *kernels, const double *solution,
  struct gkyl_array *out, int nb
)
{
  int num_cells[GKYL_MAX_CDIM];
  for (int d = 0; d < range.ndim; d++) {
    num_cells[d] = range.upper[d] - range.lower[d] + 1;
  }

  for (long cell = threadIdx.x + blockIdx.x * blockDim.x; cell < range.volume;
       cell += blockDim.x * gridDim.x) {
    int idx[GKYL_MAX_CDIM];
    gkyl_sub_range_inv_idx(&range, cell, idx);

    long loc = gkyl_range_idx(&range, idx);

    long global_idx[FEM_PROJ_MAX_NODES];
    fem_proj_global_indices(kernels, &range, num_cells, idx, global_idx);

    double nodal[FEM_PROJ_MAX_NODES];
    for (int n = 0; n < nb; n++) {
      nodal[n] = solution[global_idx[n]];
    }
    kernels->sol(nodal, (double *)gkyl_array_fetch(out, loc));
  }
}

void
fem_proj_set_rhs_cu(
  struct gkyl_fem_proj *up, const struct gkyl_array *rhs, const struct gkyl_array *boundary_values
)
{
  gkyl_culinsolver_clear_rhs(up->prob_cu, 0.0);
  double *load = gkyl_culinsolver_get_rhs_ptr(up->prob_cu, 0);
  fem_proj_load_kernel<<<up->range.nblocks, up->range.nthreads>>>(
    up->range, up->kernels_cu, up->weight_right ? up->weight_right->on_dev : NULL, rhs->on_dev,
    up->bc_values_cu, up->varying ? boundary_values->on_dev : NULL, load
  );
  up->bias_src(up, load);
}

void
fem_proj_solve_cu(struct gkyl_fem_proj *up, struct gkyl_array *out)
{
  gkyl_culinsolver_solve(up->prob_cu);
  fem_proj_solution_kernel<<<up->range.nblocks, up->range.nthreads>>>(
    up->range, up->kernels_cu, gkyl_culinsolver_get_sol_ptr(up->prob_cu, 0), out->on_dev,
    up->basis.num_basis
  );
}
