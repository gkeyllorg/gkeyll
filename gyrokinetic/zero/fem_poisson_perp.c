#include <gkyl_fem_poisson_perp.h>
#include <gkyl_fem_poisson_perp_priv.h>
#include <gkyl_array_reduce.h>

static void
fem_poisson_perp_bias_src_disabled(gkyl_fem_poisson_perp* up, struct gkyl_array *rhsin)
{
  // Do nothing.
}

static void
fem_poisson_perp_bias_src_enabled(gkyl_fem_poisson_perp* up, struct gkyl_array *rhsin)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(rhsin));

    gkyl_fem_poisson_perp_bias_src_enabled_cu(up, rhsin);
    return;
  }
#endif

  double *brhs_p = gkyl_array_fetch(up->brhs, 0);

  int idx0[GKYL_MAX_CDIM], idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
  while (gkyl_range_iter_next(&up->par_iter1d)) {
    long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

    gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
    while (gkyl_range_iter_next(&up->perp_iter2d)) {
      long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

      for (size_t d=0; d<up->ndim_perp; d++) idx1[d] = up->perp_iter2d.idx[d];
      idx1[up->pardir] = up->par_iter1d.idx[0];

      int keri = idx_to_inup_ker(up->ndim_perp, up->num_cells, up->perp_iter2d.idx);
      for (size_t d=0; d<up->ndim; d++) idx0[d] = idx1[d] - 1;
      up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

      long parProbOff = paridx*up->numnodes_global;

      for (int i=0; i<up->num_bias_line; i++) {
        // Index of the cell that abuts the line from below.
        struct gkyl_poisson_bias_line *bl = &up->bias_lines[i];
        int bl_idx_m[up->bl_ndim_perp];
        for (int d=0; d<up->bl_ndim_perp; d++) {
          int perp_dir = bl->perp_dirs[d];
          double dx = up->grid.dx[perp_dir];
          bl_idx_m[d] = (bl->perp_coords[d]-1e-3*dx - up->grid.lower[perp_dir])/dx+1;
        }

        if (
            ( idx1[bl->perp_dirs[0]] == bl_idx_m[0]   && idx1[bl->perp_dirs[1]] == bl_idx_m[1]   ) ||
            ( idx1[bl->perp_dirs[0]] == bl_idx_m[0]+1 && idx1[bl->perp_dirs[1]] == bl_idx_m[1]   ) ||
            ( idx1[bl->perp_dirs[0]] == bl_idx_m[0]   && idx1[bl->perp_dirs[1]] == bl_idx_m[1]+1 ) ||
            ( idx1[bl->perp_dirs[0]] == bl_idx_m[0]+1 && idx1[bl->perp_dirs[1]] == bl_idx_m[1]+1 )
           ) {
          int edge[2] = {
            -1+2*((bl_idx_m[0]+1)-idx1[bl->perp_dirs[0]]),
            -1+2*((bl_idx_m[1]+1)-idx1[bl->perp_dirs[1]]),
          };
          up->kernels->bias_src_ker[keri](edge, bl->perp_dirs, bl->val, parProbOff, up->globalidx, brhs_p);
        }
      }
    }

  }
}

struct gkyl_fem_poisson_perp*
gkyl_fem_poisson_perp_new(const struct gkyl_range *solve_range, const struct gkyl_rect_grid *grid,
  const struct gkyl_basis basis, struct gkyl_poisson_bc *bcs, struct gkyl_poisson_bias_line_list *bias_lines,
  struct gkyl_array *epsilon, struct gkyl_array *kSq, bool use_gpu)
{

  struct gkyl_fem_poisson_perp *up = gkyl_malloc(sizeof(struct gkyl_fem_poisson_perp));

  up->solve_range = solve_range;
  up->ndim = grid->ndim;
  up->ndim_perp = up->ndim-1;
  up->grid = *grid;
  up->num_basis =  basis.num_basis;
  up->basis_type = basis.b_type;
  up->poly_order = basis.poly_order;
  up->pardir = grid->ndim-1; // Assume parallel direction is always the last.
  up->basis = basis;
  up->use_gpu = use_gpu;
  up->epsilon = gkyl_array_acquire(epsilon);

  assert(up->ndim > 1);
  assert(up->epsilon->ncomp == (2*(up->ndim-1)-1)*basis.num_basis);

  // We assume epsilon and kSq live on the device, and we create a host-side
  // copies temporarily to compute the LHS matrix. This also works for CPU solves.
  struct gkyl_array *epsilon_ho = gkyl_array_new(GKYL_DOUBLE, up->epsilon->ncomp, up->epsilon->size);
  gkyl_array_copy(epsilon_ho, up->epsilon);
  struct gkyl_array *kSq_ho;
  if (kSq) {
    up->ishelmholtz = true;
    kSq_ho = gkyl_array_new(GKYL_DOUBLE, kSq->ncomp, kSq->size);
    gkyl_array_copy(kSq_ho, kSq);
  } else {
    up->ishelmholtz = false;
    kSq_ho = gkyl_array_new(GKYL_DOUBLE, up->num_basis, 1);
    gkyl_array_clear(kSq_ho, 0.0);

    up->kSq_null = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, kSq_ho->ncomp, kSq_ho->size)
                          : gkyl_array_acquire(kSq_ho);
    gkyl_array_clear(up->kSq_null, 0.0);
  }

  up->globalidx = gkyl_malloc(sizeof(long[up->num_basis])); // global index, one for each basis in a cell.

  for (int d=0; d<up->ndim; d++) up->num_cells[d] = up->solve_range->upper[d]-up->solve_range->lower[d]+1;

  // 2D range of perpendicular cells.
  gkyl_range_init(&up->perp_range2d, up->ndim_perp, up->solve_range->lower, up->solve_range->upper);
  // 1D range of parallel cells.
  int lower1d[] = {up->solve_range->lower[up->pardir]}, upper1d[] = {up->solve_range->upper[up->pardir]};
  gkyl_range_init(&up->par_range1d, 1, lower1d, upper1d);

  // Range of perpendicular cells at each parallel location.
  up->perp_range = (struct gkyl_range *) gkyl_malloc(up->par_range1d.volume * sizeof(struct gkyl_range));
  gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
  while (gkyl_range_iter_next(&up->par_iter1d)) {
    long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);
    int removeDim[] = {0,0,0},  loc[] = {0,0,0};
    removeDim[up->pardir] = 1;
    loc[up->pardir] = paridx;
    gkyl_range_deflate(&up->perp_range[paridx], up->solve_range, removeDim, loc);
  }
  // Range of parallel cells.
  int sublower[GKYL_MAX_CDIM], subupper[GKYL_MAX_CDIM];
  for (int d=0; d<up->ndim; d++) {
    sublower[d] = up->solve_range->lower[d];
    subupper[d] = up->solve_range->lower[d];
  }
  subupper[up->pardir] = up->solve_range->upper[up->pardir];
  gkyl_sub_range_init(&up->par_range, up->solve_range, sublower, subupper);

  // Prepare for periodic domain case.
  for (int d=0; d<up->ndim_perp; d++) {
    // Sanity check.
    if ((bcs->lo_type[d] == GKYL_POISSON_PERIODIC && bcs->up_type[d] != GKYL_POISSON_PERIODIC) ||
        (bcs->lo_type[d] != GKYL_POISSON_PERIODIC && bcs->up_type[d] == GKYL_POISSON_PERIODIC))
      assert(false);
  }
  for (int d=0; d<up->ndim_perp; d++) up->isdirperiodic[d] = bcs->lo_type[d] == GKYL_POISSON_PERIODIC;
  up->isdomperiodic = true;
  for (int d=0; d<up->ndim_perp; d++) up->isdomperiodic = up->isdomperiodic && up->isdirperiodic[d];
//  assert(up->isdomperiodic == false);  // MF 2023/06/29: there's an error in
//                                       // the periodic domain case I have not
//                                       // solved.

  if (up->isdomperiodic) {
#ifdef GKYL_HAVE_CUDA
    if (up->use_gpu) {
      up->rhs_cellavg = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, epsilon->size);
      up->rhs_avg_cu = (double*) gkyl_cu_malloc(sizeof(double));
    } else {
      up->rhs_cellavg = gkyl_array_new(GKYL_DOUBLE, 1, epsilon->size);
    }
#else
    up->rhs_cellavg = gkyl_array_new(GKYL_DOUBLE, 1, epsilon->size);
#endif
    up->rhs_avg = (double*) gkyl_malloc(sizeof(double));
    gkyl_array_clear(up->rhs_cellavg, 0.0);
    // Factor accounting for normalization when subtracting a constant from a
    // DG field and the 1/N to properly compute the volume averaged RHS.
    up->mavgfac = -pow(sqrt(2.),up->ndim)/up->perp_range2d.volume;
  }

  // Pack BC values into a single array for easier use in kernels.
  for (int d=0; d<up->ndim_perp; d++) {
    for (int k=0; k<6; k++) up->bcvals[d*2*3+k] = 0.0; // default. Not used in some cases (e.g. periodic).
    if (bcs->lo_type[d] != GKYL_POISSON_PERIODIC) {
      int vnum, voff;
      vnum = 1;
      voff = 2;
      for (int k=0; k<vnum; k++) up->bcvals[d*2*3+voff+k] = bcs->lo_value[d].v[k];

      vnum = 1;
      voff = 2;
      for (int k=0; k<vnum; k++) up->bcvals[d*2*3+voff+3+k] = bcs->up_value[d].v[k];
    }
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->bcvals_cu = (double *) gkyl_cu_malloc(sizeof(double[PERP_DIM_MAX*3*2]));
    gkyl_cu_memcpy(up->bcvals_cu, up->bcvals, sizeof(double[PERP_DIM_MAX*3*2]), GKYL_CU_MEMCPY_H2D);
  }
#endif

  // Compute the number of local and global nodes.
  up->numnodes_local = up->num_basis;
  up->numnodes_global = gkyl_fem_poisson_perp_global_num_nodes(up->ndim, up->poly_order, basis.b_type, up->num_cells, up->isdirperiodic);

  for (int d=0; d<up->ndim; d++) up->dx[d] = up->grid.dx[d];  // Cell lengths.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->dx_cu = (double *) gkyl_cu_malloc(sizeof(double[GKYL_MAX_CDIM]));
    gkyl_cu_memcpy(up->dx_cu, up->dx, sizeof(double[GKYL_MAX_CDIM]), GKYL_CU_MEMCPY_H2D);
  }
#endif

  up->brhs = gkyl_array_new(GKYL_DOUBLE, 1, up->numnodes_global*up->par_range.volume); // Global right side vector.

  up->kernels = gkyl_malloc(sizeof(struct gkyl_fem_poisson_perp_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_fem_poisson_perp_kernels));
  else
    up->kernels_cu = up->kernels;
#else
  up->kernels_cu = up->kernels;
#endif

  // Select local-to-global mapping kernels:
  fem_poisson_perp_choose_local2global_kernels(&basis, up->isdirperiodic, up->kernels->l2g);

  // Select lhs kernels:
  fem_poisson_perp_choose_lhs_kernels(&basis, bcs, up->kernels->lhsker);

  // Select rhs src kernels:
  fem_poisson_perp_choose_src_kernels(&basis, bcs, up->kernels->srcker);

  // Select sol kernel:
  up->kernels->solker = fem_poisson_perp_choose_sol_kernels(&basis);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    fem_poisson_perp_choose_kernels_cu(&basis, bcs, up->isdirperiodic, up->kernels_cu);
#endif

  // Copy the biasing line list (bias_lines) into this updater.
  up->num_bias_line = 0;
  bool *bl_in_solve_range = 0;
  if (bias_lines) {
    if (bias_lines->num_bias_line > 0) {
      // Check if any bias lines are in solve_range, and copy their info into updater.
      bl_in_solve_range = gkyl_malloc(bias_lines->num_bias_line * sizeof(bool));
      for (int i=0; i<bias_lines->num_bias_line; i++)
        bl_in_solve_range[i] = false;

      for (int i=0; i<bias_lines->num_bias_line; i++) {
        struct gkyl_poisson_bias_line *bl = &bias_lines->bl[i];

        // MF 2025/11/10: For now limit ourselves to lines perpendicular to x and z.
        up->bl_ndim_perp = 2;
        assert(bl->perp_dirs[0] == 0 && bl->perp_dirs[1] == up->ndim-1);

        double line_coords[up->ndim];
        for (int d=0; d<up->ndim; d++)
          line_coords[d] = grid->lower[d]+grid->dx[d]/2.0;

        for (int d=0; d<up->bl_ndim_perp; d++)
          line_coords[bl->perp_dirs[d]] = bl->perp_coords[d];

        // If biased line is at domain boundary, shift it minimally so it is inside the domain.
        for (int d=0; d<up->ndim; d++) {
          if (fabs(line_coords[d] - grid->lower[d]) < 1e-3*grid->dx[d]) {
            line_coords[d] += 1e-3*grid->dx[d];
          }
          if (fabs(line_coords[d] - grid->upper[d]) < 1e-3*grid->dx[d]) {
            line_coords[d] += -1e-3*grid->dx[d];
          }
        }

        bool pick_lower[3] = {true, true, true}; // If at a cell boundary, pick the cell lower than the biased line.
        int line_idx[GKYL_MAX_CDIM];
        gkyl_rect_grid_find_cell(grid, line_coords, pick_lower, (int[3]){-1,-1,-1}, line_idx);

        bl_in_solve_range[i] = gkyl_range_contains_idx(solve_range, line_idx);
        if (bl_in_solve_range[i])
          up->num_bias_line++;
      }

      if (up->num_bias_line) {
        // Copy biased lines in solve range into a temporary struct.
        size_t bl_sz = up->num_bias_line * sizeof(struct gkyl_poisson_bias_line);
        struct gkyl_poisson_bias_line *bias_lines_buff = gkyl_malloc(bl_sz);
        int blc = 0;
        for (int i=0; i<bias_lines->num_bias_line; i++) {
          if (bl_in_solve_range[i]) {
            struct gkyl_poisson_bias_line *bl = &bias_lines->bl[i];
            memcpy(&bias_lines_buff[blc], &bias_lines->bl[i], sizeof(struct gkyl_poisson_bias_line));
            blc++;
          }
        }

        if (up->use_gpu) {
          up->bias_lines = gkyl_cu_malloc(bl_sz);
          gkyl_cu_memcpy(up->bias_lines, bias_lines_buff, bl_sz, GKYL_CU_MEMCPY_H2D);
        }
        else {
          up->bias_lines = gkyl_malloc(bl_sz);
          memcpy(up->bias_lines, bias_lines_buff, bl_sz);
        }
        gkyl_free(bias_lines_buff);

        // Select biasing kernels:
        fem_poisson_perp_choose_bias_lhs_kernels(&basis, up->isdirperiodic, up->kernels->bias_lhs_ker);
        fem_poisson_perp_choose_bias_src_kernels(&basis, up->isdirperiodic, up->kernels->bias_src_ker);
      }
    }
  }

  // Create a linear Ax=B problem for each perp plane. Here A is the discrete (global)
  // matrix representation of the LHS of the perpendiculat Helmholtz equation.
  // cuSolverRF may support for A_i x_i = B_i, so we may revisit this
  // structure for the GPU solve.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->prob_cu = gkyl_culinsolver_prob_new(up->par_range.volume, up->numnodes_global, up->numnodes_global, 1);
  } else {
    up->prob = gkyl_superlu_prob_new(up->par_range.volume, up->numnodes_global, up->numnodes_global, 1);
  }
#else
  up->prob = gkyl_superlu_prob_new(up->par_range.volume, up->numnodes_global, up->numnodes_global, 1);
#endif

  up->tri = gkyl_malloc(up->par_range.volume*sizeof(struct gkyl_mat_triples *));
  for (size_t i=0; i<up->par_range.volume; i++) {
    up->tri[i] = gkyl_mat_triples_new(up->numnodes_global, up->numnodes_global);
#ifdef GKYL_HAVE_CUDA
    if (up->use_gpu) gkyl_mat_triples_set_rowmaj_order(up->tri[i]);
#endif
  }

  // Assign non-zero elements in A.
  int idx0[GKYL_MAX_CDIM],  idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
  while (gkyl_range_iter_next(&up->par_iter1d)) {
    long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

    gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
    while (gkyl_range_iter_next(&up->perp_iter2d)) {
      long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

      for (size_t d=0; d<up->ndim_perp; d++) idx1[d] = up->perp_iter2d.idx[d];
      idx1[up->pardir] = up->par_iter1d.idx[0];

      long linidx = gkyl_range_idx(up->solve_range, idx1);

      double *eps_p = gkyl_array_fetch(epsilon_ho, linidx);
      double *kSq_p = up->ishelmholtz? gkyl_array_fetch(kSq_ho, linidx) : gkyl_array_fetch(kSq_ho,0);

      int keri = idx_to_inup_ker(up->ndim_perp, up->num_cells, up->perp_iter2d.idx);
      for (size_t d=0; d<up->ndim; d++) idx0[d] = idx1[d] - 1;
      up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

      // Apply the -nabla . (epsilon*nabla_perp)-kSq stencil.
      keri = idx_to_inloup_ker(up->ndim_perp, up->num_cells, idx1);
      up->kernels->lhsker[keri](eps_p, kSq_p, up->dx, up->bcvals, up->globalidx, up->tri[paridx]);
    }
  }

  if (up->num_bias_line > 0) {
    // If biased lines are specified, replace the corresponding equation in the
    // linear system so it only has a 1.
    gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
    while (gkyl_range_iter_next(&up->par_iter1d)) {
      long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

      gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
      while (gkyl_range_iter_next(&up->perp_iter2d)) {
        long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

        for (size_t d=0; d<up->ndim_perp; d++) idx1[d] = up->perp_iter2d.idx[d];
        idx1[up->pardir] = up->par_iter1d.idx[0];

        int keri = idx_to_inup_ker(up->ndim_perp, up->num_cells, up->perp_iter2d.idx);
        for (size_t d=0; d<up->ndim; d++) idx0[d] = idx1[d] - 1;
        up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

        for (int i=0; i<bias_lines->num_bias_line; i++) {
          if (bl_in_solve_range[i]) {
            // Index of the cell that abuts the line from below.
            struct gkyl_poisson_bias_line *bl = &bias_lines->bl[i];
            int bl_idx_m[up->bl_ndim_perp];
            for (int d=0; d<up->bl_ndim_perp; d++) {
              int perp_dir = bl->perp_dirs[d];
              double dx = up->grid.dx[perp_dir];
              bl_idx_m[d] = (bl->perp_coords[d]-1e-3*dx - up->grid.lower[perp_dir])/dx+1;
            }

            if (
                ( idx1[bl->perp_dirs[0]] == bl_idx_m[0]   && idx1[bl->perp_dirs[1]] == bl_idx_m[1]   ) ||
                ( idx1[bl->perp_dirs[0]] == bl_idx_m[0]+1 && idx1[bl->perp_dirs[1]] == bl_idx_m[1]   ) ||
                ( idx1[bl->perp_dirs[0]] == bl_idx_m[0]   && idx1[bl->perp_dirs[1]] == bl_idx_m[1]+1 ) ||
                ( idx1[bl->perp_dirs[0]] == bl_idx_m[0]+1 && idx1[bl->perp_dirs[1]] == bl_idx_m[1]+1 )
               ) {
              int edge[2] = {
                -1+2*((bl_idx_m[0]+1)-idx1[bl->perp_dirs[0]]),
                -1+2*((bl_idx_m[1]+1)-idx1[bl->perp_dirs[1]]),
              };
              up->kernels->bias_lhs_ker[keri](edge, bl->perp_dirs, up->globalidx, up->tri[paridx]);
            }
          }
        }
      }
    }

    up->bias_line_src = fem_poisson_perp_bias_src_enabled;
  }
  else {
    up->bias_line_src = fem_poisson_perp_bias_src_disabled;
  }

  if (bl_in_solve_range)
    gkyl_free(bl_in_solve_range);

  if (!(up->use_gpu))
    gkyl_superlu_amat_from_triples(up->prob, up->tri);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_culinsolver_amat_from_triples(up->prob_cu, up->tri);
#endif

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    // Store offsets into csr_val in cudss_ops.cu so we can update the LHS matrix on the GPU.
    up->csr_val_idx = gkyl_array_cu_dev_new(GKYL_LONG, pow(up->basis.num_basis, 2), up->epsilon->size);
    struct gkyl_array *csr_val_idx_ho = gkyl_array_new(GKYL_LONG, up->csr_val_idx->ncomp, up->csr_val_idx->size);

    gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
    while (gkyl_range_iter_next(&up->par_iter1d)) {
      long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

      gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
      while (gkyl_range_iter_next(&up->perp_iter2d)) {
        long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

        for (size_t d=0; d<up->ndim_perp; d++) idx1[d] = up->perp_iter2d.idx[d];
        idx1[up->pardir] = up->par_iter1d.idx[0];

        long linidx = gkyl_range_idx(up->solve_range, idx1);

        long *csr_val_idx_p = gkyl_array_fetch(csr_val_idx_ho, linidx);

        int keri = idx_to_inup_ker(up->ndim_perp, up->num_cells, up->perp_iter2d.idx);
        for (size_t d=0; d<up->ndim; d++) idx0[d] = idx1[d] - 1;
        up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

        for (int k=0; k<up->basis.num_basis; k++) {
          for (int l=0; l<up->basis.num_basis; l++) {
            size_t nnz = gkyl_mat_triples_size(up->tri[paridx]); // Number of nonzero elements.
            // Given the global i,j (row-col) place in the LHS matrix, find the linear index into the mat_triples list. 
            // Here we do a brute-force search as we are unsure of the order, and this is done only once at t=0.
            long off = -1;
            gkyl_mat_triples_iter *mtt_iter = gkyl_mat_triples_iter_new(up->tri[paridx]);
            for (size_t m=0; m<nnz; ++m) {
              gkyl_mat_triples_iter_next(mtt_iter); // bump iterator.
              struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(mtt_iter);
              if ((up->globalidx[k] == mt.row) && (up->globalidx[l] == mt.col)) {
                off = m;
                break;
              }
            }
	          csr_val_idx_p[k*up->basis.num_basis+l] = paridx*nnz + off;
          }
        }
      }
    }

    gkyl_array_copy(up->csr_val_idx, csr_val_idx_ho);
    gkyl_array_release(csr_val_idx_ho);
  }
#endif

  if (up->use_gpu) {
    for (size_t i=0; i<up->par_range.volume; i++)
      gkyl_mat_triples_release(up->tri[i]);

    gkyl_free(up->tri);
  }

  gkyl_array_release(epsilon_ho);
  gkyl_array_release(kSq_ho);

  return up;
}

void
gkyl_fem_poisson_perp_set_rhs(gkyl_fem_poisson_perp *up, struct gkyl_array *rhsin)
{

  if (up->isdomperiodic && !(up->ishelmholtz)) {
    // Subtract the volume averaged RHS from the RHS.
    gkyl_array_clear(up->rhs_cellavg, 0.0);

    gkyl_dg_calc_average_range(&up->basis, 0, up->rhs_cellavg, 0, rhsin, *up->solve_range);

    gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
    while (gkyl_range_iter_next(&up->par_iter1d)) {
      long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

#ifdef GKYL_HAVE_CUDA
      if (up->use_gpu) {
        gkyl_array_reduce_range(up->rhs_avg_cu, up->rhs_cellavg, GKYL_SUM, &(up->perp_range[paridx]));
        gkyl_cu_memcpy(up->rhs_avg, up->rhs_avg_cu, sizeof(double), GKYL_CU_MEMCPY_D2H);
      } else {
        gkyl_array_reduce_range(up->rhs_avg, up->rhs_cellavg, GKYL_SUM, &(up->perp_range[paridx]));
      }
#else
      gkyl_array_reduce_range(up->rhs_avg, up->rhs_cellavg, GKYL_SUM, &(up->perp_range[paridx]));
#endif
      gkyl_array_shiftc_range(rhsin, up->mavgfac*up->rhs_avg[0], 0, &(up->perp_range[paridx]));
    }
  }


#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(rhsin));

    gkyl_fem_poisson_perp_set_rhs_cu(up, rhsin);
    return;
  }
#endif

  gkyl_array_clear(up->brhs, 0.0);
  double *brhs_p = gkyl_array_fetch(up->brhs, 0);

  int idx0[GKYL_MAX_CDIM], idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
  while (gkyl_range_iter_next(&up->par_iter1d)) {
    long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

    gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
    while (gkyl_range_iter_next(&up->perp_iter2d)) {
      long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

      for (size_t d=0; d<up->ndim_perp; d++) idx1[d] = up->perp_iter2d.idx[d];
      idx1[up->pardir] = up->par_iter1d.idx[0];

      long linidx = gkyl_range_idx(up->solve_range, idx1);

      double *eps_p = gkyl_array_fetch(up->epsilon, linidx);
      double *rhsin_p = gkyl_array_fetch(rhsin, linidx);

      int keri = idx_to_inup_ker(up->ndim_perp, up->num_cells, up->perp_iter2d.idx);
      for (size_t d=0; d<up->ndim; d++) idx0[d] = idx1[d] - 1;
      up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

      // Apply the RHS source stencil. It's mostly the mass matrix times a
      // modal-to-nodal operator times the source, modified by BCs in skin cells.
      keri = idx_to_inloup_ker(up->ndim_perp, up->num_cells, idx1);

      long parProbOff = paridx*up->numnodes_global;

      up->kernels->srcker[keri](eps_p, up->dx, rhsin_p, up->bcvals, parProbOff, up->globalidx, brhs_p);
    }

  }

  // Set the corresponding entries to the biasing potential.
  up->bias_line_src(up, rhsin);

  gkyl_superlu_brhs_from_array(up->prob, brhs_p);

}

void
gkyl_fem_poisson_perp_solve(gkyl_fem_poisson_perp *up, struct gkyl_array *phiout) {
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(phiout));
    gkyl_fem_poisson_perp_solve_cu(up, phiout);
    return;
  }
#endif

  gkyl_superlu_solve(up->prob);

  gkyl_array_clear(phiout, 0.0);

  int idx0[GKYL_MAX_CDIM], idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
  while (gkyl_range_iter_next(&up->par_iter1d)) {
    long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);
  
    gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
    while (gkyl_range_iter_next(&up->perp_iter2d)) {
      long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

      for (size_t d=0; d<up->ndim_perp; d++) idx1[d] = up->perp_iter2d.idx[d];
      idx1[up->pardir] = up->par_iter1d.idx[0];

      long linidx = gkyl_range_idx(up->solve_range, idx1);

      double *phiout_p = gkyl_array_fetch(phiout, linidx);

      int keri = idx_to_inup_ker(up->ndim_perp, up->num_cells, up->perp_iter2d.idx);
      for (size_t d=0; d<up->ndim; d++) idx0[d] = idx1[d]-1;
      up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

      long parProbOff = paridx*up->numnodes_global;

      up->kernels->solker(gkyl_superlu_get_rhs_ptr(up->prob, 0), parProbOff, up->globalidx, phiout_p);
    }
  }

}

void
gkyl_fem_poisson_perp_update_lhs(gkyl_fem_poisson_perp *up, struct gkyl_array *epsilon, struct gkyl_array *kSq)
{
  assert(up->num_bias_line == 0); // Have not accounted for bias in csr_val_idx.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(epsilon));
    if (up->ishelmholtz)
      assert(gkyl_array_is_cu_dev(kSq));

    gkyl_fem_poisson_perp_update_lhs_cu(up, epsilon, kSq);
    return;
  }
#endif

  int idx0[GKYL_MAX_CDIM],  idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
  while (gkyl_range_iter_next(&up->par_iter1d)) {
    long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

    // Clear mat triples.
    gkyl_mat_triples_clear(up->tri[paridx], 0.0);

    gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
    while (gkyl_range_iter_next(&up->perp_iter2d)) {
      long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

      for (size_t d=0; d<up->ndim_perp; d++) idx1[d] = up->perp_iter2d.idx[d];
      idx1[up->pardir] = up->par_iter1d.idx[0];

      long linidx = gkyl_range_idx(up->solve_range, idx1);

      double *eps_p = gkyl_array_fetch(epsilon, linidx);
      double *kSq_p = up->ishelmholtz? gkyl_array_fetch(kSq, linidx) : gkyl_array_fetch(up->kSq_null,0);

      int keri = idx_to_inup_ker(up->ndim_perp, up->num_cells, up->perp_iter2d.idx);
      for (size_t d=0; d<up->ndim; d++) idx0[d] = idx1[d] - 1;
      up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

      // Apply the -nabla . (epsilon*nabla_perp)-kSq stencil.
      keri = idx_to_inloup_ker(up->ndim_perp, up->num_cells, idx1);
      up->kernels->lhsker[keri](eps_p, kSq_p, up->dx, up->bcvals, up->globalidx, up->tri[paridx]);
    }
  }

  gkyl_superlu_amat_update_from_triples(up->prob, up->tri);
}

void gkyl_fem_poisson_perp_release(struct gkyl_fem_poisson_perp *up)
{
  if (up->isdomperiodic) {
    gkyl_array_release(up->rhs_cellavg);
    gkyl_free(up->rhs_avg);
  }

  if (!(up->ishelmholtz))
    gkyl_array_release(up->kSq_null);

  if (!(up->use_gpu)) {
    for (size_t i=0; i<up->par_range.volume; i++)
      gkyl_mat_triples_release(up->tri[i]);

    gkyl_free(up->tri);
    gkyl_superlu_prob_release(up->prob);
    if (up->num_bias_line > 0)
      gkyl_free(up->bias_lines);
  }

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->kernels_cu);
    gkyl_cu_free(up->dx_cu);
    if (up->isdomperiodic)
      gkyl_cu_free(up->rhs_avg_cu);

    gkyl_cu_free(up->bcvals_cu);
    gkyl_array_release(up->csr_val_idx);
    gkyl_culinsolver_prob_release(up->prob_cu);
    if (up->num_bias_line > 0)
      gkyl_cu_free(up->bias_lines);
  }
#endif

  gkyl_array_release(up->brhs);
  gkyl_free(up->kernels);
  gkyl_free(up->perp_range);
  gkyl_free(up->globalidx);
  gkyl_array_release(up->epsilon);
  gkyl_free(up);
}
