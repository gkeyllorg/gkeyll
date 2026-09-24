#include <gkyl_fem_parproj.h>
#include <gkyl_fem_parproj_priv.h>

static void
fem_parproj_bias_src_disabled(gkyl_fem_parproj* up, const struct gkyl_array *rhsin)
{
  // Do nothing.
}

static void
fem_parproj_bias_src_enabled(gkyl_fem_parproj* up, const struct gkyl_array *rhsin)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(rhsin));

    gkyl_fem_parproj_bias_src_enabled_cu(up, rhsin);
    return;
  }
#endif

  double *brhs_p = gkyl_array_fetch(up->brhs, 0);

  int idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
  while (gkyl_range_iter_next(&up->perp_iter2d)) {
    long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

    gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
    while (gkyl_range_iter_next(&up->par_iter1d)) {
      long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

      for (size_t d=0; d<up->pardir; d++) idx1[d] = up->perp_iter2d.idx[d];
      idx1[up->pardir] = up->par_iter1d.idx[0];

      int keri = up->par_iter1d.idx[0] == up->parnum_cells? 1 : 0;
      up->kernels->l2g[keri](up->parnum_cells, paridx, up->globalidx);

      long perpProbOff = perpidx*up->numnodes_global;

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
          up->kernels->bias_src_ker[keri](edge, bl->perp_dirs, bl->val, perpProbOff, up->globalidx, brhs_p);
        }
      }
    }

  }
}

struct gkyl_fem_parproj*
gkyl_fem_parproj_new(const struct gkyl_range *solve_range, const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, enum gkyl_fem_parproj_bc_type bctype, struct gkyl_poisson_bias_line_list* bias_lines,
  const struct gkyl_array *weight_left, const struct gkyl_array *weight_right, bool use_gpu)
{
  struct gkyl_fem_parproj *up = gkyl_malloc(sizeof(struct gkyl_fem_parproj));

  up->grid = *grid;
  up->solve_range = solve_range;
  up->ndim = solve_range->ndim;
  up->num_basis  = basis->num_basis;
  up->basis_type = basis->b_type;
  up->poly_order = basis->poly_order;
  up->pardir = up->ndim-1; // Assume parallel direction is always the last.
  up->isperiodic = bctype == GKYL_FEM_PARPROJ_PERIODIC;
  up->isdirichlet = bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST || GKYL_FEM_PARPROJ_DIRICHLET_SKIN;
  up->use_gpu = use_gpu;

  up->has_weight_rhs = false;
  if (weight_right) {
    up->has_weight_rhs = true;
    up->weight_rhs = gkyl_array_acquire(weight_right);
  }

  bool has_weight_lhs = false;
  struct gkyl_array *weight_left_ho;
  if (weight_left) {
    has_weight_lhs = true;
    weight_left_ho = use_gpu? gkyl_array_new(GKYL_DOUBLE, weight_left->ncomp, weight_left->size)
                            : gkyl_array_acquire(weight_left);
    gkyl_array_copy(weight_left_ho, weight_left);
  }

  up->globalidx = gkyl_malloc(sizeof(long[up->num_basis]));

  // Range of parallel cells, as a sub-range of up->solve_range.
  struct gkyl_range par_range;
  int sublower[GKYL_MAX_CDIM], subupper[GKYL_MAX_CDIM];
  for (int d=0; d<up->ndim; d++) {
    sublower[d] = up->solve_range->lower[d];
    subupper[d] = up->solve_range->lower[d];
  }
  subupper[up->pardir] = up->solve_range->upper[up->pardir];
  gkyl_sub_range_init(&par_range, up->solve_range, sublower, subupper);
  up->parnum_cells = par_range.volume;

  // Range of perpendicular cells.
  struct gkyl_range perp_range;
  gkyl_range_shorten_from_above(&perp_range, up->solve_range, up->pardir, 1);

  // 1D range of parallel cells.
  int lower1d[] = {par_range.lower[up->pardir]}, upper1d[] = {par_range.upper[up->pardir]};
  gkyl_range_init(&up->par_range1d, 1, lower1d, upper1d);
  // 2D range of perpendicular cells.
  gkyl_range_init(&up->perp_range2d, up->ndim==3 ? 2 : 1, perp_range.lower, perp_range.upper);

  // Compute the number of local and global nodes.
  up->numnodes_local = up->num_basis;
  up->numnodes_global = gkyl_fem_parproj_global_num_nodes(basis, up->isperiodic, par_range.volume);

  up->brhs = gkyl_array_new(GKYL_DOUBLE, 1, up->numnodes_global*perp_range.volume); // Global right side vector.

  // Allocate struct holding kernel pointers.
  struct gkyl_fem_parproj_kernels *kernels_ho = gkyl_malloc(sizeof(struct gkyl_fem_parproj_kernels));
  if (!use_gpu)
    up->kernels = gkyl_malloc(sizeof(struct gkyl_fem_parproj_kernels));

#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels = gkyl_cu_malloc(sizeof(struct gkyl_fem_parproj_kernels));
#endif

  // Choose kernels.
  fem_parproj_choose_kernels(basis, has_weight_lhs, up->has_weight_rhs, bctype, use_gpu, up->kernels);

  // Select kernels for building LHS matrix on host:
  fem_parproj_choose_kernels(basis, has_weight_lhs, up->has_weight_rhs, bctype, false, kernels_ho);

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

        if (!bl_in_solve_range[i]) {
          // Check if line is on the upper cell boundary. If it is, search for the cell
          // again but with pick_lower=false.
          bool on_upper_cell_boundary = false;
          double xc[GKYL_MAX_CDIM];
          gkyl_rect_grid_cell_center(grid, line_idx, xc);
          for (int d=0; d<up->ndim; d++) {
            if (fabs(line_coords[d] - (xc[d]+0.5*grid->dx[d])) < 1e-3*grid->dx[d]) {
              on_upper_cell_boundary = true;
            }
          }
          if (on_upper_cell_boundary) {
            pick_lower[0] = pick_lower[1] = pick_lower[2] = false; // If at a cell boundary, pick the cell upper than the biased line.
            gkyl_rect_grid_find_cell(grid, line_coords, pick_lower, (int[3]){-1,-1,-1}, line_idx);
            bl_in_solve_range[i] = gkyl_range_contains_idx(solve_range, line_idx);
          }
        }

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
      }
    }
  }

  // We support two cases:
  //  a) No weight, or weight is a single number so we can divide the RHS by it.
  //     Then we solve Ax=B where A is the discrete FEM projection operator,
  //     and B is a matrix with a column for each perpendicular cell.
  //  b) There's a spatially dependent weight. Then we solve A_i x_i=B_i
  //     where there's a different A_i for each perp cell and B_i is a single
  //     column matrix.
  struct gkyl_range prob_range;
  int nrhs;
  if (up->ndim == 1) {
    nrhs = 1;
    gkyl_range_init(&prob_range, 1, &((int){1}), &((int){1}));
  }
  else {
    if (has_weight_lhs || up->num_bias_line) {
      nrhs = 1;
      gkyl_range_init(&prob_range, up->perp_range2d.ndim, up->perp_range2d.lower, up->perp_range2d.upper);
    }
    else {
      nrhs = perp_range.volume;
      gkyl_range_init(&prob_range, 1, &((int){1}), &((int){1}));
    }
  }

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    up->prob_cu = gkyl_culinsolver_prob_new(prob_range.volume, up->numnodes_global, up->numnodes_global, nrhs);
  else
    up->prob = gkyl_superlu_prob_new(prob_range.volume, up->numnodes_global, up->numnodes_global, nrhs);
#else
  up->prob = gkyl_superlu_prob_new(prob_range.volume, up->numnodes_global, up->numnodes_global, nrhs);
#endif

  // Assign non-zero elements in A.
  struct gkyl_mat_triples **tri = gkyl_malloc(prob_range.volume*sizeof(struct gkyl_mat_triples *));
  for (size_t i=0; i<prob_range.volume; i++) {
    tri[i] = gkyl_mat_triples_new(up->numnodes_global, up->numnodes_global);
#ifdef GKYL_HAVE_CUDA
    if (up->use_gpu) gkyl_mat_triples_set_rowmaj_order(tri[i]);
#endif
  }

  int idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->perp_iter2d, &prob_range);
  while (gkyl_range_iter_next(&up->perp_iter2d)) {
    long perpidx = gkyl_range_idx(&prob_range, up->perp_iter2d.idx);

    gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
    while (gkyl_range_iter_next(&up->par_iter1d)) {
      long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

      const double *wgt_p = NULL;
      if (has_weight_lhs) {
        for (size_t d=0; d<up->pardir; d++) idx1[d] = up->perp_iter2d.idx[d];
        idx1[up->pardir] = up->par_iter1d.idx[0];
        long linidx = gkyl_range_idx(up->solve_range, idx1);
        wgt_p = gkyl_array_cfetch(weight_left_ho, linidx);
      }

      int keri = up->par_iter1d.idx[0] == up->parnum_cells? 1 : 0;
      kernels_ho->l2g[keri](up->parnum_cells, paridx, up->globalidx);

      // Apply the wgt*phi*basis stencil.
      keri = idx_to_inloup_ker(up->parnum_cells, up->par_iter1d.idx[0]);
      kernels_ho->lhsker[keri](wgt_p, up->globalidx, tri[perpidx]);
    }
  }

  if (up->num_bias_line > 0) {
    // If biased lines are specified, replace the corresponding equation in the
    // linear system so it only has a 1.
    gkyl_range_iter_init(&up->perp_iter2d, &prob_range);
    while (gkyl_range_iter_next(&up->perp_iter2d)) {
      long perpidx = gkyl_range_idx(&prob_range, up->perp_iter2d.idx);

      gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
      while (gkyl_range_iter_next(&up->par_iter1d)) {
        long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

        for (size_t d=0; d<up->pardir; d++) idx1[d] = up->perp_iter2d.idx[d];
        idx1[up->pardir] = up->par_iter1d.idx[0];

        int keri = up->par_iter1d.idx[0] == up->parnum_cells? 1 : 0;
        kernels_ho->l2g[keri](up->parnum_cells, paridx, up->globalidx);

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
              kernels_ho->bias_lhs_ker[keri](edge, bl->perp_dirs, up->globalidx, tri[perpidx]);
            }
          }
        }
      }
    }

    up->bias_line_src = fem_parproj_bias_src_enabled;
  }
  else {
    up->bias_line_src = fem_parproj_bias_src_disabled;
  }

  if (bl_in_solve_range)
    gkyl_free(bl_in_solve_range);

  if (!(up->use_gpu))
    gkyl_superlu_amat_from_triples(up->prob, tri);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_culinsolver_amat_from_triples(up->prob_cu, tri);
#endif

  for (size_t i=0; i<prob_range.volume; i++)
    gkyl_mat_triples_release(tri[i]);
  gkyl_free(tri);
  gkyl_free(kernels_ho);

  if (weight_left)
    gkyl_array_release(weight_left_ho);

  return up;
}

void
gkyl_fem_parproj_set_rhs(struct gkyl_fem_parproj* up, const struct gkyl_array *rhsin, const struct gkyl_array *phibc)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(rhsin));
    if (phibc)
      assert(gkyl_array_is_cu_dev(phibc));

    gkyl_fem_parproj_set_rhs_cu(up, rhsin, phibc);
    return;
  }
#endif

  gkyl_array_clear(up->brhs, 0.0);
  double *brhs_p = gkyl_array_fetch(up->brhs, 0);

  int idx1[GKYL_MAX_CDIM], dirichlet_idx[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
  while (gkyl_range_iter_next(&up->perp_iter2d)) {
    long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

    gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
    while (gkyl_range_iter_next(&up->par_iter1d)) {
      long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

      for (size_t d=0; d<up->pardir; d++) idx1[d] = up->perp_iter2d.idx[d];
      idx1[up->pardir] = up->par_iter1d.idx[0];
      long linidx = gkyl_range_idx(up->solve_range, idx1);

      const double *wgt_p = up->has_weight_rhs? gkyl_array_cfetch(up->weight_rhs, linidx) : NULL;
      const double *rhsin_p = gkyl_array_cfetch(rhsin, linidx);
      const double *phibc_p = up->kernels->get_dirichlet_value(up->pardir, up->parnum_cells, idx1, up->solve_range, phibc);

      long perpProbOff = perpidx*up->numnodes_global;

      int keri = up->par_iter1d.idx[0] == up->parnum_cells? 1 : 0;
      up->kernels->l2g[keri](up->parnum_cells, paridx, up->globalidx);

      keri = idx_to_inloup_ker(up->parnum_cells, up->par_iter1d.idx[0]);
      up->kernels->srcker[keri](wgt_p, rhsin_p, phibc_p, perpProbOff, up->globalidx, brhs_p);
    }
  }

  // Set the corresponding entries to the biasing potential.
  up->bias_line_src(up, rhsin);

  gkyl_superlu_brhs_from_array(up->prob, brhs_p);

}

void
gkyl_fem_parproj_solve(struct gkyl_fem_parproj* up, struct gkyl_array *phiout) {
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(phiout));
    gkyl_fem_parproj_solve_cu(up, phiout);
    return;
  }
#endif

  gkyl_superlu_solve(up->prob);

  int idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
  while (gkyl_range_iter_next(&up->perp_iter2d)) {
    long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

    gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
    while (gkyl_range_iter_next(&up->par_iter1d)) {
      long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

      for (size_t d=0; d<up->pardir; d++) idx1[d] = up->perp_iter2d.idx[d];
      idx1[up->pardir] = up->par_iter1d.idx[0];
      long linidx = gkyl_range_idx(up->solve_range, idx1);

      double *phiout_p = gkyl_array_fetch(phiout, linidx);

      long perpProbOff = perpidx*up->numnodes_global;

      int keri = up->par_iter1d.idx[0] == up->parnum_cells? 1 : 0;
      up->kernels->l2g[keri](up->parnum_cells, paridx, up->globalidx);

      up->kernels->solker(gkyl_superlu_get_rhs_ptr(up->prob, 0), perpProbOff, up->globalidx, phiout_p);
    }
  }

}

void gkyl_fem_parproj_release(struct gkyl_fem_parproj *up)
{
  if (up->has_weight_rhs) {
    gkyl_array_release(up->weight_rhs);
  }
  if (!up->use_gpu) {
    gkyl_free(up->kernels);
    gkyl_superlu_prob_release(up->prob);
    if (up->num_bias_line > 0)
      gkyl_free(up->bias_lines);
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->kernels);
    gkyl_culinsolver_prob_release(up->prob_cu);
    if (up->num_bias_line > 0)
      gkyl_cu_free(up->bias_lines);
  }
#endif
  gkyl_array_release(up->brhs);
  gkyl_free(up->globalidx);
  gkyl_free(up);
}
