#include <gkyl_fem_poisson.h>
#include <gkyl_fem_poisson_priv.h>
#include <gkyl_array_reduce.h>
#include <math.h>

static void
fem_poisson_bias_src_disabled(gkyl_fem_poisson* up, struct gkyl_array *rhsin)
{
}

static void
fem_poisson_bias_src_enabled(gkyl_fem_poisson* up, struct gkyl_array *rhsin)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(rhsin));

    gkyl_fem_poisson_bias_src_enabled_cu(up, rhsin);
    return;
  }
#endif

  gkyl_range_iter_init(&up->solve_iter, up->solve_range);
  int idx0[GKYL_MAX_CDIM];
  double *brhs_p = gkyl_array_fetch(up->brhs, 0);
  while (gkyl_range_iter_next(&up->solve_iter)) {
    long linidx = gkyl_range_idx(up->solve_range, up->solve_iter.idx);

    int keri = idx_to_inup_ker(up->ndim, up->num_cells, up->solve_iter.idx);
    for (size_t d=0; d<up->ndim; d++) idx0[d] = up->solve_iter.idx[d]-1;
    up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

    for (int i=0; i<up->num_bias_plane; i++) {
      // Index of the cell that abuts the plane from below.
      struct gkyl_poisson_bias_plane *bp = &up->bias_planes[i];
      double dx = up->grid.dx[bp->dir];
      int bp_idx_m = (bp->loc-1e-3*dx - up->grid.lower[bp->dir])/dx+1;
      
      if (up->solve_iter.idx[bp->dir] == bp_idx_m || up->solve_iter.idx[bp->dir] == bp_idx_m+1) {
        up->kernels->bias_src_ker[keri](-1+2*((bp_idx_m+1)-up->solve_iter.idx[bp->dir]),
          bp->dir, bp->val, up->globalidx, brhs_p);
      }
    }
  }
}

struct gkyl_fem_poisson*
gkyl_fem_poisson_new(const struct gkyl_range *solve_range, const struct gkyl_rect_grid *grid, const struct gkyl_basis basis,
  struct gkyl_poisson_bc *bcs, struct gkyl_poisson_bias_plane_list *bias_planes, struct gkyl_array *epsilon,
  struct gkyl_array *kSq, bool is_epsilon_const, bool use_gpu)
{

  struct gkyl_fem_poisson *up = gkyl_malloc(sizeof(struct gkyl_fem_poisson));

  up->kernels = gkyl_malloc(sizeof(struct gkyl_fem_poisson_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_fem_poisson_kernels));
  else
    up->kernels_cu = up->kernels;
#else
  up->kernels_cu = up->kernels;
#endif

  up->has_lhs_apply = false;
  up->bcs = *bcs;
  // Keep a host-side copy of the bias plane list (used to build the mass
  // solver in the LHS apply).
  up->bias_list_ho.num_bias_plane = 0;
  up->bias_list_ho.bp = 0;
  if (bias_planes) {
    if (bias_planes->num_bias_plane > 0) {
      up->bias_list_ho.num_bias_plane = bias_planes->num_bias_plane;
      size_t bp_sz = bias_planes->num_bias_plane * sizeof(struct gkyl_poisson_bias_plane);
      up->bias_list_ho.bp = gkyl_malloc(bp_sz);
      memcpy(up->bias_list_ho.bp, bias_planes->bp, bp_sz);
    }
  }
  up->solve_range = solve_range;
  up->ndim = grid->ndim;
  up->grid = *grid;
  up->num_basis =  basis.num_basis;
  up->basis_type = basis.b_type;
  up->poly_order = basis.poly_order;
  up->basis = basis;
  up->use_gpu = use_gpu;

  // Factor accounting for normalization when subtracting a constant from a
  // DG field and the 1/N to properly compute the volume averaged RHS.
  up->mavgfac = -pow(sqrt(2.),up->ndim)/up->solve_range->volume;

  if (!is_epsilon_const) {
    up->isvareps = true;
    up->epsilon  = epsilon;
  } else {
    up->isvareps = false;
    // Create a small gkyl_array to hold the constant epsilon value.
    double *eps_avg = gkyl_malloc(sizeof(double)); 
#ifdef GKYL_HAVE_CUDA
    up->epsilon = up->use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 1) : gkyl_array_new(GKYL_DOUBLE, 1, 1);
    struct gkyl_array *eps_cellavg = up->use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, epsilon->size)
                                                : gkyl_array_new(GKYL_DOUBLE, 1, epsilon->size);
    gkyl_dg_calc_average_range(&up->basis, 0, eps_cellavg, 0, epsilon, *up->solve_range);
    if (up->use_gpu) {
      double *eps_avg_cu = gkyl_cu_malloc(sizeof(double));
      gkyl_array_reduce_range(eps_avg_cu, eps_cellavg, GKYL_SUM, up->solve_range);
      gkyl_cu_memcpy(eps_avg, eps_avg_cu, sizeof(double), GKYL_CU_MEMCPY_D2H);
      gkyl_cu_free(eps_avg_cu);
    } else {
      gkyl_array_reduce_range(eps_avg, eps_cellavg, GKYL_SUM, up->solve_range);
    }
#else
    up->epsilon = gkyl_array_new(GKYL_DOUBLE, 1, 1);
    struct gkyl_array *eps_cellavg = gkyl_array_new(GKYL_DOUBLE, 1, epsilon->size);

    gkyl_dg_calc_average_range(&up->basis, 0, eps_cellavg, 0, epsilon, *up->solve_range);
    gkyl_array_reduce_range(eps_avg, eps_cellavg, GKYL_SUM, up->solve_range);
#endif
    gkyl_array_shiftc(up->epsilon, eps_avg[0]/up->solve_range->volume, 0);
    gkyl_array_release(eps_cellavg);
    gkyl_free(eps_avg);
  }

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
    gkyl_array_clear(kSq_ho, 0.);
  }

  up->globalidx = gkyl_malloc(sizeof(long[up->num_basis])); // global index, one for each basis in a cell.

  for (int d=0; d<up->ndim; d++) up->num_cells[d] = up->solve_range->upper[d]-up->solve_range->lower[d]+1;

  // Prepare for periodic domain case.
  for (int d=0; d<up->ndim; d++) {
    // Sanity check.
    if ((bcs->lo_type[d] == GKYL_POISSON_PERIODIC && bcs->up_type[d] != GKYL_POISSON_PERIODIC) ||
        (bcs->lo_type[d] != GKYL_POISSON_PERIODIC && bcs->up_type[d] == GKYL_POISSON_PERIODIC))
      assert(false);
  }
  for (int d=0; d<up->ndim; d++) up->isdirperiodic[d] = bcs->lo_type[d] == GKYL_POISSON_PERIODIC;
  up->isdomperiodic = true;
  for (int d=0; d<up->ndim; d++) up->isdomperiodic = up->isdomperiodic && up->isdirperiodic[d];
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
  }

  // Pack BC values into a single array for easier use in kernels.
  for (int d=0; d<up->ndim; d++) {
    for (int k=0; k<6; k++) up->bcvals[d*2*3+k] = 0.0; // default. Not used in some cases (e.g. periodic).
    if (bcs->lo_type[d] != GKYL_POISSON_PERIODIC) {
      int vnum, voff;
      vnum = bcs->lo_type[d] == GKYL_POISSON_ROBIN ? 3 : 1;
      voff = bcs->lo_type[d] == GKYL_POISSON_ROBIN ? 0 : 2;
      for (int k=0; k<vnum; k++) up->bcvals[d*2*3+voff+k] = bcs->lo_value[d].v[k];

      vnum = bcs->up_type[d] == GKYL_POISSON_ROBIN ? 3 : 1;
      voff = bcs->up_type[d] == GKYL_POISSON_ROBIN ? 0 : 2;
      for (int k=0; k<vnum; k++) up->bcvals[d*2*3+voff+3+k] = bcs->up_value[d].v[k];
    }
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->bcvals_cu = (double *) gkyl_cu_malloc(sizeof(double[GKYL_MAX_CDIM*3*2]));
    gkyl_cu_memcpy(up->bcvals_cu, up->bcvals, sizeof(double[GKYL_MAX_CDIM*3*2]), GKYL_CU_MEMCPY_H2D);
  }
#endif
  
  // Check if one of the boundaries needs a spatially varying Dirichlet BC.
  up->isdirichletvar = false;
  for (int d=0; d<up->ndim; d++) up->isdirichletvar = up->isdirichletvar ||
    (bcs->lo_type[d] == GKYL_POISSON_DIRICHLET_VARYING || bcs->up_type[d] == GKYL_POISSON_DIRICHLET_VARYING);

  // Compute the number of local and global nodes.
  up->numnodes_local = up->num_basis;
  up->numnodes_global = gkyl_fem_poisson_global_num_nodes(up->ndim, up->poly_order, basis.b_type, up->num_cells, up->isdirperiodic);

  for (int d=0; d<up->ndim; d++) up->dx[d] = up->grid.dx[d];  // Cell lengths.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->dx_cu = (double *) gkyl_cu_malloc(sizeof(double[GKYL_MAX_CDIM]));
    gkyl_cu_memcpy(up->dx_cu, up->dx, sizeof(double[GKYL_MAX_CDIM]), GKYL_CU_MEMCPY_H2D);
  }
#endif

  up->brhs = gkyl_array_new(GKYL_DOUBLE, 1, up->numnodes_global); // Global right side vector.

  // Select local-to-global mapping kernels:
  fem_poisson_choose_local2global_kernels(&basis, up->isdirperiodic, up->kernels->l2g);

  // Select lhs kernels:
  fem_poisson_choose_lhs_kernels(&basis, bcs, up->isvareps, up->kernels->lhsker);

  // Select rhs src kernels:
  fem_poisson_choose_src_kernels(&basis, bcs, up->isvareps, up->kernels->srcker);

  // Select sol kernel:
  up->kernels->solker = fem_poisson_choose_sol_kernels(&basis);

  // Copy the biasing plane list (bias_planes) to the Poisson updater
  up->num_bias_plane = 0;
  if (bias_planes) {
    if (bias_planes->num_bias_plane > 0) {
      // Select biasing kernels:
      fem_poisson_choose_bias_lhs_kernels(&basis, up->isdirperiodic, up->kernels->bias_lhs_ker);
      fem_poisson_choose_bias_src_kernels(&basis, up->isdirperiodic, up->kernels->bias_src_ker);
      // Copy biased planes' info into updater.
      up->num_bias_plane = bias_planes->num_bias_plane;
      size_t bp_sz = bias_planes->num_bias_plane * sizeof(struct gkyl_poisson_bias_plane);
      if (up->use_gpu) {
        up->bias_planes = gkyl_cu_malloc(bias_planes->num_bias_plane * sizeof(struct gkyl_poisson_bias_plane));
        gkyl_cu_memcpy(up->bias_planes, bias_planes->bp, bp_sz, GKYL_CU_MEMCPY_H2D);
      }
      else {
        up->bias_planes = gkyl_malloc(bias_planes->num_bias_plane * sizeof(struct gkyl_poisson_bias_plane));
        memcpy(up->bias_planes, bias_planes->bp, bp_sz);
      }
    }
  }

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    fem_poisson_choose_kernels_cu(&basis, bcs, up->isvareps, up->isdirperiodic, up->kernels_cu);
#endif

  // Create a linear Ax=B problem. Here A is the discrete (global) matrix
  // representation of the LHS of the Helmholtz equation.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    up->prob_cu = gkyl_culinsolver_prob_new(1, up->numnodes_global, up->numnodes_global, 1);
  else
    up->prob = gkyl_superlu_prob_new(1, up->numnodes_global, up->numnodes_global, 1);
#else
  up->prob = gkyl_superlu_prob_new(1, up->numnodes_global, up->numnodes_global, 1);
#endif

  // Assign non-zero elements in A.
  struct gkyl_mat_triples **tri = gkyl_malloc(sizeof(struct gkyl_mat_triples *));
  tri[0] = gkyl_mat_triples_new(up->numnodes_global, up->numnodes_global);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) gkyl_mat_triples_set_rowmaj_order(tri[0]);
#endif
  gkyl_range_iter_init(&up->solve_iter, up->solve_range);
  int idx0[GKYL_MAX_CDIM];
  while (gkyl_range_iter_next(&up->solve_iter)) {
    long linidx = gkyl_range_idx(up->solve_range, up->solve_iter.idx);

    double *eps_p = up->isvareps? gkyl_array_fetch(epsilon_ho, linidx) : gkyl_array_fetch(epsilon_ho,0);
    double *kSq_p = up->ishelmholtz? gkyl_array_fetch(kSq_ho, linidx) : gkyl_array_fetch(kSq_ho,0);

    int keri = idx_to_inup_ker(up->ndim, up->num_cells, up->solve_iter.idx);
    for (size_t d=0; d<up->ndim; d++) idx0[d] = up->solve_iter.idx[d]-1;
    up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

    // Apply the -nabla . (epsilon*nabla)-kSq stencil.
    keri = idx_to_inloup_ker(up->ndim, up->num_cells, up->solve_iter.idx);
    up->kernels->lhsker[keri](eps_p, kSq_p, up->dx, up->bcvals, up->globalidx, tri[0]);
  }

  if (up->num_bias_plane > 0) {
    // If biased planes are specified, replace the corresponding equation in the
    // linear system so it only has a 1.
    gkyl_range_iter_init(&up->solve_iter, up->solve_range);
    while (gkyl_range_iter_next(&up->solve_iter)) {
      long linidx = gkyl_range_idx(up->solve_range, up->solve_iter.idx);

      int keri = idx_to_inup_ker(up->ndim, up->num_cells, up->solve_iter.idx);
      for (size_t d=0; d<up->ndim; d++) idx0[d] = up->solve_iter.idx[d]-1;
      up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

      for (int i=0; i<bias_planes->num_bias_plane; i++) {
        // Index of the cell that abuts the plane from below.
        struct gkyl_poisson_bias_plane *bp = &bias_planes->bp[i];
        double dx = up->grid.dx[bp->dir];
        int bp_idx_m = (bp->loc-1e-3*dx - up->grid.lower[bp->dir])/dx+1;
        
        if (up->solve_iter.idx[bp->dir] == bp_idx_m || up->solve_iter.idx[bp->dir] == bp_idx_m+1) {
          up->kernels->bias_lhs_ker[keri](-1+2*((bp_idx_m+1)-up->solve_iter.idx[bp->dir]),
            bp->dir, up->globalidx, tri[0]);
        }
      }
    }

    up->bias_plane_src = fem_poisson_bias_src_enabled;
  }
  else {
    up->bias_plane_src = fem_poisson_bias_src_disabled;
  }
  
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_culinsolver_amat_from_triples(up->prob_cu, tri);
  } else {
    gkyl_superlu_amat_from_triples(up->prob, tri);
  }
#else
  gkyl_superlu_amat_from_triples(up->prob, tri);
#endif

  gkyl_mat_triples_release(tri[0]);
  gkyl_free(tri);
  gkyl_array_release(epsilon_ho);
  gkyl_array_release(kSq_ho);

  return up;
}

void
gkyl_fem_poisson_set_rhs(gkyl_fem_poisson* up, struct gkyl_array *rhsin, const struct gkyl_array *phibc)
{

  if (up->isdomperiodic && !(up->ishelmholtz)) {
    // Subtract the volume averaged RHS from the RHS.
    gkyl_array_clear(up->rhs_cellavg, 0.0);
    gkyl_dg_calc_average_range(&up->basis, 0, up->rhs_cellavg, 0, rhsin, *up->solve_range);
#ifdef GKYL_HAVE_CUDA
    if (up->use_gpu) {
      gkyl_array_reduce_range(up->rhs_avg_cu, up->rhs_cellavg, GKYL_SUM, up->solve_range);
      gkyl_cu_memcpy(up->rhs_avg, up->rhs_avg_cu, sizeof(double), GKYL_CU_MEMCPY_D2H);
    } else {
      gkyl_array_reduce_range(up->rhs_avg, up->rhs_cellavg, GKYL_SUM, up->solve_range);
    }
#else
    gkyl_array_reduce_range(up->rhs_avg, up->rhs_cellavg, GKYL_SUM, up->solve_range);
#endif
    gkyl_array_shiftc(rhsin, up->mavgfac*up->rhs_avg[0], 0);
  }


#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(rhsin));
    if (phibc)
      assert(gkyl_array_is_cu_dev(phibc));

    gkyl_fem_poisson_set_rhs_cu(up, rhsin, phibc);
    return;
  }
#endif

  gkyl_array_clear(up->brhs, 0.0);

  gkyl_range_iter_init(&up->solve_iter, up->solve_range);
  int idx0[GKYL_MAX_CDIM];
  double *brhs_p = gkyl_array_fetch(up->brhs, 0);
  while (gkyl_range_iter_next(&up->solve_iter)) {
    long linidx = gkyl_range_idx(up->solve_range, up->solve_iter.idx);

    double *eps_p = up->isvareps? gkyl_array_fetch(up->epsilon, linidx) : gkyl_array_fetch(up->epsilon,0);
    double *rhsin_p = gkyl_array_fetch(rhsin, linidx);
    const double *phibc_p = up->isdirichletvar? gkyl_array_cfetch(phibc, linidx) : NULL;

    int keri = idx_to_inup_ker(up->ndim, up->num_cells, up->solve_iter.idx);
    for (size_t d=0; d<up->ndim; d++) idx0[d] = up->solve_iter.idx[d]-1;
    up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

    // Apply the RHS source stencil. It's mostly the mass matrix times a
    // modal-to-nodal operator times the source, modified by BCs in skin cells.
    keri = idx_to_inloup_ker(up->ndim, up->num_cells, up->solve_iter.idx);
    up->kernels->srcker[keri](eps_p, up->dx, rhsin_p, up->bcvals, phibc_p, up->globalidx, brhs_p);
  }

  // Set the corresponding entries to the biasing potential.
  up->bias_plane_src(up, rhsin);

  gkyl_superlu_brhs_from_array(up->prob, brhs_p);

}

void
gkyl_fem_poisson_solve(gkyl_fem_poisson* up, struct gkyl_array *phiout) {
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(phiout));
    gkyl_fem_poisson_solve_cu(up, phiout);
    return;
  }
#endif

  gkyl_superlu_solve(up->prob);

  gkyl_array_clear(phiout, 0.0);

  int idx0[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->solve_iter, up->solve_range);
  while (gkyl_range_iter_next(&up->solve_iter)) {
    long linidx = gkyl_range_idx(up->solve_range, up->solve_iter.idx);

    double *phiout_p = gkyl_array_fetch(phiout, linidx);

    int keri = idx_to_inup_ker(up->ndim, up->num_cells, up->solve_iter.idx);
    for (size_t d=0; d<up->ndim; d++) idx0[d] = up->solve_iter.idx[d]-1;
    up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

    up->kernels->solker(gkyl_superlu_get_rhs_ptr(up->prob, 0), up->globalidx, phiout_p);

  }

}

static void
fem_poisson_lhs_apply_init(gkyl_fem_poisson *up)
{
  // Build a mass-matrix (M) solver: same operator with epsilon=0 and a constant
  // kSq=-1 (so the LHS stencil reduces to the mass matrix).
  struct gkyl_array *eps_zero, *kSq_mass;
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    eps_zero = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->epsilon->ncomp, up->epsilon->size);
    kSq_mass = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_basis, up->epsilon->size);
  } else {
    eps_zero = gkyl_array_new(GKYL_DOUBLE, up->epsilon->ncomp, up->epsilon->size);
    kSq_mass = gkyl_array_new(GKYL_DOUBLE, up->num_basis, up->epsilon->size);
  }
#else
  eps_zero = gkyl_array_new(GKYL_DOUBLE, up->epsilon->ncomp, up->epsilon->size);
  kSq_mass = gkyl_array_new(GKYL_DOUBLE, up->num_basis, up->epsilon->size);
#endif
  gkyl_array_clear(eps_zero, 0.0);
  gkyl_array_clear(kSq_mass, 0.0);
  gkyl_array_shiftc(kSq_mass, -pow(sqrt(2.0), up->ndim), 0);

  // Bias the mass solver like the LHS so biased values pass through the apply.
  struct gkyl_poisson_bias_plane_list *bias_list = up->bias_list_ho.num_bias_plane > 0? &up->bias_list_ho : NULL;
  up->mass = gkyl_fem_poisson_new(up->solve_range, &up->grid, up->basis, &up->bcs, bias_list,
    eps_zero, kSq_mass, true, up->use_gpu);

  up->lhs_dual = gkyl_malloc(sizeof(double[up->numnodes_global]));
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    up->lhs_dual_cu = gkyl_cu_malloc(sizeof(double[up->numnodes_global]));
#endif

  gkyl_array_release(eps_zero);
  gkyl_array_release(kSq_mass);
  up->has_lhs_apply = true;
}

void
gkyl_fem_poisson_lhs_apply(gkyl_fem_poisson* up, struct gkyl_array *xin, struct gkyl_array *xout)
{
  assert(up->ishelmholtz);
  if (!up->has_lhs_apply) fem_poisson_lhs_apply_init(up);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_fem_poisson_lhs_apply_cu(up, xin, xout);
    return;
  }
#endif

  // Recover the global nodal vector x_nodal = M^{-1}*(M_src*xin).
  gkyl_fem_poisson_set_rhs(up->mass, xin, NULL);
  gkyl_superlu_solve(up->mass->prob);
  double *x_nodal = gkyl_superlu_get_rhs_ptr(up->mass->prob, 0);

  // Dual d = (M+K)*x_nodal, then solve M*z = d (-> modal xout).
  gkyl_superlu_mat_vec(up->prob, x_nodal, up->lhs_dual);
  gkyl_superlu_brhs_from_array(up->mass->prob, up->lhs_dual);
  gkyl_fem_poisson_solve(up->mass, xout);
}

void gkyl_fem_poisson_release(gkyl_fem_poisson *up)
{
  if (up->has_lhs_apply) {
    gkyl_fem_poisson_release(up->mass);
    gkyl_free(up->lhs_dual);
#ifdef GKYL_HAVE_CUDA
    if (up->use_gpu)
      gkyl_cu_free(up->lhs_dual_cu);
#endif
  }

  if (up->bias_list_ho.num_bias_plane > 0)
    gkyl_free(up->bias_list_ho.bp);

  if (up->isdomperiodic) {
    gkyl_array_release(up->rhs_cellavg);
    gkyl_free(up->rhs_avg);
  }

  if (!up->isvareps)
    gkyl_array_release(up->epsilon);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->kernels_cu);
    gkyl_cu_free(up->dx_cu);
    if (up->isdomperiodic) gkyl_cu_free(up->rhs_avg_cu);
    gkyl_cu_free(up->bcvals_cu);
    gkyl_culinsolver_prob_release(up->prob_cu);

    if (up->num_bias_plane > 0)
      gkyl_cu_free(up->bias_planes);
  } else {
    gkyl_superlu_prob_release(up->prob);

    if (up->num_bias_plane > 0)
      gkyl_free(up->bias_planes);
  }
#else
  gkyl_superlu_prob_release(up->prob);
  if (up->num_bias_plane > 0)
    gkyl_free(up->bias_planes);
#endif

  gkyl_free(up->globalidx);
  gkyl_array_release(up->brhs);
  gkyl_free(up->kernels);
  gkyl_free(up);
}
