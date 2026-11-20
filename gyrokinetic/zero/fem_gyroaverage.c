#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_fem_gyroaverage.h>
#include <gkyl_fem_gyroaverage_priv.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_mat.h>
#include <gkyl_mat_triples.h>
#include <gkyl_superlu_ops.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static int
fem_gyroavg_extend_node(
  const struct gkyl_fem_gyroaverage *up, int dir, int node, int *ext_node, double *ext_fac
)
{
  // Nodes of the grid, and their factors, that give the value of the field at the node with
  // index 'node' in direction dir, which may be out of bounds. Periodic directions wrap around.
  // In non-periodic directions the field is extended beyond the boundary linearly about its
  // boundary value, f(bnd - s) = 2*f(bnd) - f(bnd + s), so that the gyroaverage of a linear
  // function is that function up to the boundary, or by reflection, f(bnd - s) = f(bnd + s).
  // Returns the number of nodes (1 or 2).
  int num_nodes = up->num_nodes_dir[dir];
  if (up->is_periodic[dir]) {
    ext_node[0] = ((node % num_nodes) + num_nodes) % num_nodes;
    ext_fac[0] = 1.0;
    return 1;
  }
  if (node > -1 && node < num_nodes) {
    ext_node[0] = node;
    ext_fac[0] = 1.0;
    return 1;
  }
  // Boundary node, and node mirrored about it (limited to the grid if the gyroradius is
  // larger than the domain).
  int bnd = node < 0 ? 0 : num_nodes - 1;
  int mirror = 2 * bnd - node;
  mirror = mirror < 0 ? 0 : (mirror > num_nodes - 1 ? num_nodes - 1 : mirror);
  if (up->reflect_bnd) {
    ext_node[0] = mirror;
    ext_fac[0] = 1.0;
    return 1;
  }
  ext_node[0] = bnd;
  ext_fac[0] = 2.0;
  ext_node[1] = mirror;
  ext_fac[1] = -1.0;
  return 2;
}

static long
fem_gyroavg_node_idx(const struct gkyl_fem_gyroaverage *up, const int *node)
{
  // Linear index of the node of the grid with indices 'node' in each direction (the node past
  // the last cell of a periodic direction is its first node).
  long linidx = 0;
  for (int d = 0; d < up->ndim; d++) {
    int num_nodes = up->num_nodes_dir[d];
    linidx = linidx * num_nodes + ((node[d] % num_nodes) + num_nodes) % num_nodes;
  }
  return linidx;
}

static void
fem_gyroavg_shift_weights(double shift, int *src_lo, double wgt[2][3])
{
  // Integrals over a unit cell, t in [0,1], of each linear nodal basis
  // function of the cell (l_0=1-t, l_1=t; first index of wgt) times the
  // piecewise linear nodal function of the nodes src_lo+{0,1,2} (second index
  // of wgt) evaluated at the shifted position t+shift. Node indices are
  // relative to the lower node of the cell. The integrand changes form at
  // t=1-frac, where t+shift crosses a node, so each integral has two parts.
  double shift_floor = floor(shift);
  *src_lo = (int)shift_floor;
  double frac = shift - shift_floor;
  double brk = 1.0 - frac;
  double brk_sq = brk * brk, brk_cu = brk_sq * brk;
  double frac_sq = frac * frac, frac_cu = frac_sq * frac;

  wgt[0][0] = 0.5 * brk_sq - brk_cu / 6.0;
  wgt[0][1] = (brk - brk_sq + brk_cu / 6.0) + (0.5 * frac_sq - frac_cu / 6.0);
  wgt[0][2] = frac_cu / 6.0;
  wgt[1][0] = brk_cu / 6.0;
  wgt[1][1] = (0.5 * brk_sq - brk_cu / 6.0) + (frac - frac_sq + frac_cu / 6.0);
  wgt[1][2] = 0.5 * frac_sq - frac_cu / 6.0;
}

static double
fem_gyroavg_wrap_angle(double angle)
{
  // Map an angle to [0, 2*pi).
  double wrapped = fmod(angle, 2.0 * M_PI);
  return wrapped < 0.0 ? wrapped + 2.0 * M_PI : wrapped;
}

static int
fem_gyroavg_compare_doubles(const void *lhs, const void *rhs)
{
  double lhs_val = *(const double *)lhs, rhs_val = *(const double *)rhs;
  return (lhs_val > rhs_val) - (lhs_val < rhs_val);
}

static int
fem_gyroavg_break_angles(double amp, double phase, double *angles)
{
  // Gyroangles at which the shift amp*cos(angle-phase) (in units of the cell
  // length) is an integer. The integrand of the gyroangle integral is smooth
  // between these angles. Returns the number of angles found.
  int num_angles = 0;
  if (amp > 0.0) {
    int max_cross = (int)floor(amp);
    for (int n = -max_cross; n <= max_cross; n++) {
      double cross = acos(n / amp);
      angles[num_angles++] = fem_gyroavg_wrap_angle(phase + cross);
      angles[num_angles++] = fem_gyroavg_wrap_angle(phase - cross);
    }
  }
  return num_angles;
}

static void
fem_gyroavg_csr_mat_vec(const struct fem_gyroavg_csr *mat, const double *xin, double *yout)
{
  // Compute yout = mat*xin.
  for (long i = 0; i < mat->num_rows; i++) {
    double row_sum = 0.0;
    for (long k = mat->row_ptr[i]; k < mat->row_ptr[i + 1]; k++) {
      row_sum += mat->vals[k] * xin[mat->col_idx[k]];
    }
    yout[i] = row_sum;
  }
}

static void
fem_gyroavg_csr_mat_trans_vec_accum(
  const struct fem_gyroavg_csr *mat, const double *xin, double *yout
)
{
  // Compute yout += mat^T*xin.
  for (long i = 0; i < mat->num_rows; i++) {
    for (long k = mat->row_ptr[i]; k < mat->row_ptr[i + 1]; k++) {
      yout[mat->col_idx[k]] += mat->vals[k] * xin[i];
    }
  }
}

static void
fem_gyroavg_csr_from_triples(struct fem_gyroavg_csr *mat, long num_rows, gkyl_mat_triples *tri)
{
  // Fill a CSR matrix from triples in row-major order.
  assert(gkyl_mat_triples_is_rowmaj(tri));
  mat->num_rows = num_rows;
  mat->num_nonzeros = gkyl_mat_triples_size(tri);
  mat->row_ptr = gkyl_calloc(num_rows + 1, sizeof(long));
  mat->col_idx = gkyl_malloc(sizeof(long[mat->num_nonzeros]));
  mat->vals = gkyl_malloc(sizeof(double[mat->num_nonzeros]));

  gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri);
  for (long k = 0; k < mat->num_nonzeros; k++) {
    gkyl_mat_triples_iter_next(iter);
    struct gkyl_mtriple entry = gkyl_mat_triples_iter_at(iter);
    mat->row_ptr[entry.row + 1] += 1;
    mat->col_idx[k] = entry.col;
    mat->vals[k] = entry.val;
  }
  gkyl_mat_triples_iter_release(iter);

  for (long i = 0; i < num_rows; i++) {
    mat->row_ptr[i + 1] += mat->row_ptr[i];
  }
}

static void
fem_gyroavg_mass_solve(struct gkyl_fem_gyroaverage *up, const double *rhs, double *sol)
{
  // Solve M*sol = rhs, M being the mass matrix of the nodal basis.
  gkyl_superlu_brhs_from_array(up->mass_prob, rhs);
  gkyl_superlu_solve(up->mass_prob);
  memcpy(sol, gkyl_superlu_get_rhs_ptr(up->mass_prob, 0), sizeof(double[up->num_nodes]));
}

static void
fem_gyroavg_dg_to_load(
  const struct gkyl_fem_gyroaverage *up, const struct gkyl_array *fdg, int comp_off, double *load
)
{
  // Compute the load vector int phi_i*f dX of a DG field f, phi_i being the
  // nodal basis. Uses that the modal basis is orthonormal in the unit cell.
  for (long i = 0; i < up->num_nodes; i++) {
    load[i] = 0.0;
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->local);
  long cell = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->local, iter.idx);
    const double *fdg_d = (const double *)gkyl_array_cfetch(fdg, linidx) + comp_off;
    const long *nodes = &up->cell_nodes[cell * up->num_basis];
    for (int a = 0; a < up->num_basis; a++) {
      double proj = 0.0;
      for (int k = 0; k < up->num_basis; k++) {
        proj += up->nod2mod[k][a] * fdg_d[k];
      }
      load[nodes[a]] += up->load_fac * proj;
    }
    cell++;
  }
}

static void
fem_gyroavg_nodal_to_dg(
  const struct gkyl_fem_gyroaverage *up, const double *nodal, struct gkyl_array *fdg, int comp_off
)
{
  // Compute the DG coefficients of the continuous field with nodal values 'nodal'.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->local);
  long cell = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->local, iter.idx);
    double *fdg_d = (double *)gkyl_array_fetch(fdg, linidx) + comp_off;
    const long *nodes = &up->cell_nodes[cell * up->num_basis];
    for (int k = 0; k < up->num_basis; k++) {
      double coeff = 0.0;
      for (int a = 0; a < up->num_basis; a++) {
        coeff += up->nod2mod[k][a] * nodal[nodes[a]];
      }
      fdg_d[k] = coeff;
    }
    cell++;
  }
}

static void
fem_gyroavg_mass_new(struct gkyl_fem_gyroaverage *up)
{
  // Assemble and factorize the mass matrix M_ij = int phi_i*phi_j dX.
  gkyl_mat_triples *tri = gkyl_mat_triples_new(up->num_nodes, up->num_nodes);
  gkyl_mat_triples_set_colmaj_order(tri);

  long num_cells_tot = up->local.volume;
  for (long cell = 0; cell < num_cells_tot; cell++) {
    const long *nodes = &up->cell_nodes[cell * up->num_basis];
    for (int a = 0; a < up->num_basis; a++) {
      for (int b = 0; b < up->num_basis; b++) {
        // Tensor product of 1D element mass matrices (dx/6)*[2 1; 1 2].
        double mass_elem = 1.0;
        for (int d = 0; d < up->ndim; d++) {
          bool same_node = up->node_offset[a][d] == up->node_offset[b][d];
          mass_elem *= up->grid.dx[d] * (same_node ? 2.0 : 1.0) / 6.0;
        }
        gkyl_mat_triples_accum(tri, nodes[a], nodes[b], mass_elem);
      }
    }
  }

  up->mass_prob = gkyl_superlu_prob_new(1, up->num_nodes, up->num_nodes, 1);
  gkyl_superlu_amat_from_triples(up->mass_prob, &tri);
  gkyl_superlu_ludecomp(up->mass_prob);
  gkyl_mat_triples_release(tri);
}

static void
fem_gyroavg_gyro_mat_new(
  struct gkyl_fem_gyroaverage *up, const struct gkyl_fem_gyroaverage_inp *inp, int mu_idx
)
{
  // Assemble the gyroaverage matrix of the mu_idx-th magnetic moment,
  //   B_ij = int phi_i(X) 1/(2*pi) int_0^{2*pi} phi_j(X + delta(X,alpha)) dalpha dX.
  // In each cell the shift delta is evaluated with the metric and magnetic
  // field at the cell center. For a given gyroangle the integral over the cell
  // is a product of 1D integrals with a closed form. The gyroangle integral
  // uses Gauss-Legendre quadrature between the angles at which the shift in
  // any direction is a multiple of the cell length, where the integrand is smooth.
  int ndim = up->ndim, num_basis = up->num_basis;
  double mu = inp->mu[mu_idx];
  assert(mu >= 0.0);

  gkyl_mat_triples *tri = gkyl_mat_triples_new(up->num_nodes, up->num_nodes);
  gkyl_mat_triples_set_rowmaj_order(tri);

  double cell_avg_fac = 1.0 / sqrt(pow(2.0, ndim));
  double cell_vol = up->grid.cellVolume;

  // Buffers for the break angles and for the weights of the source nodes
  // coupled to each node of a cell (grown as needed).
  long angles_cap = 0, weights_cap = 0;
  double *angles = 0, *weights = 0;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->local);
  long cell = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->local, iter.idx);
    const long *nodes = &up->cell_nodes[cell * num_basis];

    // Shift in each direction, in units of the cell length, is
    // amp[d]*cos(alpha - phase[d]).
    const double *gxx_d = gkyl_array_cfetch(inp->gxx, linidx);
    const double *bmag_d = gkyl_array_cfetch(inp->bmag, linidx);
    double gxx_c = gxx_d[0] * cell_avg_fac, bmag_c = bmag_d[0] * cell_avg_fac;
    assert(gxx_c > 0.0 && bmag_c > 0.0);
    double rho = sqrt(2.0 * inp->mass * mu / bmag_c) / fabs(inp->charge);

    double amp[FEM_GYROAVG_MAX_DIM] = {0.0}, phase[FEM_GYROAVG_MAX_DIM] = {0.0};
    amp[0] = rho * sqrt(gxx_c) / up->grid.dx[0];
    if (ndim == 2) {
      const double *gxy_d = gkyl_array_cfetch(inp->gxy, linidx);
      const double *gyy_d = gkyl_array_cfetch(inp->gyy, linidx);
      double gxy_c = gxy_d[0] * cell_avg_fac, gyy_c = gyy_d[0] * cell_avg_fac;
      double metric_det = gxx_c * gyy_c - gxy_c * gxy_c;
      assert(metric_det > 0.0);
      double amp_cos = rho * gxy_c / sqrt(gxx_c) / up->grid.dx[1];
      double amp_sin = rho * sqrt(metric_det / gxx_c) / up->grid.dx[1];
      amp[1] = hypot(amp_cos, amp_sin);
      phase[1] = atan2(amp_sin, amp_cos);
    }

    // Source nodes coupled to this cell, relative to its lower node, span
    // src_lo[d] to src_lo[d]+src_num[d]-1 in each direction.
    int src_lo[FEM_GYROAVG_MAX_DIM] = {0}, src_num[FEM_GYROAVG_MAX_DIM] = {1, 1};
    long angles_num_max = 2, src_num_tot = 1;
    for (int d = 0; d < ndim; d++) {
      int reach = (int)ceil(amp[d]) + 1;
      src_lo[d] = -reach;
      src_num[d] = 2 * reach + 3;
      src_num_tot *= src_num[d];
      angles_num_max += 2 * (2 * (long)floor(amp[d]) + 1);
    }
    if (angles_num_max > angles_cap) {
      angles_cap = angles_num_max;
      angles = gkyl_realloc(angles, sizeof(double[angles_cap]));
    }
    if (num_basis * src_num_tot > weights_cap) {
      weights_cap = num_basis * src_num_tot;
      weights = gkyl_realloc(weights, sizeof(double[weights_cap]));
    }
    for (long i = 0; i < num_basis * src_num_tot; i++) {
      weights[i] = 0.0;
    }

    // Gyroangle intervals in which the integrand is smooth.
    int num_angles = 0;
    angles[num_angles++] = 0.0;
    angles[num_angles++] = 2.0 * M_PI;
    for (int d = 0; d < ndim; d++) {
      num_angles += fem_gyroavg_break_angles(amp[d], phase[d], &angles[num_angles]);
    }
    qsort(angles, num_angles, sizeof(double), fem_gyroavg_compare_doubles);

    bool is_shifted = amp[0] > 0.0 || amp[1] > 0.0;
    int num_intervals = is_shifted ? num_angles - 1 : 1;
    for (int interval = 0; interval < num_intervals; interval++) {
      // Without a shift the integrand does not depend on the gyroangle.
      double angle_lo = is_shifted ? angles[interval] : 0.0;
      double angle_up = is_shifted ? angles[interval + 1] : 2.0 * M_PI;
      double arc = angle_up - angle_lo;
      if (arc <= 0.0) {
        continue;
      }
      int num_sub = is_shifted ? (int)ceil(arc / FEM_GYROAVG_MAX_ARC) : 1;
      int num_quad = is_shifted ? 8 : 1;
      double sub_arc = arc / num_sub;
      for (int sub = 0; sub < num_sub; sub++) {
        for (int quad = 0; quad < num_quad; quad++) {
          double ordinate = is_shifted ? gkyl_gauss_ordinates_8[quad] : 0.0;
          double quad_wgt = is_shifted ? gkyl_gauss_weights_8[quad] : 2.0;
          double angle = angle_lo + sub_arc * (sub + 0.5 * (ordinate + 1.0));
          double angle_wgt = 0.5 * sub_arc * quad_wgt / (2.0 * M_PI);

          int shift_lo[FEM_GYROAVG_MAX_DIM] = {0};
          double shift_wgt[FEM_GYROAVG_MAX_DIM][2][3] = {{{0.0}}};
          for (int d = 0; d < ndim; d++) {
            double shift = amp[d] * cos(angle - phase[d]);
            fem_gyroavg_shift_weights(shift, &shift_lo[d], shift_wgt[d]);
          }

          for (int a = 0; a < num_basis; a++) {
            double *weights_a = &weights[a * src_num_tot];
            const double *wgt_x = shift_wgt[0][up->node_offset[a][0]];
            if (ndim == 1) {
              for (int rx = 0; rx < 3; rx++) {
                weights_a[shift_lo[0] + rx - src_lo[0]] += angle_wgt * wgt_x[rx];
              }
            } else {
              const double *wgt_y = shift_wgt[1][up->node_offset[a][1]];
              for (int rx = 0; rx < 3; rx++) {
                long src_x = shift_lo[0] + rx - src_lo[0];
                for (int ry = 0; ry < 3; ry++) {
                  long src_y = shift_lo[1] + ry - src_lo[1];
                  weights_a[src_x * src_num[1] + src_y] += angle_wgt * wgt_x[rx] * wgt_y[ry];
                }
              }
            }
          }
        }
      }
    }

    // Add the contribution of this cell to the matrix.
    int cell_idx[FEM_GYROAVG_MAX_DIM] = {0};
    for (int d = 0; d < ndim; d++) {
      cell_idx[d] = iter.idx[d] - up->local.lower[d];
    }
    for (int a = 0; a < num_basis; a++) {
      const double *weights_a = &weights[a * src_num_tot];
      for (long src = 0; src < src_num_tot; src++) {
        if (weights_a[src] == 0.0) {
          continue;
        }
        int src_node[FEM_GYROAVG_MAX_DIM] = {0};
        if (ndim == 1) {
          src_node[0] = cell_idx[0] + src_lo[0] + (int)src;
        } else {
          src_node[0] = cell_idx[0] + src_lo[0] + (int)(src / src_num[1]);
          src_node[1] = cell_idx[1] + src_lo[1] + (int)(src % src_num[1]);
        }
        // Nodes of the grid that give the value at the source node (which may be out of bounds).
        int num_ext[FEM_GYROAVG_MAX_DIM] = {1, 1};
        int ext_node[FEM_GYROAVG_MAX_DIM][2] = {{0}};
        double ext_fac[FEM_GYROAVG_MAX_DIM][2] = {{1.0, 0.0}, {1.0, 0.0}};
        for (int d = 0; d < ndim; d++) {
          num_ext[d] = fem_gyroavg_extend_node(up, d, src_node[d], ext_node[d], ext_fac[d]);
        }
        for (int ex = 0; ex < num_ext[0]; ex++) {
          for (int ey = 0; ey < num_ext[1]; ey++) {
            int grid_node[FEM_GYROAVG_MAX_DIM] = {ext_node[0][ex], ext_node[1][ey]};
            gkyl_mat_triples_accum(
              tri, nodes[a], fem_gyroavg_node_idx(up, grid_node),
              cell_vol * weights_a[src] * ext_fac[0][ex] * ext_fac[1][ey]
            );
          }
        }
      }
    }
    cell++;
  }

  fem_gyroavg_csr_from_triples(&up->gyro_mat[mu_idx], up->num_nodes, tri);

  gkyl_free(angles);
  gkyl_free(weights);
  gkyl_mat_triples_release(tri);
}

static void
fem_gyroavg_csr_release(struct fem_gyroavg_csr *mat)
{
  gkyl_free(mat->row_ptr);
  gkyl_free(mat->col_idx);
  gkyl_free(mat->vals);
}

static void
fem_gyroavg_mass_inv_1d(const struct gkyl_fem_gyroaverage *up, int dir, double *mass_inv)
{
  // Compute the inverse (row-major) of the mass matrix of the 1D nodal basis in direction dir.
  // The mass matrix of the plane is the tensor product of these 1D mass matrices.
  int num_nodes = up->num_nodes_dir[dir];
  struct gkyl_mat *mass = gkyl_mat_new(num_nodes, num_nodes, 0.0);
  struct gkyl_mat *inv = gkyl_mat_new(num_nodes, num_nodes, 0.0);
  for (int cell = 0; cell < up->num_cells[dir]; cell++) {
    // 1D element mass matrix (dx/6)*[2 1; 1 2]; nodes wrap around in periodic directions.
    for (int a = 0; a < 2; a++) {
      for (int b = 0; b < 2; b++) {
        int node_a = (cell + a) % num_nodes, node_b = (cell + b) % num_nodes;
        gkyl_mat_inc(mass, node_a, node_b, up->grid.dx[dir] * (a == b ? 2.0 : 1.0) / 6.0);
      }
    }
  }
  for (int i = 0; i < num_nodes; i++) {
    gkyl_mat_set(inv, i, i, 1.0);
  }

  gkyl_mem_buff ipiv = gkyl_mem_buff_new(sizeof(long[num_nodes]));
  bool status = gkyl_mat_linsolve_lu(mass, inv, gkyl_mem_buff_data(ipiv));
  assert(status);
  for (int i = 0; i < num_nodes; i++) {
    for (int j = 0; j < num_nodes; j++) {
      mass_inv[i * num_nodes + j] = gkyl_mat_get(inv, i, j);
    }
  }

  gkyl_mem_buff_release(ipiv);
  gkyl_mat_release(mass);
  gkyl_mat_release(inv);
}

static void *
fem_gyroavg_dev_copy(const void *data, size_t size)
{
  // Copy host data to newly allocated device memory.
  void *data_dev = gkyl_cu_malloc(size);
  gkyl_cu_memcpy(data_dev, data, size, GKYL_CU_MEMCPY_H2D);
  return data_dev;
}

static void
fem_gyroavg_csr_dev_copy(struct fem_gyroavg_csr *mat_dev, const struct fem_gyroavg_csr *mat)
{
  // Copy a CSR matrix to the device.
  mat_dev->num_rows = mat->num_rows;
  mat_dev->num_nonzeros = mat->num_nonzeros;
  mat_dev->row_ptr = fem_gyroavg_dev_copy(mat->row_ptr, sizeof(long[mat->num_rows + 1]));
  mat_dev->col_idx = fem_gyroavg_dev_copy(mat->col_idx, sizeof(long[mat->num_nonzeros]));
  mat_dev->vals = fem_gyroavg_dev_copy(mat->vals, sizeof(double[mat->num_nonzeros]));
}

static void
fem_gyroavg_csr_dev_release(struct fem_gyroavg_csr *mat_dev)
{
  gkyl_cu_free(mat_dev->row_ptr);
  gkyl_cu_free(mat_dev->col_idx);
  gkyl_cu_free(mat_dev->vals);
}

struct gkyl_fem_gyroaverage_batch *
gkyl_fem_gyroaverage_batch_new(struct gkyl_fem_gyroaverage *const *ups, int num_fld)
{
  // Create the data used to apply the updaters on the GPU from their host matrices.
  struct gkyl_fem_gyroaverage_batch *batch = gkyl_malloc(sizeof(*batch));
  const struct gkyl_fem_gyroaverage *up = ups[0];
  int ndim = up->ndim, num_basis = up->num_basis, num_mu = up->num_mu;
  long num_nodes = up->num_nodes, num_cells = up->local.volume;

  batch->ndim = ndim;
  batch->num_basis = num_basis;
  batch->num_mu = num_mu;
  batch->num_fld = num_fld;
  batch->num_nodes = num_nodes;
  batch->num_cells = num_cells;
  for (int s = 0; s < num_fld; s++) {
    assert(ups[s]->num_nodes == num_nodes && ups[s]->num_mu == num_mu);
  }

  // Load matrix (see fem_gyroavg_dg_to_load) and linear index of each cell.
  long *cell_linidx = gkyl_malloc(sizeof(long[num_cells]));
  gkyl_mat_triples *tri = gkyl_mat_triples_new(num_nodes, num_cells * num_basis);
  gkyl_mat_triples_set_rowmaj_order(tri);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->local);
  long cell = 0;
  while (gkyl_range_iter_next(&iter)) {
    cell_linidx[cell] = gkyl_range_idx(&up->local, iter.idx);
    const long *nodes = &up->cell_nodes[cell * num_basis];
    for (int a = 0; a < num_basis; a++) {
      for (int k = 0; k < num_basis; k++) {
        gkyl_mat_triples_accum(
          tri, nodes[a], cell * num_basis + k, up->load_fac * up->nod2mod[k][a]
        );
      }
    }
    cell++;
  }
  struct fem_gyroavg_csr load_mat;
  fem_gyroavg_csr_from_triples(&load_mat, num_nodes, tri);
  gkyl_mat_triples_release(tri);

  // Gyroaverage matrices of all fields and magnetic moments stacked, and their transposes.
  long num_nonzeros = 0;
  for (int s = 0; s < num_fld; s++) {
    for (int n = 0; n < num_mu; n++) {
      num_nonzeros += ups[s]->gyro_mat[n].num_nonzeros;
    }
  }
  struct fem_gyroavg_csr gyro_mat = {
    .num_rows = num_fld * num_mu * num_nodes,
    .num_nonzeros = num_nonzeros,
    .row_ptr = gkyl_calloc(num_fld * num_mu * num_nodes + 1, sizeof(long)),
    .col_idx = gkyl_malloc(sizeof(long[num_nonzeros])),
    .vals = gkyl_malloc(sizeof(double[num_nonzeros])),
  };
  struct fem_gyroavg_csr gyro_mat_trans = {
    .num_rows = num_fld * num_nodes,
    .num_nonzeros = num_nonzeros,
    .row_ptr = gkyl_calloc(num_fld * num_nodes + 1, sizeof(long)),
    .col_idx = gkyl_malloc(sizeof(long[num_nonzeros])),
    .vals = gkyl_malloc(sizeof(double[num_nonzeros])),
  };
  long offset = 0;
  for (int s = 0; s < num_fld; s++) {
    for (int n = 0; n < num_mu; n++) {
      const struct fem_gyroavg_csr *mat = &ups[s]->gyro_mat[n];
      long row_offset = (s * num_mu + n) * num_nodes;
      for (long i = 0; i < num_nodes; i++) {
        gyro_mat.row_ptr[row_offset + i + 1] = offset + mat->row_ptr[i + 1];
      }
      for (long k = 0; k < mat->num_nonzeros; k++) {
        gyro_mat.col_idx[offset + k] = s * num_nodes + mat->col_idx[k];
        gyro_mat.vals[offset + k] = mat->vals[k];
        gyro_mat_trans.row_ptr[s * num_nodes + mat->col_idx[k] + 1] += 1;
      }
      offset += mat->num_nonzeros;
    }
  }
  for (long j = 0; j < num_fld * num_nodes; j++) {
    gyro_mat_trans.row_ptr[j + 1] += gyro_mat_trans.row_ptr[j];
  }
  long *num_filled = gkyl_calloc(num_fld * num_nodes, sizeof(long));
  for (int s = 0; s < num_fld; s++) {
    for (int n = 0; n < num_mu; n++) {
      const struct fem_gyroavg_csr *mat = &ups[s]->gyro_mat[n];
      for (long i = 0; i < num_nodes; i++) {
        for (long k = mat->row_ptr[i]; k < mat->row_ptr[i + 1]; k++) {
          long row = s * num_nodes + mat->col_idx[k];
          long entry = gyro_mat_trans.row_ptr[row] + num_filled[row];
          gyro_mat_trans.col_idx[entry] = (s * num_mu + n) * num_nodes + i;
          gyro_mat_trans.vals[entry] = mat->vals[k];
          num_filled[row] += 1;
        }
      }
    }
  }
  gkyl_free(num_filled);

  fem_gyroavg_csr_dev_copy(&batch->load_mat, &load_mat);
  fem_gyroavg_csr_dev_copy(&batch->gyro_mat, &gyro_mat);
  fem_gyroavg_csr_dev_copy(&batch->gyro_mat_trans, &gyro_mat_trans);
  fem_gyroavg_csr_release(&load_mat);
  fem_gyroavg_csr_release(&gyro_mat);
  fem_gyroavg_csr_release(&gyro_mat_trans);

  // Inverse of the 1D mass matrix of each direction. The node index is row-major.
  for (int d = 0; d < FEM_GYROAVG_MAX_DIM; d++) {
    batch->mass_inv[d] = 0;
    batch->node_stride[d] = 1;
    batch->num_nodes_dir[d] = 1;
  }
  for (int d = 0; d < ndim; d++) {
    batch->num_nodes_dir[d] = up->num_nodes_dir[d];
    for (int e = d + 1; e < ndim; e++) {
      batch->node_stride[d] *= up->num_nodes_dir[e];
    }
    size_t mass_inv_size = sizeof(double[up->num_nodes_dir[d] * up->num_nodes_dir[d]]);
    double *mass_inv = gkyl_malloc(mass_inv_size);
    fem_gyroavg_mass_inv_1d(up, d, mass_inv);
    batch->mass_inv[d] = fem_gyroavg_dev_copy(mass_inv, mass_inv_size);
    gkyl_free(mass_inv);
  }

  batch->cell_nodes = fem_gyroavg_dev_copy(up->cell_nodes, sizeof(long[num_cells * num_basis]));
  batch->cell_linidx = fem_gyroavg_dev_copy(cell_linidx, sizeof(long[num_cells]));
  batch->nod2mod = fem_gyroavg_dev_copy(up->nod2mod, sizeof(up->nod2mod));
  gkyl_free(cell_linidx);

  size_t nodal_size = sizeof(double[num_fld * num_mu * num_nodes]);
  batch->nodal_in = gkyl_cu_malloc(nodal_size);
  batch->nodal_buff = gkyl_cu_malloc(nodal_size);
  batch->nodal_acc = gkyl_cu_malloc(nodal_size);
  batch->nodal_mass_buff = gkyl_cu_malloc(nodal_size);

  return batch;
}

void
gkyl_fem_gyroaverage_batch_advance(
  struct gkyl_fem_gyroaverage_batch *batch, const struct gkyl_array *fin, struct gkyl_array *fout
)
{
#ifdef GKYL_HAVE_CUDA
  gkyl_fem_gyroaverage_batch_advance_cu(batch, fin, fout);
#else
  assert(false); // Batches are only available on the GPU.
#endif
}

void
gkyl_fem_gyroaverage_batch_advance_adjoint(
  struct gkyl_fem_gyroaverage_batch *batch, const struct gkyl_array *fin, struct gkyl_array *fout
)
{
#ifdef GKYL_HAVE_CUDA
  gkyl_fem_gyroaverage_batch_advance_adjoint_cu(batch, fin, fout);
#else
  assert(false); // Batches are only available on the GPU.
#endif
}

void
gkyl_fem_gyroaverage_batch_release(struct gkyl_fem_gyroaverage_batch *batch)
{
  fem_gyroavg_csr_dev_release(&batch->load_mat);
  fem_gyroavg_csr_dev_release(&batch->gyro_mat);
  fem_gyroavg_csr_dev_release(&batch->gyro_mat_trans);
  for (int d = 0; d < batch->ndim; d++) {
    gkyl_cu_free(batch->mass_inv[d]);
  }
  gkyl_cu_free(batch->cell_nodes);
  gkyl_cu_free(batch->cell_linidx);
  gkyl_cu_free(batch->nod2mod);
  gkyl_cu_free(batch->nodal_in);
  gkyl_cu_free(batch->nodal_buff);
  gkyl_cu_free(batch->nodal_acc);
  gkyl_cu_free(batch->nodal_mass_buff);
  gkyl_free(batch);
}

struct gkyl_fem_gyroaverage *
gkyl_fem_gyroaverage_inew(const struct gkyl_fem_gyroaverage_inp *inp)
{
  struct gkyl_fem_gyroaverage *up = gkyl_malloc(sizeof(*up));

  up->ndim = inp->grid->ndim;
  up->num_basis = inp->basis->num_basis;
  up->grid = *inp->grid;
  up->basis = *inp->basis;
  up->local = *inp->local;
  up->num_mu = inp->num_mu;
  up->use_gpu = inp->use_gpu;
  up->reflect_bnd = inp->reflect_bnd;

  // Only implemented for piecewise linear bases in 1D and 2D.
  assert(up->ndim == 1 || up->ndim == 2);
  assert(inp->basis->poly_order == 1 && up->num_basis == (1 << up->ndim));
  assert(up->ndim == 1 || (inp->gxy && inp->gyy));

  up->num_nodes = 1;
  for (int d = 0; d < up->ndim; d++) {
    up->is_periodic[d] = inp->is_periodic[d];
    up->num_cells[d] = up->local.upper[d] - up->local.lower[d] + 1;
    up->num_nodes_dir[d] = up->is_periodic[d] ? up->num_cells[d] : up->num_cells[d] + 1;
    up->num_nodes *= up->num_nodes_dir[d];
  }
  up->load_fac = up->grid.cellVolume / pow(2.0, up->ndim);

  // Location of the nodes in a cell, and nodal to modal transformation.
  double node_coords[FEM_GYROAVG_MAX_NODES * FEM_GYROAVG_MAX_DIM];
  up->basis.node_list(node_coords);
  for (int a = 0; a < up->num_basis; a++) {
    for (int d = 0; d < up->ndim; d++) {
      up->node_offset[a][d] = node_coords[a * up->ndim + d] > 0.0 ? 1 : 0;
    }
    double nodal[FEM_GYROAVG_MAX_NODES] = {0.0}, modal[FEM_GYROAVG_MAX_NODES];
    nodal[a] = 1.0;
    up->basis.nodal_to_modal(nodal, modal);
    for (int k = 0; k < up->num_basis; k++) {
      up->nod2mod[k][a] = modal[k];
    }
  }

  // Global index of the nodes of each cell.
  up->cell_nodes = gkyl_malloc(sizeof(long[up->local.volume * up->num_basis]));
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->local);
  long cell = 0;
  while (gkyl_range_iter_next(&iter)) {
    for (int a = 0; a < up->num_basis; a++) {
      int node[FEM_GYROAVG_MAX_DIM] = {0};
      for (int d = 0; d < up->ndim; d++) {
        node[d] = iter.idx[d] - up->local.lower[d] + up->node_offset[a][d];
      }
      up->cell_nodes[cell * up->num_basis + a] = fem_gyroavg_node_idx(up, node);
    }
    cell++;
  }

  fem_gyroavg_mass_new(up);

  up->gyro_mat = gkyl_malloc(sizeof(struct fem_gyroavg_csr[up->num_mu]));
  for (int n = 0; n < up->num_mu; n++) {
    fem_gyroavg_gyro_mat_new(up, inp, n);
  }

  up->nodal_in = gkyl_malloc(sizeof(double[up->num_nodes]));
  up->nodal_buff = gkyl_malloc(sizeof(double[up->num_nodes]));
  up->nodal_acc = gkyl_malloc(sizeof(double[up->num_nodes]));

  // The operator is assembled on the host and copied to the device.
  if (up->use_gpu) {
    up->batch = gkyl_fem_gyroaverage_batch_new(&up, 1);
  }

  return up;
}

void
gkyl_fem_gyroaverage_advance(
  struct gkyl_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_fem_gyroaverage_batch_advance_cu(up->batch, fin, fout);
  }
#endif

  // Nodal values of the projection of the input onto the continuous space.
  fem_gyroavg_dg_to_load(up, fin, 0, up->nodal_buff);
  fem_gyroavg_mass_solve(up, up->nodal_buff, up->nodal_in);

  for (int n = 0; n < up->num_mu; n++) {
    fem_gyroavg_csr_mat_vec(&up->gyro_mat[n], up->nodal_in, up->nodal_buff);
    fem_gyroavg_mass_solve(up, up->nodal_buff, up->nodal_buff);
    fem_gyroavg_nodal_to_dg(up, up->nodal_buff, fout, n * up->num_basis);
  }
}

void
gkyl_fem_gyroaverage_advance_adjoint(
  struct gkyl_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_fem_gyroaverage_batch_advance_adjoint_cu(up->batch, fin, fout);
  }
#endif

  for (long i = 0; i < up->num_nodes; i++) {
    up->nodal_acc[i] = 0.0;
  }

  for (int n = 0; n < up->num_mu; n++) {
    fem_gyroavg_dg_to_load(up, fin, n * up->num_basis, up->nodal_buff);
    fem_gyroavg_mass_solve(up, up->nodal_buff, up->nodal_in);
    fem_gyroavg_csr_mat_trans_vec_accum(&up->gyro_mat[n], up->nodal_in, up->nodal_acc);
  }

  fem_gyroavg_mass_solve(up, up->nodal_acc, up->nodal_buff);
  fem_gyroavg_nodal_to_dg(up, up->nodal_buff, fout, 0);
}

long
gkyl_fem_gyroaverage_num_nonzeros(const struct gkyl_fem_gyroaverage *up, int mu_idx)
{
  return up->gyro_mat[mu_idx].num_nonzeros;
}

void
gkyl_fem_gyroaverage_release(struct gkyl_fem_gyroaverage *up)
{
  if (up->use_gpu) {
    gkyl_fem_gyroaverage_batch_release(up->batch);
  }
  for (int n = 0; n < up->num_mu; n++) {
    fem_gyroavg_csr_release(&up->gyro_mat[n]);
  }
  gkyl_free(up->gyro_mat);
  gkyl_superlu_prob_release(up->mass_prob);
  gkyl_free(up->cell_nodes);
  gkyl_free(up->nodal_in);
  gkyl_free(up->nodal_buff);
  gkyl_free(up->nodal_acc);
  gkyl_free(up);
}
