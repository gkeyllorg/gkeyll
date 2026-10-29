#include <float.h>
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_orbit_average_gyrokinetic_priv.h>

// Frozen setup uses existing basis evaluators, rather than new generated DG
// kernels. Geometry construction runs on the host; the distribution never does
// on a GPU update.
struct gk_oa_sample {
  double energy, weight, jacobian;
  int group;
};

struct gk_oa_group {
  int mu_node, branch, minimum;
  double lower, upper;
};

struct gk_oa_geometry {
  struct gkyl_array *bmag, *phi, *jacobgeo, *mask, *wall_lo, *wall_up;
};

static struct gkyl_array *
gk_oa_host_copy(const struct gkyl_array *source)
{
  if (!source) {
    return 0;
  }
  struct gkyl_array *copy = gkyl_array_new(source->type, source->ncomp, source->size);
  gkyl_array_copy(copy, source);
  return copy;
}

static void
gk_oa_geometry_release(struct gk_oa_geometry *geo)
{
  gkyl_array_release(geo->bmag);
  gkyl_array_release(geo->phi);
  gkyl_array_release(geo->jacobgeo);
  gkyl_array_release(geo->mask);
  if (geo->wall_lo) {
    gkyl_array_release(geo->wall_lo);
  }
  if (geo->wall_up) {
    gkyl_array_release(geo->wall_up);
  }
}

static double
gk_oa_field(
  const struct gkyl_orbit_average_gyrokinetic_inp *inp, const struct gkyl_array *field, int zindex,
  double eta
)
{
  const double *modal = gkyl_array_cfetch(field, gkyl_range_idx(inp->conf_range, &zindex));
  return inp->conf_basis->eval_expand(&eta, modal);
}

static void
gk_oa_velocity(
  const struct gkyl_orbit_average_gyrokinetic_inp *inp, int vindex, int muindex, double eta_v,
  double eta_mu, double *velocity, double *mu, double *jacobian
)
{
  int index[2] = {vindex, muindex};
  const struct gkyl_velocity_map *map = inp->vel_map;
  const struct gkyl_basis *basis = &map->vmap_basis_ho;
  const double *modal = gkyl_array_cfetch(map->vmap_ho, gkyl_range_idx(&map->local_ext_vel, index));
  *velocity = basis->eval_expand(&eta_v, modal);
  *mu = basis->eval_expand(&eta_mu, modal + basis->num_basis);
  double mu_roundoff =
    64.0 * DBL_EPSILON * (fabs(modal[basis->num_basis]) + fabs(modal[basis->num_basis + 1]));
  if (*mu < 0.0 && *mu >= -mu_roundoff) {
    *mu = 0.0;
  }
  double vprime = basis->eval_grad_expand(0, &eta_v, modal) * 2.0 / inp->grid->dx[1];
  double muprime =
    basis->eval_grad_expand(0, &eta_mu, modal + basis->num_basis) * 2.0 / inp->grid->dx[2];
  *jacobian = vprime > 0.0 && muprime > 0.0 ? vprime * muprime : -1.0;
}

/** Lower-triangular Cholesky, row-major storage. */
static bool
gk_oa_cholesky(int size, double *matrix)
{
  for (int row = 0; row < size; ++row) {
    for (int col = 0; col <= row; ++col) {
      double value = matrix[row * size + col];
      for (int inner = 0; inner < col; ++inner) {
        value -= matrix[row * size + inner] * matrix[col * size + inner];
      }
      if (row == col) {
        if (!(value > 0.0)) {
          return false;
        }
        matrix[row * size + col] = sqrt(value);
      } else {
        matrix[row * size + col] = value / matrix[col * size + col];
      }
    }
  }
  return true;
}

static void
gk_oa_solve(int size, const double *factor, double *vector)
{
  for (int row = 0; row < size; ++row) {
    for (int col = 0; col < row; ++col) {
      vector[row] -= factor[row * size + col] * vector[col];
    }
    vector[row] /= factor[row * size + row];
  }
  for (int row = size - 1; row >= 0; --row) {
    for (int col = row + 1; col < size; ++col) {
      vector[row] -= factor[col * size + row] * vector[col];
    }
    vector[row] /= factor[row * size + row];
  }
}

static double
gk_oa_dot(long size, const double *left, const double *right)
{
  double value = 0.0;
  for (long index = 0; index < size; ++index) {
    value += left[index] * right[index];
  }
  return value;
}

/** Build Q and M Q by twice-reorthogonalized modified Gram-Schmidt. */
static void
gk_oa_append_column(
  int size, int *rank, double **orthogonal, double **dual, double *column, double *moment
)
{
  double original = gk_oa_dot(size, column, moment);
  if (!(original > 0.0)) {
    return;
  }
  double scale = 1.0 / sqrt(original);
  for (int index = 0; index < size; ++index) {
    column[index] *= scale;
    moment[index] *= scale;
  }
  for (int pass = 0; pass < 2; ++pass) {
    for (int mode = 0; mode < *rank; ++mode) {
      const double *qmode = *orthogonal + (long)mode * size;
      const double *dmode = *dual + (long)mode * size;
      double coefficient = gk_oa_dot(size, qmode, moment);
      for (int index = 0; index < size; ++index) {
        column[index] -= coefficient * qmode[index];
        moment[index] -= coefficient * dmode[index];
      }
    }
  }
  double residual = gk_oa_dot(size, column, moment);
  // Only remove numerically dependent columns; do not regularize the Gram
  // matrix, since regularization destroys the invariant moment identities.
  if (!(residual > 1e-24) || *rank >= size) {
    return;
  }
  *orthogonal = gkyl_realloc(*orthogonal, (long)(*rank + 1) * size * sizeof(double));
  *dual = gkyl_realloc(*dual, (long)(*rank + 1) * size * sizeof(double));
  scale = 1.0 / sqrt(residual);
  for (int index = 0; index < size; ++index) {
    (*orthogonal)[(long)*rank * size + index] = column[index] * scale;
    (*dual)[(long)*rank * size + index] = moment[index] * scale;
  }
  ++*rank;
}

/** Fill both DG traces of U on the complete parallel domain. */
static void
gk_oa_potential(
  const struct gkyl_orbit_average_gyrokinetic_inp *inp, const struct gk_oa_geometry *geo, double mu,
  double *potential
)
{
  for (int iz = inp->conf_range->lower[0]; iz <= inp->conf_range->upper[0]; ++iz) {
    int local = iz - inp->conf_range->lower[0];
    for (int side = 0; side < 2; ++side) {
      double eta = side ? 1.0 : -1.0;
      potential[2 * local + side] = mu * gk_oa_field(inp, geo->bmag, iz, eta) +
                                    inp->charge * gk_oa_field(inp, geo->phi, iz, eta);
    }
  }
}

/** Identify an accessible component and verify its physical closure. */
static bool
gk_oa_component(
  const struct gkyl_orbit_average_gyrokinetic_inp *inp, const struct gk_oa_geometry *geo,
  const double *potential, int zindex, double mu, double energy, int *minimum
)
{
  int nz = inp->conf_range->upper[0] - inp->conf_range->lower[0] + 1;
  int anchor = zindex - inp->conf_range->lower[0];
  double left_barrier = potential[2 * anchor], right_barrier = potential[2 * anchor + 1];
  for (int iz = 0; iz < anchor; ++iz) {
    left_barrier = fmax(left_barrier, fmax(potential[2 * iz], potential[2 * iz + 1]));
  }
  for (int iz = anchor + 1; iz < nz; ++iz) {
    right_barrier = fmax(right_barrier, fmax(potential[2 * iz], potential[2 * iz + 1]));
  }
  if (inp->lower_orbit == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL) {
    left_barrier = DBL_MAX;
  } else if (inp->lower_orbit == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH) {
    int edge = inp->conf_range->lower[0];
    left_barrier = fmax(
      left_barrier, mu * gk_oa_field(inp, geo->bmag, edge, -1.0) +
                      inp->charge * gk_oa_field(inp, geo->wall_lo, edge, -1.0)
    );
  }
  if (inp->upper_orbit == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL) {
    right_barrier = DBL_MAX;
  } else if (inp->upper_orbit == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH) {
    int edge = inp->conf_range->upper[0];
    right_barrier = fmax(
      right_barrier, mu * gk_oa_field(inp, geo->bmag, edge, 1.0) +
                       inp->charge * gk_oa_field(inp, geo->wall_up, edge, 1.0)
    );
  }
  double scale = fabs(energy);
  for (int iz = 0; iz < 2 * nz; ++iz) {
    scale = fmax(scale, fabs(potential[iz]));
  }
  double tolerance = 64.0 * DBL_EPSILON * fmax(scale, DBL_MIN);
  double margin = fmax(tolerance, inp->separatrix_margin);
  if (!(energy < left_barrier - margin && energy < right_barrier - margin)) {
    return false;
  }

  // Exclude the unresolved layer around every internal saddle. Ordinary
  // turning points and minima are not separatrices and are not removed.
  for (int node = 1; node < 2 * nz - 1; ++node) {
    int left = node - 1, right = node + 1;
    while (left > 0 && fabs(potential[left] - potential[node]) <= tolerance) {
      --left;
    }
    while (right < 2 * nz - 1 && fabs(potential[right] - potential[node]) <= tolerance) {
      ++right;
    }
    if (potential[node] > potential[left] + tolerance &&
        potential[node] > potential[right] + tolerance &&
        fabs(energy - potential[node]) <= margin) {
      return false;
    }
  }

  // Scan the connected allowed interval. Both traces at a DG interface must
  // be accessible before it is crossed. The lowest endpoint gives a stable
  // component label between critical energies; ties use the leftmost index.
  int lower = 2 * anchor, upper = lower + 1;
  while (lower > 0 && potential[lower] <= energy + tolerance &&
         potential[lower - 1] <= energy + tolerance) {
    --lower;
  }
  while (upper < 2 * nz - 1 && potential[upper] <= energy + tolerance &&
         potential[upper + 1] <= energy + tolerance) {
    ++upper;
  }
  *minimum = lower;
  for (int node = lower; node <= upper; ++node) {
    if (potential[node] < potential[*minimum] - tolerance) {
      *minimum = node;
    }
  }
  double kinetic = fmax(0.0, energy - potential[*minimum]);
  double speed = sqrt(2.0 * kinetic / inp->mass);
  // The whole bounce, including both signs, must fit inside velocity space.
  return -speed > inp->vel_map->vbounds[0] && speed < inp->vel_map->vbounds[2];
}

static int
gk_oa_compare_double(const void *left, const void *right)
{
  double lhs = *(const double *)left, rhs = *(const double *)right;
  return (lhs > rhs) - (lhs < rhs);
}

/** Critical levels include discontinuous extrema as separate barriers. */
static int
gk_oa_critical(int nz, const double *potential, double *levels)
{
  int count = 0;
  levels[count++] = potential[0];
  for (int node = 1; node < 2 * nz - 1; ++node) {
    // Skip a plateau to find the neighboring unequal values. This avoids
    // treating every continuous P1 interface's duplicate value as critical.
    int left = node - 1, right = node + 1;
    double tol = 64.0 * DBL_EPSILON * fmax(fabs(potential[node]), DBL_MIN);
    while (left > 0 && fabs(potential[left] - potential[node]) <= tol) {
      --left;
    }
    while (right < 2 * nz - 1 && fabs(potential[right] - potential[node]) <= tol) {
      ++right;
    }
    double slope_left = potential[node] - potential[left];
    double slope_right = potential[right] - potential[node];
    if ((slope_left > tol && slope_right < -tol) || (slope_left < -tol && slope_right > tol)) {
      levels[count++] = potential[node];
    }
  }
  levels[count++] = potential[2 * nz - 1];
  qsort(levels, count, sizeof(double), gk_oa_compare_double);
  int unique = 0;
  for (int node = 0; node < count; ++node) {
    double tol = 64.0 * DBL_EPSILON * fmax(fabs(levels[node]), DBL_MIN);
    if (unique == 0 || fabs(levels[node] - levels[unique - 1]) > tol) {
      levels[unique++] = levels[node];
    }
  }
  return unique;
}

/** Validate every monitored point, including mu interiors missed by corners. */
static bool
gk_oa_accept_cell(
  const struct gkyl_orbit_average_gyrokinetic_inp *inp, const struct gk_oa_geometry *geo,
  const int *index, const double *nodes, double *potential, double *levels
)
{
  int cell_branch[5] = {-1, -1, -1, -1, -1};
  int cell_minimum[5] = {-1, -1, -1, -1, -1};
  int nz = inp->conf_range->upper[0] - inp->conf_range->lower[0] + 1;
  // Also test the mu midpoint at (z,v) corners. Competing barriers can open
  // between both mu corners and even between the two Gauss nodes. This is an
  // additional sampled safeguard, not an exact interior-envelope proof.
  for (int node = 0; node < GKYL_ORBIT_NUM_POS + 4; ++node) {
    double midpoint[3] = {
      (node - GKYL_ORBIT_NUM_POS) / 2 ? 1.0 : -1.0, (node - GKYL_ORBIT_NUM_POS) % 2 ? 1.0 : -1.0,
      0.0
    };
    const double *eta = node < GKYL_ORBIT_NUM_POS ? nodes + 3 * node : midpoint;
    int mu_slot = node < GKYL_ORBIT_NUM_QUAD ? node % 2 :
                  node < GKYL_ORBIT_NUM_POS  ? 2 + (node - GKYL_ORBIT_NUM_QUAD) % 2 :
                                               4;
    double velocity, mu, jacvel;
    gk_oa_velocity(inp, index[1], index[2], eta[1], eta[2], &velocity, &mu, &jacvel);
    double bmag = gk_oa_field(inp, geo->bmag, index[0], eta[0]);
    double jacgeo = gk_oa_field(inp, geo->jacobgeo, index[0], eta[0]);
    double potential_local = mu * bmag + inp->charge * gk_oa_field(inp, geo->phi, index[0], eta[0]);
    double energy = 0.5 * inp->mass * velocity * velocity + potential_local;
    if (!(bmag > 0.0 && jacgeo > 0.0 && jacvel > 0.0 && mu >= 0.0)) {
      return false;
    }
    gk_oa_potential(inp, geo, mu, potential);
    int minimum;
    if (!gk_oa_component(inp, geo, potential, index[0], mu, energy, &minimum)) {
      return false;
    }
    int nlevels = gk_oa_critical(nz, potential, levels), branch = 0;
    while (branch < nlevels && energy >= levels[branch]) {
      ++branch;
    }
    // A whole-cell DG update cannot isolate different orbit families inside
    // one polynomial. Leave cells crossing a topology change unchanged.
    if (cell_branch[mu_slot] == -1) {
      cell_branch[mu_slot] = branch;
      cell_minimum[mu_slot] = minimum;
    } else if (cell_branch[mu_slot] != branch || cell_minimum[mu_slot] != minimum) {
      return false;
    }
  }
  return true;
}

/** Construct one whole-mu-cell DG projector. */
static bool
gk_oa_build_block(
  gkyl_orbit_average_gyrokinetic *up, const struct gkyl_orbit_average_gyrokinetic_inp *inp,
  const struct gk_oa_geometry *geo, int muindex, const double *nodes, const double *weights
)
{
  int nz = inp->conf_range->upper[0] - inp->conf_range->lower[0] + 1;
  long capacity = (long)nz * (inp->phase_range->upper[1] - inp->phase_range->lower[1] + 1);
  long *indices = gkyl_malloc(capacity * sizeof(long));
  int *zindices = gkyl_malloc(capacity * sizeof(int));
  int *vindices = gkyl_malloc(capacity * sizeof(int));
  double *potential = gkyl_malloc(2 * nz * sizeof(double));
  double *levels = gkyl_malloc(2 * nz * sizeof(double));
  int count = 0;
  int index[3] = {0, 0, muindex};
  for (index[1] = inp->phase_range->lower[1]; index[1] <= inp->phase_range->upper[1]; ++index[1]) {
    for (index[0] = inp->phase_range->lower[0]; index[0] <= inp->phase_range->upper[0];
         ++index[0]) {
      long location = gkyl_range_idx(inp->phase_range, index);
      const double *candidate = gkyl_array_cfetch(geo->mask, location);
      if (!(candidate[0] > 0.5)) {
        continue;
      }
      if (!gk_oa_accept_cell(inp, geo, index, nodes, potential, levels)) {
        ++up->num_rejected;
        continue;
      }
      indices[count] = location;
      zindices[count] = index[0];
      vindices[count++] = index[1];
    }
  }
  if (!count) {
    gkyl_free(indices);
    gkyl_free(zindices);
    gkyl_free(vindices);
    gkyl_free(potential);
    gkyl_free(levels);
    return true;
  }

  int nb = up->num_basis, nm = count * nb;
  long ns = (long)count * GKYL_ORBIT_NUM_QUAD;
  struct gk_oa_sample *samples = gkyl_malloc(ns * sizeof(*samples));
  struct gk_oa_group *groups = gkyl_malloc(ns * sizeof(*groups));
  double *matrices = gkyl_calloc((long)count * nb * nb, sizeof(double));
  double volume = inp->grid->dx[0] * inp->grid->dx[1] * inp->grid->dx[2] / 8.0;
  int ngroups = 0;
  double energy_min = DBL_MAX, energy_max = -DBL_MAX;
  for (int cell = 0; cell < count; ++cell) {
    for (int node = 0; node < GKYL_ORBIT_NUM_QUAD; ++node) {
      const double *eta = nodes + 3 * node;
      double velocity, mu, jacvel;
      gk_oa_velocity(inp, vindices[cell], muindex, eta[1], eta[2], &velocity, &mu, &jacvel);
      double bmag = gk_oa_field(inp, geo->bmag, zindices[cell], eta[0]);
      double jacgeo = gk_oa_field(inp, geo->jacobgeo, zindices[cell], eta[0]);
      double energy = 0.5 * inp->mass * velocity * velocity + mu * bmag +
                      inp->charge * gk_oa_field(inp, geo->phi, zindices[cell], eta[0]);
      gk_oa_potential(inp, geo, mu, potential);
      int nlevels = gk_oa_critical(nz, potential, levels), minimum;
      gk_oa_component(inp, geo, potential, zindices[cell], mu, energy, &minimum);
      int branch = 0;
      while (branch < nlevels && energy >= levels[branch]) {
        ++branch;
      }
      int mu_node = node % 2;
      int group = 0;
      while (group < ngroups &&
             !(groups[group].mu_node == mu_node && groups[group].branch == branch &&
               groups[group].minimum == minimum)) {
        ++group;
      }
      if (group == ngroups) {
        groups[group] = (struct gk_oa_group){mu_node, branch, minimum, energy, energy};
        ++ngroups;
      }
      groups[group].lower = fmin(groups[group].lower, energy);
      groups[group].upper = fmax(groups[group].upper, energy);
      energy_min = fmin(energy_min, energy);
      energy_max = fmax(energy_max, energy);
      struct gk_oa_sample *sample = samples + (long)cell * GKYL_ORBIT_NUM_QUAD + node;
      *sample =
        (struct gk_oa_sample){energy, volume * weights[node], bmag * jacgeo * jacvel, group};
      const double *basis = up->eval + node * nb;
      double weight = sample->weight / sample->jacobian;
      double *matrix = matrices + (long)cell * nb * nb;
      for (int row = 0; row < nb; ++row) {
        for (int col = 0; col < nb; ++col) {
          matrix[row * nb + col] += weight * basis[row] * basis[col];
        }
      }
    }
  }
  bool success = true;
  for (int cell = 0; cell < count; ++cell) {
    if (!gk_oa_cholesky(nb, matrices + (long)cell * nb * nb)) {
      success = false;
      break;
    }
  }
  double *orthogonal = 0, *dual = 0;
  double *column = gkyl_malloc(nm * sizeof(double));
  double *moment = gkyl_malloc(nm * sizeof(double));
  int rank = 0;
  if (success) {
    double base_spacing = (energy_max - energy_min) / (inp->num_energy_nodes - 1);
    for (int group = 0; group < ngroups; ++group) {
      double width = groups[group].upper - groups[group].lower;
      int intervals = base_spacing > 0.0 ? (int)ceil(width / base_spacing) : 1;
      intervals = intervals > 1 ? intervals : 1;
      // Constant and normalized energy first protect invariant moments before
      // any dependent hat columns are removed. Interior hats complete the
      // same piecewise-linear space as endpoint-plus-interior hat functions.
      int columns = width > 0.0 ? intervals + 1 : 1;
      for (int mode = 0; mode < columns; ++mode) {
        memset(moment, 0, nm * sizeof(double));
        for (int cell = 0; cell < count; ++cell) {
          for (int node = 0; node < GKYL_ORBIT_NUM_QUAD; ++node) {
            const struct gk_oa_sample *sample = samples + (long)cell * GKYL_ORBIT_NUM_QUAD + node;
            if (sample->group != group) {
              continue;
            }
            double coordinate = width > 0.0 ? (sample->energy - groups[group].lower) / width : 0.0;
            double value = 1.0;
            if (mode == 1) {
              value = coordinate;
            } else if (mode > 1) {
              value = fmax(0.0, 1.0 - fabs(coordinate * intervals - (mode - 1)));
            }
            const double *basis = up->eval + node * nb;
            for (int row = 0; row < nb; ++row) {
              moment[cell * nb + row] += sample->weight * value * basis[row];
            }
          }
        }
        memcpy(column, moment, nm * sizeof(double));
        for (int cell = 0; cell < count; ++cell) {
          gk_oa_solve(nb, matrices + (long)cell * nb * nb, column + cell * nb);
        }
        gk_oa_append_column(nm, &rank, &orthogonal, &dual, column, moment);
      }
    }
    success = rank > 0;
  }
  if (success) {
    int block = up->num_blocks++;
    up->blocks = gkyl_realloc(up->blocks, up->num_blocks * sizeof(*up->blocks));
    up->blocks[block] = (struct gkyl_orbit_average_gyrokinetic_block){
      up->num_cells, count, nm, rank, up->basis_size, up->num_reduced,
    };
    up->cell_indices = gkyl_realloc(up->cell_indices, (up->num_cells + count) * sizeof(long));
    memcpy(up->cell_indices + up->num_cells, indices, count * sizeof(long));
    long entries = (long)rank * nm;
    up->q = gkyl_realloc(up->q, (up->basis_size + entries) * sizeof(double));
    up->dual = gkyl_realloc(up->dual, (up->basis_size + entries) * sizeof(double));
    memcpy(up->q + up->basis_size, orthogonal, entries * sizeof(double));
    memcpy(up->dual + up->basis_size, dual, entries * sizeof(double));
    up->basis_size += entries;
    up->num_cells += count;
    up->num_reduced += rank;
  }
  gkyl_free(indices);
  gkyl_free(zindices);
  gkyl_free(vindices);
  gkyl_free(potential);
  gkyl_free(samples);
  gkyl_free(groups);
  gkyl_free(matrices);
  gkyl_free(levels);
  gkyl_free(orthogonal);
  gkyl_free(dual);
  gkyl_free(column);
  gkyl_free(moment);
  return success;
}

// Restrict rows after constructing a common global basis. Keep every block and
// reduced column, even on ranks with no selected cells, for matching collectives.
static void
gk_oa_localize(
  gkyl_orbit_average_gyrokinetic *up, const struct gkyl_orbit_average_gyrokinetic_inp *inp
)
{
  up->num_global_cells = up->num_cells;
  up->reduce = inp->reduce;
  up->reduce_ctx = inp->reduce_ctx;
  if (!inp->local_range) {
    return;
  }
  long cells = 0, entries = 0;
  for (int ib = 0; ib < up->num_blocks; ++ib) {
    struct gkyl_orbit_average_gyrokinetic_block old = up->blocks[ib];
    int *selected = gkyl_malloc((old.num_cells > 0 ? old.num_cells : 1) * sizeof(int));
    int count = 0;
    for (int cell = 0; cell < old.num_cells; ++cell) {
      int idx[GKYL_MAX_DIM];
      gkyl_range_inv_idx(inp->phase_range, up->cell_indices[old.cell_offset + cell], idx);
      if (gkyl_range_contains_idx(inp->local_range, idx)) {
        selected[count] = cell;
        up->cell_indices[cells + count++] = gkyl_range_idx(inp->local_range, idx);
      }
    }
    int modes = count * up->num_basis;
    for (int col = 0; col < old.rank; ++col) {
      for (int cell = 0; cell < count; ++cell) {
        for (int mode = 0; mode < up->num_basis; ++mode) {
          long src =
            old.basis_offset + (long)col * old.num_modes + selected[cell] * up->num_basis + mode;
          long dst = entries + (long)col * modes + cell * up->num_basis + mode;
          up->q[dst] = up->q[src];
          up->dual[dst] = up->dual[src];
        }
      }
    }
    up->blocks[ib].cell_offset = cells;
    up->blocks[ib].num_cells = count;
    up->blocks[ib].num_modes = modes;
    up->blocks[ib].basis_offset = entries;
    cells += count;
    entries += (long)modes * old.rank;
    gkyl_free(selected);
  }
  up->num_cells = cells;
  up->basis_size = entries;
  up->cell_indices = gkyl_realloc(up->cell_indices, (cells > 0 ? cells : 1) * sizeof(long));
  up->q = gkyl_realloc(up->q, (entries > 0 ? entries : 1) * sizeof(double));
  up->dual = gkyl_realloc(up->dual, (entries > 0 ? entries : 1) * sizeof(double));
}

void
gkyl_orbit_average_gyrokinetic_release(gkyl_orbit_average_gyrokinetic *up)
{
  if (!up) {
    return;
  }
#ifdef GKYL_HAVE_CUDA
  if (up->gpu) {
    gkyl_orbit_average_gyrokinetic_cu_release(up);
  }
#endif
  gkyl_free(up->blocks);
  gkyl_free(up->cell_indices);
  gkyl_free(up->q);
  gkyl_free(up->dual);
  gkyl_free(up->eval);
  gkyl_free(up->coeff);
  gkyl_free(up->correction);
  gkyl_free(up->blend);
  gkyl_free(up->reduce_buffer);
  gkyl_free(up);
}

/** Check array indexing without assuming the range is a standalone range. */
static bool
gk_oa_array_valid(const struct gkyl_array *array, const struct gkyl_range *range, int components)
{
  if (!array || array->type != GKYL_DOUBLE || array->ncomp != components || range->volume < 1) {
    return false;
  }
  long lower = gkyl_range_idx(range, range->lower), upper = gkyl_range_idx(range, range->upper);
  return lower >= 0 && upper >= lower && upper < array->size;
}

/** Invalid geometry is a setup error, not merely an ineligible orbit. */
static bool
gk_oa_geometry_valid(
  const struct gkyl_orbit_average_gyrokinetic_inp *inp, const struct gk_oa_geometry *geo
)
{
  for (int iz = inp->conf_range->lower[0]; iz <= inp->conf_range->upper[0]; ++iz) {
    for (int side = 0; side < 2; ++side) {
      double eta = side ? 1.0 : -1.0;
      double bmag = gk_oa_field(inp, geo->bmag, iz, eta);
      double jacobian = gk_oa_field(inp, geo->jacobgeo, iz, eta);
      if (!(bmag > 0.0 && jacobian > 0.0)) {
        return false;
      }
    }
  }
  for (int iv = inp->phase_range->lower[1]; iv <= inp->phase_range->upper[1]; ++iv) {
    for (int imu = inp->phase_range->lower[2]; imu <= inp->phase_range->upper[2]; ++imu) {
      for (int side = 0; side < 2; ++side) {
        double eta = side ? 1.0 : -1.0;
        double velocity, mu, jacobian;
        gk_oa_velocity(inp, iv, imu, eta, eta, &velocity, &mu, &jacobian);
        if (!(mu >= 0.0 && jacobian > 0.0)) {
          return false;
        }
      }
    }
  }
  return true;
}

gkyl_orbit_average_gyrokinetic *
gkyl_orbit_average_gyrokinetic_inew(
  const struct gkyl_orbit_average_gyrokinetic_inp *inp,
  enum gkyl_orbit_average_gyrokinetic_status *status
)
{
  enum gkyl_orbit_average_gyrokinetic_status reason = GKYL_ORBIT_AVERAGE_INVALID_INPUT;
  if (status) {
    *status = reason;
  }
  if (!inp || !inp->grid || !inp->conf_basis || !inp->phase_basis || !inp->phase_range ||
      !inp->conf_range || !inp->vel_map || !inp->bmag || !inp->phi || !inp->jacobgeo ||
      !inp->mask || !(inp->mass > 0.0) || !(inp->separatrix_margin >= 0.0) ||
      inp->num_energy_nodes < 2) {
    return 0;
  }
  if (inp->local_range) {
    if (!inp->reduce || inp->local_range->ndim != inp->phase_range->ndim) {
      return 0;
    }
    for (int d = 0; d < inp->phase_range->ndim; ++d) {
      if (inp->local_range->lower[d] < inp->phase_range->lower[d] ||
          inp->local_range->upper[d] > inp->phase_range->upper[d]) {
        return 0;
      }
    }
  }
  if (inp->conf_basis->ndim != 1 || inp->conf_basis->poly_order != 1 ||
      inp->phase_basis->ndim != 3 || inp->phase_basis->b_type != GKYL_BASIS_MODAL_GKHYBRID ||
      inp->phase_basis->num_basis != GKYL_ORBIT_NUM_BASIS ||
      inp->vel_map->vmap_basis_ho.poly_order != 1) {
    if (status) {
      *status = GKYL_ORBIT_AVERAGE_UNSUPPORTED;
    }
    return 0;
  }
#ifndef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    if (status) {
      *status = GKYL_ORBIT_AVERAGE_UNSUPPORTED;
    }
    return 0;
  }
#endif
  if (inp->grid->ndim != 3 || inp->phase_range->ndim != 3 || inp->conf_range->ndim != 1 ||
      inp->phase_range->lower[0] != inp->conf_range->lower[0] ||
      inp->phase_range->upper[0] != inp->conf_range->upper[0] ||
      inp->lower_orbit < GKYL_GK_TRAP_PASS_ORBIT_PASSING ||
      inp->lower_orbit > GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL ||
      inp->upper_orbit < GKYL_GK_TRAP_PASS_ORBIT_PASSING ||
      inp->upper_orbit > GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL ||
      (inp->lower_orbit == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH && !inp->phi_wall_lo) ||
      (inp->upper_orbit == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH && !inp->phi_wall_up)) {
    return 0;
  }
  if (!gk_oa_array_valid(inp->bmag, inp->conf_range, inp->conf_basis->num_basis) ||
      !gk_oa_array_valid(inp->phi, inp->conf_range, inp->conf_basis->num_basis) ||
      !gk_oa_array_valid(inp->jacobgeo, inp->conf_range, inp->conf_basis->num_basis) ||
      !gk_oa_array_valid(inp->mask, inp->phase_range, 1) ||
      (inp->phi_wall_lo &&
       !gk_oa_array_valid(inp->phi_wall_lo, inp->conf_range, inp->conf_basis->num_basis)) ||
      (inp->phi_wall_up &&
       !gk_oa_array_valid(inp->phi_wall_up, inp->conf_range, inp->conf_basis->num_basis))) {
    return 0;
  }
  for (int dimension = 0; dimension < 3; ++dimension) {
    if (!(inp->grid->dx[dimension] > 0.0)) {
      return 0;
    }
  }
  if (inp->conf_range->volume != inp->grid->cells[0]) {
    return 0;
  }
  for (int dimension = 0; dimension < 2; ++dimension) {
    if (inp->phase_range->lower[dimension + 1] != inp->vel_map->local_vel.lower[dimension] ||
        inp->phase_range->upper[dimension + 1] != inp->vel_map->local_vel.upper[dimension]) {
      return 0;
    }
  }

  gkyl_orbit_average_gyrokinetic *up = gkyl_calloc(1, sizeof(*up));
  up->num_basis = GKYL_ORBIT_NUM_BASIS;
  up->num_pos = GKYL_ORBIT_NUM_POS;
  up->use_gpu = inp->use_gpu;
  up->positivity = inp->positivity;
  up->frequency_is_infinite = inp->frequency_is_infinite;
  up->eval = gkyl_malloc(up->num_pos * up->num_basis * sizeof(double));
  double nodes[GKYL_ORBIT_NUM_POS * 3], weights[GKYL_ORBIT_NUM_QUAD];
  int node = 0;
  for (int iz = 0; iz < 4; ++iz) {
    for (int iv = 0; iv < 4; ++iv) {
      for (int im = 0; im < 2; ++im, ++node) {
        nodes[3 * node] = gkyl_gauss_ordinates[4][iz];
        nodes[3 * node + 1] = gkyl_gauss_ordinates[4][iv];
        nodes[3 * node + 2] = gkyl_gauss_ordinates[2][im];
        weights[node] =
          gkyl_gauss_weights[4][iz] * gkyl_gauss_weights[4][iv] * gkyl_gauss_weights[2][im];
      }
    }
  }
  for (int iz = 0; iz < 2; ++iz) {
    for (int iv = 0; iv < 2; ++iv) {
      for (int im = 0; im < 2; ++im, ++node) {
        nodes[3 * node] = iz ? 1.0 : -1.0;
        nodes[3 * node + 1] = iv ? 1.0 : -1.0;
        nodes[3 * node + 2] = im ? 1.0 : -1.0;
      }
    }
  }
  for (node = 0; node < up->num_pos; ++node) {
    inp->phase_basis->eval(nodes + 3 * node, up->eval + node * up->num_basis);
  }
  struct gk_oa_geometry geo = {gk_oa_host_copy(inp->bmag),
                               gk_oa_host_copy(inp->phi),
                               gk_oa_host_copy(inp->jacobgeo),
                               gk_oa_host_copy(inp->mask),
                               gk_oa_host_copy(inp->phi_wall_lo),
                               gk_oa_host_copy(inp->phi_wall_up)};
  if (!gk_oa_geometry_valid(inp, &geo)) {
    gk_oa_geometry_release(&geo);
    gkyl_orbit_average_gyrokinetic_release(up);
    return 0;
  }
  bool success = true;
  for (int imu = inp->phase_range->lower[2]; imu <= inp->phase_range->upper[2]; ++imu) {
    if (!gk_oa_build_block(up, inp, &geo, imu, nodes, weights)) {
      success = false;
      break;
    }
  }
  gk_oa_geometry_release(&geo);
  if (!success) {
    gkyl_orbit_average_gyrokinetic_release(up);
    if (status) {
      *status = GKYL_ORBIT_AVERAGE_SINGULAR;
    }
    return 0;
  }
  gk_oa_localize(up, inp);
  up->max_cell_index = -1;
  for (long cell = 0; cell < up->num_cells; ++cell) {
    if (up->cell_indices[cell] > up->max_cell_index) {
      up->max_cell_index = up->cell_indices[cell];
    }
  }
  long reduce_size = up->num_reduced > up->num_blocks ? up->num_reduced : up->num_blocks;
  up->reduce_buffer = gkyl_calloc(reduce_size > 0 ? reduce_size : 1, sizeof(double));
  up->coeff = gkyl_calloc(up->num_reduced > 0 ? up->num_reduced : 1, sizeof(double));
  up->correction =
    gkyl_calloc(up->num_cells > 0 ? up->num_cells * up->num_basis : 1, sizeof(double));
  up->blend = gkyl_calloc(up->num_blocks > 0 ? up->num_blocks : 1, sizeof(double));
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_orbit_average_gyrokinetic_cu_init(up);
  }
#endif
  if (status) {
    *status = GKYL_ORBIT_AVERAGE_SUCCESS;
  }
  return up;
}

struct gkyl_orbit_average_gyrokinetic_info
gkyl_orbit_average_gyrokinetic_get_info(const gkyl_orbit_average_gyrokinetic *up)
{
  return (struct gkyl_orbit_average_gyrokinetic_info){
    up->num_global_cells,
    up->num_rejected,
    up->num_reduced,
    up->num_blocks,
  };
}

struct gkyl_orbit_average_gyrokinetic_result
gkyl_orbit_average_gyrokinetic_advance(
  gkyl_orbit_average_gyrokinetic *up, double dt, double frequency, struct gkyl_array *f
)
{
  struct gkyl_orbit_average_gyrokinetic_result result = {false, 0.0, 0};
  if (!up || !f || f->type != GKYL_DOUBLE || f->ncomp != up->num_basis || !(dt >= 0.0) ||
      (!up->frequency_is_infinite && !(frequency >= 0.0)) ||
      gkyl_array_is_cu_dev(f) != up->use_gpu) {
    return result;
  }
  if (up->num_cells && up->max_cell_index >= f->size) {
    return result;
  }
  if (dt == 0.0 || (!up->frequency_is_infinite && frequency == 0.0) || up->num_global_cells == 0) {
    result.success = true;
    return result;
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_orbit_average_gyrokinetic_advance_cu(up, dt, frequency, f);
  }
#endif
  double beta = 1.0;
  if (!up->frequency_is_infinite && (dt <= 1.0 || frequency <= DBL_MAX / dt)) {
    double product = dt * frequency;
    beta = product >= 1.0 ? 1.0 - 1.0 / (1.0 + product) : product / (1.0 + product);
  }
  result.min_blend = beta;
  // Build every correction and validate every block before changing f.
  for (int iblock = 0; iblock < up->num_blocks; ++iblock) {
    const struct gkyl_orbit_average_gyrokinetic_block *block = up->blocks + iblock;
    for (int mode = 0; mode < block->rank; ++mode) {
      const double *dual = up->dual + block->basis_offset + (long)mode * block->num_modes;
      double coefficient = 0.0;
      for (int cell = 0; cell < block->num_cells; ++cell) {
        const double *modal = gkyl_array_cfetch(f, up->cell_indices[block->cell_offset + cell]);
        for (int basis = 0; basis < up->num_basis; ++basis) {
          coefficient += dual[cell * up->num_basis + basis] * modal[basis];
        }
      }
      up->coeff[block->reduce_offset + mode] = coefficient;
    }
  }
  if (up->reduce) {
    if (up->reduce(up->reduce_ctx, true, up->num_reduced, up->coeff, up->reduce_buffer)) {
      return result;
    }
    memcpy(up->coeff, up->reduce_buffer, up->num_reduced * sizeof(double));
  }
  for (int iblock = 0; iblock < up->num_blocks; ++iblock) {
    const struct gkyl_orbit_average_gyrokinetic_block *block = up->blocks + iblock;
    double bound = DBL_MAX;
    for (int cell = 0; cell < block->num_cells; ++cell) {
      const double *modal = gkyl_array_cfetch(f, up->cell_indices[block->cell_offset + cell]);
      double *correction = up->correction + (block->cell_offset + cell) * up->num_basis;
      for (int basis = 0; basis < up->num_basis; ++basis) {
        double value = -modal[basis];
        for (int mode = 0; mode < block->rank; ++mode) {
          value +=
            up->q[block->basis_offset + (long)mode * block->num_modes + cell * up->num_basis + basis] *
            up->coeff[block->reduce_offset + mode];
        }
        correction[basis] = value;
      }
      if (up->positivity) {
        for (int node = 0; node < up->num_pos; ++node) {
          const double *basis = up->eval + node * up->num_basis;
          double initial = gk_oa_dot(up->num_basis, basis, modal);
          double change = gk_oa_dot(up->num_basis, basis, correction);
          if (!(initial >= 0.0)) {
            bound = -1.0;
          }
          if (change < 0.0) {
            bound = fmin(bound, initial / -change);
          }
        }
      }
    }
    up->blend[iblock] = bound;
  }
  if (up->reduce) {
    if (up->reduce(up->reduce_ctx, false, up->num_blocks, up->blend, up->reduce_buffer)) {
      return result;
    }
    memcpy(up->blend, up->reduce_buffer, up->num_blocks * sizeof(double));
  }
  for (int iblock = 0; iblock < up->num_blocks; ++iblock) {
    double bound = up->blend[iblock];
    if (bound < 0.0) {
      return result;
    }
    up->blend[iblock] = bound < beta ? fmax(0.0, 0.99 * bound) : beta;
    if (up->blend[iblock] < beta) {
      ++result.num_limited_blocks;
    }
    result.min_blend = fmin(result.min_blend, up->blend[iblock]);
  }
  for (int iblock = 0; iblock < up->num_blocks; ++iblock) {
    const struct gkyl_orbit_average_gyrokinetic_block *block = up->blocks + iblock;
    for (int cell = 0; cell < block->num_cells; ++cell) {
      double *modal = gkyl_array_fetch(f, up->cell_indices[block->cell_offset + cell]);
      const double *correction = up->correction + (block->cell_offset + cell) * up->num_basis;
      for (int basis = 0; basis < up->num_basis; ++basis) {
        modal[basis] += up->blend[iblock] * correction[basis];
      }
    }
  }
  result.success = true;
  return result;
}
