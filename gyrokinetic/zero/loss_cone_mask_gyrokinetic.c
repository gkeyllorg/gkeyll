#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_loss_cone_mask_gyrokinetic_priv.h>
#include <gkyl_range.h>

// The loss-cone mask is built from the lower- and upper-z escape barriers.
// A sheath-trapped orbit includes the material-wall potential as an endpoint.

/** Allocate host or device array storage for precomputed corner data. */
static struct gkyl_array*
gk_lcm_mkarr(long nc, long size, bool use_gpu)
{
  return use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
    : gkyl_array_new(GKYL_DOUBLE, nc, size);
}

/**
 * Precompute the configuration-space basis at every tensor-product P1 cell
 * corner. These are interpolation/cell-classification corners, not Gaussian
 * quadrature, positivity, or modal-basis nodes.
 */
static void
gk_lcm_init_corner_values(const struct gkyl_range *corner_range,
  const struct gkyl_basis *basis, struct gkyl_array **basis_at_corners,
  bool use_gpu)
{
  struct gkyl_array *basis_at_corners_ho = gkyl_array_new(GKYL_DOUBLE,
    basis->num_basis, corner_range->volume);
  *basis_at_corners = gk_lcm_mkarr(basis->num_basis, corner_range->volume, use_gpu);

  double eta[GKYL_MAX_DIM] = { 0.0 };
  struct gkyl_range_iter corner_iter;
  gkyl_range_iter_init(&corner_iter, corner_range);
  while (gkyl_range_iter_next(&corner_iter)) {
    for (int d = 0; d < corner_range->ndim; ++d) {
      eta[d] = corner_iter.idx[d] ? 1.0 : -1.0;
    }
    long corner = gkyl_range_idx(corner_range, corner_iter.idx);
    basis->eval(eta, gkyl_array_fetch(basis_at_corners_ho, corner));
  }

  gkyl_array_copy(*basis_at_corners, basis_at_corners_ho);
  gkyl_array_release(basis_at_corners_ho);
}

struct gkyl_loss_cone_mask_gyrokinetic*
gkyl_loss_cone_mask_gyrokinetic_inew(const struct gkyl_loss_cone_mask_gyrokinetic_inp *inp)
{
  gkyl_loss_cone_mask_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  up->vel_map = gkyl_velocity_map_acquire(inp->vel_map);
  up->mass = inp->mass;
  up->charge = inp->charge;
  up->use_gpu = inp->use_gpu;
  up->lower_orbit = inp->lower_orbit;
  up->upper_orbit = inp->upper_orbit;

  assert(up->lower_orbit >= GKYL_GK_TRAP_PASS_ORBIT_PASSING &&
    up->lower_orbit <= GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL);
  assert(up->upper_orbit >= GKYL_GK_TRAP_PASS_ORBIT_PASSING &&
    up->upper_orbit <= GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL);

  up->cdim = inp->conf_basis->ndim;
  up->num_basis_conf = inp->conf_basis->num_basis;
  int corner_shape[GKYL_MAX_CDIM];
  for (int d = 0; d < up->cdim; ++d) {
    corner_shape[d] = 2;
  }
  gkyl_range_init_from_shape(&up->conf_corner_range, up->cdim, corner_shape);

  gk_lcm_init_corner_values(&up->conf_corner_range, inp->conf_basis,
    &up->basis_at_corners_conf, inp->use_gpu);

  return up;
}

void
gkyl_loss_cone_mask_gyrokinetic_advance(gkyl_loss_cone_mask_gyrokinetic *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *bmag, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall_lo, const struct gkyl_array *phi_wall_up,
  struct gkyl_array *mask_out)
{
  assert(up->lower_orbit != GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH || phi_wall_lo);
  assert(up->upper_orbit != GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH || phi_wall_up);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_loss_cone_mask_gyrokinetic_advance_cu(up, phase_range, conf_range, bmag, phi,
      phi_wall_lo, phi_wall_up, mask_out);
    return;
  }
#endif

  int cdim = up->cdim;
  int pdim = phase_range->ndim;
  int vdim = pdim - cdim;
  assert(vdim == 2);
  int num_basis_conf = up->num_basis_conf;

  // Tensor range over the phase-space P1 cell corners. The first cdim
  // indices identify the configuration corner and the remaining indices the
  // velocity corner, without assuming a flattened corner ordering.
  int corner_shape[GKYL_MAX_DIM];
  for (int d = 0; d < pdim; ++d) {
    corner_shape[d] = 2;
  }
  struct gkyl_range phase_corner_range;
  gkyl_range_init_from_shape(&phase_corner_range, pdim, corner_shape);

  // Outer loop over all phase-space cells.
  struct gkyl_range_iter phase_iter;
  gkyl_range_iter_init(&phase_iter, phase_range);
  while (gkyl_range_iter_next(&phase_iter)) {
    int conf_idx[GKYL_MAX_CDIM] = { 0 };
    int vel_idx[GKYL_MAX_VDIM] = { 0 };
    for (int d = 0; d < cdim; ++d) {
      conf_idx[d] = phase_iter.idx[d];
    }
    for (int d = 0; d < vdim; ++d) {
      vel_idx[d] = phase_iter.idx[cdim + d];
    }

    const struct gkyl_velocity_map *gvm = up->vel_map;
    long linidx_vel = gkyl_range_idx(&gvm->local_ext_vel, vel_idx);
    const double *vmap_d = gkyl_array_cfetch(gvm->vmap, linidx_vel);

    // A cell is trapped only if all of its P1 corners are trapped.
    bool cell_trapped = true;
    struct gkyl_range_iter corner_iter;
    gkyl_range_iter_init(&corner_iter, &phase_corner_range);
    while (cell_trapped && gkyl_range_iter_next(&corner_iter)) {
      int conf_corner = gkyl_range_idx(&up->conf_corner_range, corner_iter.idx);

      double vpar_eta[1] = { corner_iter.idx[cdim] ? 1.0 : -1.0 };
      double mu_eta[1] = { corner_iter.idx[cdim + 1] ? 1.0 : -1.0 };
      double vpar = gvm->vmap_basis->eval_expand(vpar_eta, vmap_d);
      double mu = gvm->vmap_basis->eval_expand(mu_eta,
        vmap_d + gvm->vmap_basis->num_basis);

      // Evaluate the Hamiltonian at this phase-space corner.
      long linidx_conf = gkyl_range_idx(conf_range, conf_idx);
      double bmag_curr = gk_lcm_field_corner_value(bmag, up->basis_at_corners_conf,
        num_basis_conf, linidx_conf, conf_corner);
      double phi_curr = gk_lcm_field_corner_value(phi, up->basis_at_corners_conf,
        num_basis_conf, linidx_conf, conf_corner);
      double kinetic_energy = 0.5 * up->mass * vpar * vpar;
      double magnetic_energy = mu * bmag_curr;
      double electric_energy = up->charge * phi_curr;

      // Determine both directional escape barriers at this corner.
      int zdim = cdim - 1;
      double barrier_left, barrier_right;
      gk_lcm_escape_barriers(cdim, num_basis_conf, conf_range,
        &up->conf_corner_range, up->basis_at_corners_conf, phi, bmag,
        phi_wall_lo, phi_wall_up, conf_idx, conf_idx[zdim], conf_corner,
        mu, up->charge, up->lower_orbit, up->upper_orbit,
        &barrier_left, &barrier_right);

      // Equality means that the particle reaches the loss endpoint with zero
      // parallel kinetic energy, so it is not trapped.
      cell_trapped = gk_lcm_hamiltonian_below_barrier(kinetic_energy,
        magnetic_energy, electric_energy, barrier_left)
        && gk_lcm_hamiltonian_below_barrier(kinetic_energy, magnetic_energy,
          electric_energy, barrier_right);
    }

    long linidx_phase = gkyl_range_idx(phase_range, phase_iter.idx);
    double *mask_d = gkyl_array_fetch(mask_out, linidx_phase);
    mask_d[0] = cell_trapped ? 1.0 : 0.0;
  }
}

void
gkyl_loss_cone_mask_gyrokinetic_release(gkyl_loss_cone_mask_gyrokinetic *up)
{
  gkyl_velocity_map_release(up->vel_map);
  gkyl_array_release(up->basis_at_corners_conf);
  gkyl_free(up);
}
