#include <float.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_const.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_loss_cone_mask_gyrokinetic_priv.h>
#include <gkyl_range.h>

// The loss-cone mask is built from the left and right escape barriers. A
// sheath wall is represented by a virtual endpoint with potential phi_wall;
// this is distinct from phi at the plasma/sheath entrance.

static struct gkyl_array *mkarr(long nc, long size, bool use_gpu)
{
  return use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size) :
                   gkyl_array_new(GKYL_DOUBLE, nc, size);
}

static int init_node_values(
  int ndim, const struct gkyl_basis *basis, struct gkyl_array **basis_at_nodes, bool use_gpu
)
{
  int num_nodes = (int)pow(2.0, (double)ndim);
  struct gkyl_array *basis_at_nodes_ho = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, num_nodes);
  *basis_at_nodes = mkarr(basis->num_basis, num_nodes, use_gpu);

  double eta[GKYL_MAX_DIM] = {0.0};
  for (int node = 0; node < num_nodes; ++node) {
    nodal_coords(ndim, node, eta);
    basis->eval(eta, gkyl_array_fetch(basis_at_nodes_ho, node));
  }

  gkyl_array_copy(*basis_at_nodes, basis_at_nodes_ho);
  gkyl_array_release(basis_at_nodes_ho);

  return num_nodes;
}

struct gkyl_loss_cone_mask_gyrokinetic *
gkyl_loss_cone_mask_gyrokinetic_inew(const struct gkyl_loss_cone_mask_gyrokinetic_inp *inp)
{
  gkyl_loss_cone_mask_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  up->vel_map = gkyl_velocity_map_acquire(inp->vel_map);
  up->mass = inp->mass;
  up->charge = inp->charge;
  up->use_gpu = inp->use_gpu;
  up->lower_boundary = inp->lower_boundary;
  up->upper_boundary = inp->upper_boundary;

  assert(
    up->lower_boundary >= GKYL_LOSS_CONE_BC_OPEN && up->lower_boundary <= GKYL_LOSS_CONE_BC_CLOSED
  );
  assert(
    up->upper_boundary >= GKYL_LOSS_CONE_BC_OPEN && up->upper_boundary <= GKYL_LOSS_CONE_BC_CLOSED
  );

  up->cdim = inp->conf_basis->ndim;
  up->num_basis_conf = inp->conf_basis->num_basis;
  up->num_nodes_conf = (int)pow(2.0, (double)up->cdim);

  init_node_values(up->cdim, inp->conf_basis, &up->basis_at_nodes_conf, inp->use_gpu);

  return up;
}

void gkyl_loss_cone_mask_gyrokinetic_advance(
  gkyl_loss_cone_mask_gyrokinetic *up, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range, const struct gkyl_array *bmag, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall_lo, const struct gkyl_array *phi_wall_up,
  struct gkyl_array *mask_out
)
{
  assert(up->lower_boundary != GKYL_LOSS_CONE_BC_SHEATH || phi_wall_lo);
  assert(up->upper_boundary != GKYL_LOSS_CONE_BC_SHEATH || phi_wall_up);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_loss_cone_mask_gyrokinetic_advance_cu(
      up, phase_range, conf_range, bmag, phi, phi_wall_lo, phi_wall_up, mask_out
    );
    return;
  }
#endif

  int cdim = up->cdim;
  int pdim = phase_range->ndim;
  int vdim = pdim - cdim;
  int num_basis_conf = up->num_basis_conf;
  int num_phase_nodes = (int)pow(2.0, (double)pdim);
  int num_vel_nodes = (int)pow(2.0, (double)vdim);

  // Outer loop over all phase space cells
  struct gkyl_range_iter phase_iter;
  gkyl_range_iter_init(&phase_iter, phase_range);
  while (gkyl_range_iter_next(&phase_iter)) {
    int conf_idx[GKYL_MAX_DIM] = {0};
    int vel_idx[GKYL_MAX_DIM] = {0};
    for (int d = 0; d < cdim; ++d) {
      conf_idx[d] = phase_iter.idx[d];
    }
    for (int d = 0; d < vdim; ++d) {
      vel_idx[d] = phase_iter.idx[cdim + d];
    }

    const struct gkyl_velocity_map *gvm = up->vel_map;
    long linidx_vel = gkyl_range_idx(&gvm->local_ext_vel, vel_idx);
    const double *vmap_d = gkyl_array_cfetch(gvm->vmap, linidx_vel);

    // At each phase space cell, loop over all nodes and check if they are trapped or not. If any node is not trapped, the whole cell is not trapped, so we can break out of the loop early in that case.
    bool cell_trapped = true;
    for (int node = 0; node < num_phase_nodes && cell_trapped; ++node) {
      int conf_node = node / num_vel_nodes;
      int vel_node = node % num_vel_nodes;

      // Evaluate velocity coordinates at the node
      double xmu[GKYL_MAX_DIM] = {0.0};
      for (int vd = 0; vd < vdim; ++vd) {
        double vel_eta[GKYL_MAX_DIM] = {0.0};
        nodal_coords(vdim, vel_node, vel_eta);
        double xcomp[1] = {vel_eta[vd]};
        xmu[cdim + vd] =
          gvm->vmap_basis->eval_expand(xcomp, vmap_d + vd * gvm->vmap_basis->num_basis);
      }
      double mu = xmu[cdim + 1];
      double vpar = xmu[cdim];

      // Evaluate Hamiltonian at the node
      long linidx_conf = gkyl_range_idx(conf_range, conf_idx);
      double bmag_curr =
        field_node_val(bmag, up->basis_at_nodes_conf, num_basis_conf, linidx_conf, conf_node);
      double phi_curr =
        field_node_val(phi, up->basis_at_nodes_conf, num_basis_conf, linidx_conf, conf_node);
      double kinetic_energy = 0.5 * up->mass * vpar * vpar;
      double magnetic_energy = mu * bmag_curr;
      double electric_energy = up->charge * phi_curr;

      // Determine escape barriers at the node
      int zdim = cdim - 1;
      double barrier_left, barrier_right;
      escape_barriers(
        cdim, num_basis_conf, conf_range, up->basis_at_nodes_conf, phi, bmag, phi_wall_lo,
        phi_wall_up, conf_idx, conf_idx[zdim], conf_node, mu, up->charge, up->lower_boundary,
        up->upper_boundary, &barrier_left, &barrier_right
      );

      // Equality means that the particle reaches the wall with zero parallel
      // kinetic energy, matching the sheath boundary's absorbing cutoff.
      cell_trapped =
        hamiltonian_below_barrier(kinetic_energy, magnetic_energy, electric_energy, barrier_left) &&
        hamiltonian_below_barrier(kinetic_energy, magnetic_energy, electric_energy, barrier_right);
    }

    long linidx_phase = gkyl_range_idx(phase_range, phase_iter.idx);
    double *mask_d = gkyl_array_fetch(mask_out, linidx_phase);
    mask_d[0] = cell_trapped ? 1.0 : 0.0;
  }
}

void gkyl_loss_cone_mask_gyrokinetic_release(gkyl_loss_cone_mask_gyrokinetic *up)
{
  gkyl_velocity_map_release(up->vel_map);
  gkyl_array_release(up->basis_at_nodes_conf);
  gkyl_free(up);
}
