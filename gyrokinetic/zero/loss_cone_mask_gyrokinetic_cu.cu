/* -*- c++ -*- */
#include <float.h>
#include <math.h>

extern "C" {
#include <gkyl_array.h>
#include <gkyl_const.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_range.h>
#include <gkyl_velocity_map.h>
}

#include <gkyl_loss_cone_mask_gyrokinetic_priv.h>

__global__ static void gkyl_loss_cone_mask_gyrokinetic_advance_cu_ker(
  int cdim, int num_basis_conf, int num_phase_nodes, int num_vel_nodes, double mass, double charge,
  struct gkyl_range phase_range, struct gkyl_range conf_range,
  const struct gkyl_array *basis_at_nodes_conf, const struct gkyl_array *bmag,
  const struct gkyl_array *phi, const struct gkyl_array *phi_wall_lo,
  const struct gkyl_array *phi_wall_up, enum gkyl_loss_cone_boundary_type lower_boundary,
  enum gkyl_loss_cone_boundary_type upper_boundary, const struct gkyl_velocity_map *gvm,
  struct gkyl_array *mask_out
)
{
  int pdim = phase_range.ndim;
  int vdim = pdim - cdim;

  int phase_idx[GKYL_MAX_DIM];
  int conf_idx[GKYL_MAX_DIM] = {0};
  int vel_idx[GKYL_MAX_DIM] = {0};

  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < phase_range.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, phase_idx);

    for (int d = 0; d < cdim; ++d) {
      conf_idx[d] = phase_idx[d];
    }
    for (int d = 0; d < vdim; ++d) {
      vel_idx[d] = phase_idx[cdim + d];
    }

    long linidx_vel = gkyl_range_idx(&gvm->local_ext_vel, vel_idx);
    const double *vmap_d = (const double *)gkyl_array_cfetch(gvm->vmap, linidx_vel);

    bool cell_trapped = true;
    for (int node = 0; node < num_phase_nodes && cell_trapped; ++node) {
      int conf_node = node / num_vel_nodes;
      int vel_node = node % num_vel_nodes;

      double xmu[GKYL_MAX_DIM] = {0.0};
      double vel_eta[GKYL_MAX_DIM] = {0.0};
      nodal_coords(vdim, vel_node, vel_eta);

      for (int vd = 0; vd < vdim; ++vd) {
        double xcomp[1] = {vel_eta[vd]};
        xmu[cdim + vd] =
          gvm->vmap_basis->eval_expand(xcomp, vmap_d + vd * gvm->vmap_basis->num_basis);
      }

      double mu = xmu[cdim + 1];
      double vpar = xmu[cdim];

      long linidx_conf = gkyl_range_idx(&conf_range, conf_idx);
      double bmag_curr =
        field_node_val(bmag, basis_at_nodes_conf, num_basis_conf, linidx_conf, conf_node);
      double phi_curr =
        field_node_val(phi, basis_at_nodes_conf, num_basis_conf, linidx_conf, conf_node);
      double kinetic_energy = 0.5 * mass * vpar * vpar;
      double magnetic_energy = mu * bmag_curr;
      double electric_energy = charge * phi_curr;

      int zdim = cdim - 1;
      double barrier_left, barrier_right;
      escape_barriers(
        cdim, num_basis_conf, &conf_range, basis_at_nodes_conf, phi, bmag, phi_wall_lo, phi_wall_up,
        conf_idx, conf_idx[zdim], conf_node, mu, charge, lower_boundary, upper_boundary,
        &barrier_left, &barrier_right
      );

      cell_trapped =
        hamiltonian_below_barrier(kinetic_energy, magnetic_energy, electric_energy, barrier_left) &&
        hamiltonian_below_barrier(kinetic_energy, magnetic_energy, electric_energy, barrier_right);
    }

    long linidx_phase = gkyl_range_idx(&phase_range, phase_idx);
    double *mask_d = (double *)gkyl_array_fetch(mask_out, linidx_phase);
    mask_d[0] = cell_trapped ? 1.0 : 0.0;
  }
}

extern "C" void gkyl_loss_cone_mask_gyrokinetic_advance_cu(
  gkyl_loss_cone_mask_gyrokinetic *up, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range, const struct gkyl_array *bmag, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall_lo, const struct gkyl_array *phi_wall_up,
  struct gkyl_array *mask_out
)
{
  int pdim = phase_range->ndim;
  int vdim = pdim - up->cdim;
  int num_phase_nodes = (int)pow(2.0, (double)pdim);
  int num_vel_nodes = (int)pow(2.0, (double)vdim);

  const struct gkyl_velocity_map *gvm = up->vel_map;

  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;

  const struct gkyl_array *phi_wall_lo_dev = phi_wall_lo ? phi_wall_lo->on_dev : 0;
  const struct gkyl_array *phi_wall_up_dev = phi_wall_up ? phi_wall_up->on_dev : 0;

  gkyl_loss_cone_mask_gyrokinetic_advance_cu_ker<<<nblocks, nthreads>>>(
    up->cdim, up->num_basis_conf, num_phase_nodes, num_vel_nodes, up->mass, up->charge,
    *phase_range, *conf_range, up->basis_at_nodes_conf->on_dev, bmag->on_dev, phi->on_dev,
    phi_wall_lo_dev, phi_wall_up_dev, up->lower_boundary, up->upper_boundary, gvm->on_dev,
    mask_out->on_dev
  );
}
