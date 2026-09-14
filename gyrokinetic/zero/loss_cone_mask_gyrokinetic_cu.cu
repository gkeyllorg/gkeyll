/* -*- c++ -*- */
#include <assert.h>

extern "C" {
#include <gkyl_array.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_loss_cone_mask_gyrokinetic_priv.h>
#include <gkyl_range.h>
#include <gkyl_velocity_map.h>
}

__global__ static void
gkyl_loss_cone_mask_gyrokinetic_advance_cu_ker(int cdim, int num_basis_conf,
  struct gkyl_range conf_corner_range, struct gkyl_range phase_corner_range,
  double mass, double charge,
  struct gkyl_range phase_range, struct gkyl_range conf_range,
  const struct gkyl_array *basis_at_corners_conf,
  const struct gkyl_array *bmag, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall_lo, const struct gkyl_array *phi_wall_up,
  enum gkyl_gk_trapped_passing_orbit_type lower_orbit,
  enum gkyl_gk_trapped_passing_orbit_type upper_orbit,
  const struct gkyl_velocity_map *gvm,
  struct gkyl_array *mask_out)
{
  int pdim = phase_range.ndim;
  int vdim = pdim - cdim;

  int phase_idx[GKYL_MAX_DIM];
  int conf_idx[GKYL_MAX_CDIM] = { 0 };
  int vel_idx[GKYL_MAX_VDIM] = { 0 };

  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x;
    tid < phase_range.volume; tid += blockDim.x * gridDim.x) {

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
    for (long corner = 0; corner < phase_corner_range.volume && cell_trapped; ++corner) {
      int corner_idx[GKYL_MAX_DIM] = { 0 };
      gkyl_range_inv_idx(&phase_corner_range, corner, corner_idx);
      int conf_corner = gkyl_range_idx(&conf_corner_range, corner_idx);

      double vpar_eta[1] = { corner_idx[cdim] ? 1.0 : -1.0 };
      double mu_eta[1] = { corner_idx[cdim + 1] ? 1.0 : -1.0 };
      double vpar = gvm->vmap_basis->eval_expand(vpar_eta, vmap_d);
      double mu = gvm->vmap_basis->eval_expand(mu_eta,
        vmap_d + gvm->vmap_basis->num_basis);

      long linidx_conf = gkyl_range_idx(&conf_range, conf_idx);
      double bmag_curr = gk_lcm_field_corner_value(bmag, basis_at_corners_conf,
        num_basis_conf, linidx_conf, conf_corner);
      double phi_curr = gk_lcm_field_corner_value(phi, basis_at_corners_conf,
        num_basis_conf, linidx_conf, conf_corner);
      double kinetic_energy = 0.5 * mass * vpar * vpar;
      double magnetic_energy = mu * bmag_curr;
      double electric_energy = charge * phi_curr;

      int zdim = cdim - 1;
      double barrier_left, barrier_right;
      gk_lcm_escape_barriers(cdim, num_basis_conf, &conf_range, &conf_corner_range,
        basis_at_corners_conf, phi, bmag, phi_wall_lo, phi_wall_up, conf_idx,
        conf_idx[zdim], conf_corner, mu, charge, lower_orbit, upper_orbit,
        &barrier_left, &barrier_right);

      cell_trapped = gk_lcm_hamiltonian_below_barrier(kinetic_energy,
        magnetic_energy, electric_energy, barrier_left)
        && gk_lcm_hamiltonian_below_barrier(kinetic_energy, magnetic_energy,
          electric_energy, barrier_right);
    }

    long linidx_phase = gkyl_range_idx(&phase_range, phase_idx);
    double *mask_d = (double *)gkyl_array_fetch(mask_out, linidx_phase);
    mask_d[0] = cell_trapped ? 1.0 : 0.0;
  }
}

void
gkyl_loss_cone_mask_gyrokinetic_advance_cu(gkyl_loss_cone_mask_gyrokinetic *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *bmag, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall_lo, const struct gkyl_array *phi_wall_up,
  struct gkyl_array *mask_out)
{
  int pdim = phase_range->ndim;
  int vdim = pdim - up->cdim;
  assert(vdim == 2);

  int corner_shape[GKYL_MAX_DIM];
  for (int d = 0; d < pdim; ++d) {
    corner_shape[d] = 2;
  }
  struct gkyl_range phase_corner_range;
  gkyl_range_init_from_shape(&phase_corner_range, pdim, corner_shape);

  const struct gkyl_velocity_map *gvm = up->vel_map;

  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;

  const struct gkyl_array *phi_wall_lo_dev = phi_wall_lo ? phi_wall_lo->on_dev : 0;
  const struct gkyl_array *phi_wall_up_dev = phi_wall_up ? phi_wall_up->on_dev : 0;

  gkyl_loss_cone_mask_gyrokinetic_advance_cu_ker<<<nblocks, nthreads>>>(
    up->cdim, up->num_basis_conf, up->conf_corner_range, phase_corner_range,
    up->mass, up->charge, *phase_range, *conf_range,
    up->basis_at_corners_conf->on_dev, bmag->on_dev, phi->on_dev,
    phi_wall_lo_dev, phi_wall_up_dev, up->lower_orbit, up->upper_orbit,
    gvm->on_dev, mask_out->on_dev);
}
