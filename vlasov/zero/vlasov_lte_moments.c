#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_calc_canonical_pb_vars.h>
#include <gkyl_mom_type.h> 
#include <gkyl_mom_vlasov.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_moments_priv.h>

#include <assert.h>

struct gkyl_vlasov_lte_moments*
gkyl_vlasov_lte_moments_inew(const struct gkyl_vlasov_lte_moments_inp *inp)
{
  gkyl_vlasov_lte_moments *up = gkyl_malloc(sizeof(*up));

  up->conf_basis = *inp->conf_basis;
  up->phase_basis = *inp->phase_basis;
  up->num_conf_basis = inp->conf_basis->num_basis;
  up->vdim = up->phase_basis.ndim - up->conf_basis.ndim;
  up->model_id = inp->model_id;
  up->use_extended_hamil_def = inp->use_extended_hamil_def;

  long conf_local_ncells = inp->conf_range->volume;
  long conf_local_ext_ncells = inp->conf_range_ext->volume;

  if (inp->use_gpu) {
    up->M0 = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->M1i = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->vdim*up->num_conf_basis, conf_local_ext_ncells);
    up->V_drift = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->vdim*up->num_conf_basis, conf_local_ext_ncells);
    up->V_drift_dot_M1i = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->pressure = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->temperature = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->mem = gkyl_dg_bin_op_mem_cu_dev_new(conf_local_ncells, up->num_conf_basis);
  }
  else {
    up->M0 = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->M1i = gkyl_array_new(GKYL_DOUBLE, up->vdim*up->num_conf_basis, conf_local_ext_ncells);
    up->V_drift = gkyl_array_new(GKYL_DOUBLE, up->vdim*up->num_conf_basis, conf_local_ext_ncells);
    up->V_drift_dot_M1i = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->pressure = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->temperature = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->mem = gkyl_dg_bin_op_mem_new(conf_local_ncells, up->num_conf_basis);
  }

  // Create moment updaters for computing LTE moments. 
  struct gkyl_mom_vlasov_inp inp_mom = {
    .conf_basis = inp->conf_basis,
    .phase_basis = inp->phase_basis,
    .vel_range = inp->vel_range,
    .vel_map = inp->vel_map,
    .hamil_range = inp->hamil_range,
    .hamil = inp->hamil,
    .model_id = inp->model_id,
    .hamil_id = inp->hamil_id,
    .use_gpu = inp->use_gpu,
  };  
  inp_mom.mom_type = GKYL_F_MOMENT_M0;
  struct gkyl_mom_type *M0_type = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M1;
  struct gkyl_mom_type *M1_type = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M2;
  struct gkyl_mom_type *M2_type = gkyl_mom_vlasov_inew(&inp_mom);

  up->M0_calc = gkyl_mom_calc_new(inp->phase_grid, M0_type, inp->use_gpu);
  up->M1i_calc = gkyl_mom_calc_new(inp->phase_grid, M1_type, inp->use_gpu);
  up->Pcalc = gkyl_mom_calc_new(inp->phase_grid, M2_type, inp->use_gpu);

  gkyl_mom_type_release(M0_type); 
  gkyl_mom_type_release(M1_type); 
  gkyl_mom_type_release(M2_type); 

  // Allocate model-specific data for obtaining LTE moments from base moments. 
  if (up->model_id == GKYL_MODEL_SR || up->model_id == GKYL_MODEL_TRIAD_GR) {
    // spatial components of the four-velocity squared, u_i^2 = (GammaV*V_drift)^2
    // and bulk four-velocity Lorentz boost factor GammaV = sqrt(1 + |u_i|^2) and its square
    if (inp->use_gpu) {
      up->V_drift_sq = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->vdim*up->num_conf_basis, conf_local_ext_ncells);
      up->GammaV = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      up->GammaV_sq = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    }
    else {
      up->V_drift_sq = gkyl_array_new(GKYL_DOUBLE, up->vdim*up->num_conf_basis, conf_local_ext_ncells);
      up->GammaV = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      up->GammaV_sq = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    }
    up->hamil = gkyl_array_acquire(inp->hamil); 
    up->gamma_inv = gkyl_array_acquire(inp->gamma_inv); 
    up->sr_vars = gkyl_dg_calc_sr_vars_new(inp->phase_grid, inp->vel_grid,
      inp->conf_basis, inp->vel_basis, inp->conf_range, inp->vel_range,
      inp->vel_map, inp->use_gpu);
  }
  else if (up->model_id == GKYL_MODEL_CANONICAL_PB) {
    up->h_ij = gkyl_array_acquire(inp->h_ij);
    up->h_ij_inv = gkyl_array_acquire(inp->h_ij_inv);
    up->det_h = gkyl_array_acquire(inp->det_h);
    if (inp->use_gpu) {
      up->energy_moment = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      up->M1i_cov = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->vdim*up->num_conf_basis, conf_local_ext_ncells);
      up->V_drift_cov = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->vdim*up->num_conf_basis, conf_local_ext_ncells);
    }
    else {
      up->energy_moment = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      up->M1i_cov = gkyl_array_new(GKYL_DOUBLE, up->vdim*up->num_conf_basis, conf_local_ext_ncells);
      up->V_drift_cov = gkyl_array_new(GKYL_DOUBLE, up->vdim*up->num_conf_basis, conf_local_ext_ncells);
    }

    // Variables needed for extended hamiltonians 
    if (up->use_extended_hamil_def) {
      up->effective_potential = gkyl_array_acquire(inp->effective_potential);
      if (inp->use_gpu) {
        up->M0_V = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      }
      else {
        up->M0_V = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      }
    }

    // Temperature moment is modified by can-pb, requires computing g^{ij}w_iw_j kernel
    up->can_pb_vars = gkyl_dg_calc_canonical_pb_vars_new(inp->phase_grid, inp->conf_basis, inp->phase_basis, inp->use_gpu);
  }
  return up;
}

void 
gkyl_vlasov_lte_density_moment_advance(struct gkyl_vlasov_lte_moments *lte_moms, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local, 
  const struct gkyl_array *fin, struct gkyl_array *density_out)
{
  // compute lab frame moment M0
  gkyl_mom_calc_advance(lte_moms->M0_calc, phase_local, conf_local, fin, lte_moms->M0); 

  // If we are relativistic, compute M1i and find the rest-frame density 
  // n = Gamma_inv*M0 where Gamma_inv = sqrt(1 - |V_drift|^2) and V_drift = M1i/M0
  if (lte_moms->model_id == GKYL_MODEL_SR || lte_moms->model_id == GKYL_MODEL_TRIAD_GR) {
    gkyl_mom_calc_advance(lte_moms->M1i_calc, phase_local, conf_local, fin, lte_moms->M1i); 
    gkyl_dg_calc_sr_vars_n(lte_moms->sr_vars, 
      lte_moms->M0, lte_moms->M1i, density_out);
  }
  else {
    gkyl_array_set_range(density_out, 1.0, lte_moms->M0, conf_local);
  }
}

void 
gkyl_vlasov_lte_moments_advance(struct gkyl_vlasov_lte_moments *lte_moms, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local, 
  const struct gkyl_array *fin, struct gkyl_array *moms_out)
{
  int vdim = lte_moms->vdim;
  int num_conf_basis = lte_moms->num_conf_basis;

  // compute lab frame moments M0 and M1i
  gkyl_mom_calc_advance(lte_moms->M0_calc, phase_local, conf_local, fin, lte_moms->M0); 
  gkyl_mom_calc_advance(lte_moms->M1i_calc, phase_local, conf_local, fin, lte_moms->M1i); 

  if (lte_moms->model_id == GKYL_MODEL_SR || lte_moms->model_id == GKYL_MODEL_TRIAD_GR) {
    // If we are relativistic, first compute rest-frame density n = Gamma_inv*M0,
    // Gamma_inv = sqrt(1 - |V_drift|^2), and V_drift = M1i/M0 (using weak division).
    // Done as a separate operator for robustness checks which insure V_drift < c.
    gkyl_dg_calc_sr_vars_n(lte_moms->sr_vars, 
      lte_moms->M0, lte_moms->M1i, moms_out); 

    // Isolate spatial component of the bulk four-velocity u_i = M1i/n = GammaV*V_drift
    // We store the output in the common V_drift array since we return the spatial
    // component of the bulk four-velocity in the LTE moms array in relativity. 
    for (int d = 0; d < vdim; ++d) {
      gkyl_dg_div_op_range(lte_moms->mem, &lte_moms->conf_basis, 
        d, lte_moms->V_drift,
        d, lte_moms->M1i, 0, moms_out, conf_local);
    }

    // Compute needed quantities for pressure velocity moment including 
    // bulk four-velocity Lorentz boost factor GammaV = sqrt(1 + |u_i|^2)
    gkyl_dg_calc_sr_vars_GammaV(lte_moms->sr_vars, conf_local, 
      lte_moms->V_drift, lte_moms->V_drift_sq, 
      lte_moms->GammaV, lte_moms->GammaV_sq);     

    // Compute the pressure moment.
    // This moment is computed *in the stationary frame* with the appropriate weight.
    // We find computing this moment *in the stationary frame* to be more accurate than 
    // computing the lab frame moment and then Lorentz transforming to the stationary frame. 
    gkyl_dg_calc_sr_vars_pressure(lte_moms->sr_vars, 
      conf_local, phase_local, 
      lte_moms->hamil, lte_moms->gamma_inv, 
      lte_moms->V_drift, lte_moms->V_drift_sq, 
      lte_moms->GammaV, lte_moms->GammaV_sq, 
      fin, lte_moms->pressure);  
  }
  else {
    // Isolate drift velocity by dividing M1i by M0
    // (For Canonical-pb only: This actually computes Jv*nv/(Jv*n) eliminating Jn)
    for (int d = 0; d < vdim; ++d) {
      gkyl_dg_div_op_range(lte_moms->mem, &lte_moms->conf_basis, 
        d, lte_moms->V_drift,
        d, lte_moms->M1i, 0, lte_moms->M0, conf_local);
    }

    // For non-relativistic Can-pb only
    if (lte_moms->model_id == GKYL_MODEL_CANONICAL_PB) {

      // M1i computes contravariant u^i, covert here to covariant u_i
      // u_i = h_{ij}u^i (where V_drift = u^i, M1i = Jnu^i)
      // This functions converts both V_drift and M1i to covariant components
      gkyl_canonical_pb_contra_to_covariant_m1i(lte_moms->can_pb_vars, conf_local, lte_moms->h_ij, 
        lte_moms->V_drift, lte_moms->M1i, lte_moms->V_drift_cov, lte_moms->M1i_cov);

      // Compute MEnergy
      gkyl_mom_calc_advance(lte_moms->Pcalc, phase_local, conf_local, fin, lte_moms->energy_moment);
      // Take covaraint momentum components and M1i components: V = 0 for non-extended can-pb cases
      // Solve for d*P*Jv: d*P*Jv + 2(Jv*n*V)  = 2*E - Jv*n*h^{ij}*u_i*u_j 
      //                                       = 2*E - h^{ij}*M1_i*V_drift_j
      gkyl_canonical_pb_pressure(lte_moms->can_pb_vars, conf_local, lte_moms->h_ij_inv, lte_moms->energy_moment,
        lte_moms->V_drift_cov, lte_moms->M1i_cov, lte_moms->pressure);

      // Subtract off the effective potential in place
      //d*P*Jv = (d*P*Jv + 2(Jv*n*V)) - 2(Jv*n*V)
      if (lte_moms->use_extended_hamil_def) {
        gkyl_dg_mul_op_range(&lte_moms->conf_basis, 0, lte_moms->M0_V, 0, lte_moms->M0, 0, lte_moms->effective_potential, conf_local);
        gkyl_array_accumulate_range(lte_moms->pressure, -2.0, lte_moms->M0_V, conf_local);
      }

      // Rescale pressure by 2.0/vdim and set the first component of moms_out to be the density. 
      gkyl_array_scale(lte_moms->pressure, 1.0/vdim);
    }
    else if (lte_moms->model_id == GKYL_MODEL_CANONICAL_PB_GR) {

      // Current no pressure/temperature or energy moments are defined for the
      // GR can-pb species
      assert(false);
    }
    else {
      // Compute the lab frame M2 = 0.5*(vdim*P/m + V_drift dot M1i).
      gkyl_mom_calc_advance(lte_moms->Pcalc, phase_local, conf_local, fin, lte_moms->pressure);
      // Subtract off V_drift dot M1i from total M2
      gkyl_array_clear(lte_moms->V_drift_dot_M1i, 0.0);
      gkyl_dg_dot_product_op_range(&lte_moms->conf_basis, 
        lte_moms->V_drift_dot_M1i, lte_moms->V_drift, lte_moms->M1i, conf_local); 
      gkyl_array_accumulate_range(lte_moms->pressure, -0.5, 
        lte_moms->V_drift_dot_M1i, conf_local); 
      // Rescale pressure by 2.0/vdim and set the first component of moms_out to be the density. 
      gkyl_array_scale(lte_moms->pressure, 2.0/vdim);
    }
    gkyl_array_set_range(moms_out, 1.0, lte_moms->M0, conf_local);
  }
  // ( T/m = P/(mn) ) 
  gkyl_dg_div_op_range(lte_moms->mem, &lte_moms->conf_basis, 
    0, lte_moms->temperature,
    0, lte_moms->pressure, 0, moms_out, conf_local);

  // Save the outputs to moms_out (n, V_drift, T/m):
  if (lte_moms->model_id == GKYL_MODEL_CANONICAL_PB || lte_moms->model_id == GKYL_MODEL_CANONICAL_PB_GR) {
    gkyl_array_set_offset_range(moms_out, 1.0, lte_moms->V_drift_cov, 1*num_conf_basis, conf_local);
  }
  else {
    gkyl_array_set_offset_range(moms_out, 1.0, lte_moms->V_drift, 1*num_conf_basis, conf_local);
  }
  gkyl_array_set_offset_range(moms_out, 1.0, lte_moms->temperature, (vdim+1)*num_conf_basis, conf_local);
}

void 
gkyl_vlasov_lte_moments_release(gkyl_vlasov_lte_moments *lte_moms)
{
  gkyl_array_release(lte_moms->M0);
  gkyl_array_release(lte_moms->M1i);
  gkyl_array_release(lte_moms->V_drift);
  gkyl_array_release(lte_moms->V_drift_dot_M1i);
  gkyl_array_release(lte_moms->pressure);
  gkyl_array_release(lte_moms->temperature);
  gkyl_dg_bin_op_mem_release(lte_moms->mem);
  if (lte_moms->use_extended_hamil_def) {
    gkyl_array_release(lte_moms->M0_V);
    gkyl_array_release(lte_moms->effective_potential);
  }
  if (lte_moms->model_id == GKYL_MODEL_SR || lte_moms->model_id == GKYL_MODEL_TRIAD_GR) {
    gkyl_dg_calc_sr_vars_release(lte_moms->sr_vars);
    gkyl_array_release(lte_moms->V_drift_sq);
    gkyl_array_release(lte_moms->GammaV);
    gkyl_array_release(lte_moms->GammaV_sq);
    gkyl_array_release(lte_moms->hamil);
    gkyl_array_release(lte_moms->gamma_inv);
  }
  else if (lte_moms->model_id == GKYL_MODEL_CANONICAL_PB) {
    // Only the non-relativistic canonical-PB model allocates the metric and
    // covariant-moment work arrays in the constructor (CANONICAL_PB_GR does not).
    gkyl_array_release(lte_moms->h_ij);
    gkyl_array_release(lte_moms->h_ij_inv);
    gkyl_array_release(lte_moms->det_h);
    gkyl_array_release(lte_moms->energy_moment);
    gkyl_dg_calc_canonical_pb_vars_release(lte_moms->can_pb_vars);
    gkyl_array_release(lte_moms->M1i_cov);
    gkyl_array_release(lte_moms->V_drift_cov);
  } 

  gkyl_mom_calc_release(lte_moms->M0_calc); 
  gkyl_mom_calc_release(lte_moms->M1i_calc); 
  gkyl_mom_calc_release(lte_moms->Pcalc); 
  gkyl_free(lte_moms);
}
