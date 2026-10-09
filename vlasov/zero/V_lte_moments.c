#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_util.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_vlasov_dist_moments.h>
#include <gkyl_vlasov_dist_moments_priv.h>
#include <gkyl_dg_vlasov_pressure.h>

#include <assert.h>


struct gkyl_vlasov_lte_moments
{
  struct gkyl_vlasov_dist_moments dist_moms;
  struct gkyl_array *M0;
  struct gkyl_array *M1i;
  struct gkyl_array *M2;
  struct gkyl_array *V_drift;
  struct gkyl_array *pressure;
  struct gkyl_array *temperature;
  struct gkyl_dg_bin_op_mem *mem;
  struct gkyl_dg_updater_moment *M0_calc;
  struct gkyl_dg_updater_moment *M1i_calc;
  struct gkyl_dg_updater_moment *M2_calc;
  struct gkyl_dg_vlasov_pressure *pressure_calc;
};

static void
vlasov_lte_density_moment_advance(struct gkyl_vlasov_dist_moments *dist_moms,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local,
  const struct gkyl_array *fin, struct gkyl_array *density_out)
{
  struct gkyl_vlasov_lte_moments *lte_moms = (struct gkyl_vlasov_lte_moments *) dist_moms;

  // Compute lab frame moment M0.
  gkyl_dg_updater_moment_advance(lte_moms->M0_calc, phase_local, conf_local, 
    fin, lte_moms->M0);
  gkyl_array_set_range(density_out, 1.0, lte_moms->M0, conf_local);
}


static void
vlasov_lte_moments_advance(struct gkyl_vlasov_dist_moments *dist_moms,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local,
  const struct gkyl_array *fin, struct gkyl_array *moms_out)
{
  struct gkyl_vlasov_lte_moments *lte_moms = (struct gkyl_vlasov_lte_moments *) dist_moms;

  int vdim = dist_moms->vdim;
  int num_conf_basis = dist_moms->num_conf_basis;

  // Compute lab frame moments M0 and M1i.
  gkyl_dg_updater_moment_advance(lte_moms->M0_calc, phase_local, conf_local, 
    fin, lte_moms->M0);
  gkyl_dg_updater_moment_advance(lte_moms->M1i_calc, phase_local, conf_local, 
    fin, lte_moms->M1i);

  // Isolate drift velocity by dividing M1i by M0.
  for (int d = 0; d < vdim; ++d) {
    gkyl_dg_div_op_range(lte_moms->mem, &dist_moms->conf_basis,
      d, lte_moms->V_drift,
      d, lte_moms->M1i, 0, lte_moms->M0, conf_local);
  }

  // Compute the lab frame M2.
  gkyl_dg_updater_moment_advance(lte_moms->M2_calc, phase_local, conf_local,
    fin, lte_moms->M2);

  // Compute scalar pressure.
  gkyl_dg_vlasov_pressure_scalar(lte_moms->pressure_calc, conf_local, lte_moms->M2, lte_moms->M1i,
  lte_moms->V_drift, lte_moms->pressure);
  gkyl_array_set_range(moms_out, 1.0, lte_moms->M0, conf_local);

  // T = P/n.
  gkyl_dg_div_op_range(lte_moms->mem, &dist_moms->conf_basis,
    0, lte_moms->temperature,
    0, lte_moms->pressure, 0, lte_moms->M0, conf_local);

  // T/m = P/(mn).
  gkyl_array_scale_range(lte_moms->temperature, 1.0/dist_moms->mass, conf_local);
  // Save the outputs to moms_out (n, V_drift, T/m).
  gkyl_array_set_offset_range(moms_out, 1.0, lte_moms->V_drift, 1*num_conf_basis, conf_local);
  gkyl_array_set_offset_range(moms_out, 1.0, lte_moms->temperature, (vdim+1)*num_conf_basis, conf_local);
}

static void
vlasov_lte_moments_release(
  struct gkyl_vlasov_dist_moments *dist_moms)
{
  struct gkyl_vlasov_lte_moments *lte_moms = (struct gkyl_vlasov_lte_moments *) dist_moms;

  gkyl_array_release(lte_moms->M0);
  gkyl_array_release(lte_moms->M1i);
  gkyl_dg_updater_moment_release(lte_moms->M2_calc);
  gkyl_array_release(lte_moms->V_drift);
  gkyl_array_release(lte_moms->pressure);
  gkyl_array_release(lte_moms->temperature);
  gkyl_dg_bin_op_mem_release(lte_moms->mem);
  gkyl_dg_updater_moment_release(lte_moms->M0_calc);
  gkyl_dg_updater_moment_release(lte_moms->M1i_calc);
  gkyl_dg_vlasov_pressure_release(lte_moms->pressure_calc);
  gkyl_free(lte_moms);
}


struct gkyl_vlasov_dist_moments*
gkyl_vlasov_lte_moments_new(
  const struct gkyl_vlasov_dist_moments_inp *inp)
{
  struct gkyl_vlasov_lte_moments *lte_moms = gkyl_malloc(sizeof(*lte_moms));
  struct gkyl_vlasov_dist_moments *dist_moms = &lte_moms->dist_moms;

  dist_moms->conf_basis = *inp->conf_basis;
  dist_moms->phase_basis = *inp->phase_basis;
  dist_moms->num_conf_basis = inp->conf_basis->num_basis;
  dist_moms->vdim = inp->phase_basis->ndim - inp->conf_basis->ndim;
  dist_moms->num_mom = dist_moms->vdim + 2;
  dist_moms->dist_id = GKYL_DIST_TYPE_LTE;
  dist_moms->model_id = inp->model_id;
  dist_moms->mass = inp->mass;
  dist_moms->use_gpu = inp->use_gpu;
  dist_moms->density_moment = vlasov_lte_density_moment_advance;
  dist_moms->moments = vlasov_lte_moments_advance;
  dist_moms->release = vlasov_lte_moments_release;
  long conf_local_ncells = inp->conf_range->volume;
  long conf_local_ext_ncells = inp->conf_range_ext->volume;

  if (inp->use_gpu) {
    lte_moms->M0 = gkyl_array_cu_dev_new(GKYL_DOUBLE, dist_moms->num_conf_basis, conf_local_ext_ncells);
    lte_moms->M1i = gkyl_array_cu_dev_new(GKYL_DOUBLE, dist_moms->vdim*dist_moms->num_conf_basis, conf_local_ext_ncells);
    lte_moms->M2 = gkyl_array_cu_dev_new(GKYL_DOUBLE, dist_moms->num_conf_basis, conf_local_ext_ncells);
    lte_moms->V_drift = gkyl_array_cu_dev_new(GKYL_DOUBLE, dist_moms->vdim*dist_moms->num_conf_basis, conf_local_ext_ncells);
    lte_moms->pressure = gkyl_array_cu_dev_new(GKYL_DOUBLE, dist_moms->num_conf_basis, conf_local_ext_ncells);
    lte_moms->temperature = gkyl_array_cu_dev_new(GKYL_DOUBLE, dist_moms->num_conf_basis, conf_local_ext_ncells);
    lte_moms->mem = gkyl_dg_bin_op_mem_cu_dev_new(conf_local_ncells, dist_moms->num_conf_basis);
  }
  else {
    lte_moms->M0 = gkyl_array_new(GKYL_DOUBLE, dist_moms->num_conf_basis, conf_local_ext_ncells);
    lte_moms->M1i = gkyl_array_new(GKYL_DOUBLE, dist_moms->vdim*dist_moms->num_conf_basis, conf_local_ext_ncells);
    lte_moms->M2 = gkyl_array_new(GKYL_DOUBLE, dist_moms->num_conf_basis, conf_local_ext_ncells);
    lte_moms->V_drift = gkyl_array_new(GKYL_DOUBLE, dist_moms->vdim*dist_moms->num_conf_basis, conf_local_ext_ncells);
    lte_moms->pressure = gkyl_array_new(GKYL_DOUBLE, dist_moms->num_conf_basis, conf_local_ext_ncells);
    lte_moms->temperature = gkyl_array_new(GKYL_DOUBLE, dist_moms->num_conf_basis, conf_local_ext_ncells);
    lte_moms->mem = gkyl_dg_bin_op_mem_new(conf_local_ncells, dist_moms->num_conf_basis);
  }

  // Moment calculators for needed moments (M0, M1i, and M2 for non-relativistic).
  lte_moms->M0_calc = gkyl_dg_updater_moment_new(inp->phase_grid, inp->conf_basis,
      inp->phase_basis, inp->conf_range, inp->vel_range, inp->phase_range,
      dist_moms->model_id, 0, GKYL_F_MOMENT_M0, false, inp->use_gpu);

  lte_moms->M1i_calc = gkyl_dg_updater_moment_new(inp->phase_grid, inp->conf_basis,
      inp->phase_basis, inp->conf_range, inp->vel_range, inp->phase_range,
      dist_moms->model_id, 0, GKYL_F_MOMENT_M1, false, inp->use_gpu);

  lte_moms->M2_calc = gkyl_dg_updater_moment_new(inp->phase_grid, inp->conf_basis,
      inp->phase_basis, inp->conf_range, inp->vel_range, inp->phase_range,
      dist_moms->model_id, 0, GKYL_F_MOMENT_M2, false, inp->use_gpu);

  lte_moms->pressure_calc = gkyl_dg_vlasov_pressure_new(inp->conf_basis, inp->conf_range,
      inp->conf_range_ext, dist_moms->vdim, dist_moms->mass, inp->use_gpu);

  return dist_moms;
}