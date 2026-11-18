#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_util.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_vlasov_dist_moments.h>
#include <gkyl_vlasov_dist_moments_priv.h>
#include <gkyl_dg_vlasov_pressure.h>

#include <assert.h>

struct gkyl_vlasov_multi_temp_max_moments
{
  struct gkyl_vlasov_dist_moments dist_moms;
  struct gkyl_array *M0;
  struct gkyl_array *M1i;
  struct gkyl_array *M2ij;
  struct gkyl_array *V_drift;
  struct gkyl_array *pressure_tensor;
  struct gkyl_array *temperature;
  struct gkyl_dg_bin_op_mem *mem;
  struct gkyl_mom_calc *M0_calc;
  struct gkyl_mom_calc *M1i_calc;
  struct gkyl_mom_calc *M2ij_calc;
  struct gkyl_dg_vlasov_pressure *pressure_calc;
};

// Index of the symmetric tensor component (i,j).
static inline int
tensor_idx(int i, int j, int vdim)
{
  return i*vdim - i*(i-1)/2 + (j-i);
}

static void
vlasov_multi_temp_max_density_moment_advance(struct gkyl_vlasov_dist_moments *dist_moms,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local,
  const struct gkyl_array *fin, struct gkyl_array *density_out)
{
  struct gkyl_vlasov_multi_temp_max_moments *mtm_moms = (struct gkyl_vlasov_multi_temp_max_moments *) dist_moms;

  // Compute lab frame moment M0.
  gkyl_mom_calc_advance(mtm_moms->M0_calc, phase_local, conf_local, 
    fin, mtm_moms->M0);

  gkyl_array_set_range(density_out, 1.0, mtm_moms->M0, conf_local);
}

static void
vlasov_multi_temp_max_moments_advance(struct gkyl_vlasov_dist_moments *dist_moms,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local,
  const struct gkyl_array *fin, struct gkyl_array *moms_out)
{
  struct gkyl_vlasov_multi_temp_max_moments *mtm_moms = (struct gkyl_vlasov_multi_temp_max_moments *) dist_moms;

  int vdim = dist_moms->vdim;
  int num_conf_basis = dist_moms->num_conf_basis;

  // Compute lab frame moments M0 and M1i.
  gkyl_mom_calc_advance(mtm_moms->M0_calc,phase_local, conf_local, 
    fin, mtm_moms->M0);
  gkyl_mom_calc_advance(mtm_moms->M1i_calc, phase_local, conf_local, 
    fin, mtm_moms->M1i);

  // Isolate drift velocity by dividing M1i by M0.
  for (int d = 0; d < vdim; ++d) {
    gkyl_dg_div_op_range(mtm_moms->mem, &dist_moms->conf_basis,
      d, mtm_moms->V_drift,
      d, mtm_moms->M1i, 0, mtm_moms->M0, conf_local);
  }

  // Compute lab frame  M2ij.
  gkyl_mom_calc_advance(mtm_moms->M2ij_calc, phase_local, conf_local, 
    fin, mtm_moms->M2ij);

  // Compute pressure tensor:
  gkyl_dg_vlasov_pressure_tensor(mtm_moms->pressure_calc, conf_local, mtm_moms->M2ij, mtm_moms->M1i,
    mtm_moms->V_drift, mtm_moms->pressure_tensor);

  // Extract the diagonal pressure components and compute
  // T_i/m = P_ii/(m*M0).
  for (int d = 0; d < vdim; ++d) {
    int comp = tensor_idx(d, d, vdim);

    gkyl_dg_div_op_range(mtm_moms->mem, &dist_moms->conf_basis, d, mtm_moms->temperature,
      comp, mtm_moms->pressure_tensor, 0, mtm_moms->M0, conf_local);
  }

  gkyl_array_scale_range(mtm_moms->temperature, 1.0/dist_moms->mass, conf_local);

  // Save outputs to moms_out (n, V_1, ..., V_vdim, T_1/m, ..., T_vdim/m).
  gkyl_array_set_range(moms_out, 1.0, mtm_moms->M0, conf_local);
  gkyl_array_set_offset_range(moms_out, 1.0, mtm_moms->V_drift, num_conf_basis, conf_local);
  gkyl_array_set_offset_range(moms_out, 1.0, mtm_moms->temperature, (vdim+1)*num_conf_basis, conf_local);
}

static void
vlasov_multi_temp_max_moments_release(struct gkyl_vlasov_dist_moments *dist_moms)
{
  struct gkyl_vlasov_multi_temp_max_moments *mtm_moms = (struct gkyl_vlasov_multi_temp_max_moments *) dist_moms;

  gkyl_array_release(mtm_moms->M0);
  gkyl_array_release(mtm_moms->M1i);
  gkyl_array_release(mtm_moms->M2ij);
  gkyl_array_release(mtm_moms->V_drift);
  gkyl_array_release(mtm_moms->pressure_tensor);
  gkyl_array_release(mtm_moms->temperature);
  gkyl_dg_bin_op_mem_release(mtm_moms->mem);
  gkyl_mom_calc_release(mtm_moms->M0_calc);
  gkyl_mom_calc_release(mtm_moms->M1i_calc);
  gkyl_mom_calc_release(mtm_moms->M2ij_calc);
  gkyl_dg_vlasov_pressure_release(mtm_moms->pressure_calc);
  gkyl_free(mtm_moms);
}

struct gkyl_vlasov_dist_moments *gkyl_vlasov_multi_temp_max_moments_new(
  const struct gkyl_vlasov_dist_moments_inp *inp)
{
  struct gkyl_vlasov_multi_temp_max_moments *mtm_moms = gkyl_malloc(sizeof(*mtm_moms));
  struct gkyl_vlasov_dist_moments *dist_moms = &mtm_moms->dist_moms;

  dist_moms->conf_basis = *inp->conf_basis;
  dist_moms->phase_basis = *inp->phase_basis;
  dist_moms->num_conf_basis = inp->conf_basis->num_basis;
  dist_moms->vdim = inp->phase_basis->ndim - inp->conf_basis->ndim;
  dist_moms->num_mom = 2*dist_moms->vdim + 1;
  dist_moms->dist_id = GKYL_DIST_TYPE_MULTI_TEMP_MAX;
  dist_moms->model_id = inp->model_id;
  dist_moms->mass = inp->mass;
  dist_moms->use_gpu = inp->use_gpu;
  dist_moms->density_moment = vlasov_multi_temp_max_density_moment_advance;
  dist_moms->moments = vlasov_multi_temp_max_moments_advance;
  dist_moms->release = vlasov_multi_temp_max_moments_release;

  int vdim = dist_moms->vdim;
  int num_conf_basis = dist_moms->num_conf_basis;
  int num_tensor_comp = vdim*(vdim+1)/2;

  long conf_local_ncells = inp->conf_range->volume;
  long conf_local_ext_ncells = inp->conf_range_ext->volume;

  if (inp->use_gpu) {
    mtm_moms->M0 = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_conf_basis, conf_local_ext_ncells);
    mtm_moms->M1i = gkyl_array_cu_dev_new(GKYL_DOUBLE, vdim*num_conf_basis, conf_local_ext_ncells);
    mtm_moms->M2ij = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_tensor_comp*num_conf_basis, conf_local_ext_ncells);
    mtm_moms->V_drift = gkyl_array_cu_dev_new(GKYL_DOUBLE, vdim*num_conf_basis, conf_local_ext_ncells);
    mtm_moms->pressure_tensor = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_tensor_comp*num_conf_basis, conf_local_ext_ncells);
    mtm_moms->temperature = gkyl_array_cu_dev_new(GKYL_DOUBLE, vdim*num_conf_basis, conf_local_ext_ncells);
    mtm_moms->mem = gkyl_dg_bin_op_mem_cu_dev_new(conf_local_ncells, num_conf_basis);
  }
  else {
    mtm_moms->M0 = gkyl_array_new(GKYL_DOUBLE, num_conf_basis, conf_local_ext_ncells);
    mtm_moms->M1i = gkyl_array_new(GKYL_DOUBLE, vdim*num_conf_basis, conf_local_ext_ncells);
    mtm_moms->M2ij = gkyl_array_new(GKYL_DOUBLE, num_tensor_comp*num_conf_basis, conf_local_ext_ncells);
    mtm_moms->V_drift = gkyl_array_new(GKYL_DOUBLE, vdim*num_conf_basis, conf_local_ext_ncells);
    mtm_moms->pressure_tensor = gkyl_array_new(GKYL_DOUBLE, num_tensor_comp*num_conf_basis, conf_local_ext_ncells);
    mtm_moms->temperature = gkyl_array_new(GKYL_DOUBLE, vdim*num_conf_basis, conf_local_ext_ncells);
    mtm_moms->mem = gkyl_dg_bin_op_mem_new(conf_local_ncells, num_conf_basis);
  }

  // Moment calculators.
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
  inp_mom.mom_type = GKYL_F_MOMENT_M2IJ;
  struct gkyl_mom_type *M2ij_type = gkyl_mom_vlasov_inew(&inp_mom);

  mtm_moms->M0_calc = gkyl_mom_calc_new(inp->phase_grid, M0_type, inp->use_gpu);
  mtm_moms->M1i_calc = gkyl_mom_calc_new(inp->phase_grid, M1_type, inp->use_gpu);
  mtm_moms->M2ij_calc = gkyl_mom_calc_new(inp->phase_grid, M2ij_type, inp->use_gpu);
  mtm_moms->pressure_calc = gkyl_dg_vlasov_pressure_new(inp->conf_basis, inp->conf_range,
    inp->conf_range_ext, vdim,dist_moms->mass, inp->use_gpu);

  gkyl_mom_type_release(M0_type);
  gkyl_mom_type_release(M1_type);
  gkyl_mom_type_release(M2ij_type);

  return dist_moms;
}
