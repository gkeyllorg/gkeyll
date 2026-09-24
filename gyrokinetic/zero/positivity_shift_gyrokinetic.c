#include <gkyl_positivity_shift_gyrokinetic.h>
#include <gkyl_positivity_shift_gyrokinetic_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <float.h>

struct gkyl_positivity_shift_gyrokinetic*
gkyl_positivity_shift_gyrokinetic_new(struct gkyl_basis cbasis, struct gkyl_basis pbasis,
  struct gkyl_rect_grid grid, double mass,
  const struct gk_geometry *gk_geom,
  const struct gkyl_velocity_map *vel_map, const struct gkyl_range *conf_rng_ext, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_positivity_shift_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  assert(pbasis.poly_order == 1); // Because of the way a rescale/division is
                                  // done in advance.

  up->ffloor_fac = 0.0; // ffloor will be set to max(f)*ffloor_fac.
  up->grid = grid;
  up->num_cbasis = cbasis.num_basis;
  up->mass = mass;

  up->gk_geom = gkyl_gk_geometry_acquire(gk_geom);
  up->vel_map = gkyl_velocity_map_acquire(vel_map);
  up->use_gpu = use_gpu;
  up->cellav_fac = 1./pow(sqrt(2.),pbasis.ndim);

  if (!use_gpu) {
    up->kernels = gkyl_malloc(sizeof(struct gkyl_positivity_shift_gyrokinetic_kernels));

    up->ffloor = gkyl_malloc(sizeof(double[1]));
    up->ffloor[0] = 0.0;  // Gets updated after 1st call to _advance.
  }
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->kernels = gkyl_cu_malloc(sizeof(struct gkyl_positivity_shift_gyrokinetic_kernels));

    up->ffloor = gkyl_cu_malloc(sizeof(double[1]));
    double ffloor_zero[] = {0.};  // Gets updated after 1st call to _advance.
    gkyl_cu_memcpy(up->ffloor, ffloor_zero, sizeof(double[1]), GKYL_CU_MEMCPY_H2D);

    up->shiftedf = gkyl_array_cu_dev_new(GKYL_INT, 1, conf_rng_ext->volume);
  }
#endif

  // Choose kernels that shift f and compute int moms of Deltaf.
  
  enum gkyl_positivity_shift_type shift_type = GKYL_POSITIVITY_SHIFT_TYPE_SHIFT_ONLY;
  pos_shift_gk_choose_shift_kernel(up->kernels, cbasis, pbasis, shift_type, use_gpu);

  return up;
}

void
gkyl_positivity_shift_gyrokinetic_advance(gkyl_positivity_shift_gyrokinetic* up,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  struct gkyl_array *GKYL_RESTRICT distf, struct gkyl_array *GKYL_RESTRICT m0,
  struct gkyl_array *GKYL_RESTRICT delta_m0)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_positivity_shift_gyrokinetic_advance_cu(up, conf_rng, phase_rng,
      distf, m0, delta_m0);
    return;
  }
#endif

  double distf_max = -DBL_MAX;

  int num_cbasis = up->num_cbasis;

  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  int rem_dir[GKYL_MAX_DIM] = {0};
  for (int d=0; d<conf_rng->ndim; ++d) rem_dir[d] = 1;

  gkyl_range_iter_init(&conf_iter, conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long clinidx = gkyl_range_idx(conf_rng, conf_iter.idx);

    const double *bmag_c = gkyl_array_cfetch(up->gk_geom->geo_int.bmag, clinidx);
    const double *jacobtot_c = gkyl_array_cfetch(up->gk_geom->geo_int.jacobtot, clinidx);
    const double *jacobtot_inv_c = gkyl_array_cfetch(up->gk_geom->geo_int.jacobtot_inv, clinidx);

    double *m0_c = gkyl_array_fetch(m0, clinidx);
    double *delta_m0_c = gkyl_array_fetch(delta_m0, clinidx);
    double m0in_c[num_cbasis];
    for (int k=0; k<num_cbasis; k++) {
      m0_c[k] = 0.0;
      delta_m0_c[k] = 0.0;
      m0in_c[k] = 0.0;
    }

    bool shiftedf = false;

    gkyl_range_deflate(&vel_rng, phase_rng, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);

    while (gkyl_range_iter_next(&vel_iter)) {
      long plinidx = gkyl_range_idx(&vel_rng, vel_iter.idx);
      double *distf_c = gkyl_array_fetch(distf, plinidx);
      double *jacobvel_c = gkyl_array_fetch(up->vel_map->jacobvel, plinidx);

      long vlinidx = gkyl_range_idx(&up->vel_map->local_vel, vel_iter.idx);
      const double *vmap_c = gkyl_array_cfetch(up->vel_map->vmap, vlinidx);

      // Contribution to the old number density from this v-space cell.
      double m0phase_in_c[num_cbasis];
      for (int k=0; k<num_cbasis; k++)
        m0phase_in_c[k] = 0.0;
      up->kernels->m0(up->grid.dx, vmap_c, up->mass, bmag_c, distf_c, m0phase_in_c);

      // Add to the old number density.
      for (int k=0; k<num_cbasis; k++)
        m0in_c[k] += m0phase_in_c[k];

      // Shift f if needed.
      bool shifted_node = false;

      // Divide by jacobtot and jacobvel so that we are shifting just f.
      up->kernels->conf_phase_mul_op(jacobtot_inv_c, distf_c, distf_c);
      for (int k=0; k<distf->ncomp; k++)
        distf_c[k] /= jacobvel_c[0];
      // Shift f to enforce positivity if needed.
      shifted_node = up->kernels->shift(up->ffloor[0], distf_c);
      // Multiply by jacobtot and jacobvel to compute M0.
      up->kernels->conf_phase_mul_op(jacobtot_c, distf_c, distf_c);
      for (int k=0; k<distf->ncomp; k++)
        distf_c[k] *= jacobvel_c[0];

      if (shifted_node) {
        // Compute the new number density in this phase-space cell.
        double m0phase_out_c[num_cbasis];
        for (int k=0; k<num_cbasis; k++)
          m0phase_out_c[k] = 0.0;
        up->kernels->m0(up->grid.dx, vmap_c, up->mass, bmag_c, distf_c, m0phase_out_c);

        if (m0phase_in_c[0] > 0.0 && m0phase_out_c[0] > 0.0) {
          // Rescale f in this cell so it keeps the same cell-averaged density.
          double m0ratio = m0phase_in_c[0]/m0phase_out_c[0];

          for (unsigned int k=0; k<distf->ncomp; ++k)
            distf_c[k] *= m0ratio;

          // Add contribution from this phase-space cell to the new number density.
          for (unsigned int k = 0; k < m0->ncomp; ++k)
            m0_c[k] += m0ratio*m0phase_out_c[k];
        }
        else {
          // Add contribution from this phase-space cell to the new number density.
          for (int k=0; k<num_cbasis; k++)
            m0_c[k] += m0phase_out_c[k];

          shiftedf = true;
        }
      }
      else {
        // Add contribution from this phase-space cell to the new number density.
        for (int k=0; k<num_cbasis; k++)
          m0_c[k] += m0phase_in_c[k];
      }

      distf_max = GKYL_MAX2(distf_max, distf_c[0]);
    }

    if (shiftedf) {
      if (up->kernels->is_m0_positive(m0in_c)) {
        // Rescale f so it has the same m0 at this conf-space cell.
        double m0ratio_c[num_cbasis];
        up->kernels->conf_inv_op(m0_c, m0ratio_c);
        up->kernels->conf_mul_op(m0in_c, m0ratio_c, m0ratio_c);

        gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
        while (gkyl_range_iter_next(&vel_iter)) {
          long plinidx = gkyl_range_idx(&vel_rng, vel_iter.idx);
          double *distf_c = gkyl_array_fetch(distf, plinidx);
          up->kernels->conf_phase_mul_op(m0ratio_c, distf_c, distf_c);
        }

        for (int k=0; k<num_cbasis; k++)
          m0_c[k] = m0in_c[k];
      }
      else {
        for (int k=0; k<num_cbasis; k++)
          delta_m0_c[k] = m0_c[k] - m0in_c[k];
      }
    }
  }

  up->ffloor[0] = up->ffloor_fac * distf_max * up->cellav_fac;
}

void
gkyl_positivity_shift_gyrokinetic_quasineutrality_scale(gkyl_positivity_shift_gyrokinetic* up,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  const struct gkyl_array *GKYL_RESTRICT delta_m0s, const struct gkyl_array *GKYL_RESTRICT delta_m0s_tot,
  const struct gkyl_array *GKYL_RESTRICT delta_m0r_tot, const struct gkyl_array *GKYL_RESTRICT m0s,
  struct gkyl_array *GKYL_RESTRICT fs)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_positivity_shift_gyrokinetic_quasineutrality_scale_cu(up, conf_rng, phase_rng,
      delta_m0s, delta_m0s_tot, delta_m0r_tot, m0s, fs);
    return;
  }
#endif

  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  int rem_dir[GKYL_MAX_DIM] = {0};
  for (int d=0; d<conf_rng->ndim; ++d) rem_dir[d] = 1;

  int num_cbasis = up->num_cbasis;

  gkyl_range_iter_init(&conf_iter, conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long clinidx = gkyl_range_idx(conf_rng, conf_iter.idx);

    const double *delta_m0r_tot_c = gkyl_array_cfetch(delta_m0r_tot, clinidx);
    const double *delta_m0s_tot_c = gkyl_array_cfetch(delta_m0s_tot, clinidx);

    // First condition in this if-statement is equivalent to
    //   max((delta_m0r_tot_c[0]-delta_m0s_tot_c[0])/abs(delta_m0r_tot_c[0]-delta_m0s_tot_c[0]))
    // and the second condition is to avoid division by 0.
    if (delta_m0r_tot_c[0] > delta_m0s_tot_c[0] &&
        delta_m0s_tot_c[0] > 0.0) {

      const double *delta_m0s_c = gkyl_array_cfetch(delta_m0s, clinidx);
      const double *m0s_c = gkyl_array_cfetch(m0s, clinidx);

      // Compute the scaling factor (n_s+h)/n_s with h=(Delta n_r,tot - Delta_n_s,tot) * Delta n_s/Delta_n_s,tot
      //   - Delta n_s = shift in density of this species due to positivity shift.
      //   - Delta n_s,tot = sum of Delta n for species with same charge sign.
      //   - Delta n_r,tot = sum of Delta n for species with opposite charge sign.
      double delta_m0fac_c[num_cbasis];
      for (int k=0; k<num_cbasis; k++)
        delta_m0fac_c[k] = delta_m0r_tot_c[k] - delta_m0s_tot_c[k];

      up->kernels->conf_mul_op(delta_m0fac_c, delta_m0s_c, delta_m0fac_c);

      double delta_m0s_tot_inv_c[num_cbasis];
      up->kernels->conf_inv_op(delta_m0s_tot_c, delta_m0s_tot_inv_c);

      up->kernels->conf_mul_op(delta_m0fac_c, delta_m0s_tot_inv_c, delta_m0fac_c);

      for (int k=0; k<num_cbasis; k++)
        delta_m0fac_c[k] += m0s_c[k];

      double m0s_inv_c[num_cbasis];
      up->kernels->conf_inv_op(m0s_c, m0s_inv_c);

      up->kernels->conf_mul_op(delta_m0fac_c, m0s_inv_c, delta_m0fac_c);

      // Scale the distribution function.
      gkyl_range_deflate(&vel_rng, phase_rng, rem_dir, conf_iter.idx);
      gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);

      while (gkyl_range_iter_next(&vel_iter)) {
        long plinidx = gkyl_range_idx(&vel_rng, vel_iter.idx);
        double *fs_c = gkyl_array_fetch(fs, plinidx);

        up->kernels->conf_phase_mul_op(delta_m0fac_c, fs_c, fs_c);
      }
    }

  }
}

void
gkyl_positivity_shift_gyrokinetic_release(gkyl_positivity_shift_gyrokinetic* up)
{
  // Release memory associated with this updater.
  gkyl_gk_geometry_release(up->gk_geom);
  gkyl_velocity_map_release(up->vel_map);
  if (!up->use_gpu) {
    gkyl_free(up->ffloor);
    gkyl_free(up->kernels);
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->ffloor);
    gkyl_cu_free(up->kernels);
    gkyl_array_release(up->shiftedf);
  }
#endif
  gkyl_free(up);
}
