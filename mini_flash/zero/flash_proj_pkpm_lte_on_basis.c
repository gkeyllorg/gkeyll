#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_flash_proj_pkpm_lte_on_basis.h>
#include <gkyl_flash_proj_pkpm_lte_on_basis_priv.h>

struct gkyl_flash_proj_pkpm_lte_on_basis*
gkyl_flash_proj_pkpm_lte_on_basis_new(struct gkyl_flash_proj_pkpm_lte_on_basis_inp *inp)
{
  if (inp->use_gpu)
    gkyl_exit("gkyl_flash_proj_pkpm_lte_on_basis does not support GPUs.");

  struct gkyl_flash_proj_pkpm_lte_on_basis *up = gkyl_malloc(sizeof(*up));

  up->phase_grid = *inp->phase_grid;
  up->conf_basis = *inp->conf_basis;
  up->phase_basis = *inp->phase_basis;

  up->cdim = up->conf_basis.ndim;
  up->pdim = up->phase_basis.ndim;
  up->vdim = up->pdim - up->cdim;
  up->num_conf_basis = up->conf_basis.num_basis;
  up->num_phase_basis = up->phase_basis.num_basis;

  int num_quad = up->conf_basis.poly_order+3;
  bool is_vdim_p2[3];
  int num_quad_v = vel_quad_shape(&up->phase_basis, up->vdim, num_quad, is_vdim_p2);

  // initialize data needed for conf-space quadrature; only the basis at the
  // ordinates is needed.
  struct gkyl_array *conf_ordinates, *conf_weights;
  up->tot_conf_quad = init_quad_values(up->cdim, &up->conf_basis, 
    inp->quad_type, num_quad, num_quad_v, is_vdim_p2,
    &conf_ordinates, &conf_weights, &up->conf_basis_at_ords);
  gkyl_array_release(conf_ordinates);
  gkyl_array_release(conf_weights);

  // initialize data needed for phase-space quadrature 
  up->tot_quad = init_quad_values(up->cdim, &up->phase_basis, 
    inp->quad_type, num_quad, num_quad_v, is_vdim_p2,
    &up->ordinates, &up->weights, &up->basis_at_ords);

  up->f0_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);
  up->g_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);

  up->conf_qrange = get_qrange(up->cdim, up->cdim, num_quad, num_quad_v, is_vdim_p2);
  up->phase_qrange = get_qrange(up->cdim, up->pdim, num_quad, num_quad_v, is_vdim_p2);

  return up;
}

void
gkyl_flash_proj_pkpm_lte_on_basis_advance(gkyl_flash_proj_pkpm_lte_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *moms, struct gkyl_array *f_lte)
{
  double f_floor = 1.e-40; 
  int cdim = up->cdim, pdim = up->pdim, vdim = up->vdim;
  int tot_quad = up->tot_quad;
  int num_phase_basis = up->num_phase_basis;  

  int tot_conf_quad = up->tot_conf_quad;
  int num_conf_basis = up->num_conf_basis;

  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_range->ndim; ++d) rem_dir[d] = 1;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  double n_quad[tot_conf_quad], T_over_m_quad[tot_conf_quad];
  double expamp_quad[tot_conf_quad];

  // outer loop over configuration space cells; for each
  // config-space cell inner loop walks over velocity space
  gkyl_range_iter_init(&conf_iter, conf_range);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_range, conf_iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, midx);
    const double *n_d = moms_d;
    const double *T_over_m_d = &moms_d[num_conf_basis*(GKYL_MAX_VDIM+1)];

    // Evaluate n and T/m at the conf-space ordinates. (Note PKPM is in moving reference frame).
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = gkyl_array_cfetch(up->conf_basis_at_ords, n);

      // Zero out quadrature values
      n_quad[n] = 0.0;
      T_over_m_quad[n] = 0.0;

      // Compute the configuration-space quadrature
      for (int k=0; k<num_conf_basis; ++k) {
        n_quad[n] += n_d[k]*b_ord[k];
        T_over_m_quad[n] += T_over_m_d[k]*b_ord[k];
      }

      if ((n_quad[n] > 0.0) && (T_over_m_quad[n] > 0.0)) {
        expamp_quad[n] = n_quad[n]/sqrt(pow(2.0*GKYL_PI*T_over_m_quad[n], vdim));
      }
      else {
        expamp_quad[n] = 0.0;
      }
    }

    // inner loop over velocity space
    gkyl_range_deflate(&vel_rng, phase_range, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
    while (gkyl_range_iter_next(&vel_iter)) {
      
      copy_idx_arrays(conf_range->ndim, phase_range->ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(&up->phase_grid, pidx, xc);

      struct gkyl_range_iter qiter;
      // compute LTE distribution function at phase-space quadrature nodes
      gkyl_range_iter_init(&qiter, &up->phase_qrange);
      while (gkyl_range_iter_next(&qiter)) {

        int cqidx = gkyl_range_idx(&up->conf_qrange, qiter.idx);
        int pqidx = gkyl_range_idx(&up->phase_qrange, qiter.idx);

        const double *xcomp_d = gkyl_array_cfetch(up->ordinates, pqidx);

        comp_to_phys(pdim, xcomp_d, up->phase_grid.dx, xc, xmu);
  
        double *f0_q = gkyl_array_fetch(up->f0_at_ords, pqidx);
        double *g_q = gkyl_array_fetch(up->g_at_ords, pqidx);
        f0_q[0] = f_floor;
        g_q[0] = f_floor;
        
        if (T_over_m_quad[cqidx] > 0.0) {
          double v_sq = 0.0;
          for (int d=0; d<vdim; ++d) {
            v_sq += pow(xmu[cdim+d], 2.0);
          }

          double f0_tmp = expamp_quad[cqidx]*exp(-v_sq/(2.0*T_over_m_quad[cqidx])); 
          f0_q[0] += f0_tmp;
          g_q[0] += T_over_m_quad[cqidx]*f0_tmp;
        }
      }
      // compute expansion coefficients of LTE distribution functions on basis
      long lidx = gkyl_range_idx(&vel_rng, vel_iter.idx);
      double *f_lte_d = (double *)gkyl_array_fetch(f_lte, lidx);
      proj_on_basis(up, up->f0_at_ords, f_lte_d);
      proj_on_basis(up, up->g_at_ords, &f_lte_d[num_phase_basis]);
    }
  }
}

void
gkyl_flash_proj_pkpm_lte_on_basis_release(gkyl_flash_proj_pkpm_lte_on_basis *up)
{
  gkyl_array_release(up->ordinates);
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);
  gkyl_array_release(up->conf_basis_at_ords);
  gkyl_array_release(up->f0_at_ords);
  gkyl_array_release(up->g_at_ords);

  gkyl_free(up);
}