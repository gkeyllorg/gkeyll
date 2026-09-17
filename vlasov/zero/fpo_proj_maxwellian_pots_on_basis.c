#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis_priv.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_mat.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];
  for (int i=cdim; i<pdim; ++i)
    out[i] = vidx[i-cdim];
}

struct gkyl_proj_maxwellian_pots_on_basis*
gkyl_proj_maxwellian_pots_on_basis_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, int num_quad, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_proj_maxwellian_pots_on_basis_cu_ho_new(grid, conf_range,
      phase_range, conf_basis, phase_basis, num_quad);
  }
#endif

  gkyl_proj_maxwellian_pots_on_basis *up = gkyl_malloc(sizeof(gkyl_proj_maxwellian_pots_on_basis));
  up->grid = *grid;
  up->cdim = conf_basis->ndim;
  up->pdim = phase_basis->ndim;
  int vdim = up->pdim - up->cdim;
  up->num_quad = num_quad;

  up->conf_basis = *conf_basis;

  up->use_gpu = use_gpu;

  if (phase_basis->poly_order == 1) {
    gkyl_cart_modal_hybrid(&up->surf_basis, up->cdim, vdim-1);
  } else {
    gkyl_cart_modal_serendip(&up->surf_basis, up->pdim-1, phase_basis->poly_order);
  }

  up->num_conf_basis = conf_basis->num_basis;
  up->num_phase_basis = phase_basis->num_basis;
  up->num_surf_basis = up->surf_basis.num_basis;
 
  // Quadrature points for phase space expansion
  up->tot_quad = init_quad_values(up->cdim, phase_basis, num_quad,
    &up->ordinates, &up->weights, &up->basis_at_ords, false);

  // Initialize quadrature for surface expansion
  up->tot_surf_quad = init_quad_values(up->cdim, &up->surf_basis, num_quad,
    &up->surf_ordinates, &up->surf_weights, &up->surf_basis_at_ords, false);

  // Hybrid basis support: uses p=2 in velocity space
  int num_quad_v = num_quad;
  bool is_vdim_p2[] = {false, false, false};  // 3 is the max vdim.
  if (phase_basis->b_type == GKYL_BASIS_MODAL_HYBRID) {
    num_quad_v = num_quad+1;
    // Maxwellian potentials are always 3V
    is_vdim_p2[0] = true, is_vdim_p2[1] = true, is_vdim_p2[2] = true;
  }

  up->phase_qrange = get_qrange(up->cdim, up->pdim, num_quad, num_quad_v, is_vdim_p2);
  up->surf_qrange = get_qrange(up->cdim, up->pdim-1, num_quad, num_quad_v, is_vdim_p2);

  // Nodes for nodal surface expansions
  struct gkyl_array *surf_nodes_ho = gkyl_array_new(GKYL_DOUBLE, grid->ndim-1, up->surf_basis.num_basis);
  up->surf_basis.node_list(gkyl_array_fetch(surf_nodes_ho, 0));

  up->fpo_h_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);
  up->fpo_g_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);
  up->fpo_dhdv_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);
  up->fpo_d2gdv2_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);

  up->fpo_h_at_surf_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_surf_quad);
  up->fpo_g_at_surf_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_surf_quad);
  up->fpo_dhdv_at_surf_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_surf_quad);
  up->fpo_dgdv_at_surf_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_surf_quad);
  up->fpo_d2gdv2_at_surf_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_surf_quad);

  return up;
}

void
gkyl_proj_maxwellian_pots_on_basis_advance(const gkyl_proj_maxwellian_pots_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array* prim_moms,
  struct gkyl_array *fpo_h, struct gkyl_array *fpo_g,
  struct gkyl_array *fpo_h_surf, struct gkyl_array *fpo_g_surf,
  struct gkyl_array *fpo_dhdv_surf, struct gkyl_array *fpo_dgdv_surf,
  struct gkyl_array *fpo_d2gdv2_surf)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_proj_maxwellian_pots_on_basis_advance_cu(up, phase_range, conf_range, prim_moms,
      fpo_h, fpo_g, fpo_h_surf, fpo_g_surf, fpo_dhdv_surf, fpo_dgdv_surf, fpo_d2gdv2_surf);
  }
#endif

  // Calculate Maxwellian potentials using primitive moments
  int cdim = up->cdim, pdim = up->pdim;
  int vdim = pdim-cdim;
  int tot_quad = up->tot_quad;
  int num_phase_basis = up->num_phase_basis;
  int num_conf_basis = up->num_conf_basis;
  int num_surf_basis = up->num_surf_basis;

  struct gkyl_range_iter phase_iter, qiter, surf_qiter;

  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<cdim; ++d) rem_dir[d] = 1;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  // Loop over configuration space cells
  gkyl_range_iter_init(&phase_iter, phase_range);
  while (gkyl_range_iter_next(&phase_iter)) {
    long linc = gkyl_range_idx(conf_range, phase_iter.idx);
    long linp = gkyl_range_idx(phase_range, phase_iter.idx);
    gkyl_copy_int_arr(pdim, phase_iter.idx, pidx);
    gkyl_rect_grid_cell_center(&up->grid, pidx, xc);

    const double *prim_moms_d = gkyl_array_cfetch(prim_moms, linc);
    const double *m0_d = &prim_moms_d[0];
    const double *u_drift_d = &prim_moms_d[num_conf_basis];
    const double *vtsq_d = &prim_moms_d[num_conf_basis*(vdim+1)];

    // Iterate over quadrature points and compute compute potentials H and G
    gkyl_range_iter_init(&qiter, &up->phase_qrange);
    while (gkyl_range_iter_next(&qiter)) {
      long linq = gkyl_range_idx(&up->phase_qrange, qiter.idx);
      const double* quad_node = gkyl_array_cfetch(up->ordinates, linq);

      comp_to_phys(pdim, quad_node, up->grid.dx, xc, xmu);

      // Evaluate moments at quadrature node
      double den_q = up->conf_basis.eval_expand(quad_node, m0_d);
      double vtsq_q = up->conf_basis.eval_expand(quad_node, vtsq_d);

      double rel_speedsq_q = 0.0;
      for (int d=0; d<vdim; ++d) {
        double udrift_q = up->conf_basis.eval_expand(quad_node, &u_drift_d[d*num_conf_basis]);
        rel_speedsq_q += pow(xmu[cdim+d]-udrift_q,2);
      }
      double rel_speed_q = sqrt(rel_speedsq_q);

      double* fpo_h_q = gkyl_array_fetch(up->fpo_h_at_ords, linq);
      double* fpo_g_q = gkyl_array_fetch(up->fpo_g_at_ords, linq);

      fpo_h_q[0] = eval_fpo_h(den_q, rel_speed_q, vtsq_q);
      fpo_g_q[0] = eval_fpo_g(den_q, rel_speed_q, vtsq_q);
    }

    // Project potentials onto basis
    proj_on_basis(up, up->fpo_h_at_ords, gkyl_array_fetch(fpo_h, linp));
    proj_on_basis(up, up->fpo_g_at_ords, gkyl_array_fetch(fpo_g, linp));

    // Iterate over velocity directions for surface expansions at domain boundaries
    for (int d1=0; d1<vdim; ++d1) {
      int dir1 = d1 + cdim;
      bool is_edge_in_dir1 = 0;
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range->lower[dir1]);
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range->upper[dir1]);

      if (is_edge_in_dir1) {
        // Compute H, dH/dv, G, d2G/dv2 at dir1 boundaries.
        // Velocity value at dir1 boundary
        double vbound = pidx[dir1] == phase_range->lower[dir1] ?
          up->grid.lower[dir1] : up->grid.upper[dir1];

        // Calculate surface expansions using Gauss-Legendre quadrature,
        // similar to proj_on_basis routine.
        gkyl_range_iter_init(&surf_qiter, &up->surf_qrange);
        while (gkyl_range_iter_next(&surf_qiter)) {
          int surf_qidx = gkyl_range_idx(&up->surf_qrange, surf_qiter.idx);
          const double* surf_ord = gkyl_array_cfetch(up->surf_ordinates, surf_qidx);
          surf_comp_to_phys(dir1, pdim, surf_ord, up->grid.dx, xc, xmu);
          xmu[dir1] = vbound;

          // Have to map pdim-1 surface quadrature index to pdim phase quadrature index
          // to get correct phase space variables.
          int edge_idx = pidx[dir1] == phase_range->lower[dir1] ? 0 : up->num_quad-1;
          int phase_idx[GKYL_MAX_DIM];
          edge_idx_to_phase_idx(pdim, dir1, surf_qiter.idx, edge_idx, phase_idx);
          int phase_qidx = gkyl_range_idx(&up->phase_qrange, phase_idx);
          comp_to_phys(pdim, gkyl_array_cfetch(up->ordinates, phase_qidx), up->grid.dx, xc, xmu);
          xmu[dir1] = vbound;

          double den_q = up->conf_basis.eval_expand(surf_ord, m0_d);
          double vtsq_q = up->conf_basis.eval_expand(surf_ord, vtsq_d);

          double rel_speedsq_q = 0.0;
          double rel_vel_in_dir_q = 0.0;
          for (int d=0; d<vdim; ++d) {
            double udrift_q = up->conf_basis.eval_expand(surf_ord, &u_drift_d[d*num_conf_basis]);
              rel_speedsq_q += pow(xmu[cdim+d]-udrift_q, 2);
              if (d == d1)
                rel_vel_in_dir_q += xmu[dir1]-udrift_q;
          }
          double rel_speed_q = sqrt(rel_speedsq_q);

          double* fpo_h_at_surf_ords_q = gkyl_array_fetch(up->fpo_h_at_surf_ords, surf_qidx);
          double* fpo_g_at_surf_ords_q = gkyl_array_fetch(up->fpo_g_at_surf_ords, surf_qidx);
          double* fpo_dhdv_at_surf_ords_q = gkyl_array_fetch(up->fpo_dhdv_at_surf_ords, surf_qidx);
          double* fpo_d2gdv2_at_surf_ords_q = gkyl_array_fetch(up->fpo_d2gdv2_at_surf_ords, surf_qidx);

          fpo_h_at_surf_ords_q[0] = eval_fpo_h(den_q, rel_speed_q, vtsq_q);
          fpo_g_at_surf_ords_q[0] = eval_fpo_g(den_q, rel_speed_q, vtsq_q);
          fpo_dhdv_at_surf_ords_q[0] = eval_fpo_dhdv(den_q, rel_vel_in_dir_q, vtsq_q, rel_speed_q);
          fpo_d2gdv2_at_surf_ords_q[0] = eval_fpo_d2gdv2(den_q, rel_vel_in_dir_q, vtsq_q, rel_speed_q);
        }

        proj_on_surf_basis(up, d1, up->fpo_h_at_surf_ords, gkyl_array_fetch(fpo_h_surf, linp));
        proj_on_surf_basis(up, d1, up->fpo_g_at_surf_ords, gkyl_array_fetch(fpo_g_surf, linp));
        proj_on_surf_basis(up, d1, up->fpo_dhdv_at_surf_ords, gkyl_array_fetch(fpo_dhdv_surf, linp));
        proj_on_surf_basis(up, d1, up->fpo_d2gdv2_at_surf_ords, gkyl_array_fetch(fpo_d2gdv2_surf, linp));

        // Iterate over transverse directions at dir1 boundary for dG/dv
        for (int d2=0; d2<vdim; ++d2) {
          if (d1 == d2) continue;
          int dir2 = d2+cdim;

          // Calculate surface expansions using Gauss-Legendre quadrature,
          // similar to proj_on_basis routine.
          gkyl_range_iter_init(&surf_qiter, &up->surf_qrange);
          while (gkyl_range_iter_next(&surf_qiter)) {
            int surf_qidx = gkyl_range_idx(&up->surf_qrange, surf_qiter.idx);
            const double* surf_ord = gkyl_array_cfetch(up->surf_ordinates, surf_qidx);
            surf_comp_to_phys(dir1, pdim, surf_ord, up->grid.dx, xc, xmu);
            xmu[dir1] = vbound;

            double den_q = up->conf_basis.eval_expand(surf_ord, m0_d);
            double vtsq_q = up->conf_basis.eval_expand(surf_ord, vtsq_d);

            double rel_speedsq_q = 0.0;
            double rel_vel_in_dir2_q = 0.0;
            for (int d=0; d<vdim; ++d) {
              double udrift_q = up->conf_basis.eval_expand(surf_ord, &u_drift_d[d*num_conf_basis]);
                rel_speedsq_q += pow(xmu[cdim+d]-udrift_q, 2);
                if (d == d2)
                  rel_vel_in_dir2_q += xmu[dir2]-udrift_q;
            }
            double rel_speed_q = sqrt(rel_speedsq_q);

            // Compute dg/dv2 at v1 boundary
            double *fpo_dgdv_at_surf_ords_q = gkyl_array_fetch(up->fpo_dgdv_at_surf_ords, surf_qidx);
            fpo_dgdv_at_surf_ords_q[0] = eval_fpo_dgdv(den_q, rel_vel_in_dir2_q,
              vtsq_q, rel_speedsq_q);
          }
          int dgdv_index = (d1 < d2) ? d2-1 : d2;
          dgdv_index += (vdim-1)*d1;
          proj_on_surf_basis(up, dgdv_index, up->fpo_dgdv_at_surf_ords, gkyl_array_fetch(fpo_dgdv_surf, linp));
        }
      }
    }
  }
}

void
gkyl_proj_maxwellian_pots_deriv_on_basis_advance(const gkyl_proj_maxwellian_pots_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *prim_moms,
  struct gkyl_array *fpo_dhdv, struct gkyl_array *fpo_d2gdv2)
{
  // Calculate Maxwellian potentials using primitive moments
  int cdim = up->cdim, pdim = up->pdim;
  int vdim = pdim-cdim;
  int tot_quad = up->tot_quad;
  int num_phase_basis = up->num_phase_basis;
  int num_conf_basis = up->num_conf_basis;
  int num_surf_basis = up->num_surf_basis;

  struct gkyl_range vel_range;
  struct gkyl_range_iter conf_iter, vel_iter;

  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_range->ndim; ++d) rem_dir[d] = 1;

  // struct gkyl_array *fpo_dgdv_at_surf_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_surf_basis);

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  // Loop over configuration space cells for quad integration of moments
  gkyl_range_iter_init(&conf_iter, conf_range);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_range, conf_iter.idx);

    const double *prim_moms_d = gkyl_array_cfetch(prim_moms, midx);
    const double *m0_d = &prim_moms_d[0];
    const double *u_drift_d = &prim_moms_d[num_conf_basis];
    const double *vtsq_d = &prim_moms_d[num_conf_basis*(vdim+1)];

    // Inner loop over velocity space
    gkyl_range_deflate(&vel_range, phase_range, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_range);
    while (gkyl_range_iter_next(&vel_iter)) {
      copy_idx_arrays(conf_range->ndim, phase_range->ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(&up->grid, pidx, xc);
 
      long lidx = gkyl_range_idx(&vel_range, vel_iter.idx);

      // Project potential onto basis
      proj_on_basis(up, up->fpo_dhdv_at_ords, gkyl_array_fetch(fpo_dhdv, lidx));

      // Iterate over directions to calculate derivatives of H and G
      for (int d1=0; d1<vdim; ++d1) {
        for (int d2=0; d2<vdim; ++d2) {
          struct gkyl_range_iter qiter;
          gkyl_range_iter_init(&qiter, &up->phase_qrange);
 
          // Iterate over quadrature points and compute compute potentials
          while (gkyl_range_iter_next(&qiter)) {
            int pqidx = gkyl_range_idx(&up->phase_qrange, qiter.idx);
            const double* quad_node = gkyl_array_cfetch(up->ordinates, pqidx);
 
            comp_to_phys(pdim, quad_node, up->grid.dx, xc, xmu);
 
            // Evaluate moments at quadrature node
            double den_q = up->conf_basis.eval_expand(quad_node, m0_d);
            double u_drift_q = up->conf_basis.eval_expand(quad_node, u_drift_d);
            double vtsq_q = up->conf_basis.eval_expand(quad_node, vtsq_d);
 
            double rel_speedsq_q = 0.0;
            double rel_vel_in_dir1_q = 0.0;
            double rel_vel_in_dir2_q = 0.0;
            for (int d=0; d<vdim; ++d) {
              double udrift_q = up->conf_basis.eval_expand(quad_node, &u_drift_d[d*num_conf_basis]);
              rel_speedsq_q += pow(xmu[cdim+d]-udrift_q,2);
              if (d == d1)
                rel_vel_in_dir1_q += xmu[d1+cdim]-udrift_q;
              if (d == d2)
                rel_vel_in_dir2_q += xmu[d2+cdim]-udrift_q;
            }
            double rel_speed_q = sqrt(rel_speedsq_q);
 
            double* fpo_dhdv_q = gkyl_array_fetch(up->fpo_dhdv_at_ords, pqidx);
            double* fpo_d2gdv2_q = gkyl_array_fetch(up->fpo_d2gdv2_at_ords, pqidx);

            if (d1 == d2) {
              fpo_dhdv_q[0] = eval_fpo_dhdv(den_q, rel_vel_in_dir1_q, vtsq_q, rel_speed_q);

              fpo_d2gdv2_q[0] = eval_fpo_d2gdv2(den_q,
                rel_vel_in_dir1_q, vtsq_q, rel_speed_q);
            }
            else {
              fpo_d2gdv2_q[0] = eval_fpo_d2gdv2_cross(den_q,
                rel_vel_in_dir1_q, rel_vel_in_dir2_q, rel_speed_q, vtsq_q);
            }
          }
          double *fpo_d2gdv2_d = gkyl_array_fetch(fpo_d2gdv2, lidx);
          proj_on_basis(up, up->fpo_d2gdv2_at_ords, &fpo_d2gdv2_d[(d1*vdim+d2)*up->num_phase_basis]);

          if (d1 == d2) {
            double *fpo_dhdv_d = gkyl_array_fetch(fpo_dhdv, lidx);
            proj_on_basis(up, up->fpo_dhdv_at_ords, &fpo_dhdv_d[d1*up->num_phase_basis]);
          }
        }
      }
    }
  }
}

void
gkyl_proj_maxwellian_pots_on_basis_release(gkyl_proj_maxwellian_pots_on_basis *up)
{
  gkyl_array_release(up->ordinates);
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);

  gkyl_array_release(up->surf_ordinates);
  gkyl_array_release(up->surf_weights);
  gkyl_array_release(up->surf_basis_at_ords);

  if (up->use_gpu) {
    gkyl_array_release(up->conf_basis_at_nodes);

    gkyl_array_release(up->prim_moms_conf_quad);
    gkyl_array_release(up->prim_moms_conf_nodes);
 
    gkyl_array_release(up->pot_phase_quad);
    gkyl_array_release(up->pot_surf_quad);
    gkyl_array_release(up->pot_deriv_surf_quad);
    gkyl_array_release(up->sol_pot_surf_modal);
    gkyl_array_release(up->fpo_dgdv_at_surf_ords);
 
    gkyl_mat_mm_array_mem_release(up->phase_quad_nodal_to_modal_mem);
    gkyl_mat_mm_array_mem_release(up->surf_quad_nodal_to_modal_mem);
    gkyl_mat_mm_array_mem_release(up->dgdv_surf_quad_nodal_to_modal_mem);

    gkyl_cart_modal_basis_release_cu(up->surf_basis_dev);

    gkyl_array_release(up->conf_ordinates);
    gkyl_array_release(up->conf_weights);
    gkyl_array_release(up->conf_basis_at_ords);
  }
  else {
    gkyl_array_release(up->fpo_h_at_ords);
    gkyl_array_release(up->fpo_g_at_ords);
    gkyl_array_release(up->fpo_dhdv_at_ords);
    gkyl_array_release(up->fpo_d2gdv2_at_ords);

    gkyl_array_release(up->fpo_h_at_surf_ords);
    gkyl_array_release(up->fpo_g_at_surf_ords);
    gkyl_array_release(up->fpo_dhdv_at_surf_ords);
    gkyl_array_release(up->fpo_dgdv_at_surf_ords);
    gkyl_array_release(up->fpo_d2gdv2_at_surf_ords);
  }

  gkyl_free(up);
}

