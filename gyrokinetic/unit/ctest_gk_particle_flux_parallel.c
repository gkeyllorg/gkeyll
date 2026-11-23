#include <acutest.h>
#include <gkyl_array_ops.h>
#include <gkyl_gk_particle_flux_parallel_priv.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_gk_collisionless_flux_priv.h>

// Affine nonidentity map with a zero crossing inside a velocity cell.
static void
velocity_map(double t, const double *x, double *v, void *ctx)
{
  v[0] = 2.0 * x[0] + 0.2;
  if (*(int *)ctx == 2) {
    v[1] = 1.5 * x[1];
  }
}

// A polynomial with variation in every configuration and velocity direction.
static void
trace_profile(double t, const double *x, double *out, void *ctx)
{
  int *dims = ctx, cdim = dims[0], vdim = dims[1];
  out[0] = 1.0;
  for (int d = 0; d < cdim; ++d) {
    out[0] *= 1.0 + 0.1 * (d + 1) * x[d];
  }
  out[0] *= 1.0 + 0.15 * x[cdim] + 0.1 * x[cdim] * x[cdim];
  if (vdim == 2) {
    out[0] *= 1.0 + 0.2 * x[cdim + 1];
  }
}

// Check constant skin distributions with nonzero drift geometry, unequal end
// areas, both field orientations, and poisoned distribution ghost cells.
static void
check_flux(int cdim, int vdim, bool gpu, bool mapped)
{
  struct gkyl_basis cb, pb;
  gkyl_cart_modal_serendip(&cb, cdim, 1);
  gkyl_cart_modal_gkhybrid(&pb, cdim, vdim);
  double lo[5] = {-1, -1, -1, -1, -1}, hi[5] = {1, 1, 1, 1, 1};
  int cells[5] = {2, 2, 2, 2, 2}, ghosts[5] = {0};
  for (int d = 0; d < cdim; ++d) {
    ghosts[d] = 1;
  }
  if (vdim == 2) {
    lo[cdim + 1] = 0.2;
    hi[cdim + 1] = 1.8;
  }
  if (mapped) {
    cells[cdim] = 3;
  }
  struct gkyl_rect_grid grid, cgrid, vgrid;
  gkyl_rect_grid_init(&grid, cdim + vdim, lo, hi, cells);
  gkyl_rect_grid_init(&cgrid, cdim, lo, hi, cells);
  gkyl_rect_grid_init(&vgrid, vdim, lo + cdim, hi + cdim, cells + cdim);
  struct gkyl_range phase, phase_ext, conf, conf_ext, vel, vel_ext;
  gkyl_create_grid_ranges(&grid, ghosts, &phase_ext, &phase);
  gkyl_create_grid_ranges(&cgrid, ghosts, &conf_ext, &conf);
  gkyl_create_grid_ranges(&vgrid, (int[]){0, 0}, &vel_ext, &vel);
  struct gkyl_gk_geometry_inp gi = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .grid = cgrid,
    .basis = cb,
    .local = conf,
    .local_ext = conf_ext,
    .global = conf,
    .global_ext = conf_ext,
  };
  struct gk_geometry *geom = gkyl_gk_geometry_new(NULL, &gi, gpu);
  struct gkyl_array *bmag = gkyl_array_new(GKYL_DOUBLE, cb.num_basis, conf_ext.volume);
  gkyl_array_clear(bmag, 0.0);
  gkyl_array_shiftc(bmag, 2.0 * pow(sqrt(2.0), cdim), 0);
  gkyl_array_copy(geom->geo_corn.bmag, bmag);
  struct gkyl_velocity_map *vm = gkyl_velocity_map_new(
    (struct gkyl_mapc2p_inp){.mapping = mapped ? velocity_map : NULL, .ctx = &vdim}, grid, vgrid,
    phase, phase_ext, vel, vel_ext, gpu
  );
  struct gkyl_dg_geom_inp di = {
    .grid = &cgrid,
    .range = &conf_ext,
    .skip_geometry_creation = true,
    .nquad = 2,
  };
  struct gkyl_gk_dg_geom_inp ki = {
    .grid = &cgrid,
    .range = &conf_ext,
    .skip_geometry_creation = true,
    .nquad = 2,
  };
  struct gkyl_dg_geom *dg = gkyl_dg_geom_new(&di);
  struct gkyl_gk_dg_geom *gk = gkyl_gk_dg_geom_new(&ki);
  int nq = 1 << (cdim - 1);
  struct gkyl_array *fh = gkyl_array_new(GKYL_DOUBLE, pb.num_basis, phase_ext.volume);
  gkyl_array_clear(fh, 731.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &phase);
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(fh, gkyl_range_idx(&phase, iter.idx));
    for (int k = 0; k < pb.num_basis; ++k) {
      f[k] = 0.0;
    }
    f[0] = pow(sqrt(2.0), cdim + vdim);
  }
  int dims[] = {cdim, vdim};
  if (mapped) {
    struct gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &pb, 3, 1, trace_profile, dims);
    gkyl_proj_on_basis_advance(proj, 0.0, &phase, fh);
    gkyl_proj_on_basis_release(proj);
  }
  struct gkyl_array *fin = gpu ?
                             gkyl_array_cu_dev_new(GKYL_DOUBLE, pb.num_basis, phase_ext.volume) :
                             gkyl_array_acquire(fh);
  gkyl_array_copy(fin, fh);
  for (int orientation = -1; orientation <= 1; orientation += 2) {
    for (long cell = 0; cell < conf_ext.volume; ++cell) {
      struct gkyl_dg_surf_geom *ds = gkyl_array_fetch(dg->surf_geom[cdim - 1], cell);
      struct gkyl_gk_dg_surf_geom *gs = gkyl_array_fetch(gk->surf_geom[cdim - 1], cell);
      for (int q = 0; q < nq; ++q) {
        ds[q].area_elem = 3.0 + 0.1 * q;
        gs[q] = (struct gkyl_gk_dg_surf_geom){
          .B3 = orientation * 1.0,
          .bmag = 2.0,
          .Jc = 1.0,
          .normcurlbhat = 9.0,
          .bhat = {.x = {2.0, 3.0, 4.0}},
        };
      }
    }
    struct gkyl_dg_geom *dg_run = gkyl_dg_geom_new_from_host(&di, dg, gpu);
    struct gkyl_gk_dg_geom *gk_run = gkyl_gk_dg_geom_new_from_host(&ki, gk, gpu);
    for (int edge = 0; edge < 2; ++edge) {
      enum gkyl_edge_loc loc = edge ? GKYL_UPPER_EDGE : GKYL_LOWER_EDGE;
      struct gkyl_range cs, cg, ps, pg;
      gkyl_skin_ghost_ranges(&cs, &cg, cdim - 1, loc, &conf_ext, ghosts);
      gkyl_skin_ghost_ranges(&ps, &pg, cdim - 1, loc, &phase_ext, ghosts);
      // Exclude transverse ghosts from the physical boundary.
      int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
      gkyl_copy_int_arr(cdim + vdim, phase.lower, lower);
      gkyl_copy_int_arr(cdim + vdim, phase.upper, upper);
      lower[cdim - 1] = upper[cdim - 1] = edge ? conf.upper[cdim - 1] : conf.lower[cdim - 1];
      gkyl_sub_range_init(&ps, &phase_ext, lower, upper);
      gkyl_sub_range_init(&cs, &conf_ext, lower, upper);
      struct gkyl_gk_particle_flux_parallel_inp inp = {
        .edge = loc,
        .cbasis = &cb,
        .pbasis = &pb,
        .phase_grid = &grid,
        .conf_skin = &cs,
        .phase_skin = &ps,
        .charge = 1,
        .mass = 1,
        .gk_geom = geom,
        .dg_geom = dg_run,
        .gk_dg_geom = gk_run,
        .vel_map = vm,
        .use_gpu = gpu,
      };
      struct gkyl_gk_particle_flux_parallel *up = gkyl_gk_particle_flux_parallel_inew(&inp);
      struct gkyl_array *out = gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nq, cs.volume) :
                                     gkyl_array_new(GKYL_DOUBLE, nq, cs.volume);
      struct gkyl_array *oh = gkyl_array_new(GKYL_DOUBLE, nq, cs.volume);
      gkyl_gk_particle_flux_parallel_advance(up, fin, out);
      gkyl_array_copy(oh, out);
      if (!gpu && cdim > 1) {
        // Compare the reused general kernel with the existing parallel no_by
        // kernel using the original, unmasked drift geometry.
        struct gkyl_gk_particle_flux_parallel ref = *up;
        ref.kernel = edge ? choose_gk_collisionless_flux_no_by_edge_surf_conf_kern(
                              cdim - 1, cdim, vdim, 1, GKYL_BC_GK_SPECIES_ABSORB
                            ) :
                            choose_gk_collisionless_flux_no_by_surf_conf_kern(
                              cdim - 1, cdim, vdim, 1, GKYL_BC_GK_SPECIES_ABSORB
                            );
        ref.surf_geo_fields =
          gkyl_array_new(GKYL_USER, nq * sizeof(struct gkyl_gk_dg_surf_geom), cs.volume);
        for (long cell = 0; cell < cs.volume; ++cell) {
          int idx[GKYL_MAX_DIM];
          gkyl_range_inv_idx(&ref.compact, cell, idx);
          if (edge) {
            idx[cdim - 1]++;
          }
          memcpy(
            gkyl_array_fetch(ref.surf_geo_fields, cell),
            gkyl_gk_dg_geom_get_surf(gk, cdim - 1, idx), nq * sizeof(struct gkyl_gk_dg_surf_geom)
          );
        }
        struct gkyl_array *expected = gkyl_array_new(GKYL_DOUBLE, nq, cs.volume);
        gkyl_array_clear(expected, 0.0);
        struct gkyl_range_iter pit;
        gkyl_range_iter_init(&pit, &ps);
        while (gkyl_range_iter_next(&pit)) {
          double val[4];
          particle_flux_parallel_cell(&ref, pit.idx, fh, val);
          double *sum = gkyl_array_fetch(expected, gkyl_range_idx(&ref.compact, pit.idx));
          for (int q = 0; q < nq; ++q) {
            sum[q] += val[q];
          }
        }
        for (long cell = 0; cell < cs.volume; ++cell) {
          const double *a = gkyl_array_cfetch(oh, cell), *b = gkyl_array_cfetch(expected, cell);
          for (int q = 0; q < nq; ++q) {
            TEST_CHECK(gkyl_compare(a[q], b[q], 1e-12));
          }
        }
        gkyl_array_release(expected);
        gkyl_array_release(ref.surf_geo_fields);
      }

      struct gkyl_range compact;
      gkyl_range_init(&compact, cdim, cs.lower, cs.upper);
      double nodes[8], qw[4];
      if (cdim > 1) {
        gkyl_ndim_ordinates_weights(cdim - 1, nodes, qw, 2);
      }
      for (long cell = 0; cell < cs.volume; ++cell) {
        const double *val = gkyl_array_cfetch(oh, cell);
        for (int q = 0; q < nq; ++q) {
          double expected = 0.5 * (3.0 + 0.1 * q) / 2.0 * (vdim == 2 ? 1.6 * 2.0 * M_PI : 1.0);
          if (mapped) {
            int idx[GKYL_MAX_DIM];
            gkyl_range_inv_idx(&compact, cell, idx);
            double x[GKYL_MAX_DIM] = {0.0};
            for (int d = 0; d < cdim - 1; ++d) {
              x[d] =
                lo[d] + (idx[d] - 0.5) * grid.dx[d] + 0.5 * grid.dx[d] * nodes[q * (cdim - 1) + d];
            }
            x[cdim - 1] = edge ? hi[cdim - 1] : lo[cdim - 1];
            double sum = 0.0;
            for (int iv = 0; iv < cells[cdim]; ++iv) {
              for (int j = 0; j < 3; ++j) {
                x[cdim] = lo[cdim] + (iv + 0.5) * grid.dx[cdim] +
                          0.5 * grid.dx[cdim] * gkyl_gauss_ordinates[3][j];
                double velocity = (edge ? 1.0 : -1.0) * orientation * (2.0 * x[cdim] + 0.2);
                if (velocity <= 0.0) {
                  continue;
                }
                int nm = vdim == 2 ? cells[cdim + 1] : 1;
                for (int im = 0; im < nm; ++im) {
                  x[cdim + 1] = vdim == 2 ? lo[cdim + 1] + (im + 0.5) * grid.dx[cdim + 1] : 0.0;
                  double f;
                  trace_profile(0.0, x, &f, dims);
                  // The input is Jf; the affine velocity Jacobian is already in it.
                  sum += gkyl_gauss_weights[3][j] * grid.dx[cdim] / 2.0 * velocity * f *
                         (vdim == 2 ? grid.dx[cdim + 1] * 2.0 * M_PI : 1.0);
                }
              }
            }
            expected = (3.0 + 0.1 * q) / 2.0 * sum;
          }
          TEST_CHECK(gkyl_compare(val[q], expected, 1e-12));
          TEST_MSG(
            "cdim=%d vdim=%d edge=%d orientation=%d q=%d: %.16g != %.16g", cdim, vdim, edge,
            orientation, q, val[q], expected
          );
        }
      }
      gkyl_array_release(oh);
      gkyl_array_release(out);
      gkyl_gk_particle_flux_parallel_release(up);
    }
    if (gpu) {
      gkyl_dg_geom_release(dg_run);
      gkyl_gk_dg_geom_release(gk_run);
    }
  }
  gkyl_array_release(fin);
  gkyl_array_release(fh);
  gkyl_array_release(bmag);
  gkyl_dg_geom_release(dg);
  gkyl_gk_dg_geom_release(gk);
  gkyl_velocity_map_release(vm);
  gkyl_gk_geometry_release(geom);
}

// Exercise each supported phase-space dimensionality on the CPU.
static void
test_particle_flux_parallel(void)
{
  for (int mapped = 0; mapped < 2; ++mapped) {
    check_flux(1, 1, false, mapped);
    check_flux(1, 2, false, mapped);
    check_flux(2, 2, false, mapped);
    check_flux(3, 2, false, mapped);
  }
}
#ifdef GKYL_HAVE_CUDA
// Exercise the same numerical checks on the device.
static void
test_particle_flux_parallel_dev(void)
{
  for (int mapped = 0; mapped < 2; ++mapped) {
    check_flux(1, 1, true, mapped);
    check_flux(1, 2, true, mapped);
    check_flux(2, 2, true, mapped);
    check_flux(3, 2, true, mapped);
  }
}
#endif
TEST_LIST = {
  {"particle_flux_parallel", test_particle_flux_parallel},
#ifdef GKYL_HAVE_CUDA
  {"particle_flux_parallel_dev", test_particle_flux_parallel_dev},
#endif
  {NULL, NULL}
};
