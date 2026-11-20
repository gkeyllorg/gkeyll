#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_mom_weighted_gyrokinetic.h>
#include <gkyl_rect_decomp.h>

/**
 * Allocate a double array on the selected backend.
 *
 * @param gpu Whether to use device storage.
 * @param ncomp Components per cell.
 * @param size Number of cells.
 * @return Array owned by the caller.
 */
static struct gkyl_array *
mkarr(bool gpu, int ncomp, long size)
{
  return gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp, size) :
               gkyl_array_new(GKYL_DOUBLE, ncomp, size);
}

/**
 * Evaluate a modal expansion at a reference-cell coordinate.
 *
 * @param basis Modal basis.
 * @param coord Reference coordinate.
 * @param coeff Modal coefficients.
 * @return Value of the represented polynomial.
 */
static double
eval(const struct gkyl_basis *basis, const double *coord, const double *coeff)
{
  double values[112], result = 0.0;
  basis->eval(coord, values);
  for (int k = 0; k < basis->num_basis; ++k) {
    result += values[k] * coeff[k];
  }
  return result;
}

/**
 * Supply a nonlinear map so moments also exercise projected velocity coordinates.
 *
 * @param time Unused time.
 * @param coord Computational velocity coordinates.
 * @param mapped Output physical velocity coordinates.
 * @param ctx Pointer to the velocity dimension.
 */
static void
map_velocity(double time, const double *coord, double *mapped, void *ctx)
{
  int vdim = *(int *)ctx;
  mapped[0] = coord[0] + 0.1 * coord[0] * coord[0];
  if (vdim == 2) {
    mapped[1] = coord[1] + 0.05 * coord[1] * coord[1];
  }
}

/**
 * Check kinetic weighted moments against independent tensor Gauss quadrature.
 *
 * @param cdim Configuration dimension.
 * @param vdim Velocity dimension.
 * @param order Polynomial order.
 * @param gpu Whether to exercise CUDA updaters.
 * @param mapped Whether to use nonlinear velocity coordinates.
 */
static void
check_weighted(int cdim, int vdim, int order, bool gpu, bool mapped)
{
  struct gkyl_basis cbasis, pbasis, vbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, order);

  if (order == 1) {
    gkyl_cart_modal_gkhybrid(&pbasis, cdim, vdim);
    gkyl_cart_modal_gkhybrid_vel(&vbasis, vdim);
  } else {
    gkyl_cart_modal_serendip(&pbasis, cdim + vdim, order);
    gkyl_cart_modal_serendip(&vbasis, vdim, order);
  }

  double lower[5] = {-1, -1, -1, -1, -1}, upper[5] = {1, 1, 1, 1, 1};
  int cells[5] = {2, 2, 2, 2, 2}, ghosts[5] = {0};

  if (vdim == 2) {
    lower[cdim + 1] = 0.2;
    upper[cdim + 1] = 1.8;
  }
  for (int d = 0; d < cdim; ++d) {
    ghosts[d] = 1;
  }

  struct gkyl_rect_grid grid, cgrid, vgrid;
  gkyl_rect_grid_init(&grid, cdim + vdim, lower, upper, cells);
  gkyl_rect_grid_init(&cgrid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&vgrid, vdim, lower + cdim, upper + cdim, cells + cdim);

  struct gkyl_range phase, phase_ext, conf, conf_ext, vel, vel_ext;
  gkyl_create_grid_ranges(&grid, ghosts, &phase_ext, &phase);
  gkyl_create_grid_ranges(&cgrid, ghosts, &conf_ext, &conf);
  gkyl_create_grid_ranges(&vgrid, (int[]){0, 0}, &vel_ext, &vel);

  struct gkyl_gk_geometry_inp geom_in = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .grid = cgrid,
    .basis = cbasis,
    .local = conf,
    .local_ext = conf_ext,
    .global = conf,
    .global_ext = conf_ext,
  };
  struct gk_geometry *geom = gkyl_gk_geometry_new(NULL, &geom_in, gpu);

  struct gkyl_array *bmag = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, conf_ext.volume);
  for (long cell = 0; cell < bmag->size; ++cell) {
    double *coeff = gkyl_array_fetch(bmag, cell);
    coeff[0] = 1.3 * pow(sqrt(2.0), cdim);
    coeff[1] = 0.2;
  }
  gkyl_array_copy(geom->geo_corn.bmag, bmag);

  struct gkyl_velocity_map *vmap = gkyl_velocity_map_new(
    (struct gkyl_mapc2p_inp){.mapping = mapped ? map_velocity : NULL, .ctx = &vdim}, grid, vgrid,
    phase, phase_ext, vel, vel_ext, gpu
  );
  struct gkyl_array *vmap_ho = gkyl_array_new(GKYL_DOUBLE, vmap->vmap->ncomp, vmap->vmap->size);
  gkyl_array_copy(vmap_ho, vmap->vmap);

  struct gkyl_array *phi_ho = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, conf_ext.volume);
  for (long cell = 0; cell < phi_ho->size; ++cell) {
    double *coeff = gkyl_array_fetch(phi_ho, cell);
    coeff[0] = 2 * pow(sqrt(2.0), cdim);
    coeff[1] = 0.3;
  }

  struct gkyl_array *phi = mkarr(gpu, cbasis.num_basis, conf_ext.volume);
  gkyl_array_copy(phi, phi_ho);

  struct gkyl_array *fin_ho = gkyl_array_new(GKYL_DOUBLE, pbasis.num_basis, phase_ext.volume);
  for (long cell = 0; cell < fin_ho->size; ++cell) {
    double *coeff = gkyl_array_fetch(fin_ho, cell);
    for (int k = 0; k < pbasis.num_basis; ++k) {
      coeff[k] = (k % 2 ? -1 : 1) * (0.1 + 0.003 * k + 0.001 * cell);
    }
  }

  struct gkyl_array *fin = mkarr(gpu, pbasis.num_basis, phase_ext.volume);
  gkyl_array_copy(fin, fin_ho);

  const enum gkyl_distribution_moments moments[] = {
    GKYL_F_MOMENT_M0,
    GKYL_F_MOMENT_M1,
    GKYL_F_MOMENT_M2,
    GKYL_F_MOMENT_M2PAR,
    GKYL_F_MOMENT_M2PERP,
    GKYL_F_MOMENT_M0M1M2,
    GKYL_F_MOMENT_M0M1M2PARM2PERP,
    GKYL_F_MOMENT_M3,
    GKYL_F_MOMENT_M0M1M2M3,
    GKYL_F_MOMENT_M0M1M2PARM2PERPM3PARM3PERP,
    GKYL_F_MOMENT_HAMILTONIAN
  };

  const int mapping[11][6] = {{0},       {1},          {2}, {3},          {4},
                              {0, 1, 2}, {0, 1, 3, 4}, {5}, {0, 1, 2, 5}, {0, 1, 3, 4, 6, 7},
                              {0, 1, 2}};

  for (int space = 0; space < 5; ++space) {
    const struct gkyl_basis *wbasis = space == 1 || space == 4 ? &vbasis : &cbasis;
    int weight_dim = wbasis->ndim;
    int wlower[3] = {-2, -2, -2}, wupper[3] = {4, 4, 4};

    struct gkyl_range wrange;
    gkyl_range_init(&wrange, weight_dim, wlower, wupper);

    struct gkyl_array *weight_ho = gkyl_array_new(GKYL_DOUBLE, wbasis->num_basis, wrange.volume);
    for (long cell = 0; cell < weight_ho->size; ++cell) {
      double *coeff = gkyl_array_fetch(weight_ho, cell);
      for (int k = 0; k < wbasis->num_basis; ++k) {
        coeff[k] = space >= 2 ? (k == 0 ? pow(sqrt(2.0), weight_dim) : 0) :
                                (k % 2 ? -0.13 : 0.17) + 0.002 * cell;
      }
    }

    struct gkyl_array *weight = mkarr(gpu, wbasis->num_basis, wrange.volume);
    gkyl_array_copy(weight, weight_ho);

    for (int integrated = 0; integrated < 2; ++integrated) {
      int ncoeff = integrated ? 1 : cbasis.num_basis;
      struct gkyl_array *reference = gkyl_array_new(GKYL_DOUBLE, 9 * ncoeff, conf_ext.volume);

      // Five points integrate all products here exactly, including velocity weights.
      int nquad = 1;
      for (int d = 0; d < grid.ndim; ++d) {
        nquad *= 5;
      }

      struct gkyl_range_iter iter;
      gkyl_range_iter_init(&iter, &phase);
      while (gkyl_range_iter_next(&iter)) {
        long pcell = gkyl_range_idx(&phase, iter.idx), ccell = gkyl_range_idx(&conf, iter.idx);
        long wcell = gkyl_range_idx(&wrange, space == 1 || space == 4 ? iter.idx + cdim : iter.idx);

        const double *fc = gkyl_array_cfetch(fin_ho, pcell),
                     *wc = gkyl_array_cfetch(weight_ho, wcell);
        const double *bc = gkyl_array_cfetch(bmag, ccell);
        double *expected = gkyl_array_fetch(reference, ccell);

        for (int quad = 0; quad < nquad; ++quad) {
          int index = quad;
          double z[5], factor = vdim == 2 ? 2 * M_PI / 2.0 : 1.0;
          for (int d = 0; d < grid.ndim; ++d) {
            int node = index % 5;
            index /= 5;
            z[d] = gkyl_gauss_ordinates[5][node];
            factor *= gkyl_gauss_weights[5][node];
            if (d >= cdim || integrated) {
              factor *= grid.dx[d] / 2;
            }
          }

          const double *vm = gkyl_array_cfetch(vmap_ho, gkyl_range_idx(&vel, iter.idx + cdim));
          double velocity = vm[0] / sqrt(2.0) + sqrt(1.5) * vm[1] * z[cdim];
          double mu = vdim == 2 ? vm[2] / sqrt(2.0) + sqrt(1.5) * vm[3] * z[cdim + 1] : 0;
          double perp = mu * eval(&cbasis, z, bc); // 2*mu*B/m, with mass=2.
          double raw[] = {
            1,
            velocity,
            velocity * velocity + perp,
            velocity * velocity,
            perp,
            velocity * (velocity * velocity + perp),
            velocity * velocity * velocity,
            velocity * perp,
            eval(&cbasis, z, gkyl_array_cfetch(phi_ho, ccell))
          };

          double basis_values[20];
          cbasis.eval(z, basis_values);

          factor *=
            eval(&pbasis, z, fc) * eval(wbasis, space == 1 || space == 4 ? z + cdim : z, wc);

          for (int mom = 0; mom < 9; ++mom) {
            for (int k = 0; k < ncoeff; ++k) {
              expected[mom * ncoeff + k] += factor * raw[mom] * (integrated ? 1 : basis_values[k]);
            }
          }
        }
      }

      for (int which = 0; which < (space == 2 ? 11 : 7); ++which) {
        if (vdim == 1 && (which == 4 || which == 9)) {
          continue;
        }

        struct gkyl_mom_weighted_gyrokinetic *up =
          space == 2 ?
            NULL :
            gkyl_mom_weighted_gyrokinetic_new(
              2.0, 1.0, &cbasis, &pbasis, &grid, vmap, geom, moments[which],
              space == 1 || space == 4 ? GKYL_F_MOMENT_WEIGHT_VEL : GKYL_F_MOMENT_WEIGHT_CONF,
              integrated, gpu
            );

        struct gkyl_mom_gyrokinetic *unweighted =
          space == 2 ?
            gkyl_mom_gyrokinetic_new(
              2.0, 1.0, &cbasis, &pbasis, &grid, vmap, geom, moments[which], integrated, gpu
            ) :
            NULL;

        int nm = unweighted ? gkyl_mom_gyrokinetic_num_mom(unweighted) :
                              gkyl_mom_weighted_gyrokinetic_num_mom(up);

        // Check supported weighted constructor options, including the 1v FourMoments alias.
        if (up) {
          int expected_num_mom = which == 5 ? 3 : which == 6 ? (vdim == 1 ? 3 : 4) : 1;
          TEST_CHECK(nm == expected_num_mom);
        }

        struct gkyl_array *out = mkarr(gpu, nm * ncoeff, conf_ext.volume),
                          *out_ho = mkarr(false, nm * ncoeff, conf_ext.volume);

        for (int repeat = 0; repeat < (space == 2 && which != 10 ? 1 : 3); ++repeat) {
          double scale = repeat == 0 ? 1 : repeat == 1 ? -0.5 : 0;

          gkyl_array_copy(weight, weight_ho);
          gkyl_array_scale(weight, scale);
          gkyl_array_copy(phi, phi_ho);
          gkyl_array_scale(phi, scale);

          struct gkyl_range selected_conf = conf;
          if (repeat == 1) {
            int selected_upper[3];
            for (int d = 0; d < cdim; ++d) {
              selected_upper[d] = conf.upper[d];
            }
            selected_upper[0] = conf.lower[0];
            gkyl_sub_range_init(&selected_conf, &conf, conf.lower, selected_upper);
          }
          gkyl_array_clear(out, 97.0);

          if (unweighted) {
            gkyl_mom_gyrokinetic_advance(unweighted, &phase, &selected_conf, phi, fin, out);
          } else {
            gkyl_mom_weighted_gyrokinetic_advance(
              up, &phase, &selected_conf, &wrange, weight, NULL, fin, out
            );
          }

          gkyl_array_copy(out_ho, out);
          gkyl_range_iter_init(&iter, &conf_ext);
          while (gkyl_range_iter_next(&iter)) {
            long cell = gkyl_range_idx(&conf_ext, iter.idx);

            const double *actual = gkyl_array_cfetch(out_ho, cell),
                         *ref = gkyl_array_cfetch(reference, cell);

            bool selected = gkyl_range_contains_idx(&selected_conf, iter.idx);
            for (int mom = 0; mom < nm; ++mom) {
              for (int k = 0; k < ncoeff; ++k) {
                double target = selected ? scale * ref[mapping[which][mom] * ncoeff + k] : 97.0;
                if (which == 10 && selected) {
                  target = ref[mapping[which][mom] * ncoeff + k];
                  if (mom == 1) {
                    target *= 2; // m*M1, with m=2.
                  }
                  if (mom == 2) {
                    target += scale * ref[8 * ncoeff + k]; // m*M2/2 + q*phi*M0.
                  }
                }
                TEST_CHECK(
                  fabs(actual[mom * ncoeff + k] - target) < 2e-11 * fmax(1.0, fabs(target))
                );
                TEST_MSG(
                  "%dx%dv p%d weight=%d integrated=%d moment=%d: %.17g != %.17g", cdim, vdim, order,
                  space, integrated, which, actual[mom * ncoeff + k], target
                );
              }
            }
          }
        }

        // An empty velocity selection clears selected outputs; empty configuration is a no-op.
        struct gkyl_range empty_phase = phase, empty_conf = conf;
        empty_phase.upper[cdim] = empty_phase.lower[cdim] - 1;
        empty_phase.volume = 0;
        empty_conf.volume = 0;
        gkyl_array_clear(out, 97);

        if (unweighted) {
          gkyl_mom_gyrokinetic_advance(unweighted, &phase, &empty_conf, phi, fin, out);
        } else {
          gkyl_mom_weighted_gyrokinetic_advance(
            up, &phase, &empty_conf, &wrange, weight, NULL, fin, out
          );
        }
        gkyl_array_copy(out_ho, out);

        TEST_CHECK(
          ((const double *)gkyl_array_cfetch(out_ho, gkyl_range_idx(&conf, conf.lower)))[0] == 97
        );

        if (unweighted) {
          gkyl_mom_gyrokinetic_advance(unweighted, &empty_phase, &conf, phi, fin, out);
        } else {
          gkyl_mom_weighted_gyrokinetic_advance(
            up, &empty_phase, &conf, &wrange, weight, NULL, fin, out
          );
        }

        gkyl_array_copy(out_ho, out);
        const double *empty_data = gkyl_array_cfetch(out_ho, gkyl_range_idx(&conf, conf.lower));
        for (int k = 0; k < nm * ncoeff; ++k) {
          TEST_CHECK(empty_data[k] == 0);
        }

        if (unweighted) {
          gkyl_mom_gyrokinetic_release(unweighted);
        } else {
          gkyl_mom_weighted_gyrokinetic_release(up);
        }
        gkyl_array_release(out);
        gkyl_array_release(out_ho);
      }
      gkyl_array_release(reference);
    }
    gkyl_array_release(weight);
    gkyl_array_release(weight_ho);
  }

  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(vmap_ho);
  gkyl_array_release(fin);
  gkyl_array_release(fin_ho);
  gkyl_array_release(bmag);
  gkyl_velocity_map_release(vmap);
  gkyl_gk_geometry_release(geom);
}

static void
test_mom_weighted_gk_ho(void)
{
  for (int order = 1; order <= 2; ++order) {
    check_weighted(1, 1, order, false, false);
    check_weighted(1, 2, order, false, false);
    check_weighted(2, 2, order, false, false);
  }
  check_weighted(3, 2, 1, false, false);
  check_weighted(1, 2, 1, false, true);
  check_weighted(1, 2, 2, false, true);
}

#ifdef GKYL_HAVE_CUDA
static void
test_mom_weighted_gk_dev(void)
{
  for (int order = 1; order <= 2; ++order) {
    check_weighted(1, 1, order, true, false);
    check_weighted(1, 2, order, true, false);
    check_weighted(2, 2, order, true, false);
  }
  check_weighted(3, 2, 1, true, false);
  check_weighted(1, 2, 1, true, true);
  check_weighted(1, 2, 2, true, true);
}
#endif

TEST_LIST = {
  {"mom_weighted_gk_ho", test_mom_weighted_gk_ho},
#ifdef GKYL_HAVE_CUDA
  {"mom_weighted_gk_dev", test_mom_weighted_gk_dev},
#endif
  {NULL, NULL}
};
