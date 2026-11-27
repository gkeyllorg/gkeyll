#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_dg_calc_gk_rad_vars.h>
#include <gkyl_rad_gyrokinetic_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <math.h>

/**
 * Allocate a double array on the host or device.
 *
 * @param on_gpu Whether to allocate on the device.
 * @param nc Number of components per element.
 * @param size Number of elements.
 * @return New array; the caller must release it.
 */
static struct gkyl_array *
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array *a;
  if (on_gpu) {
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  } else {
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  }
  return a;
}

// Compare updater dispatch, density selection, and accumulation with direct host
// kernels. Geometry and velocity maps are unused by the density-weighted advance;
// retain one reference to these stack fixtures throughout the test.
static void
test_nI_nu(int cdim, bool use_gpu)
{
  const int pdim = cdim + 2, num_dens = 3;
  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM] = {0};
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  for (int dir = 0; dir < pdim; ++dir) {
    cells[dir] = dir == 0 ? 6 : 2;
    ghost[dir] = dir < cdim ? 1 : 0;
    lower[dir] = 0.0;
    upper[dir] = cells[dir];
  }

  struct gkyl_rect_grid grid, conf_grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);

  struct gkyl_range phase, phase_ext, conf, conf_ext;
  gkyl_create_grid_ranges(&grid, ghost, &phase_ext, &phase);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_ext, &conf);

  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  gkyl_cart_modal_gkhybrid(&pbasis, cdim, 2);

  struct gk_geometry geom = {0};
  struct gkyl_velocity_map vmap = {0};
  geom.ref_count = gkyl_ref_count_init(NULL);
  vmap.ref_count = gkyl_ref_count_init(NULL);

  struct gkyl_dg_calc_gk_rad_vars *up =
    gkyl_dg_calc_gk_rad_vars_new(&grid, &cbasis, &pbasis, -1.0, 1.0, &geom, &vmap, use_gpu);

  // Electron density, ion density, and temperature, including configuration ghosts.
  struct gkyl_array *field_ho[3], *field[3];
  for (int fld = 0; fld < 3; ++fld) {
    field_ho[fld] = mkarr(false, cbasis.num_basis, conf_ext.volume);
    gkyl_array_clear(field_ho[fld], 0.0);
  }

  const double densities[] = {1.0, 4.0, 9.0};
  const double cutoffs[] = {0.5, 1.0, 1.5};

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &conf);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&conf, iter.idx);

    int fit = (iter.idx[0] - 1) % num_dens;

    double *ne = gkyl_array_fetch(field_ho[0], loc);
    double *nI = gkyl_array_fetch(field_ho[1], loc);
    double *vtsq = gkyl_array_fetch(field_ho[2], loc);

    ne[0] = (densities[fit] + 0.1) * pow(2.0, cdim / 2.0);
    for (int comp = 0; comp < cbasis.num_basis; ++comp) {
      nI[comp] = (0.2 + 0.01 * loc) / (comp + 1);
    }
    // First three cells radiate; the next three are at or below the cutoff.
    vtsq[0] = cutoffs[fit] + (iter.idx[0] <= 3 ? 0.25 : (fit == 0 ? 0.0 : -0.25));
  }

  for (int fld = 0; fld < 3; ++fld) {
    field[fld] = use_gpu ? mkarr(true, cbasis.num_basis, conf_ext.volume) :
                           gkyl_array_acquire(field_ho[fld]);
    gkyl_array_copy(field[fld], field_ho[fld]);
  }
  struct gkyl_array *ne_rad_ho = mkarr(false, 1, num_dens);
  struct gkyl_array *cutoff_ho = mkarr(false, 1, num_dens);

  for (int fit = 0; fit < num_dens; ++fit) {
    ((double *)gkyl_array_fetch(ne_rad_ho, fit))[0] = densities[fit];
    ((double *)gkyl_array_fetch(cutoff_ho, fit))[0] = cutoffs[fit];
  }

  struct gkyl_array *ne_rad = use_gpu ? mkarr(true, 1, num_dens) : gkyl_array_acquire(ne_rad_ho);
  struct gkyl_array *cutoff = use_gpu ? mkarr(true, 1, num_dens) : gkyl_array_acquire(cutoff_ho);
  gkyl_array_copy(ne_rad, ne_rad_ho);
  gkyl_array_copy(cutoff, cutoff_ho);

  int ncomp[] = {2 * cbasis.num_basis, pbasis.num_basis, 3 * cbasis.num_basis, pbasis.num_basis};
  struct gkyl_gk_rad_drag *drag[4];
  struct gkyl_array *drag_ho[4][3], *out[4], *expected[4], *actual[4];
  for (int fld = 0; fld < 4; ++fld) {
    // In 3x, drag coefficients omit y; density-weighted outputs retain it.
    int drag_ncomp = cdim == 3 ? ncomp[fld] / 2 : ncomp[fld];

    drag[fld] =
      gkyl_dg_calc_gk_rad_vars_drag_new(1, &num_dens, drag_ncomp, phase_ext.volume, use_gpu);

    expected[fld] = mkarr(false, ncomp[fld], phase_ext.volume);
    actual[fld] = mkarr(false, ncomp[fld], phase_ext.volume);
    out[fld] = mkarr(use_gpu, ncomp[fld], phase_ext.volume);

    gkyl_array_clear(expected[fld], 0.75 + fld);
    gkyl_array_copy(out[fld], expected[fld]);
    for (int fit = 0; fit < num_dens; ++fit) {
      drag_ho[fld][fit] = mkarr(false, drag_ncomp, phase_ext.volume);
      for (long loc = 0; loc < phase_ext.volume; ++loc) {
        double *coeff = gkyl_array_fetch(drag_ho[fld][fit], loc);
        for (int comp = 0; comp < drag_ncomp; ++comp) {
          coeff[comp] = 0.01 * (1 + fld) * (1 + fit) + 0.001 * comp + 0.0001 * loc;
        }
      }
      gkyl_array_copy(drag[fld]->data[fit].arr, drag_ho[fld][fit]);
    }
  }

  // Two calls must accumulate, not overwrite, the existing outputs.
  for (int repeat = 0; repeat < 2; ++repeat) {
    gkyl_dg_calc_gk_rad_vars_nI_nu_advance(
      up, &conf, &phase, drag[0], drag[1], drag[2], drag[3], ne_rad, field[0], field[1], out[0],
      out[1], out[2], out[3], cutoff, field[2]
    );

    gkyl_range_iter_init(&iter, &phase);
    while (gkyl_range_iter_next(&iter)) {
      if (iter.idx[0] > 3) {
        continue;
      }
      int fit = iter.idx[0] - 1;

      long loc = gkyl_range_idx(&phase, iter.idx);

      const double *nI = gkyl_array_cfetch(field_ho[1], gkyl_range_idx(&conf, iter.idx));

      const double *coeff[4];
      double *result[4];
      for (int fld = 0; fld < 4; ++fld) {
        coeff[fld] = gkyl_array_cfetch(drag_ho[fld][fit], loc);
        result[fld] = gkyl_array_fetch(expected[fld], loc);
      }
      switch (cdim) {
        case 1:
          rad_gyrokinetic_drag_nI_nu_1x2v_ser_p1(
            coeff[0], coeff[1], coeff[2], coeff[3], nI, result[0], result[1], result[2], result[3]
          );
          break;
        case 2:
          rad_gyrokinetic_drag_nI_nu_2x2v_ser_p1(
            coeff[0], coeff[1], coeff[2], coeff[3], nI, result[0], result[1], result[2], result[3]
          );
          break;
        case 3:
          rad_gyrokinetic_drag_nI_nu_3x2v_ser_p1(
            coeff[0], coeff[1], coeff[2], coeff[3], nI, result[0], result[1], result[2], result[3]
          );
          break;
      }
    }

    for (int fld = 0; fld < 4; ++fld) {
      gkyl_array_copy(actual[fld], out[fld]);

      for (long loc = 0; loc < phase_ext.volume; ++loc) {
        const double *want = gkyl_array_cfetch(expected[fld], loc);
        const double *got = gkyl_array_cfetch(actual[fld], loc);
        for (int comp = 0; comp < ncomp[fld]; ++comp) {
          TEST_CHECK(gkyl_compare(got[comp], want[comp], 1e-12));
          TEST_MSG(
            "cdim=%d GPU=%d repeat=%d field=%d cell=%ld component=%d: %.17g != %.17g", cdim,
            use_gpu, repeat, fld, loc, comp, got[comp], want[comp]
          );
        }
      }
    }
  }

  for (int fld = 0; fld < 4; ++fld) {
    for (int fit = 0; fit < num_dens; ++fit) {
      gkyl_array_release(drag_ho[fld][fit]);
    }
    gkyl_dg_calc_gk_rad_vars_drag_release(drag[fld], 1, use_gpu);
    gkyl_array_release(out[fld]);
    gkyl_array_release(expected[fld]);
    gkyl_array_release(actual[fld]);
  }
  for (int fld = 0; fld < 3; ++fld) {
    gkyl_array_release(field[fld]);
    gkyl_array_release(field_ho[fld]);
  }
  gkyl_array_release(ne_rad);
  gkyl_array_release(ne_rad_ho);
  gkyl_array_release(cutoff);
  gkyl_array_release(cutoff_ho);
  gkyl_dg_calc_gk_rad_vars_release(up);

  TEST_CHECK(geom.ref_count.count == 1);
  TEST_CHECK(vmap.ref_count.count == 1);
}

static void
test_nu(int cdim, bool use_gpu)
{
  // Exercise both drag callbacks with distinct scalar arguments and nonuniform
  // magnetic fields. Compare with direct host kernels, including untouched ghosts.
  int pdim = cdim + 2, cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM] = {0};
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  for (int d = 0; d < pdim; ++d) {
    cells[d] = 2;
    ghost[d] = d < cdim ? 1 : 0;
    lower[d] = d == cdim ? -1.2 : 0.2;
    upper[d] = d == cdim ? 1.7 : 1.8;
  }

  struct gkyl_rect_grid grid, cgrid, vgrid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  gkyl_rect_grid_init(&cgrid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&vgrid, 2, lower + cdim, upper + cdim, cells + cdim);

  struct gkyl_range phase, phase_ext, conf, conf_ext, vel, vel_ext;
  gkyl_create_grid_ranges(&grid, ghost, &phase_ext, &phase);
  gkyl_create_grid_ranges(&cgrid, ghost, &conf_ext, &conf);
  gkyl_create_grid_ranges(&vgrid, (int[]){0, 0}, &vel_ext, &vel);

  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  gkyl_cart_modal_gkhybrid(&pbasis, cdim, 2);

  struct gkyl_gk_geometry_inp geom_in = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .grid = cgrid,
    .basis = cbasis,
    .local = conf,
    .local_ext = conf_ext,
    .global = conf,
    .global_ext = conf_ext,
  };
  struct gk_geometry *geom = gkyl_gk_geometry_new(NULL, &geom_in, false);
  if (use_gpu) {
    struct gk_geometry *geom_dev = gkyl_gk_geometry_new(geom, &geom_in, true);
    gkyl_gk_geometry_release(geom);
    geom = geom_dev;
  }

  struct gkyl_array *bmag = mkarr(false, cbasis.num_basis, conf_ext.volume);
  for (long cell = 0; cell < conf_ext.volume; ++cell) {
    double *coeff = gkyl_array_fetch(bmag, cell);
    coeff[0] = 2.4 + 0.001 * cell;
    for (int k = 1; k < cbasis.num_basis; ++k) {
      coeff[k] = 0.01 / k;
    }
  }
  gkyl_array_copy(geom->geo_int.bmag, bmag);

  struct gkyl_velocity_map *map = gkyl_velocity_map_new(
    (struct gkyl_mapc2p_inp){0}, grid, vgrid, phase, phase_ext, vel, vel_ext, use_gpu
  );

  struct gkyl_array *vmap = mkarr(false, map->vmap->ncomp, map->vmap->size);
  struct gkyl_array *vmap_sq = mkarr(false, map->vmap_sq->ncomp, map->vmap_sq->size);
  gkyl_array_copy(vmap, map->vmap);
  gkyl_array_copy(vmap_sq, map->vmap_sq);

  struct gkyl_dg_calc_gk_rad_vars *up =
    gkyl_dg_calc_gk_rad_vars_new(&grid, &cbasis, &pbasis, -1.7, 2.3, geom, map, use_gpu);

  // The 3x drag kernels omit y; density weighting later restores that dimension.
  int nc = cdim == 3 ? cbasis.num_basis / 2 : cbasis.num_basis;
  int ncomp[] = {2 * nc, 6 * nc, 3 * nc, 6 * nc};
  struct gkyl_array *out[4], *expected[4], *actual[4];
  for (int field = 0; field < 4; ++field) {
    out[field] = mkarr(use_gpu, ncomp[field], phase_ext.volume);
    expected[field] = mkarr(false, ncomp[field], phase_ext.volume);
    actual[field] = mkarr(false, ncomp[field], phase_ext.volume);
    gkyl_array_clear(expected[field], 0.625);
    gkyl_array_copy(out[field], expected[field]);
  }

  for (int repeat = 0; repeat < 2; ++repeat) {
    double a = 0.41 + 0.2 * repeat;

    gkyl_dg_calc_gk_rad_vars_nu_advance(
      up, &conf, &phase, a, 1.2, 2.1, 0.7, 1.3, out[0], out[1], out[2], out[3]
    );

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &phase);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&phase, iter.idx);

      long vloc = gkyl_range_idx(&vel, iter.idx + cdim);

      const double *vm = gkyl_array_cfetch(vmap, vloc);
      const double *vm_sq = gkyl_array_cfetch(vmap_sq, vloc);
      const double *bm = gkyl_array_cfetch(bmag, gkyl_range_idx(&conf, iter.idx));

      double *result[4];
      for (int field = 0; field < 4; ++field) {
        result[field] = gkyl_array_fetch(expected[field], loc);
      }

      switch (cdim) {
        case 1:
          rad_gyrokinetic_drag_nuvpar_1x2v_ser_p1(
            vm, vm_sq, -1.7, 2.3, a, 1.2, 2.1, 0.7, 1.3, bm, result[0], result[1]
          );
          rad_gyrokinetic_drag_numu_1x2v_ser_p1(
            vm, vm_sq, -1.7, 2.3, a, 1.2, 2.1, 0.7, 1.3, bm, result[2], result[3]
          );
          break;
        case 2:
          rad_gyrokinetic_drag_nuvpar_2x2v_ser_p1(
            vm, vm_sq, -1.7, 2.3, a, 1.2, 2.1, 0.7, 1.3, bm, result[0], result[1]
          );
          rad_gyrokinetic_drag_numu_2x2v_ser_p1(
            vm, vm_sq, -1.7, 2.3, a, 1.2, 2.1, 0.7, 1.3, bm, result[2], result[3]
          );
          break;
        case 3:
          rad_gyrokinetic_drag_nuvpar_3x2v_ser_p1(
            vm, vm_sq, -1.7, 2.3, a, 1.2, 2.1, 0.7, 1.3, bm, result[0], result[1]
          );
          rad_gyrokinetic_drag_numu_3x2v_ser_p1(
            vm, vm_sq, -1.7, 2.3, a, 1.2, 2.1, 0.7, 1.3, bm, result[2], result[3]
          );
          break;
      }
      TEST_CHECK(isfinite(result[0][0]) && fabs(result[0][0]) > 1e-10);
      TEST_CHECK(isfinite(result[2][0]) && fabs(result[2][0]) > 1e-10);
    }

    for (int field = 0; field < 4; ++field) {
      gkyl_array_copy(actual[field], out[field]);

      for (long cell = 0; cell < phase_ext.volume; ++cell) {
        const double *want = gkyl_array_cfetch(expected[field], cell);
        const double *got = gkyl_array_cfetch(actual[field], cell);
        // Small modal coefficients result from cancellation; scale the error by
        // the largest coefficient in this cell rather than by each component.
        double scale = 0.0;
        for (int k = 0; k < ncomp[field]; ++k) {
          scale = fmax(scale, fabs(want[k]));
        }
        for (int k = 0; k < ncomp[field]; ++k) {
          TEST_CHECK(fabs(got[k] - want[k]) <= 1e-12 * scale);
          TEST_MSG(
            "cdim=%d GPU=%d repeat=%d field=%d cell=%ld comp=%d: %.17g != %.17g", cdim, use_gpu,
            repeat, field, cell, k, got[k], want[k]
          );
        }
      }
    }
  }

  for (int field = 0; field < 4; ++field) {
    gkyl_array_release(out[field]);
    gkyl_array_release(expected[field]);
    gkyl_array_release(actual[field]);
  }
  gkyl_dg_calc_gk_rad_vars_release(up);
  gkyl_array_release(vmap);
  gkyl_array_release(vmap_sq);
  gkyl_array_release(bmag);
  gkyl_velocity_map_release(map);
  gkyl_gk_geometry_release(geom);
}

static void
test_nu_ho()
{
  for (int cdim = 1; cdim <= 3; ++cdim) {
    test_nu(cdim, false);
  }
}

static void
test_1x2v_ho()
{
  test_nI_nu(1, false);
}

static void
test_2x2v_ho()
{
  test_nI_nu(2, false);
}

static void
test_3x2v_ho()
{
  test_nI_nu(3, false);
}

#ifdef GKYL_HAVE_CUDA
static void
test_nu_dev()
{
  for (int cdim = 1; cdim <= 3; ++cdim) {
    test_nu(cdim, true);
  }
}

static void
test_1x2v_dev()
{
  test_nI_nu(1, true);
}

static void
test_2x2v_dev()
{
  test_nI_nu(2, true);
}

static void
test_3x2v_dev()
{
  test_nI_nu(3, true);
}
#endif

TEST_LIST = {
  {"dg_calc_gk_rad_vars_nu_ho", test_nu_ho},
  {"dg_calc_gk_rad_vars_nI_nu_1x2v_ho", test_1x2v_ho},
  {"dg_calc_gk_rad_vars_nI_nu_2x2v_ho", test_2x2v_ho},
  {"dg_calc_gk_rad_vars_nI_nu_3x2v_ho", test_3x2v_ho},
#ifdef GKYL_HAVE_CUDA
  {"dg_calc_gk_rad_vars_nu_dev", test_nu_dev},
  {"dg_calc_gk_rad_vars_nI_nu_1x2v_dev", test_1x2v_dev},
  {"dg_calc_gk_rad_vars_nI_nu_2x2v_dev", test_2x2v_dev},
  {"dg_calc_gk_rad_vars_nI_nu_3x2v_dev", test_3x2v_dev},
#endif
  {NULL, NULL}
};
