#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_position_map.h>
#include <gkyl_rect_decomp.h>
#include <math.h>

struct quadratic_map_ctx {
  double shift, scale, curvature;
};

static void quadratic_map(double t, const double *xn, double *out, void *ctx)
{
  const struct quadratic_map_ctx *map = ctx;
  out[0] = map->shift + xn[0] * (map->scale + map->curvature * xn[0]);
}

static void quadratic_map_deriv(double t, const double *xn, double *out, void *ctx)
{
  const struct quadratic_map_ctx *map = ctx;
  out[0] = map->scale + 2.0 * map->curvature * xn[0];
}

// Independent global-index construction of P1 corners and Gauss nodes.
static double node_coordinate(double lower, double dx, int node, bool quadrature)
{
  return lower +
         dx * (quadrature ? node / 2 + 0.5 * (1.0 + (node % 2 ? 1.0 : -1.0) / sqrt(3.0)) : node);
}

static void check_mirror_mapping(
  bool analytic, enum gkyl_mirror_grid_gen_field_line_coord coord, struct quadratic_map_ctx map,
  bool map_all, bool include_axis
)
{
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(
    &grid, 3, (double[]){0.02, -0.4, -0.75}, (double[]){0.2, 0.4, 0.75}, (int[]){8, 2, 16}
  );
  struct gkyl_range global, ext;
  gkyl_create_grid_ranges(&grid, (int[]){1, 1, 1}, &ext, &global);

  // psi = R^2 (1+Z^2)/2 is represented exactly by the bicubic interpolant.
  struct gkyl_range psi_nodes;
  gkyl_range_init_from_shape(&psi_nodes, 2, (int[]){9, 17});
  struct gkyl_array *psi = gkyl_array_new(GKYL_DOUBLE, 1, psi_nodes.volume);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &psi_nodes);
  while (gkyl_range_iter_next(&iter)) {
    double radius = iter.idx[0] / 8.0, height = -1.0 + iter.idx[1] / 8.0;
    double *value = gkyl_array_fetch(psi, gkyl_range_idx(&psi_nodes, iter.idx));
    value[0] = 0.5 * radius * radius * (1.0 + height * height);
  }
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();
  pmap->maps[2] = quadratic_map;
  pmap->map_derivs[2] = quadratic_map_deriv;
  pmap->ctxs[2] = &map;
  pmap->use_map_derivs = analytic;
  if (map_all) {
    for (int dim = 0; dim < 2; ++dim) {
      pmap->maps[dim] = quadratic_map;
      pmap->map_derivs[dim] = quadratic_map_deriv;
      pmap->ctxs[dim] = &map;
    }
  }

  double slope_error = 0.0, tangent_error = 0.0, jacobian_error = 0.0, partition_error = 0.0,
         curl_error = 0.0, normal_field_error = 0.0;
  for (int kind = -2; kind < 3; ++kind) { // Interior, corners, then each surface orientation.
    for (int split_dim = 0; split_dim < 3; ++split_dim) {
      struct gkyl_mirror_grid_gen *full = 0;
      struct gkyl_range full_nodes;
      for (int cuts = 1; cuts <= 8 && cuts <= grid.cells[split_dim]; cuts *= 2) {
        for (int rank = 0; rank < cuts; ++rank) {
          int lower[3], upper[3], shape[3];
          bool quad[3];
          for (int dim = 0; dim < 3; ++dim) {
            lower[dim] = global.lower[dim];
            upper[dim] = global.upper[dim];
            if (dim == split_dim) {
              lower[dim] += rank * grid.cells[dim] / cuts;
              upper[dim] = lower[dim] + grid.cells[dim] / cuts - 1;
            }
            quad[dim] = kind == -2 || (kind >= 0 && dim != kind);
            shape[dim] = quad[dim] ? 2 * (upper[dim] - lower[dim] + 1) :
                                     upper[dim] - lower[dim] + 2;
          }
          struct gkyl_range local, nodes;
          gkyl_range_init(&local, 3, lower, upper);
          gkyl_range_init_from_shape(&nodes, 3, shape);
          struct gkyl_mirror_grid_gen_inp inp = {
            .comp_grid = &grid,
            .local = local,
            .global = global,
            .nrange = nodes,
            .R = {0.0, 1.0},
            .Z = {-1.0, 1.0},
            .nrcells = 8,
            .nzcells = 16,
            .psiRZ = psi,
            .position_map = pmap,
            .fl_coord = coord,
            .dir = kind,
            .include_axis = include_axis
          };
          struct gkyl_mirror_grid_gen *geom = kind == -2 ? gkyl_mirror_grid_gen_int_inew(&inp) :
                                              kind == -1 ? gkyl_mirror_grid_gen_inew(&inp) :
                                                           gkyl_mirror_grid_gen_surf_inew(&inp);
          TEST_ASSERT(geom != 0);
          if (cuts == 1) {
            full = geom;
            full_nodes = nodes;
          }
          gkyl_range_iter_init(&iter, &nodes);
          while (gkyl_range_iter_next(&iter)) {
            int index[3] = {iter.idx[0], iter.idx[1], iter.idx[2]};
            index[split_dim] += (quad[split_dim] ? 2 : 1) * rank * grid.cells[split_dim] / cuts;
            long loc = gkyl_range_idx(&nodes, iter.idx),
                 full_loc = gkyl_range_idx(&full_nodes, index);
            const struct gkyl_mirror_grid_gen_geom *geo = gkyl_array_cfetch(geom->nodes_geom, loc);
            const struct gkyl_mirror_grid_gen_geom *ref =
              gkyl_array_cfetch(full->nodes_geom, full_loc);
            const double *rz = gkyl_array_cfetch(geom->nodes_rza, loc);
            double comp = node_coordinate(grid.lower[2], grid.dx[2], index[2], quad[2]), slope;
            quadratic_map_deriv(0.0, &comp, &slope, &map);
            slope_error = fmax(slope_error, fabs(geo->tang[2].x[2] - slope));
            double expected_jac = slope / (1.0 + rz[1] * rz[1]);
            if (coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
              expected_jac *= 2.0 * sqrt(0.5 * rz[0] * rz[0] * (1.0 + rz[1] * rz[1]));
            }
            if (map_all) {
              for (int dim = 0; dim < 2; ++dim) {
                double coord_lo = grid.lower[dim], coord_dx = grid.dx[dim];
                if (dim == 0 && coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
                  coord_lo = sqrt(coord_lo);
                  coord_dx = (sqrt(grid.upper[0]) - coord_lo) / grid.cells[0];
                }
                double xc = node_coordinate(coord_lo, coord_dx, index[dim], quad[dim]), ds;
                quadratic_map_deriv(0.0, &xc, &ds, &map);
                expected_jac *= ds;
              }
            }
            // B_R = R Z, B_Z = -(1+Z^2), B_phi = 0.
            double radius = rz[0], height = rz[1], bz = -(1.0 + height * height),
                   br = radius * height;
            double bmag = sqrt(br * br + bz * bz), dbdr = radius * height * height / bmag;
            double dbdz = (radius * radius * height - 2.0 * height * bz) / bmag;
            double curl_phi = radius / bmag + (dbdr * bz - dbdz * br) / (bmag * bmag);
            curl_error = fmax(curl_error, fabs(radius * geo->curlbhat.x[1] - curl_phi));
            double normal_field = geo->dual[2].x[2] * geo->B.x[2] / fabs(geo->dual[2].x[2]);
            normal_field_error = fmax(normal_field_error, fabs(normal_field - bz));
            jacobian_error = fmax(jacobian_error, fabs(geo->Jc - expected_jac));
            partition_error = fmax(partition_error, fabs(geo->Jc - ref->Jc));
            // Reciprocal bases, including the constant-psi field-line tangent.
            for (int idir = 0; idir < 3; ++idir) {
              for (int jdir = 0; jdir < 3; ++jdir) {
                double dot = geo->tang[idir].x[0] * geo->dual[jdir].x[0] +
                             rz[0] * rz[0] * geo->tang[idir].x[1] * geo->dual[jdir].x[1] +
                             geo->tang[idir].x[2] * geo->dual[jdir].x[2];
                if (radius > 0.0) {
                  tangent_error = fmax(tangent_error, fabs(dot - (idir == jdir)));
                }
              }
            }
          }
          if (cuts != 1) {
            gkyl_mirror_grid_gen_release(geom);
          }
        }
      }
      gkyl_mirror_grid_gen_release(full);
    }
  }
  TEST_CHECK(slope_error < 2e-11);
  TEST_MSG("mapping slope error %.17g", slope_error);
  TEST_CHECK(tangent_error < 2e-13);
  TEST_MSG("reciprocal basis error %.17g", tangent_error);
  TEST_CHECK(jacobian_error < 2e-10);
  TEST_MSG("Jacobian error %.17g", jacobian_error);
  TEST_CHECK(curl_error < 2e-10);
  TEST_MSG("curl error %.17g", curl_error);
  TEST_CHECK(normal_field_error < 2e-10);
  TEST_MSG("normal B error %.17g", normal_field_error);
  TEST_CHECK(partition_error < 2e-12);
  TEST_MSG("partition error %.17g", partition_error);
  gkyl_position_map_release(pmap);
  gkyl_array_release(psi);
}

static void test_mirror_identity(void)
{
  check_mirror_mapping(
    false, GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z, (struct quadratic_map_ctx){0.0, 1.0, 0.0},
    false, false
  );
}
static void test_mirror_affine(void)
{
  check_mirror_mapping(
    true, GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z, (struct quadratic_map_ctx){0.05, 0.75, 0.0},
    false, false
  );
}
static void test_mirror_nonlinear_numeric(void)
{
  check_mirror_mapping(
    false, GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z, (struct quadratic_map_ctx){0.05, 0.75, 0.12},
    false, false
  );
}
static void test_mirror_nonlinear_analytic(void)
{
  check_mirror_mapping(
    true, GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z, (struct quadratic_map_ctx){0.05, 0.75, 0.12},
    false, false
  );
}
static void test_mirror_sqrt_psi(void)
{
  check_mirror_mapping(
    false, GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z,
    (struct quadratic_map_ctx){0.05, 0.75, 0.12}, false, false
  );
}

static void curved_cartesian_map(double t, const double *xn, double *out, void *ctx)
{
  for (int dim = 0; dim < 3; ++dim) {
    out[dim] = xn[dim];
  }
  out[0] += 0.1 * xn[2] * xn[2];
}

static void curved_field(double t, const double *xn, double *out, void *ctx)
{
  out[0] = 0.4 * xn[2];
  out[1] = 0.0;
  out[2] = 2.0;
}

static void compare_nodes(
  const struct gkyl_array *actual, const struct gkyl_array *expected,
  const struct gkyl_range *nodes, const struct gkyl_range *full_nodes, int split_dim,
  int node_offset, double *error, const char *name
)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, nodes);
  while (gkyl_range_iter_next(&iter)) {
    int index[3] = {iter.idx[0], iter.idx[1], iter.idx[2]};
    index[split_dim] += node_offset;
    const double *value = gkyl_array_cfetch(actual, gkyl_range_idx(nodes, iter.idx));
    const double *ref = gkyl_array_cfetch(expected, gkyl_range_idx(full_nodes, index));
    for (int comp = 0; comp < actual->ncomp; ++comp) {
      union {
        double value;
        uint64_t bits;
      } actual_bits = {.value = value[comp]}, ref_bits = {.value = ref[comp]};
      bool finite = (actual_bits.bits & UINT64_C(0x7ff0000000000000)) !=
                      UINT64_C(0x7ff0000000000000) &&
                    (ref_bits.bits & UINT64_C(0x7ff0000000000000)) != UINT64_C(0x7ff0000000000000);
      TEST_CHECK(finite);
      TEST_MSG(
        "%s at %d %d %d component %d: actual %g reference %g", name, iter.idx[0], iter.idx[1],
        iter.idx[2], comp, value[comp], ref[comp]
      );
      if (!finite) {
        return;
      }
      *error = fmax(*error, fabs(value[comp] - ref[comp]) / fmax(1.0, fabs(ref[comp])));
    }
  }
}

static void check_geometry_partitions(bool tokamak, bool analytic)
{
  struct gkyl_rect_grid grid;
  double pi = M_PI;
  gkyl_rect_grid_init(
    &grid, 3, (double[]){0.1, -0.4, -pi + 1e-14}, (double[]){0.2, 0.4, pi - 1e-14},
    (int[]){4, 2, 16}
  );
  struct gkyl_range global, global_ext;
  gkyl_create_grid_ranges(&grid, (int[]){1, 1, 1}, &global_ext, &global);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 3, 1);
  struct quadratic_map_ctx maps[] = {{0.02, 0.8, 0.2}, {0.05, 0.8, 0.2}, {0.0, 0.7, 0.03}};
  struct gkyl_position_map_inp map_inp = {
    .id = analytic ? GKYL_PMAP_USER_INPUT_W_DERIVATIVE : GKYL_PMAP_USER_INPUT,
    .maps = {quadratic_map, tokamak ? 0 : quadratic_map, quadratic_map},
    .map_derivs = {quadratic_map_deriv, quadratic_map_deriv, quadratic_map_deriv},
    .ctxs = {&maps[0], &maps[1], &maps[2]}
  };
  struct gkyl_position_map *pmap =
    gkyl_position_map_new(map_inp, grid, global, global_ext, global, global_ext, basis);
  struct gkyl_gk_geometry_inp inp = {
    .geometry_id = tokamak ? GKYL_GEOMETRY_TOKAMAK : GKYL_GEOMETRY_MAPC2P,
    .mapc2p = curved_cartesian_map,
    .bfield_func = curved_field,
    .position_map = pmap,
    .efit_info =
      {.filepath = "gyrokinetic/data/eqdsk/straight_cylinder.geqdsk",
       .rz_poly_order = 2,
       .flux_poly_order = 1,
       .reflect = true},
    .tok_grid_info =
      {.rclose = 0.5, .zmin = -1.0, .zmax = 1.0, .rleft = 0.001, .rmax = 1.0, .rright = 1.0},
    .grid = grid,
    .basis = basis,
    .global = global,
    .global_ext = global_ext,
    .geo_grid = grid,
    .geo_basis = basis,
    .geo_global = global,
    .geo_global_ext = global_ext,
    .local = global,
    .local_ext = global_ext,
    .geo_local = global,
    .geo_local_ext = global_ext
  };
  struct gk_geometry *full = tokamak ? gkyl_gk_geometry_tok_new(&inp) :
                                       gkyl_gk_geometry_mapc2p_new(&inp);
  double error = 0.0, slope_error = 0.0, curl_error = 0.0;
  for (int split_dim = 0; split_dim < 3; ++split_dim) {
    for (int cuts = 2; cuts <= grid.cells[split_dim] && cuts <= 8; cuts *= 2) {
      for (int rank = 0; rank < cuts; ++rank) {
        int lower[3], upper[3], elo[3], eup[3];
        for (int dim = 0; dim < 3; ++dim) {
          lower[dim] = global.lower[dim];
          upper[dim] = global.upper[dim];
          if (dim == split_dim) {
            lower[dim] += rank * grid.cells[dim] / cuts;
            upper[dim] = lower[dim] + grid.cells[dim] / cuts - 1;
          }
          elo[dim] = lower[dim] - 1;
          eup[dim] = upper[dim] + 1;
        }
        struct gkyl_range local, ext;
        gkyl_range_init(&ext, 3, elo, eup);
        gkyl_sub_range_init(&local, &ext, lower, upper);
        inp.local = inp.geo_local = local;
        inp.local_ext = inp.geo_local_ext = ext;
        struct gk_geometry *geom = tokamak ? gkyl_gk_geometry_tok_new(&inp) :
                                             gkyl_gk_geometry_mapc2p_new(&inp);
        int offset = rank * grid.cells[split_dim] / cuts;
        compare_nodes(
          geom->geo_corn.mc2p_nodal, full->geo_corn.mc2p_nodal, &geom->nrange_corn,
          &full->nrange_corn, split_dim, offset, &error, "corners"
        );
        compare_nodes(
          geom->geo_corn.mc2nu_pos_nodal, full->geo_corn.mc2nu_pos_nodal, &geom->nrange_corn,
          &full->nrange_corn, split_dim, offset, &error, "corners"
        );
#define COMPARE_INT(field)                                                                     \
  compare_nodes(                                                                               \
    geom->geo_int.field, full->geo_int.field, &geom->nrange_int, &full->nrange_int, split_dim, \
    2 * offset, &error, "interior " #field                                                     \
  )
        COMPARE_INT(mc2p_nodal);
        COMPARE_INT(jacobgeo_nodal);
        COMPARE_INT(g_ij_nodal);
        COMPARE_INT(dxdz_nodal);
        COMPARE_INT(dzdx_nodal);
        COMPARE_INT(bmag_nodal);
        COMPARE_INT(curlbhat_nodal);
        COMPARE_INT(B3_nodal);
#undef COMPARE_INT
        compare_nodes(
          geom->geo_int.bcart, full->geo_int.bcart, &geom->local, &full->local, split_dim, 0,
          &error, "modal Cartesian b"
        );
        for (int dir = 0; dir < 3; ++dir) {
#define COMPARE_SURF(field)                                                          \
  compare_nodes(                                                                     \
    geom->geo_surf[dir].field, full->geo_surf[dir].field, &geom->nrange_surf[dir],   \
    &full->nrange_surf[dir], split_dim, (dir == split_dim ? 1 : 2) * offset, &error, \
    "surface " #field                                                                \
  )
          COMPARE_SURF(jacobgeo_nodal);
          COMPARE_SURF(g_ij_nodal);
          COMPARE_SURF(dxdz_nodal);
          COMPARE_SURF(dzdx_nodal);
          COMPARE_SURF(bmag_nodal);
          COMPARE_SURF(curlbhat_nodal);
          COMPARE_SURF(B3_nodal);
          COMPARE_SURF(lenr_nodal);
#undef COMPARE_SURF
          // Independent axial derivative and physical B projection in a
          // Curved Cartesian map or a straight cylinder.
          struct gkyl_range_iter iter;
          gkyl_range_iter_init(&iter, &geom->nrange_surf[dir]);
          while (gkyl_range_iter_next(&iter)) {
            int znode = iter.idx[2] + (dir == 2 ? 1 : 2) * (lower[2] - global.lower[2]);
            double zcomp = node_coordinate(grid.lower[2], grid.dx[2], znode, dir != 2), slope;
            quadratic_map_deriv(0.0, &zcomp, &slope, &maps[2]);
            const double *tangent = gkyl_array_cfetch(
              geom->geo_surf[dir].dxdz_nodal, gkyl_range_idx(&geom->nrange_surf[dir], iter.idx)
            );
            slope_error = fmax(slope_error, fabs(tangent[8] - slope / (tokamak ? pi : 1.0)));
            const double *curl = gkyl_array_cfetch(
              geom->geo_surf[dir].curlbhat_nodal, gkyl_range_idx(&geom->nrange_surf[dir], iter.idx)
            );
            double zphys;
            quadratic_map(0.0, &zcomp, &zphys, &maps[2]);
            double curl_y = tokamak ? 0.0 : 0.2 / pow(1.0 + 0.04 * zphys * zphys, 1.5);
            for (int dim = 0; dim < 3; ++dim) {
              curl_error = fmax(curl_error, fabs(curl[dim] - (dim == 1 ? curl_y : 0.0)));
            }
          }
        }
        gkyl_gk_geometry_release(geom);
      }
    }
  }
  TEST_CHECK(curl_error < 2e-6);
  TEST_MSG("analytic curl error %.17g", curl_error);
  TEST_CHECK(error < 2e-10);
  TEST_MSG("partition error %.17g", error);
  TEST_CHECK(slope_error < 2e-10);
  TEST_MSG("axial tangent error %.17g", slope_error);
  gkyl_gk_geometry_release(full);
  gkyl_position_map_release(pmap);
}

static void test_mapc2p_partitions(void)
{
  check_geometry_partitions(false, false);
}
static void test_tokamak_partitions_numeric(void)
{
  check_geometry_partitions(true, false);
}
static void test_tokamak_partitions_analytic(void)
{
  check_geometry_partitions(true, true);
}

static void test_mirror_all_maps(void)
{
  check_mirror_mapping(
    false, GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z, (struct quadratic_map_ctx){0.05, 0.75, 0.12},
    true, false
  );
}
static void test_mirror_nonuniform_axis(void)
{
  check_mirror_mapping(
    false, GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z,
    (struct quadratic_map_ctx){0.05, 0.75, 0.12}, false, true
  );
}

TEST_LIST = {
  {"mirror_all_maps", test_mirror_all_maps},
  {"mirror_nonuniform_axis", test_mirror_nonuniform_axis},
  {"mapc2p_partitions", test_mapc2p_partitions},
  {"tokamak_partitions_numeric", test_tokamak_partitions_numeric},
  {"tokamak_partitions_analytic", test_tokamak_partitions_analytic},
  {"mirror_identity", test_mirror_identity},
  {"mirror_affine", test_mirror_affine},
  {"mirror_nonlinear_numeric", test_mirror_nonlinear_numeric},
  {"mirror_nonlinear_analytic", test_mirror_nonlinear_analytic},
  {"mirror_sqrt_psi", test_mirror_sqrt_psi},
  {NULL, NULL}
};
