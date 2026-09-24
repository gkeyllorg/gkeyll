#include <acutest.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_math.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_position_map.h>

static inline double
SQ(double x)
{
  return x * x;
};

static void
test_wham(bool include_axis, enum gkyl_mirror_grid_gen_field_line_coord fl_coord)
{
  double clower[] = {2.0e-6, 0.0, -2.0};
  double cupper[] = {3.0e-3, 2 * M_PI, 2.0};
  int cells[] = {10, 16, 32};
  int cdim = 3;

  const char *fname = "gyrokinetic/data/unit/wham_hires.geqdsk_psi.gkyl";

  // computational grid
  struct gkyl_rect_grid comp_grid;
  gkyl_rect_grid_init(&comp_grid, cdim, clower, cupper, cells);

  TEST_ASSERT(gkyl_check_file_exists(fname));

  // read psi(R,Z) from file
  struct gkyl_rect_grid psi_grid;
  struct gkyl_array *psi = gkyl_grid_array_new_from_file(&psi_grid, fname);

  struct gkyl_range node_range;
  gkyl_range_init_from_shape(&node_range, 3, (int[3]){cells[0] + 1, cells[1] + 1, cells[2] + 1});

  struct gkyl_range ext_range, range;
  int nghost[3] = {1, 1, 1};
  gkyl_create_grid_ranges(&comp_grid, nghost, &ext_range, &range);

  struct gkyl_position_map *gpm = gkyl_position_map_null_new();

  // create mirror geometry
  struct gkyl_mirror_grid_gen *geom = gkyl_mirror_grid_gen_inew(&(struct gkyl_mirror_grid_gen_inp){
    .comp_grid = &comp_grid,
    .nrange = node_range,
    .local = range,
    .global = range,

    .R = {psi_grid.lower[0], psi_grid.upper[0]},
    .Z = {psi_grid.lower[1], psi_grid.upper[1]},

    // psi(R,Z) grid size
    .nrcells = psi_grid.cells[0] - 1, // cells and not nodes
    .nzcells = psi_grid.cells[1] - 1, // cells and not nodes

    .psiRZ = psi,
    .fl_coord = fl_coord,
    .include_axis = include_axis,
    .write_psi_cubic = false,

    .position_map = gpm,
  });

  TEST_ASSERT(geom != NULL);
  TEST_CHECK(include_axis == gkyl_mirror_grid_gen_is_include_axis(geom));
  TEST_CHECK(fl_coord == gkyl_mirror_grid_gen_fl_coord(geom));

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &node_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&node_range, iter.idx);

    const double *rz = gkyl_array_cfetch(geom->nodes_rza, loc);

    const struct gkyl_mirror_grid_gen_geom *g = gkyl_array_cfetch(geom->nodes_geom, loc);

    // check Jacobian
    double Jac = gkyl_vec3_triple(
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[0]),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[1]),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[2])
    );

    if (rz[0] > 0) {
      TEST_CHECK(Jac > 0.0);
      TEST_CHECK(gkyl_compare_double(Jac, g->Jc, 1e-14));
    }

    // check C = Jc*Bmag/sqrt(g33)
    double g33 = gkyl_vec3_dot(
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[2]),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[2])
    );
    double Bmag = gkyl_vec3_len(gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->B));

    if (fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z) {
      TEST_CHECK(gkyl_compare_double(Bmag * g->Jc / sqrt(g33), 1.0, 1e-14));
    }

    // check B only points in the parallel direction
    // ... B^1 = 0
    double B1 = gkyl_vec3_dot(
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->B),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->dual[0])
    );
    TEST_CHECK(gkyl_compare_double(B1, 0.0, 1e-14));

    // ... B^2 = 0
    double B2 = gkyl_vec3_dot(
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->B),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->dual[1])
    );
    TEST_CHECK(gkyl_compare_double(B2, 0.0, 1e-14));

    // B = grad(psi) x grad(phi) follows the third coordinate in the
    // right-handed (psi, phi, Z) basis.
    double B3 = gkyl_vec3_dot(
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->B),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->dual[2])
    );
    TEST_CHECK(B3 > 0.0);

    // check relationship between tangents and duals
    for (int i = 0; i < 3; ++i) {
      struct gkyl_vec3 tcart = gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[i]);

      for (int j = 0; j < 3; ++j) {
        struct gkyl_vec3 dcart = gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->dual[j]);

        // NOTE: the tangent/dual relations only hold off-axis as at
        // r=0 the coordinate system is singular
        if (rz[0] > 0) {
          // tang[i] dot dual[j] = delta_{i,j}
          double tdotd = gkyl_vec3_dot(tcart, dcart);
          if (i == j) {
            TEST_CHECK(gkyl_compare_double(tdotd, 1.0, 1e-14));
          } else {
            TEST_CHECK(gkyl_compare_double(tdotd, 0.0, 1e-14));
          }
        }
      }
    }
  }

  gkyl_mirror_grid_gen_release(geom);
  gkyl_array_release(psi);
  gkyl_position_map_release(gpm);

  return;
}

static void
test_mirror_grid_gen_wham_no_axis_psi_ho(void)
{
  test_wham(false, GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z);
}

static void
test_mirror_grid_gen_wham_with_axis_psi_ho(void)
{
  test_wham(true, GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z);
}

static void
test_mirror_grid_gen_wham_no_axis_sqrt_psi_ho(void)
{
  test_wham(false, GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z);
}

static void
test_mirror_grid_gen_wham_with_axis_sqrt_psi_ho(void)
{
  test_wham(true, GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z);
}

static void
test_quad_geom(bool include_axis, enum gkyl_mirror_grid_gen_field_line_coord fl_coord)
{
  double clower[] = {1.0e-3, 0.0, -0.75};
  double cupper[] = {0.5, 2 * M_PI, 0.75};
  int cells[] = {3, 3, 3};
  int cdim = 3;

  // computational grid
  struct gkyl_rect_grid comp_grid;
  gkyl_rect_grid_init(&comp_grid, cdim, clower, cupper, cells);

  // construct analytical psi(R,Z) on a nodal grid:
  int psi_nodes[] = {9, 17};
  struct gkyl_range psi_nodes_range;
  gkyl_range_init_from_shape(&psi_nodes_range, 2, psi_nodes);

  struct gkyl_rect_grid psi_grid;
  gkyl_rect_grid_init(&psi_grid, 2, (double[]){0.0, -1.0}, (double[]){1.0, 1.0}, psi_nodes);

  struct gkyl_array *psi = gkyl_array_new(GKYL_DOUBLE, 1, psi_nodes_range.volume);
  double dnodes[] = {1.0 / (psi_nodes[0] - 1), 2.0 / (psi_nodes[1] - 1)};

  struct gkyl_range_iter psi_nodes_iter;
  gkyl_range_iter_init(&psi_nodes_iter, &psi_nodes_range);
  while (gkyl_range_iter_next(&psi_nodes_iter)) {
    double R = 0.0 + psi_nodes_iter.idx[0] * dnodes[0];
    double Z = -1.0 + psi_nodes_iter.idx[1] * dnodes[1];

    double *pn = gkyl_array_fetch(psi, gkyl_range_idx(&psi_nodes_range, psi_nodes_iter.idx));
    pn[0] = 0.5 * (R * R) * (Z * Z + 1.0); // psi(R,Z) = 1/2*R^2*(Z^2+1.0)
  }

  struct gkyl_range node_range;
  gkyl_range_init_from_shape(&node_range, 3, (int[3]){cells[0] + 1, cells[1] + 1, cells[2] + 1});

  struct gkyl_range ext_range, range;
  int nghost[3] = {1, 1, 1};
  gkyl_create_grid_ranges(&comp_grid, nghost, &ext_range, &range);

  struct gkyl_position_map *gpm = gkyl_position_map_null_new();

  // create mirror geometry
  struct gkyl_mirror_grid_gen *geom = gkyl_mirror_grid_gen_inew(&(struct gkyl_mirror_grid_gen_inp){
    .comp_grid = &comp_grid,
    .nrange = node_range,
    .local = range,
    .global = range,

    .R = {psi_grid.lower[0], psi_grid.upper[0]},
    .Z = {psi_grid.lower[1], psi_grid.upper[1]},

    // psi(R,Z) grid size
    .nrcells = psi_grid.cells[0] - 1, // cells and not nodes
    .nzcells = psi_grid.cells[1] - 1, // cells and not nodes

    .psiRZ = psi,
    .fl_coord = fl_coord,
    .include_axis = include_axis,
    .write_psi_cubic = false,
    .psi_cubic_fname = "ctest_mirror_grid_gen_quad.gkyl",

    .position_map = gpm,
  });

  TEST_ASSERT(geom != NULL);
  TEST_CHECK(include_axis == gkyl_mirror_grid_gen_is_include_axis(geom));
  TEST_CHECK(fl_coord == gkyl_mirror_grid_gen_fl_coord(geom));

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &node_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&node_range, iter.idx);

    const double *rz = gkyl_array_cfetch(geom->nodes_rza, loc);

    double r = rz[0], z = rz[1];

    const struct gkyl_mirror_grid_gen_geom *g = gkyl_array_cfetch(geom->nodes_geom, loc);

    // construct metric tensor
    double g00 = gkyl_vec3_dot(
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[0]),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[0])
    );
    double g01 = gkyl_vec3_dot(
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[0]),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[1])
    );
    double g02 = gkyl_vec3_dot(
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[0]),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[2])
    );
    double g11 = gkyl_vec3_dot(
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[1]),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[1])
    );
    double g12 = gkyl_vec3_dot(
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[1]),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[2])
    );
    double g22 = gkyl_vec3_dot(
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[2]),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[2])
    );

    if (r > 0) {
      if (fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z) {
        // g00
        TEST_CHECK(gkyl_compare_double(1 / (SQ(r) * SQ(1 + SQ(z))), g00, 1e-14));
        // g01
        TEST_CHECK(gkyl_compare_double(0.0, g01, 1e-14));
        // g02
        TEST_CHECK(gkyl_compare_double(-SQ(r) * z / (SQ(r) * SQ(1 + SQ(z))), g02, 1e-14));

        // g11
        TEST_CHECK(gkyl_compare_double(SQ(r), g11, 1e-14));
        // g12
        TEST_CHECK(gkyl_compare_double(0.0, g12, 1e-14));

        // g22
        TEST_CHECK(gkyl_compare_double(
          (SQ(r * (1 + SQ(z))) + SQ(SQ(r) * z)) / (SQ(r) * SQ(1 + SQ(z))), g22, 1e-14
        ));
      }

      if (fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
        double psil = 0.5 * SQ(r) * (1 + SQ(z));

        // g00
        TEST_CHECK(gkyl_compare_double(4 * psil / (SQ(r) * SQ(1 + SQ(z))), g00, 1e-14));
        // g01
        TEST_CHECK(gkyl_compare_double(0.0, g01, 1e-14));
        // g02
        TEST_CHECK(
          gkyl_compare_double(-2 * sqrt(psil) * SQ(r) * z / (SQ(r) * SQ(1 + SQ(z))), g02, 1e-13)
        );

        // g11
        TEST_CHECK(gkyl_compare_double(SQ(r), g11, 1e-14));
        // g12
        TEST_CHECK(gkyl_compare_double(0.0, g12, 1e-14));

        // g22
        TEST_CHECK(gkyl_compare_double(
          (SQ(r * (1 + SQ(z))) + SQ(SQ(r) * z)) / (SQ(r) * SQ(1 + SQ(z))), g22, 1e-14
        ));
      }

      double Bmag = gkyl_vec3_len(gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->B));
      TEST_CHECK(gkyl_compare_double(sqrt(SQ(r * (1 + SQ(z))) + SQ(SQ(r) * z)) / r, Bmag, 1e-14));
    } else {
      // on-axis
      double Bmag = gkyl_vec3_len(gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->B));
      TEST_CHECK(gkyl_compare_double(1 + SQ(z), Bmag, 1e-14));
    }
  }

  gkyl_mirror_grid_gen_release(geom);
  gkyl_array_release(psi);
  gkyl_position_map_release(gpm);

  return;
}

static void
test_mirror_grid_gen_quad_geom_no_axis_psi_ho(void)
{
  test_quad_geom(false, GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z);
}

static void
test_mirror_grid_gen_quad_geom_with_axis_psi_ho(void)
{
  test_quad_geom(true, GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z);
}

static void
test_mirror_grid_gen_quad_geom_no_axis_sqrt_psi_ho(void)
{
  test_quad_geom(false, GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z);
}

static void
test_mirror_grid_gen_quad_geom_with_axis_sqrt_psi_ho(void)
{
  test_quad_geom(true, GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z);
}

struct mirror_map_ctx {
  double lower, upper, shift, scale, curvature;
};

static void mirror_map(double t, const double *xn, double *out, void *ctx)
{
  const struct mirror_map_ctx *map = ctx;
  TEST_CHECK(xn[0] >= map->lower - 1e-14 && xn[0] <= map->upper + 1e-14);
  TEST_MSG("map argument %.17g outside [%.17g, %.17g]", xn[0], map->lower, map->upper);
  out[0] = map->shift + xn[0] * (map->scale + map->curvature * xn[0]);
}

static void mirror_map_deriv(double t, const double *xn, double *out, void *ctx)
{
  const struct mirror_map_ctx *map = ctx;
  TEST_CHECK(xn[0] >= map->lower - 1e-14 && xn[0] <= map->upper + 1e-14);
  out[0] = map->scale + 2.0 * map->curvature * xn[0];
}

// Inspect bits: isfinite() can be optimized away in the normal -ffast-math build.
static void mirror_check_close(double actual, double expected, double tol, const char *name)
{
  uint64_t bits;
  memcpy(&bits, &actual, sizeof bits);
  TEST_CHECK((bits & UINT64_C(0x7ff0000000000000)) != UINT64_C(0x7ff0000000000000));
  TEST_CHECK(fabs(actual - expected) <= tol * (1.0 + fabs(expected)));
  TEST_MSG("%s: got %.17g, expected %.17g (tolerance %.3g)", name, actual, expected, tol);
}

// Quadratic in each variable, so even the boundary derivative estimates used
// to construct the bicubic interpolant reproduce this equilibrium exactly.
static struct gkyl_array *mirror_psi(double strength, double curvature)
{
  struct gkyl_range nodes;
  gkyl_range_init_from_shape(&nodes, 2, (int[]){9, 17});
  struct gkyl_array *psi = gkyl_array_new(GKYL_DOUBLE, 1, nodes.volume);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &nodes);
  while (gkyl_range_iter_next(&iter)) {
    double radius = iter.idx[0] / 8.0, height = -1.0 + iter.idx[1] / 8.0;
    double *value = gkyl_array_fetch(psi, gkyl_range_idx(&nodes, iter.idx));
    value[0] = 0.5 * strength * radius * radius * (1.0 + curvature * height * height);
  }
  return psi;
}

static int mirror_nodes_per_cell(int kind, int dim, int order)
{
  return kind == -2 || (kind >= 0 && dim != kind) ? order + 1 : kind == -1 ? order : 1;
}

static struct gkyl_range
mirror_node_range(const struct gkyl_range *local, int kind, int order, const int *lower)
{
  int upper[3];
  for (int dim = 0; dim < 3; ++dim) {
    bool quad = kind == -2 || (kind >= 0 && dim != kind);
    upper[dim] = lower[dim] +
                 mirror_nodes_per_cell(kind, dim, order) * gkyl_range_shape(local, dim) -
                 (quad ? 1 : 0);
  }
  struct gkyl_range nodes;
  gkyl_range_init(&nodes, 3, lower, upper);
  return nodes;
}

static struct gkyl_mirror_grid_gen *
mirror_generate(const struct gkyl_mirror_grid_gen_inp *inp, int kind)
{
  return kind == -2 ? gkyl_mirror_grid_gen_int_inew(inp) :
         kind == -1 ? gkyl_mirror_grid_gen_inew(inp) :
                      gkyl_mirror_grid_gen_surf_inew(inp);
}

// Independent P1/P2 nodes, without using the production coordinate helper or
// its Gauss table. 'node' is measured from the global lower cell boundary.
static double
mirror_coordinate(double lower, double upper, int cells, int node, int count, bool quad)
{
  double cell = (double)node / count;
  if (quad) {
    double gauss = count == 2 ? (node % count == 0 ? -1.0 : 1.0) / sqrt(3.0) :
                                (node % count - 1) * sqrt(3.0 / 5.0);
    cell = node / count + 0.5 * (1.0 + gauss);
  }
  return lower + (upper - lower) * cell / cells;
}

static void mirror_check_geom(
  const struct gkyl_mirror_grid_gen_geom *actual, const struct gkyl_mirror_grid_gen_geom *expected,
  double tol
)
{
  for (int dim = 0; dim < 3; ++dim) {
    for (int comp = 0; comp < 3; ++comp) {
      mirror_check_close(actual->tang[dim].x[comp], expected->tang[dim].x[comp], tol, "tangent");
      mirror_check_close(actual->dual[dim].x[comp], expected->dual[dim].x[comp], tol, "dual");
    }
    mirror_check_close(actual->B.x[dim], expected->B.x[dim], tol, "B");
    mirror_check_close(actual->curlbhat.x[dim], expected->curlbhat.x[dim], tol, "curl(bhat)");
  }
  mirror_check_close(actual->Jc, expected->Jc, tol, "Jacobian");
}

static void mirror_check_exact(
  const struct gkyl_mirror_grid_gen_inp *inp, const struct gkyl_mirror_grid_gen *geom, int kind,
  int order, const struct mirror_map_ctx *maps, double strength, double curvature
)
{
  TEST_CHECK(gkyl_mirror_grid_gen_fl_coord(geom) == inp->fl_coord);
  TEST_CHECK(gkyl_mirror_grid_gen_is_include_axis(geom) == inp->include_axis);
  bool sqrt_psi = inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &inp->nrange);
  while (gkyl_range_iter_next(&iter)) {
    double mapped[3], slope[3];
    int global_node[3];
    for (int dim = 0; dim < 3; ++dim) {
      int count = mirror_nodes_per_cell(kind, dim, order);
      bool quad = kind == -2 || (kind >= 0 && dim != kind);
      global_node[dim] = iter.idx[dim] - inp->nrange.lower[dim] +
                         count * (inp->local.lower[dim] - inp->global.lower[dim]);
      double xc = mirror_coordinate(
        maps[dim].lower, maps[dim].upper, inp->comp_grid->cells[dim], global_node[dim], count, quad
      );
      mapped[dim] = maps[dim].shift + xc * (maps[dim].scale + maps[dim].curvature * xc);
      slope[dim] = maps[dim].scale + 2.0 * maps[dim].curvature * xc;
    }
    bool axis = inp->include_axis && (kind == -1 || kind == 0) && global_node[0] == 0;
    double flux = axis ? 0.0 : sqrt_psi ? mapped[0] * mapped[0] : mapped[0];
    double height = mapped[2], factor = 1.0 + curvature * height * height;
    double radius = sqrt(2.0 * flux / (strength * factor));
    double radial_scale = slope[0] * (sqrt_psi ? 2.0 * mapped[0] : 1.0);
    long loc = gkyl_range_idx(&inp->nrange, iter.idx);
    const double *rza = gkyl_array_cfetch(geom->nodes_rza, loc);
    const double *psi = gkyl_array_cfetch(geom->nodes_psi, loc);
    const struct gkyl_mirror_grid_gen_geom *actual = gkyl_array_cfetch(geom->nodes_geom, loc);
    mirror_check_close(rza[0], radius, 2e-11, "R");
    mirror_check_close(rza[1], height, 2e-14, "Z");
    mirror_check_close(rza[2], mapped[1], 2e-14, "phi");
    mirror_check_close(psi[0], flux, 2e-14, "psi");

    struct gkyl_mirror_grid_gen_geom expected = {0};
    expected.B.x[0] = -strength * curvature * radius * height;
    expected.B.x[2] = strength * factor;
    expected.tang[1].x[1] = slope[1];
    expected.tang[2].x[2] = slope[2];
    expected.dual[2].x[2] = 1.0 / slope[2];
    expected.Jc = radial_scale * slope[1] * slope[2] / (strength * factor);
    if (axis) {
      // The public representation deliberately uses finite placeholder bases
      // at the cylindrical singularity, while B and J use their true limits.
      expected.tang[0].x[0] = slope[0];
      expected.dual[0].x[0] = 1.0 / slope[0];
      expected.dual[1].x[1] = 1.0 / slope[1];
    } else {
      expected.tang[0].x[0] = radial_scale / (strength * radius * factor);
      expected.tang[2].x[0] = -curvature * radius * height * slope[2] / factor;
      expected.dual[0].x[0] = strength * radius * factor / radial_scale;
      expected.dual[0].x[2] = strength * curvature * radius * radius * height / radial_scale;
      expected.dual[1].x[1] = 1.0 / (radius * radius * slope[1]);
      // Closed form of (d_Z b_R - d_R b_Z)/R. Independent of field strength.
      double norm2 = factor * factor + curvature * curvature * radius * radius * height * height;
      expected.curlbhat.x[1] =
        -curvature * factor * (1.0 - 2.0 * curvature * height * height) / pow(norm2, 1.5);

      // Check handedness and the Clebsch field directly in Cartesian space,
      // including the chain-rule factors for mapped psi/sqrt(psi) and phi.
      struct gkyl_vec3 tangent[3], dual[3];
      for (int dim = 0; dim < 3; ++dim) {
        tangent[dim] = gkyl_vec3_polar_con_to_cart(rza[0], rza[2], actual->tang[dim]);
        dual[dim] = gkyl_vec3_polar_con_to_cart(rza[0], rza[2], actual->dual[dim]);
      }
      double jac = gkyl_vec3_triple(tangent[0], tangent[1], tangent[2]);
      TEST_CHECK(jac > 0.0);
      mirror_check_close(actual->Jc, jac, 2e-11, "right-handed Jacobian");
      struct gkyl_vec3 field = gkyl_vec3_polar_con_to_cart(rza[0], rza[2], actual->B);
      struct gkyl_vec3 clebsch =
        gkyl_vec3_scale(radial_scale * slope[1], gkyl_vec3_cross(dual[0], dual[1]));
      for (int dim = 0; dim < 3; ++dim) {
        // Allow the radial root-finder error in the sqrt(psi) scale factor.
        mirror_check_close(field.x[dim], clebsch.x[dim], 2e-9, "grad(psi) x grad(phi)");
      }
      double B3 = gkyl_vec3_dot(field, dual[2]);
      TEST_CHECK(B3 > 0.0);
      mirror_check_close(
        B3, gkyl_vec3_len(field) / gkyl_vec3_len(tangent[2]), 2e-11, "parallel field"
      );
    }
    mirror_check_geom(actual, &expected, 2e-9);
  }
}

static void mirror_check_partition(
  const struct gkyl_mirror_grid_gen_inp *inp, const struct gkyl_mirror_grid_gen *geom,
  const struct gkyl_range *full_nodes, const struct gkyl_mirror_grid_gen *full, int kind, int order
)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &inp->nrange);
  while (gkyl_range_iter_next(&iter)) {
    int index[3];
    for (int dim = 0; dim < 3; ++dim) {
      index[dim] =
        full_nodes->lower[dim] + iter.idx[dim] - inp->nrange.lower[dim] +
        mirror_nodes_per_cell(kind, dim, order) * (inp->local.lower[dim] - inp->global.lower[dim]);
    }
    long loc = gkyl_range_idx(&inp->nrange, iter.idx), ref_loc = gkyl_range_idx(full_nodes, index);
    const double *rza = gkyl_array_cfetch(geom->nodes_rza, loc);
    const double *ref_rza = gkyl_array_cfetch(full->nodes_rza, ref_loc);
    for (int dim = 0; dim < 3; ++dim) {
      mirror_check_close(rza[dim], ref_rza[dim], 0.0, "partition coordinate");
    }
    const double *psi = gkyl_array_cfetch(geom->nodes_psi, loc);
    const double *ref_psi = gkyl_array_cfetch(full->nodes_psi, ref_loc);
    mirror_check_close(psi[0], ref_psi[0], 0.0, "partition psi");
    mirror_check_geom(
      gkyl_array_cfetch(geom->nodes_geom, loc), gkyl_array_cfetch(full->nodes_geom, ref_loc), 0.0
    );
  }
}

static void mirror_exact_cases(int order, bool sqrt_psi)
{
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(
    &grid, 3, (double[]){0.025, -0.4, -0.75}, (double[]){0.24, 0.7, 0.85}, (int[]){3, 2, 5}
  );
  // Nonstandard origins expose accidental use of absolute cell/node indices.
  struct gkyl_range global;
  gkyl_range_init(&global, 3, (int[]){4, -2, 7}, (int[]){6, -1, 11});
  struct gkyl_array *psi = mirror_psi(1.3, 0.7);
  for (int include_axis = 0; include_axis < 2; ++include_axis) {
    for (int map_mode = 0; map_mode < 4; ++map_mode) {
      struct mirror_map_ctx maps[3];
      struct gkyl_position_map *pmap = gkyl_position_map_null_new();
      pmap->use_map_derivs = map_mode != 1;
      for (int dim = 0; dim < 3; ++dim) {
        maps[dim] =
          (struct mirror_map_ctx){.lower = grid.lower[dim], .upper = grid.upper[dim], .scale = 1.0};
        if (dim == 0 && sqrt_psi) {
          maps[dim].lower = include_axis ? 0.0 : sqrt(grid.lower[0]);
          maps[dim].upper = sqrt(grid.upper[0]);
        }
        if (map_mode) {
          maps[dim].shift = dim == 0 ? (include_axis ? 0.0 : 0.015) : 0.02;
          maps[dim].scale = dim == 1 ? 1.2 : 0.7;
          maps[dim].curvature = 0.1;
        }
        pmap->maps[dim] = mirror_map;
        pmap->map_derivs[dim] = map_mode == 3 && dim == 1 ? NULL : mirror_map_deriv;
        pmap->ctxs[dim] = &maps[dim];
      }
      for (int kind = -2; kind < 3; ++kind) {
        struct gkyl_range full_nodes = mirror_node_range(&global, kind, order, (int[]){-3, 5, 2});
        struct gkyl_mirror_grid_gen *full = NULL;
        // Full grid followed by an uneven 2x2x2 decomposition, including
        // single-cell radial/angular pieces and shared corners and faces.
        for (int part = -1; part < 8; ++part) {
          int lower[3], upper[3];
          for (int dim = 0; dim < 3; ++dim) {
            int cut = global.lower[dim] + grid.cells[dim] / 2;
            lower[dim] = part < 0 || !(part & (1 << dim)) ? global.lower[dim] : cut;
            upper[dim] = part < 0 || (part & (1 << dim)) ? global.upper[dim] : cut - 1;
          }
          struct gkyl_range local;
          gkyl_sub_range_init(&local, &global, lower, upper);
          struct gkyl_range nodes =
            part < 0 ? full_nodes : mirror_node_range(&local, kind, order, (int[]){4, -3, 11});
          struct gkyl_mirror_grid_gen_inp inp = {
            .comp_grid = &grid,
            .local = local,
            .global = global,
            .nrange = nodes,
            .dir = kind,
            .position_map = pmap,
            .include_axis = include_axis,
            .fl_coord = sqrt_psi ? GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z :
                                   GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z,
            .R = {0.0, 1.0},
            .Z = {-1.0, 1.0},
            .nrcells = 8,
            .nzcells = 16,
            .psiRZ = psi
          };
          TEST_CASE_(
            "order=%d sqrt_psi=%d axis=%d map=%d kind=%d part=%d", order, sqrt_psi, include_axis,
            map_mode, kind, part
          );
          struct gkyl_mirror_grid_gen *geom = mirror_generate(&inp, kind);
          TEST_ASSERT(geom != NULL);
          mirror_check_exact(&inp, geom, kind, order, maps, 1.3, 0.7);
          if (part < 0) {
            full = geom;
          } else {
            mirror_check_partition(&inp, geom, &full_nodes, full, kind, order);
            gkyl_mirror_grid_gen_release(geom);
          }
        }
        gkyl_mirror_grid_gen_release(full);
      }
      gkyl_position_map_release(pmap);
    }
  }
  gkyl_array_release(psi);
}

static void test_mirror_exact_psi_p1_ho(void)
{
  mirror_exact_cases(1, false);
}
static void test_mirror_exact_sqrt_psi_p1_ho(void)
{
  mirror_exact_cases(1, true);
}
static void test_mirror_exact_psi_p2_ho(void)
{
  mirror_exact_cases(2, false);
}
static void test_mirror_exact_sqrt_psi_p2_ho(void)
{
  mirror_exact_cases(2, true);
}

static void test_mirror_curl_scaling_ho(void)
{
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(
    &grid, 3, (double[]){0.03, -0.3, -0.8}, (double[]){0.12, 0.5, 0.8}, (int[]){2, 1, 3}
  );
  struct gkyl_range local, ext;
  gkyl_create_grid_ranges(&grid, (int[]){0, 0, 0}, &ext, &local);
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();
  struct mirror_map_ctx maps[3];
  for (int dim = 0; dim < 3; ++dim) {
    maps[dim] =
      (struct mirror_map_ctx){.lower = grid.lower[dim], .upper = grid.upper[dim], .scale = 1.0};
  }
  // A constant axial field must have zero curl(bhat) at every radius.
  // For a curved field, multiplying psi and its flux labels by the same
  // positive constant leaves both the node locations and curl(bhat) unchanged.
  double strengths[] = {0.3, 1.0, 7.0};
  for (int curved = 0; curved < 2; ++curved) {
    for (int scale = 0; scale < 3; ++scale) {
      double strength = strengths[scale];
      struct gkyl_array *psi = mirror_psi(strength, curved);
      struct gkyl_rect_grid scaled_grid;
      gkyl_rect_grid_init(
        &scaled_grid, 3, (double[]){strength * grid.lower[0], grid.lower[1], grid.lower[2]},
        (double[]){strength * grid.upper[0], grid.upper[1], grid.upper[2]}, grid.cells
      );
      maps[0].lower = scaled_grid.lower[0];
      maps[0].upper = scaled_grid.upper[0];
      for (int kind = -2; kind < 3; ++kind) {
        struct gkyl_mirror_grid_gen_inp inp = {
          .comp_grid = &scaled_grid,
          .local = local,
          .global = local,
          .nrange = mirror_node_range(&local, kind, 1, (int[]){0, 0, 0}),
          .dir = kind,
          .position_map = pmap,
          .R = {0.0, 1.0},
          .Z = {-1.0, 1.0},
          .nrcells = 8,
          .nzcells = 16,
          .psiRZ = psi
        };
        struct gkyl_mirror_grid_gen *geom = mirror_generate(&inp, kind);
        TEST_ASSERT(geom != NULL);
        mirror_check_exact(&inp, geom, kind, 1, maps, strength, curved);
        gkyl_mirror_grid_gen_release(geom);
      }
      gkyl_array_release(psi);
    }
  }
  gkyl_position_map_release(pmap);
}

// Capture only expected rejection diagnostics; keep production error reporting.
static struct gkyl_mirror_grid_gen *
mirror_generate_rejected(const struct gkyl_mirror_grid_gen_inp *inp, int kind)
{
  FILE *diagnostics = tmpfile();
  TEST_ASSERT(diagnostics != NULL);
  fflush(stderr);
  int saved_stderr = dup(STDERR_FILENO);
  if (saved_stderr < 0) {
    fclose(diagnostics);
    TEST_ASSERT(saved_stderr >= 0);
  }
  int redirected = dup2(fileno(diagnostics), STDERR_FILENO);
  struct gkyl_mirror_grid_gen *geom = redirected >= 0 ? mirror_generate(inp, kind) : NULL;
  fflush(stderr);
  int restored = dup2(saved_stderr, STDERR_FILENO);
  close(saved_stderr);

  rewind(diagnostics);
  char message[256];
  size_t count = fread(message, 1, sizeof(message) - 1, diagnostics);
  message[count] = '\0';
  fclose(diagnostics);
  TEST_ASSERT(redirected >= 0 && restored >= 0);
  TEST_CHECK(strcmp(message, "gkyl_mirror_grid_gen failed to generate a grid\n") == 0);
  TEST_MSG("Expected grid rejection diagnostic; got: %s", message);
  return geom;
}

static void test_mirror_root_bounds_ho(void)
{
  struct gkyl_array *psi = mirror_psi(1.0, 0.0);
  struct gkyl_position_map *pmap = gkyl_position_map_null_new();
  // psi = R^2/2 on 0 <= R <= 1 admits only 0 <= psi <= 1/2.
  // Check both sides of the bracket, and failure after some valid nodes
  // have already been generated, through every public constructor.
  double lower_flux[] = {1.5, -0.2, 0.1, 0.1};
  double upper_flux[] = {2.0, -0.1, 0.7, 0.5};
  for (int trial = 0; trial < 4; ++trial) {
    struct gkyl_rect_grid grid;
    gkyl_rect_grid_init(
      &grid, 3, (double[]){lower_flux[trial], -0.1, -0.5}, (double[]){upper_flux[trial], 0.1, 0.5},
      (int[]){3, 1, 1}
    );
    struct gkyl_range local, ext;
    gkyl_create_grid_ranges(&grid, (int[]){0, 0, 0}, &ext, &local);
    for (int kind = -2; kind < 3; ++kind) {
      struct gkyl_mirror_grid_gen_inp inp = {
        .comp_grid = &grid,
        .local = local,
        .global = local,
        .nrange = mirror_node_range(&local, kind, 1, (int[]){0, 0, 0}),
        .dir = kind,
        .position_map = pmap,
        .R = {0.0, 1.0},
        .Z = {-1.0, 1.0},
        .nrcells = 8,
        .nzcells = 16,
        .psiRZ = psi
      };
      struct gkyl_mirror_grid_gen *geom = trial < 3 ? mirror_generate_rejected(&inp, kind) :
                                                      mirror_generate(&inp, kind);
      TEST_CHECK((geom != NULL) == (trial == 3));
      TEST_MSG("trial=%d kind=%d: out-of-domain flux must fail; R=1 must be accepted", trial, kind);
      if (geom) {
        if (trial == 3) {
          struct mirror_map_ctx maps[3];
          for (int dim = 0; dim < 3; ++dim) {
            maps[dim] = (struct mirror_map_ctx
            ){.lower = grid.lower[dim], .upper = grid.upper[dim], .scale = 1.0};
          }
          mirror_check_exact(&inp, geom, kind, 1, maps, 1.0, 0.0);
        }
        gkyl_mirror_grid_gen_release(geom);
      }
    }
  }
  gkyl_position_map_release(pmap);
  gkyl_array_release(psi);
}

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
            // B_R = -R Z, B_Z = 1+Z^2, B_phi = 0.
            double radius = rz[0], height = rz[1], bz = 1.0 + height * height,
                   br = -radius * height;
            double bmag = sqrt(br * br + bz * bz), dbdr = radius * height * height / bmag;
            double dbdz = (radius * radius * height + 2.0 * height * bz) / bmag;
            double curl_phi = -radius / bmag + (dbdr * bz - dbdz * br) / (bmag * bmag);
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
  {"mirror_identity", test_mirror_identity},
  {"mirror_affine", test_mirror_affine},
  {"mirror_nonlinear_numeric", test_mirror_nonlinear_numeric},
  {"mirror_nonlinear_analytic", test_mirror_nonlinear_analytic},
  {"mirror_sqrt_psi", test_mirror_sqrt_psi},
  {"mirror_exact_psi_p1_ho", test_mirror_exact_psi_p1_ho},
  {"mirror_exact_sqrt_psi_p1_ho", test_mirror_exact_sqrt_psi_p1_ho},
  {"mirror_exact_psi_p2_ho", test_mirror_exact_psi_p2_ho},
  {"mirror_exact_sqrt_psi_p2_ho", test_mirror_exact_sqrt_psi_p2_ho},
  {"mirror_curl_scaling_ho", test_mirror_curl_scaling_ho},
  {"mirror_root_bounds_ho", test_mirror_root_bounds_ho},

  {"mirror_grid_gen_wham_no_axis_psi_ho", test_mirror_grid_gen_wham_no_axis_psi_ho},
  {"mirror_grid_gen_wham_with_axis_psi_ho", test_mirror_grid_gen_wham_with_axis_psi_ho},

  {"mirror_grid_gen_wham_no_axis_sqrt_psi_ho", test_mirror_grid_gen_wham_no_axis_sqrt_psi_ho},
  {"mirror_grid_gen_wham_with_axis_sqrt_psi_ho", test_mirror_grid_gen_wham_with_axis_sqrt_psi_ho},

  {"mirror_grid_gen_quad_geom_no_axis_psi_ho", test_mirror_grid_gen_quad_geom_no_axis_psi_ho},
  {"mirror_grid_gen_quad_geom_with_axis_psi_ho", test_mirror_grid_gen_quad_geom_with_axis_psi_ho},

  {"mirror_grid_gen_quad_geom_no_axis_sqrt_psi_ho",
   test_mirror_grid_gen_quad_geom_no_axis_sqrt_psi_ho},
  {"mirror_grid_gen_quad_geom_with_axis_sqrt_psi_ho",
   test_mirror_grid_gen_quad_geom_with_axis_sqrt_psi_ho},
  {NULL, NULL}
};
