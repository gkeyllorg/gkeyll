#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_math.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_position_map.h>

static inline double SQ(double x) { return x*x; };

static void
test_wham(bool include_axis, enum gkyl_mirror_grid_gen_field_line_coord fl_coord)
{
  double clower[] = { 2.0e-6, 0.0, -2.0 };
  double cupper[] = { 3.0e-3, 2*M_PI, 2.0 };
  int cells[] = { 10, 16, 32 };
  int cdim = 3;

  const char *fname = "gyrokinetic/data/unit/wham_hires.geqdsk_psi.gkyl";
  
  // computational grid
  struct gkyl_rect_grid comp_grid;
  gkyl_rect_grid_init(&comp_grid, cdim, clower, cupper, cells);

  if (!gkyl_check_file_exists(fname)) {
    fprintf(stderr, "Unable to find file %s!\n", fname);
    goto cleanup;
  }
  
  // read psi(R,Z) from file
  struct gkyl_rect_grid psi_grid;
  struct gkyl_array *psi = gkyl_grid_array_new_from_file(&psi_grid, fname);

  struct gkyl_range node_range;
  gkyl_range_init_from_shape(&node_range, 3, (int[3]) { cells[0]+1, cells[1]+1, cells[2]+1 });

  struct gkyl_range ext_range, range;
  int nghost[3] = {1,1,1};
  gkyl_create_grid_ranges(&comp_grid, nghost, &ext_range, &range);

  struct gkyl_position_map *gpm = gkyl_position_map_null_new();

  // create mirror geometry
  struct gkyl_mirror_grid_gen *geom =
    gkyl_mirror_grid_gen_inew(&(struct gkyl_mirror_grid_gen_inp) {
        .comp_grid = &comp_grid,
        .nrange = node_range,
        .local = range,
        .global = range,

        .R = { psi_grid.lower[0], psi_grid.upper[0] },
        .Z = { psi_grid.lower[1], psi_grid.upper[1] },
        
        // psi(R,Z) grid size
        .nrcells = psi_grid.cells[0]-1, // cells and not nodes
        .nzcells = psi_grid.cells[1]-1, // cells and not nodes

        .psiRZ = psi,
        .fl_coord = fl_coord,
        .include_axis = include_axis,
        .write_psi_cubic = false,

        .position_map = gpm,
      }
    );

  TEST_CHECK( include_axis == gkyl_mirror_grid_gen_is_include_axis(geom) );
  TEST_CHECK( fl_coord == gkyl_mirror_grid_gen_fl_coord(geom) );  
  

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &node_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&node_range, iter.idx);

    const double *rz =
      gkyl_array_cfetch(geom->nodes_rza, loc);

    const struct gkyl_mirror_grid_gen_geom *g =
      gkyl_array_cfetch(geom->nodes_geom, loc);

    // check Jacobian
    double Jac = gkyl_vec3_triple(
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[0]),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[1]),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[2])
    );

    if (rz[0] > 0)
      TEST_CHECK( gkyl_compare_double(Jac, g->Jc, 1e-14) );

    // check C = Jc*Bmag/sqrt(g33)
    double g33 =
      gkyl_vec3_dot(
        gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[2]),
        gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[2])
      );
    double Bmag = gkyl_vec3_len(
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->B)
    );
    
    if (fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z)
      TEST_CHECK( gkyl_compare_double(Bmag*g->Jc/sqrt(g33), 1.0, 1e-14) );

    // check B only points in the parallel direction
    // ... B^1 = 0
    double B1 = gkyl_vec3_dot(
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->B),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->dual[0])
    );
    TEST_CHECK( gkyl_compare_double(B1, 0.0, 1e-14) );

    // ... B^2 = 0
    double B2 = gkyl_vec3_dot(
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->B),
      gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->dual[1])
    );
    TEST_CHECK( gkyl_compare_double(B1, 0.0, 1e-14) );    
    
    // check relationship between tangents and duals
    for (int i=0; i<3; ++i) {
      struct gkyl_vec3 tcart = gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->tang[i]);
      
      for (int j=0; j<3; ++j) {
        struct gkyl_vec3 dcart = gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->dual[j]);

        // NOTE: the tangent/dual relations only hold off-axis as at
        // r=0 the coordinate system is singular
        if (rz[0] > 0) {
          // tang[i] dot dual[j] = delta_{i,j}
          double tdotd = gkyl_vec3_dot(tcart, dcart);
          if (i == j)
            TEST_CHECK( gkyl_compare_double(tdotd, 1.0, 1e-14) );
          else
            TEST_CHECK( gkyl_compare_double(tdotd, 0.0, 1e-14) );
        }
      }
    }
  }

  gkyl_mirror_grid_gen_release(geom);
  gkyl_array_release(psi);
  gkyl_position_map_release(gpm);

  cleanup:

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
  double clower[] = { 1.0e-3, 0.0, -0.75 };
  double cupper[] = { 0.5, 2*M_PI, 0.75 };
  int cells[] = { 3, 3, 3 };
  int cdim = 3;
  
  // computational grid
  struct gkyl_rect_grid comp_grid;
  gkyl_rect_grid_init(&comp_grid, cdim, clower, cupper, cells);

  // construct analytical psi(R,Z) on a nodal grid:
  int psi_nodes[] = { 9, 17 };
  struct gkyl_range psi_nodes_range;
  gkyl_range_init_from_shape(&psi_nodes_range, 2, psi_nodes);
  
  struct gkyl_rect_grid psi_grid;
  gkyl_rect_grid_init(&psi_grid, 2,
    (double[]) { 0.0, -1.0 },
    (double[]) { 1.0, 1.0 },
    psi_nodes
  );
    
  struct gkyl_array *psi = gkyl_array_new(GKYL_DOUBLE, 1, psi_nodes_range.volume);
  double dnodes[] = { 1.0/(psi_nodes[0]-1), 2.0/(psi_nodes[1]-1) };
  
  struct gkyl_range_iter psi_nodes_iter;
  gkyl_range_iter_init(&psi_nodes_iter, &psi_nodes_range);
  while (gkyl_range_iter_next(&psi_nodes_iter)) {
    double R = 0.0 + psi_nodes_iter.idx[0]*dnodes[0];
    double Z = -1.0 + psi_nodes_iter.idx[1]*dnodes[1];

    double *pn = gkyl_array_fetch(psi, gkyl_range_idx(&psi_nodes_range, psi_nodes_iter.idx));
    pn[0] = 0.5*(R*R)*(Z*Z+1.0); // psi(R,Z) = 1/2*R^2*(Z^2+1.0)
  }

  struct gkyl_range node_range;
  gkyl_range_init_from_shape(&node_range, 3, (int[3]) { cells[0]+1, cells[1]+1, cells[2]+1 });

  struct gkyl_range ext_range, range;
  int nghost[3] = { 1,1,1};
  gkyl_create_grid_ranges(&comp_grid, nghost, &ext_range, &range);

  struct gkyl_position_map *gpm = gkyl_position_map_null_new();

  // create mirror geometry
  struct gkyl_mirror_grid_gen *geom =
    gkyl_mirror_grid_gen_inew(&(struct gkyl_mirror_grid_gen_inp) {
        .comp_grid = &comp_grid,
        .nrange = node_range,
        .local = range,
        .global = range,
        
        .R = { psi_grid.lower[0], psi_grid.upper[0] },
        .Z = { psi_grid.lower[1], psi_grid.upper[1] },
        
        // psi(R,Z) grid size
        .nrcells = psi_grid.cells[0]-1, // cells and not nodes
        .nzcells = psi_grid.cells[1]-1, // cells and not nodes

        .psiRZ = psi,
        .fl_coord = fl_coord,
        .include_axis = include_axis,
        .write_psi_cubic = false,
        .psi_cubic_fname = "ctest_mirror_grid_gen_quad.gkyl",

        .position_map = gpm,
      }
    );

  TEST_CHECK( include_axis == gkyl_mirror_grid_gen_is_include_axis(geom) );
  TEST_CHECK( fl_coord == gkyl_mirror_grid_gen_fl_coord(geom) );  

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &node_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&node_range, iter.idx);

    const double *rz =
      gkyl_array_cfetch(geom->nodes_rza, loc);

    double r = rz[0], z = rz[1];

    const struct gkyl_mirror_grid_gen_geom *g =
      gkyl_array_cfetch(geom->nodes_geom, loc);

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

    if (r>0) {
      if (fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z) {
        // g00
        TEST_CHECK ( gkyl_compare_double(1/(SQ(r)*SQ(1+SQ(z))), g00, 1e-14) );
        // g01
        TEST_CHECK ( gkyl_compare_double(0.0, g01, 1e-14) );
        // g02
        TEST_CHECK ( gkyl_compare_double(-SQ(r)*z/(SQ(r)*SQ(1+SQ(z))), g02, 1e-14) );
        
        // g11
        TEST_CHECK ( gkyl_compare_double(SQ(r), g11, 1e-14) );
        // g12
        TEST_CHECK ( gkyl_compare_double(0.0, g12, 1e-14) );
        
        // g22
        TEST_CHECK ( gkyl_compare_double( (SQ(r*(1+SQ(z)))+SQ(SQ(r)*z))/(SQ(r)*SQ(1+SQ(z))),
            g22, 1e-14)
        );
    }
      
      if (fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
        double psil = 0.5*SQ(r)*(1+SQ(z));
        
        // g00
        TEST_CHECK ( gkyl_compare_double(4*psil/(SQ(r)*SQ(1+SQ(z))), g00, 1e-14) );
        // g01
        TEST_CHECK ( gkyl_compare_double(0.0, g01, 1e-14) );
        // g02
        TEST_CHECK ( gkyl_compare_double(-2*sqrt(psil)*SQ(r)*z/(SQ(r)*SQ(1+SQ(z))), g02, 1e-13) );

        // g11
        TEST_CHECK ( gkyl_compare_double(SQ(r), g11, 1e-14) );
        // g12
        TEST_CHECK ( gkyl_compare_double(0.0, g12, 1e-14) );

        // g22
        TEST_CHECK ( gkyl_compare_double( (SQ(r*(1+SQ(z)))+SQ(SQ(r)*z))/(SQ(r)*SQ(1+SQ(z))),
            g22, 1e-14)
        );      
      }

      double Bmag = gkyl_vec3_len(
        gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->B)
      );    
      TEST_CHECK( gkyl_compare_double(sqrt(SQ(r*(1+SQ(z)))+SQ(SQ(r)*z))/r, Bmag, 1e-14) );
    }
    else {
      // on-axis
      double Bmag = gkyl_vec3_len(
        gkyl_vec3_polar_con_to_cart(rz[0], 0.0, g->B)
      );
      TEST_CHECK( gkyl_compare_double(1+SQ(z), Bmag, 1e-14) );
    }
    
  }

  gkyl_mirror_grid_gen_release(geom);
  gkyl_array_release(psi);
  gkyl_position_map_release(gpm);

  cleanup:

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

TEST_LIST = {
  { "mirror_grid_gen_wham_no_axis_psi_ho", test_mirror_grid_gen_wham_no_axis_psi_ho },
  { "mirror_grid_gen_wham_with_axis_psi_ho", test_mirror_grid_gen_wham_with_axis_psi_ho },

  { "mirror_grid_gen_wham_no_axis_sqrt_psi_ho", test_mirror_grid_gen_wham_no_axis_sqrt_psi_ho },
  { "mirror_grid_gen_wham_with_axis_sqrt_psi_ho", test_mirror_grid_gen_wham_with_axis_sqrt_psi_ho },

  { "mirror_grid_gen_quad_geom_no_axis_psi_ho", test_mirror_grid_gen_quad_geom_no_axis_psi_ho},
  { "mirror_grid_gen_quad_geom_with_axis_psi_ho", test_mirror_grid_gen_quad_geom_with_axis_psi_ho},

  { "mirror_grid_gen_quad_geom_no_axis_sqrt_psi_ho", test_mirror_grid_gen_quad_geom_no_axis_sqrt_psi_ho },
  { "mirror_grid_gen_quad_geom_with_axis_sqrt_psi_ho", test_mirror_grid_gen_quad_geom_with_axis_sqrt_psi_ho },
  { NULL, NULL },
};
