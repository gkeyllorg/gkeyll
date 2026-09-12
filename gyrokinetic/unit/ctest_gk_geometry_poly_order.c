// Gate 3 regressions: polynomial order in the tokamak nodal construction.
//
// These tests were written BEFORE the fixes they pin, which is gate 3's rule:
// each one FAILED on the unfixed library and passes once its defect is fixed.
//
// READ THIS FIRST -- the p2 tests do not currently run.
//
// Measured 2026-09-09: p2 tokamak geometry cannot be completed by this library
// at all, because construction forms bmag_inv = 1/|B| through
// gkyl_dg_inv_op_range and there is no ser_3x_p2_inv kernel anywhere in the
// tree. That is generated code in core, not geometry code, so no geometry fix
// removes it. The library therefore REFUSES p2 up front, by name.
//
// Consequence, stated plainly: the dx_fact fix (defect 3.1) is real, and was
// observed failing before the fix -- but it has NO live test coverage while
// that kernel is missing, because the only way to observe it is to build a p2
// geometry. The three p2 tests below ask the library the same question the
// library asks itself, print a loud SKIPPED line naming what is uncovered, and
// begin enforcing automatically the day the kernel lands. Nothing has to be
// remembered and re-enabled.
//
// The fixture is straight_cylinder.geqdsk, an analytic equilibrium with
// psi = 0.25*R^2, so every node position has a closed form and the assertions
// are exact rather than golden values. It carries no vessel outline in the
// file, so the absence is declared -- silence never disables the wall.
//
// What each test pins:
//
//   p1_corner_node_spacing   Positive control. p1 must keep working; it is the
//                            only order this campaign has ever exercised.
//   p2_corner_node_spacing   Defect 3.1. `dx_fact` is written
//                            `poly_order == 1.0/poly_order` -- a comparison,
//                            not a division. At p1 it yields 1.0, which is the
//                            correct factor, so p1 hides it. At p2 it yields
//                            0.0, so dpsi and dalpha become zero and every
//                            node collapses onto the block's lower bound.
//   p2_fd_steps_nonzero      Defect 3.1, second consequence. The same zero
//                            factor reaches dzc[], the finite-difference steps
//                            the metric calculation divides by.
//   p2_interior_quadrature   Defect 3.1b. Interior nodes are each cell's
//                            Gauss-Legendre points, and p2 carries three per
//                            direction. calc_running_coord walks a hard-coded
//                            two-step cycle, which IS the 2-point rule, so it
//                            is right at p1 and cannot express three points at
//                            any stride. NOT fixed by correcting dx_fact.
//   unsupported_poly_order   Defect 3.3. The nodal shape is computed for
//                            poly_order 1 and 2 with no else and no assert,
//                            over an uninitialised array. An unsupported order
//                            must be refused explicitly, not sized from stack
//                            garbage.
//
// Note on theta: theta node positions are NOT affected by defect 3.1. The
// theta loop derives its own spacing as arcL_tot/(poly_order*cells), which is
// already order-correct. Only psi and alpha collapse. The assertions below
// reflect that rather than the whole-grid collapse first predicted.
//
// The four tests over the corner grid are independent of the interior one: the
// corner nodes are UNIFORM at poly_order*cells+1, while the interior nodes are
// each cell's quadrature points. Those are different constructions, and at p2
// they fail for different reasons.

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_position_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_tok_geo.h>

enum { PSI_IDX, AL_IDX, TH_IDX };

static const double psi_min = 0.1, psi_max = 0.2;
static const double alpha_min = -1.0, alpha_max = 1.0;
static const int cells[3] = { 4, 1, 4 };

// psi = 0.25*R^2 on the straight cylinder, so R = 2*sqrt(psi).
static double
cylinder_R_of_psi(double psi)
{
  return 2.0*sqrt(psi);
}

static struct gk_geometry *
make_cylinder_geometry(int poly_order, struct gkyl_position_map **pmap_out)
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 3, poly_order);

  double lower[3] = { psi_min, alpha_min, -M_PI+1e-14 };
  double upper[3] = { psi_max, alpha_max,  M_PI-1e-14 };

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, (int[3]){ cells[0], cells[1], cells[2] });

  struct gkyl_range ext_range, range;
  int nghost[3] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  struct gkyl_efit_inp einp = {
    .filepath = "gyrokinetic/data/eqdsk/straight_cylinder.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };
  struct gkyl_tok_geo_grid_inp ginp = {
    // straight_cylinder.geqdsk is an analytic equilibrium with no vessel
    // outline in the file, so the absence is declared here.
    .no_vessel_outline = true,
    .rclose = 0.5,
    .zmin = -1., .zmax = 1.,
    .rleft = 0.001, .rmax = 1.0, .rright = 1.0,
  };
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_TOKAMAK,
    .efit_info = einp,
    .tok_grid_info = ginp,
    .position_map = pmap,
    .grid = grid, .local = range, .local_ext = ext_range,
    .global = range, .global_ext = ext_range, .basis = basis,
    .geo_grid = grid, .geo_local = range, .geo_local_ext = ext_range,
    .geo_global = range, .geo_global_ext = ext_range, .geo_basis = basis,
  };

  struct gk_geometry *geom = gkyl_gk_geometry_tok_new(&geometry_input);
  *pmap_out = pmap;
  return geom;
}

// Can a geometry of this order actually be completed by this build?
//
// Construction forms bmag_inv = 1/|B| through gkyl_dg_inv_op_range, whose
// kernel table is sparse: measured 2026-09-09, there is no ser_3x_p2_inv
// anywhere in the tree, so p2 tokamak geometry cannot complete regardless of
// the gridding. The tests below ask the library the same question the library
// asks itself, so they start enforcing automatically the day the kernel lands
// rather than needing to be remembered and re-enabled.
static bool
geometry_order_supported(int poly_order)
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 3, poly_order);
  return gkyl_dg_inv_op_supported(&basis);
}

// Announce a skip loudly. A silently-green test is worse than a missing one:
// the dx_fact fix is NOT covered while this prints.
static bool
skip_unless_supported(int poly_order, const char *what)
{
  if (geometry_order_supported(poly_order)) return false;
  fprintf(stderr,
    "SKIPPED %s: p%d geometry cannot be built by this library "
    "(no DG inverse kernel for bmag_inv). This test does NOT currently cover "
    "the dx_fact fix; it activates when ser_3x_p%d_inv exists.\n",
    what, poly_order, poly_order);
  return true;
}

// The shared body of the corner-node check. At order p the corner nodal range
// is p*cells+1 in every direction and the nodes are uniformly spaced, so node
// ip sits at psi_min + ip*(dx_psi/p). Both the stored computational coordinate
// and the physical R it produced are checked, so a collapse cannot hide behind
// a coincidentally-correct mapping.
static void
check_corner_nodes(int poly_order)
{
  struct gkyl_position_map *pmap = NULL;
  struct gk_geometry *geom = make_cylinder_geometry(poly_order, &pmap);
  TEST_ASSERT(geom != NULL);

  const double dpsi_expected = (psi_max-psi_min)/cells[PSI_IDX]/poly_order;
  const double dalpha_expected = (alpha_max-alpha_min)/cells[AL_IDX]/poly_order;

  // The nodal range itself must be order-correct: p*cells+1 per direction.
  for (int d=0; d<3; ++d) {
    int expected = poly_order*cells[d] + 1;
    TEST_CHECK(gkyl_range_shape(&geom->nrange_corn, d) == expected);
    TEST_MSG("dir %d: corner nodal shape %d, expected %d (p%d)",
      d, gkyl_range_shape(&geom->nrange_corn, d), expected, poly_order);
  }

  // A tolerance on psi of 1e-12 is ~1e-11 relative here; the construction is
  // arithmetic on the bounds, not an iterative solve, so this is generous.
  const double tol_psi = 1.0e-12;
  const double tol_R = 1.0e-9;   // R comes from a root find in the EFIT map.

  int idx[3];
  double psi_seen_lo = 0.0, psi_seen_hi = 0.0;
  int nfail_psi = 0, nfail_R = 0, nfail_alpha = 0;

  for (int ip=geom->nrange_corn.lower[PSI_IDX]; ip<=geom->nrange_corn.upper[PSI_IDX]; ++ip) {
    for (int ia=geom->nrange_corn.lower[AL_IDX]; ia<=geom->nrange_corn.upper[AL_IDX]; ++ia) {
      for (int it=geom->nrange_corn.lower[TH_IDX]; it<=geom->nrange_corn.upper[TH_IDX]; ++it) {
        idx[PSI_IDX] = ip; idx[AL_IDX] = ia; idx[TH_IDX] = it;
        long lin = gkyl_range_idx(&geom->nrange_corn, idx);
        const double *nu = gkyl_array_cfetch(geom->geo_corn.mc2nu_pos_nodal, lin);
        const double *xp = gkyl_array_cfetch(geom->geo_corn.mc2p_nodal, lin);

        double psi_expected = psi_min + ip*dpsi_expected;
        double alpha_expected = alpha_min + ia*dalpha_expected;

        if (fabs(nu[PSI_IDX]-psi_expected) > tol_psi) {
          if (nfail_psi++ == 0)
            TEST_MSG("psi node (%d,%d,%d): got %.17g, expected %.17g (dpsi=%.17g)",
              ip, ia, it, nu[PSI_IDX], psi_expected, dpsi_expected);
        }
        if (fabs(nu[AL_IDX]-alpha_expected) > tol_psi) {
          if (nfail_alpha++ == 0)
            TEST_MSG("alpha node (%d,%d,%d): got %.17g, expected %.17g",
              ip, ia, it, nu[AL_IDX], alpha_expected);
        }
        // mc2p_nodal holds [R, Z, phi] at corner nodes.
        double R_expected = cylinder_R_of_psi(psi_expected);
        if (fabs(xp[0]-R_expected) > tol_R) {
          if (nfail_R++ == 0)
            TEST_MSG("R node (%d,%d,%d): got %.17g, expected %.17g",
              ip, ia, it, xp[0], R_expected);
        }
        if (ip == geom->nrange_corn.lower[PSI_IDX]) psi_seen_lo = nu[PSI_IDX];
        if (ip == geom->nrange_corn.upper[PSI_IDX]) psi_seen_hi = nu[PSI_IDX];
      }
    }
  }

  TEST_CHECK(nfail_psi == 0);
  TEST_MSG("p%d: %d corner nodes have the wrong psi", poly_order, nfail_psi);
  TEST_CHECK(nfail_alpha == 0);
  TEST_MSG("p%d: %d corner nodes have the wrong alpha", poly_order, nfail_alpha);
  TEST_CHECK(nfail_R == 0);
  TEST_MSG("p%d: %d corner nodes have the wrong R", poly_order, nfail_R);

  // The blunt non-degeneracy statement, kept separate so a collapsed grid is
  // unmistakable in the log rather than being one more tolerance failure.
  TEST_CHECK(fabs((psi_seen_hi-psi_seen_lo) - (psi_max-psi_min)) < tol_psi);
  TEST_MSG("p%d: psi spans %.17g across the block, expected %.17g%s",
    poly_order, psi_seen_hi-psi_seen_lo, psi_max-psi_min,
    fabs(psi_seen_hi-psi_seen_lo) < tol_psi
      ? "  -- the nodal grid has COLLAPSED to a single psi" : "");

  gkyl_gk_geometry_release(geom);
  gkyl_position_map_release(pmap);
}

// Positive control: p1 is the only order the campaign has ever run, and the
// defect is invisible there. If this fails, the fixture is wrong, not the code.
static void
test_p1_corner_node_spacing(void)
{
  check_corner_nodes(1);
}

// Defect 3.1. Fails on the unfixed library: dpsi and dalpha are multiplied by
// zero, so every corner node reports psi = psi_min and alpha = alpha_min.
static void
test_p2_corner_node_spacing(void)
{
  if (skip_unless_supported(2, "p2_corner_node_spacing")) return;
  check_corner_nodes(2);
}

// Defect 3.1, second consequence. dzc[] carries the finite-difference steps
// used by the metric calculation; they are dx*1e-2 scaled by the same factor.
// At p2 the unfixed library sets all three to exactly zero, which the metric
// code then divides by.
static void
test_p2_fd_steps_nonzero(void)
{
  if (skip_unless_supported(2, "p2_fd_steps_nonzero")) return;
  struct gkyl_position_map *pmap = NULL;
  struct gk_geometry *geom = make_cylinder_geometry(2, &pmap);
  TEST_ASSERT(geom != NULL);

  const char *name[3] = { "dzc[0] (psi)", "dzc[1] (alpha)", "dzc[2] (theta)" };
  for (int d=0; d<3; ++d) {
    TEST_CHECK(geom->dzc[d] > 0.0 && isfinite(geom->dzc[d]));
    TEST_MSG("%s = %.17g; the metric finite differences divide by this",
      name[d], geom->dzc[d]);
  }

  // The steps are 1e-2 of the nodal spacing, so state the expected values too
  // rather than only that they are nonzero.
  double dpsi_node = (psi_max-psi_min)/cells[PSI_IDX]/2.0;
  TEST_CHECK(fabs(geom->dzc[0] - dpsi_node*1e-2) < 1.0e-15);
  TEST_MSG("dzc[0] = %.17g, expected %.17g", geom->dzc[0], dpsi_node*1e-2);

  gkyl_gk_geometry_release(geom);
  gkyl_position_map_release(pmap);
}

// Defect 3.1b. Interior nodes are the Gauss-Legendre points of each cell, and
// a basis of order p carries p+1 of them per direction -- three at p2, not two.
// They are placed by calc_running_coord, which walks a hard-coded two-step
// cycle {1/sqrt3, 1-1/sqrt3}. That cycle IS the 2-point Gauss rule, so it is
// exactly right at p1 and cannot express three points at any stride.
//
// This is independent of the dx_fact typo and is NOT fixed by correcting it:
// with dx_fact = 1/2 the walk puts the three nodes at 0.2113h, 0.5h, 0.7113h
// instead of 0.1127h, 0.5h, 0.8873h, and then steps onto the cell boundary.
// Only the centre node is right.
static void
test_p2_interior_quadrature_nodes(void)
{
  if (skip_unless_supported(2, "p2_interior_quadrature_nodes")) return;
  struct gkyl_position_map *pmap = NULL;
  struct gk_geometry *geom = make_cylinder_geometry(2, &pmap);
  TEST_ASSERT(geom != NULL);

  const int nq = 3;  // poly_order+1 quadrature points per direction at p2

  // Gauss-Legendre 3-point ordinates on [-1,1], written out independently
  // rather than read from the library's own table, so that a wrong table
  // cannot make this test agree with the code it is checking.
  const double X[3] = { -0.7745966692414833770359, 0.0, 0.7745966692414833770359 };

  const double h = (psi_max-psi_min)/cells[PSI_IDX];   // cell width in psi
  // The defect displaces a node by ~0.0986*h, so a tolerance three orders of
  // magnitude below that separates it from the EFIT interpolation error
  // without being a tuned constant.
  const double tol = 1.0e-3*h;

  TEST_CHECK(gkyl_range_shape(&geom->nrange_int, PSI_IDX) == cells[PSI_IDX]*nq);
  TEST_MSG("interior nodal shape in psi is %d, expected %d",
    gkyl_range_shape(&geom->nrange_int, PSI_IDX), cells[PSI_IDX]*nq);

  int idx[3];
  int nfail = 0;
  for (int ip=geom->nrange_int.lower[PSI_IDX]; ip<=geom->nrange_int.upper[PSI_IDX]; ++ip) {
    int i = ip - geom->nrange_int.lower[PSI_IDX];
    int cell = i/nq, k = i%nq;
    double psi_expected = psi_min + cell*h + h*(1.0 + X[k])/2.0;

    idx[PSI_IDX] = ip;
    idx[AL_IDX] = geom->nrange_int.lower[AL_IDX];
    idx[TH_IDX] = geom->nrange_int.lower[TH_IDX];
    const double *xp = gkyl_array_cfetch(geom->geo_int.mc2p_nodal,
      gkyl_range_idx(&geom->nrange_int, idx));

    // psi = 0.25 R^2 on this equilibrium, so the node's psi is recoverable
    // from the physical R it produced.
    double psi_actual = 0.25*xp[0]*xp[0];
    if (fabs(psi_actual-psi_expected) > tol) {
      if (nfail++ < 4)
        TEST_MSG("interior node %d (cell %d, gauss point %d): psi %.17g, "
          "expected %.17g, off by %.3g (cell width %.17g)",
          i, cell, k, psi_actual, psi_expected,
          psi_actual-psi_expected, h);
    }
  }
  TEST_CHECK(nfail == 0);
  TEST_MSG("%d of %d interior psi nodes are not at the cell's Gauss-3 points",
    nfail, cells[PSI_IDX]*nq);

  gkyl_gk_geometry_release(geom);
  gkyl_position_map_release(pmap);
}

// The p2 blocker itself, asserted as a contract: a geometry the library cannot
// finish must be refused up front, by name, leaving nothing behind -- not
// discovered as a NULL function pointer thousands of nodes into construction.
static void
test_p2_refused_when_no_inverse_kernel(void)
{
  if (geometry_order_supported(2)) {
    // The kernel now exists: the refusal must be gone, not merely unreachable.
    TEST_CHECK(true);
    TEST_MSG("p2 is now supported; the refusal below no longer applies and the "
      "p2 node-placement tests should be enforcing.");
    return;
  }

  FILE *capture = tmpfile();
  TEST_ASSERT(capture != NULL);
  int errfd = fileno(capture);
  fflush(NULL);
  pid_t pid = fork();
  TEST_ASSERT(pid >= 0);
  if (pid == 0) {
    dup2(errfd, STDERR_FILENO);
    struct gkyl_position_map *pmap = NULL;
    struct gk_geometry *geom = make_cylinder_geometry(2, &pmap);
    fprintf(stderr, "RETURNED_WITHOUT_REFUSING geom=%p\n", (void *) geom);
    fflush(stderr);
    _exit(0);
  }
  int status = 0;
  TEST_ASSERT(waitpid(pid, &status, 0) == pid);

  char report[8192] = { 0 };
  rewind(capture);
  size_t n = fread(report, 1, sizeof(report)-1, capture);
  report[n] = '\0';
  fclose(capture);

  bool exited_cleanly = WIFEXITED(status) && WEXITSTATUS(status) == 0;
  TEST_CHECK(!exited_cleanly);
  TEST_MSG("p2 was accepted though it cannot be completed; child said: %s",
    report[0] ? report : "(nothing)");

  TEST_CHECK(strstr(report, "GKYL_GEOMETRY_UNSUPPORTED_POLY_ORDER") != NULL);
  TEST_MSG("no named refusal. stderr: %s", report[0] ? report : "(nothing)");
  TEST_CHECK(strstr(report, "reason=no_dg_inverse_kernel") != NULL);
  TEST_MSG("the refusal must name the CAUSE, not just the order. stderr: %s",
    report[0] ? report : "(nothing)");
  // It must refuse before the folded-cell guard, i.e. before building a grid.
  TEST_CHECK(strstr(report, "TOK_GEO_FOLDED_CELLS") == NULL);
  TEST_MSG("geometry was constructed before the refusal; it should refuse "
    "while nothing is allocated. stderr: %s", report);
}

// Defect 3.3. gkyl_gk_geometry_init_nodal_range sizes an uninitialised
// int nodes[GKYL_MAX_DIM] for poly_order 1 and 2 and has no else branch, so at
// p3 stack garbage reaches gkyl_range_init_from_shape. An order outside the
// implemented contract must be refused explicitly and say so.
//
// Run in a child process because the refusal is expected to abort. The check
// is on the DIAGNOSTIC, not merely on dying: an uninitialised shape can crash
// on its own, and a crash without a message would let the unfixed library pass.
static void
test_unsupported_poly_order_rejected(void)
{
  const int unsupported[] = { 3, 4 };
  for (size_t k=0; k<sizeof(unsupported)/sizeof(unsupported[0]); ++k) {
    int order = unsupported[k];

    FILE *capture = tmpfile();
    TEST_ASSERT(capture != NULL);
    int errfd = fileno(capture);
    TEST_ASSERT(errfd >= 0);

    fflush(NULL);
    pid_t pid = fork();
    TEST_ASSERT(pid >= 0);
    if (pid == 0) {
      dup2(errfd, STDERR_FILENO);
      double lower[3] = { 0.0, 0.0, 0.0 }, upper[3] = { 1.0, 1.0, 1.0 };
      struct gkyl_rect_grid grid;
      gkyl_rect_grid_init(&grid, 3, lower, upper, (int[3]){ 2, 2, 2 });
      struct gkyl_range ext_range, range;
      gkyl_create_grid_ranges(&grid, (int[3]){ 1, 1, 1 }, &ext_range, &range);
      struct gkyl_range nrange;
      gkyl_gk_geometry_init_nodal_range(&nrange, &range, order);
      // Reached only if the unsupported order was silently accepted.
      fprintf(stderr, "RETURNED_WITHOUT_REJECTING volume=%ld\n",
        (long) nrange.volume);
      fflush(stderr);
      _exit(0);
    }

    int status = 0;
    TEST_ASSERT(waitpid(pid, &status, 0) == pid);

    char report[4096] = { 0 };
    rewind(capture);
    size_t n = fread(report, 1, sizeof(report)-1, capture);
    report[n] = '\0';
    fclose(capture);

    bool exited_cleanly = WIFEXITED(status) && WEXITSTATUS(status) == 0;
    bool named_the_order = strstr(report, "GKYL_GEOMETRY_UNSUPPORTED_POLY_ORDER") != NULL;

    TEST_CHECK(!exited_cleanly);
    TEST_MSG("poly_order %d was accepted rather than refused; child said: %s",
      order, report[0] ? report : "(nothing)");
    TEST_CHECK(named_the_order);
    TEST_MSG("poly_order %d: no GKYL_GEOMETRY_UNSUPPORTED_POLY_ORDER diagnostic. "
      "Child %s. stderr: %s", order,
      WIFEXITED(status) ? "exited" : "was killed by a signal",
      report[0] ? report : "(nothing)");
  }
}

TEST_LIST = {
  { "p1_corner_node_spacing", test_p1_corner_node_spacing },
  { "p2_corner_node_spacing", test_p2_corner_node_spacing },
  { "p2_fd_steps_nonzero", test_p2_fd_steps_nonzero },
  { "p2_interior_quadrature_nodes", test_p2_interior_quadrature_nodes },
  { "p2_refused_when_no_inverse_kernel", test_p2_refused_when_no_inverse_kernel },
  { "unsupported_poly_order_rejected", test_unsupported_poly_order_rejected },
  { NULL, NULL },
};
