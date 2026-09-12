/* Positive-coverage fixture for GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN.
 *
 * DN_SOL_IN is declared in the enum and implemented in tok_geo.c -- handled
 * symmetrically with DN_SOL_OUT at the arc-length, alpha and midplane sites --
 * but as of 2026-09-11 no fixture anywhere in the tree exercises it, so it is
 * the one declared type with no positive evidence (USN_UP aside, which is an
 * accepted exclusion).
 *
 * Every geometry parameter below is COPIED from the inboard blocks of
 * rt_gk_multib_step_2x2v_p1.c, not invented. That fixture's DN_SOL_IN_LO,
 * DN_SOL_IN_MID and DN_SOL_IN_UP blocks share one identical tok_grid_info and
 * differ only in ftype and block bounds, and their union is the full inboard
 * SOL that DN_SOL_IN describes. The psi interval is that fixture's
 * [psisep - win, psisep] with win = 0.05 and npsi_inner_sol = 2.
 *
 * This is geometry-only: no species, no field, no stepping. It exists to prove
 * the type constructs, not to be a physics regression test.
 *
 * step.geqdsk has limitr=1 (limiter_status=2, a degenerate outline), so the
 * vessel-outline acknowledgement is required -- see commits 158df0bd5 and
 * f48fcc610.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_efit.h>
#include <gkyl_gk_block_geom.h>
#include <gkyl_gyrokinetic_comms.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_tok_geo.h>
#include <rt_arg_parse.h>

/* Inboard plate shapes, copied verbatim from rt_gk_multib_step_2x2v_p1.c. */
static void
shaped_pfunc_upper_inner(double s, double *RZ)
{
  RZ[0] = 1.651 + (1.8 - 1.651)*s;
  RZ[1] = 6.331 + (6.777 - 6.331)*s;
}

static void
shaped_pfunc_lower_inner(double s, double *RZ)
{
  RZ[0] = 1.651 + (1.8 - 1.651)*s;
  RZ[1] = -(6.331 + (6.777 - 6.331)*s);
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  struct gkyl_efit_inp efit_inp = {
    .filepath = "gyrokinetic/data/eqdsk/step.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };

  /* The inboard SOL psi interval, as declared in the multiblock STEP fixture. */
  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  double psisep = efit->psisep;
  gkyl_efit_release(efit);

  const double win = 0.05;                 /* inner-SOL width  */
  const double zinner = 6.34;              /* inboard z extent */
  const double rright_out = 5.2;
  double psi_lo_inner_sol = psisep - win;
  double psi_up_inner_sol = psisep;

  struct gkyl_gk_block_geom *bgeom = gkyl_gk_block_geom_new(2, 1);

  gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
    .lower = { psi_lo_inner_sol, -M_PI },
    .upper = { psi_up_inner_sol,  M_PI },
    .cells = { 2, 16 },                    /* npsi_inner_sol = 2 */
    .cuts  = { 1, 1 },
    .geometry = {
      .world = { 0.0 },
      .geometry_id = GKYL_GEOMETRY_TOKAMAK,
      .efit_info = efit_inp,
      .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
        .ftype = GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN,
        .no_vessel_outline = true,
        .rleft = 2.0,
        .rright = rright_out,
        .rmin = 0.0,
        .rmax = 6.2,
        .zmin = -zinner,
        .zmax = zinner,
        .plate_spec = true,
        .plate_func_upper = shaped_pfunc_upper_inner,
        .plate_func_lower = shaped_pfunc_lower_inner,
      },
    },
    /* One isolated block: every edge is a physical boundary, so no reciprocal
     * connection has to be invented for a type whose neighbours are absent. */
    .connections[0] = {
      { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL },
      { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL },
    },
    .connections[1] = {
      { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL },
      { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL },
    },
  });

  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi,
    app_args.use_gpu, stderr);

  struct gkyl_gyrokinetic_multib app_inp = {
    .name = "rt_gk_multib_step_dn_sol_in_geom",
    .cdim = 2,
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .use_gpu = app_args.use_gpu,
    .gk_block_geom = bgeom,
    .num_species = 0,
    .num_neut_species = 0,
    .skip_field = true,
    .comm = comm,
  };

  struct gkyl_gyrokinetic_multib_app *app =
    gkyl_gyrokinetic_multib_app_new_geom(&app_inp);
  if (app == 0) {
    fprintf(stderr, "GEOM_FAIL: could not construct DN_SOL_IN geometry\n");
    gkyl_gk_block_geom_release(bgeom);
    gkyl_gyrokinetic_comms_release(comm);
    return 1;
  }

  /* A geometry-only app is torn down with _release_geom; the general
   _release path assumes solver state this app never built. */
  gkyl_gyrokinetic_multib_app_release_geom(app);
  gkyl_gk_block_geom_release(bgeom);
  gkyl_gyrokinetic_comms_release(comm);
  return 0;
}
