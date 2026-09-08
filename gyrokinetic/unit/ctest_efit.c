#include <stdio.h>
#include <math.h>
#include <string.h>
#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>

#include <gkyl_efit.h>

void test_efit_solovev_ho()
{
  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/solovev.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "solovev_psi.gkyl");
  gkyl_grid_sub_array_write(
    &efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "solovev_fpol.gkyl"
  );
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "solovev_q.gkyl");

  gkyl_efit_release(efit);
}

void test_efit_step_ho()
{
  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/step.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "step_psi.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "step_fpol.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "step_q.gkyl");

  gkyl_efit_release(efit);
}

void test_efit_nstxu_ho()
{
  // Uses DN configuration by default, but one can switch to SN by changing the filepath if desired.
  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/nstxu_DN.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  // printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g psisep=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry, efit->psisep);
  gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "nstxu_DN_psi.gkyl");
  gkyl_grid_sub_array_write(
    &efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "nstxu_DN_fpol.gkyl"
  );
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "nstxu_DN_q.gkyl");

  gkyl_efit_release(efit);
}

void test_efit_asdex_ho()
{
  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk", .rz_poly_order = 2, .flux_poly_order = 1
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "asdex_psi.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "asdex_fpol.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "asdex_q.gkyl");

  gkyl_efit_release(efit);
}

void test_efit_cerfon_ho()
{
  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/cerfon.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "cerfon_psi.gkyl");
  gkyl_grid_sub_array_write(
    &efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "cerfon_fpol.gkyl"
  );
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "cerfon_q.gkyl");

  gkyl_efit_release(efit);
}

void test_efit_elliptical_ho()
{
  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/elliptical.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "elliptical_psi.gkyl");
  gkyl_grid_sub_array_write(
    &efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "elliptical_fpol.gkyl"
  );
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "elliptical_q.gkyl");

  gkyl_efit_release(efit);
}

void test_efit_wham_ho()
{
  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/wham.geqdsk", .rz_poly_order = 2, .flux_poly_order = 1
    //.reflect =  true,
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "wham_psi.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "wham_fpol.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "wham_q.gkyl");

  gkyl_efit_release(efit);
}

void test_efit_tcv_ho()
{
  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/tcv.geqdsk", .rz_poly_order = 2, .flux_poly_order = 1
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "tcv_psi.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "tcv_fpol.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "tcv_q.gkyl");

  gkyl_efit_release(efit);
}

void test_efit_mast_ho()
{
  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/mast.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "mast_psi.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "mast_fpol.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "mast_q.gkyl");

  gkyl_efit_release(efit);
}

void test_efit_ltx_ho()
{
  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/LTX_103955_03.eqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "ltx_psi.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "ltx_fpol.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "ltx_q.gkyl");

  gkyl_efit_release(efit);
}

TEST_LIST = {
  {"test_efit_solovev_ho", test_efit_solovev_ho},
  {"test_efit_step_ho", test_efit_step_ho},
  {"test_efit_asdex_ho", test_efit_asdex_ho},
  {"test_efit_nstxu_ho", test_efit_nstxu_ho},
  {"test_efit_cerfon_ho", test_efit_cerfon_ho},
  {"test_efit_elliptical_ho", test_efit_elliptical_ho},
  {"test_efit_wham_ho", test_efit_wham_ho},
  {"test_efit_tcv_ho", test_efit_tcv_ho},
  {"test_efit_mast_ho", test_efit_mast_ho},
  {"test_efit_ltx_ho", test_efit_ltx_ho},
  {NULL, NULL}
};
