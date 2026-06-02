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

void test_solovev(){

  struct gkyl_efit_inp inp  = {
    .filepath = "gyrokinetic/data/eqdsk/solovev.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect =  true,
  };
  struct gkyl_efit* efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);

  gkyl_efit_release(efit);

}

void test_step(){
  struct gkyl_efit_inp inp  = {
    .filepath = "gyrokinetic/data/eqdsk/step.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect =  true,
  };
  struct gkyl_efit* efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);

  gkyl_efit_release(efit);

}

void test_nstxu(){
  // Uses DN configuration by default, but one can switch to SN by changing the filepath if desired.
  struct gkyl_efit_inp inp  = {
    .filepath = "gyrokinetic/data/eqdsk/nstxu_DN.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect =  true,
  };
  struct gkyl_efit* efit = gkyl_efit_new(&inp);

  // printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g psisep=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry, efit->psisep);

  gkyl_efit_release(efit);

}

void test_asdex(){
  struct gkyl_efit_inp inp  = {
    .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };
  struct gkyl_efit* efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);

  gkyl_efit_release(efit);

}

void test_cerfon(){
  
  struct gkyl_efit_inp inp  = {
    .filepath = "gyrokinetic/data/eqdsk/cerfon.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect =  true,
  };
  struct gkyl_efit* efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);

  gkyl_efit_release(efit);

}

void test_elliptical(){
  struct gkyl_efit_inp inp  = {
    .filepath = "gyrokinetic/data/eqdsk/elliptical.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect =  true,
  };
  struct gkyl_efit* efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);

  gkyl_efit_release(efit);

}

void test_wham(){
  struct gkyl_efit_inp inp  = {
    .filepath = "gyrokinetic/data/eqdsk/wham.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    //.reflect =  true,
  };
  struct gkyl_efit* efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);

  gkyl_efit_release(efit);

}


void test_tcv(){
  struct gkyl_efit_inp inp  = {
    .filepath = "gyrokinetic/data/eqdsk/tcv.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };
  struct gkyl_efit* efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);

  gkyl_efit_release(efit);

}

void test_mast(){
  struct gkyl_efit_inp inp  = {
    .filepath = "gyrokinetic/data/eqdsk/mast.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect =  true,
  };
  struct gkyl_efit* efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);

  gkyl_efit_release(efit);

}

void test_ltx(){
  struct gkyl_efit_inp inp  = {
    .filepath = "gyrokinetic/data/eqdsk/LTX_103955_03.eqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect =  true,
  };
  struct gkyl_efit* efit = gkyl_efit_new(&inp);

  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);

  gkyl_efit_release(efit);

}

TEST_LIST = {
  { "test_solovev", test_solovev},
  { "test_step", test_step},
  { "test_asdex", test_asdex},
  { "test_nstxu", test_nstxu},
  { "test_cerfon", test_cerfon},
  { "test_elliptical", test_elliptical},
  { "test_wham", test_wham},
  { "test_tcv", test_tcv},
  { "test_mast", test_mast},
  { "test_ltx", test_ltx},
  { NULL, NULL },
};
