#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
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

static void
check_same_array(
  const struct gkyl_array *reference, const struct gkyl_array *actual,
  const struct gkyl_range *range
)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);
    const double *expected = gkyl_array_cfetch(reference, loc);
    const double *found = gkyl_array_cfetch(actual, loc);
    for (int k = 0; k < reference->ncomp; ++k) {
      TEST_CHECK(fabs(expected[k] - found[k]) < 1e-8 * (1.0 + fabs(expected[k])));
    }
  }
}

static void
write_cocos_variant(
  const char *source, const char *target, int nr, int nz, double psi_factor, double f_factor,
  double q_factor
)
{
  FILE *input = fopen(source, "r"), *output = fopen(target, "w");
  char header[256];
  TEST_ASSERT(input && output && fgets(header, sizeof(header), input));
  fputs(header, output);
  int psi_start = 20 + 4 * nr, q_start = psi_start + nr * nz;
  int count = q_start + nr;
  for (int i = 0; i < count; ++i) {
    double value;
    TEST_ASSERT(fscanf(input, "%lf", &value) == 1);
    if (i == 7 || i == 8 || i == 11 || i == 17 || (i >= psi_start && i < q_start)) {
      value /= psi_factor;
    } else if (i == 9 || i == 10 || (i >= 20 && i < 20 + nr)) {
      value /= f_factor;
    } else if (i >= 20 + 2 * nr && i < 20 + 3 * nr) {
      value *= psi_factor;
    } else if (i >= q_start) {
      value /= q_factor;
    }
    fprintf(output, "%.17e\n", value);
  }
  fclose(input);
  fclose(output);
}

void
test_efit_cocos_conversion(void)
{
  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };

  gkyl_efit *reference = gkyl_efit_new(&inp);
  TEST_ASSERT(reference && reference->eqdsk_cocos_index == 3);

  const int indices[] = {1, 2, 3, 7, 11, 17};
  const double psi_factors[] = {-1, 1, 1, 1, -1 / (2 * M_PI), 1 / (2 * M_PI)};
  const double f_factors[] = {1, -1, 1, 1, 1, 1};
  const double q_factors[] = {-1, -1, 1, -1, -1, -1};
  for (int i = 0; i < 6; ++i) {
    char path[] = "gkyl_efit_cocos_XXXXXX";
    int fd = mkstemp(path);
    TEST_ASSERT(fd >= 0);
    close(fd);

    write_cocos_variant(
      inp.filepath, path, reference->nr, reference->nz, psi_factors[i], f_factors[i], q_factors[i]
    );

    struct gkyl_efit_inp variant_inp = inp;
    snprintf(variant_inp.filepath, sizeof variant_inp.filepath, "%s", path);
    variant_inp.eqdsk_cocos_index = indices[i];

    gkyl_efit *variant = gkyl_efit_new(&variant_inp);

    TEST_ASSERT(variant && variant->eqdsk_cocos_index == indices[i]);
    TEST_CHECK(fabs(variant->simag - reference->simag) < 1e-10);
    TEST_CHECK(fabs(variant->sibry - reference->sibry) < 1e-10);
    TEST_CHECK(fabs(variant->bcentr - reference->bcentr) < 1e-10);
    TEST_CHECK(fabs(variant->current - reference->current) < 1e-8);
    check_same_array(reference->psizr, variant->psizr, &reference->rzlocal);
    check_same_array(reference->bmagzr, variant->bmagzr, &reference->rzlocal);
    check_same_array(reference->fpolflux, variant->fpolflux, &reference->fluxlocal);
    check_same_array(reference->fpolprimeflux, variant->fpolprimeflux, &reference->fluxlocal);
    check_same_array(reference->qflux, variant->qflux, &reference->fluxlocal);

    gkyl_efit_release(variant);
    unlink(path);
  }

  inp.eqdsk_cocos_index = 4;
  TEST_CHECK(gkyl_efit_new(&inp) == NULL);

  gkyl_efit_release(reference);
}

void
test_efit_solovev_ho()
{
  bool write_files = false;

  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/solovev.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  if (write_files) {
    printf(
      "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n",
      efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis,
      efit->simag, efit->sibry, efit->bcentr, efit->current, efit->simag, efit->rmaxis,
      efit->zmaxis, efit->sibry
    );
    gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "solovev_psi.gkyl");
    gkyl_grid_sub_array_write(
      &efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "solovev_fpol.gkyl"
    );
    gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "solovev_q.gkyl");
  }

  gkyl_efit_release(efit);
}

void
test_efit_step_ho()
{
  bool write_files = false;

  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/step.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  if (write_files) {
    printf(
      "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n",
      efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis,
      efit->simag, efit->sibry, efit->bcentr, efit->current, efit->simag, efit->rmaxis,
      efit->zmaxis, efit->sibry
    );
    gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "step_psi.gkyl");
    gkyl_grid_sub_array_write(
      &efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "step_fpol.gkyl"
    );
    gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "step_q.gkyl");
  }

  gkyl_efit_release(efit);
}

void
test_efit_nstxu_ho()
{
  bool write_files = false;

  // Uses DN configuration by default, but one can switch to SN by changing the filepath if desired.
  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/nstxu_DN.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  if (write_files) {
    printf(
      "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g psisep=%g \n",
      efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis,
      efit->simag, efit->sibry, efit->bcentr, efit->current, efit->simag, efit->rmaxis,
      efit->zmaxis, efit->sibry, efit->psisep
    );
    gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "nstxu_DN_psi.gkyl");
    gkyl_grid_sub_array_write(
      &efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "nstxu_DN_fpol.gkyl"
    );
    gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "nstxu_DN_q.gkyl");
  }

  gkyl_efit_release(efit);
}

void
test_efit_asdex_ho()
{
  bool write_files = false;

  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  if (write_files) {
    printf(
      "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n",
      efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis,
      efit->simag, efit->sibry, efit->bcentr, efit->current, efit->simag, efit->rmaxis,
      efit->zmaxis, efit->sibry
    );
    gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "asdex_psi.gkyl");
    gkyl_grid_sub_array_write(
      &efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "asdex_fpol.gkyl"
    );
    gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "asdex_q.gkyl");
  }

  gkyl_efit_release(efit);
}

void
test_efit_cerfon_ho()
{
  bool write_files = false;

  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/cerfon.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  if (write_files) {
    printf(
      "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n",
      efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis,
      efit->simag, efit->sibry, efit->bcentr, efit->current, efit->simag, efit->rmaxis,
      efit->zmaxis, efit->sibry
    );
    gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "cerfon_psi.gkyl");
    gkyl_grid_sub_array_write(
      &efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "cerfon_fpol.gkyl"
    );
    gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "cerfon_q.gkyl");
  }

  gkyl_efit_release(efit);
}

void
test_efit_elliptical_ho()
{
  bool write_files = false;

  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/elliptical.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  if (write_files) {
    printf(
      "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n",
      efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis,
      efit->simag, efit->sibry, efit->bcentr, efit->current, efit->simag, efit->rmaxis,
      efit->zmaxis, efit->sibry
    );
    gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "elliptical_psi.gkyl");
    gkyl_grid_sub_array_write(
      &efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "elliptical_fpol.gkyl"
    );
    gkyl_grid_sub_array_write(
      &efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "elliptical_q.gkyl"
    );
  }

  gkyl_efit_release(efit);
}

void
test_efit_wham_ho()
{
  bool write_files = false;

  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/wham.geqdsk", .rz_poly_order = 2, .flux_poly_order = 1,
    //.reflect =  true,
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  if (write_files) {
    printf(
      "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n",
      efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis,
      efit->simag, efit->sibry, efit->bcentr, efit->current, efit->simag, efit->rmaxis,
      efit->zmaxis, efit->sibry
    );
    gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "wham_psi.gkyl");
    gkyl_grid_sub_array_write(
      &efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "wham_fpol.gkyl"
    );
    gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "wham_q.gkyl");
  }

  gkyl_efit_release(efit);
}

void
test_efit_tcv_ho()
{
  bool write_files = false;

  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/tcv.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  if (write_files) {
    printf(
      "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n",
      efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis,
      efit->simag, efit->sibry, efit->bcentr, efit->current, efit->simag, efit->rmaxis,
      efit->zmaxis, efit->sibry
    );
    gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "tcv_psi.gkyl");
    gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "tcv_fpol.gkyl");
    gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "tcv_q.gkyl");
  }

  gkyl_efit_release(efit);
}

void
test_efit_mast_ho()
{
  bool write_files = false;

  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/mast.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  if (write_files) {
    printf(
      "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n",
      efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis,
      efit->simag, efit->sibry, efit->bcentr, efit->current, efit->simag, efit->rmaxis,
      efit->zmaxis, efit->sibry
    );
    gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "mast_psi.gkyl");
    gkyl_grid_sub_array_write(
      &efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "mast_fpol.gkyl"
    );
    gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "mast_q.gkyl");
  }

  gkyl_efit_release(efit);
}

void
test_efit_ltx_ho()
{
  bool write_files = false;

  struct gkyl_efit_inp inp = {
    .filepath = "gyrokinetic/data/eqdsk/LTX_103955_03.eqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };
  struct gkyl_efit *efit = gkyl_efit_new(&inp);

  if (write_files) {
    printf(
      "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n",
      efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis,
      efit->simag, efit->sibry, efit->bcentr, efit->current, efit->simag, efit->rmaxis,
      efit->zmaxis, efit->sibry
    );
    gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "ltx_psi.gkyl");
    gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "ltx_fpol.gkyl");
    gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "ltx_q.gkyl");
  }

  gkyl_efit_release(efit);
}

TEST_LIST = {
  {"test_efit_cocos_conversion", test_efit_cocos_conversion},
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
