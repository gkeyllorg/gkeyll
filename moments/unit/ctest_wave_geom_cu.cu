/* -*- c -*- */

#include <stdio.h>

extern "C" {
#include <float.h>
#include <gkyl_array.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_geom_priv.h>
int cu_wave_geom_test(const struct gkyl_wave_geom *wg);
}

// FIXME: duplicate of gkyl_compare_double in util.c
GKYL_CU_D static int compare(double a, double b, double eps)
{
  double absa = fabs(a), absb = fabs(b), diff = fabs(a-b);
  if (a == b) return 1;
  if (a == 0 || b == 0 || (absa+absb < DBL_MIN)) return diff < eps;
  if (absa < eps) return diff < eps;
  if (absb < eps) return diff < eps;
  return diff/fmin(absa+absb, DBL_MAX) < eps;
}

__global__
void ker_cu_wave_geom_test(const struct gkyl_wave_geom *wg, int *nfail)
{
  *nfail = 0;

  double r_inn = 0.25, r_out = 1.25;
  double phi_max = M_PI / 2.;
  double area = 0.5*(r_out*r_out-r_inn*r_inn);
  double edge_inn = sqrtf(2) * r_inn;
  double area_c = (r_out - r_inn) * phi_max;

  int idx[] = {1, 1, 1};
  {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, idx);

    GKYL_CU_CHECK(compare(cg->kappa, area / area_c, 1e-8), nfail);
    GKYL_CU_CHECK(compare(cg->lenr[0], edge_inn / phi_max, 1e-8), nfail);
    GKYL_CU_CHECK(compare(cg->lenr[1], 1, 1e-8), nfail);
    GKYL_CU_CHECK(compare(cg->lenr[2], area / area_c, 1e-8), nfail);

    // normal to left face has phi angle 45 deg
    GKYL_CU_CHECK(compare(cg->norm[0][0], 1/sqrtf(2.), 1e-8), nfail);
    GKYL_CU_CHECK(compare(cg->norm[0][1], 1/sqrtf(2.), 1e-8), nfail);
    GKYL_CU_CHECK(compare(cg->norm[0][2], 0.0, 1e-8), nfail);

    // tangent1 to left face has phi angle 135 deg
    GKYL_CU_CHECK(compare(cg->tau1[0][0], -1/sqrtf(2.), 1e-8), nfail);
    GKYL_CU_CHECK(compare(cg->tau1[0][1], 1/sqrtf(2.), 1e-8), nfail);
    GKYL_CU_CHECK(compare(cg->tau1[0][2], 0.0, 1e-8), nfail);

    // tangent2 to left face is ez
    GKYL_CU_CHECK(compare(cg->tau2[0][0], 0.0, 1e-8), nfail);
    GKYL_CU_CHECK(compare(cg->tau2[0][1], 0.0, 1e-8), nfail);
    GKYL_CU_CHECK(compare(cg->tau2[0][2], 1.0, 1e-8), nfail);
  }
}

int cu_wave_geom_test(const struct gkyl_wave_geom *wg_dev)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));

  ker_cu_wave_geom_test<<<1, 1>>>(wg_dev, nfail_dev);
  checkCuda(cudaGetLastError());

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}

