/* -*- c++ -*- */

// Device-side functions (and getters for their device addresses) used by the
// GPU tests in ctest_proj_on_basis.c. Each function must match its host
// counterpart in that file exactly.

extern "C" {
#include <gkyl_proj_on_basis.h>
#include <gkyl_util.h>

evalf_t ctest_proj_on_basis_f_1d_cu_dev_ptr(void);
evalf_t ctest_proj_on_basis_f_2d_2c_cu_dev_ptr(void);
proj_on_basis_c2p_t ctest_proj_on_basis_c2p_1d_cu_dev_ptr(void);
}

// Context for the 2d two-component function; must match the definition in
// ctest_proj_on_basis.c.
struct ctest_proj_on_basis_2d_ctx {
  double c0, c1;
};

GKYL_CU_DH static void
ctest_pob_f_1d(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  fout[0] = x*x;
}

GKYL_CU_DH static void
ctest_pob_f_2d_2c(double t, const double *xn, double *fout, void *ctx)
{
  struct ctest_proj_on_basis_2d_ctx *tctx = (struct ctest_proj_on_basis_2d_ctx *) ctx;
  double x = xn[0], y = xn[1];
  fout[0] = tctx->c0 + x*y;
  fout[1] = tctx->c1*x*x + y;
}

GKYL_CU_DH static void
ctest_pob_c2p_1d(const double *xcomp, double *xphys, void *ctx)
{
  xphys[0] = 0.5*xcomp[0] + 0.1;
}

GKYL_DEFINE_CU_DEV_FUNC_GETTER(ctest_pob_f_1d, evalf_t, ctest_pob_f_1d_getter);
GKYL_DEFINE_CU_DEV_FUNC_GETTER(ctest_pob_f_2d_2c, evalf_t, ctest_pob_f_2d_2c_getter);
GKYL_DEFINE_CU_DEV_FUNC_GETTER(ctest_pob_c2p_1d, proj_on_basis_c2p_t, ctest_pob_c2p_1d_getter);

extern "C" evalf_t
ctest_proj_on_basis_f_1d_cu_dev_ptr(void)
{
  return ctest_pob_f_1d_getter();
}

extern "C" evalf_t
ctest_proj_on_basis_f_2d_2c_cu_dev_ptr(void)
{
  return ctest_pob_f_2d_2c_getter();
}

extern "C" proj_on_basis_c2p_t
ctest_proj_on_basis_c2p_1d_cu_dev_ptr(void)
{
  return ctest_pob_c2p_1d_getter();
}
