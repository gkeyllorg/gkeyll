// Regression test for the kernel dispatch of the tensor p=1 hybrid phase basis.
//
// Kernel tables are keyed on the configuration-basis type, and only the tensor
// tables hold the hybrid phase-space kernels (p=1 in configuration space, p=2
// in velocity space; 216 modes at 3x3v vs. 64 for the pure p=1 basis). The
// Vlasov updaters therefore dispatch through gkyl_basis_phase_kernel_type(),
// which routes a hybrid phase basis to the tensor tables no matter how the
// (p=1) configuration basis is typed. Before that routing existed, a
// Serendipity-typed configuration basis (e.g. a gyrokinetic app's global basis
// paired with a hybrid neutral distribution) selected the pure p=1 kernels,
// which mis-index the hybrid coefficients: at 3x3v the pure-p1 index of the
// xyz configuration mode lands on a hybrid velocity-space mode and M0 is wrong
// at O(1). 1x3v and 2x3v M0 happen to be unaffected, so the 3x3v case is the
// one that guards the routing.
//
// These tests project a random distribution onto *all* hybrid modes and check
// the moment against direct quadrature for both configuration-basis typings.
#include <acutest.h>
#include <math.h>
#include <stdlib.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_velocity_map.h>

static double
lcg_rand(unsigned *s)
{
  *s = (*s) * 1103515245u + 12345u;
  return ((double)((*s >> 8) & 0xffff)) / 65535.0 - 0.5;
}

// Returns max |M0 - quadrature reference| over the configuration modes when the
// moment updater is constructed with a configuration basis of type ctype.
static double
hybrid_m0_error(int cdim, enum gkyl_basis_type ctype, bool use_gpu)
{
  int vdim = 3, pdim = cdim + vdim;
  double lower[6] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0}, upper[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  int cells[6] = {1, 1, 1, 1, 1, 1}, ghost[6] = {0};
  struct gkyl_rect_grid pgrid, cgrid, vgrid;
  gkyl_rect_grid_init(&pgrid, pdim, lower, upper, cells);
  gkyl_rect_grid_init(&cgrid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&vgrid, vdim, lower + cdim, upper + cdim, cells + cdim);
  struct gkyl_range plocal, plocal_ext, clocal, clocal_ext, vlocal, vlocal_ext;
  gkyl_create_grid_ranges(&pgrid, ghost, &plocal_ext, &plocal);
  gkyl_create_grid_ranges(&cgrid, ghost, &clocal_ext, &clocal);
  gkyl_create_grid_ranges(&vgrid, ghost, &vlocal_ext, &vlocal);

  struct gkyl_basis pbasis, cbasis, cbasis_ser, vbasis;
  gkyl_cart_modal_hybrid(&pbasis, cdim, vdim);
  gkyl_cart_modal_serendip(&cbasis_ser, cdim, 1);
  if (ctype == GKYL_BASIS_MODAL_TENSOR) {
    gkyl_cart_modal_tensor(&cbasis, cdim, 1);
  } else {
    gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  }
  gkyl_cart_modal_tensor(&vbasis, vdim, 2);

  // Identity velocity map and sparse Hamiltonian H = v^2/2 on the p=2 velocity basis.
  struct gkyl_vlasov_velocity_map_inp vmap_inp[GKYL_MAX_CDIM] = {0};
  struct gkyl_vlasov_velocity_map *vmap =
    gkyl_vlasov_velocity_map_new(&vgrid, &vlocal, &vbasis, vmap_inp, false, use_gpu);
  struct gkyl_array *hamil =
    use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, vbasis.num_basis, vlocal_ext.volume) :
              gkyl_array_new(GKYL_DOUBLE, vbasis.num_basis, vlocal_ext.volume);
  struct gkyl_array *hinv =
    use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, vbasis.num_basis, vlocal_ext.volume) :
              gkyl_array_new(GKYL_DOUBLE, vbasis.num_basis, vlocal_ext.volume);
  gkyl_dg_vlasov_calc_hamil(
    &vgrid, &vbasis, &vlocal, GKYL_MODEL_DEFAULT, vmap, hamil, hinv, use_gpu
  );

  // Random coefficients on every hybrid mode.
  struct gkyl_array *f_ho = gkyl_array_new(GKYL_DOUBLE, pbasis.num_basis, plocal_ext.volume);
  unsigned seed = 7;
  double *fd = gkyl_array_fetch(f_ho, 0);
  for (int k = 0; k < pbasis.num_basis; ++k) {
    fd[k] = lcg_rand(&seed);
  }
  struct gkyl_array *f = f_ho;
  if (use_gpu) {
    f = gkyl_array_cu_dev_new(GKYL_DOUBLE, pbasis.num_basis, plocal_ext.volume);
    gkyl_array_copy(f, f_ho);
  }

  // Reference: M0_j = int f phi_j dv by 3-point Gauss quadrature per direction
  // (exact for the p<=2 integrand on the [-1,1]^pdim computational cell).
  double ref[8] = {0.0};
  const double *xq = gkyl_gauss_ordinates_3, *wq = gkyl_gauss_weights_3;
  int nq = 1;
  for (int d = 0; d < pdim; ++d) {
    nq *= 3;
  }
  for (int q = 0; q < nq; ++q) {
    double eta[6], w = 1.0;
    int qq = q;
    for (int d = 0; d < pdim; ++d) {
      int i = qq % 3;
      qq /= 3;
      eta[d] = xq[i];
      w *= wq[i];
    }
    double pb[216];
    pbasis.eval(eta, pb);
    double fv = 0.0;
    for (int k = 0; k < pbasis.num_basis; ++k) {
      fv += fd[k] * pb[k];
    }
    double cb[8];
    cbasis_ser.eval(eta, cb);
    for (int j = 0; j < cbasis_ser.num_basis; ++j) {
      ref[j] += w * fv * cb[j];
    }
  }

  struct gkyl_mom_vlasov_inp inp = {
    .conf_basis = &cbasis,
    .phase_basis = &pbasis,
    .vel_range = &vlocal,
    .vel_map = vmap,
    .hamil_range = &vlocal,
    .hamil = hamil,
    .model_id = GKYL_MODEL_TRIAD,
    .hamil_id = GKYL_HAMIL_VEL_SPARSE,
    .mom_type = GKYL_F_MOMENT_M0,
    .use_gpu = use_gpu,
  };
  struct gkyl_mom_type *mt = gkyl_mom_vlasov_inew(&inp);
  gkyl_mom_calc *mc = gkyl_mom_calc_new(&pgrid, mt, use_gpu);
  struct gkyl_array *m0 =
    use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume) :
              gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);
  gkyl_array_clear(m0, 0.0);
  gkyl_mom_calc_advance(mc, &plocal, &clocal, f, m0);
  struct gkyl_array *m0_ho = m0;
  if (use_gpu) {
    m0_ho = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, clocal_ext.volume);
    gkyl_array_copy(m0_ho, m0);
  }
  const double *m = gkyl_array_cfetch(m0_ho, 0);
  double err = 0.0;
  for (int j = 0; j < cbasis.num_basis; ++j) {
    err = fmax(err, fabs(m[j] - ref[j]));
  }

  if (use_gpu) {
    gkyl_array_release(f);
    gkyl_array_release(m0_ho);
  }
  gkyl_array_release(f_ho);
  gkyl_array_release(m0);
  gkyl_mom_calc_release(mc);
  gkyl_mom_type_release(mt);
  gkyl_array_release(hamil);
  gkyl_array_release(hinv);
  gkyl_vlasov_velocity_map_release(vmap);
  return err;
}

static void
test_hybrid_m0(int cdim, enum gkyl_basis_type ctype, bool use_gpu)
{
  double err = hybrid_m0_error(cdim, ctype, use_gpu);
  TEST_CHECK(err < 1e-12);
  TEST_MSG(
    "%dx3v %s-typed conf basis: max|M0 - quadrature| = %g", cdim,
    ctype == GKYL_BASIS_MODAL_TENSOR ? "tensor" : "Serendipity", err
  );
}

// The 2x3v and 3x3v hybrid kernels are optional build sets; only test what was built.
static void
t_ten_1x3v(void)
{
  test_hybrid_m0(1, GKYL_BASIS_MODAL_TENSOR, false);
}
static void
t_ser_1x3v(void)
{
  test_hybrid_m0(1, GKYL_BASIS_MODAL_SERENDIPITY, false);
}
#ifdef GKYL_BUILD_VLASOV_HYB_2X3V
static void
t_ten_2x3v(void)
{
  test_hybrid_m0(2, GKYL_BASIS_MODAL_TENSOR, false);
}
static void
t_ser_2x3v(void)
{
  test_hybrid_m0(2, GKYL_BASIS_MODAL_SERENDIPITY, false);
}
#ifdef GKYL_HAVE_CUDA
static void
t_ser_2x3v_cu(void)
{
  test_hybrid_m0(2, GKYL_BASIS_MODAL_SERENDIPITY, true);
}
#endif
#endif
#ifdef GKYL_BUILD_VLASOV_HYB_3X3V
static void
t_ten_3x3v(void)
{
  test_hybrid_m0(3, GKYL_BASIS_MODAL_TENSOR, false);
}
static void
t_ser_3x3v(void)
{
  test_hybrid_m0(3, GKYL_BASIS_MODAL_SERENDIPITY, false);
}
#ifdef GKYL_HAVE_CUDA
static void
t_ten_3x3v_cu(void)
{
  test_hybrid_m0(3, GKYL_BASIS_MODAL_TENSOR, true);
}
static void
t_ser_3x3v_cu(void)
{
  test_hybrid_m0(3, GKYL_BASIS_MODAL_SERENDIPITY, true);
}
#endif
#endif

TEST_LIST = {
  {"hybrid_m0_tensor_typed_conf_1x3v", t_ten_1x3v},
  {"hybrid_m0_ser_typed_conf_1x3v", t_ser_1x3v},
#ifdef GKYL_BUILD_VLASOV_HYB_2X3V
  {"hybrid_m0_tensor_typed_conf_2x3v", t_ten_2x3v},
  {"hybrid_m0_ser_typed_conf_2x3v", t_ser_2x3v},
#ifdef GKYL_HAVE_CUDA
  {"cu_hybrid_m0_ser_typed_conf_2x3v", t_ser_2x3v_cu},
#endif
#endif
#ifdef GKYL_BUILD_VLASOV_HYB_3X3V
  {"hybrid_m0_tensor_typed_conf_3x3v", t_ten_3x3v},
  {"hybrid_m0_ser_typed_conf_3x3v", t_ser_3x3v},
#ifdef GKYL_HAVE_CUDA
  {"cu_hybrid_m0_tensor_typed_conf_3x3v", t_ten_3x3v_cu},
  {"cu_hybrid_m0_ser_typed_conf_3x3v", t_ser_3x3v_cu},
#endif
#endif
  {NULL, NULL}
};
