// Systematic test of the special-relativistic auxiliary variable
// computations (gamma = sqrt(1+p^2) via the SR Hamiltonian, rest-frame
// density, bulk four-velocity Lorentz factor, and the relativistic pressure),
// scanned across polynomial order, basis type, uniform vs non-uniform
// momentum grids, and CPU vs GPU.
//
// The distribution is a drifting Maxwell-Juttner projected with the LTE
// projection and we compute:
//   M0 = GammaV*n0, M1 = n0*u  =>  sr_vars_n recovers the rest-frame n0;
//   GammaV_sq = 1 + u^2 holds weakly (exact modal identity);
//   the relativistic ideal gas obeys P = n0*T exactly.
// The gamma check is nodal: the SR Hamiltonian is built at its construction
// nodes from the stored momentum map, so hamil(eta_node) = sqrt(1+vmap^2)
// there to machine precision (geometry is known nodally; exact by definition).
#include <math.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_const.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_vlasov_velocity_map.h>

#define PMAX 16.0 // Mapped momentum extent.

// Maxwell-Juttner parameters: rest-frame density, spatial four-velocity
// (u = GammaV*V, leading vdim components used), and temperature.
#define N0 1.0
#define TEMP 0.5
static const double udrift[3] = { 0.5, -0.3, 0.2 };

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array* a = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                 : gkyl_array_new(GKYL_DOUBLE, nc, size);
#else
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
#endif
  return a;
}

static void
eval_quad_vmap(double t, const double *vc, double *vp, void *ctx)
{
  vp[0] = vc[0] < 0.0 ? -PMAX*vc[0]*vc[0] : PMAX*vc[0]*vc[0];
}

// Evaluate the per-direction momentum map (stored as a degenerate cubic) at
// logical coordinate z in [-1,1].
static inline double
vmap_eval(const struct gkyl_basis *b1, const double *vmap_c, int d, double z)
{
  double zv[1] = { z };
  return b1->eval_expand(zv, &vmap_c[4*d]);
}

static void
test_sr_vars(int vdim, int poly_order, bool use_tensor, bool use_nonuniform, bool use_gpu)
{
  int cdim = 1, pdim = cdim+vdim;

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM];
  lower[0] = 0.0; upper[0] = 1.0; cells[0] = 2;
  for (int d=cdim; d<pdim; ++d) {
    if (use_nonuniform) { lower[d] = -1.0; upper[d] = 1.0; cells[d] = 32; }
    else { lower[d] = -PMAX; upper[d] = PMAX; cells[d] = 64; }
  }

  struct gkyl_rect_grid phase_grid, conf_grid, vel_grid;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&vel_grid, vdim, &lower[cdim], &upper[cdim], &cells[cdim]);

  struct gkyl_basis pbasis, cbasis, vbasis;
  if (use_tensor) {
    gkyl_cart_modal_tensor(&pbasis, pdim, poly_order);
    gkyl_cart_modal_tensor(&cbasis, cdim, poly_order);
    gkyl_cart_modal_tensor(&vbasis, vdim, poly_order);
  }
  else {
    gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);
    gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
    gkyl_cart_modal_serendip(&vbasis, vdim, poly_order);
  }

  int conf_ghost[] = { 1 };
  struct gkyl_range conf_local, conf_local_ext;
  gkyl_create_grid_ranges(&conf_grid, conf_ghost, &conf_local_ext, &conf_local);

  int vel_ghost[] = { 0, 0, 0 };
  struct gkyl_range vel_local, vel_local_ext;
  gkyl_create_grid_ranges(&vel_grid, vel_ghost, &vel_local_ext, &vel_local);

  int phase_ghost[GKYL_MAX_DIM] = { 1, 0, 0, 0 };
  struct gkyl_range phase_local, phase_local_ext;
  gkyl_create_grid_ranges(&phase_grid, phase_ghost, &phase_local_ext, &phase_local);

  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  if (use_nonuniform)
    for (int d=0; d<vdim; ++d) inp_vmap[d].eval_vmap = eval_quad_vmap;
  struct gkyl_vlasov_velocity_map *vvm = gkyl_vlasov_velocity_map_new(&vel_grid,
    &vel_local, &vbasis, inp_vmap, false, use_gpu);

  // SR Hamiltonian: hamil = gamma = sqrt(1 + p^2), hamil_inv = 1/gamma.
  struct gkyl_array *gamma = mkarr(use_gpu, vbasis.num_basis, vel_local.volume);
  struct gkyl_array *gamma_inv = mkarr(use_gpu, vbasis.num_basis, vel_local.volume);
  gkyl_dg_vlasov_calc_hamil(&vel_grid, &vbasis, &vel_local,
    GKYL_MODEL_SR, vvm, gamma, gamma_inv, use_gpu);

  struct gkyl_array *gamma_ho = gkyl_array_new(GKYL_DOUBLE, vbasis.num_basis, vel_local.volume);
  struct gkyl_array *gamma_inv_ho = gkyl_array_new(GKYL_DOUBLE, vbasis.num_basis, vel_local.volume);
  gkyl_array_copy(gamma_ho, gamma);
  gkyl_array_copy(gamma_inv_ho, gamma_inv);

  // 1D cubic basis for evaluating the stored momentum map.
  struct gkyl_basis b1;
  gkyl_cart_modal_tensor(&b1, 1, 3);

  // (a) Nodal check of gamma at the Hamiltonian construction nodes: the SR
  // Hamiltonian kernels build gamma = sqrt(1 + |p|^2) from the stored map at
  // exactly the num_basis nodal points of the velocity basis (node_list), so
  // the construction interpolates there for every vdim, order, and basis
  // type, and the agreement is to machine precision; gamma_inv likewise.
  {
    double nodes[160*3];
    vbasis.node_list(nodes);
    struct gkyl_range_iter viter;
    gkyl_range_iter_init(&viter, &vel_local);
    while (gkyl_range_iter_next(&viter)) {
      long vloc = gkyl_range_idx(&vel_local, viter.idx);
      const double *vmap_c = gkyl_array_cfetch(vvm->vmap_host, vloc);
      const double *g_c = gkyl_array_cfetch(gamma_ho, vloc);
      const double *gi_c = gkyl_array_cfetch(gamma_inv_ho, vloc);
      for (int i=0; i<vbasis.num_basis; ++i) {
        const double *zn = &nodes[i*vdim];
        double psq = 0.0;
        for (int d=0; d<vdim; ++d) {
          double pd = vmap_eval(&b1, vmap_c, d, zn[d]);
          psq += pd*pd;
        }
        double gex = sqrt(1.0 + psq);
        TEST_CHECK( gkyl_compare_double(vbasis.eval_expand(zn, g_c), gex, 1e-13) );
        TEST_MSG("gamma %dv p%d %s %s cell=%d node=%d: %.15e vs %.15e", vdim, poly_order,
          use_tensor ? "tensor" : "ser", use_nonuniform ? "nonuni" : "uni",
          viter.idx[0], i, vbasis.eval_expand(zn, g_c), gex);
        TEST_CHECK( gkyl_compare_double(vbasis.eval_expand(zn, gi_c), 1.0/gex, 1e-13) );
      }
    }
  }

  // Project the drifting Maxwell-Juttner with the LTE updater.
  struct gkyl_array *moms_ho = gkyl_array_new(GKYL_DOUBLE, (vdim+2)*cbasis.num_basis, conf_local_ext.volume);
  gkyl_array_clear(moms_ho, 0.0);
  for (int ix=1; ix<=cells[0]; ++ix) {
    double *m = gkyl_array_fetch(moms_ho, gkyl_range_idx(&conf_local, (int[]) { ix }));
    m[0] = N0*sqrt(2.0);
    for (int d=0; d<vdim; ++d)
      m[(1+d)*cbasis.num_basis] = udrift[d]*sqrt(2.0); // u_i (spatial four-velocity)
    m[(vdim+1)*cbasis.num_basis] = TEMP*sqrt(2.0);
  }
  struct gkyl_array *moms = use_gpu ? mkarr(use_gpu, moms_ho->ncomp, moms_ho->size) : gkyl_array_acquire(moms_ho);
  gkyl_array_copy(moms, moms_ho);

  struct gkyl_vlasov_lte_proj_on_basis_inp inp_proj = {
    .phase_grid = &phase_grid,
    .vel_grid = &vel_grid,
    .conf_basis = &cbasis,
    .vel_basis = &vbasis,
    .phase_basis = &pbasis,
    .conf_range = &conf_local,
    .conf_range_ext = &conf_local_ext,
    .vel_range = &vel_local,
    .phase_range = &phase_local,
    .vel_map = vvm,
    .hamil_range = &vel_local,
    .hamil = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .hamil_id = gkyl_hamil_id_from_model_id(GKYL_MODEL_SR),
    .use_gpu = use_gpu,
  };
  struct gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_proj);
  struct gkyl_array *Jf = mkarr(use_gpu, pbasis.num_basis, phase_local_ext.volume);
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &phase_local, &conf_local, moms, Jf);

  // Lab-frame M0 and M1i from the moment kernels (SR: M1i = int (p/gamma) Jf).
  struct gkyl_mom_vlasov_inp inp_mom = {
    .conf_basis = &cbasis,
    .phase_basis = &pbasis,
    .vel_range = &vel_local,
    .hamil_range = &vel_local,
    .hamil = gamma,
    .model_id = GKYL_MODEL_SR,
    .hamil_id = gkyl_hamil_id_from_model_id(GKYL_MODEL_SR),
    .mom_type = GKYL_F_MOMENT_M0,
    .use_gpu = use_gpu,
    .vel_map = vvm,
  };
  struct gkyl_mom_type *m0_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M1;
  struct gkyl_mom_type *m1_t = gkyl_mom_vlasov_inew(&inp_mom);
  gkyl_mom_calc *m0_calc = gkyl_mom_calc_new(&phase_grid, m0_t, use_gpu);
  gkyl_mom_calc *m1_calc = gkyl_mom_calc_new(&phase_grid, m1_t, use_gpu);

  struct gkyl_array *M0 = mkarr(use_gpu, cbasis.num_basis, conf_local_ext.volume);
  struct gkyl_array *M1i = mkarr(use_gpu, vdim*cbasis.num_basis, conf_local_ext.volume);
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    gkyl_mom_calc_advance_cu(m0_calc, &phase_local, &conf_local, Jf, M0);
    gkyl_mom_calc_advance_cu(m1_calc, &phase_local, &conf_local, Jf, M1i);
  }
  else {
    gkyl_mom_calc_advance(m0_calc, &phase_local, &conf_local, Jf, M0);
    gkyl_mom_calc_advance(m1_calc, &phase_local, &conf_local, Jf, M1i);
  }
#else
  gkyl_mom_calc_advance(m0_calc, &phase_local, &conf_local, Jf, M0);
  gkyl_mom_calc_advance(m1_calc, &phase_local, &conf_local, Jf, M1i);
#endif

  // SR variables updater.
  struct gkyl_dg_calc_sr_vars *sr_vars = gkyl_dg_calc_sr_vars_new(&phase_grid, &vel_grid,
    &cbasis, &vbasis, &conf_local, &vel_local, vvm, use_gpu);

  // (b) Rest-frame density: n = sqrt(M0^2 - M1^2) must recover n0.
  struct gkyl_array *n_rest = mkarr(use_gpu, cbasis.num_basis, conf_local_ext.volume);
  gkyl_dg_calc_sr_vars_n(sr_vars, M0, M1i, n_rest);
  struct gkyl_array *n_ho = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, conf_local_ext.volume);
  gkyl_array_copy(n_ho, n_rest);
  double usq0 = 0.0;
  for (int d=0; d<vdim; ++d) usq0 += udrift[d]*udrift[d];
  double GammaV0 = sqrt(1.0 + usq0);
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, &conf_local);
  while (gkyl_range_iter_next(&citer)) {
    const double *nc = gkyl_array_cfetch(n_ho, gkyl_range_idx(&conf_local, citer.idx));
    TEST_CHECK( gkyl_compare_double(nc[0]/sqrt(2.0), N0, 1e-4) );
    TEST_MSG("n_rest p%d %s %s cell=%d: %.15e vs n0=%.15e (rel %.2e)", poly_order,
      use_tensor ? "tensor" : "ser", use_nonuniform ? "nonuni" : "uni",
      citer.idx[0], nc[0]/sqrt(2.0), N0, nc[0]/sqrt(2.0)/N0-1.0);
    for (int k=1; k<cbasis.num_basis; ++k)
      TEST_CHECK( fabs(nc[k]) < 1e-6 ); // constant moments => constant n.
  }

  // Bulk four-velocity u_i = M1i/n by weak division.
  struct gkyl_array *u_i = mkarr(use_gpu, vdim*cbasis.num_basis, conf_local_ext.volume);
  gkyl_dg_bin_op_mem *div_mem = use_gpu ? gkyl_dg_bin_op_mem_cu_dev_new(conf_local.volume, cbasis.num_basis)
                                        : gkyl_dg_bin_op_mem_new(conf_local.volume, cbasis.num_basis);
  for (int d=0; d<vdim; ++d)
    #ifdef GKYL_HAVE_CUDA
    if (use_gpu)
      gkyl_dg_div_op_range_cu(div_mem, &cbasis, d, u_i, d, M1i, 0, n_rest, &conf_local);
    else
      gkyl_dg_div_op_range(div_mem, &cbasis, d, u_i, d, M1i, 0, n_rest, &conf_local);
#else
    gkyl_dg_div_op_range(div_mem, &cbasis, d, u_i, d, M1i, 0, n_rest, &conf_local);
#endif

  // (c) GammaV algebra: u_i_sq is the weak square of u_i and
  // GammaV_sq = 1 + sum_i u_i_sq must hold as an exact modal identity.
  struct gkyl_array *u_i_sq = mkarr(use_gpu, vdim*cbasis.num_basis, conf_local_ext.volume);
  struct gkyl_array *GammaV = mkarr(use_gpu, cbasis.num_basis, conf_local_ext.volume);
  struct gkyl_array *GammaV_sq = mkarr(use_gpu, cbasis.num_basis, conf_local_ext.volume);
  gkyl_dg_calc_sr_vars_GammaV(sr_vars, &conf_local, u_i, u_i_sq, GammaV, GammaV_sq);

  struct gkyl_array *u_i_sq_ho = gkyl_array_new(GKYL_DOUBLE, u_i_sq->ncomp, u_i_sq->size);
  struct gkyl_array *GammaV_ho = gkyl_array_new(GKYL_DOUBLE, GammaV->ncomp, GammaV->size);
  struct gkyl_array *GammaV_sq_ho = gkyl_array_new(GKYL_DOUBLE, GammaV_sq->ncomp, GammaV_sq->size);
  gkyl_array_copy(u_i_sq_ho, u_i_sq);
  gkyl_array_copy(GammaV_ho, GammaV);
  gkyl_array_copy(GammaV_sq_ho, GammaV_sq);

  gkyl_range_iter_init(&citer, &conf_local);
  while (gkyl_range_iter_next(&citer)) {
    long cloc = gkyl_range_idx(&conf_local, citer.idx);
    const double *usq = gkyl_array_cfetch(u_i_sq_ho, cloc);
    const double *gv = gkyl_array_cfetch(GammaV_ho, cloc);
    const double *gvsq = gkyl_array_cfetch(GammaV_sq_ho, cloc);
    for (int k=0; k<cbasis.num_basis; ++k) {
      double expected = k == 0 ? sqrt(2.0) : 0.0; // the constant 1 in modal form (cdim=1)
      for (int d=0; d<vdim; ++d) expected += usq[d*cbasis.num_basis+k];
      TEST_CHECK( gkyl_compare_double(gvsq[k], expected, 1e-12) );
      TEST_MSG("GammaV_sq cell=%d coeff=%d: %.15e vs 1+u_sq %.15e", citer.idx[0], k, gvsq[k], expected);
    }
    // GammaV itself: constant flow => cell value sqrt(1+|u|^2) up to the
    // accuracy of the recovered u.
    TEST_CHECK( gkyl_compare_double(gv[0]/sqrt(2.0), GammaV0, poly_order == 1 ? 5e-3 : 1e-4) );
    TEST_MSG("GammaV cell=%d: %.15e vs %.15e (rel %.2e)", citer.idx[0],
      gv[0]/sqrt(2.0), GammaV0, gv[0]/sqrt(2.0)/GammaV0-1.0);
  }

  // (d) Relativistic pressure: P = n0*T exactly for the Maxwell-Juttner.
  struct gkyl_array *pressure = mkarr(use_gpu, cbasis.num_basis, conf_local_ext.volume);
  gkyl_dg_calc_sr_vars_pressure(sr_vars, &conf_local, &phase_local,
    gamma, gamma_inv, u_i, u_i_sq, GammaV, GammaV_sq, Jf, pressure);
  struct gkyl_array *pressure_ho = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, conf_local_ext.volume);
  gkyl_array_copy(pressure_ho, pressure);

  gkyl_range_iter_init(&citer, &conf_local);
  while (gkyl_range_iter_next(&citer)) {
    const double *pc = gkyl_array_cfetch(pressure_ho, gkyl_range_idx(&conf_local, citer.idx));
    TEST_CHECK( gkyl_compare_double(pc[0]/sqrt(2.0)/(N0*TEMP), 1.0, poly_order == 1 ? 5e-2 : 2e-3) );
    TEST_MSG("pressure p%d %s %s cell=%d: P=%.15e vs n0*T=%.15e (rel %.2e)", poly_order,
      use_tensor ? "tensor" : "ser", use_nonuniform ? "nonuni" : "uni",
      citer.idx[0], pc[0]/sqrt(2.0), N0*TEMP, pc[0]/sqrt(2.0)/(N0*TEMP)-1.0);
  }

  gkyl_array_release(gamma); gkyl_array_release(gamma_inv);
  gkyl_array_release(gamma_ho); gkyl_array_release(gamma_inv_ho);
  gkyl_array_release(moms_ho); gkyl_array_release(moms);
  gkyl_array_release(Jf);
  gkyl_array_release(M0); gkyl_array_release(M1i);
  gkyl_array_release(n_rest); gkyl_array_release(n_ho);
  gkyl_array_release(u_i); gkyl_array_release(u_i_sq);
  gkyl_array_release(GammaV); gkyl_array_release(GammaV_sq);
  gkyl_array_release(u_i_sq_ho); gkyl_array_release(GammaV_ho); gkyl_array_release(GammaV_sq_ho);
  gkyl_array_release(pressure); gkyl_array_release(pressure_ho);
  gkyl_dg_bin_op_mem_release(div_mem);
  gkyl_dg_calc_sr_vars_release(sr_vars);
  gkyl_mom_calc_release(m0_calc); gkyl_mom_calc_release(m1_calc);
  gkyl_mom_type_release(m0_t); gkyl_mom_type_release(m1_t);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_vlasov_velocity_map_release(vvm);
}

// Full scan: {1v,2v,3v} x {ser p1, ser p2, tensor p2} x {uniform, nonuniform}
// x {CPU, GPU}. Tensor p1 is identical to ser p1; uniform 3v is skipped (it
// needs 64^3 velocity cells to resolve the Maxwell-Juttner that the quadratic
// map handles at 32^3).
static void t_sr_1x1v_ser_p1_uni(void) { test_sr_vars(1, 1, false, false, false); }
static void t_sr_1x1v_ser_p1_non(void) { test_sr_vars(1, 1, false, true, false); }
static void t_sr_1x1v_ser_p2_uni(void) { test_sr_vars(1, 2, false, false, false); }
static void t_sr_1x1v_ser_p2_non(void) { test_sr_vars(1, 2, false, true, false); }
static void t_sr_1x1v_ten_p2_uni(void) { test_sr_vars(1, 2, true, false, false); }
static void t_sr_1x1v_ten_p2_non(void) { test_sr_vars(1, 2, true, true, false); }
static void t_sr_1x2v_ser_p1_uni(void) { test_sr_vars(2, 1, false, false, false); }
static void t_sr_1x2v_ser_p1_non(void) { test_sr_vars(2, 1, false, true, false); }
static void t_sr_1x2v_ser_p2_uni(void) { test_sr_vars(2, 2, false, false, false); }
static void t_sr_1x2v_ser_p2_non(void) { test_sr_vars(2, 2, false, true, false); }
static void t_sr_1x2v_ten_p2_uni(void) { test_sr_vars(2, 2, true, false, false); }
static void t_sr_1x2v_ten_p2_non(void) { test_sr_vars(2, 2, true, true, false); }
static void t_sr_1x3v_ser_p1_non(void) { test_sr_vars(3, 1, false, true, false); }
static void t_sr_1x3v_ser_p2_non(void) { test_sr_vars(3, 2, false, true, false); }
static void t_sr_1x3v_ten_p2_non(void) { test_sr_vars(3, 2, true, true, false); }
#ifdef GKYL_HAVE_CUDA
static void t_sr_1x1v_ser_p1_uni_gpu(void) { test_sr_vars(1, 1, false, false, true); }
static void t_sr_1x1v_ser_p1_non_gpu(void) { test_sr_vars(1, 1, false, true, true); }
static void t_sr_1x1v_ser_p2_uni_gpu(void) { test_sr_vars(1, 2, false, false, true); }
static void t_sr_1x1v_ser_p2_non_gpu(void) { test_sr_vars(1, 2, false, true, true); }
static void t_sr_1x1v_ten_p2_uni_gpu(void) { test_sr_vars(1, 2, true, false, true); }
static void t_sr_1x1v_ten_p2_non_gpu(void) { test_sr_vars(1, 2, true, true, true); }
static void t_sr_1x2v_ser_p1_uni_gpu(void) { test_sr_vars(2, 1, false, false, true); }
static void t_sr_1x2v_ser_p1_non_gpu(void) { test_sr_vars(2, 1, false, true, true); }
static void t_sr_1x2v_ser_p2_uni_gpu(void) { test_sr_vars(2, 2, false, false, true); }
static void t_sr_1x2v_ser_p2_non_gpu(void) { test_sr_vars(2, 2, false, true, true); }
static void t_sr_1x2v_ten_p2_uni_gpu(void) { test_sr_vars(2, 2, true, false, true); }
static void t_sr_1x2v_ten_p2_non_gpu(void) { test_sr_vars(2, 2, true, true, true); }
static void t_sr_1x3v_ser_p1_non_gpu(void) { test_sr_vars(3, 1, false, true, true); }
static void t_sr_1x3v_ser_p2_non_gpu(void) { test_sr_vars(3, 2, false, true, true); }
static void t_sr_1x3v_ten_p2_non_gpu(void) { test_sr_vars(3, 2, true, true, true); }
#endif

TEST_LIST = {
  { "sr_vars_1x1v_ser_p1_uniform", t_sr_1x1v_ser_p1_uni },
  { "sr_vars_1x1v_ser_p1_nonuniform", t_sr_1x1v_ser_p1_non },
  { "sr_vars_1x1v_ser_p2_uniform", t_sr_1x1v_ser_p2_uni },
  { "sr_vars_1x1v_ser_p2_nonuniform", t_sr_1x1v_ser_p2_non },
  { "sr_vars_1x1v_tensor_p2_uniform", t_sr_1x1v_ten_p2_uni },
  { "sr_vars_1x1v_tensor_p2_nonuniform", t_sr_1x1v_ten_p2_non },
  { "sr_vars_1x2v_ser_p1_uniform", t_sr_1x2v_ser_p1_uni },
  { "sr_vars_1x2v_ser_p1_nonuniform", t_sr_1x2v_ser_p1_non },
  { "sr_vars_1x2v_ser_p2_uniform", t_sr_1x2v_ser_p2_uni },
  { "sr_vars_1x2v_ser_p2_nonuniform", t_sr_1x2v_ser_p2_non },
  { "sr_vars_1x2v_tensor_p2_uniform", t_sr_1x2v_ten_p2_uni },
  { "sr_vars_1x2v_tensor_p2_nonuniform", t_sr_1x2v_ten_p2_non },
  { "sr_vars_1x3v_ser_p1_nonuniform", t_sr_1x3v_ser_p1_non },
  { "sr_vars_1x3v_ser_p2_nonuniform", t_sr_1x3v_ser_p2_non },
  { "sr_vars_1x3v_tensor_p2_nonuniform", t_sr_1x3v_ten_p2_non },
#ifdef GKYL_HAVE_CUDA
  { "sr_vars_1x1v_ser_p1_uniform_gpu", t_sr_1x1v_ser_p1_uni_gpu },
  { "sr_vars_1x1v_ser_p1_nonuniform_gpu", t_sr_1x1v_ser_p1_non_gpu },
  { "sr_vars_1x1v_ser_p2_uniform_gpu", t_sr_1x1v_ser_p2_uni_gpu },
  { "sr_vars_1x1v_ser_p2_nonuniform_gpu", t_sr_1x1v_ser_p2_non_gpu },
  { "sr_vars_1x1v_tensor_p2_uniform_gpu", t_sr_1x1v_ten_p2_uni_gpu },
  { "sr_vars_1x1v_tensor_p2_nonuniform_gpu", t_sr_1x1v_ten_p2_non_gpu },
  { "sr_vars_1x2v_ser_p1_uniform_gpu", t_sr_1x2v_ser_p1_uni_gpu },
  { "sr_vars_1x2v_ser_p1_nonuniform_gpu", t_sr_1x2v_ser_p1_non_gpu },
  { "sr_vars_1x2v_ser_p2_uniform_gpu", t_sr_1x2v_ser_p2_uni_gpu },
  { "sr_vars_1x2v_ser_p2_nonuniform_gpu", t_sr_1x2v_ser_p2_non_gpu },
  { "sr_vars_1x2v_tensor_p2_uniform_gpu", t_sr_1x2v_ten_p2_uni_gpu },
  { "sr_vars_1x2v_tensor_p2_nonuniform_gpu", t_sr_1x2v_ten_p2_non_gpu },
  { "sr_vars_1x3v_ser_p1_nonuniform_gpu", t_sr_1x3v_ser_p1_non_gpu },
  { "sr_vars_1x3v_ser_p2_nonuniform_gpu", t_sr_1x3v_ser_p2_non_gpu },
  { "sr_vars_1x3v_tensor_p2_nonuniform_gpu", t_sr_1x3v_ten_p2_non_gpu },
#endif
  { NULL, NULL },
};
