// Tests for the PKPM variables updater (gkyl_dg_calc_pkpm_vars).
//
// The source test uses inputs that are uniform in configuration space (only the cell-average
// / 0th orthonormal basis component is non-zero in each block): for a 1D p1 conf basis the
// orthonormal cell-average basis is psi_0 = 1/sqrt(2), so a field with physical value V is
// represented by setting component 0 to V*sqrt(2).
//
// Methods tested:
//   prim, grad_ln_n : the nodal quotients u = M/rho and grad(ln n) = grad(rho)/rho are the
//            weak divisions at p=1 (rho times the quotient, projected, is the numerator).
//   source : Lorentz force q/m (rho E + M x B) momentum source and M.E work.
//
#include <acutest.h>
#include <math.h>
#include <string.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

// The weak gradient: the in-cell gradient with the whole jump of every face lifted into the cell.
static void
weak_grad(
  struct gkyl_dg_calc_pkpm_vars *vars, const struct gkyl_range *range, const struct gkyl_array *q,
  struct gkyl_array *grad
)
{
  struct gkyl_array *blend = gkyl_array_new(GKYL_DOUBLE, 3, grad->size);
  gkyl_array_clear(blend, 1.0);
  gkyl_dg_calc_pkpm_vars_in_cell_grad(vars, range, blend, q, grad);
  gkyl_array_release(blend);
}

static struct gkyl_array *
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// Identity configuration-space position map (the updaters always take one).
static struct gkyl_vlasov_position_map *
identity_pos_map(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *clocal, const struct gkyl_range *clocal_ext,
  const struct gkyl_basis *cbasis
)
{
  struct gkyl_vlasov_position_map_inp inp[GKYL_MAX_CDIM] = {0};
  return gkyl_vlasov_position_map_new(cgrid, clocal, clocal_ext, cbasis, inp, false);
}

// Set the cell-average (physical) value of block `blk` (each block has
// `nb` basis components) to `val` for *every* cell.  For a 1D p1 orthonormal
// basis psi_0 = 1/sqrt(2), so component 0 must hold val*sqrt(2).
static void
set_const_block(struct gkyl_array *arr, int blk, int nb, double val)
{
  gkyl_array_shiftc(arr, val * sqrt(2.0), blk * nb);
}

// Smooth positive state [rho, rho ux, rho uy, rho uz, E] for the quotient test.
static void
eval_quot_fluid(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  double rho = 1.5 + 0.6 * sin(2.0 * x) * cos(y) + 0.3 * cos(3.0 * z);
  fout[0] = rho;
  fout[1] = rho * (0.4 + 0.5 * cos(x + y));
  fout[2] = 0.7 * sin(y - z) + 0.2 * x;
  fout[3] = -0.3 + 0.4 * x * y;
  fout[4] = 3.0 + rho;
}

// The nodal quotients are the weak divisions at p=1: for a non-uniform positive rho the p=1
// projection of rho times the quotient reproduces the numerator, rho u = M (the flow of the
// primitive variables) and rho grad(ln n) = grad(rho) (the weak gradient of rho).
static void
test_quotients(int cdim)
{
  double lower[3] = {-1.0, -1.0, -1.0}, upper[3] = {1.0, 1.0, 1.0};
  int cells[3] = {5, 3, 2}, ghost[3] = {1, 1, 1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  int nb = cbasis.num_basis;
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct gkyl_vlasov_position_map *pm = identity_pos_map(&grid, &local, &local_ext, &cbasis);
  struct gkyl_dg_calc_pkpm_vars *up = gkyl_dg_calc_pkpm_vars_new(&grid, &cbasis, pm, false);

  struct gkyl_array *fluid = mkarr(5 * nb, local_ext.volume), *u = mkarr(3 * nb, local_ext.volume);
  struct gkyl_array *xi = mkarr(3 * nb, local_ext.volume), *rho3 = mkarr(3 * nb, local_ext.volume);
  struct gkyl_array *grad_rho = mkarr(cdim * 3 * nb, local_ext.volume), *gln = mkarr(cdim * nb, local_ext.volume);
  struct gkyl_array *prod = mkarr(nb, local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &cbasis, 4, 5, eval_quot_fluid, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &local_ext, fluid);
  gkyl_proj_on_basis_release(proj);

  int num_null = gkyl_dg_calc_pkpm_vars_prim(up, &local_ext, fluid, u, xi);
  TEST_CHECK(num_null == 0);
  gkyl_array_clear(rho3, 0.0);
  gkyl_array_set_offset(rho3, 1.0, fluid, 0);
  gkyl_array_clear(grad_rho, 0.0);
  weak_grad(up, &local, rho3, grad_rho);
  num_null = gkyl_dg_calc_pkpm_vars_grad_ln_n(up, &local, fluid, grad_rho, gln);
  TEST_CHECK(num_null == 0);

  double err_u = 0.0, max_m = 0.0, err_g = 0.0, max_g = 0.0;
  struct gkyl_range_iter iter;
  for (int i = 0; i < 3 + cdim; ++i) {
    // rho (block 0 of the fluid) times u_i or grad(ln n)_d, projected onto p=1
    if (i < 3) {
      gkyl_dg_mul_op_range(&cbasis, 0, prod, 0, fluid, i, u, &local);
    }
    else {
      gkyl_dg_mul_op_range(&cbasis, 0, prod, 0, fluid, i - 3, gln, &local);
    }
    gkyl_range_iter_init(&iter, &local);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&local, iter.idx);
      const double *pr = gkyl_array_cfetch(prod, loc);
      const double *num = i < 3 ? (const double *)gkyl_array_cfetch(fluid, loc) + (1 + i) * nb :
        (const double *)gkyl_array_cfetch(grad_rho, loc) + (i - 3) * 3 * nb;
      for (int k = 0; k < nb; ++k) {
        if (i < 3) {
          err_u = fmax(err_u, fabs(pr[k] - num[k]));
          max_m = fmax(max_m, fabs(num[k]));
        }
        else {
          err_g = fmax(err_g, fabs(pr[k] - num[k]));
          max_g = fmax(max_g, fabs(num[k]));
        }
      }
    }
  }
  TEST_CHECK(err_u < 1e-14 * max_m);
  TEST_MSG("cdim=%d: rho u - M = %g (|M| %g)", cdim, err_u, max_m);
  TEST_CHECK(max_g > 0.0 && err_g < 1e-14 * max_g);
  TEST_MSG("cdim=%d: rho grad(ln n) - grad(rho) = %g (|grad rho| %g)", cdim, err_g, max_g);

  gkyl_array_release(prod);
  gkyl_array_release(gln);
  gkyl_array_release(grad_rho);
  gkyl_array_release(rho3);
  gkyl_array_release(xi);
  gkyl_array_release(u);
  gkyl_array_release(fluid);
  gkyl_dg_calc_pkpm_vars_release(up);
  gkyl_vlasov_position_map_release(pm);
}

void test_quotients_1x() { test_quotients(1); }
void test_quotients_2x() { test_quotients(2); }
void test_quotients_3x() { test_quotients(3); }

void
test_source_1x_p1()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {6};
  int cdim = 1;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  int nb = cbasis.num_basis;

  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct gkyl_vlasov_position_map *pm = identity_pos_map(&grid, &local, &local_ext, &cbasis);


  struct gkyl_dg_calc_pkpm_vars *up =
    gkyl_dg_calc_pkpm_vars_new(&grid, &cbasis, pm, false);

  // qmem = q/m*[Ex,Ey,Ez,Bx,By,Bz,phi,psi] (8 blocks).
  struct gkyl_array *qmem = mkarr(8 * nb, local_ext.volume);
  struct gkyl_array *euler = mkarr(5 * nb, local_ext.volume);
  struct gkyl_array *rhs = mkarr(5 * nb, local_ext.volume);

  // Choose B along z only, E along x only, so source = rho*E + rhou x B.
  double rho = 3.0;
  double Ex = 2.0, Ey = 0.0, Ez = 0.0;
  double Bx = 0.0, By = 0.0, Bz = 5.0;
  double rhoux = 1.0, rhouy = 4.0, rhouz = -2.0;

  gkyl_array_clear(qmem, 0.0);
  set_const_block(qmem, 0, nb, Ex);
  set_const_block(qmem, 1, nb, Ey);
  set_const_block(qmem, 2, nb, Ez);
  set_const_block(qmem, 3, nb, Bx);
  set_const_block(qmem, 4, nb, By);
  set_const_block(qmem, 5, nb, Bz);

  gkyl_array_clear(euler, 0.0);
  set_const_block(euler, 0, nb, rho);
  set_const_block(euler, 1, nb, rhoux);
  set_const_block(euler, 2, nb, rhouy);
  set_const_block(euler, 3, nb, rhouz);

  gkyl_array_clear(rhs, 0.0);
  gkyl_dg_calc_pkpm_vars_source(up, &local, qmem, euler, rhs);

  // Expected Lorentz-force cell-average source:
  //   out_x = rho*Ex + (rhouy*Bz - rhouz*By)
  //   out_y = rho*Ey + (rhouz*Bx - rhoux*Bz)
  //   out_z = rho*Ez + (rhoux*By - rhouy*Bx)
  double sx = rho * Ex + (rhouy * Bz - rhouz * By);
  double sy = rho * Ey + (rhouz * Bx - rhoux * Bz);
  double sz = rho * Ez + (rhoux * By - rhouy * Bx);
  double sE = rhoux * Ex + rhouy * Ey + rhouz * Ez; // energy source M.E
  double s2 = sqrt(2.0);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *o = gkyl_array_cfetch(rhs, loc);
    TEST_CHECK(gkyl_compare(o[0 * nb], 0.0, 1e-12));
    TEST_CHECK(gkyl_compare(o[1 * nb] / s2, sx, 1e-12));
    TEST_CHECK(gkyl_compare(o[2 * nb] / s2, sy, 1e-12));
    TEST_CHECK(gkyl_compare(o[3 * nb] / s2, sz, 1e-12));
    TEST_CHECK(gkyl_compare(o[4 * nb] / s2, sE, 1e-12));
  }

  gkyl_array_release(rhs);
  gkyl_array_release(euler);
  gkyl_array_release(qmem);
  gkyl_dg_calc_pkpm_vars_release(up);
  gkyl_vlasov_position_map_release(pm);
}

// Vertex values of a p=1 serendipity expansion.
static void
vertex_values(const struct gkyl_basis *b, const double *c, double *v)
{
  int cdim = b->ndim, nv = 1 << cdim;
  for (int k = 0; k < nv; ++k) {
    double z[3];
    for (int d = 0; d < cdim; ++d) {
      z[d] = (k >> d) & 1 ? 1.0 : -1.0;
    }
    v[k] = b->eval_expand(z, c);
  }
}

// The positivity limiter: a cell with admissible vertices is untouched; a cell whose
// vertex density and one whose vertex internal energy are inadmissible get their slopes
// scaled to the floor with the cell averages preserved; a cell with an inadmissible
// average is left as is and counted.
static void
test_positivity(int cdim)
{
  int poly_order = 1;
  double lower[3] = {-1.0, -1.0, -1.0}, upper[3] = {1.0, 1.0, 1.0};
  int cells[3] = {4, 2, 2}, ghost[3] = {1, 1, 1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  int nb = cbasis.num_basis, nv = 1 << cdim;
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct gkyl_vlasov_position_map *pm = identity_pos_map(&grid, &local, &local_ext, &cbasis);
  struct gkyl_dg_calc_pkpm_vars *vars = gkyl_dg_calc_pkpm_vars_new(&grid, &cbasis, pm, false);

  struct gkyl_array *fluid = mkarr(5 * nb, local_ext.volume), *count = mkarr(4, local_ext.volume);
  gkyl_array_clear(fluid, 0.0);
  gkyl_array_clear(count, 0.0);
  double c0 = sqrt(pow(2.0, cdim)); // constant-mode coefficient of a unit cell average
  // Every cell starts smooth and admissible; along x (other indices at their lower bound)
  // cells 2, 3 and 4 become steep density, steep energy (large momentum slope) and
  // inadmissible average.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(fluid, gkyl_range_idx(&local_ext, iter.idx));
    f[0] = c0;
    f[nb] = 0.3 * c0;
    f[4 * nb] = c0 * (0.5 + 0.5 * 0.3 * 0.3);
    f[1] = 0.1 * c0;
    f[4 * nb + 1] = 0.05 * c0;
  }
  int idx[3] = {1, 1, 1};
  double eps = 1.0e-2;
  for (int i = 2; i <= 4; ++i) {
    idx[0] = i;
    double *f = gkyl_array_fetch(fluid, gkyl_range_idx(&local_ext, idx));
    double rho = 1.0, ux = 0.3; // averages: rho, flow
    f[1] = 0.0;
    f[4 * nb + 1] = 0.0;
    if (i == 2) {
      f[1] = 1.2 * c0; // vertex density of about -1
    } else if (i == 3) {
      f[nb + 1] = 1.5 * c0; // vertex |u| of about 2.9 with the same E: negative internal energy
    } else {
      f[4 * nb] = c0 * rho * (0.5 * ux * ux - 0.1); // average internal energy negative
    }
  }
  double avg_before[4][5];
  for (int i = 1; i <= 4; ++i) {
    idx[0] = i;
    const double *f = gkyl_array_cfetch(fluid, gkyl_range_idx(&local_ext, idx));
    for (int c = 0; c < 5; ++c) {
      avg_before[i - 1][c] = f[c * nb];
    }
  }
  int num_rho = -1, num_e = -1, num_xi = -1, num_avg = -1;
  gkyl_dg_calc_pkpm_vars_positivity(vars, &local, eps, fluid, count, &num_rho, &num_e, &num_xi, &num_avg);
  TEST_CHECK(num_rho == 1);
  TEST_CHECK(num_e == 2); // the steep-density cell also needs the energy scaling after step 1
  TEST_CHECK(num_xi == 0); // states varying along x only: steps 1-2 already bound xi (test_positivity_xi)
  TEST_CHECK(num_avg == 1);

  for (int i = 1; i <= 4; ++i) {
    idx[0] = i;
    long loc = gkyl_range_idx(&local_ext, idx);
    const double *f = gkyl_array_cfetch(fluid, loc), *cnt = gkyl_array_cfetch(count, loc);
    // averages preserved
    for (int c = 0; c < 5; ++c) {
      TEST_CHECK(f[c * nb] == avg_before[i - 1][c]);
    }
    TEST_CHECK(cnt[2] == 0.0); // no thermal-speed scaling
    if (i == 4) {
      TEST_CHECK(cnt[3] == 1.0 && cnt[0] == 0.0 && cnt[1] == 0.0);
      continue;
    }
    TEST_CHECK(cnt[3] == 0.0);
    double rho_v[8], M_v[3][8], E_v[8];
    vertex_values(&cbasis, &f[0], rho_v);
    for (int d = 0; d < 3; ++d) {
      vertex_values(&cbasis, &f[(1 + d) * nb], M_v[d]);
    }
    vertex_values(&cbasis, &f[4 * nb], E_v);
    // internal energies per volume, E - |M|^2/(2 rho), of the average and of the vertices
    double rho_avg = f[0] / c0, e_avg = f[4 * nb] / c0 - 0.5 * pow(f[nb] / c0, 2) / rho_avg;
    double rho_min = 1.0, e_min = 1.0;
    for (int k = 0; k < nv; ++k) {
      double e_v = E_v[k] - 0.5 * (M_v[0][k] * M_v[0][k] + M_v[1][k] * M_v[1][k] + M_v[2][k] * M_v[2][k]) / rho_v[k];
      rho_min = fmin(rho_min, rho_v[k] / rho_avg);
      e_min = fmin(e_min, e_v / e_avg);
    }
    TEST_CHECK(rho_min >= eps * (1.0 - 1e-12));
    TEST_CHECK(e_min >= eps * (1.0 - 1e-10));
    TEST_MSG("cdim=%d cell %d: vertex rho/avg min %g, e/avg min %g", cdim, i, rho_min, e_min);
    if (i == 1) {
      TEST_CHECK(cnt[0] == 0.0 && cnt[1] == 0.0);
      TEST_CHECK(f[1] == 0.1 * c0); // untouched
    } else if (i == 2) {
      // density scaled to the floor, then all slopes scaled again for the energy
      TEST_CHECK(cnt[0] == 1.0 && cnt[1] == 1.0);
      TEST_CHECK(rho_min > eps);
      TEST_CHECK(gkyl_compare(e_min, eps, 1e-8));
    } else {
      TEST_CHECK(cnt[0] == 0.0 && cnt[1] == 1.0);
      TEST_CHECK(gkyl_compare(e_min, eps, 1e-8)); // scaled exactly to the floor
    }
  }

  gkyl_array_release(count);
  gkyl_array_release(fluid);
  gkyl_dg_calc_pkpm_vars_release(vars);
  gkyl_vlasov_position_map_release(pm);
}

void test_positivity_1x() { test_positivity(1); }
void test_positivity_2x() { test_positivity(2); }
void test_positivity_3x() { test_positivity(3); }

// Worst states found by an adversarial search (agent-scratch/gen/positivity/xi_bound_adv.py) for the vertex
// extrapolation of xi after steps 1-2 of the limiter: q(x) = q_avg + sum_s c_s prod_{d in s} x_d with
// q_avg = [rho, M, E] = [1, 0, 0, 0, 1]. In 1x the extrapolation stays at 1.2% of xi(q_avg) (above eps = 1e-2):
// steps 1-2 alone bound it. In 2x it reaches -0.18 xi(q_avg), and the 3x state is the same without z dependence.
static void
eval_xi_worst(int cdim, const double *x, double *q)
{
  const double avg[5] = {1.0, 0.0, 0.0, 0.0, 1.0};
  const double c1[5] = {1.9, 0.0, 0.0, 0.0, -0.99}; // 1x: c_x
  const double c2[3][5] = {
    {-0.645, -0.293, 0.076, 0.177, 0.053}, // c_x
    {-0.657, -0.278, -0.112, -0.034, 0.142}, // c_y
    {0.978, -0.141, -0.149, 0.066, -0.781}, // c_xy
  };
  for (int c = 0; c < 5; ++c) {
    q[c] = avg[c] + (cdim == 1 ? c1[c] * x[0] : c2[0][c] * x[0] + c2[1][c] * x[1] + c2[2][c] * x[0] * x[1]);
  }
}

// The thermal-speed step of the limiter: in every dimension the p=1 thermal speed of
// gkyl_dg_calc_pkpm_vars_prim on the limited state is >= eps xi(q_avg) at every vertex (so on every face, where
// the configuration-space flux reads it), with the cell averages preserved. In 1x the step does not act (steps
// 1-2 bound xi); in 2x and 3x it acts and puts the smallest vertex value on the floor.
static void
test_positivity_xi(int cdim)
{
  double lower[3] = {-1.0, -1.0, -1.0}, upper[3] = {1.0, 1.0, 1.0};
  int cells[3] = {1, 1, 1}, ghost[3] = {1, 1, 1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  int nb = cbasis.num_basis, nv = 1 << cdim;
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct gkyl_vlasov_position_map *pm = identity_pos_map(&grid, &local, &local_ext, &cbasis);
  struct gkyl_dg_calc_pkpm_vars *vars = gkyl_dg_calc_pkpm_vars_new(&grid, &cbasis, pm, false);
  struct gkyl_array *fluid = mkarr(5 * nb, local_ext.volume), *count = mkarr(4, local_ext.volume);
  struct gkyl_array *u = mkarr(3 * nb, local_ext.volume), *xi = mkarr(3 * nb, local_ext.volume);
  gkyl_array_clear(fluid, 0.0);
  gkyl_array_clear(count, 0.0);

  // Project the multilinear state exactly with the 2-point Gauss-Legendre rule (unit weights).
  int idx[3] = {1, 1, 1};
  long loc = gkyl_range_idx(&local_ext, idx);
  double *f = gkyl_array_fetch(fluid, loc);
  for (int n = 0; n < nv; ++n) {
    double z[3], q[5], phi[8];
    for (int d = 0; d < cdim; ++d) {
      z[d] = ((n >> d) & 1 ? 1.0 : -1.0) / sqrt(3.0);
    }
    eval_xi_worst(cdim, z, q);
    cbasis.eval(z, phi);
    for (int c = 0; c < 5; ++c) {
      for (int k = 0; k < nb; ++k) {
        f[c * nb + k] += q[c] * phi[k];
      }
    }
  }
  double avg_before[5];
  for (int c = 0; c < 5; ++c) {
    avg_before[c] = f[c * nb];
  }

  double eps = 1.0e-2;
  int num_rho = -1, num_e = -1, num_xi = -1, num_avg = -1;
  gkyl_dg_calc_pkpm_vars_positivity(vars, &local, eps, fluid, count, &num_rho, &num_e, &num_xi, &num_avg);
  TEST_CHECK(num_avg == 0);
  TEST_CHECK(num_xi == (cdim == 1 ? 0 : 1));
  TEST_MSG("cdim=%d: num_rho %d, num_e %d, num_xi %d", cdim, num_rho, num_e, num_xi);
  for (int c = 0; c < 5; ++c) {
    TEST_CHECK(f[c * nb] == avg_before[c]); // averages preserved
  }

  int num_null = gkyl_dg_calc_pkpm_vars_prim(vars, &local, fluid, u, xi);
  TEST_CHECK(num_null == 0);
  const double *xi_c = gkyl_array_cfetch(xi, loc);
  double xv[8], c0 = sqrt(pow(2.0, cdim));
  vertex_values(&cbasis, xi_c, xv);
  // xi of the cell average: q_avg = [1, 0, 0, 0, 1], e_avg = 1
  double xi_avg = sqrt(2.0 / 3.0), xi_vmin = xv[0];
  for (int k = 1; k < nv; ++k) {
    xi_vmin = fmin(xi_vmin, xv[k]);
  }
  TEST_CHECK(xi_vmin >= eps * xi_avg * (1.0 - 1e-10));
  if (cdim > 1) {
    TEST_CHECK(gkyl_compare(xi_vmin / xi_avg, eps, 1e-8)); // scaled to the floor (bisection)
  } else {
    TEST_CHECK(xi_vmin / xi_avg > eps); // above the floor without step 3 (about 1.2% of xi_avg)
  }
  TEST_MSG("cdim=%d: smallest vertex xi / xi_avg = %.10g (floor %g); rho avg %g", cdim, xi_vmin / xi_avg, eps,
    f[0] / c0);

  gkyl_array_release(xi);
  gkyl_array_release(u);
  gkyl_array_release(count);
  gkyl_array_release(fluid);
  gkyl_dg_calc_pkpm_vars_release(vars);
  gkyl_vlasov_position_map_release(pm);
}

void test_positivity_xi_1x() { test_positivity_xi(1); }
void test_positivity_xi_2x() { test_positivity_xi(2); }
void test_positivity_xi_3x() { test_positivity_xi(3); }

TEST_LIST = {
  {"test_positivity_1x", test_positivity_1x},
  {"test_positivity_2x", test_positivity_2x},
  {"test_positivity_3x", test_positivity_3x},
  {"test_positivity_xi_1x", test_positivity_xi_1x},
  {"test_positivity_xi_2x", test_positivity_xi_2x},
  {"test_positivity_xi_3x", test_positivity_xi_3x},
  {"test_quotients_1x", test_quotients_1x},
  {"test_quotients_2x", test_quotients_2x},
  {"test_quotients_3x", test_quotients_3x},
  {"source_1x_p1", test_source_1x_p1},
  {NULL, NULL}
};
