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
//   positivity : the limiter scales the slopes of a cell with an inadmissible vertex density or
//            internal energy to the floor with the cell averages preserved; a cell with an
//            inadmissible average is left as is and counted.
//   diag : the pressure tensor, the heat flux and the pkpm variables of the output are the p=1
//            projections of their expressions in the coupling moments, b and the closures.
//
#include <acutest.h>
#include <math.h>
#include <string.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_gauss_quad_data.h>
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
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *clocal,
  const struct gkyl_range *clocal_ext, const struct gkyl_basis *cbasis
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
  struct gkyl_dg_calc_pkpm_vars *up =
    gkyl_dg_calc_pkpm_vars_inew(&(struct gkyl_dg_calc_pkpm_vars_inp){
      .conf_grid = &grid,
      .conf_basis = &cbasis,
      .pos_map = pm,
      .use_gpu = false,
    });

  struct gkyl_array *fluid = mkarr(5 * nb, local_ext.volume), *u = mkarr(3 * nb, local_ext.volume);
  struct gkyl_array *xi = mkarr(3 * nb, local_ext.volume), *rho3 = mkarr(3 * nb, local_ext.volume);
  struct gkyl_array *grad_rho = mkarr(cdim * 3 * nb, local_ext.volume),
                    *gln = mkarr(cdim * nb, local_ext.volume);
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
    } else {
      gkyl_dg_mul_op_range(&cbasis, 0, prod, 0, fluid, i - 3, gln, &local);
    }
    gkyl_range_iter_init(&iter, &local);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&local, iter.idx);
      const double *pr = gkyl_array_cfetch(prod, loc);
      const double *num = i < 3 ?
                            (const double *)gkyl_array_cfetch(fluid, loc) + (1 + i) * nb :
                            (const double *)gkyl_array_cfetch(grad_rho, loc) + (i - 3) * 3 * nb;
      for (int k = 0; k < nb; ++k) {
        if (i < 3) {
          err_u = fmax(err_u, fabs(pr[k] - num[k]));
          max_m = fmax(max_m, fabs(num[k]));
        } else {
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

void
test_quotients_1x()
{
  test_quotients(1);
}
void
test_quotients_2x()
{
  test_quotients(2);
}
void
test_quotients_3x()
{
  test_quotients(3);
}

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
    gkyl_dg_calc_pkpm_vars_inew(&(struct gkyl_dg_calc_pkpm_vars_inp){
      .conf_grid = &grid,
      .conf_basis = &cbasis,
      .pos_map = pm,
      .use_gpu = false,
    });

  // qmem = q/m*[Ex,Ey,Ez,Bx,By,Bz,phi,psi] (8 blocks).
  struct gkyl_array *qmem = mkarr(8 * nb, local_ext.volume);
  struct gkyl_array *euler = mkarr(5 * nb, local_ext.volume);
  struct gkyl_array *rhs = mkarr(5 * nb, local_ext.volume);
  struct gkyl_array *omega_p_sq = mkarr(1, local_ext.volume);
  struct gkyl_array *cflrate = mkarr(1, local_ext.volume);

  // Choose B along z only, E along x only, so source = rho*E + rhou x B.
  double rho = 3.0, energy = 10.0;
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
  set_const_block(euler, 4, nb, energy);

  // Plasma frequency below the rate of the force; the CFL rate starts below both.
  double omega_p = 2.0;
  gkyl_array_clear(omega_p_sq, omega_p * omega_p);
  gkyl_array_clear(cflrate, 1.0);
  gkyl_array_clear(rhs, 0.0);
  gkyl_dg_calc_pkpm_vars_source(up, &local, qmem, omega_p_sq, euler, cflrate, rhs);

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

  // Rate of the source: the larger of omega_p + |q/m B| and sqrt(2 |f|^2/(rho e)), and the CFL
  // rate is raised to it only where it is larger.
  double e_int = energy - 0.5 * (rhoux * rhoux + rhouy * rhouy + rhouz * rhouz) / rho;
  double rate_accel = sqrt(2.0 * (sx * sx + sy * sy + sz * sz) / (rho * e_int));
  double rate_freq = omega_p + sqrt(Bx * Bx + By * By + Bz * Bz);
  TEST_CHECK(rate_accel > rate_freq);
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *c = gkyl_array_cfetch(cflrate, loc);
    TEST_CHECK(gkyl_compare(c[0], rate_accel, 1e-12));
  }
  omega_p = 7.0;
  gkyl_array_clear(omega_p_sq, omega_p * omega_p);
  gkyl_array_clear(cflrate, 1.0);
  gkyl_dg_calc_pkpm_vars_source(up, &local, qmem, omega_p_sq, euler, cflrate, rhs);
  rate_freq = omega_p + sqrt(Bx * Bx + By * By + Bz * Bz);
  TEST_CHECK(rate_freq > rate_accel);
  gkyl_array_clear(cflrate, 100.0);
  gkyl_range_iter_init(&iter, &local);
  gkyl_range_iter_next(&iter);
  double *c_first = gkyl_array_fetch(cflrate, gkyl_range_idx(&local, iter.idx));
  c_first[0] = 1.0;
  gkyl_dg_calc_pkpm_vars_source(up, &local, qmem, omega_p_sq, euler, cflrate, rhs);
  gkyl_range_iter_init(&iter, &local);
  bool is_first = true;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *c = gkyl_array_cfetch(cflrate, loc);
    TEST_CHECK(gkyl_compare(c[0], is_first ? rate_freq : 100.0, 1e-12));
    is_first = false;
  }

  gkyl_array_release(cflrate);
  gkyl_array_release(omega_p_sq);
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
  int cdim = b->ndim;
  int vert_lower[3] = {0, 0, 0}, vert_upper[3] = {1, 1, 1};
  struct gkyl_range vert_range;
  gkyl_range_init(&vert_range, cdim, vert_lower, vert_upper);
  struct gkyl_range_iter vert_iter;
  gkyl_range_iter_init(&vert_iter, &vert_range);
  while (gkyl_range_iter_next(&vert_iter)) {
    double z[3];
    for (int d = 0; d < cdim; ++d) {
      z[d] = vert_iter.idx[d] == 1 ? 1.0 : -1.0;
    }
    v[gkyl_range_idx(&vert_range, vert_iter.idx)] = b->eval_expand(z, c);
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
  struct gkyl_dg_calc_pkpm_vars *vars =
    gkyl_dg_calc_pkpm_vars_inew(&(struct gkyl_dg_calc_pkpm_vars_inp){
      .conf_grid = &grid,
      .conf_basis = &cbasis,
      .pos_map = pm,
      .use_gpu = false,
    });

  struct gkyl_array *fluid = mkarr(5 * nb, local_ext.volume), *count = mkarr(3, local_ext.volume);
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
  int num_rho = -1, num_e = -1, num_avg = -1;
  gkyl_dg_calc_pkpm_vars_positivity(vars, &local, eps, fluid, count, &num_rho, &num_e, &num_avg);
  TEST_CHECK(num_rho == 1);
  TEST_CHECK(num_e == 2); // the steep-density cell also needs the energy scaling
  TEST_CHECK(num_avg == 1);

  for (int i = 1; i <= 4; ++i) {
    idx[0] = i;
    long loc = gkyl_range_idx(&local_ext, idx);
    const double *f = gkyl_array_cfetch(fluid, loc), *cnt = gkyl_array_cfetch(count, loc);
    // averages preserved
    for (int c = 0; c < 5; ++c) {
      TEST_CHECK(f[c * nb] == avg_before[i - 1][c]);
    }
    if (i == 4) {
      TEST_CHECK(cnt[2] == 1.0 && cnt[0] == 0.0 && cnt[1] == 0.0);
      continue;
    }
    TEST_CHECK(cnt[2] == 0.0);
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
      double e_v =
        E_v[k] -
        0.5 * (M_v[0][k] * M_v[0][k] + M_v[1][k] * M_v[1][k] + M_v[2][k] * M_v[2][k]) / rho_v[k];
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

void
test_positivity_1x()
{
  test_positivity(1);
}
void
test_positivity_2x()
{
  test_positivity(2);
}
void
test_positivity_3x()
{
  test_positivity(3);
}

// Coefficients of an arbitrary p=1 expansion in every block of every cell.
static void
fill_expansions(struct gkyl_array *arr, double shift)
{
  for (long loc = 0; loc < arr->size; ++loc) {
    double *d = gkyl_array_fetch(arr, loc);
    for (int k = 0; k < arr->ncomp; ++k) {
      d[k] = 0.3 + 0.5 * sin(shift + 1.7 * k + 0.37 * loc);
    }
  }
}

// The diagnostic arrays are the p=1 projections of the pressure tensor and heat flux of the fluid
// equations, P_ij = (p_par - p_perp) b_i b_j + p_perp delta_ij + s_i b_j + b_i s_j and
// q_i = b_i (q_par + 2 q_perp)/2 + r_i, and of the variables [u, xi^2 theta, xi alpha_0, div b,
// bb:grad u, xi^2 theta div b, xi, c]: the updater against a Gauss-Legendre projection of the same
// expressions (three nodes per direction integrate a product of four p=1 factors exactly).
static void
test_diag(int cdim)
{
  double lower[3] = {-1.0, -1.0, -1.0}, upper[3] = {1.0, 1.0, 1.0};
  int cells[3] = {3, 2, 2}, ghost[3] = {0, 0, 0};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  int nb = cbasis.num_basis;
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct gkyl_vlasov_position_map *pm = identity_pos_map(&grid, &local, &local_ext, &cbasis);
  struct gkyl_dg_calc_pkpm_vars *up =
    gkyl_dg_calc_pkpm_vars_inew(&(struct gkyl_dg_calc_pkpm_vars_inp){
      .conf_grid = &grid,
      .conf_basis = &cbasis,
      .pos_map = pm,
      .use_gpu = false,
    });

  long vol = local_ext.volume;
  struct gkyl_array *moms = mkarr(10 * nb, vol), *bvar = mkarr(3 * nb, vol);
  struct gkyl_array *grad_b = mkarr(cdim * 3 * nb, vol), *flow = mkarr(3 * nb, vol);
  struct gkyl_array *grad_u = mkarr(cdim * 3 * nb, vol), *theta = mkarr(nb, vol);
  struct gkyl_array *alpha_0 = mkarr(nb, vol), *xi = mkarr(3 * nb, vol), *c_h = mkarr(nb, vol);
  struct gkyl_array *prim_io = mkarr(9 * nb, vol), *vars_io = mkarr(10 * nb, vol);
  fill_expansions(moms, 0.1);
  fill_expansions(bvar, 0.2);
  fill_expansions(grad_b, 0.3);
  fill_expansions(flow, 0.4);
  fill_expansions(grad_u, 0.5);
  fill_expansions(theta, 0.6);
  fill_expansions(alpha_0, 0.7);
  fill_expansions(xi, 0.8);
  fill_expansions(c_h, 0.9);

  gkyl_dg_calc_pkpm_vars_diag(
    up, &local, moms, bvar, grad_b, flow, grad_u, theta, alpha_0, xi, c_h, prim_io, vars_io
  );

  static const int ij[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
  const double *ordinates = gkyl_gauss_ordinates[3], *weights = gkyl_gauss_weights[3];
  int num_quad = cdim == 1 ? 3 : (cdim == 2 ? 9 : 27);
  double err_prim = 0.0, err_vars = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *moms_c = gkyl_array_cfetch(moms, loc), *bvar_c = gkyl_array_cfetch(bvar, loc);
    const double *grad_b_c = gkyl_array_cfetch(grad_b, loc), *flow_c = gkyl_array_cfetch(flow, loc);
    const double *grad_u_c = gkyl_array_cfetch(grad_u, loc),
                 *theta_c = gkyl_array_cfetch(theta, loc);
    const double *alpha_c = gkyl_array_cfetch(alpha_0, loc), *xi_c = gkyl_array_cfetch(xi, loc);
    const double *c_h_c = gkyl_array_cfetch(c_h, loc);
    double ref_prim[9 * 8] = {0.0}, ref_vars[10 * 8] = {0.0};
    for (int quad = 0; quad < num_quad; ++quad) {
      double node[3], wgt = 1.0, basis_at_node[8];
      for (int d = 0, rem = quad; d < cdim; ++d, rem /= 3) {
        node[d] = ordinates[rem % 3];
        wgt *= weights[rem % 3];
      }
      cbasis.eval(node, basis_at_node);
      double mq[10], bq[3], uq[3], div_b = 0.0, bb_grad_u = 0.0;
      for (int c = 0; c < 10; ++c) {
        mq[c] = cbasis.eval_expand(node, moms_c + c * nb);
      }
      for (int i = 0; i < 3; ++i) {
        bq[i] = cbasis.eval_expand(node, bvar_c + i * nb);
        uq[i] = cbasis.eval_expand(node, flow_c + i * nb);
      }
      for (int d = 0; d < cdim; ++d) {
        div_b += cbasis.eval_expand(node, grad_b_c + (d * 3 + d) * nb);
        for (int j = 0; j < 3; ++j) {
          bb_grad_u += bq[d] * bq[j] * cbasis.eval_expand(node, grad_u_c + (d * 3 + j) * nb);
        }
      }
      double ppar = mq[0], pperp = mq[1], qpar = mq[2], qperp = mq[3];
      const double *sq = mq + 4, *rq = mq + 7;
      double prim[9], vars[10];
      for (int c = 0; c < 6; ++c) {
        int i = ij[c][0], j = ij[c][1];
        prim[c] =
          (ppar - pperp) * bq[i] * bq[j] + (i == j ? pperp : 0.0) + sq[i] * bq[j] + bq[i] * sq[j];
      }
      for (int i = 0; i < 3; ++i) {
        prim[6 + i] = bq[i] * (0.5 * qpar + qperp) + rq[i];
        vars[i] = uq[i];
      }
      double xi_1 = cbasis.eval_expand(node, xi_c), xi_sq = cbasis.eval_expand(node, xi_c + nb);
      double tperp = xi_sq * cbasis.eval_expand(node, theta_c);
      vars[3] = tperp;
      vars[4] = xi_1 * cbasis.eval_expand(node, alpha_c);
      vars[5] = div_b;
      vars[6] = bb_grad_u;
      vars[7] = tperp * div_b;
      vars[8] = xi_1;
      vars[9] = cbasis.eval_expand(node, c_h_c);
      for (int k = 0; k < nb; ++k) {
        for (int c = 0; c < 9; ++c) {
          ref_prim[c * nb + k] += wgt * basis_at_node[k] * prim[c];
        }
        for (int c = 0; c < 10; ++c) {
          ref_vars[c * nb + k] += wgt * basis_at_node[k] * vars[c];
        }
      }
    }
    const double *out_prim = gkyl_array_cfetch(prim_io, loc);
    const double *out_vars = gkyl_array_cfetch(vars_io, loc);
    for (int k = 0; k < 9 * nb; ++k) {
      err_prim = fmax(err_prim, fabs(out_prim[k] - ref_prim[k]));
    }
    for (int k = 0; k < 10 * nb; ++k) {
      err_vars = fmax(err_vars, fabs(out_vars[k] - ref_vars[k]));
    }
  }
  TEST_CHECK(err_prim < 1e-13);
  TEST_MSG("cdim=%d: pressure tensor and heat flux against the projection: %g", cdim, err_prim);
  TEST_CHECK(err_vars < 1e-13);
  TEST_MSG("cdim=%d: pkpm variables against the projection: %g", cdim, err_vars);

  gkyl_array_release(vars_io);
  gkyl_array_release(prim_io);
  gkyl_array_release(c_h);
  gkyl_array_release(xi);
  gkyl_array_release(alpha_0);
  gkyl_array_release(theta);
  gkyl_array_release(grad_u);
  gkyl_array_release(flow);
  gkyl_array_release(grad_b);
  gkyl_array_release(bvar);
  gkyl_array_release(moms);
  gkyl_dg_calc_pkpm_vars_release(up);
  gkyl_vlasov_position_map_release(pm);
}

void
test_diag_1x()
{
  test_diag(1);
}
void
test_diag_2x()
{
  test_diag(2);
}
void
test_diag_3x()
{
  test_diag(3);
}

TEST_LIST = {
  {"test_positivity_1x", test_positivity_1x},
  {"test_positivity_2x", test_positivity_2x},
  {"test_positivity_3x", test_positivity_3x},
  {"test_quotients_1x", test_quotients_1x},
  {"test_quotients_2x", test_quotients_2x},
  {"test_quotients_3x", test_quotients_3x},
  {"source_1x_p1", test_source_1x_p1},
  {"test_diag_1x", test_diag_1x},
  {"test_diag_2x", test_diag_2x},
  {"test_diag_3x", test_diag_3x},
  {NULL, NULL}
};
