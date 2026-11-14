// Tests of the nodal locally implicit fluid-Maxwell source update of the p=1 PKPM solver
// (gkyl_dg_calc_pkpm_em_coupling, Wang et al. JCP 2020 in the unknowns [Mbar, eps0 Ebar]):
//   * static field, B = 0: a pure electric push M^{n+1} = M^n + dt (q/m) rho E;
//   * plasma oscillation in a uniform cell: the quadratic invariant |M|^2/(2 rho) + eps0 |E|^2/2
//     is conserved to round-off (and the fluid internal energy is untouched: the fluid energy
//     gains exactly the change of |M|^2/(2 rho)), the numerical frequency is the implicit midpoint one,
//     (2/dt) atan(omega_p dt/2), and its error against omega_p is second order in dt;
//   * cyclotron rotation in a static field: |M| is conserved, the rotation per step is
//     2 atan(Omega dt/2);
//   * varying rho, E, B (1x, 2x, 3x, two species): the modal momentum increment equals the
//     exact weak projection of dt (q/m)(rho Ebar + Mbar x B) with Ebar, Mbar the modal
//     half-step fields (products of two p=1 fields are integrated exactly by the 2-point
//     rule: no aliasing), Ampere's law holds mode by mode, and the kinetic plus field energy
//     is conserved node by node.
#include <math.h>
#include <string.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_pkpm_em_coupling.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

static struct gkyl_array *
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// Set the cell-average value of block `blk` (nb basis components) to `val` in every cell.
static void
set_const_block(struct gkyl_array *arr, int blk, int nb, double val, int cdim)
{
  gkyl_array_shiftc(arr, val * pow(sqrt(2.0), cdim), blk * nb);
}

struct setup {
  int cdim, nc;
  struct gkyl_rect_grid grid;
  struct gkyl_basis basis;
  struct gkyl_range local, local_ext;
  struct gkyl_array *app_accel[2], *ext_em, *app_current, *fluid[2], *em;
};

static void
setup_init(struct setup *s, int cdim, int ns)
{
  double lower[3] = {-1.0, -1.0, -1.0}, upper[3] = {1.0, 1.0, 1.0};
  int cells[3] = {5, 3, 2}, ghost[3] = {1, 1, 1};
  s->cdim = cdim;
  gkyl_rect_grid_init(&s->grid, cdim, lower, upper, cells);
  gkyl_cart_modal_serendip(&s->basis, cdim, 1);
  gkyl_create_grid_ranges(&s->grid, ghost, &s->local_ext, &s->local);
  s->nc = s->basis.num_basis;
  long vol = s->local_ext.volume;
  for (int i = 0; i < ns; ++i) {
    s->app_accel[i] = mkarr(3 * s->nc, vol);
    s->fluid[i] = mkarr(5 * s->nc, vol); // [rho, rho ux, rho uy, rho uz, E]
    gkyl_array_clear(s->app_accel[i], 0.0);
    gkyl_array_clear(s->fluid[i], 0.0);
  }
  s->ext_em = mkarr(8 * s->nc, vol);
  s->app_current = mkarr(3 * s->nc, vol);
  s->em = mkarr(8 * s->nc, vol);
  gkyl_array_clear(s->ext_em, 0.0);
  gkyl_array_clear(s->app_current, 0.0);
  gkyl_array_clear(s->em, 0.0);
}

static void
setup_release(struct setup *s, int ns)
{
  for (int i = 0; i < ns; ++i) {
    gkyl_array_release(s->app_accel[i]);
    gkyl_array_release(s->fluid[i]);
  }
  gkyl_array_release(s->ext_em);
  gkyl_array_release(s->app_current);
  gkyl_array_release(s->em);
}

static double
avg(const struct gkyl_array *arr, long loc, int blk, int nb, int cdim)
{
  const double *d = gkyl_array_cfetch(arr, loc);
  return d[blk * nb] / pow(sqrt(2.0), cdim);
}

// Static field, B = 0: M^{n+1} = M^n + dt (q/m) rho E, E unchanged.
static void
test_static_E_push(int cdim)
{
  struct setup s;
  setup_init(&s, cdim, 1);
  int nb = s.nc;
  double qbym[1] = {-2.0}, eps0 = 1.0, dt = 0.37;
  double rho = 3.0, Ex = 2.0, Ey = -0.5, Ez = 1.0, Mx = 1.0, My = 4.0, Mz = -2.0, E0 = 7.0;
  set_const_block(s.fluid[0], 0, nb, rho, cdim);
  set_const_block(s.fluid[0], 1, nb, Mx, cdim);
  set_const_block(s.fluid[0], 2, nb, My, cdim);
  set_const_block(s.fluid[0], 3, nb, Mz, cdim);
  set_const_block(s.fluid[0], 4, nb, E0, cdim);
  set_const_block(s.em, 0, nb, Ex, cdim);
  set_const_block(s.em, 1, nb, Ey, cdim);
  set_const_block(s.em, 2, nb, Ez, cdim);

  struct gkyl_dg_calc_pkpm_em_coupling *up =
    gkyl_dg_calc_pkpm_em_coupling_new(&s.basis, &s.local, 1, qbym, eps0, true, false);
  const struct gkyl_array *aa[1] = {s.app_accel[0]};
  struct gkyl_array *ff[1] = {s.fluid[0]};
  gkyl_dg_calc_pkpm_em_coupling_advance(up, dt, aa, s.ext_em, s.app_current, ff, s.em);

  double Mxn = Mx + dt * qbym[0] * rho * Ex, Myn = My + dt * qbym[0] * rho * Ey, Mzn = Mz + dt * qbym[0] * rho * Ez;
  double En = E0 + 0.5 * (Mxn * Mxn + Myn * Myn + Mzn * Mzn - Mx * Mx - My * My - Mz * Mz) / rho;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&s.local_ext, iter.idx);
    TEST_CHECK(gkyl_compare(avg(s.fluid[0], loc, 0, nb, cdim), rho, 1e-14));
    TEST_CHECK(gkyl_compare(avg(s.fluid[0], loc, 1, nb, cdim), Mxn, 1e-12));
    TEST_CHECK(gkyl_compare(avg(s.fluid[0], loc, 2, nb, cdim), Myn, 1e-12));
    TEST_CHECK(gkyl_compare(avg(s.fluid[0], loc, 3, nb, cdim), Mzn, 1e-12));
    TEST_CHECK(gkyl_compare(avg(s.fluid[0], loc, 4, nb, cdim), En, 1e-12));
    TEST_CHECK(gkyl_compare(avg(s.em, loc, 0, nb, cdim), Ex, 1e-14));
    TEST_CHECK(gkyl_compare(avg(s.em, loc, 1, nb, cdim), Ey, 1e-14));
    for (int b = 0; b < 5; ++b) {
      for (int k = 1; k < nb; ++k) {
        const double *f = gkyl_array_cfetch(s.fluid[0], loc);
        TEST_CHECK(fabs(f[b * nb + k]) < 1e-13);
      }
    }
  }
  gkyl_dg_calc_pkpm_em_coupling_release(up);
  setup_release(&s, 1);
}

// Plasma oscillation in a uniform cell: returns the phase error of E_x after t_end against
// omega_p t_end; checks the invariant and the implicit-midpoint frequency to round-off.
static double
plasma_oscillation(double dt, double t_end)
{
  struct setup s;
  setup_init(&s, 1, 1);
  int nb = s.nc;
  double qbym[1] = {-1.0}, eps0 = 2.0, rho = 3.0, E0 = 0.7, Eint = 4.5;
  double wp = sqrt(qbym[0] * qbym[0] * rho / eps0);
  set_const_block(s.fluid[0], 0, nb, rho, 1);
  set_const_block(s.fluid[0], 4, nb, Eint, 1); // fluid energy: internal only (M = 0)
  set_const_block(s.em, 0, nb, E0, 1);
  struct gkyl_dg_calc_pkpm_em_coupling *up =
    gkyl_dg_calc_pkpm_em_coupling_new(&s.basis, &s.local, 1, qbym, eps0, false, false);
  const struct gkyl_array *aa[1] = {s.app_accel[0]};
  struct gkyl_array *ff[1] = {s.fluid[0]};
  long loc = gkyl_range_idx(&s.local_ext, s.local.lower);
  double inv0 = 0.5 * eps0 * E0 * E0, wnum = 2.0 / dt * atan(0.5 * wp * dt);
  int nsteps = (int)round(t_end / dt);
  for (int n = 1; n <= nsteps; ++n) {
    gkyl_dg_calc_pkpm_em_coupling_advance(up, dt, aa, s.ext_em, s.app_current, ff, s.em);
    double Mx = avg(s.fluid[0], loc, 1, nb, 1), Ex = avg(s.em, loc, 0, nb, 1), Ef = avg(s.fluid[0], loc, 4, nb, 1);
    double inv = 0.5 * Mx * Mx / rho + 0.5 * eps0 * Ex * Ex;
    TEST_CHECK(gkyl_compare(inv, inv0, 1e-13));
    TEST_MSG("dt=%g step %d: invariant %g vs %g", dt, n, inv, inv0);
    // the fluid energy follows the bulk kinetic energy: the internal energy is untouched
    TEST_CHECK(gkyl_compare(Ef - 0.5 * Mx * Mx / rho, Eint, 1e-13));
    TEST_MSG("dt=%g step %d: internal energy %g vs %g", dt, n, Ef - 0.5 * Mx * Mx / rho, Eint);
    // exact solution of the implicit midpoint rule
    TEST_CHECK(gkyl_compare(Ex, E0 * cos(wnum * n * dt), 1e-12));
    TEST_MSG("dt=%g step %d: Ex %g vs midpoint %g", dt, n, Ex, E0 * cos(wnum * n * dt));
  }
  double Ex = avg(s.em, loc, 0, nb, 1);
  double err = fabs(Ex - E0 * cos(wp * t_end));
  gkyl_dg_calc_pkpm_em_coupling_release(up);
  setup_release(&s, 1);
  return err;
}

static void
test_plasma_oscillation()
{
  double e1 = plasma_oscillation(0.2, 4.0), e2 = plasma_oscillation(0.1, 4.0), e3 = plasma_oscillation(0.05, 4.0);
  TEST_CHECK(e1 / e2 > 3.8 && e1 / e2 < 4.2);
  TEST_MSG("error ratio dt 0.2/0.1: %g (errors %g %g)", e1 / e2, e1, e2);
  TEST_CHECK(e2 / e3 > 3.9 && e2 / e3 < 4.1);
  TEST_MSG("error ratio dt 0.1/0.05: %g (errors %g %g)", e2 / e3, e2, e3);
}

// Cyclotron rotation in a static field: |M| conserved, rotation 2 atan(Omega dt/2) per step.
static void
test_cyclotron()
{
  struct setup s;
  setup_init(&s, 2, 1);
  int nb = s.nc, cdim = 2;
  double qbym[1] = {1.5}, eps0 = 1.0, rho = 2.0, Bz = 0.8, dt = 0.3, M0 = 1.3;
  set_const_block(s.fluid[0], 0, nb, rho, cdim);
  set_const_block(s.fluid[0], 1, nb, M0, cdim);
  set_const_block(s.ext_em, 5, nb, Bz, cdim); // external B_z
  struct gkyl_dg_calc_pkpm_em_coupling *up =
    gkyl_dg_calc_pkpm_em_coupling_new(&s.basis, &s.local, 1, qbym, eps0, true, false);
  const struct gkyl_array *aa[1] = {s.app_accel[0]};
  struct gkyl_array *ff[1] = {s.fluid[0]};
  long loc = gkyl_range_idx(&s.local_ext, s.local.lower);
  double omega = qbym[0] * Bz, dphi = 2.0 * atan(0.5 * omega * dt);
  for (int n = 1; n <= 20; ++n) {
    gkyl_dg_calc_pkpm_em_coupling_advance(up, dt, aa, s.ext_em, s.app_current, ff, s.em);
    double Mx = avg(s.fluid[0], loc, 1, nb, cdim), My = avg(s.fluid[0], loc, 2, nb, cdim);
    TEST_CHECK(gkyl_compare(sqrt(Mx * Mx + My * My), M0, 1e-13));
    TEST_CHECK(fabs(avg(s.fluid[0], loc, 4, nb, cdim)) < 1e-13); // no work: |M| is conserved
    // dM/dt = qbym M x B: (Mx, My) rotates clockwise for omega > 0
    TEST_CHECK(gkyl_compare(Mx, M0 * cos(n * dphi), 1e-12));
    TEST_CHECK(gkyl_compare(My, -M0 * sin(n * dphi), 1e-12));
    TEST_MSG("step %d: (%g, %g) vs (%g, %g)", n, Mx, My, M0 * cos(n * dphi), -M0 * sin(n * dphi));
  }
  gkyl_dg_calc_pkpm_em_coupling_release(up);
  setup_release(&s, 1);
}

// Varying fields: two species, non-uniform rho, M, E, B (and external B) projected onto p=1.
static void
eval_fields(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  // [rho0, rho1, M0(3), M1(3), E(3), B(3), Bext(3), a0(3), Japp(3)] = 20 components
  fout[0] = 2.0 + 0.5 * sin(x) * cos(y) + 0.2 * z;
  fout[1] = 1.0 + 0.3 * cos(x) + 0.1 * sin(y) * z;
  fout[2] = 0.4 * sin(y) + 0.1 * x;
  fout[3] = -0.3 + 0.2 * cos(x) * sin(z);
  fout[4] = 0.5 * x * y;
  fout[5] = 0.2 * cos(y) - 0.1 * z;
  fout[6] = 0.6 * sin(x);
  fout[7] = -0.2 + 0.3 * x * z;
  fout[8] = 0.8 * cos(x) + 0.2 * y;
  fout[9] = -0.5 * sin(y) * z;
  fout[10] = 0.3 + 0.4 * x;
  fout[11] = 1.0 + 0.4 * sin(x) * cos(y);
  fout[12] = 0.6 * cos(x) * sin(z) + 0.5 * sin(y);
  fout[13] = -0.5 + 0.3 * sin(y) * cos(z) + 0.2 * x;
  fout[14] = 0.2 * y;
  fout[15] = -0.3 * x;
  fout[16] = 0.7;
  fout[17] = 0.1 * sin(x);
  fout[18] = -0.2 * cos(y);
  fout[19] = 0.15 * z;
  fout[20] = 0.05 * x * y;
  fout[21] = -0.1;
  fout[22] = 0.2 * sin(z);
}

static const double gl4_x[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
  0.8611363115940526};
static const double gl4_w[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
  0.3478548451374538};

static void
test_varying(int cdim)
{
  struct setup s;
  setup_init(&s, cdim, 2);
  int nb = s.nc;
  double qbym[2] = {-1.0, 0.5}, eps0 = 1.3, dt = 0.21;
  struct gkyl_array *all = mkarr(23 * nb, s.local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s.grid, &s.basis, 8, 23, eval_fields, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s.local_ext, all);
  gkyl_proj_on_basis_release(proj);
  struct gkyl_array *tmp = mkarr(nb, s.local_ext.volume);
  // scatter the blocks
  for (int c = 0; c < 23; ++c) {
    gkyl_array_set_offset(tmp, 1.0, all, c * nb);
    struct gkyl_array *dst = 0;
    int blk = 0;
    if (c < 2) { dst = s.fluid[c]; blk = 0; }
    else if (c < 5) { dst = s.fluid[0]; blk = c - 2 + 1; }
    else if (c < 8) { dst = s.fluid[1]; blk = c - 5 + 1; }
    else if (c < 11) { dst = s.em; blk = c - 8; }
    else if (c < 14) { dst = s.em; blk = c - 11 + 3; }
    else if (c < 17) { dst = s.ext_em; blk = c - 14 + 3; }
    else if (c < 20) { dst = s.app_accel[0]; blk = c - 17; }
    else { dst = s.app_current; blk = c - 20; }
    gkyl_array_set_offset(dst, 1.0, tmp, blk * nb);
  }
  // fluid energies: a varying internal energy on top of the bulk kinetic energy
  gkyl_array_set_offset(tmp, 1.0, all, 0);
  gkyl_array_set_offset(s.fluid[0], 2.5, tmp, 4 * nb);
  gkyl_array_set_offset(tmp, 1.0, all, 1 * nb);
  gkyl_array_set_offset(s.fluid[1], 1.5, tmp, 4 * nb);
  struct gkyl_array *M_old[2] = {mkarr(5 * nb, s.local_ext.volume), mkarr(5 * nb, s.local_ext.volume)};
  struct gkyl_array *em_old = mkarr(8 * nb, s.local_ext.volume);
  gkyl_array_copy(M_old[0], s.fluid[0]);
  gkyl_array_copy(M_old[1], s.fluid[1]);
  gkyl_array_copy(em_old, s.em);

  struct gkyl_dg_calc_pkpm_em_coupling *up =
    gkyl_dg_calc_pkpm_em_coupling_new(&s.basis, &s.local, 2, qbym, eps0, false, false);
  const struct gkyl_array *aa[2] = {s.app_accel[0], s.app_accel[1]};
  struct gkyl_array *ff[2] = {s.fluid[0], s.fluid[1]};
  gkyl_dg_calc_pkpm_em_coupling_advance(up, dt, aa, s.ext_em, s.app_current, ff, s.em);

  int nq = 1;
  for (int d = 0; d < cdim; ++d) {
    nq *= 4;
  }
  double max_err = 0.0, max_amp = 0.0, max_ferr = 0.0, max_eerr = 0.0, max_Eferr = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&s.local_ext, iter.idx);
    const double *em_n = gkyl_array_cfetch(s.em, loc), *em_o = gkyl_array_cfetch(em_old, loc);
    const double *ext = gkyl_array_cfetch(s.ext_em, loc), *Ja = gkyl_array_cfetch(s.app_current, loc);
    // Ampere's law mode by mode: eps0 (E_new - E_old) = -dt (sum_s qbym_s Mbar_s + J_app)
    for (int c = 0; c < 3; ++c) {
      for (int k = 0; k < nb; ++k) {
        double lhs = eps0 * (em_n[c * nb + k] - em_o[c * nb + k]);
        double rhs = -dt * Ja[c * nb + k];
        for (int sp = 0; sp < 2; ++sp) {
          const double *Mn = gkyl_array_cfetch(s.fluid[sp], loc), *Mo = gkyl_array_cfetch(M_old[sp], loc);
          rhs -= dt * qbym[sp] * 0.5 * (Mn[(1 + c) * nb + k] + Mo[(1 + c) * nb + k]);
        }
        max_ferr = fmax(max_ferr, fabs(lhs - rhs));
      }
    }
    // Momentum: modal increment vs exact weak projection of dt qbym (rho Ebar + Mbar x B)
    // + dt rho (qbym E_ext + a), by 4-point quadrature of the p=1 fields' products.
    for (int sp = 0; sp < 2; ++sp) {
      const double *Mn = gkyl_array_cfetch(s.fluid[sp], loc), *Mo = gkyl_array_cfetch(M_old[sp], loc);
      const double *rho = Mn, *acc = gkyl_array_cfetch(s.app_accel[sp], loc); // density: leading block, unchanged
      double ref[24] = {0.0};
      for (int n = 0; n < nq; ++n) {
        double z[3], w = 1.0;
        int m = n;
        for (int d = cdim - 1; d >= 0; --d) {
          z[d] = gl4_x[m % 4];
          w *= gl4_w[m % 4];
          m /= 4;
        }
        double rq = s.basis.eval_expand(z, rho), Eb[3], B[3], Mb[3], Ee[3], aq[3];
        for (int c = 0; c < 3; ++c) {
          Eb[c] = 0.5 * (s.basis.eval_expand(z, &em_n[c * nb]) + s.basis.eval_expand(z, &em_o[c * nb]));
          B[c] = s.basis.eval_expand(z, &em_o[(3 + c) * nb]) + s.basis.eval_expand(z, &ext[(3 + c) * nb]);
          Mb[c] = 0.5 * (s.basis.eval_expand(z, &Mn[(1 + c) * nb]) + s.basis.eval_expand(z, &Mo[(1 + c) * nb]));
          Ee[c] = s.basis.eval_expand(z, &ext[c * nb]);
          aq[c] = s.basis.eval_expand(z, &acc[c * nb]);
        }
        double F[3] = {Mb[1] * B[2] - Mb[2] * B[1], Mb[2] * B[0] - Mb[0] * B[2], Mb[0] * B[1] - Mb[1] * B[0]};
        for (int c = 0; c < 3; ++c) {
          double src = dt * (qbym[sp] * (rq * (Eb[c] + Ee[c]) + F[c]) + rq * aq[c]);
          for (int k = 0; k < nb; ++k) {
            double unit[8] = {0.0};
            unit[k] = 1.0;
            ref[c * nb + k] += w * s.basis.eval_expand(z, unit) * src;
          }
        }
      }
      for (int k = 0; k < 3 * nb; ++k) {
        max_err = fmax(max_err, fabs(Mn[nb + k] - Mo[nb + k] - ref[k]));
        max_amp = fmax(max_amp, fabs(ref[k]));
      }
      for (int k = 0; k < nb; ++k) {
        TEST_CHECK(Mn[k] == Mo[k]); // the density is untouched
      }
    }
    // Energy at the 2-point nodes: sum_s (|M_new|^2 - |M_old|^2)/(2 rho) + eps0 (|E_new|^2 - |E_old|^2)/2
    // = dt sum_s Mbar.(qbym E_ext + a) - dt Ebar.J_app   (work of the external terms only)
    int node_lower[3] = {0, 0, 0}, node_upper[3] = {1, 1, 1};
    struct gkyl_range node_range;
    gkyl_range_init(&node_range, cdim, node_lower, node_upper);
    struct gkyl_range_iter node_iter;
    gkyl_range_iter_init(&node_iter, &node_range);
    while (gkyl_range_iter_next(&node_iter)) {
      double z[3];
      for (int d = 0; d < cdim; ++d) {
        z[d] = gkyl_gauss_ordinates_2[node_iter.idx[d]];
      }
      double de = 0.0, work = 0.0;
      for (int c = 0; c < 3; ++c) {
        double En = s.basis.eval_expand(z, &em_n[c * nb]), Eo = s.basis.eval_expand(z, &em_o[c * nb]);
        de += 0.5 * eps0 * (En * En - Eo * Eo);
        work -= dt * 0.5 * (En + Eo) * s.basis.eval_expand(z, &Ja[c * nb]);
      }
      for (int sp = 0; sp < 2; ++sp) {
        const double *Mn = gkyl_array_cfetch(s.fluid[sp], loc), *Mo = gkyl_array_cfetch(M_old[sp], loc);
        const double *acc = gkyl_array_cfetch(s.app_accel[sp], loc);
        double rq = s.basis.eval_expand(z, Mn), dke = 0.0;
        for (int c = 0; c < 3; ++c) {
          double a = s.basis.eval_expand(z, &Mn[(1 + c) * nb]), b = s.basis.eval_expand(z, &Mo[(1 + c) * nb]);
          dke += 0.5 * (a * a - b * b) / rq;
          work += dt * 0.5 * (a + b) * (qbym[sp] * s.basis.eval_expand(z, &ext[c * nb]) + s.basis.eval_expand(z, &acc[c * nb]));
        }
        de += dke;
        // the fluid energy gains exactly the nodal change of the bulk kinetic energy
        double dEf = s.basis.eval_expand(z, &Mn[4 * nb]) - s.basis.eval_expand(z, &Mo[4 * nb]);
        max_Eferr = fmax(max_Eferr, fabs(dEf - dke));
      }
      max_eerr = fmax(max_eerr, fabs(de - work));
    }
  }
  TEST_CHECK(max_amp > 0.01);
  TEST_CHECK(max_err < 1e-12 * max_amp);
  TEST_MSG("cdim=%d: momentum increment vs exact weak projection %g (amplitude %g)", cdim, max_err, max_amp);
  TEST_CHECK(max_ferr < 1e-12);
  TEST_MSG("cdim=%d: Ampere residual %g", cdim, max_ferr);
  TEST_CHECK(max_eerr < 1e-12);
  TEST_MSG("cdim=%d: nodal energy residual %g", cdim, max_eerr);
  TEST_CHECK(max_Eferr < 1e-12);
  TEST_MSG("cdim=%d: fluid energy vs bulk kinetic energy change %g", cdim, max_Eferr);

  gkyl_array_release(em_old);
  gkyl_array_release(M_old[0]);
  gkyl_array_release(M_old[1]);
  gkyl_array_release(tmp);
  gkyl_array_release(all);
  gkyl_dg_calc_pkpm_em_coupling_release(up);
  setup_release(&s, 2);
}

void test_static_E_push_1x() { test_static_E_push(1); }
void test_static_E_push_2x() { test_static_E_push(2); }
void test_varying_1x() { test_varying(1); }
void test_varying_2x() { test_varying(2); }
void test_varying_3x() { test_varying(3); }

TEST_LIST = {
  {"static_E_push_1x", test_static_E_push_1x},
  {"static_E_push_2x", test_static_E_push_2x},
  {"plasma_oscillation", test_plasma_oscillation},
  {"cyclotron", test_cyclotron},
  {"varying_1x", test_varying_1x},
  {"varying_2x", test_varying_2x},
  {"varying_3x", test_varying_3x},
  {NULL, NULL},
};
