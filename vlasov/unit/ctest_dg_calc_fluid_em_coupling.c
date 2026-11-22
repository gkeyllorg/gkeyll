// Tests of the nodal locally implicit fluid-Maxwell source update of p=1 DG fluids
// (gkyl_dg_calc_fluid_em_coupling, Wang et al. JCP 2020 in the unknowns [Mbar, eps0 Ebar]):
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
//     is conserved node by node;
//   * resistive layer: the coefficients g2, phi, phi2 are correctly rounded on both sides of the
//     series threshold and exact at sigma dt = 0; with a neutral species E alone is damped by
//     exp(-sigma dt) at each node, also for sigma dt >> 1 (no sign change), and E = -J_app/(eps0 sigma)
//     is steady; with the varying state (two species, 1x, 2x, 3x) the nodal energy changes by the
//     external work minus eps0 sigma int |E|^2 dt; sigma = 0 gives bitwise the result without a
//     layer, and sigma dt = 1e-13 agrees with it to round-off.
#include <math.h>
#include <string.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_fluid_em_coupling.h>
#include <gkyl_fluid_em_coupling_kernels.h>
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

  struct gkyl_dg_calc_fluid_em_coupling *up =
    gkyl_dg_calc_fluid_em_coupling_inew(&(struct gkyl_dg_calc_fluid_em_coupling_inp){
      .conf_basis = &s.basis,
      .mem_range = &s.local,
      .num_species = 1,
      .qbym = qbym,
      .epsilon0 = eps0,
      .is_field_static = true,
      .use_gpu = false,
    });
  const struct gkyl_array *aa[1] = {s.app_accel[0]};
  struct gkyl_array *ff[1] = {s.fluid[0]};
  gkyl_dg_calc_fluid_em_coupling_advance(up, dt, aa, s.ext_em, s.app_current, NULL, ff, s.em);

  double Mxn = Mx + dt * qbym[0] * rho * Ex, Myn = My + dt * qbym[0] * rho * Ey,
         Mzn = Mz + dt * qbym[0] * rho * Ez;
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
  gkyl_dg_calc_fluid_em_coupling_release(up);
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
  struct gkyl_dg_calc_fluid_em_coupling *up =
    gkyl_dg_calc_fluid_em_coupling_inew(&(struct gkyl_dg_calc_fluid_em_coupling_inp){
      .conf_basis = &s.basis,
      .mem_range = &s.local,
      .num_species = 1,
      .qbym = qbym,
      .epsilon0 = eps0,
      .is_field_static = false,
      .use_gpu = false,
    });
  const struct gkyl_array *aa[1] = {s.app_accel[0]};
  struct gkyl_array *ff[1] = {s.fluid[0]};
  long loc = gkyl_range_idx(&s.local_ext, s.local.lower);
  double inv0 = 0.5 * eps0 * E0 * E0, wnum = 2.0 / dt * atan(0.5 * wp * dt);
  int nsteps = (int)round(t_end / dt);
  for (int n = 1; n <= nsteps; ++n) {
    gkyl_dg_calc_fluid_em_coupling_advance(up, dt, aa, s.ext_em, s.app_current, NULL, ff, s.em);
    double Mx = avg(s.fluid[0], loc, 1, nb, 1), Ex = avg(s.em, loc, 0, nb, 1),
           Ef = avg(s.fluid[0], loc, 4, nb, 1);
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
  gkyl_dg_calc_fluid_em_coupling_release(up);
  setup_release(&s, 1);
  return err;
}

static void
test_plasma_oscillation()
{
  double e1 = plasma_oscillation(0.2, 4.0), e2 = plasma_oscillation(0.1, 4.0),
         e3 = plasma_oscillation(0.05, 4.0);
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
  struct gkyl_dg_calc_fluid_em_coupling *up =
    gkyl_dg_calc_fluid_em_coupling_inew(&(struct gkyl_dg_calc_fluid_em_coupling_inp){
      .conf_basis = &s.basis,
      .mem_range = &s.local,
      .num_species = 1,
      .qbym = qbym,
      .epsilon0 = eps0,
      .is_field_static = true,
      .use_gpu = false,
    });
  const struct gkyl_array *aa[1] = {s.app_accel[0]};
  struct gkyl_array *ff[1] = {s.fluid[0]};
  long loc = gkyl_range_idx(&s.local_ext, s.local.lower);
  double omega = qbym[0] * Bz, dphi = 2.0 * atan(0.5 * omega * dt);
  for (int n = 1; n <= 20; ++n) {
    gkyl_dg_calc_fluid_em_coupling_advance(up, dt, aa, s.ext_em, s.app_current, NULL, ff, s.em);
    double Mx = avg(s.fluid[0], loc, 1, nb, cdim), My = avg(s.fluid[0], loc, 2, nb, cdim);
    TEST_CHECK(gkyl_compare(sqrt(Mx * Mx + My * My), M0, 1e-13));
    TEST_CHECK(fabs(avg(s.fluid[0], loc, 4, nb, cdim)) < 1e-13); // no work: |M| is conserved
    // dM/dt = qbym M x B: (Mx, My) rotates clockwise for omega > 0
    TEST_CHECK(gkyl_compare(Mx, M0 * cos(n * dphi), 1e-12));
    TEST_CHECK(gkyl_compare(My, -M0 * sin(n * dphi), 1e-12));
    TEST_MSG("step %d: (%g, %g) vs (%g, %g)", n, Mx, My, M0 * cos(n * dphi), -M0 * sin(n * dphi));
  }
  gkyl_dg_calc_fluid_em_coupling_release(up);
  setup_release(&s, 1);
}

// Varying fields: two species, non-uniform rho, M, E, B (and external B) projected onto p=1.
static void
eval_fields(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  // [rho0, rho1, M0(3), M1(3), E(3), B(3), Bext(3), a0(3), Japp(3)] = 23 components
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

static const double gl4_x[4] = {
  -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526
};
static const double gl4_w[4] = {
  0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538
};

// Spatially varying state of two species and the fields: projections of eval_fields, with a
// varying internal energy on top of the bulk kinetic energy.
static void
set_varying_state(struct setup *sp, int cdim)
{
  struct setup s = *sp;
  int nb = s.nc;
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
    if (c < 2) {
      dst = s.fluid[c];
      blk = 0;
    } else if (c < 5) {
      dst = s.fluid[0];
      blk = c - 2 + 1;
    } else if (c < 8) {
      dst = s.fluid[1];
      blk = c - 5 + 1;
    } else if (c < 11) {
      dst = s.em;
      blk = c - 8;
    } else if (c < 14) {
      dst = s.em;
      blk = c - 11 + 3;
    } else if (c < 17) {
      dst = s.ext_em;
      blk = c - 14 + 3;
    } else if (c < 20) {
      dst = s.app_accel[0];
      blk = c - 17;
    } else {
      dst = s.app_current;
      blk = c - 20;
    }
    gkyl_array_set_offset(dst, 1.0, tmp, blk * nb);
  }
  // fluid energies: a varying internal energy on top of the bulk kinetic energy
  gkyl_array_set_offset(tmp, 1.0, all, 0);
  gkyl_array_set_offset(s.fluid[0], 2.5, tmp, 4 * nb);
  gkyl_array_set_offset(tmp, 1.0, all, 1 * nb);
  gkyl_array_set_offset(s.fluid[1], 1.5, tmp, 4 * nb);
  gkyl_array_release(tmp);
  gkyl_array_release(all);
}

static void
test_varying(int cdim)
{
  struct setup s;
  setup_init(&s, cdim, 2);
  int nb = s.nc;
  double qbym[2] = {-1.0, 0.5}, eps0 = 1.3, dt = 0.21;
  set_varying_state(&s, cdim);
  struct gkyl_array *M_old[2] = {
    mkarr(5 * nb, s.local_ext.volume), mkarr(5 * nb, s.local_ext.volume)
  };
  struct gkyl_array *em_old = mkarr(8 * nb, s.local_ext.volume);
  gkyl_array_copy(M_old[0], s.fluid[0]);
  gkyl_array_copy(M_old[1], s.fluid[1]);
  gkyl_array_copy(em_old, s.em);

  struct gkyl_dg_calc_fluid_em_coupling *up =
    gkyl_dg_calc_fluid_em_coupling_inew(&(struct gkyl_dg_calc_fluid_em_coupling_inp){
      .conf_basis = &s.basis,
      .mem_range = &s.local,
      .num_species = 2,
      .qbym = qbym,
      .epsilon0 = eps0,
      .is_field_static = false,
      .use_gpu = false,
    });
  const struct gkyl_array *aa[2] = {s.app_accel[0], s.app_accel[1]};
  struct gkyl_array *ff[2] = {s.fluid[0], s.fluid[1]};
  gkyl_dg_calc_fluid_em_coupling_advance(up, dt, aa, s.ext_em, s.app_current, NULL, ff, s.em);

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
    const double *ext = gkyl_array_cfetch(s.ext_em, loc),
                 *Ja = gkyl_array_cfetch(s.app_current, loc);
    // Ampere's law mode by mode: eps0 (E_new - E_old) = -dt (sum_s qbym_s Mbar_s + J_app)
    for (int c = 0; c < 3; ++c) {
      for (int k = 0; k < nb; ++k) {
        double lhs = eps0 * (em_n[c * nb + k] - em_o[c * nb + k]);
        double rhs = -dt * Ja[c * nb + k];
        for (int sp = 0; sp < 2; ++sp) {
          const double *Mn = gkyl_array_cfetch(s.fluid[sp], loc),
                       *Mo = gkyl_array_cfetch(M_old[sp], loc);
          rhs -= dt * qbym[sp] * 0.5 * (Mn[(1 + c) * nb + k] + Mo[(1 + c) * nb + k]);
        }
        max_ferr = fmax(max_ferr, fabs(lhs - rhs));
      }
    }
    // Momentum: modal increment vs exact weak projection of dt qbym (rho Ebar + Mbar x B)
    // + dt rho (qbym E_ext + a), by 4-point quadrature of the p=1 fields' products.
    for (int sp = 0; sp < 2; ++sp) {
      const double *Mn = gkyl_array_cfetch(s.fluid[sp], loc),
                   *Mo = gkyl_array_cfetch(M_old[sp], loc);
      const double *rho = Mn,
                   *acc =
                     gkyl_array_cfetch(s.app_accel[sp], loc); // density: leading block, unchanged
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
          Eb[c] =
            0.5 * (s.basis.eval_expand(z, &em_n[c * nb]) + s.basis.eval_expand(z, &em_o[c * nb]));
          B[c] = s.basis.eval_expand(z, &em_o[(3 + c) * nb]) +
                 s.basis.eval_expand(z, &ext[(3 + c) * nb]);
          Mb[c] = 0.5 * (s.basis.eval_expand(z, &Mn[(1 + c) * nb]) +
                         s.basis.eval_expand(z, &Mo[(1 + c) * nb]));
          Ee[c] = s.basis.eval_expand(z, &ext[c * nb]);
          aq[c] = s.basis.eval_expand(z, &acc[c * nb]);
        }
        double F[3] = {
          Mb[1] * B[2] - Mb[2] * B[1], Mb[2] * B[0] - Mb[0] * B[2], Mb[0] * B[1] - Mb[1] * B[0]
        };
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
        double En = s.basis.eval_expand(z, &em_n[c * nb]),
               Eo = s.basis.eval_expand(z, &em_o[c * nb]);
        de += 0.5 * eps0 * (En * En - Eo * Eo);
        work -= dt * 0.5 * (En + Eo) * s.basis.eval_expand(z, &Ja[c * nb]);
      }
      for (int sp = 0; sp < 2; ++sp) {
        const double *Mn = gkyl_array_cfetch(s.fluid[sp], loc),
                     *Mo = gkyl_array_cfetch(M_old[sp], loc);
        const double *acc = gkyl_array_cfetch(s.app_accel[sp], loc);
        double rq = s.basis.eval_expand(z, Mn), dke = 0.0;
        for (int c = 0; c < 3; ++c) {
          double a = s.basis.eval_expand(z, &Mn[(1 + c) * nb]),
                 b = s.basis.eval_expand(z, &Mo[(1 + c) * nb]);
          dke += 0.5 * (a * a - b * b) / rq;
          work += dt * 0.5 * (a + b) *
                  (qbym[sp] * s.basis.eval_expand(z, &ext[c * nb]) +
                   s.basis.eval_expand(z, &acc[c * nb]));
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
  TEST_MSG(
    "cdim=%d: momentum increment vs exact weak projection %g (amplitude %g)", cdim, max_err, max_amp
  );
  TEST_CHECK(max_ferr < 1e-12);
  TEST_MSG("cdim=%d: Ampere residual %g", cdim, max_ferr);
  TEST_CHECK(max_eerr < 1e-12);
  TEST_MSG("cdim=%d: nodal energy residual %g", cdim, max_eerr);
  TEST_CHECK(max_Eferr < 1e-12);
  TEST_MSG("cdim=%d: fluid energy vs bulk kinetic energy change %g", cdim, max_Eferr);

  gkyl_array_release(em_old);
  gkyl_array_release(M_old[0]);
  gkyl_array_release(M_old[1]);
  gkyl_dg_calc_fluid_em_coupling_release(up);
  setup_release(&s, 2);
}

// A species without an energy (isothermal fluid): the momenta, the field and the energy of the
// other species are those of the solve in which both carry one, and its fifth block is not
// written.
static void
test_no_energy(int cdim)
{
  struct setup s_ref, s_iso;
  setup_init(&s_ref, cdim, 2);
  setup_init(&s_iso, cdim, 2);
  set_varying_state(&s_ref, cdim);
  set_varying_state(&s_iso, cdim);
  int nb = s_ref.nc;
  double qbym[2] = {-1.0, 0.5}, eps0 = 1.3, dt = 0.21, sentinel = -7.5;
  // the fifth block of the isothermal species holds a sentinel
  struct gkyl_array *block = mkarr(nb, s_iso.local_ext.volume);
  gkyl_array_clear(block, 0.0);
  gkyl_array_set_offset(s_iso.fluid[0], 1.0, block, 4 * nb);
  gkyl_array_release(block);
  set_const_block(s_iso.fluid[0], 4, nb, sentinel, cdim);
  struct gkyl_array *fluid_old = mkarr(5 * nb, s_iso.local_ext.volume);
  gkyl_array_copy(fluid_old, s_iso.fluid[0]);

  bool has_energy[2] = {false, true};
  struct gkyl_dg_calc_fluid_em_coupling *up_ref =
    gkyl_dg_calc_fluid_em_coupling_inew(&(struct gkyl_dg_calc_fluid_em_coupling_inp){
      .conf_basis = &s_ref.basis,
      .mem_range = &s_ref.local,
      .num_species = 2,
      .qbym = qbym,
      .epsilon0 = eps0,
      .is_field_static = false,
      .use_gpu = false,
    });
  struct gkyl_dg_calc_fluid_em_coupling *up_iso =
    gkyl_dg_calc_fluid_em_coupling_inew(&(struct gkyl_dg_calc_fluid_em_coupling_inp){
      .conf_basis = &s_iso.basis,
      .mem_range = &s_iso.local,
      .num_species = 2,
      .qbym = qbym,
      .has_energy = has_energy,
      .epsilon0 = eps0,
      .is_field_static = false,
      .use_gpu = false,
    });
  const struct gkyl_array *aa_ref[2] = {s_ref.app_accel[0], s_ref.app_accel[1]};
  struct gkyl_array *ff_ref[2] = {s_ref.fluid[0], s_ref.fluid[1]};
  gkyl_dg_calc_fluid_em_coupling_advance(
    up_ref, dt, aa_ref, s_ref.ext_em, s_ref.app_current, NULL, ff_ref, s_ref.em
  );
  const struct gkyl_array *aa_iso[2] = {s_iso.app_accel[0], s_iso.app_accel[1]};
  struct gkyl_array *ff_iso[2] = {s_iso.fluid[0], s_iso.fluid[1]};
  gkyl_dg_calc_fluid_em_coupling_advance(
    up_iso, dt, aa_iso, s_iso.ext_em, s_iso.app_current, NULL, ff_iso, s_iso.em
  );

  double max_diff = 0.0, max_sentinel = 0.0, max_change = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s_ref.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&s_ref.local_ext, iter.idx);
    const double *f0_ref = gkyl_array_cfetch(s_ref.fluid[0], loc);
    const double *f0_iso = gkyl_array_cfetch(s_iso.fluid[0], loc);
    const double *f1_ref = gkyl_array_cfetch(s_ref.fluid[1], loc);
    const double *f1_iso = gkyl_array_cfetch(s_iso.fluid[1], loc);
    const double *em_ref = gkyl_array_cfetch(s_ref.em, loc);
    const double *em_iso = gkyl_array_cfetch(s_iso.em, loc);
    for (int k = 0; k < 4 * nb; ++k) {
      max_diff = fmax(max_diff, fabs(f0_ref[k] - f0_iso[k]));
    }
    for (int k = 0; k < 5 * nb; ++k) {
      max_diff = fmax(max_diff, fabs(f1_ref[k] - f1_iso[k]));
    }
    for (int k = 0; k < 3 * nb; ++k) {
      max_diff = fmax(max_diff, fabs(em_ref[k] - em_iso[k]));
    }
    max_sentinel = fmax(max_sentinel, fabs(f0_iso[4 * nb] - sentinel * pow(sqrt(2.0), cdim)));
    for (int k = 1; k < nb; ++k) {
      max_sentinel = fmax(max_sentinel, fabs(f0_iso[4 * nb + k]));
    }
    // the solve changes the momentum (the comparison is not of untouched states)
    const double *f0_old = gkyl_array_cfetch(fluid_old, loc);
    max_change = fmax(max_change, fabs(f0_iso[nb] - f0_old[nb]));
  }
  TEST_CHECK(max_diff == 0.0);
  TEST_MSG("largest difference from the solve with energies: %g", max_diff);
  TEST_CHECK(max_sentinel == 0.0);
  TEST_MSG("fifth block of the isothermal species changed by %g", max_sentinel);
  TEST_CHECK(max_change > 0.0);

  gkyl_array_release(fluid_old);
  gkyl_dg_calc_fluid_em_coupling_release(up_ref);
  gkyl_dg_calc_fluid_em_coupling_release(up_iso);
  setup_release(&s_ref, 2);
  setup_release(&s_iso, 2);
}

void
test_static_E_push_1x()
{
  test_static_E_push(1);
}
void
test_static_E_push_2x()
{
  test_static_E_push(2);
}
void
test_varying_1x()
{
  test_varying(1);
}
void
test_varying_2x()
{
  test_varying(2);
}
void
test_varying_3x()
{
  test_varying(3);
}
void
test_no_energy_1x()
{
  test_no_energy(1);
}
void
test_no_energy_2x()
{
  test_no_energy(2);
}
void
test_no_energy_3x()
{
  test_no_energy(3);
}

// Resistive layer, 1x: a neutral species (q/m = 0) decouples the fluid, so E alone is damped by
// exp(-sigma dt) at each node; sigma varies linearly in each cell, so the two nodes have different
// factors. B and the fluid are unchanged.
static void
test_resistive_layer(double dt)
{
  struct setup s;
  setup_init(&s, 1, 1);
  int nb = s.nc;
  double qbym[1] = {0.0}, eps0 = 1.0;
  double rho = 3.0, Mx = 1.0, E0 = 7.0, Ex = 2.0, Ey = -0.5, Ez = 1.0, Bz = 0.8;
  double s0 = 1.5, s1 = 0.6; // sigma = s0 + s1 eta in every cell, eta in [-1, 1]
  set_const_block(s.fluid[0], 0, nb, rho, 1);
  set_const_block(s.fluid[0], 1, nb, Mx, 1);
  set_const_block(s.fluid[0], 4, nb, E0, 1);
  set_const_block(s.em, 0, nb, Ex, 1);
  set_const_block(s.em, 1, nb, Ey, 1);
  set_const_block(s.em, 2, nb, Ez, 1);
  set_const_block(s.em, 5, nb, Bz, 1);
  struct gkyl_array *sigma = mkarr(2, s.local_ext.volume); // values at the 2 nodes
  gkyl_array_clear(sigma, 0.0);
  gkyl_array_shiftc(sigma, s0 - s1 / sqrt(3.0), 0);
  gkyl_array_shiftc(sigma, s0 + s1 / sqrt(3.0), 1);

  struct gkyl_dg_calc_fluid_em_coupling *up =
    gkyl_dg_calc_fluid_em_coupling_inew(&(struct gkyl_dg_calc_fluid_em_coupling_inp){
      .conf_basis = &s.basis,
      .mem_range = &s.local,
      .num_species = 1,
      .qbym = qbym,
      .epsilon0 = eps0,
      .is_field_static = false,
      .use_gpu = false,
    });
  const struct gkyl_array *aa[1] = {s.app_accel[0]};
  struct gkyl_array *ff[1] = {s.fluid[0]};
  gkyl_dg_calc_fluid_em_coupling_advance(up, dt, aa, s.ext_em, s.app_current, sigma, ff, s.em);

  double E_in[3] = {Ex, Ey, Ez};
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&s.local_ext, iter.idx);
    const double *em = gkyl_array_cfetch(s.em, loc);
    for (int nd = 0; nd < 2; ++nd) {
      double sgn = nd == 0 ? -1.0 : 1.0;
      double sig = s0 + sgn * s1 / sqrt(3.0);
      double g = exp(-sig * dt);
      for (int c = 0; c < 3; ++c) {
        // value of the p=1 expansion at the node eta = sgn/sqrt(3)
        double v = (em[c * nb] + sgn * em[c * nb + 1]) / sqrt(2.0);
        TEST_CHECK(gkyl_compare(v, g * E_in[c], 1e-14));
      }
    }
    TEST_CHECK(gkyl_compare(avg(s.em, loc, 5, nb, 1), Bz, 1e-15));
    TEST_CHECK(gkyl_compare(avg(s.fluid[0], loc, 1, nb, 1), Mx, 1e-15));
    TEST_CHECK(gkyl_compare(avg(s.fluid[0], loc, 4, nb, 1), E0, 1e-15));
  }
  gkyl_array_release(sigma);
  gkyl_dg_calc_fluid_em_coupling_release(up);
  setup_release(&s, 1);
}

// Resistive layer coefficients against references computed with 80-digit decimal arithmetic
// (series for x < 1): g2 = exp(-x), phi = (1 - g2)/x, phi2 = (1 - phi)/x.
static void
test_layer_coeffs(void)
{
  const double ref[][4] = {
    {1e-300, 1, 1, 0.5},
    {1e-12, 0.99999999999900002, 0.99999999999949996, 0.49999999999983336},
    {1e-06, 0.99999900000050002, 0.99999950000016669, 0.499999833333375},
    {0.01, 0.99004983374916811, 0.99501662508319466, 0.49833749168053576},
    {0.5, 0.60653065971263342, 0.78693868057473315, 0.4261226388505337},
    {0.999999999999, 0.36787944117181021, 0.6321205588288219, 0.36787944117154597},
    {1.0, 0.36787944117144233, 0.63212055882855767, 0.36787944117144233},
    {1.000000000001, 0.36787944117107441, 0.63212055882829343, 0.36787944117133869},
    {3.0, 0.049787068367863944, 0.31673764387737868, 0.22775411870754045},
    {40.0, 4.2483542552915889e-18, 0.025000000000000001, 0.024375000000000001},
    {800.0, 0, 0.00125, 0.0012484375}
  };
  for (int i = 0; i < (int)(sizeof ref / sizeof ref[0]); ++i) {
    double g2, phi, phi2;
    fluid_em_coupling_layer_coeffs(ref[i][0], &g2, &phi, &phi2);
    TEST_CHECK(gkyl_compare(g2, ref[i][1], 1e-15));
    TEST_CHECK(gkyl_compare(phi, ref[i][2], 1e-15));
    TEST_CHECK(gkyl_compare(phi2, ref[i][3], 1e-15));
    TEST_MSG("x = %g: g2 %.17g phi %.17g phi2 %.17g", ref[i][0], g2, phi, phi2);
  }
  double g2, phi, phi2;
  fluid_em_coupling_layer_coeffs(0.0, &g2, &phi, &phi2);
  TEST_CHECK(g2 == 1.0 && phi == 1.0 && phi2 == 0.5);
}

// Resistive layer balance, 1x: a neutral species and an applied current J leave E = -J/(eps0 sigma)
// steady (sigma uniform), for sigma dt below and well above 1.
static void
test_resistive_balance(double dt)
{
  struct setup s;
  setup_init(&s, 1, 1);
  int nb = s.nc;
  double qbym[1] = {0.0}, eps0 = 1.3, sig = 2.5;
  double J[3] = {0.7, -1.1, 0.4};
  set_const_block(s.fluid[0], 0, nb, 1.0, 1);
  for (int c = 0; c < 3; ++c) {
    set_const_block(s.app_current, c, nb, J[c], 1);
    set_const_block(s.em, c, nb, -J[c] / (eps0 * sig), 1);
  }
  struct gkyl_array *sigma = mkarr(nb, s.local_ext.volume); // values at the nb = 2 nodes
  gkyl_array_clear(sigma, 0.0);
  for (int n = 0; n < nb; ++n) {
    gkyl_array_shiftc(sigma, sig, n);
  }

  struct gkyl_dg_calc_fluid_em_coupling *up =
    gkyl_dg_calc_fluid_em_coupling_inew(&(struct gkyl_dg_calc_fluid_em_coupling_inp){
      .conf_basis = &s.basis,
      .mem_range = &s.local,
      .num_species = 1,
      .qbym = qbym,
      .epsilon0 = eps0,
      .is_field_static = false,
      .use_gpu = false,
    });
  const struct gkyl_array *aa[1] = {s.app_accel[0]};
  struct gkyl_array *ff[1] = {s.fluid[0]};
  gkyl_dg_calc_fluid_em_coupling_advance(up, dt, aa, s.ext_em, s.app_current, sigma, ff, s.em);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&s.local_ext, iter.idx);
    for (int c = 0; c < 3; ++c) {
      TEST_CHECK(gkyl_compare(avg(s.em, loc, c, nb, 1), -J[c] / (eps0 * sig), 1e-14));
    }
  }
  gkyl_array_release(sigma);
  gkyl_dg_calc_fluid_em_coupling_release(up);
  setup_release(&s, 1);
}

static void
eval_sigma(double t, const double *xn, double *restrict fout, void *ctx)
{
  const double *c = ctx; // {cdim, scale}
  double y = c[0] > 1 ? xn[1] : 0.3, z = c[0] > 2 ? xn[2] : -0.2;
  fout[0] = c[1] * (2.0 + sin(3.0 * xn[0]) + 0.5 * cos(2.0 * y) * cos(z));
}

// Values of f at the 2^cdim Gauss-Legendre nodes of each cell (row-major, the last direction
// fastest), the form in which the source solve takes sigma.
static void
sample_at_nodes(const struct setup *s, evalf_t f, void *ctx, struct gkyl_array *out)
{
  int cdim = s->cdim, node_lower[3] = {0, 0, 0}, node_upper[3] = {1, 1, 1};
  struct gkyl_range node_range;
  gkyl_range_init(&node_range, cdim, node_lower, node_upper);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s->local_ext);
  while (gkyl_range_iter_next(&iter)) {
    double xc[3], x[3];
    gkyl_rect_grid_cell_center(&s->grid, iter.idx, xc);
    double *o = gkyl_array_fetch(out, gkyl_range_idx(&s->local_ext, iter.idx));
    struct gkyl_range_iter node_iter;
    gkyl_range_iter_init(&node_iter, &node_range);
    while (gkyl_range_iter_next(&node_iter)) {
      for (int d = 0; d < cdim; ++d) {
        x[d] = xc[d] + 0.5 * s->grid.dx[d] * gkyl_gauss_ordinates_2[node_iter.idx[d]];
      }
      f(0.0, x, &o[gkyl_range_idx(&node_range, node_iter.idx)], ctx);
    }
  }
}

static double
layer_phi(double x)
{
  return x > 0.0 ? -expm1(-x) / x : 1.0;
}

// Resistive layer with the varying state (two species, E, B, external fields, applied
// acceleration and current) and a varying sigma: at each node the change of
// sum_s |M_s|^2/(2 rho_s) + eps0 |E|^2/2 equals dt sum_s Mbar_s.(qbym_s E_ext + a_s) - dt J_app.Ebar
// - eps0 sigma int |E|^2 dt, with E(t) = E^n exp(-sigma t) + S (1 - exp(-sigma t))/sigma the exact
// solution for the midpoint current S, and Ebar its time average.
static void
test_resistive_energy(int cdim, double scale)
{
  struct setup s;
  setup_init(&s, cdim, 2);
  int nb = s.nc;
  double qbym[2] = {-1.0, 0.5}, eps0 = 1.3, dt = 0.21;
  set_varying_state(&s, cdim);
  struct gkyl_array *sigma = mkarr(nb, s.local_ext.volume); // values at the nb = 2^cdim nodes
  double sig_ctx[2] = {cdim, scale};
  sample_at_nodes(&s, eval_sigma, sig_ctx, sigma);
  struct gkyl_array *M_old[2] = {
    mkarr(5 * nb, s.local_ext.volume), mkarr(5 * nb, s.local_ext.volume)
  };
  struct gkyl_array *em_old = mkarr(8 * nb, s.local_ext.volume);
  gkyl_array_copy(M_old[0], s.fluid[0]);
  gkyl_array_copy(M_old[1], s.fluid[1]);
  gkyl_array_copy(em_old, s.em);

  struct gkyl_dg_calc_fluid_em_coupling *up =
    gkyl_dg_calc_fluid_em_coupling_inew(&(struct gkyl_dg_calc_fluid_em_coupling_inp){
      .conf_basis = &s.basis,
      .mem_range = &s.local,
      .num_species = 2,
      .qbym = qbym,
      .epsilon0 = eps0,
      .is_field_static = false,
      .use_gpu = false,
    });
  const struct gkyl_array *aa[2] = {s.app_accel[0], s.app_accel[1]};
  struct gkyl_array *ff[2] = {s.fluid[0], s.fluid[1]};
  gkyl_dg_calc_fluid_em_coupling_advance(up, dt, aa, s.ext_em, s.app_current, sigma, ff, s.em);

  double max_eerr = 0.0, max_loss = 0.0, min_x = 1e300, max_x = 0.0;
  int node_lower[3] = {0, 0, 0}, node_upper[3] = {1, 1, 1};
  struct gkyl_range node_range;
  gkyl_range_init(&node_range, cdim, node_lower, node_upper);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&s.local_ext, iter.idx);
    const double *em_n = gkyl_array_cfetch(s.em, loc), *em_o = gkyl_array_cfetch(em_old, loc);
    const double *ext = gkyl_array_cfetch(s.ext_em, loc);
    const double *Ja = gkyl_array_cfetch(s.app_current, loc);
    const double *sg = gkyl_array_cfetch(sigma, loc);
    struct gkyl_range_iter node_iter;
    gkyl_range_iter_init(&node_iter, &node_range);
    while (gkyl_range_iter_next(&node_iter)) {
      double z[3];
      for (int d = 0; d < cdim; ++d) {
        z[d] = gkyl_gauss_ordinates_2[node_iter.idx[d]];
      }
      double sig = sg[gkyl_range_idx(&node_range, node_iter.idx)], x = sig * dt;
      min_x = fmin(min_x, x);
      max_x = fmax(max_x, x);
      double g2 = exp(-x), phi = layer_phi(x), phi2x = layer_phi(2.0 * x);
      double I1 = dt * phi2x, I2 = dt * dt * (phi - phi2x) / x,
             I3 = dt * dt * dt * (1.0 - 2.0 * phi + phi2x) / (x * x);
      double de = 0.0, work = 0.0, intE2 = 0.0;
      for (int c = 0; c < 3; ++c) {
        double En = s.basis.eval_expand(z, &em_n[c * nb]),
               Eo = s.basis.eval_expand(z, &em_o[c * nb]);
        double S = (En - g2 * Eo) / (dt * phi); // midpoint current of the exact solution
        double Ebar = phi * Eo + (1.0 - phi) / x * dt * S;
        de += 0.5 * eps0 * (En * En - Eo * Eo);
        work -= dt * Ebar * s.basis.eval_expand(z, &Ja[c * nb]);
        intE2 += Eo * Eo * I1 + 2.0 * Eo * S * I2 + S * S * I3;
      }
      for (int sp = 0; sp < 2; ++sp) {
        const double *Mn = gkyl_array_cfetch(s.fluid[sp], loc),
                     *Mo = gkyl_array_cfetch(M_old[sp], loc);
        const double *acc = gkyl_array_cfetch(s.app_accel[sp], loc);
        double rq = s.basis.eval_expand(z, Mn);
        for (int c = 0; c < 3; ++c) {
          double a = s.basis.eval_expand(z, &Mn[(1 + c) * nb]),
                 b = s.basis.eval_expand(z, &Mo[(1 + c) * nb]);
          de += 0.5 * (a * a - b * b) / rq;
          work += dt * 0.5 * (a + b) *
                  (qbym[sp] * s.basis.eval_expand(z, &ext[c * nb]) +
                   s.basis.eval_expand(z, &acc[c * nb]));
        }
      }
      double loss = eps0 * sig * intE2;
      max_loss = fmax(max_loss, loss);
      max_eerr = fmax(max_eerr, fabs(de - (work - loss)));
    }
  }
  TEST_CHECK(max_loss > 1e-3);
  TEST_CHECK(max_eerr < 1e-12);
  TEST_MSG(
    "cdim=%d, sigma dt in [%g, %g]: nodal energy residual %g (largest loss %g)", cdim, min_x, max_x,
    max_eerr, max_loss
  );
  gkyl_array_release(sigma);
  gkyl_array_release(em_old);
  gkyl_array_release(M_old[0]);
  gkyl_array_release(M_old[1]);
  gkyl_dg_calc_fluid_em_coupling_release(up);
  setup_release(&s, 2);
}

// Resistive layer and no layer, 2x varying state: sigma = 0 gives bitwise the result of
// sigma = NULL, and sigma dt = 1e-13 agrees with it to round-off.
static void
test_resistive_null(void)
{
  int cdim = 2;
  double qbym[2] = {-1.0, 0.5}, eps0 = 1.3, dt = 0.21;
  double sig_vals[2] = {0.0, 1e-13 / 0.21};
  struct setup s0;
  setup_init(&s0, cdim, 2);
  set_varying_state(&s0, cdim);
  int nb = s0.nc;
  struct gkyl_dg_calc_fluid_em_coupling *up =
    gkyl_dg_calc_fluid_em_coupling_inew(&(struct gkyl_dg_calc_fluid_em_coupling_inp){
      .conf_basis = &s0.basis,
      .mem_range = &s0.local,
      .num_species = 2,
      .qbym = qbym,
      .epsilon0 = eps0,
      .is_field_static = false,
      .use_gpu = false,
    });
  const struct gkyl_array *aa0[2] = {s0.app_accel[0], s0.app_accel[1]};
  struct gkyl_array *ff0[2] = {s0.fluid[0], s0.fluid[1]};
  gkyl_dg_calc_fluid_em_coupling_advance(up, dt, aa0, s0.ext_em, s0.app_current, NULL, ff0, s0.em);

  for (int k = 0; k < 2; ++k) {
    struct setup s1;
    setup_init(&s1, cdim, 2);
    set_varying_state(&s1, cdim);
    struct gkyl_array *sigma = mkarr(nb, s1.local_ext.volume); // values at the nb = 4 nodes
    gkyl_array_clear(sigma, 0.0);
    for (int n = 0; n < nb; ++n) {
      gkyl_array_shiftc(sigma, sig_vals[k], n);
    }
    const struct gkyl_array *aa1[2] = {s1.app_accel[0], s1.app_accel[1]};
    struct gkyl_array *ff1[2] = {s1.fluid[0], s1.fluid[1]};
    gkyl_dg_calc_fluid_em_coupling_advance(
      up, dt, aa1, s1.ext_em, s1.app_current, sigma, ff1, s1.em
    );
    long n_em = s0.em->size * s0.em->ncomp, n_fl = s0.fluid[0]->size * s0.fluid[0]->ncomp;
    if (k == 0) {
      TEST_CHECK(memcmp(s0.em->data, s1.em->data, n_em * sizeof(double)) == 0);
      for (int sp = 0; sp < 2; ++sp) {
        TEST_CHECK(memcmp(s0.fluid[sp]->data, s1.fluid[sp]->data, n_fl * sizeof(double)) == 0);
      }
    } else {
      double max_diff = 0.0;
      const double *a = s0.em->data, *b = s1.em->data;
      for (long i = 0; i < n_em; ++i) {
        max_diff = fmax(max_diff, fabs(a[i] - b[i]));
      }
      for (int sp = 0; sp < 2; ++sp) {
        a = s0.fluid[sp]->data;
        b = s1.fluid[sp]->data;
        for (long i = 0; i < n_fl; ++i) {
          max_diff = fmax(max_diff, fabs(a[i] - b[i]));
        }
      }
      TEST_CHECK(max_diff < 1e-12);
      TEST_MSG("sigma dt = 1e-13: largest difference to no layer %g", max_diff);
    }
    gkyl_array_release(sigma);
    setup_release(&s1, 2);
  }
  gkyl_dg_calc_fluid_em_coupling_release(up);
  setup_release(&s0, 2);
}

void
test_resistive_layer_1x()
{
  test_resistive_layer(0.37);
}
void
test_resistive_layer_stiff()
{
  test_resistive_layer(4.0);
}

void
test_resistive_balance_1x()
{
  test_resistive_balance(0.3);
}
void
test_resistive_balance_stiff()
{
  test_resistive_balance(7.0);
}
void
test_resistive_energy_1x()
{
  test_resistive_energy(1, 1.0);
}
void
test_resistive_energy_2x()
{
  test_resistive_energy(2, 1.0);
}
void
test_resistive_energy_3x()
{
  test_resistive_energy(3, 1.0);
}
void
test_resistive_energy_stiff()
{
  test_resistive_energy(2, 20.0);
}

TEST_LIST = {
  {"static_E_push_1x", test_static_E_push_1x},
  {"static_E_push_2x", test_static_E_push_2x},
  {"plasma_oscillation", test_plasma_oscillation},
  {"cyclotron", test_cyclotron},
  {"varying_1x", test_varying_1x},
  {"varying_2x", test_varying_2x},
  {"varying_3x", test_varying_3x},
  {"no_energy_1x", test_no_energy_1x},
  {"no_energy_2x", test_no_energy_2x},
  {"no_energy_3x", test_no_energy_3x},
  {"resistive_layer_1x", test_resistive_layer_1x},
  {"resistive_layer_stiff", test_resistive_layer_stiff},
  {"layer_coeffs", test_layer_coeffs},
  {"resistive_balance_1x", test_resistive_balance_1x},
  {"resistive_balance_stiff", test_resistive_balance_stiff},
  {"resistive_energy_1x", test_resistive_energy_1x},
  {"resistive_energy_2x", test_resistive_energy_2x},
  {"resistive_energy_3x", test_resistive_energy_3x},
  {"resistive_energy_stiff", test_resistive_energy_stiff},
  {"resistive_null", test_resistive_null},
  {NULL, NULL}
};
