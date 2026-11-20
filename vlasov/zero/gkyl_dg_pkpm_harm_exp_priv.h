// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_dg_conf_nodal_priv.h>
#include <gkyl_dg_pkpm_harm_exp.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_pkpm_layout.h>
#include <gkyl_util.h>

// Coefficients of one configuration node: (q/m) B_n, whose magnitude sets the mixing rate and whose
// direction b_n the plane the mixing damps (the component along b_n is not damped).
enum {
  PKPM_HX_OM = 0, // (q/m) B_n (3)
  PKPM_HX_NC = 3 // per configuration node
};

// Largest number of phase-space nodes of a configuration cell (w nodes per configuration node).
#define PKPM_HX_MAX_NODES PKPM_MAX_PHASE_BASIS

// A transform between the expansions of a phase-space cell and values at its nodes.
typedef void (*pkpm_harm_exp_kern_t)(const double *in, double *GKYL_RESTRICT out);

// Nodal transforms (p=1 serendipity only; indexed [cdim-1]): the first-harmonic blocks at the
// phase nodes and back, F and G at the phase nodes, and the w integrals at the configuration nodes.
GKYL_CU_D static const pkpm_harm_exp_kern_t ser_pkpm_harm_exp_to_nodes_kernels[] = {
  pkpm_harm_exp_to_nodes_1x1v_ser_p1, pkpm_harm_exp_to_nodes_2x1v_ser_p1,
  pkpm_harm_exp_to_nodes_3x1v_ser_p1
};

GKYL_CU_D static const pkpm_harm_exp_kern_t ser_pkpm_harm_exp_to_modal_kernels[] = {
  pkpm_harm_exp_to_modal_1x1v_ser_p1, pkpm_harm_exp_to_modal_2x1v_ser_p1,
  pkpm_harm_exp_to_modal_3x1v_ser_p1
};

GKYL_CU_D static const pkpm_harm_exp_kern_t ser_pkpm_harm_exp_fg_to_nodes_kernels[] = {
  pkpm_harm_exp_fg_to_nodes_1x1v_ser_p1, pkpm_harm_exp_fg_to_nodes_2x1v_ser_p1,
  pkpm_harm_exp_fg_to_nodes_3x1v_ser_p1
};

GKYL_CU_D static const pkpm_harm_exp_kern_t ser_pkpm_harm_exp_wint_kernels[] = {
  pkpm_harm_exp_wint_1x1v_ser_p1, pkpm_harm_exp_wint_2x1v_ser_p1, pkpm_harm_exp_wint_3x1v_ser_p1
};

struct gkyl_dg_pkpm_harm_exp {
  int cdim; // Configuration-space dimensions.
  int num_conf_nodes; // Nodes of a configuration cell.
  int num_phase_nodes; // Nodes of a phase cell: three w nodes per configuration node.
  pkpm_harm_exp_kern_t to_nodes; // First-harmonic blocks at the phase nodes.
  pkpm_harm_exp_kern_t to_modal; // First-harmonic blocks from their values at the phase nodes.
  pkpm_harm_exp_kern_t fg_to_nodes; // F and G at the phase nodes.
  pkpm_harm_exp_kern_t wint; // w integrals of F and M0_perp at the configuration nodes.
  dg_conf_nodal_kern_t conf_to_nodal; // Configuration-space expansion at the configuration nodes.
  double qbym; // Charge to mass ratio.
  double mix; // Gyrophase-mixing rate in units of |Omega| (0: off).
  double half_dw; // Half the velocity cell width.
  bool use_gpu; // Whether the updater runs on the GPU.
  uint32_t flags;
  struct gkyl_dg_pkpm_harm_exp *on_dev; // pointer to itself or device data
};

// Set the sizes and parameters of the updater (everything but the kernels and the device data).
void gkyl_dg_pkpm_harm_exp_setup(
  struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_dg_pkpm_harm_exp_inp *inp
);

// Choose the kernels of the updater for its number of configuration-space dimensions.
GKYL_CU_D static void
pkpm_hx_choose_kernels(struct gkyl_dg_pkpm_harm_exp *up)
{
  int cdim = up->cdim;
  up->to_nodes = ser_pkpm_harm_exp_to_nodes_kernels[cdim - 1];
  up->to_modal = ser_pkpm_harm_exp_to_modal_kernels[cdim - 1];
  up->fg_to_nodes = ser_pkpm_harm_exp_fg_to_nodes_kernels[cdim - 1];
  up->wint = ser_pkpm_harm_exp_wint_kernels[cdim - 1];
  up->conf_to_nodal = dg_conf_nodal_conf_to_nodal_kernels.kernels[cdim - 1];
}

// |(q/m) B_n| and its direction (zero where B_n = 0).
GKYL_CU_DH static inline double
pkpm_hx_axis(const double *om, double *oh)
{
  double mag = sqrt(om[0] * om[0] + om[1] * om[1] + om[2] * om[2]);
  for (int a = 0; a < 3; ++a) {
    oh[a] = mag > 0.0 ? om[a] / mag : 0.0;
  }
  return mag;
}

// out = d v_perp + v_par for the plane perpendicular to oh (the integrating factor e^(sL) of a node
// damped at the rate nu has d = e^(-s nu)).
GKYL_CU_DH static inline void
pkpm_hx_damp(const double *oh, double d, const double *v, double *out)
{
  double vb = (1.0 - d) * (v[0] * oh[0] + v[1] * oh[1] + v[2] * oh[2]);
  for (int a = 0; a < 3; ++a) {
    out[a] = d * v[a] + vb * oh[a];
  }
}

// Gyrophase-mixing rate of one phase node: mix |Omega| (1 - 1/r)_+, with r the larger of the
// realizability ratios r_M = sqrt2 |M0_perp|/sqrt(F G) and r_N = |N_perp| sqrt(F)/G^(3/2), which are
// <= 1 when the one-harmonic reconstruction is non-negative (r = infinity where F or G <= 0). Capped at
// ln(r)/h, the rate that takes the unforced state to the bound in one step.
GKYL_CU_DH static inline double
pkpm_hx_rate(double mix, double mag, double h, double F, double G, const double *m0, const double *n)
{
  double rate = mix * mag;
  if (rate <= 0.0) {
    return 0.0;
  }
  if (F <= 0.0 || G <= 0.0) {
    return rate;
  }
  double am = sqrt(m0[0] * m0[0] + m0[1] * m0[1] + m0[2] * m0[2]);
  double an = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  double r = fmax(sqrt(2.0) * am / sqrt(F * G), an * sqrt(F) / (G * sqrt(G)));
  if (r <= 1.0) {
    return 0.0;
  }
  return fmin(rate * (1.0 - 1.0 / r), log(r) / h);
}

// Coefficients of one configuration cell, from the magnetic field of em = [E (3), B (3), ...].
GKYL_CU_DH static inline void
pkpm_hx_coef_cell(const struct gkyl_dg_pkpm_harm_exp *up, const double *em, double *coef)
{
  int num_nodes = up->num_conf_nodes;
  for (int a = 0; a < 3; ++a) {
    double Bn[PKPM_MAX_CONF_BASIS];
    up->conf_to_nodal(&em[(3 + a) * num_nodes], Bn);
    for (int n = 0; n < num_nodes; ++n) {
      coef[n * PKPM_HX_NC + PKPM_HX_OM + a] = up->qbym * Bn[n];
    }
  }
}

// Stage s of one phase cell (see gkyl_dg_pkpm_harm_exp_advance). At stage 1 (f0 the stage state) it
// first sets the mixing rate nu of every node, and act to 1 when a node of the cell is damped (act may
// be shared by the cells of a column: every writer writes 1). A cell without a damped node is not
// touched; at a damped node, with Y = out, X_1 = e^(hL) Y, X_2 = 3 (e^(hL/2) - 1) M + e^(-hL/2) Y and
// X_3 = (e^(hL) - 1) M / 2 + e^(hL/2) Y. The factors 3 and 1/2 are the ratios c_M / c_Y of the
// weights of the combination; with lam (node values at the configuration nodes, or 0) the
// combination has the weights w_M / w_Y = c_M / (c_Y lam), and the factors are divided by lam.
GKYL_CU_DH static inline void
pkpm_hx_advance_cell(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, const double *coef, const double *f0,
  const double *lam, double *nu, double *act, double *out
)
{
  double Mn[PKPM_HX_MAX_NODES * 6];
  if (stage == 1) {
    up->to_nodes(f0, Mn);
    double fgn[PKPM_HX_MAX_NODES * 2];
    up->fg_to_nodes(f0, fgn);
    for (int p = 0; p < up->num_phase_nodes; ++p) {
      double oh[3];
      double mag = pkpm_hx_axis(&coef[(p / PKPM_NUM_VEL_MODES) * PKPM_HX_NC + PKPM_HX_OM], oh);
      nu[p] = pkpm_hx_rate(up->mix, mag, h, fgn[2 * p], fgn[2 * p + 1], &Mn[6 * p], &Mn[6 * p + 3]);
      if (nu[p] > 0.0) {
        act[0] = 1.0;
      }
    }
  }
  bool damped = false;
  for (int p = 0; p < up->num_phase_nodes; ++p) {
    damped = damped || nu[p] > 0.0;
  }
  if (!damped) {
    return;
  }
  if (stage > 1) {
    up->to_nodes(f0, Mn);
  }
  // exponents of the factor of Y and of M (e^(sL) = e^(-s nu)), and the weight cm of e^(.) M - M
  double sy = stage == 1 ? 1.0 : (stage == 2 ? -0.5 : 0.5);
  double sm = stage == 2 ? 0.5 : 1.0;
  double cm = stage == 2 ? 3.0 : 0.5;
  double Yn[PKPM_HX_MAX_NODES * 6];
  up->to_nodes(out, Yn);
  for (int p = 0; p < up->num_phase_nodes; ++p) {
    if (nu[p] <= 0.0) {
      continue;
    }
    double oh[3];
    pkpm_hx_axis(&coef[(p / PKPM_NUM_VEL_MODES) * PKPM_HX_NC + PKPM_HX_OM], oh);
    double dy = exp(-sy * h * nu[p]), dm = exp(-sm * h * nu[p]);
    for (int v = 0; v < 2; ++v) {
      double *y = &Yn[6 * p + 3 * v];
      double x[3];
      pkpm_hx_damp(oh, dy, y, x);
      if (stage > 1) {
        const double *m = &Mn[6 * p + 3 * v];
        double mf[3];
        pkpm_hx_damp(oh, dm, m, mf);
        double cm_node = lam ? cm / lam[p / PKPM_NUM_VEL_MODES] : cm;
        for (int i = 0; i < 3; ++i) {
          x[i] += cm_node * (mf[i] - m[i]);
        }
      }
      for (int i = 0; i < 3; ++i) {
        y[i] = x[i];
      }
    }
  }
  up->to_modal(Yn, out);
}

// Add the integrals over w of F and M0_perp over one phase cell f, at the configuration nodes, to
// col[node*4 + (0: F, 1-3: M0_perp)]; summed over the cells of a column, the integrals over w.
GKYL_CU_DH static inline void
pkpm_hx_col_add(const struct gkyl_dg_pkpm_harm_exp *up, const double *f, double *col)
{
  double wi[PKPM_MAX_CONF_BASIS * 4];
  up->wint(f, wi);
  for (int i = 0; i < 4 * up->num_conf_nodes; ++i) {
    col[i] += up->half_dw * wi[i];
  }
}

// Momentum restoration of one phase cell: with m = int M0_perp dw / int F dw at the configuration
// nodes, M0_perp -= F m and N_perp -= 2 G m (the first-harmonic change of a perpendicular frame shift).
GKYL_CU_DH static inline void
pkpm_hx_restore_cell(const struct gkyl_dg_pkpm_harm_exp *up, const double *col, double *out)
{
  double Xn[PKPM_HX_MAX_NODES * 6], fgn[PKPM_HX_MAX_NODES * 2];
  up->to_nodes(out, Xn);
  up->fg_to_nodes(out, fgn);
  for (int p = 0; p < up->num_phase_nodes; ++p) {
    const double *cn = &col[(p / PKPM_NUM_VEL_MODES) * 4];
    if (cn[0] <= 0.0) {
      continue;
    }
    for (int a = 0; a < 3; ++a) {
      double m = cn[1 + a] / cn[0];
      Xn[6 * p + a] -= fgn[2 * p] * m;
      Xn[6 * p + 3 + a] -= 2.0 * fgn[2 * p + 1] * m;
    }
  }
  up->to_modal(Xn, out);
}
