#pragma once

#include <gkyl_util.h>

// Node layout of the w flux of the eight blocks [F, G, M0_perp (x, y, z), N_perp (x, y, z)] at a
// configuration point, read by the w-face and closure kernels and the node functions of
// gkyl_pkpm_harm_coeffs.h: value s < PKPM_HS_W_END is the raw value at the point of the w-flux
// coefficient whose projection is the entry pkpm_ho_wslot[s] of the stored operator
// (gkyl_pkpm_harm_ops.h; the build kernels pkpm_harm_build_* write both), followed by theta, b,
// alpha_0 and c, and the w-independent part of the bounds of the coupling's spectrum
// (pkpm_harm_node_spec).
// The w flux is alpha_k U + T(U), alpha_k = -K w - (b.grad xi) w^2 on every block and
//   T: F <- t.M0 + beta G, G <- t.N + beta H, M0 <- G t + R_B1 N, N <- H t + R_B1 O,
// with H = 4 theta G - 2 theta^2 F, O = 6 theta N - 6 theta^2 M0 and t = t0 + w t1.

// w flux: shared diagonal alpha_k
#define PKPM_HS_KW 0 // -K (w^1)
#define PKPM_HS_BG 1 // -b.grad xi (w^2)
// w flux: first-harmonic couplings T (a, c = 0, 1, 2; j the power of w)
#define PKPM_HS_TF(a, j) (2 + 2 * (a) + (j)) // F <- M0_a: t_a (j = 0, 1)
#define PKPM_HS_TG(a, j) (8 + 2 * (a) + (j)) // G <- N_a: t_a
#define PKPM_HS_TM(a, j) (14 + 2 * (a) + (j)) // M0_a <- G: t_a
#define PKPM_HS_TNG(a, j) (20 + 2 * (a) + (j)) // N_a <- G: 4 theta t_a
#define PKPM_HS_TNF(a, j) (26 + 2 * (a) + (j)) // N_a <- F: -2 theta^2 t_a
#define PKPM_HS_RMN(a, c) (32 + 3 * (a) + (c)) // M0_a <- N_c: (R_B1)_ac (w^0)
#define PKPM_HS_RNN(a, c) (41 + 3 * (a) + (c)) // N_a <- N_c: 6 theta (R_B1)_ac
#define PKPM_HS_RNM(a, c) (50 + 3 * (a) + (c)) // N_a <- M0_c: -6 theta^2 (R_B1)_ac
// w flux: the two-field couplings of F and G (w^0)
#define PKPM_HS_BF 59 // F <- G: beta = xi div b
#define PKPM_HS_BGG 60 // G <- G: 4 theta beta
#define PKPM_HS_BGF 61 // G <- F: -2 theta^2 beta
// the w-flux values are [0, 62)
#define PKPM_HS_W_END 62

// then theta, b (three components), alpha_0 and c
#define PKPM_HS_NODE_THETA 62
#define PKPM_HS_NODE_B 63
#define PKPM_HS_NODE_A0 66
#define PKPM_HS_NODE_CH 67
// then the w-independent part of the spectrum bounds of T (pkpm_harm_node_spec)
#define PKPM_HS_NODE_SLO 68 // smallest eigenvalue candidate
#define PKPM_HS_NODE_SHI 69 // largest eigenvalue candidate
#define PKPM_HS_NODE_ST 70 // sqrt((2 + sqrt 2) theta), the factor of |t|
#define PKPM_HS_NODE_NV 71

