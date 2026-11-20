# EM twist-shift z-boundary instability (CBC, beta_e >= 0.5%)

Investigation of the dt collapse of the kinetic-electron EM+FLR Cyclone runs
(`/pscratch/sd/a/ah1032/gkeyll_adiabe/cbc/const_grad/wk_ke_em_flr_beta0p5`) that
remains after the phi double-twist-shift fix of Oct 16, 2026. Experiments live in
`.../const_grad/exp_zbnd/` (beta_e = 0.5%, Ny = 8, 48x8x16 x 8x4, one GPU each).

## Symptoms

- dt steady at 3.4e-9 s until t ~ 3 R/c_s, then collapse within ~600 steps.
- A single ky mode (m = 1, ky rho_s = 0.052) with kx ~ 4-6 (the ky*shat kx of one
  twist-shift pass), confined to the first and last z cells (decay x4 per cell inward),
  growing at ~6 c_s/R identically in phi, dApar/dt and n_e from ~1e-3 V.
- In the growth phase the boundary cells carry a sub-cell parallel slope of phi 65x the
  interior and electron u_par 5x, while the cell-averaged n_e is still at seed level;
  |E_par|/|d_z phi| ~ 0.76 (partial Alfvenic cancellation).

## Discriminating runs (all beta_e = 0.5% unless noted)

| run | outcome | conclusion |
| --- | --- | --- |
| base | unstable (boundary/interior phi 0.24 -> 460, dt collapse 4.3 R/c_s) | reproduces the full grid |
| beta_e = 0.1% | stable to 10.5 R/c_s | beta threshold between 0.1% and 0.5% |
| no shear / no twist-shift | stable, boundary ratio ~1 | needs the TS cut |
| TS filter off (-1) | unstable, same timing | not the filter/upsampling |
| cfl_frac_omegaH 0.5 | unstable, same timing | not a time-step issue |
| no initial ion flow | unstable, ~2x later | flow only accelerates it |
| FLR off | unstable, same timing | not the FLR/Ohm fixed-point code |
| Apar not z-smoothed each step | same linear growth, saturates | not the smoothing accumulation |
| `is_apar_static` (also sets dApar/dt = 0) | stable | needs dApar/dt |
| electrons ignore dApar/dt (Apar still evolves) | exponential mode gone | the kinetic dApar/dt term drives it |
| Apar frozen (Ohm's law and kinetic dApar/dt kept) | grows at ~half the rate, no dt collapse | flutter dynamics only enhance it |
| m_i/m_e = 25 | stable to 16 R/c_s | amplification ~ beta_e (m_i/m_e)/(k_perp^2 rho_s^2) |
| Ly = 60 rho_s (ky_min doubled) | stable to 6 R/c_s | idem (~7 unstable, <~2 stable) |

## Mechanism

The parallel FEM smoother of phi and Apar uses `GKYL_FEM_PARPROJ_DIRICHLET_GHOST`:
the end nodes are set from one-sided DG face values, the lower one to T_LU[phi_up(z_max)]
alone and the upper one to phi_up's own face value, whereas interior nodes are
least-squares averages of both adjacent cells. The lower skin cell is therefore pinned to
a foreign value and carries a spurious sub-cell d_z phi every step. Electrostatically this
is damped; with the inductive response it is amplified by the cancellation-type factor
beta_e (m_i/m_e)/(k_perp^2 rho_s^2), largest for the lowest ky, and becomes a growing mode
at the cut.

## Fix (local tree, `gyrokinetic/apps/gk_field_2x3x.c`)

`gk_field_2x3x_par_ghosts_ts_sym`: each parallel ghost gets the average of the
twist-shifted opposite skin and the boundary value of its own skin, so both end nodes are
the mean of the two face values meeting at the cut and ghost_L = T_LU(ghost_U). Used by
the 3x twist-shift phi and Apar projections without LCFS (the IWL branch keeps the old
one-sided treatment for now). Built in `installs/gkeyll_em_flr_exp3`; the production
install is untouched.

Verification: `symts` stable to 8.3 R/c_s with boundary ratio ~0.45 and the boundary
d_z phi equal to the interior; `symts_nofilt` stable; `symts_b0p1` vs `b0p1` field energy
ratio 0.95-1.03 over 2-25 us.

Side effects: the tiny 3x twist-shift regression cases (`rt_gk_cbc_3x2v_p1`,
`rt_gk_cbc_passive_3x2v_p1`, `rt_gk_passive_3x2v_p1`, `rt_gk_tcv_core_3x2v_p1`) change by
30-66% because they are initial-transient tests whose t = 0 charge density is not
field-line continuous; their reference outputs need regenerating. Also fixed the latent
precedence bug `fem_parproj.c:89` (`isdirichlet` is unused).

## Oct 19 update: the end-node fixes are only partial

- Full grid (48x32x16) beta_e = 1.3% with the symmetric averaging (branch f184bcc6e): crashed at
  1.7 R/c_s with the same cut mode (one-sided treatment: 3 R/c_s).
- Ny = 8, beta_e = 1.35%: one-sided crashes at 1.0e-5 s; averaging delays the mode (boundary ratio
  52 at 2.1e-5 s); a second variant - both z ghosts twist-shifted, phi/Apar projected on the
  ghost-extended range with free ends (`gk_field_2x3x_par_ghosts_ts_both`, local tree) - delays it a
  bit more (7 at 2.0e-5 s) but it still grows. Both variants are stable at beta_e = 0.5%.
- Ohm's law and Poisson use `fem_poisson_perp` (one 2D solve per z cell, z kept DG), so there is no
  z-node sharing in those solves; the deflated solver only serves the FLR operators (excluded).
- Remaining suspect: with least-squares end nodes the two ends are twist-shift images of each other
  only to within the DG shift-interpolation error; the EM coupling amplifies that residual. Tests of
  the interpolation accuracy (upsample 8, plain TS) at beta_e = 1.35% are queued (exp_zbnd/ext_b1p3_*).

## Oct 19, later: root cause = twist-shift sign (PR #1177)

PR #1177 (`tsbc_shift_sign_check`) states that the twist-shift sign used on `main` (and by the
input shift functions of the CBC inputs: lower +Cy[alpha(z_min)-alpha(z_max)], upper the opposite)
is inconsistent with the mapping: a blob advected through theta = +-pi is not continuous. Test with
the unchanged one-sided production library at beta_e = 1.35%, Ny = 8, with only the two input shift
functions' signs flipped (`exp_zbnd/flipsign_b1p3`): stable to 5.7 R/c_s with a flat boundary ratio
(0.5-0.9), where the same code with the inputs' sign crashed at 1.0e-5 s (ratio 441 by 7.6e-6 s).

Interpretation: the cut connects each (x, y) to the wrong field line; electrostatic ITG tolerates the
misconnection, the Alfvenic E_par dynamics does not, and the error is amplified by
beta_e (m_i/m_e)/(k_perp^2 rho_s^2). The end-node variants above (averaging, extended projection,
pinned) only delayed the symptom. A pinned-projection restart (`exp_zbnd/pin_b1p3`) and the full-grid
beta_e = 1.3% run with the flipped sign (`wk_ke_em_flr_beta1p3_flipsign`, interactive) are running.

Consequences: all CBC runs so far (ES and EM) used the misconnecting sign; growth rates and fluxes
should be rechecked once #1177 (internal shift from mapc2p, corrected sign) is merged. #1136
(one-sided surface geometry at the periodic z boundary) is a second cut inconsistency worth taking.

## Test plan

1. Rebuild the production install with the fix; rerun beta_e = 0.5% and 1.3% on the
   48x32x16 grid (expected: no boundary mode, dt steady).
2. Production-resolution ES CBC (`wk_ke_const_grad` setup) with and without the fix:
   growth rates and saturated fluxes should agree within noise.
3. Regenerate the 3x twist-shift regression references; keep the Rosenbluth-Hinton and
   2x cases (untouched code paths) as controls.
4. Apply the same symmetric end treatment in the IWL branch (TCV EM runs) and repeat 1-2.

## Update Oct 20: full-grid check of the sign flip

The Ny=32 production-grid β=1.35% run with the flipped (PR #1177) shift sign and the one-sided
TS ghosts (`wk_ke_em_flr_beta1p3_flipsign`) also failed, at t = 1.25e-5 s (3.3 R/c_s, step 8600),
essentially where the old-sign run failed (1.15e-5 s). So the shift sign is not the whole story at
full resolution; the Ny=8 test (stable to 5.7 R/c_s) was not representative.

The failure differs in character from the old-sign run, though:

| | old sign | flipped sign |
| --- | --- | --- |
| φ boundary/interior energy, frames 1/2/3 | 0.58 / 177 / 445 (monotonic growth) | 90 / 1.5 / 0.53 (burst, then decays) |
| dt history | smooth decline to collapse | bounces 1.6e-9 ↔ 0.7e-9 from step 3900, then collapses within 100 steps |
| crash frame | φ, A∥, ∂tA∥, n_e all 1e2–1e3× interior | only ∂tA∥ (1.3e3, TS jump 3.9e3); φ, A∥ ratios ~1 |

Interpretation: with the correct connection the clean exponential cut mode is gone, but an
intermittent (bursty) boundary excitation remains, and the final event is a sudden ∂tA∥
blow-up at the cut rather than slow growth. Candidates: the end-node treatment of the parallel
projection (symptom treatment earlier delayed but never removed growth), or a time-step /
ωH-bound issue when the burst locally raises k∥ (dt already bouncing). Next data: the 5 R/c_s
reruns of all four β values with the corrected sign (jobs 60033782/83/89/90).

## Update Oct 20 (later): the input sign flip is retracted

5 R/c_s reruns of all four cases (ES β=0, β=0.1%, 0.5%, 1.3%) with the flipped input shift sign
on the production grid (Ny=32, 2 GPUs each, z split in 2):

| case | result with flipped sign | same case with main's sign |
| --- | --- | --- |
| ES β=0 | blows up at 1.05–1.3 R/c_s: field energy 1e-2 → 11.6 (γ ≈ 43 c_s/R), max in the last z cell, n_e boundary/interior 33 with TS jump 3.6 | calm for 47 R/c_s |
| β=0.1% | burst at 1.05–1.3 R/c_s (φ boundary/interior 8), saturates at field energy ~34 | field energy ~1.5e-2 at 2.6 R/c_s, no burst |
| β=0.5% | burst (ratio 27 at 0.75 R/c_s), then ratio ~1, dt bouncy, alive at 3.75 R/c_s | crashed at ~3.2 R/c_s (cut mode) |
| β=1.3% | burst (ratio 12), alive at 3.8 R/c_s with z split in 2; crashed at 3.3 R/c_s with z split in 4 | crashed at 3 R/c_s (cut mode) |

The burst spectrum at the cut is dominated by m = 1, 2 and 6 (k_y ρ_s ≈ 0.3). The Ny=8 test that
motivated the flip resolves only m ≤ 4 and never showed the burst (its field energy decays
smoothly 47 → 2.8), so it was a false positive.

Conclusion: with the current `bc_twistshift`, main's input sign is the self-consistent connection
(the ES case is the clean control). The flipped sign misconnects the cut and is violently unstable
even electrostatically. The EM-only cut instability with main's sign is still unexplained. The
convention change of PR #1177 can only be tested with that PR's code, i.e. the merged install.

## Update Oct 20 (evening): the correct connection is the unstable one

Rerun of the four cases on the geometry-fixes branch (`get_em_back-flr-geom_fixes` a446f1cd8: main +
PR #1136 + PR #1177, shift computed internally from mapc2p). All four die at the cut almost at once
(β=1.3% within 100 steps, β=0.5% 200, β=0.1% 300 steps, ES β=0 at 0.6 R/c_s), boundary/interior
field ratios 60–840. The internal shift is exactly the negative of main's input shift (same magnitude).

Geometry check, independent of any plotting tool. The twist-shift updater fills
ghost(y) = donor(y − S(x)) (`ts_find_intersect`: yTar − yShift(x) = yDo, donor = periodic copy of the
opposite skin). The input's mapc2p is φ = y/Cy + α(r, θ) with b̂ ∥ e_3, so a field line at label y
leaving z_up re-enters z_lo at y + Cy[α(z_max) − α(z_min)]; hence lower ghost(y) = upper skin(y − Cy Δα)
and S_lo = Cy[α(z_max) − α(z_min)], i.e. PR #1177's sign. PR #1177 is geometrically right; main's
inputs, and therefore every CBC run so far, misconnect the cut by a y-reflection of the shift.

Consequence: the misconnection scrambles the phase of every mode across the cut (2 k_y S mod 2π,
varying with x through the shear), which acts as damping. With the correct closure a numerical
amplifier at the cut can feed back coherently, and it does: kinetic-electron ES blows up with
γ ≈ 43 c_s/R at Ny=32 (m = 1, 2, 6 at the cut), EM faster with β. Ny=8 (m ≤ 4) never shows it, so
the amplified modes are near the grid scale. Suspects: the twist-shift interpolation with its
anti-aliasing filter and upsampling, or the parallel scheme at the cut (end-node projection for EM).
Proposed tests on the geometry-fixes install (ES β=0, Ny=32): (a) filter off, (b) adiabatic electrons,
(c) Ny=16.
