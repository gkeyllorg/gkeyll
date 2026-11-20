# CBC kinetic-electron runs on the geometry-fixes branch: review of the crashes (Oct 22, 2026)

Scope: `/pscratch/sd/a/ah1032/gkeyll_adiabe/cbc/const_grad` on Perlmutter, Gkeyll branch
`get_em_back-flr-geom_fixes` (install `gkylsoft_em_flr_geom`, frames stamped `a446f1cd8-geomfix-1177`),
compared with the earlier EM+FLR install (`2682b2cef-emflr`, Oct 12) and the Sep 22 ES build.

## Summary

1. The beta = 0 crash on the geometry-fixes branch is not an electromagnetic bug. The pure electrostatic
   run (`geomfix/wk_ke_em_flr_es`: `GKYL_GK_COLLISIONLESS_ES`, beta_e = 0, no Ampere solve) dies at the
   same step and time as the `beta0` run: step ~2200, t = 2.398e-6 s = 0.63 R/c_s, dt collapse. The EM
   code path is never exercised ("Apar eqn solve: 0.00 %").
2. What changed relative to every earlier CBC run is the twist-shift (TS) connection: PR #1177 computes
   the shift from `mapc2p` and its sign is the negative of the shift functions all previous inputs
   provided (|S_lo| 1.39-2.72 m, difference 4.7e-4 m = 1.4e-3 Ly, slope within 2.6 %).
3. I verified independently that #1177's sign is the geometrically correct one:
   - the updater `gkyl_bc_twistshift` fills ghost(y) = opposite-skin(y - S) for both edges (numerical test
     with a Gaussian donor, constant and sheared shift);
   - with the CBC `mapc2p` (phi = y/Cy + alpha(r, theta), alpha(z_max) - alpha(z_min) = -2 pi q), the lower
     ghost at label y is the physical point of the upper skin at y + 2 pi q Cy, hence S_lo = -2 pi q(x) Cy,
     which is #1177's value. With it the effective k_x = k_x + k_y Cy q' theta is continuous through the cut;
     with the old input sign it jumps by 4 pi k_y Cy q'.
   So all earlier CBC runs (ES adiabatic, ES kinetic electrons, EM with the input shift) were run with a
   misconnected flux tube. In those runs the cut acts as an absorber (fluctuation amplitude in the end
   z cells is 1/7 of the interior at 5 R/c_s in `wk_ke_Ny32_const_grad`), which is why they were calm.
4. The crash is a cut-localized numerical mode of the correctly connected kinetic-electron system, and it
   is reproduced by the OLD code when the input shift sign is flipped (Oct 20 runs): same frame-0 and
   frame-1 fields to three digits, same growth rate (13-15 c_s/R), same k_y band. It is therefore not a
   regression introduced by #1136, #1133 or the main merges. The geometry used by the two builds is
   identical to 1e-12 in every coefficient the kernels read (interior and all z faces, cut included);
   #1136 is inactive under twist-shift because z is not flagged periodic.
5. The unstable band is the k_y range where the TS maps modes beyond the x Nyquist of the grid. The chart
   k_x jump per crossing is k_y |dS/dx|; with |dS/dx| = 2.8 (x_min), 5.1 (center), 7.6 (x_max) and
   dx = 2 rho_s this reaches pi/dx at m = 10.6, 5.9, 4.0 respectively (m = k_y Ly/2 pi, k_y rho_s = 0.052 m).
   Observed: geomfix ES, m = 5-7 at x cells 34-40 (k_y|dS/dx| dx/pi = 1.3-1.6); old code flipped input,
   m = 7-9 at x cells 5-24 (0.9-1.07); Ny = 8 runs (m <= 4, ratio <= 1.0) never showed it. The band
   m = 6-9 is k_y rho_s = 0.31-0.47, i.e. the ITG peak, so it cannot be dropped from the simulation.
   The anti-aliasing filter in use (upsample 2, half-width 1) attenuates but does not remove this content.

## Run inventory (const_grad)

| directory | build | TS sign | type | outcome |
| --- | --- | --- | --- | --- |
| wk_Ny32_const_grad, wk_Ny16/48/64, *_noflr, *_api* | Sep 22 (56d884cdbf2c) | input (old) | ES adiabatic | calm, 100 R/c_s; misconnected cut |
| wk_ke_Ny32_const_grad | Sep 22 | input (old) | ES kinetic e | calm, 100 R/c_s; misconnected cut |
| wk_ke_em_flr_beta0p1/wrong_ts_sign_run | emflr Oct 12 | input (old) | EM 0.1 % | completed 100 R/c_s after the phi double-TS fix |
| wk_ke_em_flr_beta{0p5,1p3}/crashed_wrong_ts_sign | emflr | input (old) | EM 0.5/1.3 % | cut mode, crash ~3 R/c_s (m = 1-2 at the cut) |
| wk_ke_em_flr_beta{0,0p1,0p5,1p3} (Oct 20) | emflr | input flipped | ES / EM | burst at 1.0-1.3 R/c_s; ES blows up (m = 9 at the cut) |
| geomfix/wk_ke_em_flr_beta{0p1,0p5,1p3} (Oct 20) | geomfix a446f1cd8 | internal (#1177) | EM | dt collapse within 300/200/100 steps |
| geomfix/wk_ke_em_flr_beta0, geomfix/wk_ke_em_flr_es (Oct 22) | geomfix | internal | ES | identical crash at 0.63 R/c_s, m = 5-7 at the cut |

## Evidence

- Spectra at the cut (|phi_m| rms over x, z_up cell): geomfix ES frame 1 (0.25 R/c_s) m = 6: 0.44 V vs
  interior 1e-2 V; frame 2: 18 V (gamma ~ 15 c_s/R); frame 3 everything > 1e6 V. Old code flipped ES:
  m = 9 0.13 -> 3.7 -> 90 V over frames 3-5 (gamma ~ 13 c_s/R). Calm main-sign run at 5 R/c_s: cut
  m = 6-9 at 6e-3 V, interior 5e-2 V.
- Frame 0 and frame 1 y-averaged phi(x, z) agree between geomfix and the old code with flipped input to
  3 digits in every z cell (including the end cells), so no geometry/field-solver difference is visible.
- Surface geometry files (`geo_surf2_*`) differ between the builds only in components 4-5 (and the
  corresponding padding of the vector arrays) at z faces 3 and 11: that is stale padding beyond the
  4-coefficient 2D surface basis, never read by the kernels. Cosmetic.
- Twist-shift updater accuracy (local build, `ts_crossing_test.c` in the session scratchpad): plain TS
  reproduces the exact shift of the DG donor to 1e-6 (quadrature limit) in every cell, including cells where
  S(x) crosses an integer multiple of Ly (shifted edges coincide with the periodic boundary, as in the CBC
  grid where Ly = 32 dy and S/Ly crosses 5, 6, 7, 8 at x cells 14/15, 27/28, 37/38, 45/46). With the
  production filter (upsample 2, half-width 1) the ghost deviates smoothly by 2-12 % (filter design), with
  no anomaly at the crossings. GPU kernels were not tested here.
- Resolution requirement for the correct connection: k_y,max |dS/dx|_max <= pi/dx needs dx <= 0.5 rho_s
  (Nx >= 194) for m <= 16, or dx <= 1.0 rho_s (Nx >= 97) for m <= 8 (k_y rho_s <= 0.42).

## Minor code notes

- `gyrokinetic/zero/fem_parproj.c:89`: `isdirichlet` is computed but never used.
- Written surface-geometry arrays carry uninitialised padding in components >= num_surf_basis.
- The internal shift from `deltats` differs from the analytic input shift by 1.4e-3 Ly (projection from
  Gauss nodes vs node values, as noted in the PR); irrelevant for the dynamics.

## Recommendations

1. Keep #1177's sign. Treat the earlier CBC validation (adiabatic and kinetic-electron ES) as obtained
   with a misconnected cut; re-validate once the connected setup is stable.
2. Tests to confirm the mechanism, cheapest first, all on the geomfix install, ES, 5 R/c_s:
   a. adiabatic electrons, Nx = 48 (does the correct connection alone misbehave without fast electrons?);
   b. kinetic electrons, Ny = 32, Nx = 96 and Nx = 144 (prediction: the m = 5-9 cut mode disappears or
      moves to m >= 12);
   c. kinetic electrons, Nx = 48 with `ts_upsample_factor = 4`, `ts_filter_half_width = 3-4`
      (prediction: delayed/damped if aliasing feeds it);
   d. production compromise: Ny = 16, Nx = 96 (k_y rho_s <= 0.42 well resolved, TS resolved everywhere).
3. Only then redo the EM beta scan; beta amplifies the same cut mode (crash in 100-300 steps at 0.1-1.3 %).
4. Code-side option: have the TS drop (not alias) content above pi/dx, or at least warn at start-up when
   k_y,max max|dS/dx| > pi/dx.
