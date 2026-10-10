-- Circular geodesic of a phase-space blob around a Schwarzschild black hole (a = 0) in ingoing
-- Kerr-Schild coordinates (equatorial r-phi preset geometry) through the GR triad Vlasov solver
-- (2x2v, collisionless). The blob starts at r0 = 8 M with the tetrad momenta of the prograde
-- circular orbit (Boyer-Lindquist E and L converted to the Kerr-Schild covariant momenta, which
-- have a radial component because g_tr is nonzero), so its density centroid must stay at r0 and
-- advance in phi at the orbital frequency Omega = sqrt(M) / (r0^(3/2) + a sqrt(M)) over t = 25 M;
-- the Kerr and Schwarzschild frequencies differ by 3.5 percent. Figures of merit
-- (tensor p1, 12 x 24 x 6 x 6 cells): at t = 25 M the density centroid is at (r, phi) =
-- (8.024, 4.259) against the exact (8, 4.246), i.e. within 0.07 cells in r and 0.05 cells in phi;
-- mass changes by 5.7e-4.
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
mass_bh = 1.0 -- Black hole mass M.
spin_bh = 0.0 -- Black hole spin a = J/M.
r0 = 8.0 -- Orbit radius.
phi0 = math.pi -- Initial azimuth of the blob.

v_c = math.sqrt(mass_bh / r0) -- Orbital velocity parameter sqrt(M / r0).
-- Specific energy of the prograde circular orbit (Boyer-Lindquist).
E_c = (1.0 - (2.0 * v_c * v_c) + (spin_bh * v_c * v_c * v_c)) /
  math.sqrt(1.0 - (3.0 * v_c * v_c) + (2.0 * spin_bh * v_c * v_c * v_c))
-- Specific angular momentum of the circular orbit.
L_c = math.sqrt(mass_bh * r0) *
  (1.0 - (2.0 * spin_bh * v_c * v_c * v_c) + (spin_bh * spin_bh * v_c * v_c * v_c * v_c)) /
  math.sqrt(1.0 - (3.0 * v_c * v_c) + (2.0 * spin_bh * v_c * v_c * v_c))
g_tt = -(1.0 - (2.0 * mass_bh / r0)) -- Boyer-Lindquist metric component g_tt at r0 (equator).
g_tp = -2.0 * mass_bh * spin_bh / r0 -- Boyer-Lindquist metric component g_tphi at r0.
-- Boyer-Lindquist metric component g_phiphi at r0.
g_pp = (r0 * r0) + (spin_bh * spin_bh) + (2.0 * mass_bh * spin_bh * spin_bh / r0)
det_g = (g_tt * g_pp) - (g_tp * g_tp) -- Determinant of the (t, phi) metric block.
p_t_up = ((g_pp * (-E_c)) - (g_tp * L_c)) / det_g -- Contravariant p^t of the orbit.
p_phi_up = ((-g_tp * (-E_c)) + (g_tt * L_c)) / det_g -- Contravariant p^phi of the orbit.
h_rr = 1.0 + (2.0 * mass_bh / r0) -- Kerr-Schild spatial metric h_rr at r0.
h_rp = -spin_bh * (1.0 + (2.0 * mass_bh / r0)) -- Kerr-Schild spatial metric h_rphi at r0.
-- Kerr-Schild spatial metric h_phiphi at r0.
h_pp = (r0 * r0) + (spin_bh * spin_bh * (1.0 + (2.0 * mass_bh / r0)))
-- Covariant Kerr-Schild p_r of the orbit (nonzero since g_tr != 0).
p_r_cov = ((2.0 * mass_bh / r0) * p_t_up) + (h_rp * p_phi_up)
pr0 = p_r_cov / math.sqrt(h_rr) -- Tetrad momentum p_r of the orbit.
-- Tetrad momentum p_phi of the orbit.
pphi0 = (L_c - ((h_rp / math.sqrt(h_rr)) * pr0)) / math.sqrt(h_pp - (h_rp * h_rp / h_rr))
-- Orbital angular frequency dphi/dt of the circular orbit.
omega_orb = math.sqrt(mass_bh) / ((r0 * math.sqrt(r0)) + (spin_bh * math.sqrt(mass_bh)))

Nr = 12 -- Cell count (configuration space: radial direction).
Nphi = 24 -- Cell count (configuration space: azimuthal direction).
Nvr = 6 -- Cell count (velocity space: radial momentum).
Nvphi = 6 -- Cell count (velocity space: azimuthal momentum).
r_min = 6.0 -- Inner radius of the domain.
r_max = 10.0 -- Outer radius of the domain.
pr_lo = 0.1 -- Lower boundary of the radial momentum domain.
pr_hi = 0.5 -- Upper boundary of the radial momentum domain.
pphi_lo = 0.3 -- Lower boundary of the azimuthal momentum domain.
pphi_hi = 0.6 -- Upper boundary of the azimuthal momentum domain.
sigma_r = (r_max - r_min) / (2.0 * Nr) -- Blob width in r (half a cell).
sigma_phi = 2.0 * math.pi / (2.0 * Nphi) -- Blob width in phi (half a cell).
sigma_pr = (pr_hi - pr_lo) / (4.0 * Nvr) -- Blob width in p_r (a quarter cell).
sigma_pphi = (pphi_hi - pphi_lo) / (4.0 * Nvphi) -- Blob width in p_phi (a quarter cell).
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 25.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { r_min, 0.0 },
  upper = { r_max, 2.0 * math.pi },
  cells = { Nr, Nphi },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = { 2 }, -- Periodic directions.

  -- Preset geometry.
  geom = Vlasov.Geom.new {
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.GR_KS_rphi,
    massBH = mass_bh,
    spinBH = spin_bh,
  },

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.TriadGR,
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { pr_lo, pphi_lo },
    upper = { pr_hi, pphi_hi },
    cells = { Nvr, Nvphi },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,
        init = function (t, xn)
          local r, phi, p_r, p_phi = xn[1], xn[2], xn[3], xn[4]

          -- Gaussian blob in phase space centred on the circular orbit (r0, phi0, pr0, pphi0).
          local dr, dphi = (r - r0) / sigma_r, (phi - phi0) / sigma_phi
          local dpr, dpphi = (p_r - pr0) / sigma_pr, (p_phi - pphi0) / sigma_pphi
          return math.exp(-0.5 * ((dr * dr) + (dphi * dphi) + (dpr * dpr) + (dpphi * dpphi)))
        end
      }
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcAbsorb
      },
      upper = {
        type = G0.SpeciesBc.bcAbsorb
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1 }
  },

  skipField = true,
}

vlasovApp:run()
