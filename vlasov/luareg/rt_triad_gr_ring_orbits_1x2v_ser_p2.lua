-- Ring of bound orbits around a Schwarzschild black hole in ingoing Kerr-Schild coordinates
-- (equatorial r-phi preset geometry) through the GR triad Vlasov solver (1x2v, collisionless). f =
-- sqrt(h) n0 exp(-(E - E_c(L))/T) exp(-(L - L_c)^2 / 2 sigma_L^2), with sqrt(h) the metric
-- determinant that the triad distribution carries as a density in the tetrad momenta, E = -p_t and
-- L = p_phi the integrals of motion, E_c(L) the energy of the circular orbit with angular momentum
-- L and L_c that of the circular orbit at r = 10 M, T = 1e-3 (the epicyclic frequency there is 0.63
-- of the orbital one, so this energy spread already makes the ring two M wide) and sigma_L = 0.1,
-- is an exact stationary state, the relativistic analogue of the Newtonian ring test; its density
-- profile must be unchanged after t = 50 M, a quarter of the orbital period, with absorbing walls
-- at r = 5 M and 15 M, beyond the pericentres and apocentres of the ring's eccentric edge orbits.
-- Figures of merit (serendipity p2, 25 x 16 x 36 cells, 1117 steps): at t = 50 M the density
-- profile differs from its initial projection by 6.6e-3 in relative L1 (1.4e-2 at the peak); mass
-- and energy change by 7.7e-3 and 7.5e-3 as the eccentric edge orbits reach the absorbing walls.
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
mass_bh = 1.0 -- Black hole mass M.
spin_bh = 0.0 -- Black hole spin a = J/M.
n0 = 1.0 -- Amplitude of the distribution function.
r_c = 10.0 -- Radius of the central circular orbit.
T0 = 1.0e-3 -- Energy spread of the ring about the circular-orbit energy.
sigma_L = 0.1 -- Width of the Gaussian angular-momentum window.

-- Specific angular momentum of the circular orbit at r_c.
L_c = math.sqrt(mass_bh * r_c * r_c / (r_c - (3.0 * mass_bh)))

Nr = 25 -- Cell count (configuration space: radial direction).
Nvr = 16 -- Cell count (velocity space: radial momentum).
Nvphi = 36 -- Cell count (velocity space: azimuthal momentum).
r_min = 5.0 -- Inner radius of the domain (absorbing).
r_max = 15.0 -- Outer radius of the domain (absorbing).
pr_lo = 0.0 -- Lower boundary of the radial momentum domain (orbits have p_r = 0.18 to 0.28).
pr_hi = 0.5 -- Upper boundary of the radial momentum domain.
pphi_lo = 0.23 -- Lower boundary of the azimuthal momentum domain (L / r spans 0.27 to 0.7).
pphi_hi = 0.72 -- Upper boundary of the azimuthal momentum domain.
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 50.0 -- Final simulation time.
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
  lower = { r_min },
  upper = { r_max },
  cells = { Nr },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = {  }, -- Periodic directions.

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
          local r, p_r, p_phi = xn[1], xn[2], xn[3]
          local M = mass_bh

          -- Kerr-Schild Hamiltonian of the equatorial r-phi preset (a = 0), E = -p_t, and the
          -- angular momentum L = p_phi (covariant) = r p_phi (tetrad).
          local alpha = 1.0 / math.sqrt(1.0 + (2.0 * M / r))
          local E = alpha * (math.sqrt(1.0 + (p_r * p_r) + (p_phi * p_phi)) - ((2.0 * M / r) * p_r))
          local L = r * p_phi

          -- Energy of the stable circular orbit with this angular momentum, E_c(L) at
          -- r_c(L) = (L^2 / 2M) (1 + sqrt(1 - 12 M^2 / L^2)): the minimum energy of bound orbits
          -- with L.
          local disc = 1.0 - (12.0 * M * M / (L * L))
          if disc < 0.0 then
            return 0.0 -- no circular orbit below L = 2 sqrt(3) M: outside the ring
          end
          local rc = (L * L / (2.0 * M)) * (1.0 + math.sqrt(disc))
          local E_c = (1.0 - (2.0 * M / rc)) / math.sqrt(1.0 - (3.0 * M / rc))

          -- Ring F = n0 exp(-(E - E_c(L))/T0) exp(-(L - L_c)^2 / 2 sigma_L^2): a function of the
          -- integrals of motion only, hence an exact stationary state, peaked on circular orbits.
          -- The triad distribution is a density in the tetrad momenta, so it carries the metric
          -- determinant sqrt(h) = r sqrt(1 + 2M/r) of the r-phi slice (Liouville measure
          -- dx dp = sqrt(h) dx dp_hat).
          local dL = (L - L_c) / sigma_L
          local sqrt_h = r * math.sqrt(1.0 + (2.0 * M / r))
          return sqrt_h * n0 * math.exp(-(E - E_c) / T0) * math.exp(-0.5 * dL * dL)
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
