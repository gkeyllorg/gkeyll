-- Generalized Brio-Wu Riemann problem with two DG 5-moment (Euler) fluids and the Maxwell field of
-- the Vlasov app (1x), the counterpart of moments/creg/rt_5m_riem.c. The initial conditions follow
-- entry JE4 of Ammar's Simulation Journal
-- (https://ammar-hakim.org/sj/je/je4/je4-twofluid-shock.html), adapted from Section 7.1 of A.
-- Hakim, J. Loverich and U. Shumlak (2006), "A high resolution wave propagation scheme for ideal
-- Two-Fluid plasma equations", Journal of Computational Physics 219 (1): 418-442: m_i / m_e =
-- 1836.2, c = 1, and the ion skin depth equal to the domain length. The fluid-field coupling is the
-- staggered operator split with the implicit fluid-Maxwell source step.
--
-- The electrons become supersonic in the current sheet: they reach Mach 2.7 at t = 0.09, with their
-- internal energy down to 20% of their total energy (cell averages). The fluids use the HLLC flux
-- of the Riemann solver (rp_type) with the vertex positivity limiter, and the field its slope
-- limiter (limit_em).
--
-- Figures of merit at t = 10 (256 cells, serendipity p1, 8534 steps): relative L1 differences of
-- the cell averages from the moments finite-volume solution on 1024 cells, 1.1% (ion density), 1.2%
-- of the maximum (ion u_x), 1.6% (ion pressure), 1.5% (electron pressure), 1.2% and 1.6% of the
-- maximum (B_y, B_z). The electron pressure lies between 4.0e-6 and 5.6e-5 and the ion density
-- between 0.118 and 1.07.

local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_ion = 1.0 -- Proton mass.
charge_ion = 1.0 -- Proton charge.
mass_elc = 1.0 / 1836.2 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

rhol_ion = 1.0 -- Left ion mass density.
rhor_ion = 0.125 -- Right ion mass density.
pl = 5.0e-5 -- Left electron/ion pressure.
pr = 5.0e-6 -- Right electron/ion pressure.

Bx = 0.75e-2 -- Total magnetic field (x-direction).
Bzl = 1.0e-2 -- Left total magnetic field (z-direction).
Bzr = -1.0e-2 -- Right total magnetic field (z-direction).

-- Derived physical quantities (using normalized code units).
rhol_elc = rhol_ion * mass_elc / mass_ion -- Left electron mass density.
rhor_elc = rhor_ion * mass_elc / mass_ion -- Right electron mass density.

-- Simulation parameters.
Nx = 256 -- Cell count (configuration space: x-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 10.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = 100 -- Number of times to calculate field energy.
integrated_mom_calcs = 100 -- Number of times to calculate integrated moments.
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
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions.

  elc = Vlasov.FluidSpecies.new {
    charge = charge_elc, mass = mass_elc,
    equation = Euler.new { gasGamma = gas_gamma },

    init = function (t, xn)
      local x = xn[1]

      local rho, p
      if x < 0.5 then
        rho, p = rhol_elc, pl -- Left state.
      else
        rho, p = rhor_elc, pr -- Right state.
      end

      return rho, 0.0, 0.0, 0.0, p / (gas_gamma - 1.0)
    end,

    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
  },

  ion = Vlasov.FluidSpecies.new {
    charge = charge_ion, mass = mass_ion,
    equation = Euler.new { gasGamma = gas_gamma },

    init = function (t, xn)
      local x = xn[1]

      local rho, p
      if x < 0.5 then
        rho, p = rhol_ion, pl -- Left state.
      else
        rho, p = rhor_ion, pr -- Right state.
      end

      return rho, 0.0, 0.0, 0.0, p / (gas_gamma - 1.0)
    end,

    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
  },

  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,
    limitField = true,

    init = function (t, xn)
      local x = xn[1]

      local Bz = x < 0.5 and Bzl or Bzr

      return 0.0, 0.0, 0.0, Bx, 0.0, Bz, 0.0, 0.0
    end,

    bcx = { G0.FieldBc.bcCopy, G0.FieldBc.bcCopy },
  },
}

vlasovApp:run()
