-- Electrostatic shock problem with kinetic ions and a DG isothermal Euler fluid for the electrons
-- (1x1v), the variant of rt_vlasov_hybrid_es_shock_1x1v with isothermal electrons (p = n Te at the
-- inflow temperature, Roe flux, vertex positivity limiter of the density) and the hybrid twin of
-- vlasov/luareg/rt_vlasov_es_shock_1x1v_tensor_p1.lua. Half domain: a plasma flowing toward a
-- reflecting wall at x = 0 (Mach 2 relative to c_s = sqrt(Te / m_i), Te / Ti = 4, m_i / m_e = 100)
-- piles up against its own reflection; the inflow boundary at x = Lx holds the upstream state
-- fixed. The electrons are coupled to the field by the implicit source step of the staggered
-- operator split, in which the ions deposit their current as a known current.
--
-- Figures of merit at t = 600 / omega_pe (tensor p1, hybrid basis of the ions, 192 x 32 cells, 9000
-- steps, omega_pe dt = 0.067, the CFL limit of the light wave of the field), measured as in
-- rt_vlasov_hybrid_es_shock_1x1v, with those of the adiabatic electron fluid of that test and of
-- the kinetic test in parentheses. The front is at x = 91.8 lambda_D and moves at 1.54 c_s between
-- t = 300 and t = 600 (25.2 lambda_D and 0.43 c_s; 66.7 lambda_D and 1.11 c_s). The ion density is
-- 1.74 n0 behind and 1.56 n0 ahead of it (2.62 and 1.77; 2.20 and 1.27). The potential jump across
-- it is e dphi / Te = 0.11 (1.09; 3.20), the logarithm of the density ratio, the momentum balance
-- of a massless isothermal fluid. Ions moving away from the wall ahead of the front: 0.54 n0 at
-- 2.42 c_s (0.67 n0 at 2.15 c_s; 0.24 n0 at 3.46 c_s). Neither electron fluid gives the potential
-- jump of the kinetic test.
local Vlasov = G0.Vlasov
local IsoEuler = G0.Vlasov.Eq.IsoEuler

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 100.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
Te_over_Ti = 4.0 -- Electron temperature / ion temperature.
mach_num = 2.0 -- Mach number of the inflow (relative to the ion sound speed).

-- Derived physical quantities (using normalized code units).
Te = 1.0 * charge_ion -- Electron temperature.
Ti = Te / Te_over_Ti -- Ion temperature.

vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.
cs = math.sqrt(Te / mass_ion) -- Ion sound speed.
Vx_drift = -mach_num * cs -- Drift velocity of the inflow (x-direction).

lambda_D = math.sqrt(epsilon0 * Te / (n0 * charge_ion * charge_ion)) -- Electron Debye length.
omega_pe = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)) -- Electron plasma frequency.

mu0 = 1.0 / 25.0 -- Permeability of free space (speed of light = 5 vte).

-- Simulation parameters.
Nx = 192 -- Cell count (configuration space: x-direction).
Nvx_ion = 32 -- Cell count (ion velocity space: vx-direction).
Lx = 192.0 * lambda_D -- Domain size (configuration space: x-direction).
vx_max_ion = 6.0 * cs -- Domain boundary (ion velocity space: vx-direction).
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set (hybrid phase-space basis of the kinetic species).
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 600.0 / omega_pe -- Final simulation time.
num_frames = 2 -- Number of output frames.
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
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Electrons (fluid): the closure is the equation object (isothermal Euler at the inflow
  -- temperature, Roe fluxes by default). The vertex positivity limiter of the density is on by
  -- default.
  elc = Vlasov.FluidSpecies.new {
    charge = charge_elc, mass = mass_elc,
    equation = IsoEuler.new { vThermal = math.sqrt(Te / mass_elc) },

    init = function (t, xn)
      local rho = mass_elc * n0 -- Electron mass density.
      local ux = Vx_drift -- Electron velocity (x-direction).

      -- Mass density and momentum density.
      return rho, rho * ux, 0.0, 0.0
    end,

    bcx = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcFixedFunc }
  },

  -- Ions (kinetic).
  ion = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,

    -- Velocity space grid.
    lower = { -vx_max_ion },
    upper = { vx_max_ion },
    cells = { Nvx_ion },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return n0 -- Ion total number density.
        end,
        temperatureInit = function (t, xn)
          return Ti -- Ion isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift -- Ion drift velocity.
        end
      }
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcReflect
      },
      upper = {
        type = G0.SpeciesBc.bcFixedFunc
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0,

    bcx = { G0.FieldBc.bcSymWall, G0.FieldBc.bcCopy }
  }
}

vlasovApp:run()
