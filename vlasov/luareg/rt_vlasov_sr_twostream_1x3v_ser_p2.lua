-- Relativistic two-stream instability with the Vlasov-Maxwell system (1x3v).
-- Two counter-streaming Maxwell-Juttner electron beams (drift +-0.9 c, T = 0.04 mc^2, n0/2 each) on
-- a neutralizing ion background are perturbed in density at k = 0.3 omega_pe/c with amplitude 1e-5;
-- the electric field follows from Gauss's law. The 3V beams keep py and pz on 4 cells each.
-- Kinetic linear theory for the 1D Maxwell-Juttner beams at k = 0.3: a purely growing mode with
-- gamma = 0.1561 omega_pe, omega_pe = sqrt(n0 e^2 / (epsilon0 m)) = 1 with n0 the rest-frame
-- density.
-- Figures of merit (serendipity p2, 8 x 16 x 4 x 4 cells): the field energy grows at gamma = 0.1478
-- over t = 30-60 and saturates at t = 76 with energy 2.79 at t = 80; electron number conserved to
-- 1e-14.

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

n0 = 1.0 -- Reference number density.
T = 0.04 -- Temperature (units of mc^2).
Vx_drift = 0.9 -- Drift velocity (x-direction).

alpha = 1.0e-5 -- Applied perturbation amplitude.
kx = 0.3 -- Perturbed wave number (x-direction); close to the fastest-growing mode.

-- Derived physical quantities (using normalized code units).
gamma = 1.0 / math.sqrt(1.0 - (Vx_drift * Vx_drift)) -- Gamma factor.
Vx_drift_SR = gamma * Vx_drift -- Relativistic drift velocity (x-direction).

-- Simulation parameters.
Nx = 8 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: vx-direction).
Nvy = 4 -- Cell count (velocity space: vy-direction).
Nvz = 4 -- Cell count (velocity space: vz-direction).
Lx = 2.0 * pi / kx -- Domain size (configuration space: x-direction).
vx_max = 5.0 -- Domain boundary (velocity space: vx-direction).
vy_max = 1.2 -- Domain boundary (velocity space: vy-direction).
vz_max = 1.2 -- Domain boundary (velocity space: vz-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 80.0 -- Final simulation time.
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
  lower = { -0.5 * Lx },
  upper = { 0.5 * Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Electrons.
  elc = Vlasov.Species.new {
    modelID = G0.Model.SR,
    charge = charge_elc, mass = mass_elc,
    
    -- Velocity space grid.
    lower = { -vx_max, -vy_max, -vz_max },
    upper = { vx_max, vy_max, vz_max },
    cells = { Nvx, Nvy, Nvz },

    -- Initial conditions.
    numInit = 2,
    projections = {
      -- Two counter-streaming Maxwellians.
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]

          local n = 0.5 * (1.0 + alpha * math.cos(kx * x)) * n0 -- Total number density.
          return n
        end,
        temperatureInit = function (t, xn)
          return T -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift_SR, 0.0, 0.0 -- Total left-going relativistic drift velocity.
        end,

        correctAllMoments = true,
        useLastConverged = true
      },
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]

          local n = 0.5 * (1.0 + alpha * math.cos(kx * x)) * n0 -- Total number density.
          return n
        end,
        temperatureInit = function (t, xn)
          return T -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return -Vx_drift_SR, 0.0, 0.0 -- Total right-going relativistic drift velocity.
        end,

        correctAllMoments = true,
        useLastConverged = true
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1 }
  },

  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Ex = -alpha * gamma * math.sin(kx * x) / kx -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,
    
    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()