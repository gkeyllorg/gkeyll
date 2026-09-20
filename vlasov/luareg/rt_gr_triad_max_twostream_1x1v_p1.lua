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

-- Derived physical quantities (using normalized code units).
gamma = 1.0 / math.sqrt(1.0 - (Vx_drift * Vx_drift)) -- Gamma factor.
Vx_drift_SR = gamma * Vx_drift -- Relativistic drift velocity (x-direction).

alpha = 1.0e-8 -- Applied perturbation amplitude.
kx = 0.5 -- Perturbed wave number (x-direction).

-- Simulation parameters.
Nx = 64 -- Cell count (configuration space: x-direction).
Nvx = 64 -- Cell count (velocity space: vx-direction).
Lx = 2.0 * pi / kx -- Domain size (configuration space: x-direction).
vx_max = 8.0 -- Domain boundary (velocity space: vx-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 500.0 -- Final simulation time.
num_frames = 100 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

vlasovApp = Vlasov.App.new {

  -- Flat preset geometry makes the GR triad/GR-Maxwell system reduce to the
  -- standard Cartesian two-stream setup.
  geom = Vlasov.Geom.new {
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.Flat,
    massBH = 0.0,
    spinBH = 0.0,
  },

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
  decompCuts = { 1 }, -- Cuts in each coordinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic direction (x-direction only).

  -- Electrons.
  elc = Vlasov.Species.new {
    modelID = G0.Model.TriadGR,
    charge = charge_elc, mass = mass_elc,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 2,
    projections = {
      -- Two counter-streaming relativistic LTE distributions.
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]

          local metric_det = 1.0
          local n = 0.5 * (1.0 + alpha * math.cos(kx * x)) * n0
          return metric_det * n
        end,
        temperatureInit = function (t, xn)
          return T -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift_SR -- Total left-going relativistic drift velocity.
        end,

        correctAllMoments = true,
        useLastConverged = true
      },
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]

          local metric_det = 1.0
          local n = 0.5 * (1.0 + alpha * math.cos(kx * x)) * n0
          return metric_det * n
        end,
        temperatureInit = function (t, xn)
          return T -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return -Vx_drift_SR -- Total right-going relativistic drift velocity.
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

    -- Use GR-Maxwell field model in the flat limit.
    fieldID = G0.FieldModel.GR,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Dx = -alpha * gamma * math.sin(kx * x) / kx -- Electric displacement (x-direction).
      local Dy = 0.0 -- Electric displacement (y-direction).
      local Dz = 0.0 -- Electric displacement (z-direction).

      local Bx = 0.0 -- Magnetic field (x-direction).
      local By = 0.0 -- Magnetic field (y-direction).
      local Bz = 0.0 -- Magnetic field (z-direction).

      local metric_det = 1.0

      return metric_det * Dx, metric_det * Dy, metric_det * Dz,
        metric_det * Bx, metric_det * By, metric_det * Bz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()
