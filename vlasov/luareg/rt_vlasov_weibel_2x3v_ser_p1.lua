-- Weibel instability of two counter-streaming electron beams (non-relativistic Vlasov-Maxwell).
-- 2x3v: the oblique mode at 45 degrees (k0 = 1) is seeded through a density perturbation with the
-- electric field from Gauss's law in the ratio E_y/E_x of the growing eigenmode, as in the 2x2v
-- input, with a passive third velocity dimension. Ions are a neutralizing background.
-- Growth rate of the magnetic energy from linear theory: gamma = 0.1828 (measured 0.173 over
-- t = 30-60 with four serendipity p = 1 cells per wavelength); the mode saturates at t = 80.

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

n0 = 1.0 -- Reference number density.
uy_drift = 0.3 -- Drift velocity of the beams (y-direction).
T_elc = 0.01 -- Electron temperature.

theta = (45.0 / 180.0) * pi -- Perturbation angle.
k0 = 1.0 -- Perturbed wave number.
alpha = 1.6202 -- Ratio E_y / E_x of the growing eigenmode (linear theory).
perturb_n = 2.0e-6 -- Perturbation density.

-- Derived physical quantities (using normalized code units).
vte = math.sqrt(T_elc / mass_elc) -- Electron thermal velocity.
kx = k0 * math.cos(theta) -- Perturbed wave number (x-direction).
ky = k0 * math.sin(theta) -- Perturbed wave number (y-direction).

-- Simulation parameters.
Nx = 4 -- Cell count (configuration space: x-direction).
Ny = 4 -- Cell count (configuration space: y-direction).
Nvx = 12 -- Cell count (velocity space: vx-direction).
Nvy = 12 -- Cell count (velocity space: vy-direction).
Nvz = 4 -- Cell count (velocity space: vz-direction).
Lx = 2.0 * pi / kx -- Domain size (configuration space: x-direction).
Ly = 2.0 * pi / ky -- Domain size (configuration space: y-direction).
vx_max = 1.0 -- Domain boundary (velocity space: vx-direction).
vy_max = 1.0 -- Domain boundary (velocity space: vy-direction).
vz_max = 1.0 -- Domain boundary (velocity space: vz-direction).
poly_order = 1 -- Polynomial order.
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
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x- and y-directions).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1, 2 }, -- Periodic directions (x- and y-directions).

  -- Electrons.
  elc = Vlasov.Species.new {
    modelID = G0.Model.Default,
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
          local x, y = xn[1], xn[2]

          local n = 0.5 * (1.0 + perturb_n * math.cos((kx * x) + (ky * y))) * n0 -- Beam number density.
          return n
        end,
        temperatureInit = function (t, xn)
          return T_elc -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return 0.0, uy_drift, 0.0 -- Drift velocity of the first beam.
        end,

        correctAllMoments = true,
        useLastConverged = true
      },
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x, y = xn[1], xn[2]

          local n = 0.5 * (1.0 + perturb_n * math.cos((kx * x) + (ky * y))) * n0 -- Beam number density.
          return n
        end,
        temperatureInit = function (t, xn)
          return T_elc -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return 0.0, -uy_drift, 0.0 -- Drift velocity of the second beam.
        end,

        correctAllMoments = true,
        useLastConverged = true
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
      local x, y = xn[1], xn[2]

      local Ex = -perturb_n * math.sin((kx * x) + (ky * y)) / (kx + (ky * alpha)) -- Total electric field (x-direction).
      local Ey = alpha * Ex -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = (kx * Ey) - (ky * Ex) -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()
