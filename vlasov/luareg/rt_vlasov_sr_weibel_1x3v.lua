local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

n_elc1 = 0.5 -- First electron number density.
n_elc2 = 0.5 -- Second electron number density.
ux_elc1 = 0.0 -- First electron velocity (x-direction).
ux_elc2 = 0.0 -- Second electron velocity (x-direction).
uy_elc1 = 0.9 -- First electron velocity (y-direction).
uy_elc2 = -0.9 -- Second electron velocity (y-direction).
uz_elc1 = 0.0 -- First electron velocity (z-direction).
uz_elc2 = 0.0 -- Second electron velocity (z-direction).
T_elc1 = 0.04 -- First electron temperature (units of mc^2).
T_elc2 = 0.04 -- Second electron temperature (units of mc^2).

alpha = 1.0e-3 -- Applied perturbation amplitude.
kx = 0.4 -- Perturbed wave number (x-direction).

-- Derived physical quantities (using normalized code units).
gamma_elc1 = 1.0 / math.sqrt(1.0 - (ux_elc1 * ux_elc1) - (uy_elc1 * uy_elc1) - (uz_elc1 * uz_elc1)) -- First electron gamma factor.
gamma_elc2 = 1.0 / math.sqrt(1.0 - (ux_elc2 * ux_elc2) - (uy_elc2 * uy_elc2) - (uz_elc2 * uz_elc2)) -- Second electron gamma factor.

ux_elc1_sr = gamma_elc1 * ux_elc1 -- First electron relativistic velocity (x-direction).
ux_elc2_sr = gamma_elc2 * ux_elc2 -- Second electron relativistic velocity (x-direction).
uy_elc1_sr = gamma_elc1 * uy_elc1 -- First electron relativistic velocity (y-direction).
uy_elc2_sr = gamma_elc2 * uy_elc2 -- Second electron relativistic velocity (y-direction).
uz_elc1_sr = gamma_elc1 * uz_elc1 -- First electron relativistic velocity (z-direction).
uz_elc2_sr = gamma_elc2 * uz_elc2 -- Second electron relativistic velocity (z-direction).

-- Simulation parameters.
Nx = 8 -- Cell count (configuration space: x-direction).
Nvx = 8 -- Cell count (velocity space: vx-direction).
Nvy = 8 -- Cell count (velocity space: vy-direction).
Nvz = 8 -- Cell count (velocity space: vz-direction).
Lx = 2.0 * pi / kx -- Domain size (configuration space: x-direction).
vx_max = 8.0 -- Domain boundary (velocity space: vx-direction).
vy_max = 8.0 -- Domain boundary (velocity space: vy-direction).
vz_max = 8.0 -- Domain boundary (velocity space: vz-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 10.0 -- Final simulation time.
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
          return n_elc1 -- Total left-going number density.
        end,
        temperatureInit = function (t, xn)
          return T_elc1 -- Total left-going temperature.
        end,
        driftVelocityInit = function (t, xn)
          return ux_elc1_sr, uy_elc1_sr, uz_elc1_sr -- Total left-going relativistic drift velocity.
        end,

        correctAllMoments = true,
        useLastConverged = true
      },
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return n_elc2 -- Total right-going number density.
        end,
        temperatureInit = function (t, xn)
          return T_elc2 -- Total right-going temperature.
        end,
        driftVelocityInit = function (t, xn)
          return ux_elc2_sr, uy_elc2_sr, uz_elc2_sr -- Total right-going relativistic drift velocity.
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

      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = alpha * math.sin(kx * x) -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,
    
    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()