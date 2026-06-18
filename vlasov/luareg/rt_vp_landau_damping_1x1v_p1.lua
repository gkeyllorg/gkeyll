local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

n0 = 1.0 -- Reference number density.
Te = 1.0 -- Electron temperature.

alpha = 1.0e-4 -- Applied perturbation amplitude.

-- Derived physical quantities (using normalized code units).
vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
omega_pe = math.sqrt((charge_elc * charge_elc) * n0 / (epsilon0 * mass_elc)) -- Electron plasma frequency.
lambda_D = vte / omega_pe -- Electron Debye length.

k0 = 0.5 / lambda_D -- Perturbed wave number.

-- Simulation parameters.
Nx = 64 -- Cell count (configuration space: x-direction).
Nvx = 64 -- Cell count (velocity space: vx-direction).
Lx = 2.0 * pi / k0 -- Domain size (configuration space: x-direction).
vx_max = 6.0 * vte -- Domain boundary (velocity space: vx-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.6 -- CFL coefficient.

t_end = 15.0 / omega_pe -- Final simulation time.
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
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,
    
    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local x, vx = xn[1], xn[2]

          local n = (1.0 + alpha * math.cos(k0 * x)) *
            (1.0 / math.sqrt(2.0 * pi * vte * vte)) * (math.exp(-(vx * vx) / (2.0 * vte * vte))) -- Distribution function.

          return n
        end
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  isElectrostatic = true,

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0,

    poissonBcs = {
      lowerType = {
        G0.PoissonBc.bcPeriodic
      },
      upperType = {
        G0.PoissonBc.bcPeriodic
      }
    }
  }
}

vlasovApp:run()