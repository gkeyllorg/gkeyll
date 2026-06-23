local Vlasov = G0.Vlasov

-- Physical constants (using non-normalized physical units).
epsilon0 = 8.854187817620389850536563031710750260608e-12 -- Permittivity of free space.
mass_elc = 9.10938215e-31 -- Electron mass.
charge_elc = -1.602176487e-19 -- Electron charge.
mass_ion = 1.672621637e-27 -- Ion mass.
charge_ion = 1.602176487e-19 -- Ion charge.

n0 = 1.0e18 -- Reference number density.
Vx_drift_elc = 0.0 -- Electron drift velocity (x-direction).
Vx_drift_ion = 0.0 -- Ion drift velocity (x-direction).
Vy_drift_elc = 0.0 -- Electron drift velocity (y-direction).
Vy_drift_ion = 0.0 -- Ion drift velocity (y-direction).

B0 = 0.001 -- Reference magnetic field strength.

-- Derived physical quantities (using non-normalized physical units).
Te = 10.0 * charge_ion -- Electron temperature.
Ti = 10.0 * charge_ion -- Ion temperature.

vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.

lambda_D = math.sqrt(epsilon0 * Te / (n0 * charge_ion * charge_ion)) -- Electron Debye length.
omega_pe = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)) -- Electron plasma frequency.

-- Simulation parameters.
Nx = 64 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: vx-direction).
Nvy = 16 -- Cell count (velocity space: vy-direction).
Lx = 128.0 * lambda_D -- Domain size (configuration space: x-direction).
vx_max_elc = 6.0 * vte -- Domain boundary (electron velocity space: vx-direction).
vx_max_ion = 6.0 * vti -- Domain boundary (ion velocity space: vx-direction).
vy_max_elc = 6.0 * vte -- Domain boundary (electron velocity space: vy-direction).
vy_max_ion = 6.0 * vti -- Domain boundary (ion velocity space: vy-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 100.0 / omega_pe -- Final simulation time.
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

  -- Electrons.
  elc = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,
    
    -- Velocity space grid.
    lower = { -vx_max_elc, -vy_max_elc },
    upper = { vx_max_elc, vy_max_elc },
    cells = { Nvx, Nvy },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return n0 -- Electron total number density.
        end,
        temperatureInit = function (t, xn)
          return Te -- Electron isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift_elc, Vy_drift_elc -- Electron drift velocity.
        end
      }
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcAbsorb
      },
      upper = {
        type = G0.SpeciesBc.bcReflect
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.LTEMoments }
  },

  -- Ions.
  ion = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,
    
    -- Velocity space grid.
    lower = { -vx_max_ion, -vy_max_ion },
    upper = { vx_max_ion, vy_max_ion },
    cells = { Nvx, Nvy },

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
          return Vx_drift_ion, Vy_drift_ion -- Ion drift velocity.
        end
      }
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcAbsorb
      },
      upper = {
        type = G0.SpeciesBc.bcReflect
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.LTEMoments }
  },

  isElectrostatic = true,

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0,

    poissonBcs = {
      lowerType = {
        G0.PoissonBc.bcDirichlet
      },
      upperType = {
        G0.PoissonBc.bcNeumann
      },
      lowerValue = {
        0.0
      },
      upperValue = {
        0.0
      }
    },

    externalPotentialInit = function (t, xn)
      local x = xn[1]

      local phi = 0.0 -- External electric scalar potential.
      local Ax = 0.0 -- External magnetic vector potential (x-direction).
      local Ay = B0 * x -- External magnetic vector potential (y-direction).
      local Az = 0.0 -- External magnetic vector potential (z-direction).

      return phi, Ax, Ay, Az
    end,
    evolveExternalPotential = false
  }
}

vlasovApp:run()