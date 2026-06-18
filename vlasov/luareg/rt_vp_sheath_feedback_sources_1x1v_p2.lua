local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 1836.153 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
Vx_drift_elc = 0.0 -- Electron drift velocity (x-direction).
Vx_drift_ion = 0.0 -- Ion drift velocity (x-direction).

-- Derived physical quantities (using normalized code units).
Te = 1.0 * charge_ion -- Electron temperature.
Ti = 1.0 * charge_ion -- Ion temperature.

vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.

lambda_D = math.sqrt(epsilon0 * Te / (n0 * charge_ion * charge_ion)) -- Electron Debye length.
omega_pe = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)) -- Electron plasma frequency.

-- Simulation parameters.
Nx = 256 -- Cell count (configuration space: x-direction).
Nvx = 64 -- Cell count (velocity space: vx-direction).
Lx = 256.0 * lambda_D -- Domain size (configuration space: x-direction).
Ls = 100.0 * lambda_D -- Domain size (source).
vx_max_elc = 6.0 * vte -- Domain boundary (electron velocity space: vx-direction).
vx_max_ion = 6.0 * vti -- Domain boundary (ion velocity space: vx-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 20.0 / omega_pe -- Final simulation time.
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
  periodicDirs = { }, -- Periodic directions (none).

  -- Electrons.
  elc = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,
    
    -- Velocity space grid.
    lower = { -vx_max_elc },
    upper = { vx_max_elc },
    cells = { Nvx },

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
          return Vx_drift_elc -- Electron drift velocity.
        end
      }
    },

    source = {
      sourceID = G0.Source.BoundaryFlux,
      sourceLength = Ls,
      sourceSpecies = "ion",

      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.LTE,

          densityInit = function (t, xn)
            local x = xn[1]

            local n = 0.0

            if math.abs(x) < Ls then
              n = (Ls - math.abs(x)) / Ls -- Electron source total number density (left).
            else
              n = 0.0 -- Electron source total number density (right).
            end

            return n
          end,
          temperatureInit = function (t, xn)
            return Te -- Electron source isotropic temperature.
          end,
          driftVelocityInit = function (t, xn)
            return Vx_drift_elc -- Electron source drift velocity.
          end
        }
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
    diagnostics = { G0.Moment.LTEMoments }
  },

  -- Ions.
  ion = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,
    
    -- Velocity space grid.
    lower = { -vx_max_ion },
    upper = { vx_max_ion },
    cells = { Nvx },

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
          return Vx_drift_ion -- Ion drift velocity.
        end
      }
    },

    source = {
      sourceID = G0.Source.BoundaryFlux,
      sourceLength = Ls,
      sourceSpecies = "ion",

      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.LTE,

          densityInit = function (t, xn)
            local x = xn[1]

            local n = 0.0

            if math.abs(x) < Ls then
              n = (Ls - math.abs(x)) / Ls -- Ion source total number density (left).
            else
              n = 0.0 -- Ion source total number density (right).
            end

            return n
          end,
          temperatureInit = function (t, xn)
            return Ti -- Ion source isotropic temperature.
          end,
          driftVelocityInit = function (t, xn)
            return Vx_drift_ion -- Ion source drift velocity.
          end
        }
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
        G0.PoissonBc.bcDirichlet
      },
      lowerValue = {
        0.0
      },
      upperValue = {
        0.0
      }
    }
  }
}

vlasovApp:run()