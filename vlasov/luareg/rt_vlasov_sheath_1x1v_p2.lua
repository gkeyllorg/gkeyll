local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using non-normalized physical units).
epsilon0 = 8.854e-12 -- Permittivity of free space.
mu0 = 1.257e-6 -- Permeability of free space.
mass_elc = 9.109e-31 -- Electron mass.
charge_elc = -1.602e-19 -- Electron charge.
mass_ion = 1836.153 * 9.109e-31 -- Ion mass.
charge_ion = 1.602e-19 -- Ion charge.

n0 = 1.0e17 -- Reference number density.

-- Derived physical quantities (using non-normalized physical units).
Te = 10.0 * charge_ion -- Electron temperature.
Ti = 10.0 * charge_ion -- Ion temperature.

vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.

lambda_D = math.sqrt(epsilon0 * Te / (n0 * charge_ion * charge_ion)) -- Electron Debye length.
omega_pe = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)) -- Electron plasma frequency.

-- Simulation parameters.
Nx = 128 -- Cell count (configuration space: x-direction).
Nvx = 32 -- Cell count (velocity space: vx-direction).
Lx = 128.0 * lambda_D -- Domain size (configuration space: x-direction).
Ls = 100.0 * lambda_D -- Domain size (source).
vx_max = 4.0 * vte -- Domain boundary (velocity space: vx-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 10.0 / omega_pe -- Final simulation time.
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
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local vx = xn[2]

          local n = n0 / math.sqrt(2.0 * pi * (vte * vte)) * (math.exp(-(vx * vx) / (2.0 * (vte * vte)))) -- Electron distribution function.

          return n
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
          projectionID = G0.Projection.Func,

          init = function (t, xn)
            local x, vx = xn[1], xn[2]

            local n = 0.0

            if math.abs(x) < Ls then
              n = 2.0 * (Ls - math.abs(x)) / Ls * (1.0 / math.sqrt(2.0 * pi * (vte * vte)) * (math.exp(-(vx * vx) / (2.0 * (vte * vte))))) -- Electron source distribution function (left).
            else
              n = 0.0 -- Electron source distribution function (right).
            end

            return n
          end
        }
      }
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcWall
      },
      upper = {
        type = G0.SpeciesBc.bcAbsorb
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Ions.
  ion = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,
    
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
          local vx = xn[2]

          local n = n0 / math.sqrt(2.0 * pi * (vti * vti)) * (math.exp(-(vx * vx) / (2.0 * (vti * vti)))) -- Ion distribution function.

          return n
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
          projectionID = G0.Projection.Func,

          init = function (t, xn)
            local x, vx = xn[1], xn[2]

            local n = 0.0

            if math.abs(x) < Ls then
              n = 2.0 * (Ls - math.abs(x)) / Ls * (1.0 / math.sqrt(2.0 * pi * (vti * vti)) * (math.exp(-(vx * vx) / (2.0 * (vti * vti))))) -- Ion source distribution function (left).
            else
              n = 0.0 -- Ion source distribution function (right).
            end

            return n
          end
        }
      }
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcWall
      },
      upper = {
        type = G0.SpeciesBc.bcAbsorb
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
      local x = xn[1]

      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    bcx = {
      lower = {
        type = G0.FieldBc.bcSymWall
      },
      upper = {
        type = G0.FieldBc.bcWall
      }
    },

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()