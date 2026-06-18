local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 1836.153 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

Te_over_Ti = 4.0 -- Electron temperature / ion temperature.

-- Derived physical quantities (using normalized code units).
vte = 1.0 -- Electron thermal velocity.
vti = vte / math.sqrt(Te_over_Ti * mass_ion) -- Ion thermal velocity.
cs = vte / math.sqrt(mass_ion) -- Sound speed.

Vx_drift = 2.0 * cs -- Drift velocity (x-direction).

nu_elc = 1.0e-4 -- Electron collision frequency.
nu_ion = (nu_elc / math.sqrt(mass_ion)) * (Te_over_Ti * math.sqrt(Te_over_Ti)) -- Ion collision frequency.

-- Simulation parameters.
Nx = 256 -- Cell count (configuration space: x-direction).
Nvx = 64 -- Cell count (velocity space: vx-direction).
Lx = 256.0 -- Domain size (configuration space: x-direction).
vx_max_elc = 6.0 * vte -- Domain boundary (electron velocity space: vx-direction).
vx_max_ion = 16.0 * vti -- Domain boundary (ion velocity space: vx-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 20.0 -- Final simulation time.
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
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local x, vx = xn[1], xn[2]

          local v_sq_m = (vx - Vx_drift) * (vx - Vx_drift)
          local v_sq_p = (vx + Vx_drift) * (vx + Vx_drift)

          local n = 0.0
          if x < 0.0 then
            n = (1.0 / math.sqrt(2.0 * pi * vte * vte)) * (math.exp(-v_sq_m / (2.0 * vte * vte))) -- Distribution function (left).
          else
            n = (1.0 / math.sqrt(2.0 * pi * vte * vte)) * (math.exp(-v_sq_p / (2.0 * vte * vte))) -- Distribution function (right).
          end

          return n
        end
      },
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_elc -- Collision frequency.
      end,
          
      numCrossCollisions = 1,
      collideWith = { "ion" }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
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
        projectionID = G0.Projection.Func,
  
        init = function (t, xn)
          local x, vx = xn[1], xn[2]
  
          local v_sq_m = (vx - Vx_drift) * (vx - Vx_drift)
          local v_sq_p = (vx + Vx_drift) * (vx + Vx_drift)
  
          local n = 0.0
          if x < 0.0 then
            n = (1.0 / math.sqrt(2.0 * pi * vti * vti)) * (math.exp(-v_sq_m / (2.0 * vti * vti))) -- Distribution function (left).
          else
            n = (1.0 / math.sqrt(2.0 * pi * vti * vti)) * (math.exp(-v_sq_p / (2.0 * vti * vti))) -- Distribution function (right).
          end

          return n
        end
      },
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ion -- Collision frequency.
      end,
    
      numCrossCollisions = 1,
      collideWith = { "elc" }
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
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()