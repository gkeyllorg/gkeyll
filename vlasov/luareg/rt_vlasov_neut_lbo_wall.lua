local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
mass_neut = 1.0 -- Neutral mass.
charge_neut = 0.0 -- Neutral charge.

n0 = 1.0 -- Reference number density.
u0 = 1.0 -- Reference velocity.
vt = 0.5 -- Neutral thermal velocity.
nu = 10.0 -- Collision frequency.

-- Simulation parameters.
Nx = 128 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: vx-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
vx_max = 12.0 * vt -- Domain boundary (velocity space: vx-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.3 -- Final simulation time.
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

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_neut, mass = mass_neut,
    
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

          local v_sq = (vx - u0) * (vx - u0)
          local n = (n0 / math.sqrt(2.0 * pi * vt * vt)) * math.exp(-v_sq / (2.0 * vt *  vt)) -- Distribution function.

          return n
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu -- Collision frequency.
      end,
      
      correctAllMoments = true
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcWall
      },
      upper = {
        type = G0.SpeciesBc.bcWall
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  skipField = true,

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = 1.0, mu0 = 1.0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0,

    isStatic = true
  }
}

vlasovApp:run()