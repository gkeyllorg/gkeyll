local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
mass_neut1 = 1.0 -- First neutral mass.
mass_neut2 = 0.05 -- Second neutral mass.
charge_neut1 = 0.0 -- First neutral charge.
charge_neut2 = 0.0 -- Second neutral charge.

p_neut1 = 1.0 -- First neutral pressure.
p_neut2 = 0.5 -- Second neutral pressure.
n0_neut1 = 1.0 -- First neutral reference number density.
n0_neut2 = 1.0 -- Second netural reference number density.

ux0_neut1 = 0.1 -- First neutral reference velocity (x-direction).
ux0_neut2 = 2.5 -- Second neutral reference velocity (x-direction).
uy0_neut1 = 0.0 -- First neutral reference velocity (y-direction).
uy0_neut2 = 0.0 -- Second neutral reference velocity (y-direction).
nu_neut1 = 1.0 / 0.01 -- First neutral collision frequency.
nu_neut2 = math.sqrt(0.5 / 0.05) / 0.01 -- Second neutral collision frequency.

-- Derived physical quantities (using normalized code units).
vt_neut1 = math.sqrt(p_neut1 / (n0_neut1 * mass_neut1)) -- First neutral thermal velocity.
vt_neut2 = math.sqrt(p_neut2 / (n0_neut2 * mass_neut2)) -- Second neutral thermal velocity.

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Nvx = 32 -- Cell count (velocity space: vx-direction).
Nvy = 32 -- Cell count (velocity space: vy-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
vx_max_neut1 = 6.0 * vt_neut1 -- First neutral domain boundary (velocity space: vx-direction).
vx_max_neut2 = 6.0 * vt_neut2 -- Second neutral domain boundary (velocity space: vx-direction).
vy_max_neut1 = 6.0 * vt_neut1 -- First neutral domain boundary (velocity space: vy-direction).
vy_max_neut2 = 6.0 * vt_neut2 -- Second neutral domain boundary (velocity space: vy-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.0025 -- Final simulation time.
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
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- First neutral species.
  neut1 = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_neut1, mass = mass_neut1,
    
    -- Velocity space grid.
    lower = { -vx_max_neut1, -vy_max_neut1 },
    upper = { vx_max_neut1, vy_max_neut1 },
    cells = { Nvx, Nvy },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local vx, vy = xn[2], xn[3]

          local v_sq = ((vx - ux0_neut1) * (vx - ux0_neut1)) + ((vy - uy0_neut1) * (vy - uy0_neut1))
          local n = (n0_neut1 / math.sqrt(2.0 * pi * vt_neut1 * vt_neut1)) * math.exp(-v_sq / (2.0 * vt_neut1 * vt_neut1)) -- Distribution function.

          return n
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_neut1 -- Collision frequency.
      end,
      
      numCrossCollisions = 1,
      collideWith = { "neut2" }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Second neutral species.
  neut2 = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_neut2, mass = mass_neut2,
    
    -- Velocity space grid.
    lower = { -vx_max_neut2, -vy_max_neut2 },
    upper = { vx_max_neut2, vy_max_neut2 },
    cells = { Nvx, Nvy },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local vx, vy = xn[2], xn[3]

          local v_sq = ((vx - ux0_neut2) * (vx - ux0_neut2)) + ((vy - uy0_neut2) * (vy - uy0_neut2))
          local n = (n0_neut2 / math.sqrt(2.0 * pi * vt_neut2 * vt_neut2)) * math.exp(-v_sq / (2.0 * vt_neut2 * vt_neut2)) -- Distribution function.

          return n
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_neut2 -- Collision frequency.
      end,
      
      numCrossCollisions = 1,
      collideWith = { "neut1" }
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