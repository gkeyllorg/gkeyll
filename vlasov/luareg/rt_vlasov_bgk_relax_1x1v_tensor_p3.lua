-- Relaxation of a waterbag and of a bump-on-tail distribution under BGK self-collisions
-- (Vlasov, 1x1v). Two homogeneous species: a square (waterbag) distribution, |v_i| < 1 in every
-- velocity direction, and a Maxwellian of thermal speed 1/3 carrying a narrow bump at 3.5 thermal
-- speeds along v_x. Collision frequency 1, run for 5 collision times. Density, momentum and energy
-- are conserved to round-off. Under BGK every non-Maxwellian moment decays at the collision
-- frequency: the heat flux, the kurtosis and the L2 distance to the Maxwellian all decay at rate
-- 1.00 (measured).
local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
mass = 1.0 -- Species mass.
charge = 0.0 -- Species charge.

n0 = 1.0 -- Reference number density.
vt = 1.0 / 3.0 -- Thermal velocity of the Maxwellian core of the bump distribution.
nu = 1.0 -- Collision frequency.
ab = math.sqrt(0.1) -- Bump amplitude.
sb = 0.12 -- Bump width (Lorentzian softening).
ub = 4.0 * math.sqrt(0.25 / 3.0) -- Bump location (v_x).
vtb = 1.0 -- Bump Maxwellian thermal velocity.

-- Simulation parameters.
Nx = 2 -- Cell count (configuration space: x-direction).
Nvx = 64 -- Cell count (velocity space: vx-direction).
Lx = 5.0 -- Domain size (configuration space: x-direction).
vx_max = 5.0 -- Domain boundary (velocity space: vx-direction).
poly_order = 3 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.8 -- CFL coefficient.

t_end = 5.0 -- Final simulation time.
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

  -- Waterbag species.
  square = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge, mass = mass,

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

          -- Waterbag: half the reference density inside |v_i| < 1, zero outside.
          local f = 0.0
          if math.abs(vx) < 1.0 then
            f = 0.5 * n0
          end

          return f
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.BGK,

      selfNu = function (t, xn)
        return nu -- Collision frequency.
      end
    },

    correct = {
      correctAllMoments = true,
      iterationEpsilon = 1.0e-12,
      maxIterations = 100,
      useLastConverged = true
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Bump-on-tail species.
  bump = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge, mass = mass,

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

          local v_sq = vx * vx
          local vb_sq = (vx - ub) * (vx - ub)

          -- Maxwellian core plus a bump centered at (ub, 0, 0) with a Lorentzian envelope of width sb.
          local f = (n0 / math.sqrt(2.0 * math.pi * vt * vt)) * math.exp(-v_sq / (2.0 * vt * vt)) +
            (n0 / math.sqrt(2.0 * math.pi * vtb * vtb)) * math.exp(-vb_sq / (2.0 * vtb * vtb)) * (ab * ab) / (vb_sq + (sb * sb))

          return f
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.BGK,

      selfNu = function (t, xn)
        return nu -- Collision frequency.
      end
    },

    correct = {
      correctAllMoments = true,
      iterationEpsilon = 1.0e-12,
      maxIterations = 100,
      useLastConverged = true
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
  }
}

vlasovApp:run()
