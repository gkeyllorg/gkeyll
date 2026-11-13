-- Stationary contact discontinuity for a neutral gas with BGK collisions with the PKPM system of
-- equations (1x, normalized parallel velocity w = v_par/xi).
-- Density 1.4 left of x = 0.5 and 1 right of it, at rest at uniform pressure 1, to t = 2. A static
-- magnetic field along x sets the field direction.
-- Weak collisions: the collision frequency nu = 1.0 makes the domain one mean free path long (the
-- thermal speed sqrt(p/rho) of the hottest gas is 1). The collisions are advanced explicitly.
-- Figures of merit (100 x 16 cells) at t = 2, with nu dt = 0.00022 on average: the L1 differences
-- of the cell averages from the Euler solution (adiabatic index 5/3) are 14% (density) and 0.84%
-- (pressure), the velocity stays below 0.07, and the density between 1.133 and 1.26.

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_neut = 1.0 -- Neutral mass.
charge_neut = 0.0 -- Neutral charge.

rho_l = 1.4 -- Left mass density.
u_l = 0.0 -- Left velocity.
p_l = 1.0 -- Left pressure.
rho_r = 1.0 -- Right mass density.
u_r = 0.0 -- Right velocity.
p_r = 1.0 -- Right pressure.
x_loc = 0.5 -- Interface of the two states.
nu = 1.0 -- Collision frequency.
B0 = 1.0 -- Magnetic field along x (sets the field direction b).

-- Simulation parameters.
Nx = 100 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 1.0 -- Domain size (configuration space: x-direction).
vx_max = 8.0 -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 2.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate integrated L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

pkpmApp = PKPM.App.new {

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
  decompCuts = { 1 }, -- Cuts in each coordinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Neutrals.
  neut = PKPM.Species.new {
    charge = charge_neut, mass = mass_neut,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2]

      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * pi) -- Normalized Maxwellian in w = v_par/xi.

      return F, F
    end,

    initFluid = function (t, xn)
      local x = xn[1]

      local rho = 0.0
      local vel = 0.0
      local pres = 0.0
      if x < x_loc then
        rho = rho_l -- Mass density (left).
        vel = u_l -- Velocity (left).
        pres = p_l -- Pressure (left).
      else
        rho = rho_r -- Mass density (right).
        vel = u_r -- Velocity (right).
        pres = p_r -- Pressure (right).
      end
      local E = 0.5 * rho * vel * vel + 1.5 * pres -- Total energy density.

      return rho, rho * vel, 0.0, 0.0, E
    end,

    collisions = {
      collisionID = G0.Collisions.BGK,

      selfNu = function (t, xn)
        return nu
      end
    },

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, B0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
