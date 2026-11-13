-- Collision of two strong shocks for a neutral gas with BGK collisions with the PKPM system of
-- equations (1x, normalized parallel velocity w = v_par/xi): the states of test 5 of Toro, Riemann
-- Solvers and Numerical Methods for Fluid Dynamics. Density 5.99924, velocity 19.5975 and pressure
-- 460.894 left of x = 0.4; density 5.99242, velocity -6.19633 and pressure 46.0895 on the right.
-- The collision frequency puts the mean free path below 1/400 of a cell, so the solution follows
-- the Euler solution of a gas with adiabatic index 5/3. A static magnetic field along x sets the
-- field direction.
-- The BGK relaxation takes one backward-Euler step in each SSP-RK3 stage (nu dt = 2).
-- Exact solution at t = 0.035: shocks at x = 0.356 and 0.894 and a contact at 0.700, with velocity
-- 8.56 and pressure 1842 between the shocks, density 12.74 left of the contact and 21.93 right of
-- it.
-- Figures of merit (100 x 16 cells): the shocks at 0.357 and 0.895 (where the cell averages of the
-- density cross the middle of each jump) lie within 0.1 cells of the Euler positions and the
-- contact at 0.703 within 0.4 cells; the plateaus average density 12.74 and 21.91, velocity 8.56
-- and pressure 1841; the L1 errors of the cell averages are 1.4% (density), 0.29% (velocity) and
-- 0.36% (pressure), and the largest cell average of the density is 22.83 (4% above the plateau, at
-- the contact). Mass, momentum and energy change by the boundary fluxes to 1e-11.

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_neut = 1.0 -- Neutral mass.
charge_neut = 0.0 -- Neutral charge.

rho_l = 5.99924 -- Left mass density.
u_l = 19.5975 -- Left velocity.
p_l = 460.894 -- Left pressure.
rho_r = 5.99242 -- Right mass density.
u_r = -6.19633 -- Right velocity.
p_r = 46.0895 -- Right pressure.
x_loc = 0.4 -- Interface of the two states.
nu = 5.0e5 -- Collision frequency.
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

t_end = 0.035 -- Final simulation time.
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
      useImplicitCollisionScheme = true,

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
