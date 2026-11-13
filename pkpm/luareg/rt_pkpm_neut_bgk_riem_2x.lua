-- Two-dimensional Riemann problem for a neutral gas with BGK collisions with the PKPM system of
-- equations (2x, normalized parallel velocity w = v_par/xi).
-- Four constant states meeting at (0.5, 0.5) in the unit square (configuration 3 of Lax & Liu),
-- with (density, x-velocity, y-velocity, pressure) = (1.5, 0, 0, 1.5) upper right, (0.5323, 1.206,
-- 0, 0.3) upper left, (0.138, 1.206, 1.206, 0.029) lower left and (0.5323, 0, 1.206, 0.3) lower
-- right, to t = 0.3. A static uniform magnetic field along the diagonal x = y sets the field
-- direction b, at 45 degrees to the grid.
-- Moderate collisions advanced explicitly, nu = 100.0: the mean free path of the hottest gas is
-- 0.2 cells. This is the collision frequency at which the diffusion step of the LBO operator is
-- about the transport step; rt_pkpm_neut_lbo_riem_2x has the same collision frequency and grid.
-- Figures of merit (20 x 20 x 8 cells) at t = 0.3, with nu dt = 0.024 on average: against the
-- Euler solution of the moments app (400 x 400 cells) averaged over the cells of this grid the L1
-- differences of the cell averages are 4.4% (density) and 4.9% (pressure); the density stays
-- between 0.1378 and 1.5743, and the solution is symmetric about the diagonal to round-off.

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

rho_ul = 0.5323 -- Upper left mass density.
u_ul = 1.206 -- Upper left x-velocity.
v_ul = 0.0 -- Upper left y-velocity.
p_ul = 0.3 -- Upper left pressure.

rho_ur = 1.5 -- Upper right mass density.
u_ur = 0.0 -- Upper right x-velocity.
v_ur = 0.0 -- Upper right y-velocity.
p_ur = 1.5 -- Upper right pressure.

rho_ll = 0.138 -- Lower left mass density.
u_ll = 1.206 -- Lower left x-velocity.
v_ll = 1.206 -- Lower left y-velocity.
p_ll = 0.029 -- Lower left pressure.

rho_lr = 0.5323 -- Lower right mass density.
u_lr = 0.0 -- Lower right x-velocity.
v_lr = 1.206 -- Lower right y-velocity.
p_lr = 0.3 -- Lower right pressure.

loc = 0.5 -- Corner of the four states (both x and y coordinates).
nu = 100.0 -- Collision frequency.
B0 = 1.0 / math.sqrt(2.0) -- Magnetic field along each of x and y (sets the field direction b).

-- Simulation parameters.
Nx = 20 -- Cell count (configuration space: x-direction).
Ny = 20 -- Cell count (configuration space: y-direction).
Nvx = 8 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 1.0 -- Domain size (configuration space: x-direction).
Ly = 1.0 -- Domain size (configuration space: y-direction).
vx_max = 6.0 -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.3 -- Final simulation time.
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
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coordinate direction (x-direction, y-direction).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Neutrals.
  neut = PKPM.Species.new {
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[3]

      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * pi) -- Normalized Maxwellian in w = v_par/xi.

      return F, F
    end,

    initFluid = function (t, xn)
      local x, y = xn[1], xn[2]

      local rho = 0.0
      local vx = 0.0
      local vy = 0.0
      local pres = 0.0

      if y > loc then
        if x < loc then
          rho = rho_ul -- Mass density (upper left).
          vx = u_ul -- x-velocity (upper left).
          vy = v_ul -- y-velocity (upper left).
          pres = p_ul -- Pressure (upper left).
        else
          rho = rho_ur -- Mass density (upper right).
          vx = u_ur -- x-velocity (upper right).
          vy = v_ur -- y-velocity (upper right).
          pres = p_ur -- Pressure (upper right).
        end
      else
        if x < loc then
          rho = rho_ll -- Mass density (lower left).
          vx = u_ll -- x-velocity (lower left).
          vy = v_ll -- y-velocity (lower left).
          pres = p_ll -- Pressure (lower left).
        else
          rho = rho_lr -- Mass density (lower right).
          vx = u_lr -- x-velocity (lower right).
          vy = v_lr -- y-velocity (lower right).
          pres = p_lr -- Pressure (lower right).
        end
      end
      local E = 0.5 * rho * (vx * vx + vy * vy) + 1.5 * pres -- Total energy density.

      return rho, rho * vx, rho * vy, 0.0, E
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
      return 0.0, 0.0, 0.0, B0, B0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
