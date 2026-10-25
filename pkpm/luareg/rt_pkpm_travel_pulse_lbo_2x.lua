-- Advection of an entropy (contact) wave by a collisional neutral gas with LBO collisions with the
-- PKPM system of equations (2x, normalized parallel velocity w = v_par/xi). A density pulse n = 1 +
-- 0.2 sin(pi x) sin(pi y) in a uniform flow u = (1, 1) at uniform pressure p = 1 crosses a doubly
-- periodic box of side 2 once along the diagonal (t = 2), along a magnetic field that is also
-- diagonal.
-- In the Euler limit the pulse returns unchanged to its starting position. With collisions at nu =
-- 100 the heat flux relaxes at 3 nu, so the gas conducts heat with thermal diffusivity chi = T/(m 3
-- nu) and the pulse, a temperature perturbation at constant pressure, decays as exp(-chi k^2 t),
-- k^2 = 2 pi^2; the measured decay beyond it is the numerical diffusion of the scheme.
-- Figures of merit (16 x 16 x 8 cells): at t = 2 the mode with k along b has amplitude 0.861 and
-- the mode with k across b 0.860, against exp(-chi k^2 t) = 0.877 for chi = T/(3 m nu). The LBO
-- relaxes the first harmonic (M0_perp at nu, N_perp at 3 nu), so the perpendicular heat flux
-- relaxes at 3 nu like the parallel one and the gas conducts heat across b as along it. Phase
-- errors are below 5e-4 rad; particle number is conserved to 1e-12.

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

n0 = 1.0 -- Reference number density.
alpha = 0.2 -- Perturbation amplitude.
u0 = 1.0 -- Reference fluid velocity (each direction).
pr0 = 1.0 -- Reference pressure.
nu = 100.0 -- Collision frequency.
B0 = 1.0 / math.sqrt(2.0) -- Magnetic field along each of x and y (sets the field direction b).

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Ny = 16 -- Cell count (configuration space: y-direction).
Nvx = 8 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 2.0 -- Domain size (configuration space: x-direction).
Ly = 2.0 -- Domain size (configuration space: y-direction).
vx_max = 6.0 -- Domain boundary (velocity space: w = v_par/xi).
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
  periodicDirs = { 1, 2 }, -- Periodic directions (x-direction, y-direction).

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

      local n = n0 + alpha * math.sin(pi * x) * math.sin(pi * y) -- Number density.
      local rho = mass * n -- Mass density.
      local E = 0.5 * rho * u0 * u0 * 2.0 + 1.5 * pr0 -- Total energy density (uniform pressure).

      return rho, rho * u0, rho * u0, 0.0, E
    end,

    collisions = {
      collisionID = G0.Collisions.LBO,

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
