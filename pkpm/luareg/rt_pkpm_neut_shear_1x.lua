-- Sheared flow of a collisionless neutral gas across an oblique magnetic field with the PKPM system
-- of equations (1x, normalized parallel velocity w = v_par/xi). A uniform gas (n = 1, p = 1) with
-- the flow u_y = 0.5 cos(pi x) sits in the uniform static field B = (sin 0.3, -cos 0.3, 0) in a
-- periodic box of length 2. The gas has no charge, so the field only sets the direction b. The rate
-- of strain along the field, bb:grad(u) = b_x b_y du_y/dx, is not zero, and the flow exchanges
-- energy with the parallel and perpendicular pressures through the stress P_xy = (p_par - p_perp)
-- b_x b_y + s_x b_y + b_x s_y, s the first-harmonic (parallel-perpendicular) stress.
-- Figures of merit (16 x 16 cells) to t = 2: particle number, the three components of the momentum
-- and the total energy (kinetic + p_par/2 + p_perp) are conserved to 2e-13. The kinetic energy
-- falls from 0.125 to 0.0028 (t = 0.74) and is 0.0100 at t = 2, where the integrated parallel and
-- perpendicular pressures are 2.2275 and 2.0013 (both 2 initially) and the largest cell averages of
-- the two parts of P_xy are 0.051 (pressure anisotropy) and 0.155 (first harmonic). The positivity
-- limiter scales no cell.

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

n0 = 1.0 -- Reference number density.
uy0 = 0.5 -- Amplitude of the sheared flow (y-direction).
pr0 = 1.0 -- Reference pressure.
B0 = 1.0 -- Magnetic field strength (in the x-y plane).
theta = 0.3 -- Angle of the magnetic field from the y-direction (b_x = sin(theta)).

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 2.0 -- Domain size (configuration space: x-direction).
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
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Neutrals.
  neut = PKPM.Species.new {
    charge = charge, mass = mass,

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

      local rho = mass * n0 -- Mass density.
      local uy = uy0 * math.cos(pi * x) -- Sheared flow.
      local E = 0.5 * rho * uy * uy + 1.5 * pr0 -- Total energy density (uniform pressure).

      return rho, 0.0, rho * uy, 0.0, E
    end,

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, B0 * math.sin(theta), -B0 * math.cos(theta), 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
