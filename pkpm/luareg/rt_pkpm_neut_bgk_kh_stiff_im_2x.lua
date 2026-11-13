-- Kelvin-Helmholtz instability for a neutral gas with BGK collisions with the PKPM system of
-- equations (2x, normalized parallel velocity w = v_par/xi).
-- In a doubly periodic unit square a slab of density 2 moving at -0.5 along x (1/4 < y < 3/4) lies
-- between streams of density 1 moving at +0.5, at uniform pressure 1, with the smooth shear layers
-- of McNally, Lyra & Passy (2012) of width 0.05 and a velocity 0.01 sin(2 pi x) along y seeding
-- one wavelength, to t = 3. A static uniform magnetic field along x, the direction of the shear
-- flow, sets the field direction b.
-- Stiff collisions, nu = 1.0e7 (nu dt of order 1e4): the BGK relaxation takes one backward-Euler
-- step in each SSP-RK3 stage, so the transport alone sets the time step.
-- Figures of merit (16 x 16 x 8 cells) at t = 3, with nu dt = 9900 on average: the integral of rho
-- u_y^2 is 0.0196 at the end of the run (7.5e-05 at t = 0); mass is conserved to 2e-13.

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

rho_in = 2.0 -- Inner mass density (1/4 < y < 3/4).
u_in = -0.5 -- Inner x-velocity.
rho_out = 1.0 -- Outer mass density.
u_out = 0.5 -- Outer x-velocity.
p0 = 1.0 -- Pressure.
width = 0.05 -- Width of the shear layers.
alpha = 0.01 -- Amplitude of the y-velocity perturbation.
nu = 1.0e7 -- Collision frequency.
B0 = 1.0 -- Magnetic field along x (sets the field direction b, along the shear flow).

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Ny = 16 -- Cell count (configuration space: y-direction).
Nvx = 8 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 1.0 -- Domain size (configuration space: x-direction).
Ly = 1.0 -- Domain size (configuration space: y-direction).
vx_max = 6.0 -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 3.0 -- Final simulation time.
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

      -- Smooth shear layers at y = 1/4 and 3/4 (McNally, Lyra & Passy 2012).
      local rho_m = 0.5 * (rho_out - rho_in)
      local u_m = 0.5 * (u_out - u_in)
      local rho = 0.0
      local vx = 0.0
      if y < 0.25 then
        rho = rho_out - rho_m * math.exp((y - 0.25) / width)
        vx = u_out - u_m * math.exp((y - 0.25) / width)
      elseif y < 0.5 then
        rho = rho_in + rho_m * math.exp((0.25 - y) / width)
        vx = u_in + u_m * math.exp((0.25 - y) / width)
      elseif y < 0.75 then
        rho = rho_in + rho_m * math.exp((y - 0.75) / width)
        vx = u_in + u_m * math.exp((y - 0.75) / width)
      else
        rho = rho_out - rho_m * math.exp((0.75 - y) / width)
        vx = u_out - u_m * math.exp((0.75 - y) / width)
      end
      local vy = alpha * math.sin(2.0 * pi * x) -- Perturbation of the y-velocity (one wavelength in x).
      local E = 0.5 * rho * (vx * vx + vy * vy) + 1.5 * p0 -- Total energy density.

      return rho, rho * vx, rho * vy, 0.0, E
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
