-- Rayleigh-Taylor instability for a neutral gas with BGK collisions with the PKPM system of
-- equations (2x, normalized parallel velocity w = v_par/xi).
-- A gas of density 2 above one of density 1 in a gravitational acceleration 0.1 along -y, in the
-- box [0, 1/4] x [0, 1/2] with reflecting walls at the bottom and the top and periodic in x, to t
-- = 3.5. The interface is flat at y = 1/4, the pressure is hydrostatic with 0.01 at the top wall,
-- and the perturbation is the velocity of the incompressible interface mode with one wavelength. A
-- static uniform magnetic field along y, the direction of gravity, sets the field direction b.
-- Gravity is an applied acceleration (explicit sources).
-- Stiff collisions, nu = 7.5e6 (nu dt of order 1e4): the BGK relaxation takes one backward-Euler
-- step in each SSP-RK3 stage, so the transport alone sets the time step.
-- Figures of merit (12 x 24 x 8 cells) at t = 3.5, with nu dt = 8400 on average: the integral of
-- rho u_y^2 is 7.83e-05 at the end of the run (7.46e-07 at t = 0); mass is conserved to 2e-13.

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

grav = 0.1 -- Gravitational acceleration (along -y).
rho_top = 2.0 -- Top mass density.
rho_bot = 1.0 -- Bottom mass density.
p_ref = 0.01 -- Pressure at the top wall.
pert_max = 0.01 -- Amplitude of the velocity perturbation.
nu = 7.5e6 -- Collision frequency.
B0 = 1.0 -- Magnetic field along y (sets the field direction b, along gravity).

-- Simulation parameters.
Nx = 12 -- Cell count (configuration space: x-direction).
Ny = 24 -- Cell count (configuration space: y-direction).
Nvx = 8 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 0.25 -- Domain size (configuration space: x-direction).
Ly = 0.5 -- Domain size (configuration space: y-direction).
vx_max = 6.0 -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 3.5 -- Final simulation time.
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

  useExplicitSource = true, -- Gravity in the SSP-RK3 stages.

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coordinate direction (x-direction, y-direction).

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
      local w = xn[3]

      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * pi) -- Normalized Maxwellian in w = v_par/xi.

      return F, F
    end,

    initFluid = function (t, xn)
      local x, y = xn[1], xn[2]

      local kx = 2.0 * pi / Lx -- Wavenumber of the perturbation.
      local yloc = 0.5 * Ly -- Interface.

      local rho = 0.0
      local pres = 0.0
      if y > yloc then
        rho = rho_top -- Mass density (top).
        pres = p_ref + rho_top * grav * (Ly - y) -- Hydrostatic pressure (top).
      else
        rho = rho_bot -- Mass density (bottom).
        pres = p_ref + rho_top * grav * (Ly - yloc) + rho_bot * grav * (yloc - y) -- Hydrostatic pressure (bottom).
      end
      -- Velocity of the incompressible interface mode.
      local decay = math.exp(-kx * math.abs(y - yloc))
      local sgn = -1.0
      if y > yloc then
        sgn = 1.0
      end
      local vx = sgn * pert_max * math.sin(kx * x) * decay
      local vy = pert_max * math.cos(kx * x) * decay
      local E = 0.5 * rho * (vx * vx + vy * vy) + 1.5 * pres -- Total energy density.

      return rho, rho * vx, rho * vy, 0.0, E
    end,

    -- Gravity.
    appliedAcceleration = function (t, xn)
      return 0.0, -grav, 0.0
    end,

    bcy = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcReflect }, -- Reflecting walls (y-direction).

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
      return 0.0, 0.0, 0.0, 0.0, B0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
