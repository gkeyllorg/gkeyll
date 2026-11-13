-- Rayleigh-Taylor instability of a neutral gas with BGK collisions with the PKPM system of
-- equations (2x, normalized parallel velocity w = v_par/xi): a gas of density 2 above one of
-- density 1 in a gravitational acceleration 0.1 along -y, in the box [0, 1/4] x [0, 1/2] with
-- reflecting walls at the bottom and the top and periodic in x (the states and the gravity of
-- moments/creg/rt_euler_rt.c). The interface is flat, on the cell faces at y = 1/4; the pressure
-- is hydrostatic with 0.01 at the top wall; the perturbation is the velocity of the incompressible
-- interface mode with one wavelength, u_y = 0.01 cos(k x) exp(-k |y - 1/4|), k = 8 pi. A static
-- uniform magnetic field along y, the direction of gravity, sets the field direction b; the
-- charge is zero, so the field only names the direction of the kinetic velocity. Gravity is an
-- applied acceleration in the SSP-RK3 stages (explicit sources).
-- The collision frequency puts the mean free path at 1/140 of a cell at the bottom wall. The
-- BGK relaxation takes one backward-Euler step in each SSP-RK3 stage (nu dt = 2 to 4).
-- Linear theory for a sharp interface in an incompressible fluid: growth rate sqrt(A g k) = 0.915
-- with the Atwood number A = 1/3; the Euler solution of the moments app (96 x 192 cells, adiabatic
-- index 5/3) grows at 0.91 between t = 1.5 and 2.3.
-- Figures of merit (12 x 24 x 8 cells): the integral of rho u_y^2 grows at the rate 0.79 between
-- t = 1.5 and 2.3 (half the logarithmic slope), 13% below the Euler solution: the heat conduction
-- along b across the interface and the viscosity at this collision frequency, and 12 cells per
-- wavelength, reduce it. At t = 3.5 the spike of heavy gas has fallen to y = 0.18 and the bubble
-- risen to 0.30; the density stays between 0.83 and 2.34, and the solution is symmetric about the
-- middle of the box in x to round-off. The energy changes by the work of gravity, -g int M_y dt,
-- to 2e-10 of the energy; mass is conserved to 2e-13. Without the perturbation the interface
-- stays at rest: the velocity stays below 1e-3 and decays, and the heat conduction spreads the
-- density jump over two cells (density undershoot 4% at t = 0.5 and 2% at t = 3.5).

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
nu = 2.0e3 -- Collision frequency.
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
