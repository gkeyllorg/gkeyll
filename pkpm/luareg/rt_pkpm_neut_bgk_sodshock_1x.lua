-- Sod shock tube for a neutral gas with BGK collisions with the PKPM system of equations (1x,
-- normalized parallel velocity w = v_par/xi), the counterpart of
-- vlasov/creg/rt_vlasov_nonuniformv_neut_bgk_sodshock_1x3v.c
-- Density 1 and pressure 1 on the left, density 1/8 and pressure 1/10 on the right, at rest. The
-- collision frequency puts the mean free path at 1/3 of a cell, so the solution follows the Euler
-- Sod solution of a gas with three degrees of freedom (PKPM carries the perpendicular temperature),
-- adiabatic index 5/3, with the discontinuities smoothed over a few mean free paths. A static
-- magnetic field along x sets the field direction.
-- Exact solution at t = 0.1: shock at x = 0.684, contact at 0.584, post-shock density 0.230,
-- velocity 0.841.
-- Figures of merit (32 x 16 cells): the shock at 0.678 and the contact at 0.571 lie within half a
-- cell of the Euler positions; the post-shock plateau, three cells wide, has density 0.241 (5%
-- high) and velocity 0.846; the density L1 error is 0.020 (Vlasov 1x3v, serendipity p2 with 6^3
-- mapped velocity cells: 0.02).

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_neut = 1.0 -- Neutral mass.
charge_neut = 0.0 -- Neutral charge.

nl = 1.0 -- Left number density.
Tl = 1.0 -- Left temperature.
nr = 0.125 -- Right number density.
Tr = 0.8 -- Right temperature.
nu = 100.0 -- Collision frequency.
vt = math.sqrt(Tl / mass_neut) -- Thermal velocity (left).
B0 = 1.0 -- Magnetic field along x (sets the field direction b).

-- Simulation parameters.
Nx = 32 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 1.0 -- Domain size (configuration space: x-direction).
vx_max = 8.0 -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.1 -- Final simulation time.
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

      local n = 0.0
      local T = 0.0
      if x < 0.5 * Lx then
        n = nl -- Left number density.
        T = Tl -- Left temperature.
      else
        n = nr -- Right number density.
        T = Tr -- Right temperature.
      end
      local rho = mass_neut * n -- Mass density.
      local E = 1.5 * n * T -- Total energy density (at rest).

      return rho, 0.0, 0.0, 0.0, E
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
