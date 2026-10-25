-- Electrons in an electrostatic potential well with the PKPM system of equations (1x, normalized
-- parallel velocity w = v_par/xi), the counterpart of vlasov/creg/rt_vp_es_pot_well_1x1v.c. A
-- uniform Maxwellian (n = 1, v_t = 1) starts at rest in the field E = -sin(x) of the potential well
-- phi = -cos(x) and evolves with Maxwell's equations (omega_p = 1, Debye length 1). In one
-- dimension Ampere's law with the initial E is equivalent to the Poisson solution for the external
-- well plus the self-consistent potential of the plasma response, which is how the Vlasov-Poisson
-- test is posed. Electrons bunch at the bottom of the well and partially screen it. A uniform
-- magnetic field along x gives PKPM its field direction and does not act on the parallel motion.
-- The fluid-field coupling is integrated implicitly, split from the explicit SSP-RK3 step.
-- Figures of merit (32 x 24 cells): particle number and the total energy (kinetic, thermal and
-- field) are conserved to 2e-13. By t = 3 the density ranges from 0.65 to 1.37 and the
-- self-consistent potential reaches 0.37 of the external one (Vlasov-Poisson, serendipity p2:
-- 0.65-1.37 and 0.37).

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

vt = 1.0 -- Thermal velocity.
n0 = 1.0 -- Reference number density.
B0 = 1.0 -- Magnetic field along x (sets the field direction b).

-- Simulation parameters.
Nx = 32 -- Cell count (configuration space: x-direction).
Nvx = 24 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 2.0 * pi -- Domain size (configuration space: x-direction).
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

  -- Electrons.
  elc = PKPM.Species.new {
    charge = charge_elc, mass = mass_elc,

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
      local rho = mass_elc * n0 -- Electron mass density.
      local E = 1.5 * n0 * mass_elc * vt * vt -- Electron total energy density (at rest).

      return rho, 0.0, 0.0, 0.0, E
    end,

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Ex = -math.sin(x) -- Electric field of the well phi = -cos(x) (x-direction).

      return Ex, 0.0, 0.0, B0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
