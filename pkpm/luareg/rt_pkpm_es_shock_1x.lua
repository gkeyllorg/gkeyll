-- Electrostatic shock with the PKPM system of equations (1x, normalized parallel velocity w =
-- v_par/xi), the counterpart of vlasov/creg/rt_vlasov_es_shock_1x1v.c.
-- Half domain: a plasma (m_i/m_e = 100, Te/Ti = 4) flowing at Mach 2 toward a reflecting wall at
-- x = 0 piles up against its own reflection and launches a shock back upstream. The inflow boundary
-- at x = Lx holds the upstream state fixed (reservoir). Ions slower than the shock potential are
-- reflected ahead of it. A uniform magnetic field along x sets the field direction, so without
-- collisions the perpendicular velocities are inert and the 1x1v Vlasov run is a one-to-one
-- reference.
-- The ion velocity grid reaches w = +-16. A reflected ion enters a cold upstream cell about 11 of
-- that cell's thermal speeds from the inflow, and a population beyond the grid of the cell receiving
-- it is carried by the fluid alone, which misplaces its energy between the parallel and
-- perpendicular temperatures: with +-12 the front is 14% fast and the ion T_perp rises from 0.25 to
-- 0.41 Te, and with +-6 the reflected ions pile up behind a second front.
-- Figures of merit at t_end (96 cells, 16 electron and 32 ion velocity cells; Vlasov 1x1v,
-- serendipity p2 on the same grid, in parentheses), the front being where the potential crosses half
-- its jump: front at x = 69.8 lambda_D, speed 1.16 c_s (66.7, 1.11); downstream density 2.19 n0
-- (2.20); potential jump across the front e*dphi/Te = 3.27 (3.37); reflected ions ahead of the front
-- 0.235 n0 at 3.52 c_s (0.242 at 3.48); downstream parallel temperatures 3.59 Te for the electrons
-- (3.65) and 1.04 Te for the ions (0.97). The perpendicular temperatures stay at their upstream
-- values, 1.00 Te and 0.26 Te (0.25).

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 / 25.0 -- Permeability of free space (speed of light = 5 vte).
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 100.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
Te_over_Ti = 4.0 -- Electron temperature / ion temperature.
mach_num = 2.0 -- Mach number of the inflow (relative to the ion sound speed).
B0 = 1.0 -- Magnetic field along x (sets the field direction b).

-- Derived physical quantities (using normalized code units).
Te = 1.0 * charge_ion -- Electron temperature.
Ti = Te / Te_over_Ti -- Ion temperature.

vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.
cs = math.sqrt(Te / mass_ion) -- Ion sound speed.
Vx_drift = -mach_num * cs -- Drift velocity of the inflow (x-direction).

lambda_D = math.sqrt(epsilon0 * Te / (n0 * charge_ion * charge_ion)) -- Electron Debye length.
omega_pe = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)) -- Electron plasma frequency.

-- Simulation parameters.
Nx = 96 -- Cell count (configuration space: x-direction).
Nvx_elc = 16 -- Cell count (electron velocity space: normalized parallel velocity w).
Nvx_ion = 32 -- Cell count (ion velocity space: normalized parallel velocity w).
Lx = 192.0 * lambda_D -- Domain size (configuration space: x-direction).
vx_max_elc = 6.0 -- Domain boundary (electron velocity space: w = v_par/xi).
vx_max_ion = 16.0 -- Domain boundary (ion velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 600.0 / omega_pe -- Final simulation time.
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

  -- Electrons.
  elc = PKPM.Species.new {
    charge = charge_elc, mass = mass_elc,

    -- Velocity space grid.
    lower = { -vx_max_elc },
    upper = { vx_max_elc },
    cells = { Nvx_elc },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2]

      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * pi) -- Normalized Maxwellian in w = v_par/xi.

      return F, F
    end,

    initFluid = function (t, xn)
      local rho = mass_elc * n0 -- Electron mass density.
      local ux = Vx_drift -- Electron drift.
      local E = 0.5 * rho * ux * ux + 1.5 * rho * vte * vte -- Electron total energy density.

      return rho, rho * ux, 0.0, 0.0, E
    end,

    -- Reflecting wall at x = 0, the upstream state held fixed at x = Lx.
    bcx = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcFixedFunc },

    evolve = true -- Evolve species?
  },

  -- Ions.
  ion = PKPM.Species.new {
    charge = charge_ion, mass = mass_ion,

    -- Velocity space grid.
    lower = { -vx_max_ion },
    upper = { vx_max_ion },
    cells = { Nvx_ion },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2]

      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * pi) -- Normalized Maxwellian in w = v_par/xi.

      return F, F
    end,

    initFluid = function (t, xn)
      local rho = mass_ion * n0 -- Ion mass density.
      local ux = Vx_drift -- Ion drift.
      local E = 0.5 * rho * ux * ux + 1.5 * rho * vti * vti -- Ion total energy density.

      return rho, rho * ux, 0.0, 0.0, E
    end,

    -- Reflecting wall at x = 0, the upstream state held fixed at x = Lx.
    bcx = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcFixedFunc },

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    externalFieldInit = function (t, xn)
      return 0.0, 0.0, 0.0, B0, 0.0, 0.0
    end,

    bcx = { G0.FieldBc.bcSymWall, G0.FieldBc.bcCopy }, -- Symmetry boundary and open boundary (x-direction).

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
