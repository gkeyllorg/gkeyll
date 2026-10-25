-- Relaxation of a waterbag and of a bump-on-tail distribution under LBO self-collisions with the
-- PKPM system of equations (1x, normalized parallel velocity w = v_par/xi), the counterpart of
-- vlasov/creg/rt_vlasov_lbo_relax_1x1v.c.
-- Two homogeneous neutral species carry the parallel distributions of that test: a waterbag,
-- |v_par| < 1, and a Maxwellian of thermal speed 1/3 carrying a narrow bump at 3.5 thermal speeds.
-- Both are Maxwellian in the perpendicular directions with T_perp = T_par, so the collisions keep
-- the temperature isotropic and the parallel distribution relaxes as under the one-dimensional
-- operator. A static magnetic field along x sets the field direction. Collision frequency 1, run
-- for 5 collision times.
-- Density, momentum and energy are conserved to 1e-10. The LBO relaxes the Hermite modes of the
-- parallel distribution at k nu: fitted while above the discretization floor, the heat flux decays
-- at 3.000 nu and the kurtosis at 4.000 nu (Vlasov 1x1v: 3.00 and 4.00). The discontinuous waterbag
-- projects with p_par 1.6% below p_perp; the collisions remove that anisotropy.

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass = 1.0 -- Species mass.
charge = 0.0 -- Species charge.

n0 = 1.0 -- Reference number density.
vt = 1.0 / 3.0 -- Thermal velocity of the Maxwellian core of the bump distribution.
nu = 1.0 -- Collision frequency.
ab = math.sqrt(0.1) -- Bump amplitude.
sb = 0.12 -- Bump width (Lorentzian softening).
ub = 4.0 * math.sqrt(0.25 / 3.0) -- Bump location (v_par).
vtb = 1.0 -- Bump Maxwellian thermal velocity.
B0 = 1.0 -- Magnetic field along x (sets the field direction b).

n_sq = n0 -- Waterbag number density.
T_sq = 1.0 / 3.0 -- Waterbag temperature (T/m, isotropic).
n_b = 1.951494084012749 * n0 -- Bump distribution number density (quadrature).
u_b = 0.5629997754412951 -- Bump distribution mean velocity (quadrature).
T_b = 0.4342855799839833 -- Bump distribution temperature T/m (quadrature, isotropic).

-- Simulation parameters.
Nx = 2 -- Cell count (configuration space: x-direction).
Nvx = 72 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 5.0 -- Domain size (configuration space: x-direction).
vx_max = 5.0 * math.sqrt(3.0) -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 5.0 -- Final simulation time.
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

  -- Waterbag species.
  square = PKPM.Species.new {
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2]

      -- Waterbag of half the reference density on |v_par| < 1, in w = v_par/xi with xi = sqrt(T_sq):
      -- F = xi f/n on |w| < 1/xi, and G = (T_perp/T) F = F.
      local xi = math.sqrt(T_sq) -- Normalizing thermal velocity.
      local F = 0.0
      if math.abs(w) < 1.0 / xi then
        F = 0.5 * n0 * xi / n_sq
      end

      return F, F
    end,

    initFluid = function (t, xn)
      local rho = mass * n_sq -- Mass density.
      local E = 1.5 * rho * T_sq -- Total energy density (at rest).

      return rho, 0.0, 0.0, 0.0, E
    end,

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu
      end
    },

    evolve = true -- Evolve species?
  },

  -- Bump-on-tail species.
  bump = PKPM.Species.new {
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2]

      -- Maxwellian core plus a bump with a Lorentzian envelope, at v_par = u_b + xi w, xi = sqrt(T_b):
      -- F = xi f/n and G = (T_perp/T) F = F.
      local xi = math.sqrt(T_b) -- Normalizing thermal velocity.
      local v = u_b + xi * w -- Parallel velocity.
      local vb_sq = (v - ub) * (v - ub) -- Squared distance to the bump.
      local f = (n0 / math.sqrt(2.0 * pi * vt * vt)) * math.exp(-v * v / (2.0 * vt * vt)) + (n0 / math.sqrt(2.0 * pi * vtb * vtb)) * math.exp(-vb_sq / (2.0 * vtb * vtb)) * (ab * ab) / (vb_sq + (sb * sb)) -- Parallel distribution.
      local F = xi * f / n_b -- Normalized parallel distribution.

      return F, F
    end,

    initFluid = function (t, xn)
      local rho = mass * n_b -- Mass density.
      local mom_x = rho * u_b -- Momentum density (x-direction).
      local E = 0.5 * rho * u_b * u_b + 1.5 * rho * T_b -- Total energy density.

      return rho, mom_x, 0.0, 0.0, E
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
      return 0.0, 0.0, 0.0, B0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
