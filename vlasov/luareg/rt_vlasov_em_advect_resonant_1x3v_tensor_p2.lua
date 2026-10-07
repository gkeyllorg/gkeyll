-- Advection in specified electromagnetic fields (Vlasov, 1x3v). The initial conditions follow entry
-- JE32 of Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je32/je32-vlasov-test-
-- ptcl.html), rotated so that the oscillating electric field E_z = cos(omega t) is along z and the
-- background magnetic field B_0 = 1 along x; omega = 1.0 Omega_c with Omega_c = q B_0 / m the
-- cyclotron frequency, the resonant case. The fields are external and the species does not feed
-- back on them; the 3V setup mirrors the PKPM test of the same name. The resonant drive grows the
-- bulk velocity secularly, u ~ (q / m) E_0 t / 2: at t = 10, u_y = -2.72010 and u_z = 4.46736,
-- exact to 1e-5, and the temperature stays 1 to 1e-4.
local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

n0 = 1.0 -- Reference density. 
vt = 1.0 -- Thermal velocity.

-- External EM field parameters.
omega = 1.0 -- Oscillating electric field frequency normalized to cyclotron frequency.
B0 = 1.0 -- Reference magnetic field strength.

-- Simulation parameters.
Nx = 2 -- Cell count (configuration space: x-direction).
Nvx = 8 -- Cell count (velocity space: vx-direction).
Nvy = 8 -- Cell count (velocity space: vy-direction).
Nvz = 16 -- Cell count (velocity space: vz-direction).
Lx = 4.0 * pi -- Domain size (configuration space: x-direction).
vx_max = 8.0 -- Domain boundary (velocity space: vx-direction).
vy_max = 8.0 -- Domain boundary (velocity space: vy-direction).
vz_max = 8.0 -- Domain boundary (velocity space: vz-direction).
poly_order = 2 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 10.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { -0.5 * Lx },
  upper = { 0.5 * Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Electrons.
  elc = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,
    
    -- Velocity space grid.
    lower = { -vx_max, -vy_max, -vz_max },
    upper = { vx_max, vy_max, vz_max },
    cells = { Nvx, Nvy, Nvz },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]

          local n = n0

          return n
        end,
        temperatureInit = function (t, xn)
          local x = xn[1]

          local T = vt*vt*mass_elc

          return T
        end,
        driftVelocityInit = function (t, xn)
          return 0.0, 0.0, 0.0 -- Total drift velocity.
        end,

        correctAllMoments = true
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments }
  },

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    externalFieldInit = function (t, xn)
      local Ex = 0.0 -- External electric field (x-direction).
      local Ey = 0.0 -- External electric field (y-direction).
      local Ez = math.cos(omega * t) -- External electric field (z-direction).

      local Bx = B0 -- External magnetic field (x-direction).
      local By = 0.0 -- External magnetic field (y-direction).
      local Bz = 0.0 -- External magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz
    end,

    evolve = false, -- Evolve field?
    evolveExternalField = true, 
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()