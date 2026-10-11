-- Buneman instability with the Vlasov-Maxwell system (1x1v).
-- Parameters follow entry JE33 of Ammar's Simulation Journal
-- (https://ammar-hakim.org/sj/je/je33/je33-buneman.html): m_i/m_e = 25, electron drift 0.159 so
-- that k v_d = omega_pe for the box mode k = 2 pi, vte = 0.02, vti = 0.001, density perturbation
-- 1e-5 with the electric field from Gauss's law. The ghost current cancels the net current so that
-- Ampere's law reduces to the periodic Poisson problem (see JE33).
-- Kinetic linear theory at k = 2 pi: gamma = 0.1976 omega_pe, omega_r = 0.1425 omega_pe (cold-beam
-- maximum 0.235).
-- Figures of merit (serendipity p3, 16 x 64 cells): the field energy grows at gamma = 0.1950 over t
-- = 20-65 and saturates at t = 72.1 with energy 1.5e-03 at t = 80; electron number conserved to
-- 1e-12.

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
mass_ion = 25.0 -- Ion mass. 
charge_elc = -1.0 -- Electron charge.
charge_ion = 1.0 -- Ion charge. 

n0 = 1.0 -- Reference number density.
Vx_drift_elc = 0.159 -- Electron drift velocity.
Vx_drift_ion = 0.0 -- Ion drift velocity
vte = 0.02 -- Electron thermal velocity.
vti = 0.001 -- Ion thermal velocity. 

alpha = 1.0e-5 -- Applied perturbation amplitude.

-- Derived physical quantities (using normalized code units).
Te = (vte * vte) * mass_elc -- Electron temperature.
Ti = (vti * vti) * mass_ion -- Ion temperature.

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Nvx = 64 -- Cell count (velocity space: vx-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
kx = 2.0 * pi / Lx -- Perturbed wave number (x-direction); kx * Vx_drift_elc is close to the electron plasma frequency.
vx_max_elc = 6.0 * Vx_drift_elc -- Domain boundary (velocity space: vx-direction).
vx_max_ion = 128.0 * vti -- Domain boundary (velocity space: vx-direction).
poly_order = 3 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.6 -- CFL coefficient.

t_end = 80.0 -- Final simulation time.
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
  lower = { 0.0 },
  upper = { Lx },
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
    lower = { -vx_max_elc },
    upper = { vx_max_elc },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]
          return n0 * (1.0 + alpha * math.cos(kx * x)) -- Electron total number density.
        end,
        temperatureInit = function (t, xn)
          return Te -- Electron isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift_elc -- Electron drift velocity.
        end
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Ions.
  ion = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,
    
    -- Velocity space grid.
    lower = { -vx_max_ion },
    upper = { vx_max_ion },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return n0 -- Ion total number density.
        end,
        temperatureInit = function (t, xn)
          return Ti -- Ion isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift_ion -- Ion drift velocity.
        end
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },


  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,
    useGhostCurrent = true, 

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Ex = -alpha * math.sin(kx * x) / kx -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,
    
    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()