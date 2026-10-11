-- Relativistic two-stream instability with the Vlasov-Maxwell system (1x1v).
-- Two counter-streaming Maxwell-Juttner electron beams (drift +-0.9 c, T = 0.04 mc^2, n0/2 each) on
-- a neutralizing ion background are perturbed in density at k = 0.5 omega_pe/c with amplitude 1e-5;
-- the electric field follows from Gauss's law. The momentum grid is uniform near p = 0 (cells of
-- width 1/32 over the inner half of the cells) and continues quadratically to |p| = 128.
-- Kinetic linear theory for the 1D Maxwell-Juttner beams at k = 0.5: a purely growing mode with
-- gamma = 0.0416 omega_pe, omega_pe = sqrt(n0 e^2 / (epsilon0 m)) = 1 with n0 the rest-frame
-- density.
-- Figures of merit (serendipity p3, 64 x 64 cells): the field energy grows at gamma = 0.0422 over t
-- = 110-220 and saturates at t = 251 with energy 2.37 at t = 265; electron number conserved to
-- 7e-13.

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

n0 = 1.0 -- Reference number density.
T = 0.04 -- Temperature (units of mc^2).
Vx_drift = 0.9 -- Drift velocity (x-direction).

alpha = 1.0e-5 -- Applied perturbation amplitude.
kx = 0.5 -- Perturbed wave number (x-direction).

-- Derived physical quantities (using normalized code units).
gamma = 1.0 / math.sqrt(1.0 - (Vx_drift * Vx_drift)) -- Gamma factor.
Vx_drift_SR = gamma * Vx_drift -- Relativistic drift velocity (x-direction).

-- Simulation parameters.
Nx = 64 -- Cell count (configuration space: x-direction).
Nvx = 64 -- Cell count (velocity space: vx-direction).
Lx = 2.0 * pi / kx -- Domain size (configuration space: x-direction).
vx_max = 128.0 -- Domain boundary (velocity space: vx-direction).
nonuniform_v_pow = 2.0 -- Quadratic velocity map. 
vx_linear_res = 1.0/32.0 -- Transition from linear to quadratic velocity map. 
poly_order = 3 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 265.0 -- Final simulation time.
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
    modelID = G0.Model.SR,
    charge = charge_elc, mass = mass_elc,
    
    -- Velocity space grid.
    lower = { -1.0 },
    upper = { 1.0 },
    cells = { Nvx },

    mapc2pVel = { 
      -- vx mapping 
      { 
        vmap = function (t, xn)
          local vc = xn[1]
          local vp = 0.0
          local ncells_linear = Nvx/2

          if (vc < 0.0) then 
            vp = vx_linear_res*ncells_linear*vc - vx_max*vc^nonuniform_v_pow
          else
            vp = vx_linear_res*ncells_linear*vc + vx_max*vc^nonuniform_v_pow
          end
          return vp
        end
      },
    },

    -- Initial conditions.
    numInit = 2,
    projections = {
      -- Two counter-streaming Maxwellians.
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]

          local n = 0.5 * (1.0 + alpha * math.cos(kx * x)) * n0 -- Total number density.
          return n
        end,
        temperatureInit = function (t, xn)
          return T -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift_SR -- Total left-going relativistic drift velocity.
        end,

        correctAllMoments = true,
        useLastConverged = true
      },
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]

          local n = 0.5 * (1.0 + alpha * math.cos(kx * x)) * n0 -- Total number density.
          return n
        end,
        temperatureInit = function (t, xn)
          return T -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return -Vx_drift_SR -- Total right-going relativistic drift velocity.
        end,

        correctAllMoments = true,
        useLastConverged = true
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1 }
  },

  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Ex = -alpha * gamma * math.sin(kx * x) / kx -- Total electric field (x-direction).
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