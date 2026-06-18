local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

n_elc1 = 0.5 -- First electron number density.
n_elc2 = 0.5 -- Second electron number density.
ux_elc1 = 0.0 -- First electron velocity (x-direction).
ux_elc2 = 0.0 -- Second electron velocity (x-direction).
uy_elc1 = 0.3 -- First electron velocity (y-direction).
uy_elc2 = -0.3 -- Second electron velocity (y-direction).
T_elc1 = 0.01 -- First electron temperature.
T_elc2 = 0.01 -- Second electron temperature.

alpha = 1.0e-3 -- Applied perturbation amplitude.
kx = 0.4 -- Perturbed wave number (x-direction).

-- Derived physical quantities (using normalized code units).
vt_elc1 = math.sqrt(T_elc1 / mass_elc) -- First electron thermal velocity.
vt_elc2 = math.sqrt(T_elc2 / mass_elc) -- Second electron thermal velocity.

-- Simulation parameters.
Nx = 24 -- Cell count (configuration space: x-direction).
Nvx = 12 -- Cell count (velocity space: vx-direction).
Nvy = 12 -- Cell count (velocity space: vy-direction).
Lx = 2.0 * pi / kx -- Domain size (configuration space: x-direction).
vx_max = 1.0 -- Domain boundary (velocity space: vx-direction).
vy_max = 1.0 -- Domain boundary (velocity space: vy-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

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
    lower = { -vx_max, -vy_max },
    upper = { vx_max, vy_max },
    cells = { Nvx, Nvy },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local vx, vy = xn[2], xn[3]

          local v_sq_elc1 = ((vx - ux_elc1) * (vx - ux_elc1)) + ((vy - uy_elc1) * (vy - uy_elc1))
          local v_sq_elc2 = ((vx - ux_elc2) * (vx - ux_elc2)) + ((vy - uy_elc2) * (vy - uy_elc2))
        
          local maxwellian1 = (n_elc1 / (2.0 * pi * vt_elc1 * vt_elc1)) * math.exp(-v_sq_elc1 / (2.0 * vt_elc1 * vt_elc1))
          local maxwellian2 = (n_elc2 / (2.0 * pi * vt_elc2 * vt_elc2)) * math.exp(-v_sq_elc2 / (2.0 * vt_elc2 * vt_elc2))
          local n = maxwellian1 + maxwellian2 -- Distribution function.

          return n
        end
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
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
      local Bz = alpha * math.sin(kx * x) -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()