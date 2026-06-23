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

theta = (45.0 / 180.0) * pi -- Perturbation angle.
R_elc = 0.333333333333333 -- Electron radius.

k0 = 1.0 -- Reference perturbed wave number.
alpha = 1.18281106421231 -- Applied perturbation amplitude.
perturb_n = 1.0e-8 -- Perturbation density.

-- Derived physical quantities (using normalized code units).
T_elc1 = mass_elc * ((R_elc * uy_elc1) * (R_elc * uy_elc1)) -- First electron temperature.
T_elc2 = mass_elc * ((R_elc * uy_elc1) * (R_elc * uy_elc1)) -- Second electron temperature.
vt_elc1 = math.sqrt(T_elc1 / mass_elc) -- First electron thermal velocity.
vt_elc2 = math.sqrt(T_elc2 / mass_elc) -- Second electron thermal velocity.

kx = k0 * math.cos(theta) -- Perturbed wave number (x-direction).
ky = k0 * math.sin(theta) -- Perturbed wave number (y-direction).

-- Simulation parameters.
Nx = 8 -- Cell count (configuration space: x-direction).
Ny = 8 -- Cell count (configuration space: y-direction).
Nvx = 16 -- Cell count (velocity space: vx-direction).
Nvy = 16 -- Cell count (velocity space: vy-direction).
Lx = 2.0 * pi / kx -- Domain size (configuration space: x-direction).
Ly = 2.0 * pi / ky -- Domain size (configuration space: y-direction).
vx_max = 0.9 -- Domain boundary (velocity space: vx-direction).
vy_max = 0.9 -- Domain boundary (velocity space: vy-direction).
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
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1, 2 }, -- Periodic directions (x- and y-directions only).

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
          local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]

          local v_sq_elc1 = ((vx - ux_elc1) * (vx - ux_elc1)) + ((vy - uy_elc1) * (vy - uy_elc1))
          local v_sq_elc2 = ((vx - ux_elc2) * (vx - ux_elc2)) + ((vy - uy_elc2) * (vy - uy_elc2))
        
          local maxwellian1 = (n_elc1 / (2.0 * pi * vt_elc1 * vt_elc1)) * math.exp(-v_sq_elc1 / (2.0 * vt_elc1 * vt_elc1))
          local maxwellian2 = (n_elc2 / (2.0 * pi * vt_elc2 * vt_elc2)) * math.exp(-v_sq_elc2 / (2.0 * vt_elc2 * vt_elc2))
          local n = (1.0 + (perturb_n * math.cos((kx * x) + (ky * y)))) * (maxwellian1 + maxwellian2) -- Distribution function.

          return n
        end
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1 }
  },

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local Ex = -perturb_n * math.sin((kx * x) + (ky * y)) / (kx + (ky * alpha)) -- Total electric field (x-direction).
      local Ey = alpha * Ex -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = (kx * Ey) - (ky * Ex) -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()