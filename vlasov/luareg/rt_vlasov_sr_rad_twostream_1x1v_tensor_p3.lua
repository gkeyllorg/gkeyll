local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 1.0 -- Positron mass.
charge_ion = 1.0 -- Positron charge.

n0 = 1.0 -- Reference number density.
T = 0.1 -- Temperature (units of mc^2).
Vx_drift = 0.5 -- Drift velocity (x-direction).
t_cool = 3000.0 -- Cooling time. 

alpha = 1.0e-3 -- Applied perturbation amplitude.
mode_init = -16 -- Initial mode to perturb.
mode_final = 16 -- Final mode to perturb.

-- Derived physical quantities (using normalized code units).
gamma = 1.0 / math.sqrt(1.0 - (Vx_drift * Vx_drift)) -- Gamma factor.
Vx_drift_SR = gamma * Vx_drift -- Relativistic drift velocity (x-direction).

-- Simulation parameters.
Nx = 32 -- Cell count (configuration space: x-direction).
Nvx = 64 -- Cell count (velocity space: vx-direction).
Lx = 4.0 * pi -- Domain size (configuration space: x-direction); two wavelengths of the fastest growing mode.
kx = 2.0*pi/Lx -- Smallest wavenumber (largest wavelength) in the domain. 
vx_max = 8.0 -- Domain boundary (velocity space: vx-direction).
poly_order = 3 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 85.0 -- Final simulation time.
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
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]
          math.randomseed(0)
          local n = n0
          for i = mode_init, mode_final do 
            if (i ~= 0) then 
              n = n + alpha*math.random()*math.cos(kx*i*x + 2.0 * pi * math.random())
            end
          end
          return n
        end,
        temperatureInit = function (t, xn)
          return T 
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift_SR 
        end,

        correctAllMoments = true,
        useLastConverged = true
      },
    },

    radiation = {
      radiationID = G0.VMRadiation.Compton, 
      coolingTime = t_cool, 
      p0 = 0.0, 
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1 }
  },

  -- Positrons.
  ion = Vlasov.Species.new {
    modelID = G0.Model.SR,
    charge = charge_ion, mass = mass_ion,
    
    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return n0
        end,
        temperatureInit = function (t, xn)
          return T 
        end,
        driftVelocityInit = function (t, xn)
          return -Vx_drift_SR 
        end,

        correctAllMoments = true,
        useLastConverged = true
      },
    },

    radiation = {
      radiationID = G0.VMRadiation.Compton, 
      coolingTime = t_cool, 
      p0 = 0.0, 
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1 }
  },

  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]
      math.randomseed(0)
      local Ex = 0.0
      -- Total electric field (x-direction). 
      for i = mode_init, mode_final do 
        if (i ~= 0) then 
          Ex = Ex - alpha*gamma*math.random()*math.sin(kx*i*x + 2.0 * pi * math.random())/(i*kx)
        end
      end
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    appliedCurrent = function (t, xn)
      local Jx = 2.0*n0*Vx_drift_SR
      return Jx, 0.0, 0.0
    end,
    
    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()