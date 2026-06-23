-- 2D Wald magnetosphere problem for a non-static (Kerr) black hole, for the general relativistic Maxwell equations.
-- Input parameters describe a uniform magnetic field surrounding a rotating black hole.
-- Based on the analytical solution for force-free electrodynamics presented in the article:
-- R. M. Wald (1974), "Black hole in a uniform magnetic field",
-- Physical Review D, Volume 10 (6): 1680.
-- https://journals.aps.org/prd/abstract/10.1103/PhysRevD.10.1680

local Moments = G0.Moments
local GRMaxwell = G0.Moments.Eq.GRMaxwell
local BlackHole = G0.Moments.Spacetime.BlackHole

-- Physical constants (using normalized code units).
light_speed = 1.0 -- Speed of light.
e_fact = 0.0 -- Factor of speed of light for electric field correction.
b_fact = 0.0 -- Factor of speed of light for magnetic field correction.

B0 = 1.0 -- Reference magnetic field strength.

-- Spacetime parameters (using geometric units).
mass = 1.0 -- Mass of the black hole.
spin = -0.9 -- Spin of the black hole.

pos_x = 0.0 -- Position of the black hole (x-direction).
pos_y = 0.0 -- Position of the black hole (y-direction).
pos_z = 0.0 -- Position of the black hole (z-direction).

-- Simulation parameters.
Nx = 256 -- Cell count (x-direction).
Ny = 256 -- Cell count (y-direction).
Lx = 10.0 -- Domain size (x-direction).
Ly = 10.0 -- Domain size (y-direction).
cfl_frac = 0.95 -- CFL coefficient.

reinit_freq = 10 -- Spacetime reinitialization frequency.

t_end = 50.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

momentApp = Moments.App.new {
  
  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { -0.5 * Lx, -0.5 * Ly },
  upper = { 0.5 * Lx, 0.5 * Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x- and y-directions).
  
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Field.
  field = Moments.Species.new {
    equation = GRMaxwell.new {
      lightSpeed = light_speed,
      elcErrorSpeedFactor = e_fact,
      mgnErrorSpeedFactor = b_fact,
      blackHoleParameters = {
        mass = mass,
        spin = spin,
        posX = pos_x,
        posY = pos_y,
        posZ = pos_z
      },
      reinitFreq = reinit_freq
    },
  
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local Dx = 0.0 -- Total electric field (x-direction).
      local Dy = 0.0 -- Total electric field (y-direction).
      local Dz = 0.0 -- Total electric field (z-direction).

      local Bx = 0.125 * x * B0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      local lapse = BlackHole.lapseFunction(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)
      local shift = BlackHole.shiftVector(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)
      local spatial_metric = BlackHole.spatialMetricTensor(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)
      local in_excision_region = BlackHole.excisionRegion(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)

      local excision = 0.0
      if in_excision_region then
        Dx, Dy, Dz, Bx, By, Bz, lapse = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        for i = 1, 3 do
          shift[i] = 0.0
          for j = 1, 3 do
            spatial_metric[i][j] = 0.0
          end  
        end

        excision = -1.0
      else
        excision = 1.0
      end
    
      return Dx, Dy, Dz, Bx, By, Bz, 0.0, 0.0,
        lapse,
        shift[1], shift[2], shift[3],
        spatial_metric[1][1], spatial_metric[1][2], spatial_metric[1][3],
        spatial_metric[2][1], spatial_metric[2][2], spatial_metric[2][3],
        spatial_metric[3][1], spatial_metric[3][2], spatial_metric[3][3],
        excision,
        0.0,
        x, y, 0.0
    end,

    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (x-direction).
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (y-direction).
  }
}

-- Run application.
momentApp:run()
