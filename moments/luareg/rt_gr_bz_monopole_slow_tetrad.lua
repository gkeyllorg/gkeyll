-- 2D Blandford-Znajek magnetosphere problem for a slowly-rotating (Kerr) black hole, for the general relativistic Maxwell equations in the tetrad basis.
-- Input parameters describe a (purely radial) monopole magnetic field surrounding a slowly-rotating black hole.
-- Based on the perturbative analytical solution for force-free electrodynamics presented in the article:
-- R. D. Blandford and R. L. Znajek (1977), "Electromagnetic extraction of energy from Kerr black holes",
-- Monthly Notices of the Royal Astronomical Society, Volume 179 (3): 433-456.
-- https://academic.oup.com/mnras/article/179/3/433/962905

local Moments = G0.Moments
local GRMaxwellTetrad = G0.Moments.Eq.GRMaxwellTetrad
local BlackHole = G0.Moments.Spacetime.BlackHole

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
light_speed = 1.0 -- Speed of light.
e_fact = 0.0 -- Factor of speed of light for electric field correction.
b_fact = 0.0 -- Factor of speed of light for magnetic field correction.

B0 = 1.0 -- Reference magnetic field strength.

-- Spacetime parameters (using geometric units).
mass = 1.0 -- Mass of the black hole.
spin = -0.1 -- Spin of the black hole.

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
    equation = GRMaxwellTetrad.new {
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

      local r = math.sqrt((x * x) + (y * y))
      local phi = 0.5 * pi
    
      local theta = 0.0;
      if math.abs(y) < math.pow(10.0, -6.0) then
        theta = 0.5 * pi;
      else
        theta = math.atan(x / y);
      end

      local lapse = BlackHole.lapseFunction(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)
      local shift = BlackHole.shiftVector(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)
      local spatial_metric = BlackHole.spatialMetricTensor(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)
      local spatial_det = BlackHole.spatialMetricDeterminant(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)
      local in_excision_region = BlackHole.excisionRegion(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)

      local B_r = B0 * math.sin(theta) / math.sqrt(spatial_det)

      local Dx = 0.0 -- Total electric field (x-direction).
      local Dy = 0.0 -- Total electric field (y-direction).
      local Dz = 0.0 -- Total electric field (z-direction).

      local Bx = math.sin(theta) * math.sin(phi) * B_r -- Total magnetic field (x-direction).
      local By = math.sin(theta) * math.cos(phi) * B_r -- Total magnetic field (y-direction).
      local Bz = math.cos(theta) * B_r -- Total magnetic field (z-direction).

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
