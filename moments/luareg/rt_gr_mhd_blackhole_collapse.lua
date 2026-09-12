local Moments = G0.Moments
local GRMHD = G0.Moments.Eq.GRMHD
local Minkowski = G0.Moments.Spacetime.Minkowski
local BlackHole = G0.Moments.Spacetime.BlackHole

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.

rhob = 1.0 -- Background fluid mass density.
omega_b = 0.2 -- Background fluid angular velocity.
pb = 1.0 -- Background fluid pressure.

rhos = 10.0 -- Star fluid mass density.
omega_s = 0.0 -- Star fluid angular velocity.
ps = 10.0 -- Star fluid pressure.

B0 = 0.05 -- Reference magnetic field strength.

light_speed = 1.0 -- Speed of light.
b_fact = 0.8 -- Factor of speed of light for magnetic field correction.

-- Predicted spacetime parameters (using geometric units).
mass = 0.3 -- Predicted mass of the black hole.
spin = 0.9 -- Predicted spin of the black hole.

pos_x = 2.5 -- Predicted position of the black hole (x-direction).
pos_y = 2.5 -- Predicted position of the black hole (y-direction).
pos_z = 0.0 -- Predicted position of the black hole (z-direction).

-- Simulation parameters.
Nx = 256 -- Cell count (x-direction).
Ny = 256 -- Cell count (y-direction).
Lx = 5.0 -- Domain size (x-direction).
Ly = 5.0 -- Domain size (y-direction).
cfl_frac = 0.9 -- CFL coefficient.

spacetime_gauge = G0.SpacetimeGauge.BlackHoleCollapse -- Spacetime gauge choice.
reinit_freq = 100 -- Spacetime reinitialization frequency.

t_end = 2.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

r_star = 1.5 -- Star radius.

momentApp = Moments.App.new {
  
  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x- and y-directions).
  
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Fluid.
  fluid = Moments.Species.new {
    equation = GRMHD.new {
      gasGamma = gas_gamma,
      lightSpeed = light_speed,
      mgnErrorSpeedFactor = b_fact,
      blackHoleParameters = {
        mass = mass,
        spin = spin,
        posX = pos_x,
        posY = pos_y,
        posZ = pos_z
      },
      spacetimeGauge = spacetime_gauge,
      reinitFreq = reinit_freq
    },

    hasGRMHD = true,
    GRMHDGasGamma = gas_gamma,
  
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]
      
      local rho = 0.0
      local u = 0.0
      local v = 0.0
      local p = 0.0

      local r = math.sqrt((x - (0.5 * Lx)) * (x - (0.5 * Lx)) + (y - (0.5 * Ly)) * (y - (0.5 * Ly)))
    
      if (r <= r_star) then
        rho = rhos -- Fluid mass density (star).
        u = -omega_s * (y - (0.5 * Ly)) -- Fluid velocity (star).
        v = omega_s * (x - (0.5 * Lx))
        p = ps -- Fluid pressure (star).
      else
        rho = rhob -- Fluid mass density (background).
        u = -omega_b * (y - (0.5 * Ly)) -- Fluid velocity (background).
        v = omega_b * (x - (0.5 * Lx))
        p = pb -- Fluid pressure (background).
      end

      local lapse = Minkowski.lapseFunction(0.0, x, y, 0.0)
      local shift = Minkowski.shiftVector(0.0, x, y, 0.0)
      local spatial_metric = Minkowski.spatialMetricTensor(0.0, x, y, 0.0)
      local spatial_det = Minkowski.spatialMetricDeterminant(0.0, x, y, 0.0)
      local extrinsic_curvature = Minkowski.extrinsicCurvatureTensor(0.0, x, y, 0.0, 1.0, 1.0, 1.0)
      local in_excision_region = Minkowski.excisionRegion(0.0, x, y, 0.0)

      local lapse_der = Minkowski.lapseFunctionDer(0.0, x, y, 0.0, 1.0, 1.0, 1.0)
      local shift_der = Minkowski.shiftVectorDer(0.0, x, y, 0.0, 1.0, 1.0, 1.0)
      local spatial_metric_der = Minkowski.spatialMetricTensorDer(0.0, x, y, 0.0, 1.0, 1.0, 1.0)

      local vel = { u, v, 0.0 }
      local v_sq = 0.0

      for i = 1, 3 do
        for j = 1, 3 do
          v_sq = v_sq + (spatial_metric[i][j] * vel[i] * vel[j])
        end
      end

      local W = 1.0 / math.sqrt(1.0 - v_sq)
      if v_sq > 1.0 - math.pow(10.0, -8.0) then
        W = 1.0 / math.sqrt(1.0 - pow(10.0, -8.0))
      end

      local mag = { 0.0, B0, 0.0 }

      local cov_mag = { 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        for j = 1, 3 do
          cov_mag[i] = cov_mag[i] + (spatial_metric[i][j] * mag[j])
        end
      end

      local cov_vel = { 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        for j = 1, 3 do
          cov_vel[i] = cov_vel[i] + (spatial_metric[i][j] * vel[j])
        end
      end

      local b0 = 0.0
      for i = 1, 3 do
        b0 = b0 + (W * mag[i] * (cov_vel[i] / lapse))
      end

      local spacetime_vel = { W / lapse, 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        spacetime_vel[i + 1] = (W * vel[i]) - (shift[i] * (W / lapse))
      end

      local b = { 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        b[i] = (mag[i] + (lapse * b0 * spacetime_vel[i + 1])) / W
      end

      local b_sq = 0.0
      for i = 1, 3 do
        b_sq = b_sq + ((mag[i] * cov_mag[i]) / (W * W))
      end
      b_sq = b_sq + (((lapse * lapse) * (b0 * b0)) / (W * W))

      local cov_b = { 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        for j = 1, 3 do
          cov_b[i] = cov_b[i] + (spatial_metric[i][j] * b[j])
        end
      end

      local h_star = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq / rho)
      local p_star = p + (0.5 * b_sq)

      local rho_rel = math.sqrt(spatial_det) * rho * W -- Fluid relativistic mass density.
      local mom_x = math.sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[1]) - (lapse * b0 * cov_b[1])) -- Fluid momentum density (x-direction).
      local mom_y = math.sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[2]) - (lapse * b0 * cov_b[2])) -- Fluid momentum density (y-direction).
      local mom_z = 0.0 -- Fluid momentum density (z-direction).
      local Etot = math.sqrt(spatial_det) * ((rho * h_star * (W * W)) - p_star - ((lapse * lapse) * (b0 * b0)) - (rho * W)) -- Fluid total energy density.

      local Bx_rel = 0.0 -- Fluid relativistic magnetic field (x-direction).
      local By_rel = math.sqrt(spatial_det) * B0 -- Fluid relativistic magnetic field (y-direction).
      local Bz_rel = 0.0 -- Fluid relativistic magnetic field (z-direction).

      local excision = 0.0
      if in_excision_region then
        rho_rel, mom_x, mom_y, mom_z, Etot, Bx_rel, By_rel, Bz_rel, lapse = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        for i = 1, 3 do
          shift[i] = 0.0
          lapse_der[i] = 0.0
          for j = 1, 3 do
            spatial_metric[i][j] = 0.0
            extrinsic_curvature[i][j] = 0.0
            shift_der[i][j] = 0.0
            for k = 1, 3 do
              spatial_metric_der[i][j][k] = 0.0
            end
          end  
        end
        
        excision = -1.0
      else
        excision = 1.0
      end
    
      return rho_rel, mom_x, mom_y, mom_z, Etot, Bx_rel, By_rel, Bz_rel, 0.0,
        lapse,
        shift[1], shift[2], shift[3],
        spatial_metric[1][1], spatial_metric[1][2], spatial_metric[1][3],
        spatial_metric[2][1], spatial_metric[2][2], spatial_metric[2][3],
        spatial_metric[3][1], spatial_metric[3][2], spatial_metric[3][3],
        extrinsic_curvature[1][1], extrinsic_curvature[1][2], extrinsic_curvature[1][3],
        extrinsic_curvature[2][1], extrinsic_curvature[2][2], extrinsic_curvature[2][3],
        extrinsic_curvature[3][1], extrinsic_curvature[3][2], extrinsic_curvature[3][3],
        excision,
        lapse_der[1], lapse_der[2], lapse_der[3],
        shift_der[1][1], shift_der[1][2], shift_der[1][3],
        shift_der[2][1], shift_der[2][2], shift_der[2][3],
        shift_der[3][1], shift_der[3][2], shift_der[3][3],
        spatial_metric_der[1][1][1], spatial_metric_der[1][1][2], spatial_metric_der[1][1][3],
        spatial_metric_der[1][2][1], spatial_metric_der[1][2][2], spatial_metric_der[1][2][3],
        spatial_metric_der[1][3][1], spatial_metric_der[1][3][2], spatial_metric_der[1][3][3],
        spatial_metric_der[2][1][1], spatial_metric_der[2][1][2], spatial_metric_der[2][1][3],
        spatial_metric_der[2][2][1], spatial_metric_der[2][2][2], spatial_metric_der[2][2][3],
        spatial_metric_der[2][3][1], spatial_metric_der[2][3][2], spatial_metric_der[2][3][3],
        spatial_metric_der[3][1][1], spatial_metric_der[3][1][2], spatial_metric_der[3][1][3],
        spatial_metric_der[3][2][1], spatial_metric_der[3][2][2], spatial_metric_der[3][2][3],
        spatial_metric_der[3][3][1], spatial_metric_der[3][3][2], spatial_metric_der[3][3][3],
        0.0,
        x, y, 0.0
    end,

    evolve = true, -- Evolve species?
    forceLowOrderFlux = true, -- Use Lax fluxes.
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (x-direction).
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (y-direction).
  }
}

-- Run application.
momentApp:run()
