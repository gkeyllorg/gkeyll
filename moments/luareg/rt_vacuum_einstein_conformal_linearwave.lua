local Moments = G0.Moments
local VacuumEinsteinConformal = G0.Moments.Eq.VacuumEinsteinConformal
local Minkowski = G0.Moments.Spacetime.Minkowski

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
amp = math.pow(10.0, -8.0) -- Wave amplitude.

-- Evolution parameters.
excision_threshold = 0.3 -- Excision threshold (lapse).
spacetime_slicing = G0.SpacetimeSlicing.Harmonic -- Spacetime slicing condition.
spacetime_evolution = G0.SpacetimeEvolution.Einstein -- Spacetime evolution system.

-- Simulation parameters.
Nx = 50 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
cfl_frac = 0.95 -- CFL coefficient.

t_end = 1000.0 -- Final simulation time.
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
  lower = { -0.5 * Lx },
  upper = { 0.5 * Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x- and y-directions).
  
  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Fluid.
  fluid = Moments.Species.new {
    equation = VacuumEinsteinConformal.new {
      excisionThreshold = excision_threshold,
      spacetimeSlicing = spacetime_slicing,
      spacetimeEvolution = spacetime_evolution
    },

    hasVacuumEinsteinConformal = true,
    vacuumEinsteinConformalExcisionThreshold = excision_threshold,
    vacuumEinsteinConformalSpacetimeSlicing = spacetime_slicing,
    vacuumEinsteinConformalSpacetimeEvolution = spacetime_evolution,
  
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local conformal_lapse = Minkowski.lapseFunction(0.0, x, 0.0, 0.0)
      local conformal_shift = Minkowski.shiftVector(0.0, x, 0.0, 0.0)
      local conformal_spatial_metric = Minkowski.spatialMetricTensor(0.0, x, 0.0, 0.0)
      local inv_conformal_spatial_metric = Minkowski.invSpatialMetricTensor(0.0, x, 0.0, 0.0)
      local conformal_spatial_det = Minkowski.spatialMetricDeterminant(0.0, x, 0.0, 0.0)
      local conformal_extrinsic_curvature = Minkowski.extrinsicCurvatureTensor(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)
      local in_excision_region = Minkowski.excisionRegion(0.0, x, 0.0, 0.0)

      local conformal_lapse_der = Minkowski.lapseFunctionDer(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)
      local conformal_shift_der = Minkowski.shiftVectorDer(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)
      local conformal_spatial_metric_der = Minkowski.spatialMetricTensorDer(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)

      local conformal_fact = Minkowski.conformalFactor(0.0, x, 0.0, 0.0)
      local bssn_conformal_fact = Minkowski.bssnConformalFactor(0.0, x, 0.0, 0.0)

      local conformal_fact_der = Minkowski.conformalFactorDer( 0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)
      local bssn_conformal_fact_der = Minkowski.bssnConformalFactorDer(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)
      local bssn_conformal_fact_der2 = Minkowski.bssnConformalFactorDer2(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)

      b = amp * math.sin(2.0 * pi * x)
      conformal_spatial_metric[2][2] = (1.0 + b) / (conformal_fact * conformal_fact * conformal_fact * conformal_fact)
      conformal_spatial_metric[3][3] = (1.0 - b) / (conformal_fact * conformal_fact * conformal_fact * conformal_fact)

      conformal_extrinsic_curvature[2][2] = -amp * pi * math.cos(2.0 * pi * x)
      conformal_extrinsic_curvature[3][3] = amp * pi * math.cos(2.0 * pi * x)

      conformal_spatial_metric_der[1][2][2] = 2.0 * amp * pi * math.cos(2.0 * pi * x)
      conformal_spatial_metric_der[1][3][3] = -2.0 * amp * pi * math.cos(2.0 * pi * x)

      conformal_spatial_det = 1.0 - (b * b)
      inv_conformal_spatial_metric[2][2] = (1.0 / (1.0 + b)) * (conformal_fact * conformal_fact * conformal_fact * conformal_fact)
      inv_conformal_spatial_metric[3][3] = (1.0 / (1.0 + b)) * (conformal_fact * conformal_fact * conformal_fact * conformal_fact)

      for i = 1, 3 do
        for j = 1, 3 do
          for k = 1, 3 do
            conformal_spatial_metric_der[i][j][k] = 0.5 * conformal_spatial_metric_der[i][j][k]
            conformal_spatial_metric_der[i][j][k] = conformal_spatial_metric_der[i][j][k] / (conformal_fact * conformal_fact * conformal_fact * conformal_fact)
            conformal_spatial_metric_der[i][j][k] = conformal_spatial_metric_der[i][j][k] - (2.0 * conformal_fact_der[i] * conformal_spatial_metric[j][k]);

          end

        conformal_shift_der[i][j] = 0.5 * conformal_shift_der[i][j]
        end
      end

      for i = 1, 3 do
        conformal_lapse_der[i] = conformal_lapse_der[i] / conformal_lapse
      end

      local conformal_spatial_metric_der_raised1 = { }
      for i = 1, 3 do
        conformal_spatial_metric_der_raised1[i] = { }
        for j = 1, 3 do
          conformal_spatial_metric_der_raised1[i][j] = { }
        end
      end
      for i = 1, 3 do
        for j = 1, 3 do
          for k = 1, 3 do
            conformal_spatial_metric_der_raised1[k][i][j] = 0.0
        
            for l = 1, 3 do
              conformal_spatial_metric_der_raised1[k][i][j] = conformal_spatial_metric_der_raised1[k][i][j] + inv_conformal_spatial_metric[k][l] * conformal_spatial_metric_der[l][i][j]
            end
          end
        end
      end
  
      local conformal_spatial_metric_der_raised3 = { }
      for i = 1, 3 do
        conformal_spatial_metric_der_raised3[i] = { }
        for j = 1, 3 do
          conformal_spatial_metric_der_raised3[i][j] = { }
        end
      end
      for i = 1, 3 do
        for j = 1, 3 do
          for k = 1, 3 do
            conformal_spatial_metric_der_raised3[i][j][k] = 0.0

            for l = 1, 3 do
              conformal_spatial_metric_der_raised3[i][j][k] = conformal_spatial_metric_der_raised3[i][j][k] + inv_conformal_spatial_metric[l][k] * conformal_spatial_metric_der[i][j][l]
            end
          end
        end
      end

      local conformal_aux_vect = { }
      for i = 1, 3 do
        conformal_aux_vect[i] = 0.0

        for s = 1, 3 do
          conformal_aux_vect[i] = conformal_aux_vect[i] + conformal_spatial_metric_der_raised3[i][s][s]
          conformal_aux_vect[i] = conformal_aux_vect[i] - conformal_spatial_metric_der_raised1[s][s][i]
        end
        
        conformal_aux_vect[i] = conformal_aux_vect[i] - 4.0 * conformal_fact_der[i]
      end

      local excision = 0.0
      if in_excision_region then
        conformal_lapse = 0.0
        conformal_fact = 0.0
        bssn_conformal_fact = 0.0
        for i = 1, 3 do
          conformal_shift[i] = 0.0
          conformal_lapse_der[i] = 0.0
          conformal_aux_vect[i] = 0.0
          conformal_fact_der[i] = 0.0
          bssn_conformal_fact_der[i] = 0.0
          for j = 1, 3 do
            conformal_spatial_metric[i][j] = 0.0
            conformal_extrinsic_curvature[i][j] = 0.0
            conformal_shift_der[i][j] = 0.0
            bssn_conformal_fact_der2[i][j] = 0.0
            for k = 1, 3 do
              conformal_spatial_metric_der[i][j][k] = 0.0
            end
          end  
        end
        
        excision = -1.0
      else
        excision = 1.0
      end
    
      return conformal_spatial_metric[1][1], conformal_spatial_metric[1][2], conformal_spatial_metric[1][3],
        conformal_spatial_metric[2][1], conformal_spatial_metric[2][2], conformal_spatial_metric[2][3],
        conformal_spatial_metric[3][1], conformal_spatial_metric[3][2], conformal_spatial_metric[3][3],
        conformal_lapse,
        conformal_extrinsic_curvature[1][1], conformal_extrinsic_curvature[1][2], conformal_extrinsic_curvature[1][3],
        conformal_extrinsic_curvature[2][1], conformal_extrinsic_curvature[2][2], conformal_extrinsic_curvature[2][3],
        conformal_extrinsic_curvature[3][1], conformal_extrinsic_curvature[3][2], conformal_extrinsic_curvature[3][3],
        conformal_spatial_metric_der[1][1][1], conformal_spatial_metric_der[1][1][2], conformal_spatial_metric_der[1][1][3],
        conformal_spatial_metric_der[1][2][1], conformal_spatial_metric_der[1][2][2], conformal_spatial_metric_der[1][2][3],
        conformal_spatial_metric_der[1][3][1], conformal_spatial_metric_der[1][3][2], conformal_spatial_metric_der[1][3][3],
        conformal_spatial_metric_der[2][1][1], conformal_spatial_metric_der[2][1][2], conformal_spatial_metric_der[2][1][3],
        conformal_spatial_metric_der[2][2][1], conformal_spatial_metric_der[2][2][2], conformal_spatial_metric_der[2][2][3],
        conformal_spatial_metric_der[2][3][1], conformal_spatial_metric_der[2][3][2], conformal_spatial_metric_der[2][3][3],
        conformal_spatial_metric_der[3][1][1], conformal_spatial_metric_der[3][1][2], conformal_spatial_metric_der[3][1][3],
        conformal_spatial_metric_der[3][2][1], conformal_spatial_metric_der[3][2][2], conformal_spatial_metric_der[3][2][3],
        conformal_spatial_metric_der[3][3][1], conformal_spatial_metric_der[3][3][2], conformal_spatial_metric_der[3][3][3],
        conformal_lapse_der[1], conformal_lapse_der[2], conformal_lapse_der[3],
        conformal_aux_vect[1], conformal_aux_vect[2], conformal_aux_vect[3],
        conformal_shift[1], conformal_shift[2], conformal_shift[3],
        conformal_shift_der[1][1], conformal_shift_der[1][2], conformal_shift_der[1][3],
        conformal_shift_der[2][1], conformal_shift_der[2][2], conformal_shift_der[2][3],
        conformal_shift_der[3][1], conformal_shift_der[3][2], conformal_shift_der[3][3],
        bssn_conformal_fact,
        bssn_conformal_fact_der[1], bssn_conformal_fact_der[2], bssn_conformal_fact_der[3],
        bssn_conformal_fact_der2[1][1], bssn_conformal_fact_der2[1][2], bssn_conformal_fact_der2[1][3],
        bssn_conformal_fact_der2[2][1], bssn_conformal_fact_der2[2][2], bssn_conformal_fact_der2[2][3],
        bssn_conformal_fact_der2[3][1], bssn_conformal_fact_der2[3][2], bssn_conformal_fact_der2[3][3]
    end,

    evolve = true, -- Evolve species?
    forceLowOrderFlux = false -- Use HLL fluxes.
  }
}

-- Run application.
momentApp:run()
