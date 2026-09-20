local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

nl = 1.0 -- Left number density.
Tl = 0.1 -- Left temperature.
V_r_drift_l = 0.0 -- Left drift velocity (radial direction).
V_phi_drift_l = 0.0 -- Left drift velocity (angular direction).

nr = 1.0 -- Right number density.
Tr = 0.1 -- Right temperature.
V_r_drift_r = 0.0 -- Right drift velocity (radial direction).
V_phi_drift_r = 0.0 -- Right drift velocity (angular direction).

nb = 1.0 -- Background number density.
Tb = 0.1 -- Background temperature.
V_r_drift_b = 0.0 -- Background drift velocity (radial direction).
V_phi_drift_b = 0.0 -- Background drift velocity (angular direction).

vt = 1.0 -- Thermal velocity.
nu = 15000.0 -- Collision frequency.

-- Simulation parameters.
Nr = 32 -- 124 -- Cell count (configuration space: radial direction).
Nvr = 48 -- 96 -- Cell count (velocity space: radial direction).
Nvphi = 48 -- 96 --Cell count (velocity space: angular direction).
Lr = 9.5 -- Domain size (configuration space: radial direction).
vr_max = 1.5 * vt -- Domain boundary (velocity space: radial direction).
vphi_max = 1.5 * vt -- Domain boundary (velocity space: angular direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 50.0 -- Final simulation time.
num_frames = 15 -- 100 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

massBH = 0.3 -- Mass of the Kerr-Schild black hole
spinBH = 0.295 -- Spin parameter, a = J/M, of the Kerr-Schild black hole

r_inner = 1.2 -- Ring inner radius.
r_outer = 10.0 -- Ring outer radius.

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.3 },
  upper = { 0.3 + Lr },
  cells = { Nr },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (phi-dir).

  -- Vlasov Geometry
  geom = Vlasov.Geom.new {
    
    -- Black hole parameters
    massBH = massBH,
    spinBH = spinBH,

    -- Use the vierbein inputs to contruct the Poisson Tensor
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.GR_KS_rphi
    
  },

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.TriadGR,
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { -vr_max, -vphi_max },
    upper = { vr_max, vphi_max },
    cells = { Nvr, Nvphi },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local r = xn[1]
          local theta = math.pi/2.0 -- r-phi plane
          local phi = 0.0 -- angle

          local n = 0.0
          if r > r_inner and r < r_outer then
            if phi > (math.pi) then
              n = nl -- Fluid-stationary frame density (left ring).
            else
              n = nr -- Fluid-stationary frame density (right ring).
            end
          else
            n = nb -- Fluid-stationary frame density (background).
          end
          
          local rho = math.sqrt(r * r + spinBH * spinBH * math.cos(theta) * math.cos(theta) )
          local metric_det = rho * math.sqrt(2*massBH*r + rho * rho) * math.sin(theta)          

          return metric_det * n
        end,
        temperatureInit = function (t, xn)
          local r = xn[1]
          local theta = math.pi/2.0 -- r-phi plane
          local phi = 0.0 -- angle

          local T = 0.0
         if r > r_inner and r < r_outer then
            if phi > (math.pi) then
              T = Tl -- Fluid-stationary frame Temperature  (left ring).
            else
              T = Tr -- Fluid-stationary frame Temperature  (right ring).
            end
          else
            T = Tb -- Fluid-stationary frame Temperature  (background).
          end

          return T
        end,
        driftVelocityInit = function (t, xn)
          local r = xn[1]
          local theta = math.pi/2.0 -- r-phi plane
          local phi = 0.0 -- angle

          local V_r_drift = 0.0
          local V_phi_drift = 0.0
          if r > r_inner and r < r_outer then
            if phi > (math.pi) then
              V_r_drift = V_r_drift_l -- Fluid velocity r (left ring).
              V_phi_drift = V_phi_drift_l ---- Fluid velocity. phi (left ring).
            else
              V_r_drift = V_r_drift_r -- Fluid velocity r (right ring).
              V_phi_drift = V_phi_drift_r ---- Fluid velocity. phi (right ring).
            end
          else
            V_r_drift = V_r_drift_b -- Fluid velocity r (background).
            V_phi_drift = V_phi_drift_b ---- Fluid velocity. phi (background).
          end

          return V_r_drift, V_phi_drift
        end,

        correctAllMoments = true,
        iterationEpsilon = 0.0,
        maxIterations = 0,
        useLastConverged = false
      }
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcAbsorb
      },
      upper = {
        type = G0.SpeciesBc.bcAbsorb
      }
    },

    correct = {
      correctAllMoments = true,
      iterationEpsilon = 1e-12,
      maxIterations = 100,
      useLastConverged = false
    },


    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments, G0.Moment.EnergyMoment }
  },

  skipField = true,

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = 1.0, mu0 = 1.0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0,

    isStatic = true
  }
}

vlasovApp:run()
