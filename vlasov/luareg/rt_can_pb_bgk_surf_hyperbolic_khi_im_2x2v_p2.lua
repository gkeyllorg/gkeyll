local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

nl = 2.0 -- Left/inner number density.
Pl = 2.5 -- Left/inner pressure.
V_theta_drift_l = -0.5 -- Left/inner drift velocity (azimuthal angular direction).
V_z_drift_l = 0.0 -- Left/inner drift velocity (z direction).

nr = 1.0 -- Right/outer number density.
Pr = 2.5 -- Right/outer pressure.
V_theta_drift_r = 0.5 -- Right/outer drift velocity (azimuthal angular direction).
V_z_drift_r = 0.0 -- Right/outer drift velocity (z direction).

vt = 1.0 -- Thermal velocity.
nu = 15000.0 -- Collision frequency.

-- Simulation parameters.
Ntheta = 32 -- Cell count (configuration space: azimuthal angular direction).
Nz = 32 -- Cell count (configuration space: z direction).
Nvtheta = 8 -- Cell count (velocity space: azimuthal angular direction).
Nvz = 8 -- Cell count (velocity space: z direction).
Ltheta = 2.0 * pi -- Domain size (configuration space: azimuthal angular direction).
Lz = 1.0 -- Domain size (configuration space: z direction).
vtheta_max = 8.0 * vt -- Domain boundary (velocity space: azimuthal angular direction).
vz_max = 8.0 * vt -- Domain boundary (velocity space: z direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.001 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

R = 1.0 -- Radius of the sphere.
z_midplane = 0.0 -- azimuthal angular z_midplane location.
z_loc = 0.5 -- azimuthal angular boundary location designating jump in quantities.

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0, -Lz },
  upper = { Ltheta, Lz },
  cells = { Ntheta, Nz },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (theta direction only).

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.CanonicalPB,
    charge = charge, mass = mass,

    hamiltonian = function (t, xn)
      local q_theta, q_z, p_theta_dot, p_z_dot = xn[1], xn[2], xn[3], xn[4]

      local inv_metric_theta_theta = 1.0 / ( R * R + q_z * q_z) -- Inverse metric tensor (azimuthal-azimuthal component).
      local inv_metric_theta_z = 0.0 -- Inverse metric tensor (azimuthal-z component).
      local inv_metric_z_z = ( R * R + q_z * q_z) / ( R * R + 2.0 * q_z * q_z ) -- Inverse metric tensor (z-z component).
   
      local hamiltonian = (0.5 * inv_metric_theta_theta * p_theta_dot * p_theta_dot) + (0.5 * (2.0 * inv_metric_theta_z * p_theta_dot * p_z_dot)) +
        (0.5 * inv_metric_z_z * p_z_dot * p_z_dot) -- Canonical Hamiltonian.

      return hamiltonian
    end,
    inverseMetric = function (t, xn)
      local q_theta = xn[1]
      local q_z = xn[2]

      local inv_metric_theta_theta = 1.0 / ( R * R + q_z * q_z) -- Inverse metric tensor (azimuthal-azimuthal component).
      local inv_metric_theta_z = 0.0 -- Inverse metric tensor (azimuthal-z component).
      local inv_metric_z_z = ( R * R + q_z * q_z) / ( R * R + 2.0 * q_z * q_z ) -- Inverse metric tensor (z-z component).

      return inv_metric_theta_theta, inv_metric_theta_z, inv_metric_z_z
    end,
    metric = function (t, xn)
      local q_theta = xn[1]
      local q_z = xn[2]

      local metric_theta_theta = R * R + q_z * q_z -- Metric tensor (azimuthal-azimuthal component).
      local metric_theta_z = 0.0 -- Metric tensor (azimuthal-z component).
      local metric_z_z = ( R * R + 2.0 * q_z * q_z ) / ( R * R + q_z * q_z ) -- Metric tensor (z-z component).

      return metric_theta_theta, metric_theta_z, metric_z_z
    end,
    metricDeterminant = function (t, xn)
      local z = xn[2]

      local metric_det = math.sqrt( R * R + 2.0 * z * z ) -- Metric tensor determinant.

      return metric_det
    end,
    
    -- Velocity space grid.
    lower = { -vtheta_max, -vz_max },
    upper = { vtheta_max, vz_max },
    cells = { Nvtheta, Nvz },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local theta = xn[1]
          local z = xn[2]

          local n = 0.0
          if math.abs(z - z_midplane) < z_loc then
            n = nl -- Total number density (left/inner).
          else
            n = nr -- Total number density (right/outer).
          end
        
          local metric_det = math.sqrt( R * R + 2.0 * z * z )

          return metric_det * n
        end,
        temperatureInit = function (t, xn)
          local theta = xn[1]
          local z = xn[2]

          local T = 0.0
          if math.abs(z - z_midplane) < z_loc then
            T = Pl/nl -- Isotropic temperature (left/inner).
          else
            T = Pr/nr -- Isotropic temperature (right/outer).
          end

          return T
        end,
        driftVelocityInit = function (t, xn)
          local theta = xn[1]
          local z = xn[2]

          local V_theta_drift = 0.0
          local V_z_drift = 0.0
        
          if math.abs(z - z_midplane) < z_loc then
            V_theta_drift = V_theta_drift_l -- azimuthal angular drift velocity (left/inner).
            V_z_drift = V_z_drift_l -- z drift velocity (left/inner).
          else
            V_theta_drift = V_theta_drift_r -- azimuthal angular drift velocity (right/outer).
            V_z_drift = V_z_drift_r -- z drift velocity (right/outer).
          end

          math.randomseed(0)

          -- Initalize noise
          local alpha = 1.0e-2
          local k_theta = 2.0 * pi / Ltheta
          local k_z = 2.0 * pi / Lz

          for i = 0, 16 do
            for j = 0, 16 do
              V_theta_drift = V_theta_drift + alpha * math.random() * math.sin(i * k_theta * theta + j * k_z * z + 2.0 * pi * math.random())
              V_z_drift = V_z_drift + alpha * math.random() * math.sin(i * k_theta * theta + j * k_z * z + 2.0 * pi * math.random())
            end
          end

          return V_theta_drift, V_z_drift
        end,

        correctAllMoments = true,
        iterationEpsilon = 0.0,
        maxIterations = 0,
        useLastConverged = false
      }
    },

    collisions = {
      collisionID = G0.Collisions.BGK,

      selfNu = function (t, xn)
        return nu -- Collision frequency.
      end,
      
      useImplicitCollisionScheme = true
    },

    correct = {
      outputfLTE = true,
      correctAllMoments = true,
      iterationEpsilon = 1e-12,
      maxIterations = 100,
      useLastConverged = false
    },

    bcy = {
      lower = {
        type = G0.SpeciesBc.bcReflect
      },
      upper = {
        type = G0.SpeciesBc.bcReflect
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments }
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