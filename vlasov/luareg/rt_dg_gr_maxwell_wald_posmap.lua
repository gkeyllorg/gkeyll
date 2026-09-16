local Vlasov = G0.Vlasov
-- Primary reference: S. S. Komissarov, Electrodynamics of black hole magnetospheres, 2004 (Eqns 101/102)
-- Axisymmetry in phi, Schwarzschild Wald solution
-- Also see: K. Parfery et al. Introducing PHAEDRA: a new spectral code for simulations of relativistic magnetospheres
-- Also see: Wald solution referenced from Section 3.2 of the Tetrad-First paper: https://arxiv.org/pdf/2410.02549

-- Mathematical constants (dimensionless).
pi = math.pi

-- Simulation parameters.
Nr = 48 -- 192 -- Cell count (r-direction).
Ntheta = 192 -- Cell count (theta-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 5.0 -- Final simulation time.
num_frames = 1 -- 30 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

massBH = 1.0 -- Mass of the black hole
spinBH = 0.0 -- Spin parameter, a = J/M, of the black hole ( Kerr-Schild, coordinates , |a| <= M).

B_0 = 1.0 -- Magnetic field amplitude 

vlasovApp = Vlasov.App.new {

  -- Vlasov Geometry
  geom = Vlasov.Geom.new {
    
    -- Black hole parameters
    massBH = massBH,
    spinBH = spinBH,

    -- Specify the geometry rtheta choice
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.GR_KS_rtheta
    
  },

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 1.5 * massBH, 0 },
  upper = { 30.0 * massBH, math.pi },
  cells = { Nr, Ntheta },
  cflFrac = cfl_frac,

  mapc2pPos = {
    -- r mapping: power-law refinement toward the horizon (finest cells at r = 1.5 M).
    {
      pmap = function (t, xn)
        local rc = xn[1]
        local r_lo, r_up = 1.5 * massBH, 30.0 * massBH
        local xi = (rc - r_lo) / (r_up - r_lo)
        return r_lo + (r_up - r_lo) * math.pow(xi, 1.5)
      end
    },
  },
    
  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (r-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = {  }, -- Periodic directions (x- and theta-directions only).

  -- Field.
  field = Vlasov.Field.new {

    -- Use GR field ID
    fieldID = G0.FieldModel.GR,

    -- Speed factors
    epsilon0 = 1.0,
    mu0 = 1.0,
    elcErrorSpeedFactor = 1.0, -- chi = c*elcErrorSpeedFactor = 1.
    mgnErrorSpeedFactor = 1.0, -- gamma = c*mgnErrorSpeedFactor = 1.
    K_phi = 1.0, -- Damping Constant (electric field).
    K_psi = 1.0, -- Damping Constant (magnetic field).

    -- Initial conditions function.
    init = function (t, xn)
      local r, theta = xn[1], xn[2]

      -- D^i, and B^i (contravaraint components)
      -- Even though E_i is non-zero, the vacuum consituitve relations give 
      -- a non-zero D^\phi component
      local Dr = 0.0 -- Total electric field (r-direction).
      local Dtheta = 0.0 -- Total electric field (theta-direction).
      local Dphi = - 2.0 * massBH * B_0 / ( r * r * math.sqrt( 1.0 + 2.0 * massBH / r ) )  -- Total electric field (phi-direction).
    
      local Br = - B_0 * math.cos(theta) / ( math.sqrt( 1.0 + 2.0 * massBH / r ) )   -- Total magnetic field (r-direction).
      local Btheta = B_0 * math.sin(theta) / ( r * math.sqrt( 1.0 + 2.0 * massBH / r ) ) -- Total magnetic field (theta-direction).
      local Bphi = 0.0 -- Total magnetic field (phi-direction).

      -- Hand off the primative variables (U^\xi)
      return Dr, Dtheta, Dphi, Br, Btheta, Bphi, 0.0, 0.0
    end,

    -- Copy boundary conditions in r are sufficient. Theta requries theta-pole BCs
    bcx = { G0.FieldBc.bcFixedFunc, G0.FieldBc.bcFixedFunc }, -- boundary conditions (r-direction).
    bcy = { G0.FieldBc.bcThetaPole, G0.FieldBc.bcThetaPole } -- boundary conditions (theta-direction).
  }
}

-- Run application.
vlasovApp:run()