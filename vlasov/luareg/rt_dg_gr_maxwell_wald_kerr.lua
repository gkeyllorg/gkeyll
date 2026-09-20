local Vlasov = G0.Vlasov
-- Primary reference: S. S. Komissarov, Electrodynamics of black hole magnetospheres, 2004 (Eqns 101)
-- Axisymmetry in phi, Kerr Wald solution
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
spinBH = 0.95 -- Spin parameter, a = J/M, of the black hole ( Kerr-Schild, coordinates , |a| <= M).

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
  lower = { 0.9, 0 },
  upper = { 5.0, math.pi },
  cells = { Nr, Ntheta },
  cflFrac = cfl_frac,
    
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
      -- Compact geometry factors.
      local a = spinBH
      local M = massBH
      local ct = math.cos(theta)
      local st = math.sin(theta)
      local ct2 = ct*ct
      local st2 = st*st
      local a2 = a*a
      local a4 = a2*a2
      local a6 = a4*a2
      local r2 = r*r
      local r3 = r2*r
      local r4 = r2*r2
      local r5 = r4*r
      local r6 = r3*r3

      local rho_sq = r2 + a2*ct2
      local rho = math.sqrt( r2 + a2*ct2 )
      local sqrt_fac = math.sqrt(1.0 + 2.0*M*r/rho_sq)
      local denom_root = math.sqrt(2.0*M*r + rho_sq)

      -- D^i contravariant components: (Dr, Dtheta, Dphi)
      local Dr =
        -(B_0*M*a*sqrt_fac*(st2 - 2.0)*(a2*st2 - a2 + r2)*
          (a4 - a4*st2 + r4 + 2.0*a2*r2 - a2*r2*st2 + 2.0*M*a2*r*st2))/
        ((rho_sq * rho_sq * rho_sq)*(2.0*M*r + rho_sq))
        -(2.0*B_0*M*a*r*st2*sqrt_fac*
          (r5 + 2.0*a2*r3*ct2 + 2.0*M*a2*r2 - M*a4*ct2 - M*a4*ct2*ct2 + a4*r*ct2*ct2 - M*a2*r2*st2))/
        ((rho_sq * rho_sq * rho_sq)*(2.0*M*r + rho_sq))

      local Dtheta =
        -(2.0*B_0*M*a*a*a*r*ct*st*sqrt_fac*(ct2 + 1.0))/
        ((rho_sq * rho_sq)*(2.0*M*r + rho_sq))

      local Dphi =
        -(B_0*M*a2*sqrt_fac*(st2 - 2.0)*(a2*st2 - a2 + r2))/(rho_sq * rho_sq * rho_sq)
        -(2.0*B_0*M*r*sqrt_fac*
          (r5 + 2.0*a2*r3*ct2 + 2.0*M*a2*r2 - M*a4*ct2 - M*a4*ct2*ct2 + a4*r*ct2*ct2 - M*a2*r2*st2))/
        ((rho_sq * rho_sq)*(a2 - a2*st2 + r2)*(2.0*M*r + rho_sq))

      -- B^i contravariant components: (Br, Btheta, Bphi)
      local Br =
        -(B_0*ct*
          (a6*ct2*ct2 + r6 + a2*r4 + 2.0*a2*r4*ct2 + 2.0*a4*r2*ct2 + a4*r2*ct2*ct2
          - 2.0*M*a4*r - 2.0*M*a4*r*ct2*ct2 - 4.0*M*a2*r3*ct2))/
        ((rho_sq * rho_sq * rho)*denom_root)

      local Btheta =
        (B_0*st*
          (r5 + 2.0*a2*r3*ct2 + 2.0*M*a2*r2 - M*a4*ct2 - M*a4*ct2*ct2 + a4*r*ct2*ct2 - M*a2*r2*st2))/
        ((rho_sq * rho_sq * rho)*denom_root)

      local Bphi =
        -(B_0*a*ct*(a4*ct2*ct2 + 2.0*M*r3 + r4 + 2.0*a2*r2*ct2 - 2.0*M*a2*r))/
        ((rho_sq * rho_sq * rho)*denom_root)

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