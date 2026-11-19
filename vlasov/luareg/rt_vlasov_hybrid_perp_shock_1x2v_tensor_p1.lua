-- Perpendicular shock with kinetic ions and fluid electrons (hybrid Vlasov-Maxwell, 1x2v). A
-- uniform plasma with a magnetic field along z (mi / me = 25, vA / c = 0.1, beta_i = 1, beta_e =
-- 0.5; Maxwellian ions with two velocity dimensions, an adiabatic electron fluid with gamma = 5/3)
-- flows at u0 = 0.75 vA into a reflecting, perfectly conducting wall at x = 0, and the upper
-- boundary holds the upstream state. A shock forms at the wall and travels upstream, leaving the
-- plasma at rest. The fluid-EM coupling is the implicit staggered split: the source step takes the
-- ion current as known and solves the electron momentum with E, the transport step moves the ions
-- with that E.
--
-- Reference: mass, magnetic flux, momentum and energy balance across the shock in the frame of the
-- wall, v_s (Q2 - Q1) = F2 - F1 with the Maxwell stress and the Poynting flux of the motional field
-- E_y = -u0 B_z, an ion enthalpy of 2 p_i and an electron enthalpy of gamma p_e / (gamma - 1). They
-- give the compression r = 1 + u0 / v_s of the density and of B_z, but the shock speed v_s depends
-- on how the heating divides between the species, which the conservation laws do not determine.
-- With adiabatically compressed electrons (in the run p_e / n_e^gamma rises by 0.2% across the
-- shock), v_s = 1.4210 vA, r = 1.5278 (fast Mach number 1.42, Alfven Mach number 2.21), T_i2 / T_i1
-- = 1.680 and T_e2 / T_e1 = 1.327.
--
-- Figures of merit at t = 10 / Omega_ci (tensor p1 with the hybrid phase basis, 128 x 16 x 14
-- cells, 2000 steps, omega_pi dt = 0.05, set by the light wave): the shock (B_z through the mean of
-- its upstream and downstream values) moves at 1.4045 vA between t = 4 and 10 / Omega_ci (-1.2%),
-- and at 1.316, 1.346, 1.422 and 1.446 vA between successive frames. Averages over 2 < x / d_i <
-- 9.9, behind the shock: compression 1.5294 (ions), 1.5295 (electrons) and 1.5298 (B_z), +0.1%,
-- with a wave train of 6% rms in the density and 5% in B_z; T_i2 / T_i1 = 1.631 (-2.9%), T_e2 /
-- T_e1 = 1.330 (+0.3%). The total energy (field, ions and electrons, from the frames) grows by the
-- inflow through the upper boundary to 1.4e-3, and the number of ions to 4e-4.

local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 25.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Upstream number density.
vA = 0.1 -- Upstream Alfven speed (of the ions; the speed of light is 1).
B0 = vA * math.sqrt(mu0 * n0 * mass_ion) -- Upstream magnetic field (along z).
beta_ion = 1.0 -- Upstream ion plasma beta.
beta_elc = 0.5 -- Upstream electron plasma beta.
Ti = 0.5 * beta_ion * B0 * B0 / (mu0 * n0) -- Upstream ion temperature.
Te = 0.5 * beta_elc * B0 * B0 / (mu0 * n0) -- Upstream electron temperature.
gas_gamma = 5.0 / 3.0 -- Adiabatic index of the electron fluid.
vti = math.sqrt(Ti / mass_ion) -- Upstream ion thermal velocity.
u0 = 0.75 * vA -- Speed of the upstream flow toward the wall.
omega_ci = charge_ion * B0 / mass_ion -- Upstream ion cyclotron frequency.
omega_pi = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_ion)) -- Ion plasma frequency.
di = 1.0 / (math.sqrt(epsilon0 * mu0) * omega_pi) -- Upstream ion inertial length.

-- Simulation parameters.
Nx = 128 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: vx-direction).
Nvy = 14 -- Cell count (velocity space: vy-direction).
Lx = 19.2 * di -- Domain size (configuration space: x-direction).
-- The reflecting wall mirrors the vx cells, so the vx grid is symmetric about zero.
vx_max = 8.0 * vti -- Domain boundary (velocity space: vx-direction).
vy_max = 7.0 * vti -- Domain boundary (velocity space: vy-direction).
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set (hybrid phase-space basis of the kinetic species).
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 10.0 / omega_ci -- Final simulation time.
num_frames = 5 -- Number of output frames.
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
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Ions (kinetic).
  ion = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,

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
          local vx, vy = xn[2], xn[3]

          -- The uniform upstream state, a Maxwellian flowing toward the wall.
          local v_sq = (vx + u0) * (vx + u0) + vy * vy

          return n0 / (2.0 * pi * vti * vti) * math.exp(-v_sq / (2.0 * vti * vti))
        end
      }
    },

    -- Reflecting wall, and the upstream state held at the upper boundary.
    bcx = {
      lower = {
        type = G0.SpeciesBc.bcReflect
      },
      upper = {
        type = G0.SpeciesBc.bcFixedFunc
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2ij }
  },

  -- Electrons (adiabatic fluid).
  elc = Vlasov.FluidSpecies.new {
    charge = charge_elc, mass = mass_elc,
    equation = Euler.new { gasGamma = gas_gamma },

    init = function (t, xn)
      -- The uniform upstream state.
      local rho = mass_elc * n0

      return rho, -rho * u0, 0.0, 0.0, n0 * Te / (gas_gamma - 1.0) + 0.5 * rho * u0 * u0
    end,

    bcx = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcFixedFunc },
  },

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      -- The motional electric field of the upstream flow, E = -u x B.
      return 0.0, -u0 * B0, 0.0, 0.0, 0.0, B0, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?

    -- Conducting wall, and the upstream state held at the upper boundary.
    bcx = { G0.FieldBc.bcPECWall, G0.FieldBc.bcFixedFunc },
  }
}

vlasovApp:run()
