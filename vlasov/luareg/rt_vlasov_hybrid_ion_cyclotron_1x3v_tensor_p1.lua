-- Parallel ion cyclotron wave with kinetic ions and fluid electrons (hybrid Vlasov-Maxwell, 1x3v).
-- A left-hand circularly polarized wave propagates along the background magnetic field (x) at k d_i
-- = 1 in a plasma with beta_i = beta_e = 0.1, vA / c = 0.1 and mi / me = 100, and is damped by the
-- cyclotron resonance of the Maxwellian ions. The fluid-EM coupling is the implicit staggered
-- split: the source step takes the ion current as known and solves the electron momentum with E,
-- the transport step moves the ions with that E. The ions have LBO self-collisions with a constant
-- frequency nu_ii = 1e-3 Omega_ci. A kinetic species cannot collide with a fluid species, so there
-- are no ion-electron collisions; the mixed counterpart
-- vlasov/creg/rt_pkpm_vlasov_ion_cyclotron_1x3v.c (kinetic ions and PKPM electrons, at mi / me =
-- 25) has the same ion self-collisions and, in addition, electron self-collisions and the
-- collisions of the ion-electron pair.
--
-- Collisionless linear theory (Maxwell's equations with the displacement current, kinetic ions, and
-- cold electrons, since the pressure does not enter a transverse wave along B): omega^2 - c^2 k^2 +
-- omega_pi^2 (omega / (sqrt(2) k v_ti)) Z((omega - Omega_ci) / (sqrt(2) k v_ti)) - omega_pe^2 omega
-- / (omega + |Omega_ce|) = 0 gives omega = 0.51572 Omega_ci and a damping rate of 0.05099 Omega_ci.
-- Kinetic electrons give 0.51575 and 0.05101.
--
-- The wave starts with a magnetic field of 1e-4 B0, Maxwellian ions and an electron fluid with
-- transverse velocities u = -1.5953 vA dB / B0 and -0.5992 vA dB / B0, and E = 0.6031 vA (dB / B0)
-- x B0. These amplitudes solve the linear initial value problem for which the three other normal
-- modes at this wave number and polarization are not excited: the whistler that propagates the
-- opposite way (omega = -1.566 Omega_ci) and two high-frequency electron modes. The amplitude
-- relations of the wave itself start the whistler, which is not damped, at 4.5% of the field. The
-- Maxwellian is not the kinetic eigenfunction: the rest of the perturbation phase-mixes during the
-- first period, which leaves the wave 0.978 of the initial field, so the figures of merit are taken
-- after it.
--
-- Figures of merit (tensor p1 with the hybrid phase basis, 16 x 16 x 6 x 6 cells, 2327 steps,
-- omega_pi dt = 0.131, omega_pe dt = |Omega_ce| dt = 1.31, set by the light wave): between t = 1
-- and 2.5 wave periods the k d_i = 1 mode of B_y + i B_z decays at 0.05150 Omega_ci (+1.0%) with
-- omega = 0.51714 Omega_ci (+0.3%); the history of the transverse magnetic energy over the same
-- interval gives 0.05153 Omega_ci (+1.1%). The total energy changes by -4e-13 (5e-5 of the initial
-- magnetic energy of the wave).

local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 100.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
B0 = 1.0 -- Background magnetic field (Alfven speed 0.1 c).
beta_ion = 0.1 -- Ion plasma beta.
Ti = 0.5 * beta_ion * B0 * B0 / (mu0 * n0) -- Ion temperature.
Te = Ti -- Electron temperature.
gas_gamma = 5.0 / 3.0 -- Adiabatic index of the electron fluid.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.
vA = B0 / math.sqrt(mu0 * n0 * mass_ion) -- Alfven speed.
omega_ci = charge_ion * B0 / mass_ion -- Ion cyclotron frequency.
nu_ii = 1.0e-3 * omega_ci -- Ion-ion collision frequency.
omega_pi = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_ion)) -- Ion plasma frequency.
di = 1.0 / (math.sqrt(epsilon0 * mu0) * omega_pi) -- Ion inertial length.

alpha = 1.0e-4 -- Amplitude of the magnetic field of the wave, relative to B0.
k0 = 1.0 / di -- Wave number (k d_i = 1).
-- Root of the dispersion relation of kinetic ions and cold fluid electrons (see the header).
omega_r = 0.515720033889 * omega_ci -- Real frequency of the wave.
-- Amplitudes of the initial condition that starts the ion cyclotron wave alone (see the header).
ion_vel_fac = 1.595254032963 -- Ion velocity: u = -ion_vel_fac vA dB / B0.
elc_vel_fac = 0.599201362932 -- Electron velocity: u = -elc_vel_fac vA dB / B0.
E_fac = 0.603120213183 -- Electric field: E = E_fac vA (dB / B0) x B0.

-- Magnetic field of the wave relative to B0, circularly polarized in the sense the ions gyrate.
local function wave_dB(x)
  return -alpha * math.sin(k0 * x), alpha * math.cos(k0 * x)
end

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: vx-direction).
Nvy = 6 -- Cell count (velocity space: vy-direction).
Nvz = 6 -- Cell count (velocity space: vz-direction).
Lx = 2.0 * pi / k0 -- Domain size (configuration space: x-direction).
vx_max = 5.0 * vti -- Domain boundary (velocity space: vx-direction).
vy_max = 4.5 * vti -- Domain boundary (velocity space: vy-direction).
vz_max = 4.5 * vti -- Domain boundary (velocity space: vz-direction).
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set (hybrid phase-space basis of the kinetic species).
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 2.5 * (2.0 * pi / omega_r) -- Final simulation time (2.5 wave periods).
num_frames = 5 -- Number of output frames (two per wave period).
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
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Ions (kinetic).
  ion = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,

    -- Velocity space grid.
    lower = { -vx_max, -vy_max, -vz_max },
    upper = { vx_max, vy_max, vz_max },
    cells = { Nvx, Nvy, Nvz },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local x, vx, vy, vz = xn[1], xn[2], xn[3], xn[4]

          -- Maxwellian with the transverse velocity of the initial condition.
          local dBy, dBz = wave_dB(x)
          local uy = -ion_vel_fac * vA * dBy
          local uz = -ion_vel_fac * vA * dBz
          local v_sq = vx * vx + (vy - uy) * (vy - uy) + (vz - uz) * (vz - uz)

          return n0 / math.pow(2.0 * pi * vti * vti, 1.5) * math.exp(-v_sq / (2.0 * vti * vti))
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ii
      end
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Electrons (fluid).
  elc = Vlasov.FluidSpecies.new {
    charge = charge_elc, mass = mass_elc,
    equation = Euler.new { gasGamma = gas_gamma },

    init = function (t, xn)
      local x = xn[1]

      -- Transverse velocity of the initial condition (the density and the pressure are uniform).
      local dBy, dBz = wave_dB(x)
      local uy = -elc_vel_fac * vA * dBy
      local uz = -elc_vel_fac * vA * dBz
      local rho = mass_elc * n0

      local energy = n0 * Te / (gas_gamma - 1.0) + 0.5 * rho * (uy * uy + uz * uz)

      return rho, 0.0, rho * uy, rho * uz, energy
    end,
  },

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      -- Electric field E = E_fac vA (dB / B0) x B0, with B0 along x.
      local dBy, dBz = wave_dB(x)
      local Ey = E_fac * vA * B0 * dBz
      local Ez = -E_fac * vA * B0 * dBy

      return 0.0, Ey, Ez, B0, B0 * dBy, B0 * dBz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
  }
}

vlasovApp:run()
