-- Ion acoustic wave with kinetic ions and fluid electrons (hybrid Vlasov-Maxwell, 1x1v).
-- Maxwellian ions (mi/me = 100, Te/Ti = 2) and an adiabatic electron fluid (gamma = 5/3) carry a
-- travelling wave at k lambda_De = 0.05, started from the linear densities and velocities of the
-- wave with an ion density amplitude of 1e-4. The ions start as a Maxwellian with these moments,
-- which is not the kinetic eigenfunction: the rest of the perturbation phase-mixes during the first
-- period, so the figures of merit are taken after it. The ions have LBO self-collisions at nu_ii =
-- 2e-4 omega_pi, the frequency of the PKPM test of this problem
-- (vlasov/luareg/rt_pkpm_hybrid_ion_acoustic_1x.lua). Collisions between a kinetic species and a
-- fluid species are not supported, so there is no exchange with the electron fluid, which the PKPM
-- ions have (nu_ei = 2e-3 omega_pi). In one velocity dimension the collisions conserve the one
-- temperature of the ions; a PKPM species has three degrees of freedom and its collisions also
-- exchange parallel and perpendicular energy, so the figures of this test are not the values the
-- PKPM test should reproduce. The fluid-EM coupling is the implicit staggered split: the source
-- step takes the ion current as known and solves the electron momentum with E, the transport step
-- moves the ions with that E. Collisionless linear theory, 1 + chi_i + chi_e = 0 with chi_i = (1 +
-- zeta Z(zeta)) / (k lambda_Di)^2, zeta = omega / (sqrt(2) k v_ti), and chi_e = -omega_pe^2 /
-- (omega^2 - gamma k^2 Te / me): omega = 0.09470 omega_pi, damping rate 0.01333 omega_pi. Kinetic
-- electrons, which are isothermal at this frequency, give omega = 0.08590 omega_pi and a damping
-- rate of 0.02508 omega_pi. With the ion collisions (linearized kinetic equation in Hermite
-- moments, of which the LBO damps moment m >= 3 at m nu_ii, and the electron fluid) the theory
-- gives omega = 0.09466 omega_pi and a damping rate of 0.01314 omega_pi. Figures of merit (tensor
-- p1 with the hybrid phase basis, 32 x 48 cells, omega_pe dt = 1.01): between t = 1 and 3 wave
-- periods the ion density wave decays at 0.01307 omega_pi (-0.5% against the theory with
-- collisions) with omega = 0.09448 omega_pi (-0.2%). The total energy changes by -5e-13 (4e-5 of
-- twice the initial kinetic energy of the ion flow).

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
Te = 1.0 -- Electron temperature.
Ti = 0.5 -- Ion temperature.
gas_gamma = 5.0 / 3.0 -- Adiabatic index of the electron fluid.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.
omega_pi = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_ion)) -- Ion plasma frequency.
nu_ii = 2.0e-4 * omega_pi -- Ion-ion collision frequency.

alpha = 1.0e-4 -- Amplitude of the ion density of the wave, relative to n0.
k0 = 0.05 -- Wave number (k lambda_De = 0.05).
-- Root of the dispersion relation of kinetic ions and adiabatic fluid electrons (see the header),
-- in units of the electron plasma frequency.
omega_r = 0.009469560873 -- Real frequency of the wave.
omega_i = -0.001332945847 -- Growth rate of the wave.

-- Densities of the wave per unit electric field (real and imaginary parts), all perturbations
-- proportional to exp(i (k x - omega t)). Electrons (continuity and momentum of an adiabatic fluid):
-- n_e = -i q_e n0 k E / (gamma Te k^2 - m_e omega^2). Ions (Gauss's law):
-- i k eps0 E = q_i n_i + q_e n_e.
local den_r = gas_gamma * Te * k0 * k0 - mass_elc * (omega_r * omega_r - omega_i * omega_i)
local den_i = -2.0 * mass_elc * omega_r * omega_i
local den_sq = den_r * den_r + den_i * den_i
elc_resp_r = -charge_elc * n0 * k0 * den_i / den_sq
elc_resp_i = -charge_elc * n0 * k0 * den_r / den_sq
ion_resp_r = -charge_elc * elc_resp_r / charge_ion
ion_resp_i = (k0 * epsilon0 - charge_elc * elc_resp_i) / charge_ion
E0 = alpha * n0 / math.sqrt(ion_resp_r * ion_resp_r + ion_resp_i * ion_resp_i) -- Field amplitude.

-- Simulation parameters.
Nx = 32 -- Cell count (configuration space: x-direction).
Nvx = 48 -- Cell count (velocity space: vx-direction).
Lx = 2.0 * pi / k0 -- Domain size (configuration space: x-direction).
vx_max = 6.0 * vti -- Domain boundary (velocity space: vx-direction).
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set (hybrid phase-space basis of the kinetic species).
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 3.0 * (2.0 * pi / omega_r) -- Final simulation time (three wave periods).
num_frames = 3 -- Number of output frames (one per wave period).
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
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local x, vx = xn[1], xn[2]

          -- Maxwellian with the density of the wave and the velocity that moves it at the phase
          -- speed (continuity: u = omega_r n / (k n0)).
          local n1 = E0 * (ion_resp_r * math.cos(k0 * x) - ion_resp_i * math.sin(k0 * x))
          local n = n0 + n1
          local ux = omega_r * n1 / (k0 * n0)

          return n / math.sqrt(2.0 * pi * vti * vti) * math.exp(-(vx - ux) * (vx - ux) / (2.0 * vti * vti))
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ii -- Self-collision frequency.
      end
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Electrons (adiabatic fluid).
  elc = Vlasov.FluidSpecies.new {
    charge = charge_elc, mass = mass_elc,
    equation = Euler.new { gasGamma = gas_gamma },

    init = function (t, xn)
      local x = xn[1]

      -- Density, velocity (continuity: u = omega n / (k n0)) and pressure (adiabatic) of the wave.
      local cos_kx, sin_kx = math.cos(k0 * x), math.sin(k0 * x)
      local n1 = E0 * (elc_resp_r * cos_kx - elc_resp_i * sin_kx)
      local ux = E0 / (k0 * n0) * ((omega_r * elc_resp_r - omega_i * elc_resp_i) * cos_kx -
        (omega_r * elc_resp_i + omega_i * elc_resp_r) * sin_kx)
      local rho = mass_elc * (n0 + n1)
      local p = n0 * Te + gas_gamma * Te * n1

      return rho, rho * ux, 0.0, 0.0, p / (gas_gamma - 1.0) + 0.5 * rho * ux * ux
    end,
  },

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      return E0 * math.cos(k0 * x), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
  }
}

vlasovApp:run()
