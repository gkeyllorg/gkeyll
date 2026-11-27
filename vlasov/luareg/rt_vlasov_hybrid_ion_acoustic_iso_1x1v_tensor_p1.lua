-- Ion acoustic wave with kinetic ions and isothermal fluid electrons (hybrid Vlasov-Maxwell, 1x1v),
-- the variant of rt_vlasov_hybrid_ion_acoustic_1x1v with an isothermal Euler fluid (p = n Te, Roe
-- flux) for the electrons. Maxwellian ions (mi/me = 100, Te/Ti = 2) and the electron fluid carry a
-- travelling wave at k lambda_De = 0.05, started from the linear densities and velocities of the
-- wave with an ion density amplitude of 1e-4. The ions start as a Maxwellian with these moments,
-- which is not the kinetic eigenfunction: the rest of the perturbation phase-mixes during the first
-- period, so the figures of merit are taken after it. The ions have LBO self-collisions at nu_ii =
-- 2e-4 omega_pi, as in rt_vlasov_hybrid_ion_acoustic_1x1v; there are no collisions with the
-- electron fluid. The fluid-EM coupling is the implicit staggered split: the source step takes the
-- ion current as known and solves the electron momentum with E, the transport step moves the ions
-- with that E. Collisionless linear theory, 1 + chi_i + chi_e = 0 with chi_i = (1 + zeta Z(zeta)) /
-- (k lambda_Di)^2, zeta = omega / (sqrt(2) k v_ti), and chi_e = -omega_pe^2 / (omega^2 - k^2 Te /
-- me): omega = 0.08339 omega_pi, damping rate 0.02024 omega_pi. Kinetic electrons give omega =
-- 0.08590 omega_pi and a damping rate of 0.02508 omega_pi; adiabatic fluid electrons (gamma = 5/3,
-- rt_vlasov_hybrid_ion_acoustic_1x1v) 0.09470 and 0.01333. With the ion collisions (linearized
-- kinetic equation in Hermite moments, of which the LBO damps moment m >= 3 at m nu_ii, and the
-- isothermal electron fluid) the theory gives omega = 0.08347 omega_pi and a damping rate of
-- 0.01994 omega_pi. Figures of merit (tensor p1 with the hybrid phase basis, 32 x 48 cells, 1727
-- steps, omega_pe dt = 1.31): between t = 1 and 3 wave periods the ion density wave decays at
-- 0.01986 omega_pi (-0.4% against the theory with collisions) with omega = 0.08306 omega_pi
-- (-0.5%). Without collisions the run gives 0.02009 and 0.08295: the collisions change the rate by
-- -1.2% (theory -1.5%) and omega by +0.13% (theory +0.10%). The total energy (field, ions, and
-- (1/2) rho u^2 + v_t^2 rho ln rho of the electrons) changes by +3e-11 (6.1e-4 of twice the initial
-- kinetic energy of the ion flow).
local Vlasov = G0.Vlasov
local IsoEuler = G0.Vlasov.Eq.IsoEuler

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
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.
omega_pi = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_ion)) -- Ion plasma frequency.
nu_ii = 2.0e-4 * omega_pi -- Ion-ion collision frequency.

alpha = 1.0e-4 -- Amplitude of the ion density of the wave, relative to n0.
k0 = 0.05 -- Wave number (k lambda_De = 0.05).
-- Root of the dispersion relation of kinetic ions and isothermal fluid electrons (see the header),
-- in units of the electron plasma frequency.
omega_r = 0.008338515405 -- Real frequency of the wave.
omega_i = -0.002024346176 -- Growth rate of the wave.

-- Densities of the wave per unit electric field (real and imaginary parts), all perturbations
-- proportional to exp(i (k x - omega t)). Electrons (continuity and momentum of an isothermal fluid):
-- n_e = -i q_e n0 k E / (Te k^2 - m_e omega^2). Ions (Gauss's law):
-- i k eps0 E = q_i n_i + q_e n_e.
local den_r = Te * k0 * k0 - mass_elc * (omega_r * omega_r - omega_i * omega_i)
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

  -- Electrons (isothermal fluid).
  elc = Vlasov.FluidSpecies.new {
    charge = charge_elc, mass = mass_elc,
    equation = IsoEuler.new { vThermal = math.sqrt(Te / mass_elc) },

    init = function (t, xn)
      local x = xn[1]

      -- Density and velocity (continuity: u = omega n / (k n0)) of the wave.
      local cos_kx, sin_kx = math.cos(k0 * x), math.sin(k0 * x)
      local n1 = E0 * (elc_resp_r * cos_kx - elc_resp_i * sin_kx)
      local ux = E0 / (k0 * n0) * ((omega_r * elc_resp_r - omega_i * elc_resp_i) * cos_kx -
        (omega_r * elc_resp_i + omega_i * elc_resp_r) * sin_kx)
      local rho = mass_elc * (n0 + n1)

      return rho, rho * ux, 0.0, 0.0
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
