-- A held, nonuniform electron background driven by a static external scalar potential phi_ext =
-- (E0 / k) cos(k x), with the DG 5-moment (Euler) equations and the Poisson field (1x). Electrons
-- (vt = 0.01, lambda_D = 0.01, no ions) with the density n0 (1 + 0.5 cos(k x)), k = 1, are held
-- about that state at rest, and start on it. The hold is derived with no electric force of the
-- fields, so the external field, E_ext = E0 sin(k x) with E0 = 1e-4, acts on the electrons in
-- full, for every type of species (derived with it, the hold of a kinetic species would cancel it
-- for good, and that of a fluid would not). For a cold plasma every fluid element oscillates at
-- its local plasma frequency about the displacement that shields the external field: u_x = (q / m)
-- E_ext sin(omega_p(x) t) / omega_p(x), and the field of the electrons is E_x = -E_ext (1 -
-- cos(omega_p(x) t)), with omega_p(x)^2 = omega_pe^2 (1 + 0.5 cos(k x)), so the total field is
-- E_ext cos(omega_p(x) t). Thermal corrections are of the order of (k lambda_D)^2. The potential
-- of the electrons is solved from their departure from the background; the coupling is the
-- kick-drift-kick split, whose kicks take E = -grad(phi + phi_ext) (cfl_frac = 0.2). The same
-- problem with the other models and fields: rt_vp_ext_pot_background_1x1v_ser_p2.lua,
-- rt_pkpm_poisson_ext_pot_background_1x.lua, rt_dg_5m_ext_field_background_1x_ser_p1.lua,
-- rt_vlasov_ext_field_background_1x1v_ser_p2.lua.
--
-- Figures of merit at t = 4 pi / omega_pe (64 cells, serendipity p1, 77 steps, omega_pe dt =
-- 0.163): the cell averages of the field of the electrons, -grad(phi), differ from the cold local
-- solution by 9.2e-3 in L2, and by 6.0e-3 with the leapfrog frequency (2 / dt) asin(omega_p dt /
-- 2) of the split in it; those of u_x differ by 1.7e-2. The density changes by at most 3.8e-4 of
-- itself, the compression of the oscillation.


local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi -- Pi.

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

n0 = 1.0 -- Number density.
T0 = 1.0e-4 -- Electron temperature.
a_bg = 0.5 -- Relative amplitude of the density variation of the background.
E0 = 1.0e-4 -- Amplitude of the external electric field.
k = 1.0 -- Wavenumber.

-- Derived physical quantities (using normalized code units).
omega_pe = math.sqrt(n0 * charge_elc * charge_elc / (epsilon0 * mass_elc)) -- Plasma frequency.
vt = math.sqrt(T0 / mass_elc) -- Electron thermal velocity.

-- Simulation parameters.
Nx = 64 -- Cell count (configuration space: x-direction).
Lx = 2.0 * pi / k -- Domain size (configuration space: x-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.2 -- CFL coefficient.

t_end = 4.0 * pi / omega_pe -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Background of the electrons: a density that varies along x, at rest, with a uniform temperature.
function elcBackground (t, xn)
  local x = xn[1]

  local n = n0 * (1.0 + a_bg * math.cos(k * x)) -- Number density.

  return n * mass_elc, 0.0, 0.0, 0.0, n * T0 / (gas_gamma - 1.0)
end

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
  decompCuts = { 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions.

  -- Electrons.
  elc = Vlasov.FluidSpecies.new {
    charge = charge_elc, mass = mass_elc,
    equation = Euler.new { gasGamma = gas_gamma },

    init = elcBackground,
    -- The electrons are held about their background, evaluated at the nodes (continuous).
    background = elcBackground,
    evalOnNodes = true,
  },

  -- Field: the electrostatic potential of the departure of the electrons from their background, with the external potential.
  field = Vlasov.Field.new {
    fieldID = G0.FieldModel.Phi,
    epsilon0 = epsilon0,

    poissonBcs = {
      lowerType = {
        G0.PoissonBc.bcPeriodic
      },
      upperType = {
        G0.PoissonBc.bcPeriodic
      }
    },

    -- External potentials [phi, A_x, A_y, A_z]: E_ext = -grad(phi) = E0 sin(k x).
    externalPotentialInit = function (t, xn)
      return (E0 / k) * math.cos(k * xn[1]), 0.0, 0.0, 0.0
    end
  }
}

vlasovApp:run()
