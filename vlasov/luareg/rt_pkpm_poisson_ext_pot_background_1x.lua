-- A held, nonuniform electron background driven by a static external scalar potential phi_ext =
-- (E0 / k) cos(k x), with the PKPM system of equations and the Poisson field (1x, normalized
-- parallel velocity w = v_par/xi). Electrons (vt = 0.01, lambda_D = 0.01, no ions) with the
-- density n0 (1 + 0.5 cos(k x)), k = 1, are held about that state at rest, and start on it. The
-- hold is derived with no electric force of the fields, so the external field, E_ext = E0 sin(k x)
-- with E0 = 1e-4, acts on the electrons in full, for every type of species (derived with it, the
-- hold of a kinetic species would cancel it for good, and that of a fluid would not). For a cold
-- plasma every fluid element oscillates at its local plasma frequency about the displacement that
-- shields the external field: u_x = (q / m) E_ext sin(omega_p(x) t) / omega_p(x), and the field of
-- the electrons is E_x = -E_ext (1 - cos(omega_p(x) t)), with omega_p(x)^2 = omega_pe^2 (1 + 0.5
-- cos(k x)), so the total field is E_ext cos(omega_p(x) t). Thermal corrections are of the order
-- of (k lambda_D)^2. The potential of the electrons is solved from their departure from the
-- background; the coupling is the kick-drift-kick split, whose kicks take E = -grad(phi + phi_ext)
-- (cfl_frac = 0.2). A uniform external magnetic field along x gives PKPM its field direction; LBO
-- self-collisions at nu = 0.01 omega_pe act on the electrons. The same problem with the other
-- models and fields: rt_dg_5m_poisson_ext_pot_background_1x_ser_p1.lua,
-- rt_vp_ext_pot_background_1x1v_ser_p2.lua, rt_dg_5m_ext_field_background_1x_ser_p1.lua,
-- rt_vlasov_ext_field_background_1x1v_ser_p2.lua.
--
-- Figures of merit at t = 4 pi / omega_pe (64 x 16 cells, 248 steps, omega_pe dt = 0.051): the
-- cell averages of the field of the electrons, -grad(phi), differ from the cold local solution by
-- 5.9e-3 in L2 and those of u_x by 1.0e-2. The density changes by at most 3.8e-4 of itself, the
-- compression of the oscillation.


local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi -- Pi.

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

vt = 0.01 -- Thermal velocity.
B0 = 1.0 -- Magnetic field along x (sets the field direction b).
n0 = 1.0 -- Number density.
a_bg = 0.5 -- Relative amplitude of the density variation of the background.
E0 = 1.0e-4 -- Amplitude of the external electric field.
k = 1.0 -- Wavenumber.

-- Derived physical quantities (using normalized code units).
omega_pe = math.sqrt(n0 * charge_elc * charge_elc / (epsilon0 * mass_elc)) -- Plasma frequency.
nu_ee = 0.01 * omega_pe -- Electron-electron collision frequency.

-- Simulation parameters.
Nx = 64 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 2.0 * pi / k -- Domain size (configuration space: x-direction).
vx_max = 6.0 -- Domain boundary (velocity space: w = v_par/xi).
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

-- Normalized electron distributions [F, G, M0_perp (x, y, z), N_perp (x, y, z)]: a Maxwellian in
-- w = v_par/xi and no first harmonic. The initial condition and the held background.
function elcDist (t, xn)
  local w = xn[2]

  local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * pi) -- Normalized Maxwellian in w = v_par/xi.

  return F, F, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
end

-- Fluid background of the electrons: a density that varies along x, at rest, with a uniform
-- temperature.
function elcFluidBackground (t, xn)
  local x = xn[1]

  local n = n0 * (1.0 + a_bg * math.cos(k * x)) -- Electron number density.

  return mass_elc * n, 0.0, 0.0, 0.0, 1.5 * n * mass_elc * vt * vt
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
  elc = Vlasov.PKPMSpecies.new {
    charge = charge_elc, mass = mass_elc,

    kinetic = {
      -- Velocity space grid.
      lower = { -vx_max },
      upper = { vx_max },
      cells = { Nvx },

      -- Initial conditions.
      numInit = 1,
      projections = {
        {
          projectionID = G0.Projection.Func,
          init = elcDist
        }
      },
      -- The electrons are held about their background: the normalized Maxwellian here, the fluid
      -- background in the fluid block.
      background = {
        projectionID = G0.Projection.Func,
        init = elcDist
      },

      collisions = {
        collisionID = G0.Collisions.LBO,

        selfNu = function (t, xn)
          return nu_ee -- Electron-electron collision frequency.
        end
      }
    },

    fluid = {
      init = elcFluidBackground,
      background = elcFluidBackground,
      evalOnNodes = true
    },

    evolve = true -- Evolve species?
  },

  -- Field: the electrostatic potential of the departure of the electrons from their background, with the external potential,
  -- and the static external magnetic field that sets b.
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
    end,

    externalFieldInit = function (t, xn)
      return 0.0, 0.0, 0.0, B0, 0.0, 0.0
    end
  }
}

vlasovApp:run()
