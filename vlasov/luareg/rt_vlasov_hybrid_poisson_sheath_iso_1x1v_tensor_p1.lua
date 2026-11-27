-- Sheath formation with kinetic ions and an isothermal electron fluid (hybrid Vlasov-Poisson, 1x1v),
-- rt_vlasov_hybrid_poisson_sheath_1x1v.c with the isothermal closure of the electrons (thermal speed
-- vte, so no electron temperature source).
-- Half domain: symmetry plane at x = 0 (reflecting ions and electrons, zero normal field), absorbing
-- grounded wall at x = Lx (phi = 0). The absorbing wall of the electron fluid is a boundary flux:
-- that of the electrons of the Maxwellian of the wall cell (at the thermal speed vte) that move
-- toward the wall (none come back). A boundary-flux source of each species replaces the ions lost
-- to the wall, and the BGK-form temperature source of the ions at the source temperature acts
-- upstream (x < L_nu). The position map clusters the cells at the wall. Run to t = 1000/omega_pe.
-- Figures of merit at t = 1000 (16 x 24 cells, tensor p1 with the hybrid phase basis):
-- midplane-to-wall potential drop e*dphi/Te = 2.258, ion wall flux / (n_mid c_s) = 0.494, ion speed
-- at the wall / c_s = 1.754, with c_s = sqrt((Te + Ti)/mi) (adiabatic electrons,
-- rt_vlasov_hybrid_poisson_sheath_1x1v.c: 2.128, 0.495, 1.712; kinetic electrons, rt_vp_sheath_1x1v.c:
-- 1.925, 0.509, 1.629).

local Vlasov = G0.Vlasov
local IsoEuler = G0.Vlasov.Eq.IsoEuler

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 100.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
Vx_drift_ion = 0.0 -- Ion drift velocity (x-direction).

-- Derived physical quantities (using normalized code units).
Te = 1.0 * charge_ion -- Electron temperature.
Ti = 1.0 * charge_ion -- Ion temperature.

vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.

lambda_D = math.sqrt(epsilon0 * Te / (n0 * charge_ion * charge_ion)) -- Electron Debye length.
omega_pe = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)) -- Electron plasma frequency.

nu_Q_ion = vti / (5.0 * lambda_D) -- Rate of the ion temperature source (upstream).

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Nvx = 24 -- Cell count (velocity space: vx-direction).
Lx = 64.0 * lambda_D -- Domain size (configuration space: x-direction).
Ls = 40.0 * lambda_D -- Domain size (source).
L_nu = 32.0 * lambda_D -- Extent of the collisional (fixed-temperature) region.
vx_max_ion = 8.0 * vti -- Domain boundary (ion velocity space: vx-direction).
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set (hybrid phase-space basis of the kinetic species).
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 1000.0 / omega_pe -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Source profile: a density 2 (Ls - x)/Ls inside the source region (it integrates to Ls, so the
-- boundary-flux source replaces the wall losses exactly), at the plasma temperature.
source_density = function (t, xn)
  local x = xn[1]

  local n = 0.0

  if x < Ls then
    n = 2.0 * (Ls - x) / Ls -- Source total number density (inside the source region).
  else
    n = 0.0 -- Source total number density (outside the source region).
  end

  return n
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

  -- Position map clustering the cells at the wall (x = Lx): the cell size goes from
  -- (1 + A) * Lx / Nx at the symmetry boundary to (1 - A) * Lx / Nx at the wall.
  mapc2pPos = {
    {
      pmap = function (t, xn)
        local xc = xn[1]
        local A = 0.8
        return (1.0 + A) * xc - A * (xc * xc) / Lx
      end
    },
  },
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
    lower = { -vx_max_ion },
    upper = { vx_max_ion },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return n0 -- Ion total number density.
        end,
        temperatureInit = function (t, xn)
          return Ti -- Ion isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift_ion -- Ion drift velocity.
        end
      }
    },

    source = {
      sourceID = G0.Source.BoundaryFlux,
      sourceLength = Ls,
      sourceSpecies = "ion",

      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.LTE,

          densityInit = source_density,
          temperatureInit = function (t, xn)
            return Ti -- Ion source isotropic temperature.
          end,
          driftVelocityInit = function (t, xn)
            return Vx_drift_ion -- Ion source drift velocity.
          end
        }
      }
    },

    sourceBGK = {
      rateProfile = function (t, xn)
        local x = xn[1]

        local nu = nu_Q_ion / (1.0 + math.exp((x - L_nu) / (6.0 * lambda_D))) -- Ion temperature-source rate.

        return nu
      end,
      tempShape = function (t, xn)
        return Ti -- Ion source temperature.
      end
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcReflect
      },
      upper = {
        type = G0.SpeciesBc.bcAbsorb
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Electrons: an isothermal Euler fluid at the thermal speed vte (Roe flux).
  elc = Vlasov.FluidSpecies.new {
    charge = charge_elc, mass = mass_elc,
    equation = IsoEuler.new { vThermal = vte },

    -- Initial conditions.
    init = function (t, xn)
      local rho = mass_elc * n0 -- Electron mass density.

      return rho, 0.0, 0.0, 0.0
    end,

    source = {
      sourceID = G0.Source.BoundaryFlux,
      sourceLength = Ls,
      sourceSpecies = "ion",

      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.LTE,

          densityInit = source_density,
          temperatureInit = function (t, xn)
            return Te -- Electron source isotropic temperature.
          end
        }
      }
    },

    bcx = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcAbsorb } -- Symmetry boundary and absorbing wall (x-direction).
  },

  -- Field: the electrostatic potential, symmetric at x = 0 and grounded at the wall.
  field = Vlasov.Field.new {
    fieldID = G0.FieldModel.Phi,
    epsilon0 = epsilon0,

    poissonBcs = {
      lowerType = {
        G0.PoissonBc.bcNeumann
      },
      upperType = {
        G0.PoissonBc.bcDirichlet
      },
      lowerValue = {
        0.0
      },
      upperValue = {
        0.0
      }
    }
  }
}

vlasovApp:run()
