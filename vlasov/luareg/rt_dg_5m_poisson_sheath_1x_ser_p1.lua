-- Sheath formation with two DG Euler (5-moment, gas_gamma = 5/3) fluids, electrons and ions, and the
-- Poisson field: the fluid counterpart of rt_vp_sheath_1x1v.c.
-- Half domain: symmetry plane at x = 0 (reflecting fluids, zero normal field), absorbing grounded
-- wall at x = Lx (phi = 0). The absorbing wall of a fluid is a boundary flux: that of the particles
-- of the Maxwellian of the wall cell that move toward the wall (none come back). A boundary-flux
-- source replaces the particles lost to the wall, and a temperature source (the moments of the
-- BGK-form source of the kinetic tests) at the source temperature acts upstream (x < L_nu) to hold
-- the pre-sheath plasma at that temperature. The position map clusters the cells at the wall.
-- Run to t = 1000/omega_pe.
-- Figures of merit at t = 1000 (16 cells): midplane-to-wall potential drop e*dphi/Te = 2.129, ion
-- wall flux / (n_mid c_s) = 0.493, ion speed at the wall / c_s = 1.749, with c_s = sqrt((Te + Ti)/mi)
-- (rt_vp_sheath_1x1v.c: 1.925, 0.509, 1.629). Both species lose particles at the same rate, so the
-- total charge stays constant.

local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 100.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
Te_over_Ti = 1.0 -- Electron temperature / ion temperature.
gas_gamma = 5.0 / 3.0 -- Adiabatic index of the fluids.

-- Derived physical quantities (using normalized code units).
vte = 1.0 -- Electron thermal velocity.
vti = vte / math.sqrt(Te_over_Ti * mass_ion) -- Ion thermal velocity.

lambda_D = math.sqrt(epsilon0 * vte * vte * mass_elc / (n0 * charge_ion * charge_ion)) -- Electron Debye length.
omega_pe = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)) -- Electron plasma frequency.
nu_Q_elc = vte / (5.0 * lambda_D) -- Rate of the electron temperature source (upstream).
nu_Q_ion = vti / (5.0 * lambda_D) -- Rate of the ion temperature source (upstream).

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Lx = 64.0 * lambda_D -- Domain size (configuration space: x-direction).
Ls = 40.0 * lambda_D -- Extent of the source region.
L_nu = 32.0 * lambda_D -- Extent of the temperature-source region.
map_A = 0.8 -- Amplitude of the position map clustering the cells at the wall.
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 1000.0 / omega_pe -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate integrated L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Source profile: a density 2 (Ls - x)/Ls inside the source region (it integrates to Ls, so the
-- boundary-flux source replaces the wall losses exactly), at the plasma temperature and without
-- drift. The temperature source relaxes to the same temperature.
source_density = function (t, xn)
  local x = xn[1]
  local n = 0.0

  if x < Ls then
    n = 2.0 * (Ls - x) / Ls -- Source number density (inside the source region).
  end

  return n
end

-- Rates of the temperature source: nu_Q upstream (x < L_nu), vanishing over 6 lambda_D towards the
-- collisionless region next to the wall.
temp_source_rate = function (nu_Q)
  return function (t, xn)
    local x = xn[1]

    return nu_Q / (1.0 + math.exp((x - L_nu) / (6.0 * lambda_D)))
  end
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
  -- (1 + A) Lx/Nx at the symmetry boundary to (1 - A) Lx/Nx at the wall.
  mapc2pPos = {
    {
      pmap = function (t, xn)
        local xc = xn[1]

        return (1.0 + map_A) * xc - map_A * (xc * xc) / Lx
      end
    },
  },

  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coordinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Electrons: an Euler fluid with the HLLC flux.
  elc = Vlasov.FluidSpecies.new {
    charge = charge_elc, mass = mass_elc,
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions.
    init = function (t, xn)
      local rho = mass_elc * n0 -- Electron mass density.
      local E = rho * vte * vte / (gas_gamma - 1.0) -- Electron total energy density (no flow).

      return rho, 0.0, 0.0, 0.0, E
    end,

    source = {
      sourceID = G0.Source.BoundaryFlux,
      sourceLength = Ls,
      sourceSpecies = "ion",
      writeSource = true,

      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.LTE,

          densityInit = source_density,

          temperatureInit = function (t, xn)
            return mass_elc * vte * vte -- Electron source temperature.
          end,

          driftVelocityInit = function (t, xn)
            return 0.0, 0.0, 0.0 -- Source drift velocity.
          end
        }
      }
    },

    -- Temperature source at the source temperature.
    sourceBGK = {
      rateProfile = temp_source_rate(nu_Q_elc),

      tempShape = function (t, xn)
        return mass_elc * vte * vte -- Electron source temperature.
      end
    },

    bcx = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcAbsorb } -- Symmetry boundary and absorbing wall (x-direction).
  },

  -- Ions: an Euler fluid with the HLLC flux.
  ion = Vlasov.FluidSpecies.new {
    charge = charge_ion, mass = mass_ion,
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions.
    init = function (t, xn)
      local rho = mass_ion * n0 -- Ion mass density.
      local E = rho * vti * vti / (gas_gamma - 1.0) -- Ion total energy density (no flow).

      return rho, 0.0, 0.0, 0.0, E
    end,

    source = {
      sourceID = G0.Source.BoundaryFlux,
      sourceLength = Ls,
      sourceSpecies = "ion",
      writeSource = true,

      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.LTE,

          densityInit = source_density,

          temperatureInit = function (t, xn)
            return mass_ion * vti * vti -- Ion source temperature.
          end,

          driftVelocityInit = function (t, xn)
            return 0.0, 0.0, 0.0 -- Source drift velocity.
          end
        }
      }
    },

    -- Temperature source at the source temperature.
    sourceBGK = {
      rateProfile = temp_source_rate(nu_Q_ion),

      tempShape = function (t, xn)
        return mass_ion * vti * vti -- Ion source temperature.
      end
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
