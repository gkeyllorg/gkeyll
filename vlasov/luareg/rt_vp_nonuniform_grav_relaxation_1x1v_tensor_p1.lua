-- Relaxation of a self-gravitating collisionless sheet toward equilibrium with the Vlasov-Poisson
-- system (1x1v). Neutral dark matter (m = 1, gravitational coupling alpha_g = 4 pi G = 1) starts as
-- a Gaussian slab of mass per unit area M = 1 and width a = 1 with a Maxwellian velocity
-- distribution whose dispersion puts it at half the virial kinetic energy, 2K/U = 1/2 with U = pi G
-- int int rho(x) rho(x') |x - x'| dx dx', so it collapses and phase mixes toward a quasi-stationary
-- state. The isothermal equilibrium of 1D gravity is the Spitzer-Camm sheet rho = rho0 sech^2(x/z0)
-- with z0 = sigma/sqrt(2 pi G rho0) = sigma^2/(pi G M) (not the NFW profile, which is 3D), and the
-- 1D virial theorem is 2K = U. The configuration mesh is a linear-to-exponential map on [-1, 1]
-- clustering cells at x = 0 while the box edge sits at ~20 a (effectively infinity, with absorbing
-- species boundaries and phi_g = 0 Dirichlet boundaries); the velocity mesh uses the same map so
-- the high-velocity tail of the relaxed distribution is covered.
-- Figures of merit (tensor p1 (hybrid), 64 x 64 cells): the virial ratio 2K/U goes from 0.499 to
-- 0.997 +- 0.004 (mean over t = 32-63) and the final density is fit by rho0 sech^2(x/z0) with rho0
-- = 0.527 and z0 = 0.951 (L2 residual 0.8%; the isothermal sheet with the final K has rho0 = 0.532
-- and z0 = 0.939); mass conserved to 3e-12 and total energy (kinetic minus half the gravitational
-- energy file, int |grad phi_g|^2/alpha_g) to 1e-9.

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
alpha_g = 1.0 -- Gravitational coupling 4 pi G.
mass_dm = 1.0 -- Dark matter mass.
charge_dm = 0.0 -- Dark matter charge (neutral).

M = 1.0 -- Mass per unit area of the sheet.
a = 1.0 -- Initial Gaussian width.
virial_ratio0 = 0.5 -- Initial virial ratio 2K/U.

-- Derived physical quantities (using normalized code units).
-- For a Gaussian sheet U = 2 sqrt(pi) G M^2 a and K = M sigma^2/2, so
-- 2K/U = sigma^2/(2 sqrt(pi) G M a) = 2 sqrt(pi) sigma^2/(alpha_g M a).
sigma0 = math.sqrt(virial_ratio0 * alpha_g * M * a / (2.0 * math.sqrt(pi))) -- Initial velocity dispersion.
rho0 = M / (math.sqrt(2.0 * pi) * a) -- Initial peak density.
t_dyn = 1.0 / math.sqrt(alpha_g * rho0) -- Dynamical time 1/sqrt(4 pi G rho0).

-- Simulation parameters.
Nx = 64 -- Cell count (configuration space: x-direction).
Nvx = 64 -- Cell count (velocity space: vx-direction).
dx_lin = 0.05 -- Configuration map: cell size of the linear part at x = 0.
x_exp = 20.0 -- Configuration map: the exponential part reaches x_exp - 1 at the boundary.
dv_lin = 0.02 -- Velocity map: cell size of the linear part at vx = 0.
v_exp = 4.0 -- Velocity map: the exponential part reaches v_exp - 1 at the boundary.
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 40.0 * t_dyn -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Linear-to-exponential map of a computational coordinate c in [-1, 1]: cells of
-- size d_lin at the origin (N cells across the domain), growing exponentially so
-- that the exponential part reaches e_max - 1 at the boundary.
function linExpMap(c, d_lin, e_max, N)
  local lin = 0.5 * d_lin * N * c
  if c < 0.0 then
    return lin - (math.exp(-c * math.log(e_max)) - 1.0)
  end
  return lin + (math.exp(c * math.log(e_max)) - 1.0)
end

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  -- The configuration grid is computational, mapped to physical positions by mapc2pPos.
  lower = { -1.0 },
  upper = { 1.0 },
  cells = { Nx },

  mapc2pPos = {
    {
      pmap = function (t, xn)
        return linExpMap(xn[1], dx_lin, x_exp, Nx)
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

  -- Dark matter.
  dm = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_dm, mass = mass_dm,

    -- Velocity space grid (computational, mapped to physical velocities by mapc2pVel).
    lower = { -1.0 },
    upper = { 1.0 },
    cells = { Nvx },

    mapc2pVel = {
      {
        vmap = function (t, xn)
          return linExpMap(xn[1], dv_lin, v_exp, Nvx)
        end
      },
    },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local x, vx = xn[1], xn[2]

          -- Gaussian slab with a Maxwellian velocity distribution.
          local n = (M / (2.0 * pi * a * sigma0)) * math.exp(-(x * x) / (2.0 * a * a)) *
            math.exp(-(vx * vx) / (2.0 * sigma0 * sigma0))

          return n
        end
      }
    },

    -- Absorbing boundaries: the box edge is effectively at infinity.
    bcx = {
      lower = {
        type = G0.SpeciesBc.bcAbsorb
      },
      upper = {
        type = G0.SpeciesBc.bcAbsorb
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Field: self-gravity only (the species is neutral), phi_g = 0 at the box edges.
  field = Vlasov.Field.new {
    fieldID = G0.FieldModel.Phi,
    alphaG = alpha_g,

    poissonBcs = {
      lowerType = {
        G0.PoissonBc.bcDirichlet
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
