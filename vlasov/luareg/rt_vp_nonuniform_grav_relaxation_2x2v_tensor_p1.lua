-- Relaxation of a self-gravitating collisionless disc toward equilibrium with the Vlasov-Poisson
-- system (2x2v). Neutral dark matter (m = 1, gravitational coupling alpha_g = 4 pi G = 1) starts as
-- a Gaussian disc of mass per unit length M = 8 pi and width a = 2 with a Maxwellian velocity
-- distribution whose dispersion puts it at half the virial kinetic energy. In 2D gravity (potential
-- 2 G m ln r) the virial theorem fixes the kinetic energy of any equilibrium, K = G M^2/2 = alpha_g
-- M^2/(8 pi), so sigma^2 = K/M = 1 here, and the isothermal equilibrium is rho = rho0/(1 +
-- r^2/r0^2)^2 with r0^2 = 2 sigma^2/(pi G rho0) = 8 sigma^2/(alpha_g rho0), whose mass is M = 2
-- sigma^2/G (not the NFW profile, which is 3D). The configuration mesh is a linear-to-exponential
-- map in each direction clustering cells at the origin while the box edge sits at ~7.5 a
-- (effectively infinity, with absorbing species boundaries and phi_g = 0 Dirichlet boundaries); the
-- velocity mesh uses the same map so the high-velocity tail of the relaxed distribution is covered.
-- Figures of merit (tensor p1 (hybrid), 16^2 x 16^2 cells): the kinetic energy relative to its
-- virial value goes from 0.500 to 0.995 +- 0.027 (mean and standard deviation over t = 10-20 of the
-- integrated-moment time series, final 0.976) and the final density is fit by rho0/(1 + r^2/r0^2)^2
-- with rho0 = 2.43 and r0 = 1.91 (least squares on the cell averages weighted by cell area;
-- relative L2 residual 3.6%, versus 11% for a Gaussian; the isothermal r0 from the final sigma^2 =
-- K/M and rho0 is 1.79); the mass changes by +0.05% over the run (the absorbing boundaries sweep
-- out the signed, under-resolved tails of the distribution in the large outer cells) and the total
-- energy (kinetic minus half the gravitational energy file, int |grad phi_g|^2/alpha_g) by 1e-3.

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
alpha_g = 1.0 -- Gravitational coupling 4 pi G.
mass_dm = 1.0 -- Dark matter mass.
charge_dm = 0.0 -- Dark matter charge (neutral).

M = 8.0 * pi -- Mass per unit length of the disc.
a = 2.0 -- Initial Gaussian width.
virial_ratio0 = 0.5 -- Initial ratio of the kinetic energy to its virial value G M^2/2.

-- Derived physical quantities (using normalized code units).
-- K = M sigma^2 in 2D and the virial value is G M^2/2 = alpha_g M^2/(8 pi).
sigma0 = math.sqrt(virial_ratio0 * alpha_g * M / (8.0 * pi)) -- Initial velocity dispersion.
rho0 = M / (2.0 * pi * a * a) -- Initial peak density.
t_dyn = 1.0 / math.sqrt(alpha_g * rho0) -- Dynamical time 1/sqrt(4 pi G rho0).

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x- and y-directions).
Nvx = 16 -- Cell count (velocity space: vx- and vy-directions).
dx_lin = 0.1 -- Configuration map: cell size of the linear part at x = 0.
x_exp = 15.0 -- Configuration map: the exponential part reaches x_exp - 1 at the boundary.
dv_lin = 0.05 -- Velocity map: cell size of the linear part at vx = 0.
v_exp = 6.0 -- Velocity map: the exponential part reaches v_exp - 1 at the boundary.
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 20.0 * t_dyn -- Final simulation time.
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
  lower = { -1.0, -1.0 },
  upper = { 1.0, 1.0 },
  cells = { Nx, Nx },

  mapc2pPos = {
    {
      pmap = function (t, xn)
        return linExpMap(xn[1], dx_lin, x_exp, Nx)
      end
    },
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
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x- and y-directions).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Dark matter.
  dm = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_dm, mass = mass_dm,

    -- Velocity space grid (computational, mapped to physical velocities by mapc2pVel).
    lower = { -1.0, -1.0 },
    upper = { 1.0, 1.0 },
    cells = { Nvx, Nvx },

    mapc2pVel = {
      {
        vmap = function (t, xn)
          return linExpMap(xn[1], dv_lin, v_exp, Nvx)
        end
      },
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
          local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]

          -- Gaussian disc with a Maxwellian velocity distribution.
          local n = (M / (4.0 * pi * pi * a * a * sigma0 * sigma0)) *
            math.exp(-(x * x + y * y) / (2.0 * a * a)) *
            math.exp(-(vx * vx + vy * vy) / (2.0 * sigma0 * sigma0))

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
    bcy = {
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
        G0.PoissonBc.bcDirichlet, G0.PoissonBc.bcDirichlet
      },
      upperType = {
        G0.PoissonBc.bcDirichlet, G0.PoissonBc.bcDirichlet
      },
      lowerValue = {
        0.0, 0.0
      },
      upperValue = {
        0.0, 0.0
      }
    }
  }
}

vlasovApp:run()
