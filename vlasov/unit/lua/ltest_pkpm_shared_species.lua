-- Test that a PKPM species definition can be shared by several apps without the apps affecting
-- each other: a PKPM species carries two initial conditions, of its distributions (called with the
-- configuration coordinates and the normalized velocity) and of its fluid state (called with the
-- configuration coordinates), and each app copies their callback contexts with its own coordinate
-- dimensions. Building a 2D app from the same species must not change how the live 1D app calls
-- them.

local Vlasov = G0.Vlasov

local dist_ndim, fluid_ndim = 0, 0

local neut = Vlasov.PKPMSpecies.new {
  charge = 0.0, mass = 1.0,

  kinetic = {
    lower = { -6.0 },
    upper = { 6.0 },
    cells = { 8 },

    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          dist_ndim = #xn
          local w = xn[#xn]

          local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * math.pi) -- Normalized Maxwellian in w.

          return F, F, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        end
      }
    }
  },

  fluid = {
    init = function (t, xn)
      fluid_ndim = #xn
      return 1.0, 0.0, 0.0, 0.0, 1.5
    end
  }
}

local function new_field()
  return Vlasov.Field.new {
    epsilon0 = 1.0, mu0 = 1.0,

    init = function (t, xn)
      return 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = false
  }
end

local app1d = Vlasov.App.new {
  tEnd = 1.0,
  nFrame = 1,
  lower = { 0.0 },
  upper = { 1.0 },
  cells = { 4 },
  polyOrder = 1,
  basis = "serendipity",
  periodicDirs = { 1 },
  neut = neut,
  field = new_field()
}
app1d:apply_ic()
assert(dist_ndim == 2, "1D app called the distribution initial condition with " .. dist_ndim .. " coordinates")
assert(fluid_ndim == 1, "1D app called the fluid initial condition with " .. fluid_ndim .. " coordinates")

local app2d = Vlasov.App.new {
  tEnd = 1.0,
  nFrame = 1,
  lower = { 0.0, 0.0 },
  upper = { 1.0, 1.0 },
  cells = { 4, 4 },
  polyOrder = 1,
  basis = "serendipity",
  periodicDirs = { 1, 2 },
  neut = neut,
  field = new_field()
}
app2d:apply_ic()
assert(dist_ndim == 3, "2D app called the distribution initial condition with " .. dist_ndim .. " coordinates")
assert(fluid_ndim == 2, "2D app called the fluid initial condition with " .. fluid_ndim .. " coordinates")

-- The first app must be unaffected by the second.
app1d:apply_ic()
assert(dist_ndim == 2, "constructing a 2D app changed the 1D app's distribution callback to " .. dist_ndim .. "D")
assert(fluid_ndim == 1, "constructing a 2D app changed the 1D app's fluid callback to " .. fluid_ndim .. "D")
