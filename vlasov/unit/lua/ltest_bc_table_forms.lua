-- Forms of the boundary-condition tables in the Lua wrapper. The field and the fluid species take
-- the array form { lower, upper }; the kinetic species take the table form { lower = { type = ...
-- }, upper = { type = ... } } (either side may be absent). A table in the other form used to be
-- read as no boundary condition (copy on both sides) without an error: each constructor must
-- accept its own form and raise an error for the other.

local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

local function new_field(bcx)
  return Vlasov.Field.new {
    epsilon0 = 1.0, mu0 = 1.0,

    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    bcx = bcx
  }
end

local function new_fluid(bcx)
  return Vlasov.FluidSpecies.new {
    equation = Euler.new { gasGamma = 1.4 },

    init = function (t, xn)
      return 1.0, 0.0, 0.0, 0.0, 1.0
    end,

    bcx = bcx
  }
end

local function new_species(bcx)
  return Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = -1.0, mass = 1.0,

    lower = { -6.0 },
    upper = { 6.0 },
    cells = { 8 },

    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return 1.0
        end,
        temperatureInit = function (t, xn)
          return 1.0
        end,
        driftVelocityInit = function (t, xn)
          return 0.0
        end
      }
    },

    bcx = bcx
  }
end

local field_array = { G0.FieldBc.bcSymWall, G0.FieldBc.bcCopy }
local field_table = {
  lower = { type = G0.FieldBc.bcSymWall },
  upper = { type = G0.FieldBc.bcCopy }
}
local species_array = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcCopy }
local species_table = {
  lower = { type = G0.SpeciesBc.bcReflect },
  upper = { type = G0.SpeciesBc.bcCopy }
}

-- Field: the array form with both sides.
assert(pcall(new_field, field_array), "Field: the array form was rejected")
assert(pcall(new_field, nil), "Field: no boundary condition was rejected")
assert(not pcall(new_field, field_table), "Field: the table form was accepted")
assert(not pcall(new_field, { G0.FieldBc.bcSymWall }), "Field: a single side was accepted")

-- Fluid species: the array form with both sides.
assert(pcall(new_fluid, species_array), "FluidSpecies: the array form was rejected")
assert(pcall(new_fluid, nil), "FluidSpecies: no boundary condition was rejected")
assert(not pcall(new_fluid, species_table), "FluidSpecies: the table form was accepted")
assert(
  not pcall(new_fluid, { G0.SpeciesBc.bcReflect }), "FluidSpecies: a single side was accepted"
)

-- Kinetic species: the table form, either side optional.
assert(pcall(new_species, species_table), "Species: the table form was rejected")
assert(
  pcall(new_species, { upper = { type = G0.SpeciesBc.bcCopy } }), "Species: one side was rejected"
)
assert(pcall(new_species, nil), "Species: no boundary condition was rejected")
assert(not pcall(new_species, species_array), "Species: the array form was accepted")
