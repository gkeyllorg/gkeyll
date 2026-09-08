#ifdef GKYL_HAVE_LUA

#include <gkyl_alloc.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_lw.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_lua_utils.h>
#include <gkyl_lw_priv.h>
#include <gkyl_null_comm.h>
#include <gkyl_zero_lw.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
#include <limits.h>

#include <string.h>

#include <stc/coption.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif



// Gyrokinetic FEM boundary conditions -> enum map.
static const struct gkyl_str_int_pair parproj_type[] = {
  { "None", GKYL_FEM_PARPROJ_NONE },
  { "Periodic", GKYL_FEM_PARPROJ_PERIODIC },
  { "DirichletGhost", GKYL_FEM_PARPROJ_DIRICHLET_GHOST },
  { "DirichletSkin", GKYL_FEM_PARPROJ_DIRICHLET_SKIN },
  { 0, 0 }
};

// Gyrokinetic geometry type -> enum map.
static const struct gkyl_str_int_pair geometry_type[] = {
  { "Tokamak", GKYL_GEOMETRY_TOKAMAK },
  { "Mirror", GKYL_GEOMETRY_MIRROR },
  { "MapC2P", GKYL_GEOMETRY_MAPC2P },
  { "FromFile", GKYL_GEOMETRY_FROMFILE },
  { 0, 0 }
};

// Gyrokinetic position map type -> enum map.
static const struct gkyl_str_int_pair position_map_type[] = {
  { "UserInput", GKYL_PMAP_USER_INPUT },
  { "ConstantPolynomial", GKYL_PMAP_CONSTANT_DB_POLYNOMIAL },
  { "ConstantNumeric", GKYL_PMAP_CONSTANT_DB_NUMERIC },
  { 0, 0 }
};

// Gyrokinetic collisionless type -> enum map.
static const struct gkyl_str_int_pair gk_collisionless_type[] = {
  { "None", GKYL_GK_COLLISIONLESS_NONE },
  { "GKCollisionlessES", GKYL_GK_COLLISIONLESS_ES },
  { "GKCollisionlessEM", GKYL_GK_COLLISIONLESS_EM },
  { 0, 0 }
};

// Gyrokinetic field type -> enum map.
static const struct gkyl_str_int_pair gk_field_type[] = {
  { "Default", GKYL_GK_FIELD },
  { "Boltzmann", GKYL_GK_FIELD_BOLTZMANN },
  { "Adiabatic", GKYL_GK_FIELD_ADIABATIC },
  { "Electromagnetic", GKYL_GK_FIELD_EM },
  { 0, 0 }
};

// Gyrokinetic radiation type -> enum map.
static const struct gkyl_str_int_pair gk_radiation_type[] = {
  { "None", GKYL_NO_RADIATION },
  { "GKRadiation", GKYL_GK_RADIATION },
  { "VMComptonRadiation", GKYL_VM_COMPTON_RADIATION },
  { 0, 0 }
};

// Gyrokinetic radiation Te model type -> enum map.
static const struct gkyl_str_int_pair gk_radiation_te_type[] = {
  { "Conservative", GKYL_VARY_TE_CONSERVATIVE },
  { "Aggressive", GKYL_VARY_TE_AGGRESSIVE },
  { "Const", GKYL_CONST_TE },
  { 0, 0 }
};

// Gyrokinetic reaction type -> enum map.
static const struct gkyl_str_int_pair gk_react_type[] = {
  { "None", GKYL_NO_REACT },
  { "Ionization", GKYL_REACT_IZ },
  { "ChargeExchange", GKYL_REACT_CX },
  { "Recombination", GKYL_REACT_RECOMB },
  { 0, 0 }
};

// Gyrokinetic ion type -> enum map.
static const struct gkyl_str_int_pair gk_ion_type[] = {
  { "Hydrogen", GKYL_ION_H },
  { "Deuterium", GKYL_ION_D },
  { "Helium", GKYL_ION_HE },
  { "Lithium", GKYL_ION_LI },
  { "Beryllium", GKYL_ION_BE },
  { "Boron", GKYL_ION_B },
  { "Carbon", GKYL_ION_C },
  { "Nitrogen", GKYL_ION_N },
  { "Oxygen", GKYL_ION_O },
  { "Neon", GKYL_ION_NE },
  { "Argon", GKYL_ION_AR },
  { 0, 0 }
};

// Gyrokinetic self-reaction type -> enum map.
static const struct gkyl_str_int_pair gk_react_self_type[] = {
  { "Electron", GKYL_SELF_ELC },
  { "Ion", GKYL_SELF_ION },
  { "Donor", GKYL_SELF_DONOR },
  { "Receiver", GKYL_SELF_RECVR },
  { "Partner", GKYL_SELF_PARTNER },
  { 0, 0 }
};

// Species boundary conditions -> enum map.
static const struct gkyl_str_int_pair gk_bcs[] = {
  // Particle BCs.
  { "speciesSkip", GKYL_BC_GK_SKIP }, // Do not apply any BCs
  { "speciesCopy", GKYL_BC_GK_SPECIES_COPY }, // Copy skin into ghost.
  { "speciesReflect", GKYL_BC_GK_SPECIES_REFLECT }, // Reflect particles.
  { "speciesAbsorb", GKYL_BC_GK_SPECIES_ABSORB }, // Absorbing BCs.
  { "speciesFunc", GKYL_BC_GK_SPECIES_FUNC }, // Fill ghost cell using a user-function.
  { "speciesFixedFunc", GKYL_BC_GK_SPECIES_FIXED_FUNC }, // Fixed function, time-independent.
  { "speciesZeroFlux", GKYL_BC_GK_SPECIES_ZERO_FLUX }, // Zero flux.
  { "speciesSheath", GKYL_BC_GK_SPECIES_SHEATH }, // Sheath.
  { "speciesRecycle", GKYL_BC_GK_SPECIES_RECYCLE }, // Recycling.
  { "speciesPeriodic", GKYL_BC_GK_SPECIES_PERIODIC }, // Periodic.
  { "speciesTwistshift", GKYL_BC_GK_SPECIES_TWISTSHIFT }, // Twist-shift.
  // Field BCs.
  { "fieldPeriodic", GKYL_BC_GK_FIELD_PERIODIC }, // Periodic.
  { "fieldDirichlet", GKYL_BC_GK_FIELD_DIRICHLET }, // Dirichlet.
  { "fieldNeumann", GKYL_BC_GK_FIELD_NEUMANN }, // Nemann.
  { "fieldDirichletVarying", GKYL_BC_GK_FIELD_DIRICHLET_VARYING }, // Spatially varying Dirichlet.
  { "fieldBoundaryValue", GKYL_BC_GK_FIELD_BOUNDARY_VALUE }, // Skin value at the boundary.
  { "fieldTwistshift", GKYL_BC_GK_FIELD_TWISTSHIFT }, // Twist-shift.
  { 0, 0 }
};


void
gkyl_register_gyrokinetic_fem_bc_types(lua_State *L)
{
  register_types(L, parproj_type, "ParProjBc");
}

void
gkyl_register_gyrokinetic_geometry_types(lua_State *L)
{
  register_types(L, geometry_type, "Geometry");
}

void
gkyl_register_gyrokinetic_position_map_types(lua_State *L)
{
  register_types(L, position_map_type, "PositionMap");
}

void
gkyl_register_gyrokinetic_collisionless_types(lua_State *L)
{
  register_types(L, gk_collisionless_type, "GKCollisionless");
}

void
gkyl_register_gyrokinetic_field_types(lua_State *L)
{
  register_types(L, gk_field_type, "GKField");
}

void
gkyl_register_gyrokinetic_radiation_types(lua_State *L)
{
  register_types(L, gk_radiation_type, "Radiation");
}

void
gkyl_register_gyrokinetic_radiation_Te_types(lua_State *L)
{
  register_types(L, gk_radiation_te_type, "TeMinModel");
}

void
gkyl_register_gyrokinetic_reaction_types(lua_State *L)
{
  register_types(L, gk_react_type, "Reaction");
}

void
gkyl_register_gyrokinetic_ion_types(lua_State *L)
{
  register_types(L, gk_ion_type, "Ion");
}

void
gkyl_register_gyrokinetic_self_reaction_types(lua_State *L)
{
  register_types(L, gk_react_self_type, "Self");
}

void
gkyl_register_gyrokinetic_bc_types(lua_State *L)
{
  register_types(L, gk_bcs, "GyrokineticBc");
}

// Magic IDs for use in distinguishing various species and field types.
enum gyrokinetic_magic_ids {
  GYROKINETIC_SPECIES_DEFAULT = 100, // Standard gyrokinetic species.
  GYROKINETIC_NEUTRAL_SPECIES_DEFAULT, // Neutral gyrokinetic species.
  GYROKINETIC_FIELD_DEFAULT, // Gyrokinetic Poisson equation.
};

/* *************** */
/* Species methods */
/* *************** */

// Metatable name for species input struct.
#define GYROKINETIC_SPECIES_METATABLE_NM "GkeyllZero.App.Gyrokinetic.Species"

// Lua userdata object for constructing species input.
struct gyrokinetic_species_lw {
  int magic; // This must be first element in the struct.
  
  struct gkyl_gyrokinetic_species gk_species; // Input struct to construct species.
  int vdim; // Velocity space dimensions.
  bool evolve; // Is this species evolved?

  bool has_mapc2p_mapping_func; // Is there a non-uniform velocity space mapping function?
  struct lua_func_ctx mapc2p_mapping_func_ref; // Lua registry reference to non-uniofrm velocity space mapping function.

  enum gkyl_projection_id proj_id; // Projection type.

  bool has_init_func; // Is there an initialization function?
  struct lua_func_ctx init_func_ref; // Lua registry reference to initialization function.

  bool has_density_init_func; // Is there a density initialization function?
  struct lua_func_ctx density_init_func_ref; // Lua registry reference to density initialization function.

  bool has_Upar_init_func; // Is there a parallel velocity initialiation function?
  struct lua_func_ctx Upar_init_func_ref; // Lua registry reference to parallel velocity initialization function.

  bool has_temp_init_func; // Is there a temperature initialization function?
  struct lua_func_ctx temp_init_func_ref; // Lua registry reference to temperature initialization function.

  bool has_par_temp_init_func; // Is there a parallel temperature initialization function?
  struct lua_func_ctx par_temp_init_func_ref; // Lua registry reference to parallel temperature initialization function.

  bool has_perp_temp_init_func; // Is there a perpendicular temperature initialization function?
  struct lua_func_ctx perp_temp_init_func_ref; // Lua registry reference to perpendicular temperature initialization function.

  bool proj_correct_all_moms; // Are we correcting all moments in projection, or only density?

  enum gkyl_gk_collisionless_type collisionless_type; // Collisionless type.

  enum gkyl_collision_id collision_id; // Collision type.
  
  bool has_self_nu_func; // Is there a self-collision frequency function?
  struct lua_func_ctx self_nu_func_ref; // Lua registry reference to self-collision frequency function.

  int num_cross_collisions; // Number of species that we cross-collide with.
  char collide_with[GKYL_MAX_SPECIES][128]; // Names of species that we cross-collide with.

  bool collision_norm_nu; // Are we rescaling the collision frequency?
  double collision_den_ref; // Density used to calculate Coulomb logarithm for collision frequency.
  double collision_temp_ref; // Temperature used to calculate Coulomb logarithm for collision frequency.
  double collision_hbar; // Reduced Planck's constant for calculating collision frequency.
  double collision_eps0; // Vacuum permittivity for calculating collision frequency.
  double collision_eV; // Elementary charge for calculating collision frequency.

  bool correct_all_moms; // Are we correcting all moments in collisions, or only density?
  double iter_eps; // Error tolerance for moment fixes in collisions (density is always exact).
  int max_iter; // Maximum number of iterations for moment fixes in collisions.
  bool use_last_converged; // Use last iteration value in collisions regardless of convergence?

  enum gkyl_source_id source_id; // Source type.

  int num_sources; // Number of projection objects in source.
  enum gkyl_projection_id source_proj_id[GKYL_MAX_PROJ]; // Projection type in source.

  bool source_has_init_func[GKYL_MAX_PROJ]; // Is there an initialization function in source?
  struct lua_func_ctx source_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to initialization function in source.

  bool source_has_density_init_func[GKYL_MAX_PROJ]; // Is there a density initialization function in source?
  struct lua_func_ctx source_density_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to density initialization function in source.

  bool source_has_Upar_init_func[GKYL_MAX_PROJ]; // Is there a parallel velocity initialization function in source?
  struct lua_func_ctx source_Upar_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to parallel velocity initialization function in source.

  bool source_has_temp_init_func[GKYL_MAX_PROJ]; // Is there a temperature initialization function in source?
  struct lua_func_ctx source_temp_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to temperature initialization function in source.

  enum gkyl_radiation_id radiation_id; // Radiation type.

  int radiation_num_cross_collisions; // Number of radiation species that we cross-collide with.
  char radiation_collide_with[GKYL_MAX_SPECIES][128]; // Names of radiation species that we cross-collide with.

  int radiation_z[GKYL_MAX_SPECIES]; // Atomic Z of radiation species that we are colliding with.
  int radiation_charge_state[GKYL_MAX_SPECIES]; // Charge state of radiation species that we are colliding with.
  int radiation_num_of_densities[GKYL_MAX_SPECIES]; // Maximum number of densities to use per charge state of radiation species that we are colliding with.

  enum gkyl_te_min_model radiation_te_min_model; // How is the radiation turned off (constant, or with varying electron temperature)?
  double radiation_Te_min; // Minimum temperature (in J) at which to stop radiating.

  int num_react; // Number of reaction types.

  enum gkyl_react_id react_id[GKYL_MAX_REACT]; // What type of reaction (ionization, charge exchange, recombination)?
  enum gkyl_react_self_type react_type_self[GKYL_MAX_REACT]; // What is the role of the species in this reaction?
  enum gkyl_ion_type react_ion_id[GKYL_MAX_REACT]; // What type of ion is reacting?

  char react_elc_nm[GKYL_MAX_REACT][128]; // Name of electron species in the reaction.
  char react_ion_nm[GKYL_MAX_REACT][128]; // Name of ion species in the reaction.
  char react_donor_nm[GKYL_MAX_REACT][128]; // Name of donor species in the reaction.
  char react_recvr_nm[GKYL_MAX_REACT][128]; // Name of receiver species in the reaction.

  int react_charge_state[GKYL_MAX_REACT]; // Charge state of species in the reaction.
  double react_ion_mass[GKYL_MAX_REACT]; // Mass of ion species in the reaction.
  double react_elc_mass[GKYL_MAX_REACT]; // Mass of electron species in the reaction.

  int num_neut_react; // Number of neutral reaction types.

  enum gkyl_react_id neut_react_id[GKYL_MAX_REACT]; // What type of neutral reaction (ionization, charge exchange, recombination)?
  enum gkyl_react_self_type neut_react_type_self[GKYL_MAX_REACT]; // What is the role of the species in this neutral reaction?
  enum gkyl_ion_type neut_react_ion_id[GKYL_MAX_REACT]; // What type of ion in the neutral reaction?

  char neut_react_elc_nm[GKYL_MAX_REACT][128]; // Name of electron species in the neutral reaction.
  char neut_react_ion_nm[GKYL_MAX_REACT][128]; // Name of ion species in the neutral reaction.
  char neut_react_donor_nm[GKYL_MAX_REACT][128]; // Name of donor species in the neutral reaction.
  char neut_react_recvr_nm[GKYL_MAX_REACT][128]; // Name of receiver species in the neutral reaction.
  char neut_react_partner_nm[GKYL_MAX_REACT][128]; // Name of partner species in the neutral reaction.

  int neut_react_charge_state[GKYL_MAX_REACT]; // Charge state of species in the neutral reaction.
  double neut_react_ion_mass[GKYL_MAX_REACT]; // Mass of ion species in the neutral reaction.
  double neut_react_elc_mass[GKYL_MAX_REACT]; // Mass of electron species in the neutral reaction.
};

static int
gyrokinetic_species_lw_new(lua_State *L)
{
  int vdim = 0;
  struct gkyl_gyrokinetic_species gk_species = { };
  
  gk_species.charge = glua_tbl_get_number(L, "charge", 0.0);
  gk_species.mass = glua_tbl_get_number(L, "mass", 1.0);
  gk_species.polarization_density = glua_tbl_get_number(L, "polarizationDensity", 0.0);

  with_lua_tbl_tbl(L, "cells") {
    vdim = glua_objlen(L);

    for (int d = 0; d < vdim; d++) {
      gk_species.cells[d] = glua_tbl_iget_integer(L, d + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "lower") {
    for (int d = 0; d < vdim; d++) {
      gk_species.lower[d] = glua_tbl_iget_number(L, d + 1, 0.0);
    }
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d = 0; d < vdim; d++) {
      gk_species.upper[d] = glua_tbl_iget_number(L, d + 1, 0.0);
    }
  }

  bool evolve = glua_tbl_get_integer(L, "evolve", true);

  with_lua_tbl_tbl(L, "diagnostics") {
    int num_diag_moments = glua_objlen(L);

    for (int i = 0; i < num_diag_moments; i ++) {
      gk_species.diag_moments[i] = glua_tbl_iget_integer(L, i+1, 0);
    }

    gk_species.num_diag_moments = num_diag_moments;
  }

  with_lua_tbl_tbl(L, "bcs") {
    int num_bcs = glua_objlen(L);
    for (int i = 0; i < num_bcs; i++) {
      gk_species.bcs[i].dir = glua_tbl_get_integer(L, "dir", 0);
      gk_species.bcs[i].edge = glua_tbl_get_integer(L, "edge", 0);
      gk_species.bcs[i].type = glua_tbl_get_integer(L, "type", 0);
      with_lua_tbl_tbl(L, "value") {
        int num_vals = glua_objlen(L);
        for (int k = 0; k < num_vals; k++) {
          gk_species.bcs[i].value[k] = glua_tbl_iget_number(L, k + 1, 0.0);;
        }
      }
    }
  }

  bool has_mapc2p_mapping_func = false;
  int mapc2p_mapping_func_ref = LUA_NOREF;

  with_lua_tbl_tbl(L, "mapc2p") {
    if (glua_tbl_get_func(L, "mapping")) {
      mapc2p_mapping_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_mapc2p_mapping_func = true;
    }
  };
  
  enum gkyl_projection_id proj_id = GKYL_PROJ_FUNC;

  bool has_init_func = false;
  int init_func_ref = LUA_NOREF;

  bool has_density_init_func = false;
  int density_init_func_ref = LUA_NOREF;

  bool has_Upar_init_func = false;
  int Upar_init_func_ref = LUA_NOREF;

  bool has_temp_init_func = false;
  int temp_init_func_ref = LUA_NOREF;

  bool has_par_temp_init_func = false;
  int par_temp_init_func_ref = LUA_NOREF;

  bool has_perp_temp_init_func = false;
  int perp_temp_init_func_ref = LUA_NOREF;

  bool proj_correct_all_moms = false;

  with_lua_tbl_tbl(L, "projection") {
    proj_id = glua_tbl_get_integer(L, "projectionID", 0);

    if (glua_tbl_get_func(L, "init")) {
      init_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_init_func = true;
    }

    if (glua_tbl_get_func(L, "densityInit")) {
      density_init_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_density_init_func = true;
    }

    if (glua_tbl_get_func(L, "parallelVelocityInit")) {
      Upar_init_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_Upar_init_func = true;
    }

    if (glua_tbl_get_func(L, "temperatureInit")) {
      temp_init_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_temp_init_func = true;
    }

    if (glua_tbl_get_func(L, "parallelTemperatureInit")) {
      par_temp_init_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_par_temp_init_func = true;
    }

    if (glua_tbl_get_func(L, "perpendicularTemperatureInit")) {
      perp_temp_init_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_perp_temp_init_func = true;
    }

    proj_correct_all_moms = glua_tbl_get_bool(L, "correctAllMoments", false);
  }

  bool correct_all_moms = false;
  double iter_eps = 0.0;
  int max_iter = 0; 
  bool use_last_converged = true; 
  with_lua_tbl_tbl(L, "correct") {
    correct_all_moms = glua_tbl_get_bool(L, "correctAllMoments", false);
    iter_eps = glua_tbl_get_number(L, "iterationEpsilon", 0.0);
    max_iter = glua_tbl_get_integer(L, "maxIterations", 0);
    use_last_converged = glua_tbl_get_bool(L, "useLastConverged", true);
  }

  enum gkyl_gk_collisionless_type collisionless_type = GKYL_GK_COLLISIONLESS_NONE;

  with_lua_tbl_tbl(L, "collisionless") {
    collisionless_type = glua_tbl_get_integer(L, "type", 0);
  }

  enum gkyl_collision_id collision_id = GKYL_NO_COLLISIONS;

  bool has_self_nu_func = false;
  int self_nu_func_ref = LUA_NOREF;

  int num_cross_collisions = 0;
  char collide_with[GKYL_MAX_SPECIES][128];

  bool collision_norm_nu = false;
  double collision_den_ref = 0.0;
  double collision_temp_ref = 0.0;
  double collision_hbar = 0.0;
  double collision_eps0 = 0.0;
  double collision_eV = 0.0;

  with_lua_tbl_tbl(L, "collisions") {
    collision_id = glua_tbl_get_integer(L, "collisionID", 0);

    if (glua_tbl_get_func(L, "selfNu")) {
      self_nu_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_self_nu_func = true;
    }

    num_cross_collisions = glua_tbl_get_integer(L, "numCrossCollisions", 0);
    with_lua_tbl_tbl(L, "collideWith") {
      for (int i = 0; i < num_cross_collisions; i++) {
        const char* collide_with_char = glua_tbl_iget_string(L, i + 1, "");
        strcpy(collide_with[i], collide_with_char);
      }
    }

    collision_norm_nu = glua_tbl_get_bool(L, "normalizeNu", false);
    collision_den_ref = glua_tbl_get_number(L, "referenceDensity", 0.0);
    collision_temp_ref = glua_tbl_get_number(L, "referenceTemperature", 0.0);
    collision_hbar = glua_tbl_get_number(L, "hbar", 0.0);
    collision_eps0 = glua_tbl_get_number(L, "epsilon0", 0.0);
    collision_eV = glua_tbl_get_number(L, "eV", 0.0);
  }

  enum gkyl_source_id source_id = GKYL_NO_SOURCE;

  int num_sources = 0;
  enum gkyl_projection_id source_proj_id[GKYL_MAX_PROJ];

  bool source_has_init_func[GKYL_MAX_PROJ];
  int source_init_func_ref[GKYL_MAX_PROJ];

  bool source_has_density_init_func[GKYL_MAX_PROJ];
  int source_density_init_func_ref[GKYL_MAX_PROJ];

  bool source_has_Upar_init_func[GKYL_MAX_PROJ];
  int source_Upar_init_func_ref[GKYL_MAX_PROJ];

  bool source_has_temp_init_func[GKYL_MAX_PROJ];
  int source_temp_init_func_ref[GKYL_MAX_PROJ];

  with_lua_tbl_tbl(L, "source") {
    source_id = glua_tbl_get_integer(L, "sourceID", 0);
    num_sources = glua_tbl_get_integer(L, "numSources", 0);

    with_lua_tbl_tbl(L, "projections") {
      for (int i = 0; i < num_sources; i++) {
        if (glua_tbl_iget_tbl(L, i + 1)) {
          source_proj_id[i] = glua_tbl_get_integer(L, "projectionID", 0);

          source_init_func_ref[i] = LUA_NOREF;
          source_has_init_func[i] = false;
          if (glua_tbl_get_func(L, "init")) {
            source_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
            source_has_init_func[i] = true;
          }

          source_density_init_func_ref[i] = LUA_NOREF;
          source_has_density_init_func[i] = false;
          if (glua_tbl_get_func(L, "densityInit")) {
            source_density_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
            source_has_density_init_func[i] = true;
          }

          source_Upar_init_func_ref[i] = LUA_NOREF;
          source_has_Upar_init_func[i] = false;
          if (glua_tbl_get_func(L, "parallelVelocityInit")) {
            source_Upar_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
            source_has_Upar_init_func[i] = true;
          }

          source_temp_init_func_ref[i] = LUA_NOREF;
          source_has_temp_init_func[i] = false;
          if (glua_tbl_get_func(L, "temperatureInit")) {
            source_temp_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
            source_has_temp_init_func[i] = true;
          }

          lua_pop(L, 1);
        }
      }
    }
  }

  enum gkyl_radiation_id radiation_id = GKYL_NO_RADIATION;

  int radiation_num_cross_collisions = 0;
  char radiation_collide_with[GKYL_MAX_SPECIES][128];

  int radiation_z[GKYL_MAX_SPECIES];
  int radiation_charge_state[GKYL_MAX_SPECIES];
  int radiation_num_of_densities[GKYL_MAX_SPECIES];

  int radiation_te_min_model = GKYL_VARY_TE_CONSERVATIVE;
  double radiation_Te_min = 0.0;

  with_lua_tbl_tbl(L, "radiation") {
    radiation_id = glua_tbl_get_integer(L, "radiationID", GKYL_NO_RADIATION);

    radiation_num_cross_collisions = glua_tbl_get_integer(L, "numCrossCollisions", 0);
    with_lua_tbl_tbl(L, "collideWith") {
      for (int i = 0; i < radiation_num_cross_collisions; i++) {
        const char* radiation_collide_with_char = glua_tbl_iget_string(L, i + 1, "");
        strcpy(radiation_collide_with[i], radiation_collide_with_char);
      }
    }
    
    with_lua_tbl_tbl(L, "atomicZ") {
      for (int i = 0; i < radiation_num_cross_collisions; i++) {
        radiation_z[i] = glua_tbl_iget_integer(L, i + 1, 0);
      }
    }
    with_lua_tbl_tbl(L, "chargeState") {
      for (int i = 0; i < radiation_num_cross_collisions; i++) {
        radiation_charge_state[i] = glua_tbl_iget_integer(L, i + 1, 0);
      }
    }
    with_lua_tbl_tbl(L, "numDensities") {
      for (int i = 0; i < radiation_num_cross_collisions; i++) {
        radiation_num_of_densities[i] = glua_tbl_iget_integer(L, i + 1, 0);
      }
    }

    radiation_te_min_model = glua_tbl_get_integer(L, "TeMinModel", 0);
    radiation_Te_min = glua_tbl_get_number(L, "TeMin", 0.0);
  }

  int num_react = 0;
  enum gkyl_react_id react_id[GKYL_MAX_REACT];
  enum gkyl_react_self_type react_type_self[GKYL_MAX_REACT];
  enum gkyl_ion_type react_ion_id[GKYL_MAX_REACT];

  char react_elc_nm[GKYL_MAX_REACT][128];
  char react_ion_nm[GKYL_MAX_REACT][128];
  char react_donor_nm[GKYL_MAX_REACT][128];
  char react_recvr_nm[GKYL_MAX_REACT][128];

  int react_charge_state[GKYL_MAX_REACT];
  double react_ion_mass[GKYL_MAX_REACT];
  double react_elc_mass[GKYL_MAX_REACT];

  with_lua_tbl_tbl(L, "reaction") {
    num_react = glua_tbl_get_integer(L, "numReactions", 0);

    with_lua_tbl_tbl(L, "reactionTypes") {
      for (int i = 0; i < num_react; i++) {
        if (glua_tbl_iget_tbl(L, i + 1)) {
          react_id[i] = glua_tbl_get_integer(L, "reactionID", GKYL_NO_REACT);
          react_type_self[i] = glua_tbl_get_integer(L, "selfType", GKYL_SELF_ELC);
          react_ion_id[i] = glua_tbl_get_integer(L, "ionType", GKYL_ION_H);

          const char* react_elc_nm_char = glua_tbl_get_string(L, "electronName", "");
          strcpy(react_elc_nm[i], react_elc_nm_char);

          const char* react_ion_nm_char = glua_tbl_get_string(L, "ionName", "");
          strcpy(react_ion_nm[i], react_ion_nm_char);

          const char* react_donor_nm_char = glua_tbl_get_string(L, "donorName", "");
          strcpy(react_donor_nm[i], react_donor_nm_char);

          const char* react_recvr_nm_char = glua_tbl_get_string(L, "receiverName", "");
          strcpy(react_recvr_nm[i], react_recvr_nm_char);

          react_charge_state[i] = glua_tbl_get_integer(L, "chargeState", 0);
          react_ion_mass[i] = glua_tbl_get_number(L, "ionMass", 0.0);
          react_elc_mass[i] = glua_tbl_get_number(L, "electronMass", 0.0);

          lua_pop(L, 1);
        }
      }
    }
  }

  int num_neut_react = 0;
  enum gkyl_react_id neut_react_id[GKYL_MAX_REACT];
  enum gkyl_react_self_type neut_react_type_self[GKYL_MAX_REACT];
  enum gkyl_ion_type neut_react_ion_id[GKYL_MAX_REACT];

  char neut_react_elc_nm[GKYL_MAX_REACT][128];
  char neut_react_ion_nm[GKYL_MAX_REACT][128];
  char neut_react_donor_nm[GKYL_MAX_REACT][128];
  char neut_react_recvr_nm[GKYL_MAX_REACT][128];
  char neut_react_partner_nm[GKYL_MAX_REACT][128];

  int neut_react_charge_state[GKYL_MAX_REACT];
  double neut_react_ion_mass[GKYL_MAX_REACT];
  double neut_react_elc_mass[GKYL_MAX_REACT];

  with_lua_tbl_tbl(L, "neutralReaction") {
    num_neut_react = glua_tbl_get_integer(L, "numReactions", 0);

    with_lua_tbl_tbl(L, "reactionTypes") {
      for (int i = 0; i < num_neut_react; i++) {
        if (glua_tbl_iget_tbl(L, i + 1)) {
          neut_react_id[i] = glua_tbl_get_integer(L, "reactionID", GKYL_NO_REACT);
          neut_react_type_self[i] = glua_tbl_get_integer(L, "selfType", GKYL_SELF_ELC);
          neut_react_ion_id[i] = glua_tbl_get_integer(L, "ionType", GKYL_ION_H);

          const char* neut_react_elc_nm_char = glua_tbl_get_string(L, "electronName", "");
          strcpy(neut_react_elc_nm[i], neut_react_elc_nm_char);

          const char* neut_react_ion_nm_char = glua_tbl_get_string(L, "ionName", "");
          strcpy(neut_react_ion_nm[i], neut_react_ion_nm_char);

          const char* neut_react_donor_nm_char = glua_tbl_get_string(L, "donorName", "");
          strcpy(neut_react_donor_nm[i], neut_react_donor_nm_char);

          const char* neut_react_recvr_nm_char = glua_tbl_get_string(L, "receiverName", "");
          strcpy(neut_react_recvr_nm[i], neut_react_recvr_nm_char);

          const char* neut_react_partner_nm_char = glua_tbl_get_string(L, "partnerName", "");
          strcpy(neut_react_partner_nm[i], neut_react_partner_nm_char);

          neut_react_charge_state[i] = glua_tbl_get_integer(L, "chargeState", 0);
          neut_react_ion_mass[i] = glua_tbl_get_number(L, "ionMass", 0.0);
          neut_react_elc_mass[i] = glua_tbl_get_number(L, "electronMass", 0.0);

          lua_pop(L, 1);
        }
      }
    }
  }
  
  struct gyrokinetic_species_lw *gks_lw = lua_newuserdata(L, sizeof(*gks_lw));
  gks_lw->magic = GYROKINETIC_SPECIES_DEFAULT;
  gks_lw->vdim = vdim;
  gks_lw->evolve = evolve;
  gks_lw->gk_species = gk_species;

  gks_lw->has_mapc2p_mapping_func = has_mapc2p_mapping_func;
  gks_lw->mapc2p_mapping_func_ref = (struct lua_func_ctx) {
    .func_ref = mapc2p_mapping_func_ref,
    .ndim = 2,
    .nret = 2,
    .L = L,
  };

  gks_lw->proj_id = proj_id;

  gks_lw->has_init_func = has_init_func;
  gks_lw->init_func_ref = (struct lua_func_ctx) {
    .func_ref = init_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1,
    .L = L,
  };

  gks_lw->has_density_init_func = has_density_init_func;
  gks_lw->density_init_func_ref = (struct lua_func_ctx) {
    .func_ref = density_init_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1,
    .L = L,
  };

  gks_lw->has_Upar_init_func = has_Upar_init_func;
  gks_lw->Upar_init_func_ref = (struct lua_func_ctx) {
    .func_ref = Upar_init_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1,
    .L = L,
  };

  gks_lw->has_temp_init_func = has_temp_init_func;
  gks_lw->temp_init_func_ref = (struct lua_func_ctx) {
    .func_ref = temp_init_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1,
    .L = L,
  };

  gks_lw->has_par_temp_init_func = has_par_temp_init_func;
  gks_lw->par_temp_init_func_ref = (struct lua_func_ctx) {
    .func_ref = par_temp_init_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1,
    .L = L,
  };

  gks_lw->has_perp_temp_init_func = has_perp_temp_init_func;
  gks_lw->perp_temp_init_func_ref = (struct lua_func_ctx) {
    .func_ref = perp_temp_init_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1,
    .L = L,
  };

  gks_lw->proj_correct_all_moms = proj_correct_all_moms;

  gks_lw->source_id = source_id;
  gks_lw->num_sources = num_sources;
  
  for (int i = 0; i < num_sources; i++) {
    gks_lw->source_proj_id[i] = source_proj_id[i];

    gks_lw->source_has_init_func[i] = source_has_init_func[i];
    gks_lw->source_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = source_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };

    gks_lw->source_has_density_init_func[i] = source_has_density_init_func[i];
    gks_lw->source_density_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = source_density_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };

    gks_lw->source_has_Upar_init_func[i] = source_has_Upar_init_func[i];
    gks_lw->source_Upar_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = source_Upar_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };

    gks_lw->source_has_temp_init_func[i] = source_has_temp_init_func[i];
    gks_lw->source_temp_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = source_temp_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };
  }

  gks_lw->correct_all_moms = correct_all_moms;
  gks_lw->iter_eps = iter_eps;
  gks_lw->max_iter = max_iter;
  gks_lw->use_last_converged = use_last_converged;

  gks_lw->collisionless_type = collisionless_type;

  gks_lw->collision_id = collision_id;

  gks_lw->has_self_nu_func = has_self_nu_func;
  gks_lw->self_nu_func_ref = (struct lua_func_ctx) {
    .func_ref = self_nu_func_ref,
    .ndim = 0,
    .nret = 1,
    .L = L,
  };

  gks_lw->num_cross_collisions = num_cross_collisions;
  for (int i = 0; i < num_cross_collisions; i++) {
    strcpy(gks_lw->collide_with[i], collide_with[i]);
  }

  gks_lw->collision_norm_nu = collision_norm_nu;
  gks_lw->collision_den_ref = collision_den_ref;
  gks_lw->collision_temp_ref = collision_temp_ref;
  gks_lw->collision_hbar = collision_hbar;
  gks_lw->collision_eps0 = collision_eps0;
  gks_lw->collision_eV = collision_eV;

  gks_lw->radiation_id = radiation_id;

  gks_lw->radiation_num_cross_collisions = radiation_num_cross_collisions;
  for (int i = 0; i < radiation_num_cross_collisions; i++) {
    strcpy(gks_lw->radiation_collide_with[i], radiation_collide_with[i]);

    gks_lw->radiation_z[i] = radiation_z[i];
    gks_lw->radiation_charge_state[i] = radiation_charge_state[i];
    gks_lw->radiation_num_of_densities[i] = radiation_num_of_densities[i];
  }

  gks_lw->radiation_te_min_model = radiation_te_min_model;
  gks_lw->radiation_Te_min = radiation_Te_min;

  gks_lw->num_react = num_react;

  for (int i = 0; i < num_react; i++) {
    gks_lw->react_id[i] = react_id[i];
    gks_lw->react_type_self[i] = react_type_self[i];
    gks_lw->react_ion_id[i] = react_ion_id[i];

    strcpy(gks_lw->react_elc_nm[i], react_elc_nm[i]);
    strcpy(gks_lw->react_ion_nm[i], react_ion_nm[i]);
    strcpy(gks_lw->react_donor_nm[i], react_donor_nm[i]);
    strcpy(gks_lw->react_recvr_nm[i], react_recvr_nm[i]);

    gks_lw->react_charge_state[i] = react_charge_state[i];
    gks_lw->react_ion_mass[i] = react_ion_mass[i];
    gks_lw->react_elc_mass[i] = react_elc_mass[i];
  }

  gks_lw->num_neut_react = num_neut_react;

  for (int i = 0; i < num_neut_react; i++) {
    gks_lw->neut_react_id[i] = neut_react_id[i];
    gks_lw->neut_react_type_self[i] = neut_react_type_self[i];
    gks_lw->neut_react_ion_id[i] = neut_react_ion_id[i];

    strcpy(gks_lw->neut_react_elc_nm[i], neut_react_elc_nm[i]);
    strcpy(gks_lw->neut_react_ion_nm[i], neut_react_ion_nm[i]);
    strcpy(gks_lw->neut_react_donor_nm[i], neut_react_donor_nm[i]);
    strcpy(gks_lw->neut_react_recvr_nm[i], neut_react_recvr_nm[i]);
    strcpy(gks_lw->neut_react_partner_nm[i], neut_react_partner_nm[i]);

    gks_lw->neut_react_charge_state[i] = neut_react_charge_state[i];
    gks_lw->neut_react_ion_mass[i] = neut_react_ion_mass[i];
    gks_lw->neut_react_elc_mass[i] = neut_react_elc_mass[i];
  }
  
  // Set metatable.
  luaL_getmetatable(L, GYROKINETIC_SPECIES_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor.
static struct luaL_Reg gk_species_ctor[] = {
  { "new", gyrokinetic_species_lw_new },
  { 0, 0 }
};

/* *********************** */
/* Neutral Species methods */
/* *********************** */

// Metatable name for species input struct.
#define GYROKINETIC_NEUTRAL_SPECIES_METATABLE_NM "GkeyllZero.App.Gyrokinetic.NeutralSpecies"

// Lua userdata object for constructing neutral species input.
struct gyrokinetic_neutral_species_lw {
  int magic; // This must be first element in the struct.
  
  struct gkyl_gyrokinetic_neut_species gk_neut_species; // Input struct to construct neutral species.

  enum gkyl_projection_id proj_id; // Projection type.

  bool has_density_init_func; // Is there a density initialization function?
  struct lua_func_ctx density_init_func_ref; // Lua registry reference to density initialization function.

  bool has_Udrift_init_func; // Is there a drift velocity initialiation function?
  struct lua_func_ctx Udrift_init_func_ref; // Lua registry reference to drift velocity initialization function.

  bool has_temp_init_func; // Is there a temperature initialization function?
  struct lua_func_ctx temp_init_func_ref; // Lua registry reference to temperature initialization function.
};

static int
gyrokinetic_neutral_species_lw_new(lua_State *L)
{
  struct gkyl_gyrokinetic_neut_species gk_neut_species = { };
  
  gk_neut_species.mass = glua_tbl_get_number(L, "mass", 1.0);
  gk_neut_species.is_static = glua_tbl_get_bool(L, "isStatic", false);

  with_lua_tbl_tbl(L, "cells") {
    for (int d = 0; d < 3; d++) {
      gk_neut_species.cells[d] = glua_tbl_iget_integer(L, d + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "lower") {
    for (int d = 0; d < 3; d++) {
      gk_neut_species.lower[d] = glua_tbl_iget_number(L, d + 1, 0.0);
    }
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d = 0; d < 3; d++) {
      gk_neut_species.upper[d] = glua_tbl_iget_number(L, d + 1, 0.0);
    }
  }

  with_lua_tbl_tbl(L, "diagnostics") {
    int num_diag_moments = glua_objlen(L);

    for (int i = 0; i < num_diag_moments; i ++) {
      gk_neut_species.diag_moments[i] = glua_tbl_iget_integer(L, i+1, 0);
    }

    gk_neut_species.num_diag_moments = num_diag_moments;
  }

  with_lua_tbl_tbl(L, "bcs") {
    int num_bcs = glua_objlen(L);
    for (int i = 0; i < num_bcs; i++) {
      gk_neut_species.bcs[i].dir = glua_tbl_get_integer(L, "dir", 0);
      gk_neut_species.bcs[i].edge = glua_tbl_get_integer(L, "edge", 0);
      gk_neut_species.bcs[i].type = glua_tbl_get_integer(L, "type", 0);
      with_lua_tbl_tbl(L, "value") {
        int num_vals = glua_objlen(L);
        for (int k = 0; k < num_vals; k++) {
          gk_neut_species.bcs[i].value[k] = glua_tbl_iget_number(L, k + 1, 0.0);;
        }
      }
    }
  }


  enum gkyl_projection_id proj_id = GKYL_PROJ_FUNC;

  bool has_density_init_func = false;
  int density_init_func_ref = LUA_NOREF;

  bool has_Udrift_init_func = false;
  int Udrift_init_func_ref = LUA_NOREF;

  bool has_temp_init_func = false;
  int temp_init_func_ref = LUA_NOREF;

  with_lua_tbl_tbl(L, "projection") {
    proj_id = glua_tbl_get_integer(L, "projectionID", 0);

    if (glua_tbl_get_func(L, "densityInit")) {
      density_init_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_density_init_func = true;
    }

    if (glua_tbl_get_func(L, "driftVelocityInit")) {
      Udrift_init_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_Udrift_init_func = true;
    }

    if (glua_tbl_get_func(L, "temperatureInit")) {
      temp_init_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_temp_init_func = true;
    }
  }
  
  struct gyrokinetic_neutral_species_lw *gkns_lw = lua_newuserdata(L, sizeof(*gkns_lw));
  gkns_lw->magic = GYROKINETIC_NEUTRAL_SPECIES_DEFAULT;
  gkns_lw->gk_neut_species = gk_neut_species;

  gkns_lw->proj_id = proj_id;

  gkns_lw->has_density_init_func = has_density_init_func;
  gkns_lw->density_init_func_ref = (struct lua_func_ctx) {
    .func_ref = density_init_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1,
    .L = L,
  };

  gkns_lw->has_Udrift_init_func = has_Udrift_init_func;
  gkns_lw->Udrift_init_func_ref = (struct lua_func_ctx) {
    .func_ref = Udrift_init_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 3,
    .L = L,
  };

  gkns_lw->has_temp_init_func = has_temp_init_func;
  gkns_lw->temp_init_func_ref = (struct lua_func_ctx) {
    .func_ref = temp_init_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1,
    .L = L,
  };
  
  // Set metatable.
  luaL_getmetatable(L, GYROKINETIC_NEUTRAL_SPECIES_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor.
static struct luaL_Reg gk_neutral_species_ctor[] = {
  { "new", gyrokinetic_neutral_species_lw_new },
  { 0, 0 }
};

/* ************* */
/* Field methods */
/* ************* */

// Metatable name for field input struct.
#define GYROKINETIC_FIELD_METATABLE_NM "GkeyllZero.App.Gyrokinetic.Field"

// Lua userdata object for constructing field input.
struct gyrokinetic_field_lw {
  int magic; // This must be first element in the struct.
  
  struct gkyl_gyrokinetic_field gk_field; // Input struct to construct field.
};

static int
gyrokinetic_field_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_gyrokinetic_field gk_field = { };

  gk_field.gkfield_id = glua_tbl_get_integer(L, "fieldID", 0);
  gk_field.electron_mass = glua_tbl_get_number(L, "electronMass", 0.0);
  gk_field.electron_charge = glua_tbl_get_number(L, "electronCharge", 0.0);
  gk_field.electron_density = glua_tbl_get_number(L, "electronDensity", 0.0);
  gk_field.electron_temp = glua_tbl_get_number(L, "electronTemperature", 0.0);

  gk_field.kperpSq = glua_tbl_get_number(L, "kPerpSq", 0.0);

  gk_field.zero_init_field = glua_tbl_get_bool(L, "zeroInitField", false);
  gk_field.is_static = glua_tbl_get_bool(L, "isStatic", false);

  with_lua_tbl_tbl(L, "poissonBcs") {
    int num_bcs = glua_objlen(L);
    for (int i = 0; i < num_bcs; i++) {
      gk_field.poisson_bcs[i].dir = glua_tbl_get_integer(L, "dir", 0);
      gk_field.poisson_bcs[i].edge = glua_tbl_get_integer(L, "edge", 0);
      gk_field.poisson_bcs[i].type = glua_tbl_get_integer(L, "type", 0);
      with_lua_tbl_tbl(L, "value") {
        int num_vals = glua_objlen(L);
        for (int k = 0; k < num_vals; k++) {
          gk_field.poisson_bcs[i].value[k] = glua_tbl_iget_number(L, k + 1, 0.0);;
        }
      }
    }
  }

  struct gyrokinetic_field_lw *gkf_lw = lua_newuserdata(L, sizeof(*gkf_lw));

  gkf_lw->magic = GYROKINETIC_FIELD_DEFAULT;
  gkf_lw->gk_field = gk_field;
  
  // Set metatable.
  luaL_getmetatable(L, GYROKINETIC_FIELD_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Field constructor.
static struct luaL_Reg gk_field_ctor[] = {
  { "new",  gyrokinetic_field_lw_new },
  { 0, 0 }
};

/* *********** */
/* App methods */
/* *********** */

// Metatable name for top-level Gyrokinetic App.
#define GYROKINETIC_APP_METATABLE_NM "GkeyllZero.App.Gyrokinetic"

// Lua userdata object for holding Gyrokinetic app and run parameters.
struct gyrokinetic_app_lw {
  gkyl_gyrokinetic_app *app; // Gyrokinetic app object.

  struct lua_func_ctx mapc2p_ctx; // Function context for mapc2p.
  struct lua_func_ctx bfield_ctx; // Function context for bmag.

  struct lua_func_ctx nonuniform_position_map_ctx[3]; // Function context for nonuniform position maps.

  bool has_mapc2p_mapping_func[GKYL_MAX_SPECIES]; // Is there a non-uniform velocity space mapping function?
  struct lua_func_ctx mapc2p_mapping_func_ctx[GKYL_MAX_SPECIES]; // Context for non-uniform velocity space mapping function.

  enum gkyl_projection_id proj_id[GKYL_MAX_SPECIES]; // Projection type.

  bool has_init_func[GKYL_MAX_SPECIES]; // Is there an initialization function?
  struct lua_func_ctx init_func_ctx[GKYL_MAX_SPECIES]; // Context for initialization function.

  bool has_density_init_func[GKYL_MAX_SPECIES]; // Is there a density initialization function?
  struct lua_func_ctx density_init_func_ctx[GKYL_MAX_SPECIES]; // Context for density initialization function.

  bool has_Upar_init_func[GKYL_MAX_SPECIES]; // Is there a parallel velocity initialization function?
  struct lua_func_ctx Upar_init_func_ctx[GKYL_MAX_SPECIES]; // Context for parallel velocity initialziation function.
  
  bool has_temp_init_func[GKYL_MAX_SPECIES]; // Is there a temperature initialization function?
  struct lua_func_ctx temp_init_func_ctx[GKYL_MAX_SPECIES]; // Context for temperature initialization function.

  bool has_par_temp_init_func[GKYL_MAX_SPECIES]; // Is there a parallel temperature initialization function?
  struct lua_func_ctx par_temp_init_func_ctx[GKYL_MAX_SPECIES]; // Context for parallel temperature initialization function.

  bool has_perp_temp_init_func[GKYL_MAX_SPECIES]; // Is there a perpendicular temperature initialization function?
  struct lua_func_ctx perp_temp_init_func_ctx[GKYL_MAX_SPECIES]; // Context for perpendicular temperature initialization function.

  bool proj_correct_all_moms[GKYL_MAX_SPECIES]; // Are we correcting all moments in projection, or only density?

  enum gkyl_gk_collisionless_type collisionless_type[GKYL_MAX_SPECIES]; // Collisionless type.
                                                        
  enum gkyl_collision_id collision_id[GKYL_MAX_SPECIES]; // Collision type.

  bool has_self_nu_func[GKYL_MAX_SPECIES]; // Is there a self-collision frequency function?
  struct lua_func_ctx self_nu_func_ctx[GKYL_MAX_SPECIES]; // Context for self-collision frequency function.

  int num_cross_collisions[GKYL_MAX_SPECIES]; // Number of species that we cross-collide with.
  char collide_with[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES][128]; // Names of species that we cross-collide with.

  bool collision_norm_nu[GKYL_MAX_SPECIES]; // Are we rescaling the collision frequency?
  double collision_den_ref[GKYL_MAX_SPECIES]; // Density used to calculate Coulomb logarithm for collision frequency.
  double collision_temp_ref[GKYL_MAX_SPECIES]; // Temperature used to calculate Coulomb logarithm for collision frequency.
  double collision_hbar[GKYL_MAX_SPECIES]; // Reduced Planck's constant for calculating collision frequency.
  double collision_eps0[GKYL_MAX_SPECIES]; // Vacuum permittivity for calculating collision frequency.
  double collision_eV[GKYL_MAX_SPECIES]; // Elementary charge for calculating collision frequency.

  bool correct_all_moms[GKYL_MAX_SPECIES]; // Are we correcting all moments in collisions, or only density?
  double iter_eps[GKYL_MAX_SPECIES]; // Error tolerance for moment fixes in collision (density is always exact).
  int max_iter[GKYL_MAX_SPECIES]; // Maximum number of iterations for moment fixes in collisions.
  bool use_last_converged[GKYL_MAX_SPECIES]; // Use last iteration value in collisions regardless of convergence?

  enum gkyl_source_id source_id[GKYL_MAX_SPECIES]; // Source type.

  int num_sources[GKYL_MAX_SPECIES]; // Number of projection objects in source.
  enum gkyl_projection_id source_proj_id[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Projection type in source.

  bool source_has_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there an initialization function in source?
  struct lua_func_ctx source_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for initialization function in source.

  bool source_has_density_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there a density initialization function in source?
  struct lua_func_ctx source_density_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for density initialization function in source.

  bool source_has_Upar_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there a parallel velocity initialization function in source?
  struct lua_func_ctx source_Upar_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for parallel velocity initialization function in source.

  bool source_has_temp_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there a temperature initialization function in source?
  struct lua_func_ctx source_temp_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for temperature initialization function in source.

  enum gkyl_radiation_id radiation_id[GKYL_MAX_SPECIES]; // Radiation type.

  int radiation_num_cross_collisions[GKYL_MAX_SPECIES]; // Number of radiation species that we cross-collide with.
  char radiation_collide_with[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES][128]; // Names of radiation species that we cross-collide with.

  int radiation_z[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES]; // Atomic Z of radiation species that we are colliding with.
  int radiation_charge_state[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES]; // Charge state of radiation species that we are colliding with.
  int radiation_num_of_densities[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES]; // Maximum number of densities to use per charge state of radiation species that we are colliding with.
  
  enum gkyl_te_min_model radiation_te_min_model[GKYL_MAX_SPECIES]; // How is the radiation turned off (constant, or with varying electron temperature)?
  double radiation_Te_min[GKYL_MAX_SPECIES]; // Minimum temperature (in J) at which to stop radiating.

  int num_react[GKYL_MAX_SPECIES]; // Number of reaction types.

  enum gkyl_react_id react_id[GKYL_MAX_SPECIES][GKYL_MAX_REACT]; // What type of reaction (ionization, charge exchange, recombination)?
  enum gkyl_react_self_type react_type_self[GKYL_MAX_SPECIES][GKYL_MAX_REACT]; // What is the role of the species in this reaction?
  enum gkyl_ion_type react_ion_id[GKYL_MAX_SPECIES][GKYL_MAX_REACT]; // What type of ion is reacting?

  char react_elc_nm[GKYL_MAX_SPECIES][GKYL_MAX_REACT][128]; // Name of electron species in the reaction.
  char react_ion_nm[GKYL_MAX_SPECIES][GKYL_MAX_REACT][128]; // Name of ion species in the reaction.
  char react_donor_nm[GKYL_MAX_SPECIES][GKYL_MAX_REACT][128]; // Name of donor species in the reaction.
  char react_recvr_nm[GKYL_MAX_SPECIES][GKYL_MAX_REACT][128]; // Name of receiver species in the reaction.

  int react_charge_state[GKYL_MAX_SPECIES][GKYL_MAX_REACT]; // Charge state of species in the reaction.
  double react_ion_mass[GKYL_MAX_SPECIES][GKYL_MAX_REACT]; // Mass of ion species in the reaction.
  double react_elc_mass[GKYL_MAX_SPECIES][GKYL_MAX_REACT]; // Mass of electron species in the reaction.

  int num_neut_react[GKYL_MAX_SPECIES]; // Number of neutral reaction types.

  enum gkyl_react_id neut_react_id[GKYL_MAX_SPECIES][GKYL_MAX_REACT]; // What type of neutral reaction (ionization, charge exchange, recombination)?
  enum gkyl_react_self_type neut_react_type_self[GKYL_MAX_SPECIES][GKYL_MAX_REACT]; // What is the role of the species in this neutral reaction?
  enum gkyl_ion_type neut_react_ion_id[GKYL_MAX_SPECIES][GKYL_MAX_REACT]; // What type of ion in the neutral reaction?

  char neut_react_elc_nm[GKYL_MAX_SPECIES][GKYL_MAX_REACT][128]; // Name of electron species in the neutral reaction.
  char neut_react_ion_nm[GKYL_MAX_SPECIES][GKYL_MAX_REACT][128]; // Name of ion species in the neutral reaction.
  char neut_react_donor_nm[GKYL_MAX_SPECIES][GKYL_MAX_REACT][128]; // Name of donor species in the neutral reaction.
  char neut_react_recvr_nm[GKYL_MAX_SPECIES][GKYL_MAX_REACT][128]; // Name of receiver species in the neutral reaction.
  char neut_react_partner_nm[GKYL_MAX_SPECIES][GKYL_MAX_REACT][128]; // Name of partner species in the neutral reaction.,

  int neut_react_charge_state[GKYL_MAX_SPECIES][GKYL_MAX_REACT]; // Charge state of species in the neutral reaction.
  double neut_react_ion_mass[GKYL_MAX_SPECIES][GKYL_MAX_REACT]; // Mass of ion species in the neutral reaction.
  double neut_react_elc_mass[GKYL_MAX_SPECIES][GKYL_MAX_REACT]; // Mass of electron species in the neutral reaction.

  enum gkyl_projection_id neut_proj_id[GKYL_MAX_SPECIES]; // Neutral projection type.

  bool neut_has_density_init_func[GKYL_MAX_SPECIES]; // Is there a neutral density initialization function?
  struct lua_func_ctx neut_density_init_func_ctx[GKYL_MAX_SPECIES]; // Context for neutral density initialization function.

  bool neut_has_Udrift_init_func[GKYL_MAX_SPECIES]; // Is there a neutral drift velocity initialization function?
  struct lua_func_ctx neut_Udrift_init_func_ctx[GKYL_MAX_SPECIES]; // Context for neutral drift velocity initialziation function.
  
  bool neut_has_temp_init_func[GKYL_MAX_SPECIES]; // Is there a neutral temperature initialization function?
  struct lua_func_ctx neut_temp_init_func_ctx[GKYL_MAX_SPECIES]; // Context for neutral temperature initialization function.
  
  double t_start, t_end; // Start and end times of simulation.
  int num_frames; // Number of data frames to write.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

// Gets all species objects from the App table, which must on top of
// the stack. The number of species is returned and the appropriate
// pointers set in the species pointer array.
static int
get_species_inp(lua_State *L, int cdim, struct gyrokinetic_species_lw *species[GKYL_MAX_SPECIES])
{
  enum { TKEY = -2, TVAL = -1 };
  
  int curr = 0;
  lua_pushnil(L); // Initial key is nil.
  while (lua_next(L, TKEY) != 0) {
    // Key at TKEY and value at TVAL.
    if (lua_type(L, TVAL) == LUA_TUSERDATA) {
      struct gyrokinetic_species_lw *gks = lua_touserdata(L, TVAL);

      if (gks->magic == GYROKINETIC_SPECIES_DEFAULT) {
        if (gks->has_init_func) {
          gks->init_func_ref.ndim = cdim + gks->vdim;
        }

        if (gks->has_density_init_func) {
          gks->density_init_func_ref.ndim = cdim;
        }
        
        if (gks->has_Upar_init_func) {
          gks->Upar_init_func_ref.ndim = cdim;
        }

        if (gks->has_temp_init_func) {
          gks->temp_init_func_ref.ndim = cdim;
        }

        if (gks->has_par_temp_init_func) {
          gks->par_temp_init_func_ref.ndim = cdim;
        }

        if (gks->has_perp_temp_init_func) {
          gks->perp_temp_init_func_ref.ndim = cdim;
        }

        if (gks->has_self_nu_func) {
          gks->self_nu_func_ref.ndim = cdim;
        }

        for (int i = 0; i < gks->num_sources; i++) {
          if (gks->source_has_init_func[i]) {
            gks->source_init_func_ref[i].ndim = cdim + gks->vdim;
          }

          if (gks->source_has_density_init_func[i]) {
            gks->source_density_init_func_ref[i].ndim = cdim;
          }

          if (gks->source_has_Upar_init_func[i]) {
            gks->source_Upar_init_func_ref[i].ndim = cdim;
          }

          if (gks->source_has_temp_init_func[i]) {
            gks->source_temp_init_func_ref[i].ndim = cdim;
          }
        }
        
        if (lua_type(L,TKEY) == LUA_TSTRING) {
          const char *key = lua_tolstring(L, TKEY, 0);
          strcpy(gks->gk_species.name, key);
        }
        species[curr++] = gks;
      }
    }
    lua_pop(L, 1);
  }

  return curr;
}

// Comparison method to sort species array by species name.
static int
species_compare_func(const void *a, const void *b)
{
  const struct gyrokinetic_species_lw *const *spa = a;
  const struct gyrokinetic_species_lw *const *spb = b;
  return strcmp((*spa)->gk_species.name, (*spb)->gk_species.name);
}

// Gets all neutral species objects from the App table, which must on top of
// the stack. The number of neutral species is returned and the appropriate
// pointers set in the neutral species pointer array.
static int
get_neutral_species_inp(lua_State *L, int cdim, struct gyrokinetic_neutral_species_lw *neut_species[GKYL_MAX_SPECIES])
{
  enum { TKEY = -2, TVAL = -1 };
  
  int curr = 0;
  lua_pushnil(L); // Initial key is nil.
  while (lua_next(L, TKEY) != 0) {
    // Key at TKEY and value at TVAL.
    if (lua_type(L, TVAL) == LUA_TUSERDATA) {
      struct gyrokinetic_neutral_species_lw *gkns = lua_touserdata(L, TVAL);

      if (gkns->magic == GYROKINETIC_NEUTRAL_SPECIES_DEFAULT) {
        if (gkns->has_density_init_func) {
          gkns->density_init_func_ref.ndim = cdim;
        }
        
        if (gkns->has_Udrift_init_func) {
          gkns->Udrift_init_func_ref.ndim = cdim;
        }

        if (gkns->has_temp_init_func) {
          gkns->temp_init_func_ref.ndim = cdim;
        }
        
        if (lua_type(L,TKEY) == LUA_TSTRING) {
          const char *key = lua_tolstring(L, TKEY, 0);
          strcpy(gkns->gk_neut_species.name, key);
        }
        neut_species[curr++] = gkns;
      }
    }
    lua_pop(L, 1);
  }

  return curr;
}

// Comparison method to sort neutral species array by neutral species name.
static int
neutral_species_compare_func(const void *a, const void *b)
{
  const struct gyrokinetic_neutral_species_lw *const *spa = a;
  const struct gyrokinetic_neutral_species_lw *const *spb = b;
  return strcmp((*spa)->gk_neut_species.name, (*spb)->gk_neut_species.name);
}

static struct gkyl_tool_args *
tool_args_from_argv(int optind, int argc, char *const*argv)
{
  struct gkyl_tool_args *targs = gkyl_malloc(sizeof *targs);
  
  targs->argc = argc-optind;
  targs->argv = 0;

  if (targs->argc > 0) {
    targs->argv = gkyl_malloc(targs->argc*sizeof(char *));
      for (int i = optind, j = 0; i < argc; ++i, ++j) {
        targs->argv[j] = gkyl_malloc(strlen(argv[i])+1);
        strcpy(targs->argv[j], argv[i]);
      }
  }

  return targs;
}

// CLI parser for main script.
struct script_cli {
  bool help; // Show help.
  bool step_mode; // Run for fixed number of steps? (for valgrind/cuda-memcheck)
  int num_steps; // Number of steps.
  bool use_mpi; // Should we use MPI?
  bool use_gpu; // Should this be run on GPU?
  bool trace_mem; // Should we trace memory allocation/deallocation?
  bool use_verbose; // Should we use verbose output?
  bool is_restart; // Is this a restarted simulation?
  int restart_frame; // Which frame to restart simulation from.  
  
  struct gkyl_tool_args *rest;
};

static struct script_cli
gk_parse_script_cli(struct gkyl_tool_args *acv)
{
  struct script_cli cli = {
    .help =- false,
    .step_mode = false,
    .num_steps = INT_MAX,
    .use_mpi = false,
    .use_gpu = false,
    .trace_mem = false,
    .use_verbose = false,
    .is_restart = false,
    .restart_frame = 0,
  };

#ifdef GKYL_HAVE_MPI
  cli.use_mpi = true;
#endif
#ifdef GKYL_HAVE_CUDA
  cli.use_gpu = true;
#endif
  
  coption_long longopts[] = {
    { 0 }
  };
  const char* shortopts = "+hVs:SGmr:";

  coption opt = coption_init();
  int c;
  while ((c = coption_get(&opt, acv->argc, acv->argv, shortopts, longopts)) != -1) {
    switch (c) {
      case 'h':
        cli.help = true;
        break;

      case 's':
        cli.num_steps = atoi(opt.arg);
        break;
      
      case 'S':
        cli.use_mpi = false;
        break;
      
      case 'G':
        cli.use_gpu = false;
        break;
      
      case 'm':
        cli.trace_mem = true;
        break;
      
      case 'V':
        cli.use_verbose = true;
        break;
      
      case 'r':
        cli.is_restart = true;
        cli.restart_frame = atoi(opt.arg);
        break;        
        
      case '?':
        break;
    }
  }

  cli.rest = tool_args_from_argv(opt.ind, acv->argc, acv->argv);
  
  return cli;
}

// Create top-level App object.
static int
gk_app_new(lua_State *L)
{
  struct gyrokinetic_app_lw *app_lw = gkyl_malloc(sizeof(*app_lw));

  // The output prefix to use is stored in the global
  // GKYL_OUT_PREFIX. If this is not found then "g0-gyrokinetic" is used.
  const char *sim_name = "g0-gyrokinetic";

  with_lua_global(L, "GKYL_OUT_PREFIX") {
    if (lua_isstring(L, -1)) {
      sim_name = lua_tostring(L, -1);
    }
  }
  
  // Initialize app using table inputs (table is on top of stack).

  app_lw->t_start = glua_tbl_get_number(L, "tStart", 0.0);
  app_lw->t_end = glua_tbl_get_number(L, "tEnd", 1.0);
  app_lw->num_frames = glua_tbl_get_integer(L, "nFrame", 1);
  app_lw->field_energy_calcs = glua_tbl_get_integer(L, "fieldEnergyCalcs", INT_MAX);
  app_lw->integrated_mom_calcs = glua_tbl_get_integer(L, "integratedMomentCalcs", INT_MAX);
  app_lw->dt_failure_tol = glua_tbl_get_number(L, "dtFailureTol", 1.0e-4);
  app_lw->num_failures_max = glua_tbl_get_integer(L, "numFailuresMax", 20);

  struct gkyl_gk gk = { }; // Input table for app.

  strcpy(gk.name, sim_name);
  
  int cdim = 0;
  with_lua_tbl_tbl(L, "cells") {
    gk.cdim = cdim = glua_objlen(L);

    for (int d = 0; d < cdim; d++) {
      gk.cells[d] = glua_tbl_iget_integer(L, d + 1, 0);
    }
  }

  int cuts[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; d++) {
    cuts[d] = 1;
  }
  
  with_lua_tbl_tbl(L, "decompCuts") {
    int ncuts = glua_objlen(L);

    for (int d = 0; d < ncuts; d++) {
      cuts[d] = glua_tbl_iget_integer(L, d + 1, 0);
    }
  }  

  with_lua_tbl_tbl(L, "lower") {
    for (int d = 0; d < cdim; d++) {
      gk.lower[d] = glua_tbl_iget_number(L, d + 1, 0.0);
    }
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d = 0; d < cdim; d++) {
      gk.upper[d] = glua_tbl_iget_number(L, d + 1, 0.0);
    }
  }

  gk.cfl_frac = glua_tbl_get_number(L, "cflFrac", 0.95);
  gk.poly_order = glua_tbl_get_integer(L, "polyOrder", 1);

  gk.basis_type = get_basis_type(
    glua_tbl_get_string(L, "basis", "serendipity")
  );

  gk.num_periodic_dir = 0;
  if (glua_tbl_has_key(L, "periodicDirs")) {
    with_lua_tbl_tbl(L, "periodicDirs") {
      gk.num_periodic_dir = glua_objlen(L);

      for (int d = 0; d < gk.num_periodic_dir; d++) {
        // Indices are off by 1 between Lua and C.
        gk.periodic_dirs[d] = glua_tbl_iget_integer(L, d + 1, 0) - 1;
      }
    }
  }

  with_lua_tbl_tbl(L, "geometry") {
    gk.geometry.geometry_id = glua_tbl_get_integer(L, "geometryID", 0);

    with_lua_tbl_tbl(L, "world") {
      for (int d = 0; d < 3 - cdim; d++) {
        gk.geometry.world[d] = glua_tbl_iget_number(L, d + 1, 0.0);
      }
    }

    gk.geometry.c2p_ctx = 0;
    gk.geometry.mapc2p = 0;
    bool has_mapc2p = false;
    int mapc2p_ref = LUA_NOREF;
    if (glua_tbl_get_func(L, "mapc2p")) {
      has_mapc2p = true;
      mapc2p_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    }

    gk.geometry.bfield_ctx = 0;
    gk.geometry.bfield_func = 0;
    bool has_bfield_func = false;
    int bfield_func_ref = LUA_NOREF;
    if (glua_tbl_get_func(L, "bfieldFunc")) {
      has_bfield_func = true;
      bfield_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    }

    with_lua_tbl_tbl(L, "positionMap") {
      gk.geometry.position_map_info.id = glua_tbl_get_integer(L, "ID", 0);
      bool has_nonuniform_position_map[3];
      int nonuniform_position_map_ref[3];

      with_lua_tbl_tbl(L, "maps") {
        for (int i = 0; i < 3; i++) {
          gk.geometry.position_map_info.ctxs[i] = 0;
          gk.geometry.position_map_info.maps[i] = 0;

          has_nonuniform_position_map[i] = false;
          nonuniform_position_map_ref[i] = LUA_NOREF;
          if (glua_tbl_iget_func(L, i + 1)) {
            has_nonuniform_position_map[i] = true;
            nonuniform_position_map_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        if (has_nonuniform_position_map[i]) {
          app_lw->nonuniform_position_map_ctx[i] = (struct lua_func_ctx) {
            .func_ref = nonuniform_position_map_ref[i],
            .ndim = 1,
            .nret = 1,
            .L = L,
          };
          gk.geometry.position_map_info.maps[i] = gkyl_lw_eval_cb;
          gk.geometry.position_map_info.ctxs[i] = &app_lw->nonuniform_position_map_ctx[i];
        }
      }
    }

    if (has_mapc2p) {
      app_lw->mapc2p_ctx = (struct lua_func_ctx) {
        .func_ref = mapc2p_ref,
        .ndim = 3,
        .nret = 3,
        .L = L,
      };
      gk.geometry.mapc2p = gkyl_lw_eval_cb;
      gk.geometry.c2p_ctx = &app_lw->mapc2p_ctx;
    }

    if (has_bfield_func) {
      app_lw->bfield_ctx = (struct lua_func_ctx) {
        .func_ref = bfield_func_ref,
        .ndim = 3,
        .nret = 3,
        .L = L,
      };
      gk.geometry.bfield_func = gkyl_lw_eval_cb;
      gk.geometry.bfield_ctx = &app_lw->bfield_ctx;
    }
  }

  struct gyrokinetic_species_lw *species[GKYL_MAX_SPECIES];

  // Set all species input.
  gk.num_species = get_species_inp(L, cdim, species);

  // need to sort the species[] array by name of the species before
  // proceeding as there is no way to ensure that all cores loop over
  // Lua tables in the same order
  qsort(species, gk.num_species, sizeof(struct gyrokinetic_species_lw *), species_compare_func);
  
  for (int s = 0; s < gk.num_species; s++) {
    gk.species[s] = species[s]->gk_species;
    gk.species[s].vdim = species[s]->vdim;

    app_lw->has_mapc2p_mapping_func[s] = species[s]->has_mapc2p_mapping_func;
    app_lw->mapc2p_mapping_func_ctx[s] = species[s]->mapc2p_mapping_func_ref;

    app_lw->proj_id[s] = species[s]->proj_id;

    app_lw->has_init_func[s] = species[s]->has_init_func;
    app_lw->init_func_ctx[s] = species[s]->init_func_ref;

    app_lw->has_density_init_func[s] = species[s]->has_density_init_func;
    app_lw->density_init_func_ctx[s] = species[s]->density_init_func_ref;

    app_lw->has_Upar_init_func[s] = species[s]->has_Upar_init_func;
    app_lw->Upar_init_func_ctx[s] = species[s]->Upar_init_func_ref;
    
    app_lw->has_temp_init_func[s] = species[s]->has_temp_init_func;
    app_lw->temp_init_func_ctx[s] = species[s]->temp_init_func_ref;

    app_lw->has_par_temp_init_func[s] = species[s]->has_par_temp_init_func;
    app_lw->par_temp_init_func_ctx[s] = species[s]->par_temp_init_func_ref;
    
    app_lw->has_perp_temp_init_func[s] = species[s]->has_perp_temp_init_func;
    app_lw->perp_temp_init_func_ctx[s] = species[s]->perp_temp_init_func_ref;

    app_lw->proj_correct_all_moms[s] = species[s]->proj_correct_all_moms;

    if (species[s]->has_mapc2p_mapping_func) {
      gk.species[s].mapc2p.mapping = gkyl_lw_eval_cb;
      gk.species[s].mapc2p.ctx = &app_lw->mapc2p_mapping_func_ctx[s];
    }

    gk.species[s].projection.proj_id = app_lw->proj_id[s];

    if (species[s]->has_init_func) {
      gk.species[s].projection.func = gkyl_lw_eval_cb;
      gk.species[s].projection.ctx_func = &app_lw->init_func_ctx[s];
    }

    if (species[s]->has_density_init_func) {
      gk.species[s].projection.density = gkyl_lw_eval_cb;
      gk.species[s].projection.ctx_density = &app_lw->density_init_func_ctx[s];
    }

    if (species[s]->has_Upar_init_func) {
      gk.species[s].projection.upar = gkyl_lw_eval_cb;
      gk.species[s].projection.ctx_upar = &app_lw->Upar_init_func_ctx[s];
    }

    if (species[s]->has_temp_init_func) {
      gk.species[s].projection.temp = gkyl_lw_eval_cb;
      gk.species[s].projection.ctx_temp = &app_lw->temp_init_func_ctx[s];
    }

    if (species[s]->has_par_temp_init_func) {
      gk.species[s].projection.temppar = gkyl_lw_eval_cb;
      gk.species[s].projection.ctx_temppar = &app_lw->par_temp_init_func_ctx[s];
    }

    if (species[s]->has_perp_temp_init_func) {
      gk.species[s].projection.tempperp = gkyl_lw_eval_cb;
      gk.species[s].projection.ctx_tempperp = &app_lw->perp_temp_init_func_ctx[s];
    }

    gk.species[s].projection.correct_all_moms = app_lw->proj_correct_all_moms[s];

    app_lw->correct_all_moms[s] = species[s]->correct_all_moms;
    app_lw->iter_eps[s] = species[s]->iter_eps;
    app_lw->max_iter[s] = species[s]->max_iter;
    app_lw->use_last_converged[s] = species[s]->use_last_converged;

    gk.species[s].correct.correct_all_moms = app_lw->correct_all_moms[s];
    gk.species[s].correct.iter_eps = app_lw->iter_eps[s];
    gk.species[s].correct.max_iter = app_lw->max_iter[s];
    gk.species[s].correct.use_last_converged = app_lw->use_last_converged[s];

    app_lw->collisionless_type[s] = species[s]->collisionless_type;

    app_lw->collision_id[s] = species[s]->collision_id;

    app_lw->has_self_nu_func[s] = species[s]->has_self_nu_func;
    app_lw->self_nu_func_ctx[s] = species[s]->self_nu_func_ref;

    app_lw->num_cross_collisions[s] = species[s]->num_cross_collisions;
    for (int i = 0; i < app_lw->num_cross_collisions[s]; i++) {
      strcpy(app_lw->collide_with[s][i], species[s]->collide_with[i]);
    }

    app_lw->collision_norm_nu[s] = species[s]->collision_norm_nu;
    app_lw->collision_den_ref[s] = species[s]->collision_den_ref;
    app_lw->collision_temp_ref[s] = species[s]->collision_temp_ref;
    app_lw->collision_hbar[s] = species[s]->collision_hbar;
    app_lw->collision_eps0[s] = species[s]->collision_eps0;
    app_lw->collision_eV[s] = species[s]->collision_eV;

    gk.species[s].collisionless.type = app_lw->collisionless_type[s];

    gk.species[s].collisions.collision_id = app_lw->collision_id[s];

    if (species[s]->has_self_nu_func) {
      gk.species[s].collisions.self_nu = gkyl_lw_eval_cb;
      gk.species[s].collisions.self_nu_ctx = &app_lw->self_nu_func_ctx[s];
    }

    gk.species[s].collisions.num_cross_collisions = app_lw->num_cross_collisions[s];
    for (int i = 0; i < app_lw->num_cross_collisions[s]; i++) {
      strcpy(gk.species[s].collisions.collide_with[i], app_lw->collide_with[s][i]);
    }

    gk.species[s].collisions.den_ref = app_lw->collision_den_ref[s];
    gk.species[s].collisions.temp_ref = app_lw->collision_temp_ref[s];
    gk.species[s].collisions.hbar = app_lw->collision_hbar[s];
    gk.species[s].collisions.eps0 = app_lw->collision_eps0[s];
    gk.species[s].collisions.eV = app_lw->collision_eV[s];

    app_lw->source_id[s] = species[s]->source_id;

    app_lw->num_sources[s] = species[s]->num_sources;
    for (int i = 0; i < app_lw->num_sources[s]; i++) {
      app_lw->source_proj_id[s][i] = species[s]->source_proj_id[i];

      app_lw->source_has_init_func[s][i] = species[s]->source_has_init_func[i];
      app_lw->source_init_func_ctx[s][i] = species[s]->source_init_func_ref[i];

      app_lw->source_has_density_init_func[s][i] = species[s]->source_has_density_init_func[i];
      app_lw->source_density_init_func_ctx[s][i] = species[s]->source_density_init_func_ref[i];

      app_lw->source_has_Upar_init_func[s][i] = species[s]->source_has_Upar_init_func[i];
      app_lw->source_Upar_init_func_ctx[s][i] = species[s]->source_Upar_init_func_ref[i];

      app_lw->source_has_temp_init_func[s][i] = species[s]->source_has_temp_init_func[i];
      app_lw->source_temp_init_func_ctx[s][i] = species[s]->source_temp_init_func_ref[i];
    }

    gk.species[s].source.source_id = app_lw->source_id[s];

    gk.species[s].source.num_sources = app_lw->num_sources[s];
    for (int i = 0; i < app_lw->num_sources[s]; i++) {
      gk.species[s].source.projection[i].proj_id = app_lw->source_proj_id[s][i];

      if (species[s]->source_has_init_func[i]) {
        gk.species[s].source.projection[i].func = gkyl_lw_eval_cb;
        gk.species[s].source.projection[i].ctx_func = &app_lw->source_init_func_ctx[s][i];
      }

      if (species[s]->source_has_density_init_func[i]) {
        gk.species[s].source.projection[i].density = gkyl_lw_eval_cb;
        gk.species[s].source.projection[i].ctx_density = &app_lw->source_density_init_func_ctx[s][i];
      }

      if (species[s]->source_has_Upar_init_func[i]) {
        gk.species[s].source.projection[i].upar = gkyl_lw_eval_cb;
        gk.species[s].source.projection[i].ctx_upar = &app_lw->source_Upar_init_func_ctx[s][i];
      }

      if (species[s]->source_has_temp_init_func[i]) {
        gk.species[s].source.projection[i].temp = gkyl_lw_eval_cb;
        gk.species[s].source.projection[i].ctx_temp = &app_lw->source_temp_init_func_ctx[s][i];
      }
    }

    app_lw->radiation_id[s] = species[s]->radiation_id;

    app_lw->radiation_num_cross_collisions[s] = species[s]->radiation_num_cross_collisions;
    for (int i = 0; i < app_lw->radiation_num_cross_collisions[s]; i++) {
      strcpy(app_lw->radiation_collide_with[s][i], species[s]->radiation_collide_with[i]);

      app_lw->radiation_z[s][i] = species[s]->radiation_z[i];
      app_lw->radiation_charge_state[s][i] = species[s]->radiation_charge_state[i];
      app_lw->radiation_num_of_densities[s][i] = species[s]->radiation_num_of_densities[i];
    }

    app_lw->radiation_te_min_model[s] = species[s]->radiation_te_min_model;
    app_lw->radiation_Te_min[s] = species[s]->radiation_Te_min;

    gk.species[s].radiation.radiation_id = app_lw->radiation_id[s];

    gk.species[s].radiation.num_cross_collisions = app_lw->radiation_num_cross_collisions[s];
    for (int i = 0; i < app_lw->radiation_num_cross_collisions[s]; i++) {
      strcpy(gk.species[s].radiation.collide_with[i], app_lw->radiation_collide_with[s][i]);

      gk.species[s].radiation.atomic_Z[i] = app_lw->radiation_z[s][i];
      gk.species[s].radiation.charge_state[i] = app_lw->radiation_charge_state[s][i];
      gk.species[s].radiation.num_of_densities[i] = app_lw->radiation_num_of_densities[s][i];
    }

    gk.species[s].radiation.te_min_model = app_lw->radiation_te_min_model[s];
    gk.species[s].radiation.Te_min = app_lw->radiation_Te_min[s];

    app_lw->num_react[s] = species[s]->num_react;

    for (int i = 0; i < app_lw->num_react[s]; i++) {
      app_lw->react_id[s][i] = species[s]->react_id[i];
      app_lw->react_type_self[s][i] = species[s]->react_type_self[i];
      app_lw->react_ion_id[s][i] = species[s]->react_ion_id[i];

      strcpy(app_lw->react_elc_nm[s][i], species[s]->react_elc_nm[i]);
      strcpy(app_lw->react_ion_nm[s][i], species[s]->react_ion_nm[i]);
      strcpy(app_lw->react_donor_nm[s][i], species[s]->react_donor_nm[i]);
      strcpy(app_lw->react_recvr_nm[s][i], species[s]->react_recvr_nm[i]);

      app_lw->react_charge_state[s][i] = species[s]->react_charge_state[i];
      app_lw->react_ion_mass[s][i] = species[s]->react_ion_mass[i];
      app_lw->react_elc_mass[s][i] = species[s]->react_elc_mass[i];
    }

    gk.species[s].react.num_react = app_lw->num_react[s];

    for (int i = 0; i < app_lw->num_react[s]; i++) {
      gk.species[s].react.react_type[i].react_id = app_lw->react_id[s][i];
      gk.species[s].react.react_type[i].type_self = app_lw->react_type_self[s][i];
      gk.species[s].react.react_type[i].ion_id = app_lw->react_ion_id[s][i];

      strcpy(gk.species[s].react.react_type[i].elc_nm, app_lw->react_elc_nm[s][i]);
      strcpy(gk.species[s].react.react_type[i].ion_nm, app_lw->react_ion_nm[s][i]);
      strcpy(gk.species[s].react.react_type[i].recvr_nm, app_lw->react_recvr_nm[s][i]);
      strcpy(gk.species[s].react.react_type[i].donor_nm, app_lw->react_donor_nm[s][i]);

      gk.species[s].react.react_type[i].charge_state = app_lw->react_charge_state[s][i];
      gk.species[s].react.react_type[i].ion_mass = app_lw->react_ion_mass[s][i];
      gk.species[s].react.react_type[i].elc_mass = app_lw->react_elc_mass[s][i];
    }

    app_lw->num_neut_react[s] = species[s]->num_neut_react;

    for (int i = 0; i < app_lw->num_neut_react[s]; i++) {
      app_lw->neut_react_id[s][i] = species[s]->neut_react_id[i];
      app_lw->neut_react_type_self[s][i] = species[s]->neut_react_type_self[i];
      app_lw->neut_react_ion_id[s][i] = species[s]->neut_react_ion_id[i];

      strcpy(app_lw->neut_react_elc_nm[s][i], species[s]->neut_react_elc_nm[i]);
      strcpy(app_lw->neut_react_ion_nm[s][i], species[s]->neut_react_ion_nm[i]);
      strcpy(app_lw->neut_react_donor_nm[s][i], species[s]->neut_react_donor_nm[i]);
      strcpy(app_lw->neut_react_recvr_nm[s][i], species[s]->neut_react_recvr_nm[i]);
      strcpy(app_lw->neut_react_partner_nm[s][i], species[s]->neut_react_partner_nm[i]);

      app_lw->neut_react_charge_state[s][i] = species[s]->neut_react_charge_state[i];
      app_lw->neut_react_ion_mass[s][i] = species[s]->neut_react_ion_mass[i];
      app_lw->neut_react_elc_mass[s][i] = species[s]->neut_react_elc_mass[i];
    }

    gk.species[s].react_neut.num_react = app_lw->num_neut_react[s];

    for (int i = 0; i < app_lw->num_neut_react[s]; i++) {
      gk.species[s].react_neut.react_type[i].react_id = app_lw->neut_react_id[s][i];
      gk.species[s].react_neut.react_type[i].type_self = app_lw->neut_react_type_self[s][i];
      gk.species[s].react_neut.react_type[i].ion_id = app_lw->neut_react_ion_id[s][i];

      strcpy(gk.species[s].react_neut.react_type[i].elc_nm, app_lw->neut_react_elc_nm[s][i]);
      strcpy(gk.species[s].react_neut.react_type[i].ion_nm, app_lw->neut_react_ion_nm[s][i]);
      strcpy(gk.species[s].react_neut.react_type[i].recvr_nm, app_lw->neut_react_recvr_nm[s][i]);
      strcpy(gk.species[s].react_neut.react_type[i].donor_nm, app_lw->neut_react_donor_nm[s][i]);
      strcpy(gk.species[s].react_neut.react_type[i].partner_nm, app_lw->neut_react_partner_nm[s][i]);

      gk.species[s].react_neut.react_type[i].charge_state = app_lw->neut_react_charge_state[s][i];
      gk.species[s].react_neut.react_type[i].ion_mass = app_lw->neut_react_ion_mass[s][i];
      gk.species[s].react_neut.react_type[i].elc_mass = app_lw->neut_react_elc_mass[s][i];
    }
  }

  struct gyrokinetic_neutral_species_lw *neut_species[GKYL_MAX_SPECIES];

  // Set all neutral species input.
  gk.num_neut_species = get_neutral_species_inp(L, cdim, neut_species);

  // need to sort the neut_species[] array by name of the neutral species before
  // proceeding as there is no way to ensure that all cores loop over
  // Lua tables in the same order
  qsort(neut_species, gk.num_neut_species, sizeof(struct gyrokinetic_neutral_species_lw *), neutral_species_compare_func);
  
  for (int s = 0; s < gk.num_neut_species; s++) {
    gk.neut_species[s] = neut_species[s]->gk_neut_species;

    app_lw->neut_proj_id[s] = neut_species[s]->proj_id;

    app_lw->neut_has_density_init_func[s] = neut_species[s]->has_density_init_func;
    app_lw->neut_density_init_func_ctx[s] = neut_species[s]->density_init_func_ref;

    app_lw->neut_has_Udrift_init_func[s] = neut_species[s]->has_Udrift_init_func;
    app_lw->neut_Udrift_init_func_ctx[s] = neut_species[s]->Udrift_init_func_ref;
    
    app_lw->neut_has_temp_init_func[s] = neut_species[s]->has_temp_init_func;
    app_lw->neut_temp_init_func_ctx[s] = neut_species[s]->temp_init_func_ref;

    gk.neut_species[s].projection.proj_id = app_lw->neut_proj_id[s];

    if (neut_species[s]->has_density_init_func) {
      gk.neut_species[s].projection.density = gkyl_lw_eval_cb;
      gk.neut_species[s].projection.ctx_density = &app_lw->neut_density_init_func_ctx[s];
    }

    if (neut_species[s]->has_Udrift_init_func) {
      gk.neut_species[s].projection.udrift = gkyl_lw_eval_cb;
      gk.neut_species[s].projection.ctx_udrift = &app_lw->neut_Udrift_init_func_ctx[s];
    }

    if (neut_species[s]->has_temp_init_func) {
      gk.neut_species[s].projection.temp = gkyl_lw_eval_cb;
      gk.neut_species[s].projection.ctx_temp = &app_lw->neut_temp_init_func_ctx[s];
    }
  }

  with_lua_tbl_key(L, "field") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct gyrokinetic_field_lw *gkf = lua_touserdata(L, -1);

      if (gkf->magic == GYROKINETIC_FIELD_DEFAULT) {
        gk.field = gkf->gk_field;
      }
    }
  }

  // Create parallelism.
  struct gkyl_comm *comm = 0;

  for (int d = 0; d < cdim; d++) {
    gk.parallelism.cuts[d] = cuts[d]; 
  }

  struct gkyl_tool_args *args = gkyl_tool_args_new(L);
  struct script_cli script_cli = gk_parse_script_cli(args);

#ifdef GKYL_HAVE_MPI
  if (script_cli.use_gpu && script_cli.use_mpi) {
#ifdef GKYL_HAVE_NCCL
    with_lua_global(L, "GKYL_MPI_COMM") {
      if (lua_islightuserdata(L, -1)) {
        struct { MPI_Comm comm; } *lw_mpi_comm_world = lua_touserdata(L, -1);
        MPI_Comm mpi_comm = lw_mpi_comm_world->comm;

        int nrank = 1; // Number of processors in simulation.
        MPI_Comm_size(mpi_comm, &nrank);

        comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
            .mpi_comm = mpi_comm,
          }
        );
      }
    }
#else
    printf("Using CUDA and MPI together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (script_cli.use_mpi) {
    with_lua_global(L, "GKYL_MPI_COMM") {
      if (lua_islightuserdata(L, -1)) {
        struct { MPI_Comm comm; } *lw_mpi_comm_world = lua_touserdata(L, -1);
        MPI_Comm mpi_comm = lw_mpi_comm_world->comm;

        int nrank = 1; // Number of processors in simulation.
        MPI_Comm_size(mpi_comm, &nrank);

        comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
            .mpi_comm = mpi_comm,
          }
        );
      }
    }
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = script_cli.use_gpu,
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = script_cli.use_gpu,
    }
  );
#endif

  gk.parallelism.comm = comm;
  gk.parallelism.use_gpu = script_cli.use_gpu;

  int rank;
  gkyl_comm_get_rank(comm, &rank);

  int comm_sz;
  gkyl_comm_get_size(comm, &comm_sz);

  int tot_cuts = 1;
  for (int d = 0; d < cdim; d++) {
    tot_cuts *= cuts[d];
  }

  if (tot_cuts != comm_sz) {
    printf("tot_cuts = %d (%d)\n", tot_cuts, comm_sz);
    luaL_error(L, "Number of ranks and cuts do not match!");
  }
  
  app_lw->app = gkyl_gyrokinetic_app_new(&gk);

  gkyl_comm_release(comm);

  // Create Lua userdata.
  struct gyrokinetic_app_lw **l_app_lw = lua_newuserdata(L, sizeof(struct gyrokinetic_app_lw*));
  *l_app_lw = app_lw; // Point it to the Lua app pointer.

  // Set metatable.
  luaL_getmetatable(L, GYROKINETIC_APP_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Apply initial conditions. (time) -> bool.
static int
gk_app_apply_ic(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->t_start);
  gkyl_gyrokinetic_app_apply_ic(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to species. (sidx, time) -> bool.
static int
gk_app_apply_ic_species(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double t0 = luaL_optnumber(L, 3, app_lw->t_start);
  gkyl_gyrokinetic_app_apply_ic_species(app_lw->app, sidx, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated moments. (tm) -> bool.
static int
gk_app_calc_integrated_mom(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_gyrokinetic_app_calc_integrated_mom(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated field energy (L2 norm of each field
// component). (tm) -> bool.
static int
gk_app_calc_field_energy(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_gyrokinetic_app_calc_field_energy(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Write solution (field and species) to file (time, frame) -> bool.
static int
gk_app_write(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_gyrokinetic_app_write(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write field to file (time, frame) -> bool.
static int
gk_app_write_field(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_gyrokinetic_app_write_field(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write species solution to file (sidx, time, frame) -> bool.
static int
gk_app_write_species(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double tm = luaL_checknumber(L, 3);
  int frame = luaL_checkinteger(L, 4);
  gkyl_gyrokinetic_app_write_species(app_lw->app, sidx, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write diagnostic moments to file (time, frame) -> bool.
static int
gk_app_write_mom(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_gyrokinetic_app_write_mom(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated moments to file () -> bool.
static int
gk_app_write_integrated_mom(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  gkyl_gyrokinetic_app_write_integrated_mom(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}


// Write integrated field energy to file () -> bool.
static int
gk_app_write_field_energy(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  gkyl_gyrokinetic_app_write_field_energy(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write simulation statistics to JSON. () -> bool.
static int
gk_app_stat_write(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  gkyl_gyrokinetic_app_stat_write(app_lw->app);

  lua_pushboolean(L, status);
  return 1;  
}

// Write data from simulation to file.
static void
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_gyrokinetic_app_write(app, t_curr, frame);
    gkyl_gyrokinetic_app_write_field_energy(app);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
  }
}

// Calculate and append field energy to dynvector.
static void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_gyrokinetic_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
  }
}

// Calculate and append integrated moments to dynvector.
static void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_gyrokinetic_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
  }
}

// Step message context.
struct step_message_trigs {
  int log_count; // Number of times logging called.
  int tenth, p1c; 
  struct gkyl_tm_trigger log_trig; // 10% trigger.
  struct gkyl_tm_trigger log_trig_1p; // 1% trigger.
};

// Write log message to console.
static void
write_step_message(const struct gkyl_gyrokinetic_app *app, struct step_message_trigs *trigs, int step, double t_curr, double dt_next)
{
  if (gkyl_tm_trigger_check_and_bump(&trigs->log_trig, t_curr)) {
    if (trigs->log_count > 0) {
      gkyl_gyrokinetic_app_cout(app, stdout, " Step %d at time %#11.8g.  Time-step  %.6e.  Completed %g%s\n", step, t_curr, dt_next, trigs->tenth * 10.0, "%");
    }
    else {
      trigs->log_count += 1;
    }
    
    trigs->tenth += 1;
  }
  if (gkyl_tm_trigger_check_and_bump(&trigs->log_trig_1p, t_curr)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "%d", trigs->p1c);
    trigs->p1c = (trigs->p1c+1) % 10;
  }
}

static void
show_help(const struct gkyl_gyrokinetic_app *app)
{
  gkyl_gyrokinetic_app_cout(app, stdout, "Gyrokinetic script takes the following arguments:\n");
  gkyl_gyrokinetic_app_cout(app, stdout, " -h   Print this help message and exit\n");
  gkyl_gyrokinetic_app_cout(app, stdout, " -sN  Only run N steps of simulation\n");
  gkyl_gyrokinetic_app_cout(app, stdout, " -S   Do not initialize MPI\n");
  gkyl_gyrokinetic_app_cout(app, stdout, " -G   Do not initialize CUDA\n");
  gkyl_gyrokinetic_app_cout(app, stdout, " -m   Run memory tracer\n");
  gkyl_gyrokinetic_app_cout(app, stdout, " -V   Show verbose output\n");
  gkyl_gyrokinetic_app_cout(app, stdout, " -rN  Restart simulation from frame N\n");

  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
}

// Run simulation. (num_steps) -> bool. num_steps is optional.
static int
gk_app_run(lua_State *L)
{
  bool ret_status = true;

  // Create app object.
  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;
  struct gkyl_gyrokinetic_app *app = app_lw->app;

  // Parse command lines arguments passed to input file.
  struct gkyl_tool_args *args = gkyl_tool_args_new(L);

  struct script_cli script_cli = gk_parse_script_cli(args);
  if (script_cli.help) {
    show_help(app);
    gkyl_tool_args_release(script_cli.rest);
    gkyl_tool_args_release(args);
    goto freeresources;
  }

  gkyl_tool_args_release(script_cli.rest);
  gkyl_tool_args_release(args);

  // Initial and final simulation times.
  double t_curr = app_lw->t_start, t_end = app_lw->t_end;
  long num_steps = script_cli.num_steps;

  gkyl_gyrokinetic_app_cout(app, stdout, "Initializing Gyrokinetic Simulation ...\n");

  if (script_cli.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  // Initialize simulation.
  bool is_restart = script_cli.is_restart;
  int restart_frame = script_cli.restart_frame;

  int frame_curr = 0;
  if (is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_gyrokinetic_app_apply_ic(app, t_curr);
  }

  int num_frames = app_lw->num_frames;
  int field_energy_calcs = app_lw->field_energy_calcs;
  int integrated_mom_calcs = app_lw->integrated_mom_calcs;
  // Triggers for IO and logging.
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames, .tcurr = frame_curr * (t_end / num_frames), .curr = frame_curr };
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  struct step_message_trigs m_trig = {
    .log_count = 0,
    .tenth = t_curr > 0.0 ?  (int) floor(t_curr / t_end * 10.0) : 0.0,
    .p1c = t_curr > 0.0 ?  (int) floor(t_curr / t_end * 100.0) % 10 : 0.0,
    .log_trig = { .dt = t_end / 10.0, .tcurr = t_curr },
    .log_trig_1p = { .dt = t_end / 100.0, .tcurr = t_curr },
  };

  struct timespec tm_ic0 = gkyl_wall_clock();
  // Initialize simulation.
  calc_field_energy(&fe_trig, app, t_curr, false);
  calc_integrated_mom(&im_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);

  gkyl_gyrokinetic_app_cout(app, stdout, "Initialization completed in %g sec\n\n", gkyl_time_diff_now_sec(tm_ic0));

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = app_lw->dt_failure_tol;
  int num_failures = 0, num_failures_max = app_lw->num_failures_max;

  bool use_verbose = script_cli.use_verbose;

  long step = 1;
  while ((t_curr < t_end) && (step <= num_steps)) {
    if (use_verbose) {
      gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    }
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    if (use_verbose) {
      gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    }

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_gyrokinetic_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_gyrokinetic_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        write_data(&io_trig, app, t_curr, true);

        break;
      }
    }
    else {
      num_failures = 0;
    }

    if (!use_verbose) {
      write_step_message(app, &m_trig, step, t_curr, status.dt_suggested);
    }

    step += 1;
  }

  calc_field_energy(&fe_trig, app, t_curr, false);
  calc_integrated_mom(&im_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);
  gkyl_gyrokinetic_app_stat_write(app);

  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app);

  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  // Plasma species counters. 
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of species CFL reduction calls %ld\n", stat.n_species_omega_cfl);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of moment updater calls %ld\n", stat.n_mom);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of diagnostics computed %ld\n", stat.n_diag);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls for phase space fields %ld\n", stat.n_io);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls for configuration space fields %ld\n", stat.n_diag_io);
  // Neutral species counters. 
  if (app->num_neut_species > 0) {
  
    gkyl_gyrokinetic_app_cout(app, stdout, "Number of neutral species CFL reduction calls %ld\n", stat.n_neut_species_omega_cfl);
    gkyl_gyrokinetic_app_cout(app, stdout, "Number of neutral moment updater calls %ld\n", stat.n_neut_mom);
    gkyl_gyrokinetic_app_cout(app, stdout, "Number of neutral diagnostics computed %ld\n", stat.n_neut_diag);
    gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls for neutrals' phase space fields %ld\n", stat.n_neut_io);
    gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls for neutrals' configuration space fields %ld\n", stat.n_neut_diag_io);
  }
  // Field counters. 
  if (app->field->update_field > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Number of field diagnostics computed %ld\n", stat.n_diag);
    gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls for fields %ld\n", stat.n_field_io);
    gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls for fields' diagnostics %ld\n", stat.n_field_diag_io);
  }
  // Timers.
  gkyl_gyrokinetic_app_print_timings(app, stdout);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species initial conditions took %lg secs\n", stat.init_species_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species LTE computation, including corrections, took %lg secs\n", stat.species_lte_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Neutral species initial conditions took %lg secs\n", stat.init_neut_species_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Neutral species LTE computation, including corrections, took %lg secs\n", stat.neut_species_lte_tm);

freeresources:

  lua_pushboolean(L, ret_status);
  return 1;
}

// Clean up memory allocated for simulation.
static int
gk_app_gc(lua_State *L)
{
  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  gkyl_gyrokinetic_app_release(app_lw->app);
  gkyl_free(*l_app_lw);
  
  return 0;
}

// App constructor.
static struct luaL_Reg gk_app_ctor[] = {
  { "new",  gk_app_new },
  { 0, 0 }
};

// App methods.
static struct luaL_Reg gk_app_funcs[] = {
  { "apply_ic", gk_app_apply_ic },
  { "apply_ic_species", gk_app_apply_ic_species },
  { "calc_integrated_mom", gk_app_calc_integrated_mom },
  { "calc_field_energy", gk_app_calc_field_energy },
  { "write", gk_app_write },
  { "write_field", gk_app_write_field },
  { "write_species", gk_app_write_species },
  { "write_mom", gk_app_write_mom },
  { "write_integrated_mom", gk_app_write_integrated_mom },
  { "write_field_energy", gk_app_write_field_energy },
  { "stat_write", gk_app_stat_write },
  { "run", gk_app_run },
  { 0, 0 }
};

static void
app_openlibs(lua_State *L)
{
  // Register top-level App.
  do {
    luaL_newmetatable(L, GYROKINETIC_APP_METATABLE_NM);

    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, gk_app_gc);
    lua_settable(L, -3);

    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, gk_app_funcs);
    
    luaL_register(L, "G0.Gyrokinetic.App", gk_app_ctor);
    
  }
  while (0);

  // Register Species input struct.
  do {
    luaL_newmetatable(L, GYROKINETIC_SPECIES_METATABLE_NM);
    luaL_register(L, "G0.Gyrokinetic.Species", gk_species_ctor);
  }
  while (0);

  // Register Neutral Species input struct.
  do {
    luaL_newmetatable(L, GYROKINETIC_NEUTRAL_SPECIES_METATABLE_NM);
    luaL_register(L, "G0.Gyrokinetic.NeutralSpecies", gk_neutral_species_ctor);
  }
  while (0);

  // Register Field input struct.
  do {
    luaL_newmetatable(L, GYROKINETIC_FIELD_METATABLE_NM);
    luaL_register(L, "G0.Gyrokinetic.Field", gk_field_ctor);
  }
  while (0);
}

void
gkyl_gyrokinetic_lw_openlibs(lua_State *L)
{
  // Register types for gyrokinetic FEM boundary conditions, geometry type, position map type, field ID, radiation ID, Te min model type, reaction type, ion type, and self-reaction type.
  gkyl_register_gyrokinetic_fem_bc_types(L);
  gkyl_register_gyrokinetic_geometry_types(L);
  gkyl_register_gyrokinetic_position_map_types(L);
  gkyl_register_gyrokinetic_collisionless_types(L);
  gkyl_register_gyrokinetic_field_types(L);
  gkyl_register_gyrokinetic_radiation_types(L);
  gkyl_register_gyrokinetic_radiation_Te_types(L);
  gkyl_register_gyrokinetic_reaction_types(L);
  gkyl_register_gyrokinetic_ion_types(L);
  gkyl_register_gyrokinetic_self_reaction_types(L);
  gkyl_register_gyrokinetic_bc_types(L);
  
  app_openlibs(L);
}

#endif
