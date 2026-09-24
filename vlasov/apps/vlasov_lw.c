#ifdef GKYL_HAVE_LUA

#include <gkyl_alloc.h>
#include <gkyl_eqn_type.h>
#include <gkyl_lua_utils.h>
#include <gkyl_lw_priv.h>
#include <gkyl_null_comm.h>
#include <gkyl_vlasov.h>
#include <gkyl_vlasov_lw.h>
#include <gkyl_vlasov_priv.h>
#include <gkyl_wv_advect.h>
#include <gkyl_wv_canonical_pb_fluid.h>
#include <gkyl_wv_euler.h>
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

// Poisson boundary conditions -> enum map.
static const struct gkyl_str_int_pair poisson_bcs[] = {
  { "bcPeriodic", GKYL_POISSON_PERIODIC },
  { "bcDirichlet", GKYL_POISSON_DIRICHLET },
  { "bcNeumann", GKYL_POISSON_NEUMANN },
  { "bcRobin", GKYL_POISSON_ROBIN },
  { 0, 0 }
};

// Vlasov projection type -> enum map.
static const struct gkyl_str_int_pair projection_type[] = {
  { "Func", GKYL_PROJ_FUNC },
  { "MaxwellianPrimitive", GKYL_PROJ_MAXWELLIAN_PRIM },
  { "MaxwellianLab", GKYL_PROJ_MAXWELLIAN_LAB },
  { "BiMaxwellian", GKYL_PROJ_BIMAXWELLIAN },
  { "LTE", GKYL_PROJ_VLASOV_LTE },
  { 0, 0 }
};

// Vlasov model type -> enum map.
static const struct gkyl_str_int_pair model_type[] = {
  { "Default", GKYL_MODEL_DEFAULT },
  { "SR", GKYL_MODEL_SR },
  { "CanonicalPB", GKYL_MODEL_CANONICAL_PB },
  { "CanonicalPBGR", GKYL_MODEL_CANONICAL_PB_GR },
  { 0, 0 }
};

// Vlasov collision type -> enum map.
static const struct gkyl_str_int_pair collision_type[] = {
  { "None", GKYL_NO_COLLISIONS },
  { "BGK", GKYL_BGK_COLLISIONS },
  { "LBO", GKYL_LBO_COLLISIONS },
  { "FPO", GKYL_FPO_COLLISIONS },
  { 0, 0 }
};

// Vlasov source type -> enum map.
static const struct gkyl_str_int_pair source_type[] = {
  { "None", GKYL_NO_SOURCE },
  { "Func", GKYL_FUNC_SOURCE },
  { "Proj", GKYL_PROJ_SOURCE },
  { "BoundaryFlux", GKYL_BFLUX_SOURCE },
  { 0, 0 }
};

void
gkyl_register_poisson_bc_types(lua_State *L)
{
  register_types(L, poisson_bcs, "PoissonBc");
}

void
gkyl_register_vlasov_projection_types(lua_State *L)
{
  register_types(L, projection_type, "Projection");
}

void
gkyl_register_vlasov_model_types(lua_State *L)
{
  register_types(L, model_type, "Model");
}

void
gkyl_register_vlasov_collision_types(lua_State *L)
{
  register_types(L, collision_type, "Collisions");
}

void
gkyl_register_vlasov_source_types(lua_State *L)
{
  register_types(L, source_type, "Source");
}

// Magic IDs for use in distinguishing various species and field types.
enum vlasov_magic_ids {
  VLASOV_SPECIES_DEFAULT = 100, // Non-relativistic kinetic species.
  VLASOV_FIELD_DEFAULT, // Maxwell equations.
  VLASOV_FLUID_SPECIES_DEFAULT, // Fluid species.
  VLASOV_EQN_DEFAULT, // Equation object.
};

// Metatable name for equation object input struct.
#define VLASOV_WAVE_EQN_METATABLE_NM "GkeyllZero.App.Vlasov.Eq"

// Methods for manipulating gkyl_wv_eqn objects.

// Lua userdata object for constructing wave equation objects.
struct wv_eqn_lw {
  int magic; // This must be the first element in the struct.
  struct gkyl_wv_eqn *eqn; // Equation object.
};

// Clean up memory allocated for equation object.
static int
wv_eqn_lw_gc(lua_State *L)
{
  struct wv_eqn_lw **l_wv_lw = GKYL_CHECK_UDATA(L, VLASOV_WAVE_EQN_METATABLE_NM);
  struct wv_eqn_lw *wv_lw = *l_wv_lw;

  gkyl_wv_eqn_release(wv_lw->eqn);
  gkyl_free(*l_wv_lw);
  
  return 0;
}

// Acquire equation object.
static struct gkyl_wv_eqn*
wv_eqn_get(lua_State *L)
{
  struct wv_eqn_lw **l_wv_lw = luaL_checkudata(L, -1, VLASOV_WAVE_EQN_METATABLE_NM);
  struct wv_eqn_lw *wv_lw = *l_wv_lw;

  return wv_lw->eqn;
}

/* *************** */
/* Euler Equations */
/* *************** */

// Euler.new { gasGamma = 1.4, rpType = G0.EulerRP.Roe }
// where rpType is one of G0.EulerRP.Roe, G0.EulerRP.Lax, G0.EulerRP.HLL or G0.EulerRP.HLLC.
static int
eqn_euler_lw_new(lua_State *L)
{
  struct wv_eqn_lw *euler_lw = gkyl_malloc(sizeof(*euler_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 1.4);
  enum gkyl_wv_euler_rp rp_type = glua_tbl_get_integer(L, "rpType", WV_EULER_RP_ROE);

  euler_lw->magic = VLASOV_EQN_DEFAULT;
  euler_lw->eqn = gkyl_wv_euler_inew( &(struct gkyl_wv_euler_inp) {
      .gas_gamma = gas_gamma,
      .rp_type = rp_type,
      .use_gpu = false
    }
  );

  // Create Lua userdata.
  struct wv_eqn_lw **l_euler_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_euler_lw = euler_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, VLASOV_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_euler_ctor[] = {
  { "new", eqn_euler_lw_new },
  { 0, 0 }
};

/* ****************************** */
/* Advection Equation */
/* ****************************** */

// Advect.new { }
static int
eqn_advect_lw_new(lua_State *L)
{
  struct wv_eqn_lw *advect_lw = gkyl_malloc(sizeof(*advect_lw));

  advect_lw->magic = VLASOV_EQN_DEFAULT;
  // Set a constant speed of 1.0; Advection velocity in DG advection equation
  // is handled by the app_advect function initialized in the FluidSpecies table. 
  advect_lw->eqn = gkyl_wv_advect_new(1.0, false);

  // Create Lua userdata.
  struct wv_eqn_lw **l_advect_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_advect_lw = advect_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, VLASOV_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_advect_ctor[] = {
  { "new", eqn_advect_lw_new },
  { 0, 0 }
};

/* ****************************** */
/* incompressible Euler Equations */
/* ****************************** */

// IncompressEuler.new { }
static int
eqn_incompress_euler_lw_new(lua_State *L)
{
  struct wv_eqn_lw *incompress_euler_lw = gkyl_malloc(sizeof(*incompress_euler_lw));

  incompress_euler_lw->magic = VLASOV_EQN_DEFAULT;
  incompress_euler_lw->eqn = gkyl_wv_can_pb_incompress_euler_new();

  // Create Lua userdata.
  struct wv_eqn_lw **l_incompress_euler_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_incompress_euler_lw = incompress_euler_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, VLASOV_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_incompress_euler_ctor[] = {
  { "new", eqn_incompress_euler_lw_new },
  { 0, 0 }
};

/* *********************** */
/* Hasegawa-Mima Equations */
/* *********************** */

// HasegawaMima.new { }
static int
eqn_hasegawa_mima_lw_new(lua_State *L)
{
  struct wv_eqn_lw *hasegawa_mima_lw = gkyl_malloc(sizeof(*hasegawa_mima_lw));

  hasegawa_mima_lw->magic = VLASOV_EQN_DEFAULT;
  hasegawa_mima_lw->eqn = gkyl_wv_can_pb_hasegawa_mima_new();

  // Create Lua userdata.
  struct wv_eqn_lw **l_hasegawa_mima_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_hasegawa_mima_lw = hasegawa_mima_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, VLASOV_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_hasegawa_mima_ctor[] = {
  { "new", eqn_hasegawa_mima_lw_new },
  { 0, 0 }
};

/* *************************** */
/* Hasegawa-Wakatani Equations */
/* *************************** */

// HasegawaWakatani.new { alpha = 1.0, is_modified = false }
// Set is_modified=true to utilized modified Hasegawa-Wakatani system which
// subtracts off zonal component of adiabatic coupling term. 
static int
eqn_hasegawa_wakatani_lw_new(lua_State *L)
{
  struct wv_eqn_lw *hasegawa_wakatani_lw = gkyl_malloc(sizeof(*hasegawa_wakatani_lw));

  double alpha = glua_tbl_get_number(L, "alpha", 1.0);
  bool is_modified = glua_tbl_get_bool(L, "is_modified", false);  

  hasegawa_wakatani_lw->magic = VLASOV_EQN_DEFAULT;
  hasegawa_wakatani_lw->eqn = gkyl_wv_can_pb_hasegawa_wakatani_new(alpha, is_modified);

  // Create Lua userdata.
  struct wv_eqn_lw **l_hasegawa_wakatani_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_hasegawa_wakatani_lw = hasegawa_wakatani_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, VLASOV_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_hasegawa_wakatani_ctor[] = {
  { "new", eqn_hasegawa_wakatani_lw_new },
  { 0, 0 }
};

// Register and load all wave equation objects.
static void
eqn_openlibs(lua_State *L)
{
  luaL_newmetatable(L, VLASOV_WAVE_EQN_METATABLE_NM);

  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, wv_eqn_lw_gc);
  lua_settable(L, -3);

  luaL_register(L, "G0.Vlasov.Eq.Euler", eqn_euler_ctor);
  luaL_register(L, "G0.Vlasov.Eq.LinearAdvection", eqn_advect_ctor);
  luaL_register(L, "G0.Vlasov.Eq.IncompressEuler", eqn_incompress_euler_ctor);
  luaL_register(L, "G0.Vlasov.Eq.HasegawaMima", eqn_hasegawa_mima_ctor);
  luaL_register(L, "G0.Vlasov.Eq.HasegawaWakatani", eqn_hasegawa_wakatani_ctor);
}

/* *************** */
/* Species methods */
/* *************** */

// Metatable name for species input struct.
#define VLASOV_SPECIES_METATABLE_NM "GkeyllZero.App.Vlasov.Species"

// Lua userdata object for constructing species input.
struct vlasov_species_lw {
  int magic; // This must be first element in the struct.
  
  struct gkyl_vlasov_species vm_species; // Input struct to construct species.
  int vdim; // Velocity space dimensions.
  bool evolve; // Is this species evolved?

  bool has_hamiltonian_func; // Is there a Hamiltonian function?
  struct lua_func_ctx hamiltonian_func_ref; // Lua registry reference to Hamiltonian function.

  bool has_metric_func; // Is there an metric tensor function?
  struct lua_func_ctx metric_func_ref; // Lua registry reference to metric tensor function.

  bool has_inverse_metric_func; // Is there an inverse metric tensor function?
  struct lua_func_ctx inverse_metric_func_ref; // Lua registry reference to inverse metric tensor function.

  bool has_metric_determinant_func; // Is there a metric determinant function?
  struct lua_func_ctx metric_determinant_func_ref; // Lua registry reference to metric determinant function.

  bool output_f_lte; // Should f_lte be written out (for calculating transport coefficients)?

  int num_init; // Number of projection objects.
  enum gkyl_projection_id proj_id[GKYL_MAX_PROJ]; // Projection type.

  bool has_init_func[GKYL_MAX_PROJ]; // Is there an initialization function?
  struct lua_func_ctx init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to initialization function.

  bool has_density_init_func[GKYL_MAX_PROJ]; // Is there a density initialization function?
  struct lua_func_ctx density_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to density initialization function.

  bool has_V_drift_init_func[GKYL_MAX_PROJ]; // Is there a drift velocity initialiation function?
  struct lua_func_ctx V_drift_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to drift velocity initialization function.

  bool has_temp_init_func[GKYL_MAX_PROJ]; // Is there a temperature initialization function?
  struct lua_func_ctx temp_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to temperature initialization function.

  bool correct_all_moms[GKYL_MAX_PROJ]; // Are we correcting all moments in projections, or only density?
  double iter_eps[GKYL_MAX_PROJ]; // Error tolerance for moment fixes in projections (density is always exact).
  int max_iter[GKYL_MAX_PROJ]; // Maximum number of iterations for moment fixes in projections.
  bool use_last_converged[GKYL_MAX_PROJ]; // Use last iteration value in projections regardless of convergence?

  enum gkyl_collision_id collision_id; // Collision type.
  
  bool has_self_nu_func; // Is there a self-collision frequency function?
  struct lua_func_ctx self_nu_func_ref; // Lua registry reference to self-collision frequency function.

  int num_cross_collisions; // Number of species that we cross-collide with.
  char collide_with[GKYL_MAX_SPECIES][128]; // Names of species that we cross-collide with.

  bool collision_correct_all_moms; // Are we correcting all moments in collisions, or only density?
  double collision_iter_eps; // Error tolerance for moment fixes in collisions (density is always exact).
  int collision_max_iter; // Maximum number of iterations for moment fixes in collisions.
  bool fixed_temp_relax; // Are BGK collisions relaxing to a fixed input temperature?
  bool collision_use_last_converged; // Use last iteration value in collisions regardless of convergence?
  bool has_implicit_coll_scheme; // Use implicit scheme for collisions?

  enum gkyl_source_id source_id; // Source type.

  double source_length; // Length used to scale the source function.
  char source_species[128]; // Name of species to use for the source.

  int num_sources; // Number of projection objects in source.
  enum gkyl_projection_id source_proj_id[GKYL_MAX_PROJ]; // Projection type in source.

  bool source_has_init_func[GKYL_MAX_PROJ]; // Is there an initialization function in source?
  struct lua_func_ctx source_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to initialization function in source.

  bool source_has_density_init_func[GKYL_MAX_PROJ]; // Is there a density initialization function in source?
  struct lua_func_ctx source_density_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to density initialization function in source.

  bool source_has_V_drift_init_func[GKYL_MAX_PROJ]; // Is there a drift velocity initialization function in source?
  struct lua_func_ctx source_V_drift_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to drift velocity initialization function in source.

  bool source_has_temp_init_func[GKYL_MAX_PROJ]; // Is there a temperature initialization function in source?
  struct lua_func_ctx source_temp_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to temperature initialization function in source.

  bool source_correct_all_moms[GKYL_MAX_PROJ]; // Are we correcting all moments in projections, or only density, in source?
  double source_iter_eps[GKYL_MAX_PROJ]; // Error tolerance for moment fixes in projections (density is always exact) in source.
  int source_max_iter[GKYL_MAX_PROJ]; // Maximum number of iterations for moment fixes in projections in source.
  bool source_use_last_converged[GKYL_MAX_PROJ]; // Use last iteration value in projection regardless of convergence in source?
};

static int
vlasov_species_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_vlasov_species vm_species = { };

  vm_species.model_id = glua_tbl_get_integer(L, "modelID", 0);
  
  vm_species.charge = glua_tbl_get_number(L, "charge", 0.0);
  vm_species.mass = glua_tbl_get_number(L, "mass", 1.0);

  with_lua_tbl_tbl(L, "cells") {
    vdim = glua_objlen(L);

    for (int d = 0; d < vdim; d++) {
      vm_species.cells[d] = glua_tbl_iget_integer(L, d + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "lower") {
    for (int d = 0; d < vdim; d++) {
      vm_species.lower[d] = glua_tbl_iget_number(L, d + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d = 0; d < vdim; d++) {
      vm_species.upper[d] = glua_tbl_iget_number(L, d + 1, 0);
    }
  }

  bool evolve = glua_tbl_get_bool(L, "evolve", true);

  with_lua_tbl_tbl(L, "diagnostics") {
    int num_diag_moments = glua_objlen(L);

    for (int i = 0; i < num_diag_moments; i ++) {
      vm_species.diag_moments[i] = glua_tbl_iget_integer(L, i+1, 0);
    }

    vm_species.num_diag_moments = num_diag_moments;
  }

  with_lua_tbl_tbl(L, "bcx") {
    with_lua_tbl_tbl(L, "lower") {
      vm_species.bcx.lower.type = glua_tbl_get_integer(L, "type", 0);
    }
    
    with_lua_tbl_tbl(L, "upper") {
      vm_species.bcx.upper.type = glua_tbl_get_integer(L, "type", 0);
    }
  }

  with_lua_tbl_tbl(L, "bcy") {
    with_lua_tbl_tbl(L, "lower") {
      vm_species.bcy.lower.type = glua_tbl_get_integer(L, "type", 0);
    }

    with_lua_tbl_tbl(L, "upper") {
      vm_species.bcy.upper.type = glua_tbl_get_integer(L, "type", 0);
    }
  }

  with_lua_tbl_tbl(L, "bcz") {
    with_lua_tbl_tbl(L, "lower") {
      vm_species.bcz.lower.type = glua_tbl_get_integer(L, "type", 0);
    }

    with_lua_tbl_tbl(L, "upper") {
      vm_species.bcz.upper.type = glua_tbl_get_integer(L, "type", 0);
    }
  }

  bool has_hamiltonian_func = false;
  int hamiltonian_func_ref = LUA_NOREF;

  bool has_metric_func = false;
  int metric_func_ref = LUA_NOREF;

  bool has_inverse_metric_func = false;
  int inverse_metric_func_ref = LUA_NOREF;

  bool has_metric_determinant_func = false;
  int metric_determinant_func_ref = LUA_NOREF;

  if (glua_tbl_get_func(L, "hamiltonian")) {
    hamiltonian_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_hamiltonian_func = true;
  }

  if (glua_tbl_get_func(L, "metric")) {
    metric_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_metric_func = true;
  }

  if (glua_tbl_get_func(L, "inverseMetric")) {
    inverse_metric_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_inverse_metric_func = true;
  }

  if (glua_tbl_get_func(L, "metricDeterminant")) {
    metric_determinant_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_metric_determinant_func = true;
  }

  bool output_f_lte = glua_tbl_get_bool(L, "outputfLTE", false);
  
  enum gkyl_projection_id proj_id[GKYL_MAX_PROJ];

  bool has_init_func[GKYL_MAX_PROJ];
  int init_func_ref[GKYL_MAX_PROJ];

  bool has_density_init_func[GKYL_MAX_PROJ];
  int density_init_func_ref[GKYL_MAX_PROJ];

  bool has_V_drift_init_func[GKYL_MAX_PROJ];
  int V_drift_init_func_ref[GKYL_MAX_PROJ];

  bool has_temp_init_func[GKYL_MAX_PROJ];
  int temp_init_func_ref[GKYL_MAX_PROJ];
  
  bool correct_all_moms[GKYL_MAX_PROJ];
  double iter_eps[GKYL_MAX_PROJ];
  int max_iter[GKYL_MAX_PROJ];
  bool use_last_converged[GKYL_MAX_PROJ];

  int num_init = glua_tbl_get_integer(L, "numInit", 0);

  with_lua_tbl_tbl(L, "projections") {
    for (int i = 0; i < num_init; i++) {
      if (glua_tbl_iget_tbl(L, i + 1)) {
        proj_id[i] = glua_tbl_get_integer(L, "projectionID", 0);

        init_func_ref[i] = LUA_NOREF;
        has_init_func[i] = false;
        if (glua_tbl_get_func(L, "init")) {
          init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
          has_init_func[i] = true;
        }

        density_init_func_ref[i] = LUA_NOREF;
        has_density_init_func[i] = false;
        if (glua_tbl_get_func(L, "densityInit")) {
          density_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
          has_density_init_func[i] = true;
        }

        V_drift_init_func_ref[i] = LUA_NOREF;
        has_V_drift_init_func[i] = false;
        if (glua_tbl_get_func(L, "driftVelocityInit")) {
          V_drift_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
          has_V_drift_init_func[i] = true;
        }

        temp_init_func_ref[i] = LUA_NOREF;
        has_temp_init_func[i] = false;
        if (glua_tbl_get_func(L, "temperatureInit")) {
          temp_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
          has_temp_init_func[i] = true;
        }

        correct_all_moms[i] = glua_tbl_get_bool(L, "correctAllMoments", true);
        iter_eps[i] = glua_tbl_get_number(L, "iterationEpsilon", pow(10.0, -12.0));
        max_iter[i] = glua_tbl_get_integer(L, "maxIterations", 100);
        use_last_converged[i] = glua_tbl_get_bool(L, "useLastConverged", true);

        lua_pop(L, 1);
      }
    }
  }

  enum gkyl_collision_id collision_id = GKYL_NO_COLLISIONS;

  bool has_self_nu_func = false;
  int self_nu_func_ref = LUA_NOREF;

  int num_cross_collisions = 0;
  char collide_with[GKYL_MAX_SPECIES][128];

  bool collision_correct_all_moms = false;
  double collision_iter_eps = pow(10.0, -12.0);
  int collision_max_iter = 100;
  bool fixed_temp_relax = false;
  bool collision_use_last_converged = true;
  bool has_implicit_coll_scheme = false;

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

    collision_correct_all_moms = glua_tbl_get_bool(L, "correctAllMoments", true);
    collision_iter_eps = glua_tbl_get_number(L, "iterationEpsilon", pow(10.0, -12.0));
    collision_max_iter = glua_tbl_get_integer(L, "maxIterations", 100);
    fixed_temp_relax = glua_tbl_get_bool(L, "fixedTempRelax", false);
    collision_use_last_converged = glua_tbl_get_bool(L, "useLastConverged", true);
    has_implicit_coll_scheme = glua_tbl_get_bool(L, "useImplicitCollisionScheme", false);
  }

  enum gkyl_source_id source_id = GKYL_NO_SOURCE;

  double source_length = 1.0;
  char source_species[128] = { '\0' };

  int num_sources = 0;
  enum gkyl_projection_id source_proj_id[GKYL_MAX_PROJ];

  bool source_has_init_func[GKYL_MAX_PROJ];
  int source_init_func_ref[GKYL_MAX_PROJ];

  bool source_has_density_init_func[GKYL_MAX_PROJ];
  int source_density_init_func_ref[GKYL_MAX_PROJ];

  bool source_has_V_drift_init_func[GKYL_MAX_PROJ];
  int source_V_drift_init_func_ref[GKYL_MAX_PROJ];

  bool source_has_temp_init_func[GKYL_MAX_PROJ];
  int source_temp_init_func_ref[GKYL_MAX_PROJ];

  bool source_correct_all_moms[GKYL_MAX_PROJ];
  double source_iter_eps[GKYL_MAX_PROJ];
  int source_max_iter[GKYL_MAX_PROJ];
  bool source_use_last_converged[GKYL_MAX_PROJ];

  with_lua_tbl_tbl(L, "source") {
    source_id = glua_tbl_get_integer(L, "sourceID", 0);

    source_length = glua_tbl_get_number(L, "sourceLength", 1.0);
    const char* source_species_char = glua_tbl_get_string(L, "sourceSpecies", "");
    strcpy(source_species, source_species_char);

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

          source_V_drift_init_func_ref[i] = LUA_NOREF;
          source_has_V_drift_init_func[i] = false;
          if (glua_tbl_get_func(L, "driftVelocityInit")) {
            source_V_drift_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
            source_has_V_drift_init_func[i] = true;
          }

          source_temp_init_func_ref[i] = LUA_NOREF;
          source_has_temp_init_func[i] = false;
          if (glua_tbl_get_func(L, "temperatureInit")) {
            source_temp_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
            source_has_temp_init_func[i] = true;
          }

          source_correct_all_moms[i] = glua_tbl_get_bool(L, "correctAllMoms", true);
          source_iter_eps[i] = glua_tbl_get_number(L, "iterationEpsilon", pow(10.0, -12.0));
          source_max_iter[i] = glua_tbl_get_integer(L, "maxIterations", 100);
          source_use_last_converged[i] = glua_tbl_get_bool(L, "useLastConverged", true);

          lua_pop(L, 1);
        }
      }
    }
  }
  
  struct vlasov_species_lw *vms_lw = lua_newuserdata(L, sizeof(*vms_lw));
  vms_lw->magic = VLASOV_SPECIES_DEFAULT;
  vms_lw->vdim = vdim;
  vms_lw->evolve = evolve;
  vms_lw->vm_species = vm_species;

  vms_lw->has_hamiltonian_func = has_hamiltonian_func;
  vms_lw->hamiltonian_func_ref = (struct lua_func_ctx) {
    .func_ref = hamiltonian_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1,
    .L = L,
  };

  vms_lw->has_metric_func = has_metric_func;
  vms_lw->metric_func_ref = (struct lua_func_ctx) {
    .func_ref = metric_func_ref,
    .ndim = 0, // This will be set later.
    .nret = (vdim * (vdim + 1)) / 2,
    .L = L,
  };


  vms_lw->has_inverse_metric_func = has_inverse_metric_func;
  vms_lw->inverse_metric_func_ref = (struct lua_func_ctx) {
    .func_ref = inverse_metric_func_ref,
    .ndim = 0, // This will be set later.
    .nret = (vdim * (vdim + 1)) / 2,
    .L = L,
  };

  vms_lw->has_metric_determinant_func = has_metric_determinant_func;
  vms_lw->metric_determinant_func_ref = (struct lua_func_ctx) {
    .func_ref = metric_determinant_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1,
    .L = L,
  };

  vms_lw->output_f_lte = output_f_lte;

  vms_lw->num_init = num_init;
  for (int i = 0; i < num_init; i++) {
    vms_lw->proj_id[i] = proj_id[i];

    vms_lw->has_init_func[i] = has_init_func[i];
    vms_lw->init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };

    vms_lw->has_density_init_func[i] = has_density_init_func[i];
    vms_lw->density_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = density_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };

    vms_lw->has_V_drift_init_func[i] = has_V_drift_init_func[i];
    vms_lw->V_drift_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = V_drift_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = vdim,
      .L = L,
    };

    vms_lw->has_temp_init_func[i] = has_temp_init_func[i];
    vms_lw->temp_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = temp_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };

    vms_lw->correct_all_moms[i] = correct_all_moms[i];
    vms_lw->iter_eps[i] = iter_eps[i];
    vms_lw->max_iter[i] = max_iter[i];
    vms_lw->use_last_converged[i] = use_last_converged[i];
  }

  vms_lw->source_id = source_id;
  vms_lw->num_sources = num_sources;

  strcpy(vms_lw->source_species, source_species);
  vms_lw->source_length = source_length;
  
  for (int i = 0; i < num_sources; i++) {
    vms_lw->source_proj_id[i] = source_proj_id[i];

    vms_lw->source_has_init_func[i] = source_has_init_func[i];
    vms_lw->source_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = source_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };

    vms_lw->source_has_density_init_func[i] = source_has_density_init_func[i];
    vms_lw->source_density_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = source_density_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };

    vms_lw->source_has_V_drift_init_func[i] = source_has_V_drift_init_func[i];
    vms_lw->source_V_drift_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = source_V_drift_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = vdim,
      .L = L,
    };

    vms_lw->source_has_temp_init_func[i] = source_has_temp_init_func[i];
    vms_lw->source_temp_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = source_temp_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };
  }

  vms_lw->collision_id = collision_id;

  vms_lw->has_self_nu_func = has_self_nu_func;
  vms_lw->self_nu_func_ref = (struct lua_func_ctx) {
    .func_ref = self_nu_func_ref,
    .ndim = 0,
    .nret = 1,
    .L = L,
  };

  vms_lw->num_cross_collisions = num_cross_collisions;
  for (int i = 0; i < num_cross_collisions; i++) {
    strcpy(vms_lw->collide_with[i], collide_with[i]);
  }

  vms_lw->collision_correct_all_moms = collision_correct_all_moms;
  vms_lw->collision_iter_eps = collision_iter_eps;
  vms_lw->collision_max_iter = collision_max_iter;
  vms_lw->fixed_temp_relax = fixed_temp_relax;
  vms_lw->collision_use_last_converged = collision_use_last_converged;
  vms_lw->has_implicit_coll_scheme = has_implicit_coll_scheme;
  
  // Set metatable.
  luaL_getmetatable(L, VLASOV_SPECIES_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor.
static struct luaL_Reg vm_species_ctor[] = {
  { "new", vlasov_species_lw_new },
  { 0, 0 }
};

/* ********************* */
/* Fluid Species methods */
/* ********************* */

// Metatable name for fluid species input struct.
#define VLASOV_FLUID_SPECIES_METATABLE_NM "GkeyllZero.App.Vlasov.FluidSpecies"

// Lua userdata object for constructing fluid species input.
struct vlasov_fluid_species_lw {
  int magic; // This must be first element in the struct.
  
  struct gkyl_vlasov_fluid_species vlasov_fluid_species; // Input struct to construct fluid species.
  struct lua_func_ctx init_ctx; // Lua registry reference to initialization function.

  bool has_app_advect_func; // Is there an applied advection function?
  struct lua_func_ctx app_advect_func_ref; // Lua registry reference to applied advection function.

  bool has_n0_func; // Is there a background density function?
  struct lua_func_ctx n0_func_ref; // Lua registry reference to background density function.
  
  bool has_diffusion_func; // Is there a diffusion tensor function?
  struct lua_func_ctx diffusion_func_ref; // Lua registry reference to diffusion tensor function.
};

static int
vlasov_fluid_species_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_vlasov_fluid_species vm_fluid_species = { };

  vm_fluid_species.charge = glua_tbl_get_number(L, "charge", 0.0);
  vm_fluid_species.mass = glua_tbl_get_number(L, "mass", 1.0);

  bool has_eqn = false;
  with_lua_tbl_key(L, "equation") {
    vm_fluid_species.equation = wv_eqn_get(L);
    has_eqn = true;
  }

  if (has_eqn) {
    // We need to store a reference to equation object in a global
    // table as otherwise it gets garbage-collected while the
    // simulation is being set up.
    lua_getfield(L, -1, "equation");
    int eq_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    
    lua_getglobal(L, "__vlasov_ref_table");
    int eqtbl_len = glua_objlen(L);
    lua_pushinteger(L, eqtbl_len+1);    
    lua_rawgeti(L, LUA_REGISTRYINDEX, eq_ref);
    lua_rawset(L, -3);

    lua_pop(L, 1);
  }
  else {
    return luaL_error(L, "Fluid species \"equation\" not specfied or incorrect type!");
  }

  int init_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "init")) {
    init_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }
  else {
    return luaL_error(L, "Fluid species must have an \"init\" function for initial conditions!");
  }

  bool has_app_advect_func = false;
  int app_advect_func_ref = LUA_NOREF;

  if (glua_tbl_get_func(L, "appAdvect")) {
    app_advect_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_app_advect_func = true;
  }

  bool has_n0_func = false;
  int n0_func_ref = LUA_NOREF;

  if (glua_tbl_get_func(L, "n0")) {
    n0_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_n0_func = true;
  }

  bool has_diffusion_func = false;
  int diffusion_func_ref = LUA_NOREF;

  with_lua_tbl_tbl(L, "diffusion") {
    vm_fluid_species.diffusion.D = glua_tbl_get_number(L, "diffusionCoefficient", 0.0);
    vm_fluid_species.diffusion.order = glua_tbl_get_integer(L, "diffusionOrder", 2);
    
    if (glua_tbl_get_func(L, "diffusionTensor")) {
      diffusion_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_diffusion_func = true;
    }
  }

  with_lua_tbl_tbl(L, "bcx") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      vm_fluid_species.bcx[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcy") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      vm_fluid_species.bcy[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcz") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      vm_fluid_species.bcz[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  struct vlasov_fluid_species_lw *vmfs_lw = lua_newuserdata(L, sizeof(*vmfs_lw));
  vmfs_lw->magic = VLASOV_FLUID_SPECIES_DEFAULT;
  vmfs_lw->vlasov_fluid_species = vm_fluid_species;

  vmfs_lw->init_ctx = (struct lua_func_ctx) {
    .func_ref = init_ref,
    .ndim = 0, // This will be set later.
    .nret = vm_fluid_species.equation->num_equations,
    .L = L,
  };

  vmfs_lw->has_app_advect_func = has_app_advect_func;
  vmfs_lw->app_advect_func_ref = (struct lua_func_ctx) {
    .func_ref = app_advect_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 3,
    .L = L,
  };

  vmfs_lw->has_n0_func = has_n0_func;
  vmfs_lw->n0_func_ref = (struct lua_func_ctx) {
    .func_ref = n0_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1,
    .L = L,
  };

  vmfs_lw->has_diffusion_func = has_diffusion_func;
  vmfs_lw->diffusion_func_ref = (struct lua_func_ctx) {
    .func_ref = diffusion_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1, // This will be set later.
    .L = L,
  };
  
  // Set metatable.
  luaL_getmetatable(L, VLASOV_FLUID_SPECIES_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor.
static struct luaL_Reg vm_fluid_species_ctor[] = {
  { "new", vlasov_fluid_species_lw_new },
  { 0, 0 }
};

/* ************* */
/* Field methods */
/* ************* */

// Metatable name for field input struct.
#define VLASOV_FIELD_METATABLE_NM "GkeyllZero.App.Vlasov.Field"

// Lua userdata object for constructing field input.
struct vlasov_field_lw {
  int magic; // This must be first element in the struct.
  
  struct gkyl_vlasov_field vm_field; // Input struct to construct field.
  bool evolve; // Is this field evolved?
  struct lua_func_ctx init_ref; // Lua registry reference to initilization function.

  bool has_external_potential_func; // Is there an external potential initialization function?
  struct lua_func_ctx external_potential_func_ref; // Lua registry reference to external potential initialization function.
  bool evolve_external_potential; // Is the external potential evolved?

  bool has_external_field_func; // Is there an external field initialization function?
  struct lua_func_ctx external_field_func_ref; // Lua registry reference to external field initialization function.
  bool evolve_external_field; // Is the external field evolved?

  bool has_applied_current_func; // Is there an applied current initialization function?
  struct lua_func_ctx applied_current_func_ref; // Lua registry reference to applied current initialization function.
  bool evolve_applied_current; // Is the applied current evolved?
};

static int
vlasov_field_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_vlasov_field vm_field = { };

  vm_field.field_id = GKYL_FIELD_E_B;
  
  vm_field.epsilon0 = glua_tbl_get_number(L, "epsilon0", 1.0);
  vm_field.mu0 = glua_tbl_get_number(L, "mu0", 1.0);
  vm_field.elcErrorSpeedFactor = glua_tbl_get_number(L, "elcErrorSpeedFactor", 0.0);
  vm_field.mgnErrorSpeedFactor = glua_tbl_get_number(L, "mgnErrorSpeedFactor", 0.0);
  vm_field.limit_em = glua_tbl_get_bool(L, "limitField", false);
  vm_field.use_ghost_current = glua_tbl_get_bool(L, "useGhostCurrent", false);

  bool evolve = glua_tbl_get_bool(L, "evolve", true);
  vm_field.is_static = !evolve;

  int init_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "init")) {
    init_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }

  with_lua_tbl_tbl(L, "bcx") { 
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      vm_field.bcx[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcy") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      vm_field.bcy[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcz") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      vm_field.bcz[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "poissonBcs") {
    with_lua_tbl_tbl(L, "lowerType") {
      int nbc = glua_objlen(L);
      
      for (int i = 0; i < nbc; i++) {
        vm_field.poisson_bcs.lo_type[i] = glua_tbl_iget_integer(L, i + 1, 0);
      }
    }

    with_lua_tbl_tbl(L, "upperType") {
      int nbc = glua_objlen(L);

      for (int i = 0; i < nbc; i++) {
        vm_field.poisson_bcs.up_type[i] = glua_tbl_iget_integer(L, i + 1, 0);
      }
    }

    with_lua_tbl_tbl(L, "lowerValue") {
      int nbc = glua_objlen(L);

      for (int i = 0; i < nbc; i++) {
        struct gkyl_poisson_bc_value lower_bc = {
          .v = { glua_tbl_iget_number(L, i + 1, 0.0) },
        };

        vm_field.poisson_bcs.lo_value[i] = lower_bc;
      }
    }

    with_lua_tbl_tbl(L, "upperValue") {
      int nbc = glua_objlen(L);

      for (int i = 0; i < nbc; i++) {
        struct gkyl_poisson_bc_value upper_bc = {
          .v = { glua_tbl_iget_number(L, i + 1, 0.0) },
        };

        vm_field.poisson_bcs.up_value[i] = upper_bc;
      }
    }
  }

  bool has_external_potential_func = false;
  int external_potential_func_ref = LUA_NOREF;
  bool evolve_external_potential = false;

  if (glua_tbl_get_func(L, "externalPotentialInit")) {
    external_potential_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_external_potential_func = true;

    evolve_external_potential = glua_tbl_get_bool(L, "evolveExternalPotential", false);
  }

  bool has_external_field_func = false;
  int external_field_func_ref = LUA_NOREF;
  bool evolve_external_field = false;

  if (glua_tbl_get_func(L, "externalFieldInit")) {
    external_field_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_external_field_func = true;

    evolve_external_field = glua_tbl_get_bool(L, "evolveExternalField", false);
  }

  bool has_applied_current_func = false;
  int applied_current_func_ref = LUA_NOREF;
  bool evolve_applied_current = false;

  if (glua_tbl_get_func(L, "appliedCurrent")) {
    applied_current_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_applied_current_func = true;

    evolve_applied_current = glua_tbl_get_bool(L, "evolveAppliedCurrent", false);
  }

  struct vlasov_field_lw *vmf_lw = lua_newuserdata(L, sizeof(*vmf_lw));

  vmf_lw->magic = VLASOV_FIELD_DEFAULT;
  vmf_lw->evolve = evolve;
  vmf_lw->vm_field = vm_field;
  
  vmf_lw->init_ref = (struct lua_func_ctx) {
    .func_ref = init_ref,
    .ndim = 0, // This will be set later.
    .nret = 6,
    .L = L,
  };

  vmf_lw->has_external_potential_func = has_external_potential_func;
  vmf_lw->external_potential_func_ref = (struct lua_func_ctx) {
    .func_ref = external_potential_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 4,
    .L = L,
  };
  vmf_lw->evolve_external_potential = evolve_external_potential;

  vmf_lw->has_external_field_func = has_external_field_func;
  vmf_lw->external_field_func_ref = (struct lua_func_ctx) {
    .func_ref = external_field_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 6,
    .L = L,
  };
  vmf_lw->evolve_external_field = evolve_external_field;

  vmf_lw->has_applied_current_func = has_applied_current_func;
  vmf_lw->applied_current_func_ref = (struct lua_func_ctx) {
    .func_ref = applied_current_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 3,
    .L = L,
  };
  vmf_lw->evolve_applied_current = evolve_applied_current;
  
  // Set metatable.
  luaL_getmetatable(L, VLASOV_FIELD_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Field constructor.
static struct luaL_Reg vm_field_ctor[] = {
  { "new",  vlasov_field_lw_new },
  { 0, 0 }
};

/* *********** */
/* App methods */
/* *********** */

// Metatable name for top-level Vlasov App.
#define VLASOV_APP_METATABLE_NM "GkeyllZero.App.Vlasov"

// Lua userdata object for holding Vlasov app and run parameters.
struct vlasov_app_lw {
  gkyl_vlasov_app *app; // Vlasov app object.

  bool has_hamiltonian_func[GKYL_MAX_SPECIES]; // Is there a Hamiltonian function?
  struct lua_func_ctx hamiltonian_func_ctx[GKYL_MAX_SPECIES]; // Lua registry reference to Hamiltonian function.

  bool has_metric_func[GKYL_MAX_SPECIES]; // Is there an inverse metric tensor function?
  struct lua_func_ctx metric_func_ctx[GKYL_MAX_SPECIES]; // Lua registry reference to inverse metric tensor function.

  bool has_inverse_metric_func[GKYL_MAX_SPECIES]; // Is there an inverse metric tensor function?
  struct lua_func_ctx inverse_metric_func_ctx[GKYL_MAX_SPECIES]; // Lua registry reference to inverse metric tensor function.

  bool has_metric_determinant_func[GKYL_MAX_SPECIES]; // Is there a metric determinant function?
  struct lua_func_ctx metric_determinant_func_ctx[GKYL_MAX_SPECIES]; // Lua registry reference to metric determinant function.

  bool output_f_lte[GKYL_MAX_SPECIES]; // Should f_lte be written out (for calculating transport coefficients)?

  int num_init[GKYL_MAX_SPECIES]; // Number of projection objects.
  enum gkyl_projection_id proj_id[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Projection type.

  bool has_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there an initialization function?
  struct lua_func_ctx init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for initialization function.

  bool has_density_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there a density initialization function?
  struct lua_func_ctx density_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for density initialization function.

  bool has_V_drift_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there a drift velocity initialization function?
  struct lua_func_ctx V_drift_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for drift velocity initialziation function.
  
  bool has_temp_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there a temperature initialization function?
  struct lua_func_ctx temp_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for temperature initialization function.

  bool correct_all_moms[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Are we correcting all moments in projections, or only density?
  double iter_eps[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Error tolerance for moment fixes in projections (density is always exact).
  int max_iter[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Maximum number of iterations for moment fixes in projections.
  bool use_last_converged[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Use last iteration value in projections regardless of convergence?

  enum gkyl_collision_id collision_id[GKYL_MAX_SPECIES]; // Collision type.

  bool has_self_nu_func[GKYL_MAX_SPECIES]; // Is there a self-collision frequency function?
  struct lua_func_ctx self_nu_func_ctx[GKYL_MAX_SPECIES]; // Context for self-collision frequency function.

  int num_cross_collisions[GKYL_MAX_SPECIES]; // Number of species that we cross-collide with.
  char collide_with[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES][128]; // Names of species that we cross-collide with.

  bool collision_correct_all_moms[GKYL_MAX_SPECIES]; // Are we correcting all moments in collisions, or only density?
  double collision_iter_eps[GKYL_MAX_SPECIES]; // Error tolerance for moment fixes in collision (density is always exact).
  int collision_max_iter[GKYL_MAX_SPECIES]; // Maximum number of iterations for moment fixes in collisions.
  bool fixed_temp_relax[GKYL_MAX_SPECIES]; // Are BGK collisions relaxing to a fixed input temperature?
  bool collision_use_last_converged[GKYL_MAX_SPECIES]; // Use last iteration value in collisions regardless of convergence?
  bool has_implicit_coll_scheme[GKYL_MAX_SPECIES]; // Use implicit scheme for collisions?

  enum gkyl_source_id source_id[GKYL_MAX_SPECIES]; // Source type.

  double source_length[GKYL_MAX_SPECIES]; // Length used to scale the source function.
  char source_species[GKYL_MAX_SPECIES][128]; // Name of speccies to use for the source.

  int num_sources[GKYL_MAX_SPECIES]; // Number of projection objects in source.
  enum gkyl_projection_id source_proj_id[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Projection type in source.

  bool source_has_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there an initialization function in source?
  struct lua_func_ctx source_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for initialization function in source.

  bool source_has_density_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there a density initialization function in source?
  struct lua_func_ctx source_density_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for density initialization function in source.

  bool source_has_V_drift_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there a drift velocity initialization function in source?
  struct lua_func_ctx source_V_drift_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for drift velocity initialization function in source.

  bool source_has_temp_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there a temperature initialization function in source?
  struct lua_func_ctx source_temp_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for temperature initialization function in source.

  bool source_correct_all_moms[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Are we correcting all moments in projections, or only density, in source?
  double source_iter_eps[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Error tolerance for moment fixes in projections (density is always exact) in source.
  int source_max_iter[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Maximum number of iterations for moment fixes in projections in source.
  bool source_use_last_converged[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Use last iteration value in projection regardless of convergence in source?

  struct lua_func_ctx fluid_species_init_ctx[GKYL_MAX_SPECIES]; // Function context for fluid species initial conditions.

  struct lua_func_ctx app_advect_func_ctx[GKYL_MAX_SPECIES]; // Function context for fluid species applied advection.
  struct lua_func_ctx n0_func_ctx[GKYL_MAX_SPECIES]; // Function context for fluid species background density.
  struct lua_func_ctx diffusion_func_ctx[GKYL_MAX_SPECIES]; // Function context for fluid species diffusion tensor.

  struct lua_func_ctx field_func_ctx; // Function context for field.
  struct lua_func_ctx external_potential_func_ctx; // Function context for external potential.
  struct lua_func_ctx external_field_func_ctx; // Function context for external field.
  struct lua_func_ctx applied_current_func_ctx; // Function context for applied current.
  
  double t_start, t_end; // Start and end times of simulation.
  int num_frames; // Number of data frames to write.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

// Gets all species objects from the App table, which must on top of
// the stack. The number of species is returned and the appropriate
// pointers set in the species pointer array.
static int
get_species_inp(lua_State *L, int cdim, struct vlasov_species_lw *species[GKYL_MAX_SPECIES])
{
  enum { TKEY = -2, TVAL = -1 };
  
  int curr = 0;
  lua_pushnil(L); // Initial key is nil.
  while (lua_next(L, TKEY) != 0) {
    // Key at TKEY and value at TVAL.
    if (lua_type(L, TVAL) == LUA_TUSERDATA) {
      struct vlasov_species_lw *vms = lua_touserdata(L, TVAL);

      if (vms->magic == VLASOV_SPECIES_DEFAULT) {
        if (vms->has_hamiltonian_func) {
          vms->hamiltonian_func_ref.ndim = cdim + vms->vdim;
        }

        if (vms->has_metric_func) {
          vms->metric_func_ref.ndim = cdim;
        }

        if (vms->has_inverse_metric_func) {
          vms->inverse_metric_func_ref.ndim = cdim;
        }

        if (vms->has_metric_determinant_func) {
          vms->metric_determinant_func_ref.ndim = cdim;
        }

        for (int i = 0; i < vms->num_init; i++) {
          if (vms->has_init_func[i]) {
            vms->init_func_ref[i].ndim = cdim + vms->vdim;
          }

          if (vms->has_density_init_func[i]) {
            vms->density_init_func_ref[i].ndim = cdim;
          }
          
          if (vms->has_V_drift_init_func[i]) {
            vms->V_drift_init_func_ref[i].ndim = cdim;
          }

          if (vms->has_temp_init_func[i]) {
            vms->temp_init_func_ref[i].ndim = cdim;
          }
        }

        if (vms->has_self_nu_func) {
          vms->self_nu_func_ref.ndim = cdim;
        }

        for (int i = 0; i < vms->num_sources; i++) {
          if (vms->source_has_init_func[i]) {
            vms->source_init_func_ref[i].ndim = cdim + vms->vdim;
          }

          if (vms->source_has_density_init_func[i]) {
            vms->source_density_init_func_ref[i].ndim = cdim;
          }

          if (vms->source_has_V_drift_init_func[i]) {
            vms->source_V_drift_init_func_ref[i].ndim = cdim;
          }

          if (vms->source_has_temp_init_func[i]) {
            vms->source_temp_init_func_ref[i].ndim = cdim;
          }
        }
        
        if (lua_type(L,TKEY) == LUA_TSTRING) {
          const char *key = lua_tolstring(L, TKEY, 0);
          strcpy(vms->vm_species.name, key);
        }
        species[curr++] = vms;
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
  const struct vlasov_species_lw *const *spa = a;
  const struct vlasov_species_lw *const *spb = b;
  return strcmp((*spa)->vm_species.name, (*spb)->vm_species.name);
}

// Gets all fluid species objects from the App table, which must on top of
// the stack. The number of fluid species is returned and the appropriate
// pointers set in the fluid species pointer array.
static int
get_fluid_species_inp(lua_State *L, int cdim, struct vlasov_fluid_species_lw *fluid_species[GKYL_MAX_SPECIES])
{
  enum { TKEY = -2, TVAL = -1 };
  
  int curr = 0;
  lua_pushnil(L); // Initial key is nil.
  while (lua_next(L, TKEY) != 0) {
    // Key at TKEY and value at TVAL.
    if (lua_type(L, TVAL) == LUA_TUSERDATA) {
      struct vlasov_fluid_species_lw *vmfs = lua_touserdata(L, TVAL);

      if (vmfs->magic == VLASOV_FLUID_SPECIES_DEFAULT) {
        
        vmfs->init_ctx.ndim = cdim;
        
        if (lua_type(L,TKEY) == LUA_TSTRING) {
          const char *key = lua_tolstring(L, TKEY, 0);
          strcpy(vmfs->vlasov_fluid_species.name, key);
        }
        fluid_species[curr++] = vmfs;
      }
    }
    lua_pop(L, 1);
  }

  return curr;
}

// Comparison method to sort fluid species array by fluid species name.
static int
fluid_species_compare_func(const void *a, const void *b)
{
  const struct vlasov_fluid_species_lw *const *spa = a;
  const struct vlasov_fluid_species_lw *const *spb = b;
  return strcmp((*spa)->vlasov_fluid_species.name, (*spb)->vlasov_fluid_species.name);
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
vm_parse_script_cli(struct gkyl_tool_args *acv)
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
vm_app_new(lua_State *L)
{
  struct vlasov_app_lw *app_lw = gkyl_malloc(sizeof(*app_lw));

  // The output prefix to use is stored in the global
  // GKYL_OUT_PREFIX. If this is not found then "g0-vlasov" is used.
  const char *sim_name = "g0-vlasov";

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
  app_lw->integrated_L2_f_calcs = glua_tbl_get_integer(L, "integratedL2fCalcs", INT_MAX);
  app_lw->integrated_mom_calcs = glua_tbl_get_integer(L, "integratedMomentCalcs", INT_MAX);
  app_lw->dt_failure_tol = glua_tbl_get_number(L, "dtFailureTol", 1.0e-4);
  app_lw->num_failures_max = glua_tbl_get_integer(L, "numFailuresMax", 20);

  struct gkyl_vm vm = { }; // Input table for app.

  strcpy(vm.name, sim_name);
  
  int cdim = 0;
  with_lua_tbl_tbl(L, "cells") {
    vm.cdim = cdim = glua_objlen(L);

    for (int d = 0; d < cdim; d++) {
      vm.cells[d] = glua_tbl_iget_integer(L, d + 1, 0);
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
      vm.lower[d] = glua_tbl_iget_number(L, d + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d = 0; d < cdim; d++) {
      vm.upper[d] = glua_tbl_iget_number(L, d + 1, 0);
    }
  }

  vm.cfl_frac = glua_tbl_get_number(L, "cflFrac", 0.95);
  vm.poly_order = glua_tbl_get_integer(L, "polyOrder", 1);

  vm.basis_type = get_basis_type(
    glua_tbl_get_string(L, "basis", "serendipity")
  );

  vm.num_periodic_dir = 0;
  if (glua_tbl_has_key(L, "periodicDirs")) {
    with_lua_tbl_tbl(L, "periodicDirs") {
      vm.num_periodic_dir = glua_objlen(L);

      for (int d = 0; d < vm.num_periodic_dir; d++) {
        // Indices are off by 1 between Lua and C.
        vm.periodic_dirs[d] = glua_tbl_iget_integer(L, d + 1, 0) - 1;
      }
    }
  }

  struct vlasov_species_lw *species[GKYL_MAX_SPECIES];

  // Set all species input.
  vm.num_species = get_species_inp(L, cdim, species);

  // Need to sort the species[] array by name of the species before
  // proceeding as there is no way to ensure that all cores loop over
  // Lua tables in the same order.
  qsort(species, vm.num_species, sizeof(struct vlasov_species_lw *), species_compare_func);  
  
  for (int s = 0; s < vm.num_species; s++) {
    vm.species[s] = species[s]->vm_species;
    vm.vdim = species[s]->vdim;

    app_lw->has_hamiltonian_func[s] = species[s]->has_hamiltonian_func;
    app_lw->hamiltonian_func_ctx[s] = species[s]->hamiltonian_func_ref;

    app_lw->has_metric_func[s] = species[s]->has_metric_func;
    app_lw->metric_func_ctx[s] = species[s]->metric_func_ref;

    app_lw->has_inverse_metric_func[s] = species[s]->has_inverse_metric_func;
    app_lw->inverse_metric_func_ctx[s] = species[s]->inverse_metric_func_ref;

    app_lw->has_metric_determinant_func[s] = species[s]->has_metric_determinant_func;
    app_lw->metric_determinant_func_ctx[s] = species[s]->metric_determinant_func_ref;

    if (species[s]->has_hamiltonian_func) {
      vm.species[s].hamil = gkyl_lw_eval_cb;
      vm.species[s].hamil_ctx = &app_lw->hamiltonian_func_ctx[s];
    }

    if (species[s]->has_metric_func) {
      vm.species[s].h_ij = gkyl_lw_eval_cb;
      vm.species[s].h_ij_ctx = &app_lw->metric_func_ctx[s];
    }

    if (species[s]->has_inverse_metric_func) {
      vm.species[s].h_ij_inv = gkyl_lw_eval_cb;
      vm.species[s].h_ij_inv_ctx = &app_lw->inverse_metric_func_ctx[s];
    }

    if (species[s]->has_metric_determinant_func) {
      vm.species[s].det_h = gkyl_lw_eval_cb;
      vm.species[s].det_h_ctx = &app_lw->metric_determinant_func_ctx[s];
    }

    app_lw->output_f_lte[s] = species[s]->output_f_lte;

    app_lw->num_init[s] = species[s]->num_init;
    for (int i = 0; i < app_lw->num_init[s]; i++) {
      app_lw->proj_id[s][i] = species[s]->proj_id[i];

      app_lw->has_init_func[s][i] = species[s]->has_init_func[i];
      app_lw->init_func_ctx[s][i] = species[s]->init_func_ref[i];

      app_lw->has_density_init_func[s][i] = species[s]->has_density_init_func[i];
      app_lw->density_init_func_ctx[s][i] = species[s]->density_init_func_ref[i];

      app_lw->has_V_drift_init_func[s][i] = species[s]->has_V_drift_init_func[i];
      app_lw->V_drift_init_func_ctx[s][i] = species[s]->V_drift_init_func_ref[i];
      
      app_lw->has_temp_init_func[s][i] = species[s]->has_temp_init_func[i];
      app_lw->temp_init_func_ctx[s][i] = species[s]->temp_init_func_ref[i];

      app_lw->correct_all_moms[s][i] = species[s]->correct_all_moms[i];
      app_lw->iter_eps[s][i] = species[s]->iter_eps[i];
      app_lw->max_iter[s][i] = species[s]->max_iter[i];
      app_lw->use_last_converged[s][i] = species[s]->use_last_converged[i];
    }

    vm.species[s].output_f_lte = app_lw->output_f_lte[s];

    vm.species[s].num_init = app_lw->num_init[s];
    for (int i = 0; i < app_lw->num_init[s]; i++) {
      vm.species[s].projection[i].proj_id = app_lw->proj_id[s][i];

      if (species[s]->has_init_func[i]) {
        vm.species[s].projection[i].func = gkyl_lw_eval_cb;
        vm.species[s].projection[i].ctx_func = &app_lw->init_func_ctx[s][i];
      }

      if (species[s]->has_density_init_func[i]) {
        vm.species[s].projection[i].density = gkyl_lw_eval_cb;
        vm.species[s].projection[i].ctx_density = &app_lw->density_init_func_ctx[s][i];
      }

      if (species[s]->has_V_drift_init_func[i]) {
        vm.species[s].projection[i].V_drift = gkyl_lw_eval_cb;
        vm.species[s].projection[i].ctx_V_drift = &app_lw->V_drift_init_func_ctx[s][i];
      }

      if (species[s]->has_temp_init_func[i]) {
        vm.species[s].projection[i].temp = gkyl_lw_eval_cb;
        vm.species[s].projection[i].ctx_temp = &app_lw->temp_init_func_ctx[s][i];
      }

      vm.species[s].projection[i].correct_all_moms = app_lw->correct_all_moms[s][i];
      vm.species[s].projection[i].iter_eps = app_lw->iter_eps[s][i];
      vm.species[s].projection[i].max_iter = app_lw->max_iter[s][i];
      vm.species[s].projection[i].use_last_converged = app_lw->use_last_converged[s][i];
    }

    app_lw->collision_id[s] = species[s]->collision_id;

    app_lw->has_self_nu_func[s] = species[s]->has_self_nu_func;
    app_lw->self_nu_func_ctx[s] = species[s]->self_nu_func_ref;

    app_lw->num_cross_collisions[s] = species[s]->num_cross_collisions;
    for (int i = 0; i < app_lw->num_cross_collisions[s]; i++) {
      strcpy(app_lw->collide_with[s][i], species[s]->collide_with[i]);
    }

    app_lw->collision_correct_all_moms[s] = species[s]->collision_correct_all_moms;
    app_lw->collision_iter_eps[s] = species[s]->collision_iter_eps;
    app_lw->collision_max_iter[s] = species[s]->collision_max_iter;
    app_lw->fixed_temp_relax[s] = species[s]->fixed_temp_relax;
    app_lw->collision_use_last_converged[s] = species[s]->collision_use_last_converged;
    app_lw->has_implicit_coll_scheme[s] = species[s]->has_implicit_coll_scheme;

    vm.species[s].collisions.collision_id = app_lw->collision_id[s];

    if (species[s]->has_self_nu_func) {
      vm.species[s].collisions.self_nu = gkyl_lw_eval_cb;
      vm.species[s].collisions.ctx = &app_lw->self_nu_func_ctx[s];
    }

    vm.species[s].collisions.num_cross_collisions = app_lw->num_cross_collisions[s];
    for (int i = 0; i < app_lw->num_cross_collisions[s]; i++) {
      strcpy(vm.species[s].collisions.collide_with[i], app_lw->collide_with[s][i]);
    }

    vm.species[s].collisions.correct_all_moms = app_lw->collision_correct_all_moms[s];
    vm.species[s].collisions.iter_eps = app_lw->collision_iter_eps[s];
    vm.species[s].collisions.max_iter = app_lw->collision_max_iter[s];
    vm.species[s].collisions.fixed_temp_relax = app_lw->fixed_temp_relax[s];
    vm.species[s].collisions.use_last_converged = app_lw->collision_use_last_converged[s];
    vm.species[s].collisions.has_implicit_coll_scheme = app_lw->has_implicit_coll_scheme[s];

    app_lw->source_id[s] = species[s]->source_id;

    app_lw->source_length[s] = species[s]->source_length;
    strcpy(app_lw->source_species[s], species[s]->source_species);

    app_lw->num_sources[s] = species[s]->num_sources;
    for (int i = 0; i < app_lw->num_sources[s]; i++) {
      app_lw->source_proj_id[s][i] = species[s]->source_proj_id[i];

      app_lw->source_has_init_func[s][i] = species[s]->source_has_init_func[i];
      app_lw->source_init_func_ctx[s][i] = species[s]->source_init_func_ref[i];

      app_lw->source_has_density_init_func[s][i] = species[s]->source_has_density_init_func[i];
      app_lw->source_density_init_func_ctx[s][i] = species[s]->source_density_init_func_ref[i];

      app_lw->source_has_V_drift_init_func[s][i] = species[s]->source_has_V_drift_init_func[i];
      app_lw->source_V_drift_init_func_ctx[s][i] = species[s]->source_V_drift_init_func_ref[i];

      app_lw->source_has_temp_init_func[s][i] = species[s]->source_has_temp_init_func[i];
      app_lw->source_temp_init_func_ctx[s][i] = species[s]->source_temp_init_func_ref[i];

      app_lw->source_correct_all_moms[s][i] = species[s]->source_correct_all_moms[i];
      app_lw->source_iter_eps[s][i] = species[s]->source_iter_eps[i];
      app_lw->source_max_iter[s][i] = species[s]->source_max_iter[i];
      app_lw->source_use_last_converged[s][i] = species[s]->source_use_last_converged[i];
    }

    vm.species[s].source.source_id = app_lw->source_id[s];

    vm.species[s].source.source_length = app_lw->source_length[s];
    strcpy(vm.species[s].source.source_species, app_lw->source_species[s]);

    vm.species[s].source.num_sources = app_lw->num_sources[s];
    for (int i = 0; i < app_lw->num_sources[s]; i++) {
      vm.species[s].source.projection[i].proj_id = app_lw->source_proj_id[s][i];

      if (species[s]->source_has_init_func[i]) {
        vm.species[s].source.projection[i].func = gkyl_lw_eval_cb;
        vm.species[s].source.projection[i].ctx_func = &app_lw->source_init_func_ctx[s][i];
      }

      if (species[s]->source_has_density_init_func[i]) {
        vm.species[s].source.projection[i].density = gkyl_lw_eval_cb;
        vm.species[s].source.projection[i].ctx_density = &app_lw->source_density_init_func_ctx[s][i];
      }

      if (species[s]->source_has_V_drift_init_func[i]) {
        vm.species[s].source.projection[i].V_drift = gkyl_lw_eval_cb;
        vm.species[s].source.projection[i].ctx_V_drift = &app_lw->source_V_drift_init_func_ctx[s][i];
      }

      if (species[s]->source_has_temp_init_func[i]) {
        vm.species[s].source.projection[i].temp = gkyl_lw_eval_cb;
        vm.species[s].source.projection[i].ctx_temp = &app_lw->source_temp_init_func_ctx[s][i];
      }

      vm.species[s].source.projection[i].correct_all_moms = app_lw->source_correct_all_moms[s][i];
      vm.species[s].source.projection[i].iter_eps = app_lw->source_iter_eps[s][i];
      vm.species[s].source.projection[i].max_iter = app_lw->source_max_iter[s][i];
      vm.species[s].source.projection[i].use_last_converged = app_lw->source_use_last_converged[s][i];
    }
  }

  struct vlasov_fluid_species_lw *fluid_species[GKYL_MAX_SPECIES];

  // Set all fluid species input.
  vm.num_fluid_species = get_fluid_species_inp(L, cdim, fluid_species);

  // Need to sort the fluid_species[] array by name of the fluid species before
  // proceeding as there is no way to ensure that all cores loop over
  // Lua tables in the same order.
  qsort(fluid_species, vm.num_fluid_species, sizeof(struct vlasov_fluid_species_lw *), fluid_species_compare_func);
  
  for (int s = 0; s < vm.num_fluid_species; s++) {
    vm.fluid_species[s] = fluid_species[s]->vlasov_fluid_species;

    fluid_species[s]->init_ctx.ndim = cdim; 
    app_lw->fluid_species_init_ctx[s] = fluid_species[s]->init_ctx;
    vm.fluid_species[s].init = gkyl_lw_eval_cb;
    vm.fluid_species[s].ctx = &app_lw->fluid_species_init_ctx[s];

    if (fluid_species[s]->has_app_advect_func) {
      fluid_species[s]->app_advect_func_ref.ndim = cdim; 
      app_lw->app_advect_func_ctx[s] = fluid_species[s]->app_advect_func_ref;
      vm.fluid_species[s].advection.velocity = gkyl_lw_eval_cb;
      vm.fluid_species[s].advection.velocity_ctx = &app_lw->app_advect_func_ctx[s];
    }

    if (fluid_species[s]->has_n0_func) {
      fluid_species[s]->n0_func_ref.ndim = cdim; 
      app_lw->n0_func_ctx[s] = fluid_species[s]->n0_func_ref;
      vm.fluid_species[s].can_pb_n0 = gkyl_lw_eval_cb;
      vm.fluid_species[s].can_pb_n0_ctx = &app_lw->n0_func_ctx[s];
    }

    if (fluid_species[s]->has_diffusion_func) {
      fluid_species[s]->diffusion_func_ref.ndim = cdim;
      if (cdim == 2) {
        fluid_species[s]->diffusion_func_ref.nret = 3;
      }
      else {
        fluid_species[s]->diffusion_func_ref.nret = 6;
      }
      app_lw->diffusion_func_ctx[s] = fluid_species[s]->diffusion_func_ref;
      vm.fluid_species[s].diffusion.Dij = gkyl_lw_eval_cb;
      vm.fluid_species[s].diffusion.Dij_ctx = &app_lw->diffusion_func_ctx[s];
    }
  }

  // Set field input.
  vm.skip_field = glua_tbl_get_bool(L, "skipField", false);
  vm.is_electrostatic = glua_tbl_get_bool(L, "isElectrostatic", false);

  with_lua_tbl_key(L, "field") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct vlasov_field_lw *vmf = lua_touserdata(L, -1);

      if (vmf->magic == VLASOV_FIELD_DEFAULT) {
        vmf->init_ref.ndim = cdim;

        vm.field = vmf->vm_field;

        app_lw->field_func_ctx = vmf->init_ref;
        vm.field.init = gkyl_lw_eval_cb;
        vm.field.ctx = &app_lw->field_func_ctx;

        if (vmf->has_external_potential_func) {
          vmf->external_potential_func_ref.ndim = cdim;

          app_lw->external_potential_func_ctx = vmf->external_potential_func_ref;
          vm.field.external_potentials = gkyl_lw_eval_cb;
          vm.field.external_potentials_ctx = &app_lw->external_potential_func_ctx;

          vm.field.external_potentials_evolve = vmf->evolve_external_potential;
        }

        if (vmf->has_external_field_func) {
          vmf->external_field_func_ref.ndim = cdim;

          app_lw->external_field_func_ctx = vmf->external_field_func_ref;
          vm.field.ext_em = gkyl_lw_eval_cb;
          vm.field.ext_em_ctx = &app_lw->external_field_func_ctx;

          vm.field.ext_em_evolve = vmf->evolve_external_field;
        }

        if (vmf->has_applied_current_func) {
          vmf->applied_current_func_ref.ndim = cdim;

          app_lw->applied_current_func_ctx = vmf->applied_current_func_ref;
          vm.field.app_current = gkyl_lw_eval_cb;
          vm.field.app_current_ctx = &app_lw->applied_current_func_ctx;

          vm.field.app_current_evolve = vmf->evolve_applied_current;
        }
      }
    }
  }

  // Create parallelism.
  struct gkyl_comm *comm = 0;

  for (int d = 0; d < cdim; d++) {
    vm.parallelism.cuts[d] = cuts[d]; 
  }

  struct gkyl_tool_args *args = gkyl_tool_args_new(L);
  struct script_cli script_cli = vm_parse_script_cli(args);

  script_cli.use_mpi = false;
  with_lua_global(L, "GKYL_MPI_COMM") {
    if (lua_islightuserdata(L, -1)) {
      script_cli.use_mpi = true;
    }
  }

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

  if (comm == 0)
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = script_cli.use_gpu,
      }
    );    

  vm.parallelism.comm = comm;
  vm.parallelism.use_gpu = script_cli.use_gpu;

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
  
  app_lw->app = gkyl_vlasov_app_new(&vm);

  gkyl_comm_release(comm);

  // Create Lua userdata.
  struct vlasov_app_lw **l_app_lw = lua_newuserdata(L, sizeof(struct vlasov_app_lw*));
  *l_app_lw = app_lw; // Point it to the Lua app pointer.

  // Set metatable.
  luaL_getmetatable(L, VLASOV_APP_METATABLE_NM);
  lua_setmetatable(L, -2);

  gkyl_tool_args_release(script_cli.rest);
  gkyl_tool_args_release(args);
  
  return 1;
}

// Apply initial conditions. (time) -> bool.
static int
vm_app_apply_ic(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->t_start);
  gkyl_vlasov_app_apply_ic(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to field. (time) -> bool.
static int
vm_app_apply_ic_field(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->t_start);
  gkyl_vlasov_app_apply_ic_field(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to species. (sidx, time) -> bool.
static int
vm_app_apply_ic_species(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double t0 = luaL_optnumber(L, 3, app_lw->t_start);
  gkyl_vlasov_app_apply_ic_species(app_lw->app, sidx, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute diagnostic moments. () -> bool.
static int
vm_app_calc_mom(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  gkyl_vlasov_app_calc_mom(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated moments. (tm) -> bool.
static int
vm_app_calc_integrated_mom(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_vlasov_app_calc_integrated_mom(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated L2 norm of distribution function. (tm) -> bool.
static int
vm_app_calc_integrated_L2_f(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_vlasov_app_calc_integrated_L2_f(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated field energy (L2 norm of each field
// component). (tm) -> bool.
static int
vm_app_calc_field_energy(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_vlasov_app_calc_field_energy(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Write solution (field and species) to file (time, frame) -> bool.
static int
vm_app_write(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_vlasov_app_write(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write field to file (time, frame) -> bool.
static int
vm_app_write_field(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_vlasov_app_write_field(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write species solution to file (sidx, time, frame) -> bool.
static int
vm_app_write_species(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double tm = luaL_checknumber(L, 3);
  int frame = luaL_checkinteger(L, 4);
  gkyl_vlasov_app_write_species(app_lw->app, sidx, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write diagnostic moments to file (time, frame) -> bool.
static int
vm_app_write_mom(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_vlasov_app_write_mom(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated moments to file () -> bool.
static int
vm_app_write_integrated_mom(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  gkyl_vlasov_app_write_integrated_mom(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated fluid moments to file () -> bool.
static int
vm_app_write_fluid_integrated_mom(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  gkyl_vlasov_app_write_fluid_integrated_mom(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated L2 norm of f to file () -> bool.
static int
vm_app_write_integrated_L2_f(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  gkyl_vlasov_app_write_integrated_L2_f(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated field energy to file () -> bool.
static int
vm_app_write_field_energy(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  gkyl_vlasov_app_write_field_energy(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write simulation statistics to JSON. () -> bool.
static int
vm_app_stat_write(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  gkyl_vlasov_app_stat_write(app_lw->app);

  lua_pushboolean(L, status);
  return 1;  
}

// Write data from simulation to file.
static void
write_data(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_vlasov_app_write(app, t_curr, frame);
    gkyl_vlasov_app_write_field_energy(app);
    gkyl_vlasov_app_write_integrated_mom(app);
    gkyl_vlasov_app_write_fluid_integrated_mom(app);
    gkyl_vlasov_app_write_integrated_L2_f(app);

    gkyl_vlasov_app_calc_mom(app);
    gkyl_vlasov_app_write_mom(app, t_curr, frame);
  }
}

// Calculate and append field energy to dynvector.
static void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_field_energy(app, t_curr);
  }
}

// Calculate and append integrated moments to dynvector.
static void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  }
}

// Calculate and append integrated L2 norm of distribution function to dynvector.
static void
calc_integrated_L2_f(struct gkyl_tm_trigger* l2t, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(l2t, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_integrated_L2_f(app, t_curr);
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
write_step_message(const struct gkyl_vlasov_app *app, struct step_message_trigs *trigs, int step, double t_curr, double dt_next)
{
  if (gkyl_tm_trigger_check_and_bump(&trigs->log_trig, t_curr)) {
    if (trigs->log_count > 0) {
      gkyl_vlasov_app_cout(app, stdout, " Step %6d at time %#11.8g.  Time-step  %.6e.  Completed %g%s\n", step, t_curr, dt_next, trigs->tenth * 10.0, "%");
    }
    else {
      trigs->log_count += 1;
    }
    
    trigs->tenth += 1;
  }
  if (gkyl_tm_trigger_check_and_bump(&trigs->log_trig_1p, t_curr)) {
    gkyl_vlasov_app_cout(app, stdout, "%d", trigs->p1c);
    trigs->p1c = (trigs->p1c+1) % 10;
  }
}

static void
show_help(const struct gkyl_vlasov_app *app)
{
  gkyl_vlasov_app_cout(app, stdout, "Vlasov script takes the following arguments:\n");
  gkyl_vlasov_app_cout(app, stdout, " -h   Print this help message and exit\n");
  gkyl_vlasov_app_cout(app, stdout, " -sN  Only run N steps of simulation\n");
  gkyl_vlasov_app_cout(app, stdout, " -S   Do not initialize MPI\n");
  gkyl_vlasov_app_cout(app, stdout, " -G   Do not initialize CUDA\n");
  gkyl_vlasov_app_cout(app, stdout, " -m   Run memory tracer\n");
  gkyl_vlasov_app_cout(app, stdout, " -V   Show verbose output\n");
  gkyl_vlasov_app_cout(app, stdout, " -rN  Restart simulation from frame N\n");

  gkyl_vlasov_app_cout(app, stdout, "\n");
}

// Run simulation. (num_steps) -> bool. num_steps is optional.
static int
vm_app_run(lua_State *L)
{
  bool ret_status = true;

  // Create app object.
  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;
  struct gkyl_vlasov_app *app = app_lw->app;

  // Parse command lines arguments passed to input file.
  struct gkyl_tool_args *args = gkyl_tool_args_new(L);  
  
  struct script_cli script_cli = vm_parse_script_cli(args);
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

  gkyl_vlasov_app_cout(app, stdout, "Initializing Vlasov Simulation ...\n");

  if (script_cli.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  // Initialize simulation.
  bool is_restart = script_cli.is_restart;
  int restart_frame = script_cli.restart_frame;

  int frame_curr = 0;
  if (is_restart) {
    struct gkyl_app_restart_status status = gkyl_vlasov_app_read_from_frame(app, restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_vlasov_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_vlasov_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_vlasov_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_vlasov_app_apply_ic(app, t_curr);
  }

  int num_frames = app_lw->num_frames;
  int field_energy_calcs = app_lw->field_energy_calcs;
  int integrated_mom_calcs = app_lw->integrated_mom_calcs;
  int integrated_L2_f_calcs = app_lw->integrated_L2_f_calcs;
  // Triggers for IO and logging.
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames, .tcurr = frame_curr * (t_end / num_frames), .curr = frame_curr };
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger l2f_trig = { .dt = t_end / integrated_L2_f_calcs, .tcurr = t_curr, .curr = frame_curr };

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
  calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);

  gkyl_vlasov_app_cout(app, stdout, "Initialization completed in %g sec\n\n", gkyl_time_diff_now_sec(tm_ic0));

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = app_lw->dt_failure_tol;
  int num_failures = 0, num_failures_max = app_lw->num_failures_max;

  bool use_verbose = script_cli.use_verbose;

  long step = 1;
  while ((t_curr < t_end) && (step <= num_steps)) {
    if (use_verbose) {
      gkyl_vlasov_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    }
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    if (use_verbose) {
      gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    }

    if (!status.success) {
      gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_vlasov_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_vlasov_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_vlasov_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_vlasov_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_vlasov_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        calc_integrated_L2_f(&l2f_trig, app, t_curr, true);
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
  calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);
  gkyl_vlasov_app_stat_write(app);

  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  gkyl_vlasov_app_cout(app, stdout, "\n");
  gkyl_vlasov_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_vlasov_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_vlasov_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_vlasov_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_vlasov_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_vlasov_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

freeresources:

  lua_pushboolean(L, ret_status);
  return 1;
}

// Clean up memory allocated for simulation.
static int
vm_app_gc(lua_State *L)
{
  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  gkyl_vlasov_app_release(app_lw->app);
  gkyl_free(*l_app_lw);
  
  return 0;
}

// App constructor.
static struct luaL_Reg vm_app_ctor[] = {
  { "new",  vm_app_new },
  { 0, 0 }
};

// App methods.
static struct luaL_Reg vm_app_funcs[] = {
  { "apply_ic", vm_app_apply_ic },
  { "apply_ic_field", vm_app_apply_ic_field },
  { "apply_ic_species", vm_app_apply_ic_species },
  { "calc_mom", vm_app_calc_mom },
  { "calc_integrated_mom", vm_app_calc_integrated_mom },
  { "calc_integrated_L2_f", vm_app_calc_integrated_L2_f },
  { "calc_field_energy", vm_app_calc_field_energy },
  { "write", vm_app_write },
  { "write_field", vm_app_write_field },
  { "write_species", vm_app_write_species },
  { "write_mom", vm_app_write_mom },
  { "write_integrated_mom", vm_app_write_integrated_mom },
  { "write_integrated_L2_f", vm_app_write_integrated_L2_f },
  { "write_field_energy", vm_app_write_field_energy },
  { "stat_write", vm_app_stat_write },
  { "run", vm_app_run },
  { 0, 0 }
};

static void
app_openlibs(lua_State *L)
{
  // Register top-level App.
  do {
    luaL_newmetatable(L, VLASOV_APP_METATABLE_NM);

    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, vm_app_gc);
    lua_settable(L, -3);

    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, vm_app_funcs);
    
    luaL_register(L, "G0.Vlasov.App", vm_app_ctor);
    
  }
  while (0);

  // Register Species input struct.
  do {
    luaL_newmetatable(L, VLASOV_SPECIES_METATABLE_NM);
    luaL_register(L, "G0.Vlasov.Species", vm_species_ctor);
  }
  while (0);

  // Register Fluid Species input struct.
  do {
    luaL_newmetatable(L, VLASOV_FLUID_SPECIES_METATABLE_NM);
    luaL_register(L, "G0.Vlasov.FluidSpecies", vm_fluid_species_ctor);
  }
  while (0);

  // Register Field input struct.
  do {
    luaL_newmetatable(L, VLASOV_FIELD_METATABLE_NM);
    luaL_register(L, "G0.Vlasov.Field", vm_field_ctor);
  }
  while (0);

  // Add globals and other parameters.
  do {
    // Table to store references so the objects do not get garbage-collected.
    lua_newtable(L);
    lua_setglobal(L, "__vlasov_ref_table");
  }
  while (0);
}

void
gkyl_vlasov_lw_openlibs(lua_State *L)
{
  gkyl_register_poisson_bc_types(L);

  // Register types for Vlasov projection, model ID, collision ID, and source ID initialization.
  gkyl_register_vlasov_projection_types(L);
  gkyl_register_vlasov_model_types(L);
  gkyl_register_vlasov_collision_types(L);
  gkyl_register_vlasov_source_types(L);  
  
  eqn_openlibs(L);
  app_openlibs(L);
}

#endif
