#ifdef GKYL_HAVE_LUA

#include <gkyl_alloc.h>
#include <gkyl_eqn_type.h>
#include <gkyl_lua_utils.h>
#include <gkyl_lw_priv.h>
#include <gkyl_null_comm.h>
#include <gkyl_pkpm.h>
#include <gkyl_pkpm_lw.h>
#include <gkyl_pkpm_priv.h>
#include <gkyl_zero_lw.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
#include <limits.h>

#include <string.h>

#include <stc/coption.h>
#include <kann.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

// Magic IDs for use in distinguishing various species and field types.
enum pkpm_magic_ids {
  PKPM_SPECIES_DEFAULT = 100, // Non-relativistic PKPM species.
  PKPM_FIELD_DEFAULT, // Maxwell equations.
};

/* *************** */
/* Species methods */
/* *************** */

// Metatable name for species input struct.
#define PKPM_SPECIES_METATABLE_NM "GkeyllZero.App.PKPM.Species"

// Lua userdata object for constructing species input.
struct pkpm_species_lw {
  int magic; // This must be first element in the struct.
  
  struct gkyl_pkpm_species pkpm_species; // Input struct to construct species.
  int vdim; // Velocity space dimensions.
  bool evolve; // Is this species evolved?

  bool has_dist_init_func; // Is there a distribution initialization function?
  struct lua_func_ctx dist_init_func_ref; // Lua registry reference to distribution initialization function.

  bool has_fluid_init_func; // Is there a fluid initialization function?
  struct lua_func_ctx fluid_init_func_ref; // Lua registry reference to fluid initialization function.

  bool has_applied_acceleration_func; // Is there an applied acceleration initialization function?
  struct lua_func_ctx applied_acceleration_func_ref; // Lua registry reference to applied acceleration initialization function.
  bool evolve_applied_acceleration; // Is the applied acceleration evolved?

  bool has_diffusion; // Is there a diffusion operator? 
  double D; // Diffusion coefficient. 
  double order; // Diffusion order (e.g., grad^2, grad^4, grad^6). 

  enum gkyl_collision_id collision_id; // Collision type.
  
  bool has_self_nu_func; // Is there a self-collision frequency function?
  struct lua_func_ctx self_nu_func_ref; // Lua registry reference to self-collision frequency function.

  int num_cross_collisions; // Number of species that we cross-collide with.
  char collide_with[GKYL_MAX_SPECIES][128]; // Names of species that we cross-collide with.
};

static int
pkpm_species_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_pkpm_species pkpm_species = { };

  pkpm_species.model_id = GKYL_MODEL_DEFAULT;
  
  pkpm_species.charge = glua_tbl_get_number(L, "charge", 0.0);
  pkpm_species.mass = glua_tbl_get_number(L, "mass", 1.0);

  with_lua_tbl_tbl(L, "cells") {
    vdim = glua_objlen(L);

    for (int d = 0; d < vdim; d++) {
      pkpm_species.cells[d] = glua_tbl_iget_integer(L, d + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "lower") {
    for (int d = 0; d < vdim; d++) {
      pkpm_species.lower[d] = glua_tbl_iget_number(L, d + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d = 0; d < vdim; d++) {
      pkpm_species.upper[d] = glua_tbl_iget_number(L, d + 1, 0);
    }
  }

  bool evolve = glua_tbl_get_bool(L, "evolve", true);

  with_lua_tbl_tbl(L, "bcx") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      pkpm_species.bcx[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcy") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      pkpm_species.bcy[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcz") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      pkpm_species.bcz[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  bool has_dist_init_func = false;
  int dist_init_func_ref = LUA_NOREF;

  bool has_fluid_init_func = false;
  int fluid_init_func_ref = LUA_NOREF;

  if (glua_tbl_get_func(L, "initDist")) {
    dist_init_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_dist_init_func = true;
  }

  if (glua_tbl_get_func(L, "initFluid")) {
    fluid_init_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_fluid_init_func = true;
  }

  bool has_applied_acceleration_func = false;
  int applied_acceleration_func_ref = LUA_NOREF;
  bool evolve_applied_acceleration = false;

  if (glua_tbl_get_func(L, "appliedAcceleration")) {
    applied_acceleration_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_applied_acceleration_func = true;

    evolve_applied_acceleration = glua_tbl_get_bool(L, "evolveAppliedAcceleration", false);
  }

  bool has_diffusion = false; 
  double D = 0.0; 
  int order = 0; 
  with_lua_tbl_tbl(L, "diffusion") {
    has_diffusion = true; 
    D = glua_tbl_get_number(L, "D", 0.0);
    order = glua_tbl_get_integer(L, "order", 0);
  }

  enum gkyl_collision_id collision_id = GKYL_NO_COLLISIONS;

  bool has_self_nu_func = false;
  int self_nu_func_ref = LUA_NOREF;

  int num_cross_collisions = 0;
  char collide_with[GKYL_MAX_SPECIES][128];

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
  }
  
  struct pkpm_species_lw *pkpm_s_lw = lua_newuserdata(L, sizeof(*pkpm_s_lw));
  pkpm_s_lw->magic = PKPM_SPECIES_DEFAULT;
  pkpm_s_lw->vdim = vdim;
  pkpm_s_lw->evolve = evolve;
  pkpm_s_lw->pkpm_species = pkpm_species;

  pkpm_s_lw->has_dist_init_func = has_dist_init_func;
  pkpm_s_lw->dist_init_func_ref = (struct lua_func_ctx) {
    .func_ref = dist_init_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 2,
    .L = L,
  };

  pkpm_s_lw->has_fluid_init_func = has_fluid_init_func;
  pkpm_s_lw->fluid_init_func_ref = (struct lua_func_ctx) {
    .func_ref = fluid_init_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 3,
    .L = L,
  };

  pkpm_s_lw->has_applied_acceleration_func = has_applied_acceleration_func;
  pkpm_s_lw->applied_acceleration_func_ref = (struct lua_func_ctx) {
    .func_ref = applied_acceleration_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 3,
    .L = L,
  };
  pkpm_s_lw->evolve_applied_acceleration = evolve_applied_acceleration;

  pkpm_s_lw->has_diffusion = has_diffusion; 
  pkpm_s_lw->D = D; 
  pkpm_s_lw->order = order; 

  pkpm_s_lw->collision_id = collision_id;

  pkpm_s_lw->has_self_nu_func = has_self_nu_func;
  pkpm_s_lw->self_nu_func_ref = (struct lua_func_ctx) {
    .func_ref = self_nu_func_ref,
    .ndim = 0,
    .nret = 1,
    .L = L,
  };

  pkpm_s_lw->num_cross_collisions = num_cross_collisions;
  for (int i = 0; i < num_cross_collisions; i++) {
    strcpy(pkpm_s_lw->collide_with[i], collide_with[i]);
  }
  
  // Set metatable.
  luaL_getmetatable(L, PKPM_SPECIES_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor.
static struct luaL_Reg pkpm_species_ctor[] = {
  { "new", pkpm_species_lw_new },
  { 0, 0 }
};

/* ************* */
/* Field methods */
/* ************* */

// Metatable name for field input struct.
#define PKPM_FIELD_METATABLE_NM "GkeyllZero.App.PKPM.Field"

// Lua userdata object for constructing field input.
struct pkpm_field_lw {
  int magic; // This must be first element in the struct.
  
  struct gkyl_pkpm_field pkpm_field; // Input struct to construct field.
  bool evolve; // Is this field evolved?
  struct lua_func_ctx init_ref; // Lua registry reference to initilization function.

  bool has_external_field_func; // Is there an external field initialization function?
  struct lua_func_ctx external_field_func_ref; // Lua registry reference to external field initialization function.
  bool evolve_external_field; // Is the external field evolved?

  bool has_applied_current_func; // Is there an applied current initialization function?
  struct lua_func_ctx applied_current_func_ref; // Lua registry reference to applied current initialization function.
  bool evolve_applied_current; // Is the applied current evolved?
};

static int
pkpm_field_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_pkpm_field pkpm_field = { };

  pkpm_field.field_id = GKYL_FIELD_E_B;  
  
  pkpm_field.epsilon0 = glua_tbl_get_number(L, "epsilon0", 1.0);
  pkpm_field.mu0 = glua_tbl_get_number(L, "mu0", 1.0);
  pkpm_field.elcErrorSpeedFactor = glua_tbl_get_number(L, "elcErrorSpeedFactor", 0.0);
  pkpm_field.mgnErrorSpeedFactor = glua_tbl_get_number(L, "mgnErrorSpeedFactor", 0.0);

  bool evolve = glua_tbl_get_bool(L, "evolve", true);
  pkpm_field.is_static = !evolve;

  int init_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "init")) {
    init_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }

  with_lua_tbl_tbl(L, "bcx") { 
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      pkpm_field.bcx[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcy") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      pkpm_field.bcy[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcz") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      pkpm_field.bcz[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
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

  struct pkpm_field_lw *pkpm_f_lw = lua_newuserdata(L, sizeof(*pkpm_f_lw));

  pkpm_f_lw->magic = PKPM_FIELD_DEFAULT;
  pkpm_f_lw->evolve = evolve;
  pkpm_f_lw->pkpm_field = pkpm_field;
  
  pkpm_f_lw->init_ref = (struct lua_func_ctx) {
    .func_ref = init_ref,
    .ndim = 0, // This will be set later.
    .nret = 6,
    .L = L,
  };

  pkpm_f_lw->has_external_field_func = has_external_field_func;
  pkpm_f_lw->external_field_func_ref = (struct lua_func_ctx) {
    .func_ref = external_field_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 6,
    .L = L,
  };
  pkpm_f_lw->evolve_external_field = evolve_external_field;

  pkpm_f_lw->has_applied_current_func = has_applied_current_func;
  pkpm_f_lw->applied_current_func_ref = (struct lua_func_ctx) {
    .func_ref = applied_current_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 3,
    .L = L,
  };
  pkpm_f_lw->evolve_applied_current = evolve_applied_current;
  
  // Set metatable.
  luaL_getmetatable(L, PKPM_FIELD_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor.
static struct luaL_Reg pkpm_field_ctor[] = {
  { "new", pkpm_field_lw_new },
  { 0, 0 }
};

/* *********** */
/* App methods */
/* *********** */

// Metatable name for top-level PKPM App.
#define PKPM_APP_METATABLE_NM "GkeyllZero.App.PKPM"

// Lua userdata object for holding PKPM app and run parameters.
struct pkpm_app_lw {
  gkyl_pkpm_app *app; // PKPM app object.

  struct lua_func_ctx dist_init_func_ctx[GKYL_MAX_SPECIES]; // Context for distribution initialization function.
  struct lua_func_ctx fluid_init_func_ctx[GKYL_MAX_SPECIES]; // Context for fluid initialization function.
  struct lua_func_ctx applied_acceleration_func_ctx[GKYL_MAX_SPECIES]; // Context for applied acceleration function.

  enum gkyl_collision_id collision_id[GKYL_MAX_SPECIES]; // Collision type.

  struct lua_func_ctx self_nu_func_ctx[GKYL_MAX_SPECIES]; // Context for self-collision frequency function.

  int num_cross_collisions[GKYL_MAX_SPECIES]; // Number of species that we cross-collide with.
  char collide_with[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES][128]; // Names of species that we cross-collide with.

  struct lua_func_ctx field_func_ctx; // Function context for field.
  struct lua_func_ctx external_field_func_ctx; // Function context for external field.
  struct lua_func_ctx applied_current_func_ctx; // Function context for applied current.
  
  double t_start, t_end; // Start and end times of simulation.
  int num_frames; // Number of data frames to write.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  // Training parameters.
  bool train_nn; // Train neural network on simulation data?
  bool train_ab_initio; // Train neural network ab initio?
  int nn_width; // Number of neurons to use per layer.
  int nn_depth; // Number of layers to use.
  char train_nn_file[256]; // File path of neural network to train.
  int num_trains; // Number of times to train neural network.
  int num_nn_writes; // Number of times to write out neural network.
  int num_input_moms; // Number of "input" moments to train on.
  int* input_moms; // Array of "input" moments to train on.
  int num_output_moms; // Number of "output" moments to train on.
  int* output_moms; // Array of "output" moments to train on.
  bool test_nn; // Test neural network on simulation data?
  char test_nn_file[256]; // File path of neural network to test.
  int num_tests; // Number of times to test neural network.
};

// Gets all species objects from the App table, which must on top of
// the stack. The number of species is returned and the appropriate
// pointers set in the species pointer array.
static int
get_species_inp(lua_State *L, int cdim, struct pkpm_species_lw *species[GKYL_MAX_SPECIES])
{
  enum { TKEY = -2, TVAL = -1 };
  
  int curr = 0;
  lua_pushnil(L); // Initial key is nil.
  while (lua_next(L, TKEY) != 0) {
    // Key at TKEY and value at TVAL.
    if (lua_type(L, TVAL) == LUA_TUSERDATA) {
      struct pkpm_species_lw *pkpm_s = lua_touserdata(L, TVAL);

      if (pkpm_s->magic == PKPM_SPECIES_DEFAULT) {
        if (pkpm_s->has_dist_init_func) {
          pkpm_s->dist_init_func_ref.ndim = cdim + pkpm_s->vdim;
        }

        if(pkpm_s->has_fluid_init_func) {
          pkpm_s->fluid_init_func_ref.ndim = cdim;
        }

        if(pkpm_s->has_applied_acceleration_func) {
          pkpm_s->applied_acceleration_func_ref.ndim = cdim;
        }

        if (pkpm_s->has_self_nu_func) {
          pkpm_s->self_nu_func_ref.ndim = cdim;
        }
        
        if (lua_type(L, TKEY) == LUA_TSTRING) {
          const char *key = lua_tolstring(L, TKEY, 0);
          strcpy(pkpm_s->pkpm_species.name, key);
        }
        species[curr++] = pkpm_s;
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
  const struct pkpm_species_lw *const *spa = a;
  const struct pkpm_species_lw *const *spb = b;
  return strcmp((*spa)->pkpm_species.name, (*spb)->pkpm_species.name);
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
pkpm_parse_script_cli(struct gkyl_tool_args *acv)
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
pkpm_app_new(lua_State *L)
{
  struct pkpm_app_lw *app_lw = gkyl_malloc(sizeof(*app_lw));

  // The output prefix to use is stored in the global
  // GKYL_OUT_PREFIX. If this is not found then "g0-pkpm" is used.
  const char *sim_name = "g0-pkpm";

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

  // Initialize training parameters.
  app_lw->train_nn = glua_tbl_get_bool(L, "trainNN", false);
  app_lw->train_ab_initio = glua_tbl_get_bool(L, "trainAbInitio", true);
  app_lw->nn_width = glua_tbl_get_integer(L, "NNWidth", 256);
  app_lw->nn_depth = glua_tbl_get_integer(L, "NNDepth", 5);
  const char* train_nn_file_char = glua_tbl_get_string(L, "trainNNFile", "");
  strcpy(app_lw->train_nn_file, train_nn_file_char);
  app_lw->num_trains = glua_tbl_get_integer(L, "numTrains", INT_MAX);
  app_lw->num_nn_writes = glua_tbl_get_integer(L, "numNNWrites", 1);

  app_lw->num_input_moms = 0;
  if (glua_tbl_has_key(L, "inputMoms")) {
    with_lua_tbl_tbl(L, "inputMoms") {
      app_lw->num_input_moms = glua_objlen(L);
      app_lw->input_moms = gkyl_malloc(sizeof(int[app_lw->num_input_moms]));

      for (int i = 0; i < app_lw->num_input_moms; i++) {
        // Indices are off by 1 between Lua and C.
        app_lw->input_moms[i] = glua_tbl_iget_integer(L, i + 1, 0) - 1;
      }
    }
  }

  app_lw->num_output_moms = 0;
  if (glua_tbl_has_key(L, "outputMoms")) {
    with_lua_tbl_tbl(L, "outputMoms") {
      app_lw->num_output_moms = glua_objlen(L);
      app_lw->output_moms = gkyl_malloc(sizeof(int[app_lw->num_output_moms]));

      for (int i = 0; i < app_lw->num_output_moms; i++) {
        // Indices are off by 1 between Lua and C.
        app_lw->output_moms[i] = glua_tbl_iget_integer(L, i + 1, 0) - 1;
      }
    }
  }

  app_lw->test_nn = glua_tbl_get_bool(L, "testNN", false);
  const char* test_nn_file_char = glua_tbl_get_string(L, "testNNFile", "");
  strcpy(app_lw->test_nn_file, test_nn_file_char);
  app_lw->num_tests = glua_tbl_get_integer(L, "numTests", 1);

  struct gkyl_pkpm pkpm = { }; // Input table for app.

  strcpy(pkpm.name, sim_name);
  
  int cdim = 0;
  with_lua_tbl_tbl(L, "cells") {
    pkpm.cdim = cdim = glua_objlen(L);

    for (int d = 0; d < cdim; d++) {
      pkpm.cells[d] = glua_tbl_iget_integer(L, d + 1, 0);
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
      pkpm.lower[d] = glua_tbl_iget_number(L, d + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d = 0; d < cdim; d++) {
      pkpm.upper[d] = glua_tbl_iget_number(L, d + 1, 0);
    }
  }

  pkpm.cfl_frac = glua_tbl_get_number(L, "cflFrac", 0.95);
  pkpm.poly_order = glua_tbl_get_integer(L, "polyOrder", 1);

  pkpm.use_explicit_source = glua_tbl_get_bool(L, "useExplicitSource", false);

  pkpm.basis_type = get_basis_type(
    glua_tbl_get_string(L, "basis", "serendipity")
  );

  pkpm.num_periodic_dir = 0;
  if (glua_tbl_has_key(L, "periodicDirs")) {
    with_lua_tbl_tbl(L, "periodicDirs") {
      pkpm.num_periodic_dir = glua_objlen(L);

      for (int d = 0; d < pkpm.num_periodic_dir; d++) {
        // Indices are off by 1 between Lua and C.
        pkpm.periodic_dirs[d] = glua_tbl_iget_integer(L, d + 1, 0) - 1;
      }
    }
  }

  struct pkpm_species_lw *species[GKYL_MAX_SPECIES];

  // Set all species input.
  pkpm.num_species = get_species_inp(L, cdim, species);

  // need to sort the species[] array by name of the species before
  // proceeding as there is no way to ensure that all cores loop over
  // Lua tables in the same order
  qsort(species, pkpm.num_species, sizeof(struct pkpm_species_lw *), species_compare_func);
  
  for (int s = 0; s < pkpm.num_species; s++) {
    pkpm.species[s] = species[s]->pkpm_species;
    pkpm.vdim = species[s]->vdim;

    app_lw->dist_init_func_ctx[s] = species[s]->dist_init_func_ref;
    app_lw->fluid_init_func_ctx[s] = species[s]->fluid_init_func_ref;
    if (species[s]->has_dist_init_func) {
      pkpm.species[s].init_dist = gkyl_lw_eval_cb;
      pkpm.species[s].ctx_dist = &app_lw->dist_init_func_ctx[s];
    }
    if (species[s]->has_fluid_init_func) {
      pkpm.species[s].init_fluid = gkyl_lw_eval_cb;
      pkpm.species[s].ctx_fluid = &app_lw->fluid_init_func_ctx[s];
    }

    app_lw->applied_acceleration_func_ctx[s] = species[s]->applied_acceleration_func_ref;
    if (species[s]->has_applied_acceleration_func) {
      pkpm.species[s].app_accel = gkyl_lw_eval_cb;
      pkpm.species[s].app_accel_ctx = &app_lw->applied_acceleration_func_ctx[s];
      pkpm.species[s].app_accel_evolve = species[s]->evolve_applied_acceleration;
    }

    if (species[s]->has_diffusion) {
      pkpm.species[s].diffusion.D = species[s]->D; 
      pkpm.species[s].diffusion.order = species[s]->order; 
    }

    app_lw->collision_id[s] = species[s]->collision_id;
    pkpm.species[s].collisions.collision_id = app_lw->collision_id[s];

    app_lw->self_nu_func_ctx[s] = species[s]->self_nu_func_ref;
    if (species[s]->has_self_nu_func) {
      pkpm.species[s].collisions.self_nu = gkyl_lw_eval_cb;
      pkpm.species[s].collisions.ctx = &app_lw->self_nu_func_ctx[s];
    }

    app_lw->num_cross_collisions[s] = species[s]->num_cross_collisions;
    for (int i = 0; i < app_lw->num_cross_collisions[s]; i++) {
      strcpy(app_lw->collide_with[s][i], species[s]->collide_with[i]);
    }

    pkpm.species[s].collisions.num_cross_collisions = app_lw->num_cross_collisions[s];
    for (int i = 0; i < app_lw->num_cross_collisions[s]; i++) {
      strcpy(pkpm.species[s].collisions.collide_with[i], app_lw->collide_with[s][i]);
    }
  }

  // Set field input.
  with_lua_tbl_key(L, "field") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct pkpm_field_lw *pkpm_f = lua_touserdata(L, -1);

      if (pkpm_f->magic == PKPM_FIELD_DEFAULT) {
        pkpm_f->init_ref.ndim = cdim;

        pkpm.field = pkpm_f->pkpm_field;

        app_lw->field_func_ctx = pkpm_f->init_ref;
        pkpm.field.init = gkyl_lw_eval_cb;
        pkpm.field.ctx = &app_lw->field_func_ctx;

        if (pkpm_f->has_external_field_func) {
          pkpm_f->external_field_func_ref.ndim = cdim;

          app_lw->external_field_func_ctx = pkpm_f->external_field_func_ref;
          pkpm.field.ext_em = gkyl_lw_eval_cb;
          pkpm.field.ext_em_ctx = &app_lw->external_field_func_ctx;

          pkpm.field.ext_em_evolve = pkpm_f->evolve_external_field;
        }

        if (pkpm_f->has_applied_current_func) {
          pkpm_f->applied_current_func_ref.ndim = cdim;

          app_lw->applied_current_func_ctx = pkpm_f->applied_current_func_ref;
          pkpm.field.app_current = gkyl_lw_eval_cb;
          pkpm.field.app_current_ctx = &app_lw->applied_current_func_ctx;

          pkpm.field.app_current_evolve = pkpm_f->evolve_applied_current;
        }
      }
    }
  }

  // Create parallelism.
  struct gkyl_comm *comm = 0;

  for (int d = 0; d < cdim; d++) {
    pkpm.parallelism.cuts[d] = cuts[d]; 
  }

  struct gkyl_tool_args *args = gkyl_tool_args_new(L);
  struct script_cli script_cli = pkpm_parse_script_cli(args);

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

  pkpm.parallelism.comm = comm;
  pkpm.parallelism.use_gpu = script_cli.use_gpu;

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
  
  app_lw->app = gkyl_pkpm_app_new(&pkpm);

  gkyl_comm_release(comm);

  // Create Lua userdata.
  struct pkpm_app_lw **l_app_lw = lua_newuserdata(L, sizeof(struct pkpm_app_lw*));
  *l_app_lw = app_lw; // Point it to the Lua app pointer.

  // Set metatable.
  luaL_getmetatable(L, PKPM_APP_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Apply initial conditions. (time) -> bool.
static int
pkpm_app_apply_ic(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->t_start);
  gkyl_pkpm_app_apply_ic(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to field. (time) -> bool.
static int
pkpm_app_apply_ic_field(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->t_start);
  gkyl_pkpm_app_apply_ic_field(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to species. (sidx, time) -> bool.
static int
pkpm_app_apply_ic_species(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double t0 = luaL_optnumber(L, 3, app_lw->t_start);
  gkyl_pkpm_app_apply_ic_species(app_lw->app, sidx, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated moments. (tm) -> bool.
static int
pkpm_app_calc_integrated_mom(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_pkpm_app_calc_integrated_mom(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated L2 norm of distribution function. (tm) -> bool.
static int
pkpm_app_calc_integrated_L2_f(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_pkpm_app_calc_integrated_L2_f(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated field energy (L2 norm of each field
// component). (tm) -> bool.
static int
pkpm_app_calc_field_energy(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_pkpm_app_calc_field_energy(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Write solution (field and species) to file (time, frame) -> bool.
static int
pkpm_app_write(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_pkpm_app_write(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write field to file (time, frame) -> bool.
static int
pkpm_app_write_field(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_pkpm_app_write_field(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write species solution to file (sidx, time, frame) -> bool.
static int
pkpm_app_write_species(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double tm = luaL_checknumber(L, 3);
  int frame = luaL_checkinteger(L, 4);
  gkyl_pkpm_app_write_species(app_lw->app, sidx, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated moments to file () -> bool.
static int
pkpm_app_write_integrated_mom(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  gkyl_pkpm_app_write_integrated_mom(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated L2 norm of f to file () -> bool.
static int
pkpm_app_write_integrated_L2_f(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  gkyl_pkpm_app_write_integrated_L2_f(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated field energy to file () -> bool.
static int
pkpm_app_write_field_energy(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  gkyl_pkpm_app_write_field_energy(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write simulation statistics to JSON. () -> bool.
static int
pkpm_app_stat_write(lua_State *L)
{
  bool status = true;

  struct pkpm_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  gkyl_pkpm_app_stat_write(app_lw->app);

  lua_pushboolean(L, status);
  return 1;  
}

// Write data from simulation to file.
static void
write_data(struct gkyl_tm_trigger* iot, gkyl_pkpm_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_pkpm_app_write(app, t_curr, frame);
    gkyl_pkpm_app_write_field_energy(app);
    gkyl_pkpm_app_write_integrated_mom(app);
    gkyl_pkpm_app_write_integrated_L2_f(app);
  }
}

// Calculate and append field energy to dynvector.
static void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_pkpm_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_pkpm_app_calc_field_energy(app, t_curr);
  }
}

// Calculate and append integrated moments to dynvector.
static void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_pkpm_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_pkpm_app_calc_integrated_mom(app, t_curr);
  }
}

// Calculate and append integrated L2 norm of distribution function to dynvector.
static void
calc_integrated_L2_f(struct gkyl_tm_trigger* l2t, gkyl_pkpm_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(l2t, t_curr) || force_calc) {
    gkyl_pkpm_app_calc_integrated_L2_f(app, t_curr);
  }
}

static void
train_mom(struct gkyl_tm_trigger* nn, gkyl_pkpm_app* app, double t_curr, bool force_train, struct gkyl_kann_net** ann, int num_input_moms, int* input_moms, int num_output_moms, int* output_moms,
  struct gkyl_kn_vec* input_data, struct gkyl_kn_vec* output_data)
{
  if (gkyl_tm_trigger_check_and_bump(nn, t_curr) || force_train) {
    int frame = nn->curr - 1;
    if (force_train) {
      frame = nn->curr;
    }

    gkyl_pkpm_app_train(app, t_curr, frame, ann, num_input_moms, input_moms, num_output_moms, output_moms, input_data, output_data);
  }
}

static void
write_nn(struct gkyl_tm_trigger* nnw, gkyl_pkpm_app* app, double t_curr, bool force_write, struct gkyl_kann_net** ann)
{
  if (gkyl_tm_trigger_check_and_bump(nnw, t_curr) || force_write) {
    int frame = nnw->curr - 1;
    if (force_write) {
      frame = nnw->curr;
    }

    gkyl_pkpm_app_write_nn(app, t_curr, frame, ann);
  }
}

static void
test_mom(struct gkyl_tm_trigger* nnt, gkyl_pkpm_app* app, double t_curr, bool force_test, struct gkyl_kann_net** ann, int num_input_moms, int* input_moms, int num_output_moms, int* output_moms,
  struct gkyl_kn_vec* input_data_real, struct gkyl_kn_vec* output_data_real, struct gkyl_kn_vec* output_data_predicted)
{
  if (gkyl_tm_trigger_check_and_bump(nnt, t_curr) || force_test) {
    int frame = nnt->curr - 1;
    if (force_test) {
      frame = nnt->curr;
    }

    gkyl_pkpm_app_test(app, t_curr, frame, ann, num_input_moms, input_moms, num_output_moms, output_moms, input_data_real, output_data_real, output_data_predicted);
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
write_step_message(const struct gkyl_pkpm_app *app, struct step_message_trigs *trigs, int step, double t_curr, double dt_next)
{
  if (gkyl_tm_trigger_check_and_bump(&trigs->log_trig, t_curr)) {
    if (trigs->log_count > 0) {
      gkyl_pkpm_app_cout(app, stdout, " Step %6d at time %#11.8g.  Time-step  %.6e.  Completed %g%s\n", step, t_curr, dt_next, trigs->tenth * 10.0, "%");
    }
    else {
      trigs->log_count += 1;
    }
    
    trigs->tenth += 1;
  }
  if (gkyl_tm_trigger_check_and_bump(&trigs->log_trig_1p, t_curr)) {
    gkyl_pkpm_app_cout(app, stdout, "%d", trigs->p1c);
    trigs->p1c = (trigs->p1c+1) % 10;
  }
}

static void
show_help(const struct gkyl_pkpm_app *app)
{
  gkyl_pkpm_app_cout(app, stdout, "PKPM script takes the following arguments:\n");
  gkyl_pkpm_app_cout(app, stdout, " -h   Print this help message and exit\n");
  gkyl_pkpm_app_cout(app, stdout, " -sN  Only run N steps of simulation\n");
  gkyl_pkpm_app_cout(app, stdout, " -S   Do not initialize MPI\n");
  gkyl_pkpm_app_cout(app, stdout, " -G   Do not initialize CUDA\n");
  gkyl_pkpm_app_cout(app, stdout, " -m   Run memory tracer\n");
  gkyl_pkpm_app_cout(app, stdout, " -V   Show verbose output\n");
  gkyl_pkpm_app_cout(app, stdout, " -rN  Restart simulation from frame N\n");

  gkyl_pkpm_app_cout(app, stdout, "\n");
}

// Run simulation. (num_steps) -> bool. num_steps is optional.
static int
pkpm_app_run(lua_State *L)
{
  bool ret_status = true;

  // Create app object.
  struct pkpm_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;
  struct gkyl_pkpm_app *app = app_lw->app;

  // Parse command lines arguments passed to input file.
  struct gkyl_tool_args *args = gkyl_tool_args_new(L);

  struct script_cli script_cli = pkpm_parse_script_cli(args);
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

  gkyl_pkpm_app_cout(app, stdout, "Initializing PKPM Simulation ...\n");

  if (script_cli.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  // Initialize simulation.
  bool is_restart = script_cli.is_restart;
  int restart_frame = script_cli.restart_frame;  

  int frame_curr = 0;
  if (is_restart) {
    struct gkyl_app_restart_status status = gkyl_pkpm_app_read_from_frame(app, restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_pkpm_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_pkpm_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_pkpm_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_pkpm_app_apply_ic(app, t_curr);
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

  // Create trigger for neural network training.
  int num_trains = app_lw->num_trains;
  struct gkyl_tm_trigger nn_trig = { .dt = t_end / num_trains, .tcurr = t_curr, .curr = frame_curr };

  kad_node_t **t = gkyl_malloc(sizeof(kad_node_t*) * app->num_species);
  struct gkyl_kann_net **ann = gkyl_malloc(sizeof(struct gkyl_kann_net*) * app->num_species);
  if (app_lw->train_nn) {
    if (app_lw->train_ab_initio) {
      for (int i = 0; i < app->num_species; i++ ) {
        if (app->poly_order == 1) {
          if (app->cdim == 1) {
            t[i] = kann_layer_input(app_lw->num_input_moms * 2);
          }
          else if (app->cdim == 2) {
            t[i] = kann_layer_input(app_lw->num_input_moms * 4);
          }
        }
        else if (app->poly_order == 2) {
          t[i] = kann_layer_input(app_lw->num_input_moms * 3);
        }

        for (int j = 0; j < app_lw->nn_depth; j++) {
          t[i] = kann_layer_dense(t[i], app_lw->nn_width);
          t[i] = kad_tanh(t[i]);
        }
        
        if (app->poly_order == 1) {
          if (app->cdim == 1) {
            t[i] = kann_layer_cost(t[i], app_lw->num_output_moms * 2, KANN_C_MSE);
          }
          else if (app->cdim == 2) {
            t[i] = kann_layer_cost(t[i], app_lw->num_output_moms * 4, KANN_C_MSE);
          }
        }
        else if (app->poly_order == 2) {
          t[i] = kann_layer_cost(t[i], app_lw->num_output_moms * 3, KANN_C_MSE);
        }
        ann[i] = gkyl_kann_net_new(t[i], app->use_gpu);
      }
    }
    else {
      for (int i = 0; i < app->num_species; i++) {
        const char *fmt = "%s-%s.dat";
        int sz = gkyl_calc_strlen(fmt, app_lw->train_nn_file, app->species[i].info.name);
        char fileNm[sz + 1];
        snprintf(fileNm, sizeof fileNm, fmt, app_lw->train_nn_file, app->species[i].info.name);

        FILE *file = fopen(fileNm, "r");
        if (file != NULL) {
          ann[i] = gkyl_kann_net_load(fileNm, app->use_gpu);
          fclose(file);
        }
        else {
          ann[i] = 0;
          app_lw->train_nn = false;
          fprintf(stderr, "Neural network for species %s not found! Disabling NN training.\n", app->species[i].info.name);
        }
      }
    }
  }

  int cell_count = 0;
  if (app->cdim == 1) {
    cell_count = app->grid.cells[0];
  }
  else if (app->cdim == 2) {
    cell_count = app->grid.cells[0] * app->grid.cells[1];
  }

  // Number of DG coefficients per moment carried into/out of the network.
  int nn_num_coeff = (app->poly_order == 1) ? ((app->cdim == 1) ? 2 : 4) : 3;
  int nn_input_dim = app_lw->num_input_moms * nn_num_coeff;
  int nn_output_dim = app_lw->num_output_moms * nn_num_coeff;

  struct gkyl_kn_vec *input_data = gkyl_kn_vec_new(cell_count, nn_input_dim);
  struct gkyl_kn_vec *output_data = gkyl_kn_vec_new(cell_count, nn_output_dim);

  if (app_lw->train_nn) {
    train_mom(&nn_trig, app, t_curr, false, ann, app_lw->num_input_moms, app_lw->input_moms, app_lw->num_output_moms, app_lw->output_moms, input_data, output_data);
  }

  // Create trigger for neural network writing.
  int num_nn_writes = app_lw->num_nn_writes;
  struct gkyl_tm_trigger nnw_trig = { .dt = t_end / num_nn_writes, .tcurr = t_curr, .curr = frame_curr };

  if (app_lw->train_nn) {
    write_nn(&nnw_trig, app, t_curr, false, ann);
  }

  // Create trigger for neural network testing.
  int num_tests = app_lw->num_tests;
  struct gkyl_tm_trigger nnt_trig = { .dt = t_end / num_tests, .tcurr = t_curr, .curr = frame_curr };

  struct gkyl_kann_net **ann_test = gkyl_malloc(sizeof(struct gkyl_kann_net*) * app->num_species);
  if (app_lw->test_nn) {
    for (int i = 0; i < app->num_species; i++) {
      const char *fmt = "%s-%s.dat";
      int sz = gkyl_calc_strlen(fmt, app_lw->test_nn_file, app->species[i].info.name);
      char fileNm[sz + 1];
      snprintf(fileNm, sizeof fileNm, fmt, app_lw->test_nn_file, app->species[i].info.name);

      FILE *file = fopen(fileNm, "r");
      if (file != NULL) {
        ann_test[i] = gkyl_kann_net_load(fileNm, app->use_gpu);
        fclose(file);
      }
      else {
        ann_test[i] = 0;
        app_lw->test_nn = false;
        fprintf(stderr, "Neural network for species %s not found! Disabling NN testing.\n", app->species[i].info.name);
      }
    }
  }

  struct gkyl_kn_vec *input_data_real = gkyl_kn_vec_new(cell_count, nn_input_dim);
  struct gkyl_kn_vec *output_data_real = gkyl_kn_vec_new(cell_count, nn_output_dim);
  struct gkyl_kn_vec *output_data_predicted = gkyl_kn_vec_new(cell_count, nn_output_dim);

  if (app_lw->test_nn) {
    test_mom(&nnt_trig, app, t_curr, false, ann_test, app_lw->num_input_moms, app_lw->input_moms, app_lw->num_output_moms, app_lw->output_moms,
      input_data_real, output_data_real, output_data_predicted);
  }

  gkyl_pkpm_app_cout(app, stdout, "Initialization completed in %g sec\n\n", gkyl_time_diff_now_sec(tm_ic0));
  
  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = app_lw->dt_failure_tol;
  int num_failures = 0, num_failures_max = app_lw->num_failures_max;

  bool use_verbose = script_cli.use_verbose;

  long step = 1;
  while ((t_curr < t_end) && (step <= num_steps)) {
    if (use_verbose) {
      gkyl_pkpm_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    }
    struct gkyl_update_status status = gkyl_pkpm_update(app, dt);
    if (use_verbose) {
      gkyl_pkpm_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    }

    if (!status.success) {
      gkyl_pkpm_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);
    if (app_lw->train_nn) {
      train_mom(&nn_trig, app, t_curr, false, ann, app_lw->num_input_moms, app_lw->input_moms, app_lw->num_output_moms, app_lw->output_moms, input_data, output_data);
      write_nn(&nnw_trig, app, t_curr, false, ann);
    }
    if (app_lw->test_nn) {
      test_mom(&nnt_trig, app, t_curr, false, ann_test, app_lw->num_input_moms, app_lw->input_moms, app_lw->num_output_moms, app_lw->output_moms,
        input_data_real, output_data_real, output_data_predicted);
    }

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_pkpm_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_pkpm_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_pkpm_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_pkpm_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_pkpm_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        calc_integrated_L2_f(&l2f_trig, app, t_curr, true);
        write_data(&io_trig, app, t_curr, true);
        if (app_lw->train_nn) {
          train_mom(&nn_trig, app, t_curr, true, ann, app_lw->num_input_moms, app_lw->input_moms, app_lw->num_output_moms, app_lw->output_moms, input_data, output_data);
          write_nn(&nnw_trig, app, t_curr, true, ann);
        }
        if (app_lw->test_nn) {
          test_mom(&nnt_trig, app, t_curr, true, ann_test, app_lw->num_input_moms, app_lw->input_moms, app_lw->num_output_moms, app_lw->output_moms,
            input_data_real, output_data_real, output_data_predicted);
        }

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
  if (app_lw->train_nn) {
    train_mom(&nn_trig, app, t_curr, false, ann, app_lw->num_input_moms, app_lw->input_moms, app_lw->num_output_moms, app_lw->output_moms, input_data, output_data);
    write_nn(&nnw_trig, app, t_curr, false, ann);

    for (int i = 0; i < app->num_species; i++) {
      gkyl_kann_net_release(ann[i]);
    }
  }
  if (app_lw->test_nn) {
    test_mom(&nnt_trig, app, t_curr, false, ann_test, app_lw->num_input_moms, app_lw->input_moms, app_lw->num_output_moms, app_lw->output_moms,
      input_data_real, output_data_real, output_data_predicted);

    for (int i = 0; i < app->num_species; i++) {
      gkyl_kann_net_release(ann_test[i]);
    }
  }
  gkyl_pkpm_app_stat_write(app);

  struct gkyl_pkpm_stat stat = gkyl_pkpm_app_stat(app);

  gkyl_pkpm_app_cout(app, stdout, "\n");
  gkyl_pkpm_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_pkpm_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_pkpm_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_pkpm_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_pkpm_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_pkpm_app_cout(app, stdout, "Fluid species RHS calc took %g secs\n", stat.fluid_species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species PKPM vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_pkpm_app_cout(app, stdout, "EM variables (bvar) calc took %g secs\n", stat.field_em_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm);
  gkyl_pkpm_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_pkpm_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_pkpm_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

freeresources:

  gkyl_free(app_lw->input_moms);
  gkyl_free(app_lw->output_moms);
  gkyl_free(t);
  gkyl_free(ann);
  gkyl_free(ann_test);

  gkyl_kn_vec_release(input_data);
  gkyl_kn_vec_release(output_data);
  gkyl_kn_vec_release(input_data_real);
  gkyl_kn_vec_release(output_data_real);
  gkyl_kn_vec_release(output_data_predicted);

  lua_pushboolean(L, ret_status);
  return 1;
}

// Clean up memory allocated for simulation.
static int
pkpm_app_gc(lua_State *L)
{
  struct pkpm_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, PKPM_APP_METATABLE_NM);
  struct pkpm_app_lw *app_lw = *l_app_lw;

  gkyl_pkpm_app_release(app_lw->app);
  gkyl_free(*l_app_lw);
  
  return 0;
}

// App constructor.
static struct luaL_Reg pkpm_app_ctor[] = {
  { "new", pkpm_app_new },
  { 0, 0 }
};

// App methods.
static struct luaL_Reg pkpm_app_funcs[] = {
  { "apply_ic", pkpm_app_apply_ic },
  { "apply_ic_field", pkpm_app_apply_ic_field },
  { "apply_ic_species", pkpm_app_apply_ic_species },
  { "calc_integrated_mom", pkpm_app_calc_integrated_mom },
  { "calc_integrated_L2_f", pkpm_app_calc_integrated_L2_f },
  { "calc_field_energy", pkpm_app_calc_field_energy },
  { "write", pkpm_app_write },
  { "write_field", pkpm_app_write_field },
  { "write_species", pkpm_app_write_species },
  { "write_integrated_mom", pkpm_app_write_integrated_mom },
  { "write_integrated_L2_f", pkpm_app_write_integrated_L2_f },
  { "write_field_energy", pkpm_app_write_field_energy },
  { "stat_write", pkpm_app_stat_write },
  { "run", pkpm_app_run },
  { 0, 0 }
};

static void
app_openlibs(lua_State *L)
{
  // Register top-level App.
  do {
    luaL_newmetatable(L, PKPM_APP_METATABLE_NM);

    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, pkpm_app_gc);
    lua_settable(L, -3);

    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, pkpm_app_funcs);
    
    luaL_register(L, "G0.PKPM.App", pkpm_app_ctor);
    
  }
  while (0);

  // Register Species input struct.
  do {
    luaL_newmetatable(L, PKPM_SPECIES_METATABLE_NM);
    luaL_register(L, "G0.PKPM.Species", pkpm_species_ctor);
  }
  while (0);

  // Register Field input struct.
  do {
    luaL_newmetatable(L, PKPM_FIELD_METATABLE_NM);
    luaL_register(L, "G0.PKPM.Field", pkpm_field_ctor);
  }
  while (0);
}

void
gkyl_pkpm_lw_openlibs(lua_State *L)
{
  app_openlibs(L);
}

#endif