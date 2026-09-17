// Gkyl ------------------------------------------------------------------------
//
// Top-level entry point into Gkeyll
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifdef GKYL_HAVE_LUA

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <lfs.h>
#include <whereami.h>

#include <gkyl_alloc.h>
#include <gkyl_lua_utils.h>
#include <gkyl_util.h>


#include <gkyl_zero_lw.h>
#include <gkyl_lw_priv.h>

#ifdef GKYL_HAVE_MOMENTS
#include <gkyl_moment_lw.h>
#endif

#ifdef GKYL_HAVE_VLASOV
#include <gkyl_vlasov_lw.h>
#endif

#ifdef GKYL_HAVE_PKPM
#include <gkyl_pkpm_lw.h>
#endif

#ifdef GKYL_HAVE_GYROKINETIC
#include <gkyl_gyrokinetic_lw.h>
#endif

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#endif

#define STRINGIFY_(x)  #x
#define STRINGIFY(x)   STRINGIFY_(x)



// Tool description struct
struct tool_description {
  const char *tool_name;
  const char *tool_lua;
  const char *tool_help;
};

// List of available Tools
static struct tool_description tool_list[] = {
  {"man", "man.lua", "Gkeyll online manual"},
  {"woman", "man.lua", "Gkeyll online manual (Woe without man)"},
  {"queryrdb", "queryrdb.lua", "Query/modify regression test DB"},
  {"exacteulerrp", "exacteulerrp.lua", "Exact Euler Riemann problem solver"},
  {"runregression", "runregression.lua", "Run regression/unit tests"},
  {"multimomlinear", "multimomlinear.lua",
   "Linear dispersion solver for multi-moment, multifluid equations"},  
  {"eqdskreader", "eqdskreader.lua", "Read eqdsk file, writing data to files"},
  {0, 0}
};

static int max2(int a, int b) { return a>b ? a : b; }

// Show list of available Tools
static void
show_tool_list(void)
{
  fprintf(stdout, "Following tools are available. Query tool help for more information.\n\n");

  int mlen = 0;
  for (int i=0; tool_list[i].tool_name != 0; ++i) {
    int len = strlen(tool_list[i].tool_name);
    mlen = len > mlen ? len : mlen;
  }
  
  for (int i=0; tool_list[i].tool_name != 0; ++i)
    fprintf(stdout, "%*s %s\n", mlen+2, tool_list[i].tool_name, tool_list[i].tool_help);
  fprintf(stdout, "\n");
}

// Returns tool Lua script name given tool name. Returns 0 if Tool
// does no exist
static const char *
get_tool_from_name(const char *nm)
{
  for (int i=0; tool_list[i].tool_name != 0; ++i)
    if (strcmp(tool_list[i].tool_name, nm) == 0)
      return tool_list[i].tool_lua;
  return 0;
}

static char *
find_exec_path(void)
{
  int len = wai_getExecutablePath(NULL, 0, NULL);
  char *path = gkyl_malloc(len+1);
  int dirname_len; wai_getExecutablePath(path, len, &dirname_len);
  path[dirname_len] = '\0'; // only directory path is returned
  return path;
}

static int
calc_output_prefix_len(const char *fn)
{
  const char *suff = strrchr(fn, '.');
  return strlen(fn) - (suff ? strlen(suff) : 0);
}

static const char*
get_fname(const char *fn)
{
  return strrchr(fn, '/');
}

static const char*
get_fname_with_dir(const char *fn)
{
  return fn;
}

static void
show_banner(FILE *fp)
{
  if (fp) {
    fprintf(fp, "Gkeyll built with Git changeset %s\n", STRINGIFY(GKYL_GIT_CHANGESET));
    fprintf(fp, "Gkeyll built on %s\n", STRINGIFY(GKYL_BUILD_DATE));
#ifdef GKYL_HAVE_CUDA
    fprintf(fp, "Built with CUDA\n");
#else
    fprintf(fp, "Built without CUDA\n");
#endif
#ifdef GKYL_HAVE_NCCL
    fprintf(fp, "Built with NCCL\n");
#else
    fprintf(fp, "Built without NCCL\n");
#endif
#ifdef GKYL_HAVE_MPI
    fprintf(fp, "Built with MPI\n");
#else
    fprintf(fp, "Built without MPI\n");
#endif

    int num_apps = 0;
#ifdef GKYL_HAVE_MOMENTS
    fprintf(fp, "Moments App enabled\n");
    num_apps += 1;
#endif
#ifdef GKYL_HAVE_VLASOV
    fprintf(fp, "Vlasov App enabled\n");
    num_apps += 1;
#endif
#ifdef GKYL_HAVE_GYROKINETIC
    fprintf(fp, "Gyrokinetic App enabled\n");
    num_apps += 1;
#endif
#ifdef GKYL_HAVE_PKPM
    fprintf(fp, "PKPM App enabled\n");
    num_apps += 1;
#endif
    if (0 == num_apps)
      fprintf(fp, "Built without any Apps! Core library only.\n");
    fprintf(fp, "\n");
  }
}

// show usage
static void
show_usage()
{
  fprintf(stdout, "This is the Gkeyll code. See gkeyll.rtfd.io for details.\n");
  fprintf(stdout, "Type 'gkeyll man' for help.\n\n");

  show_banner(stdout);

  fprintf(stdout, "gkeyll [OPTIONS] input-file/tool-name [APP-OPTIONS]\n");
  fprintf(stdout, "Available options are\n");
  fprintf(stdout, "  -e chunk   Execute string 'chunk'\n");
  fprintf(stdout, "  -t         Show list of registered tools\n");
  fprintf(stdout, "  -v         Show version information\n");
  fprintf(stdout, "  -d         Write output to same directory as input file\n");  
  fprintf(stdout, "  -g         Run on NVIDIA GPU (if available and built with CUDA)\n\n");
  fprintf(stdout, "  -m         Run memory tracer\n");
  fprintf(stdout, "  -S         Do not initialize MPI\n");
  fprintf(stdout, "  -V         Use verbose output\n");

  fprintf(stdout, "Most app input files take the following commands:\n");
  fprintf(stdout, "  run        Run simulation. Default if nothing specified\n");
  fprintf(stdout, "  restart    Restart simulation \n");
  fprintf(stdout, "To get help for commands type command name followed by -h\n\n");

  fprintf(stdout, "Individual tools may take other options and commands. See their specific help.\n");
}

static void
show_version()
{
  fprintf(stdout, "This is the Gkeyll code. See gkeyll.rtfd.io for details.\n");
  fprintf(stdout, "Type 'gkyl -h' for help.\n\n");
#ifdef GKYL_GIT_CHANGESET  
  fprintf(stdout, "Built with git changeset %s\n", STRINGIFY(GKYL_GIT_CHANGESET));
#endif
#ifdef GKYL_BUILD_DATE
  fprintf(stdout, "Built on %s\n", STRINGIFY(GKYL_BUILD_DATE));
#endif
#ifdef GKYL_HAVE_NCCL
  fprintf(stdout, "Built with CUDA\n");
#else
  fprintf(stdout, "Built without CUDA\n");
#endif
}

// Input arguments to gkyl executable. You must call release_app_args
// when done with this struct
struct app_args {
  bool use_gpu; // should this be run on GPU?
  bool step_mode; // run for fixed number of steps? (for valgrind/cuda-memcheck)
  bool trace_mem; // should we trace memory allocation/deallocations?
  int num_steps; // number of steps
  bool use_mpi; // should we use MPI?
  bool use_verbose; // Should we use verbose output?
  bool is_restart; // Is this a restarted simulation?
  int restart_frame; // Which frame to restart simulation from.

  bool use_dir_path; /// should we use full directory path
  
  char *echunk; // chunk of lua code to execute
  
  int num_opt_args; // number of optional arguments
  char **opt_args; // optional arguments

  char *exec_path; // location of executable
};

static void
release_opt_args(struct app_args *args)
{
  for (int i=0; i<args->num_opt_args; ++i)
    gkyl_free(args->opt_args[i]);
  gkyl_free(args->opt_args);
  
  if (args->echunk)
    gkyl_free(args->echunk);

  gkyl_free(args->exec_path);
  gkyl_free(args);
}

static struct app_args*
parse_app_args(int argc, char **argv)
{
  struct app_args *args = gkyl_malloc(sizeof(*args));

  args->use_gpu = false;
  args->step_mode = false;

#ifdef GKYL_HAVE_MPI  
  args->use_mpi = true;
#else
  args->use_mpi = false;
#endif
  
  args->trace_mem = false;
  args->num_opt_args = 0;
  args->echunk = 0;
  args->use_verbose = false;
  args->num_steps = -1;

  args->is_restart = false;
  args->restart_frame = 0;

  args->use_dir_path = false;

  int c;
  while ((c = getopt(argc, argv, "+hvtmdSe:gV")) != -1) {
    switch (c)
    {
      case 'h':
        show_usage();
        exit(-1);
        break;

      case 'v':
        show_version();
        exit(-1);
        break;        

      case 't':
        show_tool_list();
        exit(1);
        break;

      case 'e':
        args->echunk = gkyl_malloc(strlen(optarg)+1);
        strcpy(args->echunk, optarg);
        break;
        
      case 'g':
        args->use_gpu = true;
        break;

      case 'd':
        args->use_dir_path = true;
        break;        

      case 'S':
        args->use_mpi = false;
        break;
        
      case 'm':
        args->trace_mem = true;
        break;
      
      case 'V':
        args->use_verbose = true;
        break;
      
      case '?':
        break;
    }
  }

  args->num_opt_args = 0;
  // collect remaining options into a list
  for (int oind=optind; oind < argc; ++oind) args->num_opt_args += 1;
  args->opt_args = gkyl_malloc(sizeof(char*)*args->num_opt_args);

  for (int i=0, oind=optind; oind < argc; ++oind, ++i) {
    args->opt_args[i] = gkyl_malloc(strlen(argv[oind])+1);
    strcpy(args->opt_args[i], argv[oind]);
  }

  args->exec_path = find_exec_path();
  
  return args;
}



int
main(int argc, char **argv)
{
  struct app_args *app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args->use_mpi)
    MPI_Init(&argc, &argv);
#endif

  if (app_args->trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  
  lua_State *L = luaL_newstate();
  lua_gc(L, LUA_GCSTOP, 0);
  luaL_openlibs(L);
  luaopen_lfs(L); // Lua file-system library

  // G0 librararies
  gkyl_zero_lw_openlibs(L);
#ifdef GKYL_HAVE_MOMENTS 
  gkyl_moment_lw_openlibs(L);
#endif
#ifdef GKYL_HAVE_VLASOV
  gkyl_vlasov_lw_openlibs(L);
#endif
#ifdef GKYL_HAVE_GYROKINETIC  
  gkyl_gyrokinetic_lw_openlibs(L);
#endif
#ifdef GKYL_HAVE_PKPM  
  gkyl_pkpm_lw_openlibs(L);
#endif  
  lua_gc(L, LUA_GCRESTART, -1);

#ifdef GKYL_HAVE_MPI
  struct {
    MPI_Comm comm;
  } lw_mpi_comm_world = { .comm = MPI_COMM_WORLD };
#endif
  
  if (app_args->use_mpi) {
    lua_pushboolean(L, true);
    lua_setglobal(L, "GKYL_HAVE_MPI");

#ifdef GKYL_HAVE_MPI
    lua_pushlightuserdata(L, &lw_mpi_comm_world);
#else
    lua_pushlightuserdata(L, false);
#endif
    lua_setglobal(L, "GKYL_MPI_COMM");
  }
  else {
    lua_pushboolean(L, false);
    lua_setglobal(L, "GKYL_HAVE_MPI");
    
    lua_pushboolean(L, false);
    lua_setglobal(L, "GKYL_MPI_COMM");
  }

#ifdef GKYL_HAVE_NCCL
  lua_pushboolean(L, true);
  lua_setglobal(L, "GKYL_HAVE_CUDA");

  if (app_args->use_gpu) {
    lua_pushboolean(L, true);
    lua_setglobal(L, "GKYL_USE_GPU");
  }
#else
  lua_pushboolean(L, false);
  lua_setglobal(L, "GKYL_HAVE_CUDA");

  lua_pushboolean(L, false);
  lua_setglobal(L, "GKYL_USE_GPU");  
#endif

  if (app_args->use_verbose) {
    lua_pushboolean(L, true);
  }
  else {
    lua_pushboolean(L, false);
  }
  lua_setglobal(L, "GKYL_USE_VERBOSE");

  lua_pushinteger(L, app_args->num_steps);
  lua_setglobal(L, "GKYL_NUM_STEPS");

  if (app_args->is_restart) {
    lua_pushboolean(L, true);
  }
  else {
    lua_pushboolean(L, false);
  }
  lua_setglobal(L, "GKYL_IS_RESTART");

  lua_pushinteger(L, app_args->restart_frame);
  lua_setglobal(L, "GKYL_RESTART_FRAME");

  lua_pushnumber(L, DBL_MIN);
  lua_setglobal(L, "GKYL_MIN_DOUBLE");
  lua_pushnumber(L, DBL_MAX);
  lua_setglobal(L, "GKYL_MAX_DOUBLE");

  lua_pushnumber(L, FLT_MIN);
  lua_setglobal(L, "GKYL_MIN_FLOAT");
  lua_pushnumber(L, FLT_MAX);
  lua_setglobal(L, "GKYL_MAX_FLOAT");

  lua_pushinteger(L, INT_MIN);
  lua_setglobal(L, "GKYL_MIN_INT");
  lua_pushinteger(L, INT_MAX);
  lua_setglobal(L, "GKYL_MAX_INT");

  lua_pushinteger(L, LONG_MIN);
  lua_setglobal(L, "GKYL_MIN_LONG");
  lua_pushinteger(L, LONG_MAX);
  lua_setglobal(L, "GKYL_MAX_LONG");

  lua_pushinteger(L, LONG_MIN);
  lua_setglobal(L, "GKYL_MIN_LONG_LONG");
  lua_pushinteger(L, LONG_MAX);
  lua_setglobal(L, "GKYL_MAX_LONG_LONG");

  lua_pushnumber(L, DBL_EPSILON);
  lua_setglobal(L, "GKYL_EPSILON");
  
  lua_pushinteger(L, INT16_MAX);
  lua_setglobal(L, "GKYL_MAX_INT16");

  lua_pushstring(L, app_args->exec_path);
  lua_setglobal(L, "GKYL_EXEC_PATH");
  
  do {
    const char *fmt = "%s/gkyl";
    size_t len = gkyl_calc_strlen(fmt, app_args->exec_path);

    char *str = gkyl_malloc(len+1);
    snprintf(str, len+1, fmt, app_args->exec_path);

    lua_pushstring(L, str);
    lua_setglobal(L, "GKYL_EXEC");

    gkyl_free(str);
  } while (0);  

#ifdef GKYL_GIT_CHANGESET
  lua_pushstring(L, STRINGIFY(GKYL_GIT_CHANGESET));
  lua_setglobal(L, "GKYL_GIT_CHANGESET");
#endif  

#ifdef GKYL_BUILD_DATE
  lua_pushstring(L, STRINGIFY(GKYL_BUILD_DATE));
  lua_setglobal(L, "GKYL_BUILD_DATE");
#endif  

  // push extra arguments into a Lua table to Tools and App can get
  // them
  lua_newtable(L);
  if (app_args->num_opt_args > 0)  {
    lua_pushinteger(L, 1);
    lua_pushstring(L, app_args->opt_args[0]);
    lua_rawset(L, -3);    
  }
  
  for (int i=1; i<app_args->num_opt_args; ++i) {
    lua_pushinteger(L, i+1);
    lua_pushstring(L, app_args->opt_args[i]);
    lua_rawset(L, -3);
  }
  lua_setglobal(L, "GKYL_COMMANDS");

  // Pushing args without the script name into a new table: this is
  // needed for various tools CLI parsers. For some reason they do not
  // allow the script name in the args
  lua_newtable(L);
  
  for (int i=1; i<app_args->num_opt_args; ++i) {
    lua_pushinteger(L, i);
    lua_pushstring(L, app_args->opt_args[i]);
    lua_rawset(L, -3);
  }
  lua_setglobal(L, "GKYL_COMMANDS_L");

  // set package paths so we find installed libraries
  do {
    const char *fmt = "package.path = package.path .. \";%s/lua/?.lua;%s/lua/Lib/?.lua;%s/lua/?/init.lua\"";    
    size_t len = gkyl_calc_strlen(fmt, app_args->exec_path, app_args->exec_path, app_args->exec_path);

    char *str = gkyl_malloc(len+1);
    snprintf(str, len+1, fmt, app_args->exec_path, app_args->exec_path, app_args->exec_path);
    glua_run_lua(L, str, strlen(str), 0);
    
    gkyl_free(str);
  } while (0);

  // run Lua code (if it exists) before running input file
  if (app_args->echunk)
    glua_run_lua(L, app_args->echunk, strlen(app_args->echunk), 0);

  int rank = 0;
#ifdef GKYL_HAVE_MPI
  if (app_args->use_mpi)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (app_args->num_opt_args > 0) {
    bool something_run = false;

    const char *inp_name = app_args->opt_args[0];
    if (gkyl_check_file_exists(inp_name)) {

      show_banner(rank == 0 ? stdout : 0);

      const char *suff = 0, suff1 = 0;
      if (app_args->use_dir_path) {
        const char *suff = get_fname_with_dir(inp_name);
        const char *suff1 = suff;
        lua_pushlstring(L, suff1, calc_output_prefix_len(suff1));
        lua_setglobal(L, "GKYL_OUT_PREFIX");
        
      }
      else {
        const char *suff = get_fname(inp_name);
        const char *suff1 = suff ? suff+1 : inp_name;
        lua_pushlstring(L, suff1, calc_output_prefix_len(suff1));
        lua_setglobal(L, "GKYL_OUT_PREFIX");
      }
      
      int64_t sz = 0;
      char *buff = gkyl_load_file(inp_name, &sz);
      glua_run_lua(L, buff, sz, stderr);
      gkyl_free(buff);
      something_run = true;
    }
    else {
      // check if it is a Tool, and run it if so
      const char *tlua = get_tool_from_name(app_args->opt_args[0]);
      if (tlua) {
        const char *fmt = "%s/lua/Tool/%s";
        size_t len = gkyl_calc_strlen(fmt, app_args->exec_path, tlua);
        char *tool_name = gkyl_malloc(len+1);
        snprintf(tool_name, len+1, fmt, app_args->exec_path, tlua);
        
        int64_t sz = 0;
        char *buff = gkyl_load_file(tool_name, &sz);
        glua_run_lua(L, buff, sz, stderr);
        gkyl_free(buff);
        
        gkyl_free(tool_name);
        something_run = true;
      }
    }
    if (!something_run)
      fprintf(stderr, "No Lua code was run!\n");
  }
  
  lua_close(L);  

#ifdef GKYL_HAVE_MPI
  if (app_args->use_mpi)
    MPI_Finalize();
#endif

  release_opt_args(app_args);
  
}

#else

#include <stdio.h>

int
main(int argc, char **argv)
{
  fprintf(stderr, "Gkeyll built without Lua support!\n");
  return 0;
}

#endif
