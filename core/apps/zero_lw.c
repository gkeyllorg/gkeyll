#ifdef GKYL_HAVE_LUA

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_format_desc.h>
#include <gkyl_block_topo.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_lua_utils.h>
#include <gkyl_lw_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_zero_lw.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <string.h>

struct gkyl_tool_args *
gkyl_tool_args_new(lua_State *L)
{
  struct gkyl_tool_args *args = gkyl_malloc(sizeof(*args));

  with_lua_global(L, "GKYL_COMMANDS") {
    args->argc = glua_objlen(L);
    args->argv = 0;

    if (args->argc > 0) {
      args->argv = gkyl_malloc(args->argc*sizeof(char*));
    
      for (int i = 1; i <= glua_objlen(L); ++i)  {
        const char *av = glua_tbl_iget_string(L, i, "x");
        args->argv[i-1] = gkyl_malloc(1+strlen(av));
        strcpy(args->argv[i-1], av);
      }
    }
  }

  return args;
}

void
gkyl_tool_args_release(struct gkyl_tool_args* args)
{
  for (int i=0; i<args->argc; ++i)
    gkyl_free(args->argv[i]);
  gkyl_free(args->argv);
  gkyl_free(args);
}

/* *********************/
/* Rect decomp methods */
/* *********************/

// Metatable name for top-level Vlasov App
#define RECT_DECOMP_METATABLE_NM "GkeyllZero.Zero.RectDecomp"

// Lua userdata object for holding Vlasov app and run parameters
struct rect_decomp_lw {
  struct gkyl_rect_decomp *decomp; // Decomp object
};

// G0.RectDecomp.new { cells = { 100, 100}, cuts = { 2, 2 } }
static int
rect_decomp_lw_new(lua_State *L)
{
  struct rect_decomp_lw *rd_lw = gkyl_malloc(sizeof(*rd_lw));
  
  int ndim = 1;
  int cells[GKYL_MAX_DIM] = { 0 }, cuts[GKYL_MAX_DIM] = { 1 };
  
  with_lua_tbl_tbl(L, "cells") {
    ndim = glua_objlen(L);
    for (int d=0; d<ndim; ++d)
      cells[d] = glua_tbl_iget_integer(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "cuts") {
    for (int d=0; d<ndim; ++d)
      cuts[d] = glua_tbl_iget_integer(L, d+1, 0);
  }

  rd_lw->decomp = gkyl_rect_decomp_new_from_cuts_and_cells(ndim, cuts, cells);
  
  // create Lua userdata ...
  struct rect_decomp_lw **l_rd_lw = lua_newuserdata(L, sizeof(struct rect_decomp_lw*));
  *l_rd_lw = rd_lw; // ... point it to the rect decomp pointer

  // set metatable
  luaL_getmetatable(L, RECT_DECOMP_METATABLE_NM);
  lua_setmetatable(L, -2);  
  
  return 1;
}

// Clean up memory allocated for decomp
static int
rect_decomp_lw_gc(lua_State *L)
{
  struct rect_decomp_lw **l_rd_lw = GKYL_CHECK_UDATA(L, RECT_DECOMP_METATABLE_NM);
  struct rect_decomp_lw *rd_lw = *l_rd_lw;

  gkyl_rect_decomp_release(rd_lw->decomp);
  gkyl_free(*l_rd_lw);
  
  return 0;
}

// rect_decomp constructor
static struct luaL_Reg rect_decomp_ctor[] = {
  { "new",  rect_decomp_lw_new },
  { 0, 0 }
};

// rect_decomp methods
static struct luaL_Reg rect_decomp_funcs[] = {
  {0, 0}
};

static void
rect_decomp_openlibs(lua_State *L)
{
  do {
    luaL_newmetatable(L, RECT_DECOMP_METATABLE_NM);

    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, rect_decomp_lw_gc);
    lua_settable(L, -3);

    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, rect_decomp_funcs);
    
    luaL_register(L, "G0.Zero.RectDecomp", rect_decomp_ctor);
    
  } while (0);
}

/* *************************************************************** */
/* Array, RectGrid, and Range wrappers for runregression compareFiles */
/* *************************************************************** */

/* -- RectGrid userdata -- */
#define RECT_GRID_METATABLE_NM "GkeyllZero.Zero.RectGrid"

// Embeds grid by value; no C release function needed (value type).
struct rect_grid_lw {
  struct gkyl_rect_grid grid;
};

static int
rect_grid_lw_gc(lua_State *L)
{
  struct rect_grid_lw **l_g = GKYL_CHECK_UDATA(L, RECT_GRID_METATABLE_NM);
  gkyl_free(*l_g);
  return 0;
}

static struct luaL_Reg rect_grid_funcs[] = {
  { 0, 0 }
};

static void
rect_grid_openlibs(lua_State *L)
{
  luaL_newmetatable(L, RECT_GRID_METATABLE_NM);

  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, rect_grid_lw_gc);
  lua_settable(L, -3);

  lua_pushvalue(L, -1);
  lua_setfield(L, -2, "__index");
  luaL_register(L, NULL, rect_grid_funcs);

  lua_pop(L, 1); // pop metatable
}

/* -- Array userdata -- */
#define ARRAY_METATABLE_NM "GkeyllZero.Zero.Array"

// Holds a heap-allocated gkyl_array*; __gc must release it.
struct array_lw {
  struct gkyl_array *arr;
};

static int
array_lw_gc(lua_State *L)
{
  struct array_lw **l_a = GKYL_CHECK_UDATA(L, ARRAY_METATABLE_NM);
  gkyl_array_release((*l_a)->arr);
  gkyl_free(*l_a);
  return 0;
}

static struct luaL_Reg array_funcs[] = {
  { 0, 0 }
};

static void
array_openlibs(lua_State *L)
{
  luaL_newmetatable(L, ARRAY_METATABLE_NM);

  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, array_lw_gc);
  lua_settable(L, -3);

  lua_pushvalue(L, -1);
  lua_setfield(L, -2, "__index");
  luaL_register(L, NULL, array_funcs);

  lua_pop(L, 1); // pop metatable
}

/* -- Range userdata -- */
#define RANGE_METATABLE_NM "GkeyllZero.Zero.Range"

// Embeds range by value; no C release function needed (value type).
struct range_lw {
  struct gkyl_range range;
};

static int
range_lw_gc(lua_State *L)
{
  struct range_lw **l_r = GKYL_CHECK_UDATA(L, RANGE_METATABLE_NM);
  gkyl_free(*l_r);
  return 0;
}

static struct luaL_Reg range_funcs[] = {
  { 0, 0 }
};

static void
range_openlibs(lua_State *L)
{
  luaL_newmetatable(L, RANGE_METATABLE_NM);

  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, range_lw_gc);
  lua_settable(L, -3);

  lua_pushvalue(L, -1);
  lua_setfield(L, -2, "__index");
  luaL_register(L, NULL, range_funcs);

  lua_pop(L, 1); // pop metatable
}

/* -- Helper: push a rect_grid_lw onto the Lua stack -- */
static void
push_rect_grid(lua_State *L, const struct gkyl_rect_grid *grid)
{
  struct rect_grid_lw *g_lw = gkyl_malloc(sizeof(*g_lw));
  g_lw->grid = *grid;

  struct rect_grid_lw **l_g = lua_newuserdata(L, sizeof(*l_g));
  *l_g = g_lw;
  luaL_getmetatable(L, RECT_GRID_METATABLE_NM);
  lua_setmetatable(L, -2);
}

/* -- Helper: push a range_lw onto the Lua stack -- */
static void
push_range(lua_State *L, const struct gkyl_range *range)
{
  struct range_lw *r_lw = gkyl_malloc(sizeof(*r_lw));
  r_lw->range = *range;

  struct range_lw **l_r = lua_newuserdata(L, sizeof(*l_r));
  *l_r = r_lw;
  luaL_getmetatable(L, RANGE_METATABLE_NM);
  lua_setmetatable(L, -2);
}

/* -- G0.Zero.gkylFileType(fname) -> string -- */
// Returns one of: "field", "dynvector", "multi-range-field",
//                 "block-topology", "multi-block-meta", or "not-gkyl".
static int
gkyl_file_type_lw(lua_State *L)
{
  const char *fname = luaL_checkstring(L, 1);
  int ftype = gkyl_get_gkyl_file_type(fname);

  // gkyl_get_gkyl_file_type returns the raw file-header integer (1-based):
  //   1 = field, 2 = dynvector, 3 = multi-range-field,
  //   4 = block-topology, 5 = multi-block-meta
  // names[] is indexed directly by ftype so names[0] is a placeholder.
  static const char *names[] = {
    NULL,                // [0] — not a valid file type
    "field",             // [1] — GKYL_FIELD_DATA_FILE
    "dynvector",         // [2] — GKYL_DYNVEC_DATA_FILE
    "multi-range-field", // [3] — GKYL_MULTI_RANGE_DATA_FILE
    "block-topology",    // [4] — GKYL_BLOCK_TOPO_DATA_FILE
    "multi-block-meta"   // [5] — GKYL_MULTI_BLOCK_DATA_FILE
  };
  int nnames = (int)(sizeof(names) / sizeof(names[0]));

  if (ftype <= 0 || ftype >= nnames || names[ftype] == NULL)
    lua_pushstring(L, "not-gkyl");
  else
    lua_pushstring(L, names[ftype]);

  return 1;
}

/* -- G0.Zero.arrayNewFromFile(fname) -> grid_ud, array_ud  (or nil, nil on failure) -- */
static int
array_new_from_file_lw(lua_State *L)
{
  const char *fname = luaL_checkstring(L, 1);

  struct rect_grid_lw *g_lw = gkyl_malloc(sizeof(*g_lw));
  struct gkyl_array *arr = gkyl_grid_array_new_from_file(&g_lw->grid, fname);

  if (!arr) {
    gkyl_free(g_lw);
    lua_pushnil(L);
    lua_pushnil(L);
    return 2;
  }

  // push grid userdata
  push_rect_grid(L, &g_lw->grid);
  gkyl_free(g_lw); // grid was copied into the userdata by push_rect_grid

  // push array userdata
  struct array_lw *a_lw = gkyl_malloc(sizeof(*a_lw));
  a_lw->arr = arr;
  struct array_lw **l_a = lua_newuserdata(L, sizeof(*l_a));
  *l_a = a_lw;
  luaL_getmetatable(L, ARRAY_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 2; // (grid_ud, array_ud)
}

/* -- G0.Zero.rectGridCmp(g1, g2) -> bool -- */
static int
rect_grid_cmp_lw(lua_State *L)
{
  struct rect_grid_lw **l_g1 = luaL_checkudata(L, 1, RECT_GRID_METATABLE_NM);
  struct rect_grid_lw **l_g2 = luaL_checkudata(L, 2, RECT_GRID_METATABLE_NM);
  lua_pushboolean(L, gkyl_rect_grid_cmp(&(*l_g1)->grid, &(*l_g2)->grid));
  return 1;
}

/* -- G0.Zero.createGridRanges(grid_ud, nghost_table) -> range_ud, ext_range_ud -- */
static int
create_grid_ranges_lw(lua_State *L)
{
  struct rect_grid_lw **l_g = luaL_checkudata(L, 1, RECT_GRID_METATABLE_NM);
  luaL_checktype(L, 2, LUA_TTABLE);

  int ndim = (*l_g)->grid.ndim;
  int nghost[GKYL_MAX_DIM] = { 0 };
  for (int d = 0; d < ndim; ++d) {
    lua_rawgeti(L, 2, d+1);
    nghost[d] = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);
  }

  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&(*l_g)->grid, nghost, &ext_range, &range);

  push_range(L, &range);
  push_range(L, &ext_range);

  return 2; // (range_ud, ext_range_ud)
}

/* -- G0.Zero.arrayDiff(a1_ud, a2_ud, range_ud) -> table -- */
// Returned table has fields: is_compatible, max_abs_diff, min_abs_diff,
//                             max_rel_diff, min_rel_diff.
static int
array_diff_lw(lua_State *L)
{
  struct array_lw **l_a1 = luaL_checkudata(L, 1, ARRAY_METATABLE_NM);
  struct array_lw **l_a2 = luaL_checkudata(L, 2, ARRAY_METATABLE_NM);
  struct range_lw **l_r  = luaL_checkudata(L, 3, RANGE_METATABLE_NM);

  struct gkyl_array_diff diff = gkyl_array_diff(
    (*l_a1)->arr, (*l_a2)->arr, &(*l_r)->range);

  lua_newtable(L);

  lua_pushboolean(L, diff.is_compatible);
  lua_setfield(L, -2, "is_compatible");

  lua_pushnumber(L, diff.max_abs_diff);
  lua_setfield(L, -2, "max_abs_diff");

  lua_pushnumber(L, diff.min_abs_diff);
  lua_setfield(L, -2, "min_abs_diff");

  lua_pushnumber(L, diff.max_rel_diff);
  lua_setfield(L, -2, "max_rel_diff");

  lua_pushnumber(L, diff.min_rel_diff);
  lua_setfield(L, -2, "min_rel_diff");

  return 1;
}

// Returns a Lua table with all diff fields set to zero and is_compatible=false.
// Used as the early-exit return value whenever two dynvec files cannot be compared.
static int
dynvec_incompat(lua_State *L)
{
  lua_newtable(L);
  lua_pushboolean(L, 0); lua_setfield(L, -2, "is_compatible");
  lua_pushnumber(L, 0);  lua_setfield(L, -2, "max_abs_diff");
  lua_pushnumber(L, 0);  lua_setfield(L, -2, "min_abs_diff");
  lua_pushnumber(L, 0);  lua_setfield(L, -2, "max_rel_diff");
  lua_pushnumber(L, 0);  lua_setfield(L, -2, "min_rel_diff");
  lua_pushnumber(L, 0);  lua_setfield(L, -2, "tm_max_abs_diff");
  lua_pushnumber(L, 0);  lua_setfield(L, -2, "tm_min_abs_diff");
  return 1;
}

// G0.Zero.dynvecDiff(f1, f2) -> table
// Reads two dynvector files, compares their data payloads, and returns a
// table with the same fields as arrayDiff plus two informational timestamp
// fields (tm_max_abs_diff, tm_min_abs_diff).  Timestamp differences are
// computed but do not affect is_compatible: accumulated floating-point
// drift in t += dt is expected and platform-dependent.
static int
dynvec_diff_lw(lua_State *L)
{
  const char *f1 = luaL_checkstring(L, 1);
  const char *f2 = luaL_checkstring(L, 2);

  // Read metadata (fast: header only, no data allocation).
  struct gkyl_dynvec_etype_ncomp enc1 = gkyl_dynvec_read_ncomp(f1);
  struct gkyl_dynvec_etype_ncomp enc2 = gkyl_dynvec_read_ncomp(f2);

  // ncomp == 0 means the file couldn't be read or has an invalid header.
  if (enc1.ncomp == 0 || enc2.ncomp == 0) return dynvec_incompat(L);

  // gkyl_array_diff only handles GKYL_DOUBLE safely.
  if (enc1.type != GKYL_DOUBLE || enc2.type != GKYL_DOUBLE) return dynvec_incompat(L);

  // Number of components must match.
  if (enc1.ncomp != enc2.ncomp) return dynvec_incompat(L);

  // Read full dynvecs.
  gkyl_dynvec dv1 = gkyl_dynvec_new(enc1.type, enc1.ncomp);
  gkyl_dynvec dv2 = gkyl_dynvec_new(enc2.type, enc2.ncomp);

  bool ok1 = gkyl_dynvec_read(dv1, f1);
  bool ok2 = gkyl_dynvec_read(dv2, f2);

  if (!ok1 || !ok2) {
    gkyl_dynvec_release(dv1);
    gkyl_dynvec_release(dv2);
    return dynvec_incompat(L);
  }

  size_t n1 = gkyl_dynvec_size(dv1);
  size_t n2 = gkyl_dynvec_size(dv2);

  // Different number of time steps -> incompatible.
  if (n1 != n2) {
    gkyl_dynvec_release(dv1);
    gkyl_dynvec_release(dv2);
    return dynvec_incompat(L);
  }

  size_t nsteps = n1;

  // Both dynvecs are empty (no time steps recorded): trivially equal.
  // Guard required because gkyl_array_new with nelem=0 may return an array
  // whose data pointer is NULL, which then causes gkyl_aligned_free_ to
  // assert on the NULL pointer during gkyl_array_release.
  if (nsteps == 0) {
    gkyl_dynvec_release(dv1);
    gkyl_dynvec_release(dv2);
    lua_newtable(L);
    lua_pushboolean(L, 1);  lua_setfield(L, -2, "is_compatible");
    lua_pushnumber(L, 0.0); lua_setfield(L, -2, "max_abs_diff");
    lua_pushnumber(L, 0.0); lua_setfield(L, -2, "min_abs_diff");
    lua_pushnumber(L, 0.0); lua_setfield(L, -2, "max_rel_diff");
    lua_pushnumber(L, 0.0); lua_setfield(L, -2, "min_rel_diff");
    lua_pushnumber(L, 0.0); lua_setfield(L, -2, "tm_max_abs_diff");
    lua_pushnumber(L, 0.0); lua_setfield(L, -2, "tm_min_abs_diff");
    return 1;
  }

  // Allocate flat arrays for timestamps (ncomp=1) and data (ncomp=enc1.ncomp).
  struct gkyl_array *tm1 = gkyl_array_new(GKYL_DOUBLE, 1,          nsteps);
  struct gkyl_array *tm2 = gkyl_array_new(GKYL_DOUBLE, 1,          nsteps);
  struct gkyl_array *da1 = gkyl_array_new(GKYL_DOUBLE, enc1.ncomp, nsteps);
  struct gkyl_array *da2 = gkyl_array_new(GKYL_DOUBLE, enc2.ncomp, nsteps);

  gkyl_dynvec_to_array(dv1, tm1, da1);
  gkyl_dynvec_to_array(dv2, tm2, da2);

  gkyl_dynvec_release(dv1);
  gkyl_dynvec_release(dv2);

  // Build a 1-D range [0, nsteps-1] for gkyl_array_diff.
  struct gkyl_range rng;
  int shape[1] = { (int)nsteps };
  gkyl_range_init_from_shape(&rng, 1, shape);

  struct gkyl_array_diff ddiff = gkyl_array_diff(da1, da2, &rng);
  struct gkyl_array_diff tdiff = gkyl_array_diff(tm1, tm2, &rng);

  gkyl_array_release(tm1);
  gkyl_array_release(tm2);
  gkyl_array_release(da1);
  gkyl_array_release(da2);

  lua_newtable(L);

  lua_pushboolean(L, ddiff.is_compatible);
  lua_setfield(L, -2, "is_compatible");

  lua_pushnumber(L, ddiff.max_abs_diff);
  lua_setfield(L, -2, "max_abs_diff");

  lua_pushnumber(L, ddiff.min_abs_diff);
  lua_setfield(L, -2, "min_abs_diff");

  lua_pushnumber(L, ddiff.max_rel_diff);
  lua_setfield(L, -2, "max_rel_diff");

  lua_pushnumber(L, ddiff.min_rel_diff);
  lua_setfield(L, -2, "min_rel_diff");

  // Timestamp fields: informational only, never cause a test failure.
  lua_pushnumber(L, tdiff.is_compatible ? tdiff.max_abs_diff : 0.0);
  lua_setfield(L, -2, "tm_max_abs_diff");

  lua_pushnumber(L, tdiff.is_compatible ? tdiff.min_abs_diff : 0.0);
  lua_setfield(L, -2, "tm_min_abs_diff");

  return 1;
}

// G0.Zero.blockTopoCmp(f1, f2) -> bool
// Reads two block-topology files and returns true if they are identical
// (same ndim, num_blocks, and all per-block connection entries), false
// otherwise.  Block topology contains only integer/enum data, so this
// is an exact equality check rather than a numerical diff.
static int
block_topo_cmp_lw(lua_State *L)
{
  const char *f1 = luaL_checkstring(L, 1);
  const char *f2 = luaL_checkstring(L, 2);

  int st1 = 0, st2 = 0;
  struct gkyl_block_topo *bt1 = gkyl_block_topo_read(f1, &st1);
  struct gkyl_block_topo *bt2 = gkyl_block_topo_read(f2, &st2);

  if (!bt1 || !bt2) {
    if (bt1) gkyl_block_topo_release(bt1);
    if (bt2) gkyl_block_topo_release(bt2);
    lua_pushboolean(L, 0);
    return 1;
  }

  bool equal = (bt1->ndim == bt2->ndim) && (bt1->num_blocks == bt2->num_blocks);

  if (equal) {
    int ndim      = bt1->ndim;
    int nblocks   = bt1->num_blocks;
    for (int b = 0; b < nblocks && equal; ++b) {
      for (int d = 0; d < ndim && equal; ++d) {
        for (int e = 0; e < 2 && equal; ++e) {
          struct gkyl_target_edge *c1 = &bt1->conn[b].connections[d][e];
          struct gkyl_target_edge *c2 = &bt2->conn[b].connections[d][e];
          if (c1->bid != c2->bid || c1->dir != c2->dir || c1->edge != c2->edge)
            equal = false;
        }
      }
    }
  }

  gkyl_block_topo_release(bt1);
  gkyl_block_topo_release(bt2);

  lua_pushboolean(L, equal ? 1 : 0);
  return 1;
}

// Module-level functions registered under G0.Zero
static struct luaL_Reg zero_array_funcs[] = {
  { "gkylFileType",       gkyl_file_type_lw       },
  { "arrayNewFromFile",   array_new_from_file_lw  },
  { "rectGridCmp",        rect_grid_cmp_lw        },
  { "createGridRanges",   create_grid_ranges_lw   },
  { "arrayDiff",          array_diff_lw           },
  { "dynvecDiff",         dynvec_diff_lw          },
  { "blockTopoCmp",       block_topo_cmp_lw       },
  { 0, 0 }
};

void
gkyl_zero_lw_openlibs(lua_State *L)
{
  // Push empty global table called "G0".
  lua_newtable(L);
  lua_setglobal(L, "G0");

  rect_decomp_openlibs(L);

  // Register metatables for array comparison userdata types.
  rect_grid_openlibs(L);
  array_openlibs(L);
  range_openlibs(L);

  // Register G0.Zero module-level functions.
  luaL_register(L, "G0.Zero", zero_array_funcs);

  // Register types for species, field, and Poisson solver boundary conditions.
  gkyl_register_distribution_moment_types(L);
  gkyl_register_species_bc_types(L);
  gkyl_register_field_bc_types(L);
}

#endif
