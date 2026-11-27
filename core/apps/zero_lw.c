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

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

// The core is built with -ffast-math, which can optimize isfinite away.
static bool
finite_double(double value)
{
  volatile uint64_t bits;
  memcpy((void *)&bits, &value, sizeof bits);
  return (bits & UINT64_C(0x7ff0000000000000)) != UINT64_C(0x7ff0000000000000);
}

static double
infinite_double(void)
{
  // Construct infinity from its bits so fast-math cannot fold it away.
  uint64_t bits = UINT64_C(0x7ff0000000000000);
  double value;
  memcpy(&value, &bits, sizeof value);
  return value;
}

struct gkyl_tool_args *
gkyl_tool_args_new(lua_State *L)
{
  struct gkyl_tool_args *args = gkyl_malloc(sizeof(*args));

  with_lua_global(L, "GKYL_COMMANDS")
  {
    args->argc = glua_objlen(L);
    args->argv = 0;

    if (args->argc > 0) {
      args->argv = gkyl_malloc(args->argc * sizeof(char *));

      for (int i = 1; i <= glua_objlen(L); ++i) {
        const char *av = glua_tbl_iget_string(L, i, "x");
        args->argv[i - 1] = gkyl_malloc(1 + strlen(av));
        strcpy(args->argv[i - 1], av);
      }
    }
  }

  return args;
}

void
gkyl_tool_args_release(struct gkyl_tool_args *args)
{
  for (int i = 0; i < args->argc; ++i) {
    gkyl_free(args->argv[i]);
  }
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
  int cells[GKYL_MAX_DIM] = {0}, cuts[GKYL_MAX_DIM] = {1};

  with_lua_tbl_tbl(L, "cells")
  {
    ndim = glua_objlen(L);
    for (int d = 0; d < ndim; ++d) {
      cells[d] = glua_tbl_iget_integer(L, d + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "cuts")
  {
    for (int d = 0; d < ndim; ++d) {
      cuts[d] = glua_tbl_iget_integer(L, d + 1, 0);
    }
  }

  rd_lw->decomp = gkyl_rect_decomp_new_from_cuts_and_cells(ndim, cuts, cells);

  // create Lua userdata ...
  struct rect_decomp_lw **l_rd_lw = lua_newuserdata(L, sizeof(struct rect_decomp_lw *));
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
static struct luaL_Reg rect_decomp_ctor[] = {{"new", rect_decomp_lw_new}, {0, 0}};

// rect_decomp methods
static struct luaL_Reg rect_decomp_funcs[] = {{0, 0}};

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

static struct luaL_Reg rect_grid_funcs[] = {{0, 0}};

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

static struct luaL_Reg array_funcs[] = {{0, 0}};

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

static struct luaL_Reg range_funcs[] = {{0, 0}};

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
    NULL, // [0] — not a valid file type
    "field", // [1] — GKYL_FIELD_DATA_FILE
    "dynvector", // [2] — GKYL_DYNVEC_DATA_FILE
    "multi-range-field", // [3] — GKYL_MULTI_RANGE_DATA_FILE
    "block-topology", // [4] — GKYL_BLOCK_TOPO_DATA_FILE
    "multi-block-meta" // [5] — GKYL_MULTI_BLOCK_DATA_FILE
  };
  int nnames = (int)(sizeof(names) / sizeof(names[0]));

  if (ftype <= 0 || ftype >= nnames || names[ftype] == NULL) {
    lua_pushstring(L, "not-gkyl");
  } else {
    lua_pushstring(L, names[ftype]);
  }

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
  int nghost[GKYL_MAX_DIM] = {0};
  for (int d = 0; d < ndim; ++d) {
    lua_rawgeti(L, 2, d + 1);
    nghost[d] = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);
  }

  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&(*l_g)->grid, nghost, &ext_range, &range);

  push_range(L, &range);
  push_range(L, &ext_range);

  return 2; // (range_ud, ext_range_ud)
}

// Apply absolute and relative tolerances to each coefficient together.
static bool
array_within_combined_tol(
  const struct gkyl_array *a, const struct gkyl_array *b, const struct gkyl_range *range,
  double abs_tol, double rel_tol
)
{
  // Each coefficient passes if either its absolute or relative error is small.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);
    const double *av = gkyl_array_cfetch(a, loc), *bv = gkyl_array_cfetch(b, loc);
    for (int c = 0; c < a->ncomp; ++c) {
      double diff = fabs(av[c] - bv[c]);
      if (!finite_double(av[c]) || !finite_double(bv[c]) || !finite_double(diff)) {
        return false;
      }
      if (diff <= abs_tol) {
        continue;
      }
      double sum = fabs(av[c]) + fabs(bv[c]);
      double rel_diff = av[c] == 0 || bv[c] == 0 || sum < DBL_MIN ? diff :
                                                                    diff / fmin(sum, DBL_MAX);
      if (rel_diff > rel_tol) {
        return false;
      }
    }
  }
  return true;
}

/* -- G0.Zero.arrayDiff(a1_ud, a2_ud, range_ud, abs_tol?, rel_tol?) -> table -- */
// Returned table has fields: is_compatible, max_abs_diff, min_abs_diff,
//                             max_rel_diff, min_rel_diff, within_combined_tol.
static int
array_diff_lw(lua_State *L)
{
  struct array_lw **l_a1 = luaL_checkudata(L, 1, ARRAY_METATABLE_NM);
  struct array_lw **l_a2 = luaL_checkudata(L, 2, ARRAY_METATABLE_NM);
  struct range_lw **l_r = luaL_checkudata(L, 3, RANGE_METATABLE_NM);

  double abs_tol = luaL_optnumber(L, 4, 1e-12);
  double rel_tol = luaL_optnumber(L, 5, 1e-12);
  luaL_argcheck(L, finite_double(abs_tol) && abs_tol >= 0, 4, "invalid absolute tolerance");
  luaL_argcheck(L, finite_double(rel_tol) && rel_tol >= 0, 5, "invalid relative tolerance");
  struct gkyl_array_diff diff = gkyl_array_diff((*l_a1)->arr, (*l_a2)->arr, &(*l_r)->range);
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

  if (!lua_isnoneornil(L, 4) && !lua_isnoneornil(L, 5)) {
    bool passed = diff.is_compatible && array_within_combined_tol(
                                          (*l_a1)->arr, (*l_a2)->arr, &(*l_r)->range,
                                          luaL_checknumber(L, 4), luaL_checknumber(L, 5)
                                        );
    lua_pushboolean(L, passed);
    lua_setfield(L, -2, "within_combined_tol");
  }

  return 1;
}

// Compare a selected gyrokinetic moment using baseline-only physical scales.
// nil means the required companion data cannot establish a valid scale.
static int
scaled_moment_diff_lw(lua_State *L)
{
  struct array_lw **bl = luaL_checkudata(L, 1, ARRAY_METATABLE_NM);
  struct array_lw **cl = luaL_checkudata(L, 2, ARRAY_METATABLE_NM);
  struct gkyl_array *b = (*bl)->arr, *c = (*cl)->arr;
  const char *kind = luaL_checkstring(L, 3);
  double tol = luaL_checknumber(L, 4);
  struct gkyl_array *m0 = NULL, *mp = NULL, *mt = NULL;
  if (!lua_isnoneornil(L, 5)) {
    m0 = (*(struct array_lw **)luaL_checkudata(L, 5, ARRAY_METATABLE_NM))->arr;
  }
  if (!lua_isnoneornil(L, 6)) {
    mp = (*(struct array_lw **)luaL_checkudata(L, 6, ARRAY_METATABLE_NM))->arr;
  }
  if (!lua_isnoneornil(L, 7)) {
    mt = (*(struct array_lw **)luaL_checkudata(L, 7, ARRAY_METATABLE_NM))->arr;
  }
  if (b->type != GKYL_DOUBLE || c->type != GKYL_DOUBLE || b->ncomp == 0 || b->ncomp != c->ncomp ||
      b->size != c->size) {
    lua_pushnil(L);
    return 1;
  }
  int ncomp = b->ncomp, offset = 0, nb = ncomp;
  // Maxwellian files pack several moments; compare only their drift block here.
  bool maxwell = strcmp(kind, "Maxwellian") == 0, bimax = strcmp(kind, "BiMaxwellian") == 0;
  if (maxwell || bimax) {
    int nmom = bimax ? 4 : 3;
    if (ncomp % nmom) {
      lua_pushnil(L);
      return 1;
    }
    nb = ncomp / nmom;
    offset = nb;
  } else if (!m0 || !mp || m0->type != GKYL_DOUBLE || mp->type != GKYL_DOUBLE || m0->ncomp == 0 ||
             mp->ncomp == 0 || m0->size != b->size || mp->size != b->size ||
             (strcmp(kind, "M3perp") == 0 &&
              (!mt || mt->type != GKYL_DOUBLE || mt->ncomp == 0 || mt->size != b->size))) {
    lua_pushnil(L);
    return 1;
  }
  // Use baseline moments over the whole file to set one physical scale.
  double n = 0, p = 0, t = 0;
  for (long i = 0; i < b->size; ++i) {
    if (maxwell || bimax) {
      const double *v = gkyl_array_cfetch(b, i);
      if (!finite_double(v[2 * nb])) {
        lua_pushnil(L);
        return 1;
      }
      t = fmax(t, v[2 * nb]);
    } else {
      double n0 = *(const double *)gkyl_array_cfetch(m0, i);
      double p0 = *(const double *)gkyl_array_cfetch(mp, i);
      if (!finite_double(n0) || !finite_double(p0)) {
        lua_pushnil(L);
        return 1;
      }
      n = fmax(n, fabs(n0));
      p = fmax(p, fabs(p0));
      if (mt) {
        double t0 = *(const double *)gkyl_array_cfetch(mt, i);
        if (!finite_double(t0)) {
          lua_pushnil(L);
          return 1;
        }
        t = fmax(t, fabs(t0));
      }
    }
  }
  double scale = 0;
  if (maxwell || bimax) {
    struct rect_grid_lw **grid = luaL_checkudata(L, 8, RECT_GRID_METATABLE_NM);
    int ndim = (*grid)->grid.ndim;
    double basis_zero = pow(2.0, ndim / 2.0);
    scale = basis_zero * sqrt(t / basis_zero);
  } else if (n > 0 && p > 0) {
    if (strcmp(kind, "M1") == 0) {
      scale = sqrt(n) * sqrt(p);
    }
    if (strcmp(kind, "M3par") == 0) {
      scale = p * sqrt(p / n);
    }
    if (strcmp(kind, "M3perp") == 0) {
      scale = t * sqrt(p / n);
    }
  }
  if (!(scale > 0) || !finite_double(scale)) {
    lua_pushnil(L);
    return 1;
  }
  // Check the mean coefficient and the combined higher-mode error separately.
  double worst = 0;
  bool passed = true;
  for (long i = 0; i < b->size; ++i) {
    const double *bv = gkyl_array_cfetch(b, i), *cv = gkyl_array_cfetch(c, i);
    double high_sq = 0;
    for (int k = 0; k < nb; ++k) {
      double diff = fabs(cv[offset + k] - bv[offset + k]);
      if (!finite_double(bv[offset + k]) || !finite_double(cv[offset + k]) ||
          !finite_double(diff)) {
        lua_pushboolean(L, 0);
        lua_pushnumber(L, infinite_double());
        return 2;
      }
      if (k == 0) {
        double ratio = diff / fmax(fabs(bv[offset]), scale);
        worst = fmax(worst, ratio);
        if (ratio > tol) {
          passed = false;
        }
      } else {
        high_sq += (diff / scale) * (diff / scale);
      }
    }
    double high_ratio = sqrt(high_sq);
    worst = fmax(worst, high_ratio);
    if (high_ratio > tol) {
      passed = false;
    }
    if (maxwell || bimax) {
      for (int k = 0; k < ncomp; ++k) {
        if (k >= offset && k < offset + nb) {
          continue;
        }
        double diff = fabs(cv[k] - bv[k]);
        if (!finite_double(diff) || !finite_double(cv[k]) || !finite_double(bv[k]) ||
            (diff > tol && diff > tol * fabs(bv[k]))) {
          passed = false;
        }
      }
    }
  }
  lua_pushboolean(L, passed);
  lua_pushnumber(L, worst);
  return 2;
}

// Compare Vlasov vector moments against a global baseline density/thermal scale.
static int
scaled_vlasov_moment_diff_lw(lua_State *L)
{
  struct gkyl_array *b = (*(struct array_lw **)luaL_checkudata(L, 1, ARRAY_METATABLE_NM))->arr;
  struct gkyl_array *c = (*(struct array_lw **)luaL_checkudata(L, 2, ARRAY_METATABLE_NM))->arr;
  const char *kind = luaL_checkstring(L, 3);
  double tol = luaL_checknumber(L, 4), abs_tol = luaL_checknumber(L, 5);
  struct gkyl_array *m0 = (*(struct array_lw **)luaL_checkudata(L, 6, ARRAY_METATABLE_NM))->arr;
  struct gkyl_array *energy =
    lua_isnoneornil(L, 7) ? NULL :
                            (*(struct array_lw **)luaL_checkudata(L, 7, ARRAY_METATABLE_NM))->arr;
  struct gkyl_array *lte = lua_isnoneornil(L, 8) ?
                             NULL :
                             (*(struct array_lw **)luaL_checkudata(L, 8, ARRAY_METATABLE_NM))->arr;
  struct rect_grid_lw **grid = luaL_checkudata(L, 9, RECT_GRID_METATABLE_NM);
  bool is_lte = strcmp(kind, "LTE") == 0;
  if ((!is_lte && strcmp(kind, "M1") != 0) || b->type != GKYL_DOUBLE || c->type != GKYL_DOUBLE ||
      m0->type != GKYL_DOUBLE || b->size != c->size || b->size != m0->size ||
      b->ncomp != c->ncomp || m0->ncomp < 1) {
    lua_pushnil(L);
    return 1;
  }
  int nb = m0->ncomp, vdim = b->ncomp / nb - (is_lte ? 2 : 0);
  if (vdim < 1 || vdim > 3 || b->ncomp != (vdim + (is_lte ? 2 : 0)) * nb ||
      (energy && (energy->type != GKYL_DOUBLE || energy->size != b->size || energy->ncomp != nb)) ||
      (lte && (lte->type != GKYL_DOUBLE || lte->size != b->size || lte->ncomp != (vdim + 2) * nb)
      ) ||
      (!is_lte && !energy && !lte)) {
    lua_pushnil(L);
    return 1;
  }
  // Extract a global density and thermal scale from baseline diagnostics.
  double n = 0, p = 0, t = 0;
  for (long i = 0; i < b->size; ++i) {
    double n0 = *(const double *)gkyl_array_cfetch(m0, i);
    if (!finite_double(n0)) {
      lua_pushnil(L);
      return 1;
    }
    n = fmax(n, fabs(n0));
    if (energy) {
      double e0 = *(const double *)gkyl_array_cfetch(energy, i);
      if (!finite_double(e0)) {
        lua_pushnil(L);
        return 1;
      }
      p = fmax(p, e0);
    }
    if (lte || is_lte) {
      const double *lv = gkyl_array_cfetch(is_lte ? b : lte, i);
      double t0 = lv[(vdim + 1) * nb];
      if (!finite_double(t0)) {
        lua_pushnil(L);
        return 1;
      }
      t = fmax(t, t0);
    }
  }
  double basis_zero = pow(2.0, (*grid)->grid.ndim / 2.0);
  double scale = is_lte ? sqrt(basis_zero * t) : (energy ? sqrt(n * p) : n * sqrt(t / basis_zero));
  if (!(scale > 0) || !finite_double(scale)) {
    lua_pushnil(L);
    return 1;
  }
  // LTE stores density and temperature around its drift-velocity blocks.
  bool passed = true;
  double worst = 0;
  for (long i = 0; i < b->size; ++i) {
    const double *bv = gkyl_array_cfetch(b, i), *cv = gkyl_array_cfetch(c, i);
    for (int d = 0; d < vdim; ++d) {
      int offset = (d + (is_lte ? 1 : 0)) * nb;
      double high_sq = 0;
      for (int k = 0; k < nb; ++k) {
        double diff = fabs(cv[offset + k] - bv[offset + k]);
        if (!finite_double(bv[offset + k]) || !finite_double(cv[offset + k]) ||
            !finite_double(diff)) {
          lua_pushboolean(L, 0);
          lua_pushnumber(L, infinite_double());
          return 2;
        }
        if (k == 0) {
          double ratio = diff / fmax(fabs(bv[offset]), scale);
          worst = fmax(worst, ratio);
          if (ratio > tol) {
            passed = false;
          }
        } else {
          high_sq += (diff / scale) * (diff / scale);
        }
      }
      double ratio = sqrt(high_sq);
      worst = fmax(worst, ratio);
      if (ratio > tol) {
        passed = false;
      }
    }
    if (is_lte) {
      for (int block = 0; block <= vdim + 1; block += vdim + 1) {
        for (int k = 0; k < nb; ++k) {
          int j = block * nb + k;
          double diff = fabs(cv[j] - bv[j]);
          if (!finite_double(bv[j]) || !finite_double(cv[j]) || !finite_double(diff)) {
            lua_pushboolean(L, 0);
            lua_pushnumber(L, infinite_double());
            return 2;
          }
          if (diff > abs_tol && diff > tol * (fabs(bv[j]) + fabs(cv[j]))) {
            passed = false;
          }
        }
      }
    }
  }
  lua_pushboolean(L, passed);
  lua_pushnumber(L, worst);
  return 2;
}

// Reject nonfinite coefficients before a tolerance or scale-aware comparison.
static int
array_finite_lw(lua_State *L)
{
  struct array_lw **a = luaL_checkudata(L, 1, ARRAY_METATABLE_NM);
  if ((*a)->arr->type != GKYL_DOUBLE) {
    lua_pushboolean(L, 0);
    return 1;
  }
  for (long i = 0; i < (*a)->arr->size; ++i) {
    const double *v = gkyl_array_cfetch((*a)->arr, i);
    for (int k = 0; k < (*a)->arr->ncomp; ++k) {
      if (!finite_double(v[k])) {
        lua_pushboolean(L, 0);
        return 1;
      }
    }
  }
  lua_pushboolean(L, 1);
  return 1;
}

// Returns a Lua table with all diff fields set to zero and is_compatible=false.
// Used as the early-exit return value whenever two dynvec files cannot be compared.
static int
dynvec_incompat(lua_State *L)
{
  lua_newtable(L);
  lua_pushboolean(L, 0);
  lua_setfield(L, -2, "is_compatible");
  lua_pushnumber(L, 0);
  lua_setfield(L, -2, "max_abs_diff");
  lua_pushnumber(L, 0);
  lua_setfield(L, -2, "min_abs_diff");
  lua_pushnumber(L, 0);
  lua_setfield(L, -2, "max_rel_diff");
  lua_pushnumber(L, 0);
  lua_setfield(L, -2, "min_rel_diff");
  lua_pushnumber(L, 0);
  lua_setfield(L, -2, "tm_max_abs_diff");
  lua_pushnumber(L, 0);
  lua_setfield(L, -2, "tm_min_abs_diff");
  return 1;
}

// G0.Zero.dynvecDiff(f1, f2, abs_tol?, rel_tol?) -> table
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
  double abs_tol = luaL_optnumber(L, 3, 1e-12);
  double rel_tol = luaL_optnumber(L, 4, 1e-12);
  luaL_argcheck(L, finite_double(abs_tol) && abs_tol >= 0, 3, "invalid absolute tolerance");
  luaL_argcheck(L, finite_double(rel_tol) && rel_tol >= 0, 4, "invalid relative tolerance");

  // Read metadata (fast: header only, no data allocation).
  struct gkyl_dynvec_etype_ncomp enc1 = gkyl_dynvec_read_ncomp(f1);
  struct gkyl_dynvec_etype_ncomp enc2 = gkyl_dynvec_read_ncomp(f2);

  // ncomp == 0 means the file couldn't be read or has an invalid header.
  if (enc1.ncomp == 0 || enc2.ncomp == 0) {
    return dynvec_incompat(L);
  }

  // gkyl_array_diff only handles GKYL_DOUBLE safely.
  if (enc1.type != GKYL_DOUBLE || enc2.type != GKYL_DOUBLE) {
    return dynvec_incompat(L);
  }

  // Number of components must match.
  if (enc1.ncomp != enc2.ncomp) {
    return dynvec_incompat(L);
  }

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
    lua_pushboolean(L, 1);
    lua_setfield(L, -2, "is_compatible");
    lua_pushnumber(L, 0.0);
    lua_setfield(L, -2, "max_abs_diff");
    lua_pushnumber(L, 0.0);
    lua_setfield(L, -2, "min_abs_diff");
    lua_pushnumber(L, 0.0);
    lua_setfield(L, -2, "max_rel_diff");
    lua_pushnumber(L, 0.0);
    lua_setfield(L, -2, "min_rel_diff");
    lua_pushnumber(L, 0.0);
    lua_setfield(L, -2, "tm_max_abs_diff");
    lua_pushnumber(L, 0.0);
    lua_setfield(L, -2, "tm_min_abs_diff");
    if (!lua_isnoneornil(L, 3) && !lua_isnoneornil(L, 4)) {
      lua_pushboolean(L, 1);
      lua_setfield(L, -2, "within_combined_tol");
    }
    return 1;
  }

  // Allocate flat arrays for timestamps (ncomp=1) and data (ncomp=enc1.ncomp).
  struct gkyl_array *tm1 = gkyl_array_new(GKYL_DOUBLE, 1, nsteps);
  struct gkyl_array *tm2 = gkyl_array_new(GKYL_DOUBLE, 1, nsteps);
  struct gkyl_array *da1 = gkyl_array_new(GKYL_DOUBLE, enc1.ncomp, nsteps);
  struct gkyl_array *da2 = gkyl_array_new(GKYL_DOUBLE, enc2.ncomp, nsteps);

  gkyl_dynvec_to_array(dv1, tm1, da1);
  gkyl_dynvec_to_array(dv2, tm2, da2);

  gkyl_dynvec_release(dv1);
  gkyl_dynvec_release(dv2);

  // Build a 1-D range [0, nsteps-1] for gkyl_array_diff.
  struct gkyl_range rng;
  int shape[1] = {(int)nsteps};
  gkyl_range_init_from_shape(&rng, 1, shape);

  struct gkyl_array_diff ddiff = gkyl_array_diff(da1, da2, &rng);
  struct gkyl_array_diff tdiff = gkyl_array_diff(tm1, tm2, &rng);
  bool check_tol = !lua_isnoneornil(L, 3) && !lua_isnoneornil(L, 4);
  bool passed =
    check_tol && ddiff.is_compatible &&
    array_within_combined_tol(da1, da2, &rng, luaL_checknumber(L, 3), luaL_checknumber(L, 4));

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

  if (check_tol) {
    lua_pushboolean(L, passed);
    lua_setfield(L, -2, "within_combined_tol");
  }

  // Timestamp fields: informational only, never cause a test failure.
  lua_pushnumber(L, tdiff.is_compatible ? tdiff.max_abs_diff : 0.0);
  lua_setfield(L, -2, "tm_max_abs_diff");

  lua_pushnumber(L, tdiff.is_compatible ? tdiff.min_abs_diff : 0.0);
  lua_setfield(L, -2, "tm_min_abs_diff");

  return 1;
}

// Retry the default four-component integrated moments with a baseline M1 scale.
static int
scaled_integrated_m1_diff_lw(lua_State *L)
{
  const char *bf = luaL_checkstring(L, 1), *cf = luaL_checkstring(L, 2);
  double tol = luaL_checknumber(L, 3);
  struct gkyl_dynvec_etype_ncomp be = gkyl_dynvec_read_ncomp(bf);
  struct gkyl_dynvec_etype_ncomp ce = gkyl_dynvec_read_ncomp(cf);
  if (be.type != GKYL_DOUBLE || ce.type != GKYL_DOUBLE || be.ncomp != 4 || ce.ncomp != 4) {
    lua_pushnil(L);
    return 1;
  }
  gkyl_dynvec b = gkyl_dynvec_new(GKYL_DOUBLE, 4), c = gkyl_dynvec_new(GKYL_DOUBLE, 4);
  bool valid = gkyl_dynvec_read(b, bf) && gkyl_dynvec_read(c, cf) &&
               gkyl_dynvec_size(b) == gkyl_dynvec_size(c);
  if (!valid) {
    gkyl_dynvec_release(b);
    gkyl_dynvec_release(c);
    lua_pushnil(L);
    return 1;
  }
  size_t count = gkyl_dynvec_size(b);
  if (count == 0) {
    gkyl_dynvec_release(b);
    gkyl_dynvec_release(c);
    lua_pushboolean(L, 1);
    lua_pushnumber(L, 0);
    return 2;
  }
  // Convert records to arrays so each time sample has the same layout.
  struct gkyl_array *ba = gkyl_array_new(GKYL_DOUBLE, 4, count);
  struct gkyl_array *ca = gkyl_array_new(GKYL_DOUBLE, 4, count);
  struct gkyl_array *bt = gkyl_array_new(GKYL_DOUBLE, 1, count);
  struct gkyl_array *ct = gkyl_array_new(GKYL_DOUBLE, 1, count);
  gkyl_dynvec_to_array(b, bt, ba);
  gkyl_dynvec_to_array(c, ct, ca);
  bool passed = true;
  double worst = 0;
  for (size_t i = 0; i < count; ++i) {
    const double *bv = gkyl_array_cfetch(ba, i), *cv = gkyl_array_cfetch(ca, i);
    for (int k = 0; k < 4; ++k) {
      if (!finite_double(bv[k]) || !finite_double(cv[k])) {
        passed = false;
        worst = infinite_double();
        break;
      }
    }
    if (!passed) {
      break;
    }
    // M0 and M2 set the physical scale for the near-zero M1 component.
    double scale = bv[0] > 0 && bv[2] > 0 ? sqrt(bv[0]) * sqrt(bv[2]) : 0;
    if (!(scale > 0) || !finite_double(scale)) {
      valid = false;
      break;
    }
    for (int k = 0; k < 4; ++k) {
      double diff = fabs(cv[k] - bv[k]);
      if (!finite_double(diff) || !finite_double(bv[k]) || !finite_double(cv[k])) {
        passed = false;
        worst = infinite_double();
        break;
      }
      double limit = tol * (k == 1 ? fmax(fabs(bv[k]), scale) : fmax(1, fabs(bv[k])));
      if (diff > limit) {
        passed = false;
      }
      if (k == 1) {
        worst = fmax(worst, diff / fmax(fabs(bv[k]), scale));
      }
    }
  }
  gkyl_array_release(ba);
  gkyl_array_release(ca);
  gkyl_array_release(bt);
  gkyl_array_release(ct);
  gkyl_dynvec_release(b);
  gkyl_dynvec_release(c);
  if (!valid) {
    lua_pushnil(L);
    return 1;
  }
  lua_pushboolean(L, passed);
  lua_pushnumber(L, worst);
  return 2;
}

// Compare Vlasov integrated M0, vector M1, M2 records with a baseline M1 scale.
static int
scaled_vlasov_integrated_diff_lw(lua_State *L)
{
  const char *bf = luaL_checkstring(L, 1), *cf = luaL_checkstring(L, 2);
  double tol = luaL_checknumber(L, 3), abs_tol = luaL_checknumber(L, 4);
  struct gkyl_dynvec_etype_ncomp be = gkyl_dynvec_read_ncomp(bf);
  struct gkyl_dynvec_etype_ncomp ce = gkyl_dynvec_read_ncomp(cf);
  if (be.type != GKYL_DOUBLE || ce.type != GKYL_DOUBLE || be.ncomp != ce.ncomp || be.ncomp < 3 ||
      be.ncomp > 5) {
    lua_pushnil(L);
    return 1;
  }
  int nc = be.ncomp;
  gkyl_dynvec b = gkyl_dynvec_new(GKYL_DOUBLE, nc), c = gkyl_dynvec_new(GKYL_DOUBLE, nc);
  bool valid = gkyl_dynvec_read(b, bf) && gkyl_dynvec_read(c, cf) &&
               gkyl_dynvec_size(b) == gkyl_dynvec_size(c);
  if (!valid) {
    gkyl_dynvec_release(b);
    gkyl_dynvec_release(c);
    lua_pushnil(L);
    return 1;
  }
  size_t count = gkyl_dynvec_size(b);
  // Empty vectors use one allocated slot, but have no records to compare.
  struct gkyl_array *ba = gkyl_array_new(GKYL_DOUBLE, nc, count ? count : 1);
  struct gkyl_array *ca = gkyl_array_new(GKYL_DOUBLE, nc, count ? count : 1);
  struct gkyl_array *bt = gkyl_array_new(GKYL_DOUBLE, 1, count ? count : 1);
  struct gkyl_array *ct = gkyl_array_new(GKYL_DOUBLE, 1, count ? count : 1);
  if (count) {
    gkyl_dynvec_to_array(b, bt, ba);
    gkyl_dynvec_to_array(c, ct, ca);
  }
  bool passed = true;
  double worst = 0;
  for (size_t i = 0; i < count; ++i) {
    const double *bv = gkyl_array_cfetch(ba, i), *cv = gkyl_array_cfetch(ca, i);
    if (!finite_double(bv[0]) || !finite_double(bv[nc - 1]) || !(bv[0] > 0) || !(bv[nc - 1] > 0)) {
      valid = false;
      break;
    }
    // Baseline density and energy normalize every momentum component.
    double scale = sqrt(bv[0]) * sqrt(bv[nc - 1]);
    if (!(scale > 0) || !finite_double(scale)) {
      valid = false;
      break;
    }
    for (int k = 0; k < nc; ++k) {
      double diff = fabs(cv[k] - bv[k]);
      if (!finite_double(bv[k]) || !finite_double(cv[k]) || !finite_double(diff)) {
        passed = false;
        worst = infinite_double();
        break;
      }
      if (k > 0 && k < nc - 1) {
        double ratio = diff / fmax(fabs(bv[k]), scale);
        worst = fmax(worst, ratio);
        if (ratio > tol) {
          passed = false;
        }
      } else if (diff > abs_tol && diff > tol * (fabs(bv[k]) + fabs(cv[k]))) {
        passed = false;
      }
    }
  }
  gkyl_array_release(ba);
  gkyl_array_release(ca);
  gkyl_array_release(bt);
  gkyl_array_release(ct);
  gkyl_dynvec_release(b);
  gkyl_dynvec_release(c);
  if (!valid) {
    lua_pushnil(L);
    return 1;
  }
  lua_pushboolean(L, passed);
  lua_pushnumber(L, worst);
  return 2;
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
    if (bt1) {
      gkyl_block_topo_release(bt1);
    }
    if (bt2) {
      gkyl_block_topo_release(bt2);
    }
    lua_pushboolean(L, 0);
    return 1;
  }

  bool equal = (bt1->ndim == bt2->ndim) && (bt1->num_blocks == bt2->num_blocks);

  if (equal) {
    int ndim = bt1->ndim;
    int nblocks = bt1->num_blocks;
    for (int b = 0; b < nblocks && equal; ++b) {
      for (int d = 0; d < ndim && equal; ++d) {
        for (int e = 0; e < 2 && equal; ++e) {
          struct gkyl_target_edge *c1 = &bt1->conn[b].connections[d][e];
          struct gkyl_target_edge *c2 = &bt2->conn[b].connections[d][e];
          if (c1->bid != c2->bid || c1->dir != c2->dir || c1->edge != c2->edge) {
            equal = false;
          }
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
  {"gkylFileType", gkyl_file_type_lw},
  {"arrayNewFromFile", array_new_from_file_lw},
  {"rectGridCmp", rect_grid_cmp_lw},
  {"createGridRanges", create_grid_ranges_lw},
  {"arrayDiff", array_diff_lw},
  {"dynvecDiff", dynvec_diff_lw},
  {"scaledMomentDiff", scaled_moment_diff_lw},
  {"scaledIntegratedM1Diff", scaled_integrated_m1_diff_lw},
  {"scaledVlasovMomentDiff", scaled_vlasov_moment_diff_lw},
  {"scaledVlasovIntegratedDiff", scaled_vlasov_integrated_diff_lw},
  {"arrayFinite", array_finite_lw},
  {"blockTopoCmp", block_topo_cmp_lw},
  {0, 0}
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
