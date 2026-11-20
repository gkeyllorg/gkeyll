#pragma once

#ifdef GKYL_HAVE_LUA

#include <lua.h>

/**
 * Load the PKPM Lua wrapper module into the Lua
 * interpreter.
 *
 * @param L Lua state
 */
void gkyl_pkpm_lw_openlibs(lua_State *L);

#endif
