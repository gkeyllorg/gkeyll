#pragma once

#ifdef GKYL_HAVE_LUA

#include <lua.h>

/**
 * Load the Vlasov-Maxwell Lua wrapper module into the Lua
 * interpreter.
 *
 * @param L Lua state
 */
void gkyl_vlasov_lw_openlibs(lua_State *L);

/**
 * Add boundary condition flags for Poisson solver into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_poisson_bc_types(lua_State *L);

/**
 * Add projection type flags for Vlasov species initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_vlasov_projection_types(lua_State *L);

/**
 * Add model type flags for Vlasov species initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_vlasov_model_types(lua_State *L);

/**
 * Add collision type flags for Vlasov species initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_vlasov_collision_types(lua_State *L);

/**
 * Add source type flags for Vlasov species initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_vlasov_source_types(lua_State *L);


#endif
