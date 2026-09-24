#pragma once

#ifdef GKYL_HAVE_LUA

#include <lua.h>

/**
 * Load the gyrokinetic Lua wrapper module into the Lua
 * interpreter.
 *
 * @param L Lua state
 */
void gkyl_gyrokinetic_lw_openlibs(lua_State *L);


/**
 * Add FEM boundary condition flags for gyrokinetic field initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_gyrokinetic_fem_bc_types(lua_State *L);

/**
 * Add geometry type flags for gyrokinetic app initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_gyrokinetic_geometry_types(lua_State *L);

/**
 * Add position map type flags for gyrokinetic app initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_gyrokinetic_position_map_types(lua_State *L);

/**
 * Add field type flags for gyrokinetic field initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_gyrokinetic_field_types(lua_State *L);

/**
 * Add radiation type flags for gyrokinetic species initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_gyrokinetic_radiation_types(lua_State *L);

/**
 * Add Te model type flags for gyrokinetic radiation initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_gyrokinetic_radiation_Te_types(lua_State *L);

/**
 * Add reaction type flags for gyrokinetic species initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_gyrokinetic_reaction_types(lua_State *L);

/**
 * Add ion type flags for gyrokinetic species initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_gyrokinetic_ion_types(lua_State *L);

/**
 * Add self-reaction type flags for gyrokinetic species initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_gyrokinetic_self_reaction_types(lua_State *L);

#endif
