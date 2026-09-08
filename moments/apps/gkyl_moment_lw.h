#pragma once

#ifdef GKYL_HAVE_LUA

#include <lua.h>

/**
 * Load the Moments Lua wrapper module into the Lua interpreter.
 *
 * @param L Lua state
 */
void gkyl_moment_lw_openlibs(lua_State *L);

/**
* Add moment reconstruction type flags for fluid solvers into Lua interpreter.
*
* @param L Lua state to use.
 */
void gkyl_register_moment_scheme_types(lua_State *L);

/**
* Add wave limiter type flags for fluid solvers into Lua interpreter.
*
* @param L Lua state to use.
 */
void gkyl_register_wave_limiter_types(lua_State *L);

/**
* Add Riemann problem type flags for Euler equations into Lua interpreter.
*
* @param L Lua state to use.
 */
void gkyl_register_euler_rp_types(lua_State *L);

/**
* Add Riemann problem type flags for MHD equations into Lua interpreter.
*
* @param L Lua state to use.
 */
void gkyl_register_mhd_rp_types(lua_State *L);

/**
* Add divergence correction type flags for MHD equations into Lua interpreter.
*
* @param L Lua state to use.
 */
void gkyl_register_mhd_divb_types(lua_State *L);

/**
* Add Braginskii type flags for moment equations into Lua interpreter.
*
* @param L Lua state to use.
 */
void gkyl_register_braginskii_types(lua_State *L);

/**
* Add spacetime gauge type flags for moment equations into Lua interpreter.
*
* @param L Lua state to use.
 */
void gkyl_register_spacetime_gauge_types(lua_State *L);

/**
* Add spacetime slicing type flags for Einstein equations into Lua interpreter.
*
* @param L Lua state to use.
 */
void gkyl_register_spacetime_slicing_types(lua_State *L);

/**
* Add spacetime evolution type flags for Einstein equations into Lua interpreter.
*
* @param L Lua state to use.
 */
void gkyl_register_spacetime_evolution_types(lua_State *L);

#endif
