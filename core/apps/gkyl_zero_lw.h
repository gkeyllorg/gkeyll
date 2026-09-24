#pragma once

#ifdef GKYL_HAVE_LUA

#include <lua.h>

// struct to construct argc and argv from GKYL_COMMAND global table
struct gkyl_tool_args {
  int argc; // number of arguments
  char **argv; // arguments
};

/**
 * Construct Tool or script CLI args. Must be freed by calling the
 * "release" method. This method assumes that the commands are in a
 * global Lua table called GKYL_COMMANDS.
 *
 * @param L Lua state
 * @return New argument struct.
 */
struct gkyl_tool_args *gkyl_tool_args_new(lua_State *L);

/**
 * Free tools list
 *
 * @param args Argument struct to free
 */
void gkyl_tool_args_release(struct gkyl_tool_args* args);

/**
 * Load the zero-level Lua interface into the Lua interpreter.
 *
 * @param L Lua state
 */
void gkyl_zero_lw_openlibs(lua_State *L);

#endif
