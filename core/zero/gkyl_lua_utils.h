#pragma once

#ifdef GKYL_HAVE_LUA

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <stdio.h>
#include <stdbool.h>

// Gets a table with given name from a table, pushing it on the
// stack. Table is popped when the scope ends
#define with_lua_tbl_tbl(L, key)                                      \
    for (bool _break = (lua_getfield(L, -1, key), (lua_isnil(L,-1) || !lua_istable(L, -1) ? (lua_pop(L, 1), false) : true)); \
         _break;                                                        \
         _break = false, lua_pop(L, 1))

// Pushes the value associated with key on stack, popping it when the scope exits
#define with_lua_tbl_key(L, key)                                        \
    for (bool _break = (lua_getfield(L, -1, key), (lua_isnil(L,-1) ? (lua_pop(L, 1), false) : true)); \
         _break;                                                        \
         _break = false, lua_pop(L, 1))

// This macro pushes the named global on the stack and restores the
// stack when the scope is complete. If the global is a function it is
// not popped on scope exit. This is not a problem if the function is
// used inside the scope. If you do not use it, then you must call an
// explict pop yourself.
#define with_lua_global(L, name)                                        \
    for (bool _break = (lua_getglobal(L, name), (lua_isnil(L,-1) ? (lua_pop(L, 1), false) : true)), _isfun = lua_isfunction(L,-1); \
         _break;                                                        \
         _break = false, _isfun ? 0 : lua_pop(L, 1))

// Check and fetch user-data based on metatable name
#define GKYL_CHECK_UDATA(L, mnm) luaL_checkudata(L, 1, mnm)

// For debugging top of stack
#define gkyl_lua_trace_stack_top(L, fnm) do { \
      fprintf(stdout, "Inside function %s\n", fnm);                              \
      fprintf(stdout, "--> Top of stack is %s\n", lua_typename(L, lua_type(L, -1))); \
    } while (0);

/**
 * Get length of object on top of stack. (Table size, string length
 * and memory allocated)
 *
 * @return Length of object on top of stack.
 */
static inline size_t glua_objlen(lua_State *L) { return lua_objlen(L, -1); }

/**
 * Check if table has specifed key in it. Table must be on top of the
 * stack
 *
 * @param L Lua state
 * @param key Key to check
 * @return true if key is present, false otherwise
 */
bool glua_tbl_has_key(lua_State *L, const char *key);

// In the following table access methods, the get_XYZ methods take a
// string as a key while the iget_XYZ takes a integer index as key.

/**
 * Return number from table, keyed by @a key. Table must be on top of
 * the stack. If number does not exist, @a def is returned.
 *
 * @param L Lua state
 * @param key Key
 * @param def Default value if key is not present in table
 * @return number corresponding to key, or def
 */
double glua_tbl_get_number(lua_State *L, const char *key, double def);
double glua_tbl_iget_number(lua_State *L, long key, double def);

/**
 * Return integer from table, keyed by @a key. Table must be on top of
 * the stack. If number does not exist, @a def is returned.
 *
 * @param L Lua state
 * @param key Key
 * @param def Default value if key is not present in table
 * @return integer corresponding to key, or def
 */
long glua_tbl_get_integer(lua_State *L, const char *key, long def);
long glua_tbl_iget_integer(lua_State *L, long key, long def);

/**
 * Return boolean from table, keyed by @a key. Table must be on top of
 * the stack. If bool does not exist, @a def is returned. (In Lua, all
 * values except false and nil are true).
 *
 * @param L Lua state
 * @param key Key
 * @param def Default value if key is not present in table
 * @return integer corresponding to key, or def
 */
int glua_tbl_get_bool(lua_State *L, const char *key, int def);
int glua_tbl_iget_bool(lua_State *L, long key, int def);

/**
 * Return string from table, keyed by @a key. Table must be on top of
 * the stack. If string does not exist, @a def is returned.
 *
 * @param L Lua state
 * @param key Key
 * @param def Default value if key is not present in table
 * @return string corresponding to key, or def
 */
const char *glua_tbl_get_string(lua_State *L, const char *key, const char *def);
const char *glua_tbl_iget_string(lua_State *L, long key, const char *def);

/**
 * Fetches table named @a key from table on top of stack and pushes
 * that on stack. Returns false if no such table was found. You must
 * pop the table youself once you are done with it.
 *
 * @param L Lua state.
 * @param key Name of table to fetch
 * @return true if table exists, false otherwise
 */
bool glua_tbl_get_tbl(lua_State *L, const char *key);
bool glua_tbl_iget_tbl(lua_State *L, long key);

/**
 * Fetches function named @a key from table on top of stack and pushes
 * that on stack. Returns false if no such function was found.
 *
 * @param L Lua state.
 * @param key Name of function to fetch
 * @return true if function exists, false otherwise 
 */
bool glua_tbl_get_func(lua_State *L, const char *key);
bool glua_tbl_iget_func(lua_State *L, long key);

/**
 * Run Lua code stored in @a str buffer. The size of the buffer is
 * sz. This buffer can either be Lua byte-code or plain-text Lua code.
 * On error the error message is written to the @a err file
 * stream. Pass NULL to not print any error.
 *
 * @param L Lua state
 * @param str Buffer with Lua code (text or bytecode)
 * @param sz Size of buffer
 * @param err Error stream to write errors (if any)
 * @return 0 on success.
 */
int glua_run_lua(lua_State *L, const char *str, long sz, FILE *err);

#endif
