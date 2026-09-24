#include <acutest.h>

#ifdef GKYL_HAVE_LUA

#include <gkyl_lua_utils.h>

static lua_State *
new_lua_State(void)
{
  lua_State *L = luaL_newstate();
  lua_gc(L, LUA_GCSTOP, 0);
  luaL_openlibs(L);
  lua_gc(L, LUA_GCRESTART, -1);
  return L;
}

void
test_0(void)
{
  lua_State *L = new_lua_State();

  const char *lcode1 = "tbl = { x = 1 }";
  glua_run_lua(L, lcode1, strlen(lcode1), stderr);

  with_lua_global(L, "tbl") {

    with_lua_tbl_tbl(L, "nope") {
      TEST_CHECK(  false );
    }

    with_lua_tbl_tbl(L, "x") {
      TEST_CHECK(  false );
    }

    with_lua_tbl_key(L, "nope") {
      TEST_CHECK(  false );
    }

    TEST_CHECK( 1 == glua_tbl_get_integer(L, "x", 0) );
  }
  
  lua_close(L);  
}

void
test_1(void)
{
  lua_State *L = new_lua_State();
  
  const char *lcode1 =
    "kvpairs = { x1 = 25.5, n1 = 2345, x2 = 32.5, n2 = 1234, bc = \"periodic\" }";
  glua_run_lua(L, lcode1, strlen(lcode1), stderr);

  const char *lcode2 =
    "nums = { 0.5, 10.5, 20.5, 30.5 }";
  glua_run_lua(L, lcode2, strlen(lcode2), stderr);

  TEST_CHECK( 0 == lua_gettop(L) );
  with_lua_global(L, "does-not-exist") {
    TEST_CHECK( false );
  }
  TEST_CHECK( 0 == lua_gettop(L) );  

  with_lua_global(L, "kvpairs") {
    TEST_CHECK( 1 == lua_gettop(L) );
  }
  TEST_CHECK( 0 == lua_gettop(L) );

  // push table on top of stack
  with_lua_global(L, "kvpairs") {

    TEST_CHECK( 0 == glua_objlen(L) );

    TEST_CHECK( glua_tbl_has_key(L, "x1") );
    TEST_CHECK( glua_tbl_has_key(L, "x3") == false );
  
    TEST_CHECK( 25.5 == glua_tbl_get_number(L, "x1", 0.0) );
    TEST_CHECK( 32.5 == glua_tbl_get_number(L, "x2", 0.0) );
    TEST_CHECK( 12.5 == glua_tbl_get_number(L, "x3", 12.5) );

    TEST_CHECK( glua_tbl_has_key(L, "n1") );
    TEST_CHECK( glua_tbl_has_key(L, "n3") == false );  
    
    TEST_CHECK( glua_tbl_has_key(L, "bc") );
    TEST_CHECK( strcmp("periodic", glua_tbl_get_string(L, "bc", "N")) == 0 );

    with_lua_tbl_key(L, "bc") {
      const char *val = lua_tostring(L, -1);
      TEST_CHECK( strcmp("periodic", val) == 0 );
    }

    with_lua_tbl_key(L, "does-not-exist") {
      TEST_CHECK( false );
    }

    TEST_CHECK( 2345 == glua_tbl_get_integer(L, "n1", 0) );
    TEST_CHECK( 1234 == glua_tbl_get_integer(L, "n2", 0) );
    TEST_CHECK( 7890 == glua_tbl_get_integer(L, "n3", 7890) );
  }

  TEST_CHECK( 0 == lua_gettop(L) );

  lua_close(L);
}

void
test_2(void)
{
  lua_State *L = new_lua_State();
  
  const char *lcode1 =
    "kvpairs = { tEnd = 101.1, nums = { 11, 12, 13, v = 2222 }, names = { \"Vlasov\", \"Maxwell\"}  }";
  glua_run_lua(L, lcode1, strlen(lcode1), stderr);

  with_lua_global(L, "kvpairs") {

    with_lua_tbl_tbl(L, "nope") {
      TEST_CHECK( false );
    }

    with_lua_tbl_tbl(L, "nums") {
      TEST_CHECK(glua_tbl_has_key(L, "v"));
      TEST_CHECK(2222 == glua_tbl_get_integer(L, "v", 0));

      TEST_CHECK(3 == glua_objlen(L));

      for (int i = 1; i <= glua_objlen(L); ++i)
        TEST_CHECK(10 + i == glua_tbl_iget_integer(L, i, 0));

      for (int i = 1; i <= glua_objlen(L); ++i)
        TEST_CHECK(10 + i == glua_tbl_iget_number(L, i, 0));
    }

    TEST_CHECK(glua_tbl_has_key(L, "names"));

    with_lua_tbl_tbl(L, "names") {
      TEST_CHECK(2 == glua_objlen(L));
      TEST_CHECK(strcmp("Vlasov", glua_tbl_iget_string(L, 1, "Einstein")) == 0);
      TEST_CHECK(strcmp("Maxwell", glua_tbl_iget_string(L, 2, "Einstein")) == 0);
    }

    TEST_CHECK(glua_tbl_has_key(L, "tEnd"));
    TEST_CHECK(101.1 == glua_tbl_get_number(L, "tEnd", 0.0));
  }

  lua_close(L);
}

void
test_3(void)
{
  lua_State *L = new_lua_State();

  const char *lcode1 =
    "function mysq(x) return x*x end; function twov(x) return x, 2*x end";
  glua_run_lua(L, lcode1, strlen(lcode1), stderr);

  TEST_CHECK( 0 == lua_gettop(L) );
  
  with_lua_global(L, "mysq") {
    lua_pushnumber(L, 2.5);
    if (lua_pcall(L, 1, 1, 0)) {
      // signal error condition
    }
    else {
      double res = lua_tonumber(L, -1);
      TEST_CHECK( 2.5*2.5 == res );
      lua_pop(L, 1);
    }
  }

  with_lua_global(L, "twov") {
    lua_pushnumber(L, 2.5);
    if (lua_pcall(L, 1, 2, 0)) {
      // signal error condition
    }
    else {
      // returned values are accessed in reverse order
      double r2 = lua_tonumber(L, -1);
      TEST_CHECK( 2*2.5 == r2 );
      lua_pop(L, 1);

      double r1 = lua_tonumber(L, -1);
      TEST_CHECK( 2.5 == r1 );
      lua_pop(L, 1);
    }
  }  

  TEST_CHECK( 0 == lua_gettop(L) );

  const char *lcode2 = "kvpairs = { sq = function(x) return x*x end, x = 10.5 } ";
  glua_run_lua(L, lcode2, strlen(lcode2), stderr);

  with_lua_global(L, "kvpairs") {

    if (glua_tbl_get_func(L, "no_sq")) {
      TEST_CHECK(false);
    }

    if (glua_tbl_get_func(L, "sq")) {
      lua_pushnumber(L, 3.5);
      if (lua_pcall(L, 1, 1, 0)) {
        // signal error condition
      }
      else {
        double res = lua_tonumber(L, -1);
        TEST_CHECK( 3.5*3.5 == res );
        lua_pop(L, 1);
      }
    }
  }

  TEST_CHECK( 0 == lua_gettop(L) );  

  lua_close(L);
}

void
test_4(void)
{
  /* lua_State *L = new_lua_State(); */
  
  /* const char *lcode1 = */
  /*   "kvpairs = { tEnd = 101.1, nums = { 11, 12, 13, v = 2222 }, names = { \"Vlasov\", \"Maxwell\"}  }"; */
  /* glua_run_lua(L, lcode1, strlen(lcode1), stderr); */

  /* TEST_CHECK( 0 == lua_gettop(L) ); */

  /* with_lua_global(L, "kvpairs") { */
  /*   printf("Stack top: %d\n", lua_gettop(L)); */

  /*   lua_pushnil(L); // initial key is nil */
    
  /*   while (lua_next(L, -2) != 0) { */
  /*     // key at -2 and value at -1 */
      
  /*     printf("-> Stack top: %d\n", lua_gettop(L)); */
  /*     printf("%s - %s\n", */
  /*       lua_typename(L, lua_type(L,-2)), */
  /*       lua_typename(L, lua_type(L,-1))); */
      
  /*     if (lua_type(L,-2) == LUA_TSTRING) { */
  /*       const char *key = lua_tolstring(L, -2, 0); */
  /*       printf("--> key is '%s'\n", key); */
  /*     } */
  /*     lua_pop(L, 1); */
  /*   } */
  /* } */

  /* TEST_CHECK( 0 == lua_gettop(L) ); */

  /* lua_close(L); */

}

TEST_LIST = {
  { "test_0", test_0 },
  { "test_1", test_1 },
  { "test_2", test_2 },
  { "test_3", test_3 },
  {NULL, NULL},
};

#else

// nothing to test if not building with MPI
TEST_LIST = {
  {NULL, NULL},
};

#endif
