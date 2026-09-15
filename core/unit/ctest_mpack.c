#include <acutest.h>
#include <mpack.h>
#include <gkyl_util.h>

void
test_mpack_map_1_ho(void)
{
  char *data;
  mpack_writer_t writer;
  size_t data_length = 0;
  mpack_writer_init_growable(&writer, &data, &data_length);

  // add some data to mpack
  mpack_build_map(&writer);
  
  mpack_write_cstr(&writer, "time");
  mpack_write_double(&writer, 1.5);

  mpack_write_cstr(&writer, "frame");
  mpack_write_i64(&writer, 100);

  mpack_write_cstr(&writer, "basis");
  mpack_write_cstr(&writer, "ms");

  mpack_write_cstr(&writer, "values");
  mpack_build_array(&writer);
  mpack_write_double(&writer, 0.1);
  mpack_write_double(&writer, 0.2);
  mpack_write_double(&writer, 0.3);
  mpack_complete_array(&writer);

  mpack_complete_map(&writer);

  // finish writing
  int status = mpack_writer_destroy(&writer);
  TEST_CHECK( mpack_ok == status );

  // read it back
  mpack_tree_t tree;
  mpack_tree_init_data(&tree, data, data_length);
  mpack_tree_parse(&tree);

  mpack_node_t root = mpack_tree_root(&tree);
  TEST_CHECK(mpack_node_type(root) == mpack_type_map);

  TEST_CHECK( 4 == mpack_node_map_count(root) );

  TEST_CHECK( mpack_node_map_contains_cstr(root, "time") );
  TEST_CHECK( mpack_node_map_contains_cstr(root, "bogus") == false );

  mpack_node_t tm_node = mpack_node_map_cstr(root, "time");
  TEST_CHECK( mpack_node_double(tm_node) == 1.5 );

  mpack_node_t fr_node = mpack_node_map_cstr(root, "frame");
  TEST_CHECK( mpack_node_i64(fr_node) == 100 );

  mpack_node_t bs_node = mpack_node_map_cstr(root, "basis");
  char *basis_str = mpack_node_cstr_alloc(bs_node, mpack_node_strlen(bs_node)+1);
  TEST_CHECK( strcmp("ms", basis_str) == 0 );
  MPACK_FREE(basis_str);

  mpack_node_t vl_node = mpack_node_map_cstr(root, "values");
  TEST_CHECK( mpack_node_type(vl_node) == mpack_type_array );
  TEST_CHECK( mpack_node_array_length(vl_node) == 3 );

  double v = 0.1;
  for (int i=0; i<3; ++i) {
    mpack_node_t dn = mpack_node_array_at(vl_node, i);
    TEST_CHECK( gkyl_compare_double(v, mpack_node_double(dn), 1e-15) );
    v += 0.1;
  }
    
  status = mpack_tree_destroy(&tree);
  TEST_CHECK( mpack_ok == status );

  free(data);
}

void
test_msgpack_1_ho(void)
{
  struct gkyl_msgpack_map_elem elist[] = {
    GKYL_MSGPACK_MAP_ELEM("bool", true), // 0
    GKYL_MSGPACK_MAP_ELEM("int", -10), // 1
    GKYL_MSGPACK_MAP_ELEM("unsigned", 10u), // 2
    GKYL_MSGPACK_MAP_ELEM("float", 1.0f), // 3
    GKYL_MSGPACK_MAP_ELEM("double", 1.0), // 4
    GKYL_MSGPACK_MAP_ELEM("string", "Hello, World!") // 5
  };

  TEST_CHECK( elist[0].elem_type == GKYL_MP_INT ); // bools are int?
  TEST_CHECK( elist[0].bval == true );  
  
  TEST_CHECK( elist[1].elem_type == GKYL_MP_INT );
  TEST_CHECK( elist[1].ival == -10 );    

  TEST_CHECK( elist[2].elem_type == GKYL_MP_UNSIGNED_INT );
  TEST_CHECK( elist[2].uval == 10u );  

  TEST_CHECK( elist[3].elem_type == GKYL_MP_FLOAT );
  TEST_CHECK( elist[3].fval == 1.0f );  

  TEST_CHECK( elist[4].elem_type == GKYL_MP_DOUBLE );
  TEST_CHECK( elist[4].dval == 1.0 );  

  TEST_CHECK( elist[5].elem_type == GKYL_MP_STRING );
  TEST_CHECK( strcmp(elist[5].cval, "Hello, World!") == 0 );
}

TEST_LIST = {
  { "mpack_map_1_ho", test_mpack_map_1_ho },
  { "msgpack_1_ho", test_msgpack_1_ho },  
  { 0, 0 },
};
