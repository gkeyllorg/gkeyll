#include <acutest.h>
#include <gkyl_block_topo.h>
#include <gkyl_array_rio.h>

static struct gkyl_block_topo *
create_L_domain(void)
{
  // 2D with 3 blocks
  struct gkyl_block_topo *btopo = gkyl_block_topo_new(2, 3);

  TEST_CHECK( 3 == btopo->num_blocks );
  // topology is inconsistent at this point!
  TEST_CHECK( 0 == gkyl_block_topo_check_consistency(btopo) );

  /* Block layout

     +------+
     |0     |
     |      |
     +------+-----+
     |1     |2    |
     |      |     |
     +------+-----+
    
  */

  // block 0
  btopo->conn[0] = (struct gkyl_block_connections) {
    .connections[0] = { // x-direction connections
      { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, // physical boundary
      { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }  // physical boundary
    },
    .connections[1] = { // y-direction connections
      { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE },
      { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } // physical boundary
    }
  };
  // topology is inconsistent at this point!
  TEST_CHECK( 0 == gkyl_block_topo_check_consistency(btopo) );
  
  // block 1
  btopo->conn[1] = (struct gkyl_block_connections) {
    .connections[0] = { // x-direction connections
      { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, // physical boundary
      { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE }
    },
    .connections[1] = { // y-direction connections
      { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
      { .bid = 0, .dir = 1, .edge = GKYL_LOWER_POSITIVE }
    }
  };
  // topology is inconsistent at this point!
  TEST_CHECK( 0 == gkyl_block_topo_check_consistency(btopo) );
  
  // block 2
  btopo->conn[2] = (struct gkyl_block_connections) {
    .connections[0] = { // x-direction connections
      { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
      { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL } // physical boundary
    },
    .connections[1] = { // y-direction connections
      { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
      { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } // physical boundary
    }
  };

  return btopo;
}

static void
test_L_domain()
{
  struct gkyl_block_topo *btopo = create_L_domain();

  // should be fully consistent as all blocks properly specified
  TEST_CHECK( 1 == gkyl_block_topo_check_consistency(btopo) );

  gkyl_block_topo_release(btopo);
}

static void
test_mobius_domain()
{
  // 2D with 1 block
  struct gkyl_block_topo *btopo = gkyl_block_topo_new(2, 1);

  TEST_CHECK( 1 == btopo->num_blocks );
  // topology is inconsistent at this point!
  TEST_CHECK( 0 == gkyl_block_topo_check_consistency(btopo) );

  /* Block layout

     periodic
     +------+
     |0     | periodic with twist
     |      |
     +------+
    
  */

  // block 0
  btopo->conn[0] = (struct gkyl_block_connections) {
    .connections[0] = { // x-direction connections
      { .bid = 0, .dir = 0, .edge = GKYL_UPPER_NEGATIVE }, // note twist
      { .bid = 0, .dir = 0, .edge = GKYL_LOWER_NEGATIVE }
    },
    .connections[1] = { // y-direction connections
      { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE },
      { .bid = 0, .dir = 1, .edge = GKYL_LOWER_POSITIVE }
    }
  };

  TEST_CHECK( 1 == gkyl_block_topo_check_consistency(btopo) );

  gkyl_block_topo_release(btopo);
}

static void
test_topo_io()
{
  struct gkyl_block_topo *btopo = create_L_domain();
  int status_out = gkyl_block_topo_write(btopo, "ctest_block_topo_L_domain.gkyl");
  TEST_CHECK( GKYL_ARRAY_RIO_SUCCESS == status_out );

  int status_inp;
  struct gkyl_block_topo *btopo_inp = gkyl_block_topo_read("ctest_block_topo_L_domain.gkyl",
    &status_inp);
  TEST_CHECK( GKYL_ARRAY_RIO_SUCCESS == status_out );

  if (status_out == GKYL_ARRAY_RIO_SUCCESS) {
    TEST_CHECK( 2 == btopo_inp->ndim );
    TEST_CHECK( 3 == btopo_inp->num_blocks );
    TEST_CHECK( 1 == gkyl_block_topo_check_consistency(btopo_inp) );
    
    gkyl_block_topo_release(btopo_inp);
  }
  
  gkyl_block_topo_release(btopo);
}

TEST_LIST = {
  { "mobius_domain", test_mobius_domain },
  { "L_domain", test_L_domain },
  { "topo_io", test_topo_io },  
  { NULL, NULL },
};
