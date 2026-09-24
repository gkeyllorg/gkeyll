#include <acutest.h>
#include <gkyl_block_geom.h>

static void
test_L_domain()
{
  // 2D with 3 blocks
  struct gkyl_block_geom *bgeom = gkyl_block_geom_new(2, 3);

  TEST_CHECK( 0 == gkyl_block_geom_check_consistency(bgeom) );
  TEST_CHECK( 3 == gkyl_block_geom_num_blocks(bgeom) );

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
  gkyl_block_geom_set_block(bgeom, 0, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 10, 10 },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }  // physical boundary
      },
      .connections[1] = { // y-direction connections
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } // physical boundary
      }
    }
  );
  // topology is inconsistent at this point!
  TEST_CHECK( 0 == gkyl_block_geom_check_consistency(bgeom) );
  
  // block 1
  gkyl_block_geom_set_block(bgeom, 1, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 10, 10 },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE }
      },
      .connections[1] = { // y-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 0, .dir = 1, .edge = GKYL_LOWER_POSITIVE }
      }
    }
  );
  // geomlogy is inconsistent at this point!
  TEST_CHECK( 0 == gkyl_block_geom_check_consistency(bgeom) );

  // block 2
  gkyl_block_geom_set_block(bgeom, 2, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 10, 10 },
      
      .connections[0] = { // x-direction connections
        { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL } // physical boundary
      },
      .connections[1] = { // y-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } // physical boundary
      }
    }
  );

  // should be fully consistent as all blocks properly specified
  TEST_CHECK( 1 == gkyl_block_geom_check_consistency(bgeom) );

  struct gkyl_block_topo *btopo = gkyl_block_geom_topo(bgeom);

  TEST_CHECK( 3 == btopo->num_blocks );
  TEST_CHECK( 1 == gkyl_block_topo_check_consistency(btopo) );

  int nblocks = gkyl_block_geom_num_blocks(bgeom);
  for (int i=0; i<nblocks; ++i) {
    const struct gkyl_block_geom_info *ginfo = gkyl_block_geom_get_block(bgeom, i);
    TEST_CHECK( 10 == ginfo->cells[0] );
    TEST_CHECK( 10 == ginfo->cells[1] );

    TEST_CHECK( 1 == ginfo->cuts[0] );
    TEST_CHECK( 1 == ginfo->cuts[1] );
  }
    

  gkyl_block_geom_release(bgeom);
  gkyl_block_topo_release(btopo);
}

static void
test_mobius_domain()
{
  // 2D with 1 block
  struct gkyl_block_geom *bgeom = gkyl_block_geom_new(2, 1);

  /* Block layout

     periodic
     +------+
     |0     | periodic with twist
     |      |
     +------+
    
  */

  // block 0
  gkyl_block_geom_set_block(bgeom, 0, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 10, 10 },      
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_UPPER_NEGATIVE }, // note twist
        { .bid = 0, .dir = 0, .edge = GKYL_LOWER_NEGATIVE }
      },
      .connections[1] = { // y-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 0, .dir = 1, .edge = GKYL_LOWER_POSITIVE }
      }
    }
  );

  TEST_CHECK( 1 == gkyl_block_geom_check_consistency(bgeom) );
  
  gkyl_block_geom_release(bgeom);
}



TEST_LIST = {
  { "mobius_domain", test_mobius_domain },
  { "L_domain", test_L_domain },
  { NULL, NULL },
};
