#include <acutest.h>

#include <gkyl_rrobin_decomp.h>

static void
test_rrobin_decomp_1_ho(void)
{
  const struct gkyl_rrobin_decomp *rr = gkyl_rrobin_decomp_new(1, 3,
    (int[]) { 1, 1, 1 });

  TEST_CHECK( 1 == rr->total_ranks );
  TEST_CHECK( 3 == rr->nblocks );

  int b0_ranks[1];
  gkyl_rrobin_decomp_getranks(rr, 0, b0_ranks);
  TEST_CHECK( 0 == b0_ranks[0] );

  int b1_ranks[1];
  gkyl_rrobin_decomp_getranks(rr, 1, b1_ranks);
  TEST_CHECK( 0 == b1_ranks[0] );

  int b2_ranks[1];
  gkyl_rrobin_decomp_getranks(rr, 2, b2_ranks);
  TEST_CHECK( 0 == b2_ranks[0] );
  
  gkyl_rrobin_decomp_release(rr);
}

static void
test_rrobin_decomp_2_ho(void)
{
  const struct gkyl_rrobin_decomp *rr = gkyl_rrobin_decomp_new(4, 3,
    (int[]) { 4, 1, 1 });

  TEST_CHECK( 4 == rr->total_ranks );
  TEST_CHECK( 3 == rr->nblocks );

  int b0_ranks[4];
  gkyl_rrobin_decomp_getranks(rr, 0, b0_ranks);
  TEST_CHECK( 0 == b0_ranks[0] );
  TEST_CHECK( 1 == b0_ranks[1] );
  TEST_CHECK( 2 == b0_ranks[2] );
  TEST_CHECK( 3 == b0_ranks[3] );

  int b1_ranks[1];
  gkyl_rrobin_decomp_getranks(rr, 1, b1_ranks);
  TEST_CHECK( 0 == b1_ranks[0] );

  int b2_ranks[1];
  gkyl_rrobin_decomp_getranks(rr, 2, b2_ranks);
  TEST_CHECK( 1 == b2_ranks[0] );
  
  gkyl_rrobin_decomp_release(rr);
}

static void
test_rrobin_decomp_3_ho(void)
{
  const struct gkyl_rrobin_decomp *rr = gkyl_rrobin_decomp_new(6, 3,
    (int[]) { 4, 1, 1 });

  TEST_CHECK( 6 == rr->total_ranks );
  TEST_CHECK( 3 == rr->nblocks );

  int b0_ranks[4];
  gkyl_rrobin_decomp_getranks(rr, 0, b0_ranks);
  TEST_CHECK( 0 == b0_ranks[0] );
  TEST_CHECK( 1 == b0_ranks[1] );
  TEST_CHECK( 2 == b0_ranks[2] );
  TEST_CHECK( 3 == b0_ranks[3] );

  int b1_ranks[1];
  gkyl_rrobin_decomp_getranks(rr, 1, b1_ranks);
  TEST_CHECK( 4 == b1_ranks[0] );

  int b2_ranks[1];
  gkyl_rrobin_decomp_getranks(rr, 2, b2_ranks);
  TEST_CHECK( 5 == b2_ranks[0] );
  
  gkyl_rrobin_decomp_release(rr);
}

static void
test_rrobin_decomp_4_ho(void)
{
  const struct gkyl_rrobin_decomp *rr = gkyl_rrobin_decomp_new(6, 3,
    (int[]) { 4, 2, 3 });

  TEST_CHECK( 6 == rr->total_ranks );
  TEST_CHECK( 3 == rr->nblocks );

  int b0_ranks[4];
  gkyl_rrobin_decomp_getranks(rr, 0, b0_ranks);
  TEST_CHECK( 0 == b0_ranks[0] );
  TEST_CHECK( 1 == b0_ranks[1] );
  TEST_CHECK( 2 == b0_ranks[2] );
  TEST_CHECK( 3 == b0_ranks[3] );

  int b1_ranks[2];
  gkyl_rrobin_decomp_getranks(rr, 1, b1_ranks);
  TEST_CHECK( 4 == b1_ranks[0] );
  TEST_CHECK( 5 == b1_ranks[1] );

  int b2_ranks[3];
  gkyl_rrobin_decomp_getranks(rr, 2, b2_ranks);
  TEST_CHECK( 0 == b2_ranks[0] );
  TEST_CHECK( 1 == b2_ranks[1] );
  TEST_CHECK( 2 == b2_ranks[2] ); 
  
  gkyl_rrobin_decomp_release(rr);
}

static void
test_rrobin_decomp_5_ho(void)
{
  const struct gkyl_rrobin_decomp *rr = gkyl_rrobin_decomp_new(9, 3,
    (int[]) { 4, 2, 3 });

  TEST_CHECK( 9 == rr->total_ranks );
  TEST_CHECK( 3 == rr->nblocks );

  int b0_ranks[4];
  gkyl_rrobin_decomp_getranks(rr, 0, b0_ranks);
  TEST_CHECK( 0 == b0_ranks[0] );
  TEST_CHECK( 1 == b0_ranks[1] );
  TEST_CHECK( 2 == b0_ranks[2] );
  TEST_CHECK( 3 == b0_ranks[3] );

  int b1_ranks[2];
  gkyl_rrobin_decomp_getranks(rr, 1, b1_ranks);
  TEST_CHECK( 4 == b1_ranks[0] );
  TEST_CHECK( 5 == b1_ranks[1] );

  int b2_ranks[3];
  gkyl_rrobin_decomp_getranks(rr, 2, b2_ranks);
  TEST_CHECK( 6 == b2_ranks[0] );
  TEST_CHECK( 7 == b2_ranks[1] );
  TEST_CHECK( 8 == b2_ranks[2] ); 
  
  gkyl_rrobin_decomp_release(rr);
}

static void
test_rrobin_decomp_6_ho(void)
{
  const struct gkyl_rrobin_decomp *rr = gkyl_rrobin_decomp_new(6, 3,
    (int[]) { 2, 4, 3 });

  TEST_CHECK( 6 == rr->total_ranks );
  TEST_CHECK( 3 == rr->nblocks );

  int b0_ranks[2];
  gkyl_rrobin_decomp_getranks(rr, 0, b0_ranks);
  TEST_CHECK( 0 == b0_ranks[0] );
  TEST_CHECK( 1 == b0_ranks[1] );

  int b1_ranks[4];
  gkyl_rrobin_decomp_getranks(rr, 1, b1_ranks);
  TEST_CHECK( 2 == b1_ranks[0] );
  TEST_CHECK( 3 == b1_ranks[1] );
  TEST_CHECK( 4 == b1_ranks[2] );
  TEST_CHECK( 5 == b1_ranks[3] );

  int b2_ranks[3];
  gkyl_rrobin_decomp_getranks(rr, 2, b2_ranks);
  TEST_CHECK( 0 == b2_ranks[0] );
  TEST_CHECK( 1 == b2_ranks[1] );
  TEST_CHECK( 2 == b2_ranks[2] ); 
  
  gkyl_rrobin_decomp_release(rr);
}

static void
test_rrobin_decomp_7_ho(void)
{
  // this is a rather strange decomposition: only 2 ranks while the
  // largest block needs 4 ranks for full concurrency
  
  const struct gkyl_rrobin_decomp *rr = gkyl_rrobin_decomp_new(2, 3,
    (int[]) { 2, 4, 3 });

  TEST_CHECK( 2 == rr->total_ranks );
  TEST_CHECK( 3 == rr->nblocks );

  int b0_ranks[2];
  gkyl_rrobin_decomp_getranks(rr, 0, b0_ranks);
  TEST_CHECK( 0 == b0_ranks[0] );
  TEST_CHECK( 1 == b0_ranks[1] );

  int b1_ranks[4];
  gkyl_rrobin_decomp_getranks(rr, 1, b1_ranks);
  TEST_CHECK( 0 == b1_ranks[0] );
  TEST_CHECK( 1 == b1_ranks[1] );
  TEST_CHECK( 0 == b1_ranks[2] );
  TEST_CHECK( 1 == b1_ranks[3] );

  int b2_ranks[3];
  gkyl_rrobin_decomp_getranks(rr, 2, b2_ranks);
  TEST_CHECK( 0 == b2_ranks[0] );
  TEST_CHECK( 1 == b2_ranks[1] );
  TEST_CHECK( 0 == b2_ranks[2] ); 
  
  gkyl_rrobin_decomp_release(rr);
}

TEST_LIST = {
  { "test_rrobin_decomp_1_ho", test_rrobin_decomp_1_ho },
  { "test_rrobin_decomp_2_ho", test_rrobin_decomp_2_ho },
  { "test_rrobin_decomp_3_ho", test_rrobin_decomp_3_ho },
  { "test_rrobin_decomp_4_ho", test_rrobin_decomp_4_ho },
  { "test_rrobin_decomp_5_ho", test_rrobin_decomp_5_ho },
  { "test_rrobin_decomp_6_ho", test_rrobin_decomp_6_ho },
  { "test_rrobin_decomp_7_ho", test_rrobin_decomp_7_ho },
  { NULL, NULL },
};
