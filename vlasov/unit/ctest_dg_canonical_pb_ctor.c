// Tests for the canonical-pb (Poisson bracket) DG equation object constructor.
// Serendipity basis requires poly_order == 2.
//
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_dg_canonical_pb.h>
#include <gkyl_dg_canonical_pb_priv.h>

void
test_canonical_pb_1x1v_p2()
{
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, 2);
  gkyl_cart_modal_serendip(&pbasis, 2, 2);

  struct gkyl_range prange;
  gkyl_range_init_from_shape(&prange, 2, (int[]) { 8, 16 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_canonical_pb_new(&cbasis, &pbasis, &prange, false);

  // canonical-pb is a scalar equation
  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->vol_term != 0 );
  TEST_CHECK( eqn->surf_term != 0 );
  TEST_CHECK( eqn->boundary_surf_term != 0 );
  TEST_CHECK( eqn->on_dev == eqn );

  struct dg_canonical_pb *can = container_of(eqn, struct dg_canonical_pb, eqn);
  TEST_CHECK( can->cdim == 1 );
  TEST_CHECK( can->pdim == 2 );
  TEST_CHECK( can->phase_range.volume == 8*16 );
  TEST_CHECK( can->phase_range.ndim == 2 );

  // streaming + acceleration surface kernels for the single config dim
  TEST_CHECK( can->stream_surf[0] != 0 );
  TEST_CHECK( can->stream_boundary_surf[0] != 0 );
  TEST_CHECK( can->accel_surf[0] != 0 );
  TEST_CHECK( can->accel_boundary_surf[0] != 0 );

  // aux fields start NULL
  TEST_CHECK( can->auxfields.hamil == 0 );
  TEST_CHECK( can->auxfields.alpha_surf == 0 );

  gkyl_dg_eqn_release(eqn);
}

void
test_canonical_pb_2x2v_p2()
{
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 2, 2);
  gkyl_cart_modal_serendip(&pbasis, 4, 2);

  struct gkyl_range prange;
  gkyl_range_init_from_shape(&prange, 4, (int[]) { 4, 4, 8, 8 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_canonical_pb_new(&cbasis, &pbasis, &prange, false);

  TEST_CHECK( eqn->num_equations == 1 );

  struct dg_canonical_pb *can = container_of(eqn, struct dg_canonical_pb, eqn);
  TEST_CHECK( can->cdim == 2 );
  TEST_CHECK( can->pdim == 4 );
  TEST_CHECK( can->phase_range.volume == 4*4*8*8 );

  // both config dims have streaming kernels
  TEST_CHECK( can->stream_surf[0] != 0 );
  TEST_CHECK( can->stream_surf[1] != 0 );
  TEST_CHECK( can->accel_surf[0] != 0 );
  TEST_CHECK( can->accel_surf[1] != 0 );

  gkyl_dg_eqn_release(eqn);
}

void
test_canonical_pb_1x2v_p2()
{
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, 2);
  gkyl_cart_modal_serendip(&pbasis, 3, 2);

  struct gkyl_range prange;
  gkyl_range_init_from_shape(&prange, 3, (int[]) { 6, 6, 6 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_canonical_pb_new(&cbasis, &pbasis, &prange, false);

  TEST_CHECK( eqn->num_equations == 1 );

  struct dg_canonical_pb *can = container_of(eqn, struct dg_canonical_pb, eqn);
  TEST_CHECK( can->cdim == 1 );
  TEST_CHECK( can->pdim == 3 );
  TEST_CHECK( can->phase_range.volume == 216 );
  TEST_CHECK( can->stream_surf[0] != 0 );

  gkyl_dg_eqn_release(eqn);
}

TEST_LIST = {
  { "canonical_pb_1x1v_p2", test_canonical_pb_1x1v_p2 },
  { "canonical_pb_2x2v_p2", test_canonical_pb_2x2v_p2 },
  { "canonical_pb_1x2v_p2", test_canonical_pb_1x2v_p2 },
  { NULL, NULL },
};
