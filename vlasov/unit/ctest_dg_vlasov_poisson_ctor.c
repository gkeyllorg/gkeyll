// Tests for the Vlasov-Poisson DG equation object constructor.
//
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_dg_vlasov_poisson.h>
#include <gkyl_dg_vlasov_poisson_priv.h>

static void
check_common(struct gkyl_dg_eqn *eqn, int cdim, int pdim, long cvol, long pvol)
{
  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->vol_term != 0 );
  TEST_CHECK( eqn->surf_term != 0 );
  TEST_CHECK( eqn->boundary_surf_term != 0 );
  TEST_CHECK( eqn->on_dev == eqn );

  struct dg_vlasov_poisson *vp = container_of(eqn, struct dg_vlasov_poisson, eqn);
  TEST_CHECK( vp->cdim == cdim );
  TEST_CHECK( vp->pdim == pdim );
  TEST_CHECK( vp->conf_range.volume == cvol );
  TEST_CHECK( vp->phase_range.volume == pvol );
  TEST_CHECK( vp->auxfields.potentials == 0 );
  TEST_CHECK( vp->auxfields.fields_ext == 0 );
}

void
test_vp_phi()
{
  int cdim = 1, vdim = 1, pdim = cdim+vdim;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 2);
  gkyl_cart_modal_serendip(&pbasis, pdim, 2);

  struct gkyl_range crange, prange;
  gkyl_range_init_from_shape(&crange, cdim, (int[]) { 10 });
  gkyl_range_init_from_shape(&prange, pdim, (int[]) { 10, 8 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_vlasov_poisson_new(&cbasis, &pbasis, &crange, &prange,
    GKYL_MODEL_DEFAULT, GKYL_FIELD_PHI, false);

  check_common(eqn, cdim, pdim, 10, 80);

  gkyl_dg_eqn_release(eqn);
}

void
test_vp_phi_ext_potentials()
{
  int cdim = 1, vdim = 1, pdim = cdim+vdim;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 2);
  gkyl_cart_modal_serendip(&pbasis, pdim, 2);

  struct gkyl_range crange, prange;
  gkyl_range_init_from_shape(&crange, cdim, (int[]) { 6 });
  gkyl_range_init_from_shape(&prange, pdim, (int[]) { 6, 6 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_vlasov_poisson_new(&cbasis, &pbasis, &crange, &prange,
    GKYL_MODEL_DEFAULT, GKYL_FIELD_PHI_EXT_POTENTIALS, false);

  check_common(eqn, cdim, pdim, 6, 36);

  gkyl_dg_eqn_release(eqn);
}

void
test_vp_phi_ext_fields()
{
  int cdim = 1, vdim = 1, pdim = cdim+vdim;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 2);
  gkyl_cart_modal_serendip(&pbasis, pdim, 2);

  struct gkyl_range crange, prange;
  gkyl_range_init_from_shape(&crange, cdim, (int[]) { 4 });
  gkyl_range_init_from_shape(&prange, pdim, (int[]) { 4, 4 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_vlasov_poisson_new(&cbasis, &pbasis, &crange, &prange,
    GKYL_MODEL_DEFAULT, GKYL_FIELD_PHI_EXT_FIELDS, false);

  check_common(eqn, cdim, pdim, 4, 16);

  gkyl_dg_eqn_release(eqn);
}

TEST_LIST = {
  { "vp_phi", test_vp_phi },
  { "vp_phi_ext_potentials", test_vp_phi_ext_potentials },
  { "vp_phi_ext_fields", test_vp_phi_ext_fields },
  { NULL, NULL },
};
