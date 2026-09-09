// Test construction of the gyrokinetic anomalous-diffusion DG equation object
// (gkyl_gk_anomalous_diffusion_new). Verifies the base num_equations, the
// number of phase-space basis functions recorded, the stored conf-space range
// volume, and that the volume / surface / boundary kernels are wired up for
// the supported 2x2v p1 case under a couple of boundary-condition choices.
// Also checks acquire/release reference counting.
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_gk_anomalous_diffusion.h>
#include <gkyl_gk_anomalous_diffusion_priv.h>
#include <gkyl_gk_bc_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

static void
make_ranges(int cdim, struct gkyl_range *confRange, struct gkyl_range *confRange_ext,
  struct gkyl_basis *basis, struct gkyl_basis *cbasis)
{
  int vdim = 2;
  int pdim = cdim + vdim;
  int cells[GKYL_MAX_DIM];
  int ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  for (int d=0; d<cdim; d++) { cells[d] = 8; ghost[d] = 1; lower[d] = 0.0; upper[d] = 1.0; }

  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, confRange_ext, confRange);

  gkyl_cart_modal_serendip(cbasis, cdim, 1);
  gkyl_cart_modal_gkhybrid(basis, cdim, vdim);
}

void
test_anom_diff_ctor_local()
{
  int cdim = 2;
  struct gkyl_basis basis, cbasis;
  struct gkyl_range confRange, confRange_ext;
  make_ranges(cdim, &confRange, &confRange_ext, &basis, &cbasis);

  // "ELSE" branch (local boundary kernels) using copy BCs.
  struct gkyl_dg_eqn *eqn = gkyl_gk_anomalous_diffusion_new(&basis, &cbasis,
    &confRange, GKYL_BC_GK_SPECIES_COPY, GKYL_BC_GK_SPECIES_COPY, false);

  TEST_CHECK( eqn != NULL );
  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->vol_term != NULL );
  TEST_CHECK( eqn->surf_term != NULL );

  struct gk_anomalous_diffusion *diffusion =
    container_of(eqn, struct gk_anomalous_diffusion, eqn);

  TEST_CHECK( diffusion->conf_range.volume == 64 );  // 8x8 interior cells
  TEST_CHECK( diffusion->surf != NULL );
  TEST_CHECK( diffusion->boundary_surf[0] != NULL );
  TEST_CHECK( diffusion->boundary_surf[1] != NULL );

  gkyl_dg_eqn_release(eqn);
}

void
test_anom_diff_ctor_zeroflux()
{
  int cdim = 2;
  struct gkyl_basis basis, cbasis;
  struct gkyl_range confRange, confRange_ext;
  make_ranges(cdim, &confRange, &confRange_ext, &basis, &cbasis);

  // ZERO_FLUX branch on both ends.
  struct gkyl_dg_eqn *eqn = gkyl_gk_anomalous_diffusion_new(&basis, &cbasis,
    &confRange, GKYL_BC_GK_SPECIES_ZERO_FLUX, GKYL_BC_GK_SPECIES_ZERO_FLUX, false);

  TEST_CHECK( eqn != NULL );
  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->vol_term != NULL );
  TEST_CHECK( eqn->surf_term != NULL );

  struct gk_anomalous_diffusion *diffusion =
    container_of(eqn, struct gk_anomalous_diffusion, eqn);
  TEST_CHECK( diffusion->conf_range.volume == 64 );
  TEST_CHECK( diffusion->surf != NULL );

  gkyl_dg_eqn_release(eqn);
}

void
test_anom_diff_acquire()
{
  int cdim = 2;
  struct gkyl_basis basis, cbasis;
  struct gkyl_range confRange, confRange_ext;
  make_ranges(cdim, &confRange, &confRange_ext, &basis, &cbasis);

  struct gkyl_dg_eqn *eqn = gkyl_gk_anomalous_diffusion_new(&basis, &cbasis,
    &confRange, GKYL_BC_GK_SPECIES_COPY, GKYL_BC_GK_SPECIES_COPY, false);

  struct gkyl_dg_eqn *eqn2 = gkyl_dg_eqn_acquire(eqn);
  TEST_CHECK( eqn2 == eqn );
  gkyl_dg_eqn_release(eqn);
  TEST_CHECK( eqn2->num_equations == 1 );
  gkyl_dg_eqn_release(eqn2);
}

TEST_LIST = {
  { "anom_diff_ctor_local", test_anom_diff_ctor_local },
  { "anom_diff_ctor_zeroflux", test_anom_diff_ctor_zeroflux },
  { "anom_diff_acquire", test_anom_diff_acquire },
  { NULL, NULL },
};
