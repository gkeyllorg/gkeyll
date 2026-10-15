#include <gkyl_gpython_gyrokinetic.h>
#include <gkyl_gpython_core_priv.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_mom_calc.h>
#include <gkyl_velocity_map.h>
#include <gkyl_array_ops.h>
#include <math.h>
#include <string.h>

int
gpython_gyrokinetic_api_version(void)
{
  return GPYTHON_GYROKINETIC_API_VERSION;
}

gpython_array *
gpython_mom_gyrokinetic(
  const gpython_basis *cbasis, const gpython_basis *pbasis, const double *lower,
  const double *upper, const int *cells, double mass, double charge,
  const gpython_array *distribution, const gpython_array *bmag, const gpython_array *phi,
  const char *moment
)
{
  const struct gkyl_basis *cb = CBAS(cbasis), *pb = CBAS(pbasis);
  int cdim = cb->ndim, pdim = pb->ndim, vdim = pdim - cdim, order = cb->poly_order;
  if (cdim < 1 || cdim > 3 || vdim < 1 || vdim > 2 || (vdim == 1 && cdim != 1) ||
      cb->b_type != GKYL_BASIS_MODAL_SERENDIPITY || pb->poly_order != order ||
      (order != 1 && order != 2) || (cdim == 3 && order != 1) ||
      (order == 1 && pb->b_type != GKYL_BASIS_MODAL_GKHYBRID) ||
      (order == 2 && pb->b_type != GKYL_BASIS_MODAL_SERENDIPITY) || !isfinite(mass) || mass <= 0 ||
      !isfinite(charge)) {
    return NULL;
  }

  static const struct {
    const char *name;
    enum gkyl_distribution_moments id;
  } moments[] = {
    {"M0", GKYL_F_MOMENT_M0},
    {"M1", GKYL_F_MOMENT_M1},
    {"M2", GKYL_F_MOMENT_M2},
    {"M2par", GKYL_F_MOMENT_M2PAR},
    {"M2perp", GKYL_F_MOMENT_M2PERP},
    {"M3par", GKYL_F_MOMENT_M3PAR},
    {"M3perp", GKYL_F_MOMENT_M3PERP},
    {"ThreeMoments", GKYL_F_MOMENT_M0M1M2},
    {"FourMoments", GKYL_F_MOMENT_M0M1M2PARM2PERP},
    {"HamiltonianMoments", GKYL_F_MOMENT_HAMILTONIAN}
  };
  size_t choice = 0;
  for (; choice < sizeof(moments) / sizeof(moments[0]); ++choice) {
    if (strcmp(moment, moments[choice].name) == 0) {
      break;
    }
  }
  if (choice == sizeof(moments) / sizeof(moments[0])) {
    return NULL;
  }
  enum gkyl_distribution_moments id = moments[choice].id;
  if ((vdim == 1 && (id == GKYL_F_MOMENT_M2PERP || id == GKYL_F_MOMENT_M3PERP)) ||
      (id == GKYL_F_MOMENT_HAMILTONIAN && !phi)) {
    return NULL;
  }
  for (int d = 0; d < pdim; ++d) {
    if (cells[d] < 1 || !isfinite(lower[d]) || !isfinite(upper[d]) || lower[d] >= upper[d]) {
      return NULL;
    }
  }
  if (vdim == 2 && lower[cdim + 1] < 0) {
    return NULL;
  }

  struct gkyl_range crange, prange, vrange;
  gkyl_range_init_from_shape1(&crange, cdim, cells);
  gkyl_range_init_from_shape1(&prange, pdim, cells);
  gkyl_range_init_from_shape1(&vrange, vdim, cells + cdim);
  if (CARR(distribution)->ncomp != pb->num_basis || CARR(distribution)->size != prange.volume ||
      CARR(bmag)->ncomp != cb->num_basis || CARR(bmag)->size != crange.volume ||
      (phi && (CARR(phi)->ncomp != cb->num_basis || CARR(phi)->size != crange.volume))) {
    return NULL;
  }

  struct gkyl_rect_grid pgrid, cgrid, vgrid;
  gkyl_rect_grid_init(&pgrid, pdim, lower, upper, cells);
  gkyl_rect_grid_init(&cgrid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&vgrid, vdim, lower + cdim, upper + cdim, cells + cdim);
  struct gkyl_gk_geometry_inp inp = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .grid = cgrid,
    .basis = *cb,
    .local = crange,
    .local_ext = crange,
    .global = crange,
    .global_ext = crange,
  };
  struct gk_geometry *geom = gkyl_gk_geometry_new(NULL, &inp, false);
  gkyl_array_copy(geom->geo_corn.bmag, CARR(bmag));
  struct gkyl_mapc2p_inp identity = {0};
  struct gkyl_velocity_map *vmap =
    gkyl_velocity_map_new(identity, pgrid, vgrid, prange, prange, vrange, vrange, false);
  struct gkyl_mom_type *momt = gkyl_mom_gyrokinetic_new(
    cb, pb, &crange, mass, charge, vmap, geom, phi ? ARR(phi) : NULL, id, false
  );
  struct gkyl_array *out =
    gkyl_array_new(GKYL_DOUBLE, cb->num_basis * momt->num_mom, crange.volume);
  struct gkyl_mom_calc *calc = gkyl_mom_calc_new(&pgrid, momt, false);
  gkyl_mom_calc_advance(calc, &prange, &crange, CARR(distribution), out);
  gkyl_mom_calc_release(calc);
  gkyl_mom_type_release(momt);
  gkyl_velocity_map_release(vmap);
  gkyl_gk_geometry_release(geom);
  return (gpython_array *)out;
}
