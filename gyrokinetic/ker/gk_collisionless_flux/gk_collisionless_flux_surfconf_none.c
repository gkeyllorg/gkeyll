#include <gkyl_gk_collisionless_flux_kernels.h> 
GKYL_CU_DH double gk_collisionless_flux_surfconf_none(
    const double *w, const double *dxv,
    const double *vmap, const double *vmapSq, const double q_, const double m_,
    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR, 
    const double *phi, const double *apar, const double *apardot, 
    const double *JfL, const double *JfR, double* GKYL_RESTRICT flux_surf) 
{ 
  return 0.0; 
}
