#include <gkyl_gk_collisionless_flux_kernels.h> 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_none(
    const double *w, const double *dxv,
    const double *vmap_prime_l, const double *vmap_prime_r,
    const double *vmap, const double *vmapSq, const double q_, const double m_, 
    const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
    const double *bmag, const double *phi, const double *apar, const double *apardot, 
    const double *JfL, const double *JfR, double* GKYL_RESTRICT flux_surf) 
{ 
  return 0.0; 
}
