#pragma once

#include <math.h>
#include <gkyl_util.h>
#include <gkyl_dg_geom.h>
#include <gkyl_gk_dg_geom.h>

#define gsign(x) (((x) > 0) - ((x) < 0))

EXTERN_C_BEG

GKYL_CU_DH double gk_collisionless_flux_surfx_1x1v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_periodic_1x1v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_periodic_1x1v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_nonperiodic_1x1v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_nonperiodic_1x1v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_multib_1x1v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_multib_1x1v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_1x1v_ser_p1(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *yfield, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 

GKYL_CU_DH double gk_collisionless_flux_surfx_1x1v_ser_p2(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_periodic_1x1v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_periodic_1x1v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_nonperiodic_1x1v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_nonperiodic_1x1v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_multib_1x1v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_multib_1x1v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_1x1v_ser_p2(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *yfield, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 

GKYL_CU_DH double gk_collisionless_flux_surfx_1x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_periodic_1x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_periodic_1x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_nonperiodic_1x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_nonperiodic_1x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_multib_1x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_multib_1x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_1x2v_ser_p1(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *yfield, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 

GKYL_CU_DH double gk_collisionless_flux_surfx_1x2v_ser_p2(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_periodic_1x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_periodic_1x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_nonperiodic_1x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_nonperiodic_1x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_multib_1x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_multib_1x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_1x2v_ser_p2(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *yfield, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 

GKYL_CU_DH double gk_collisionless_flux_surfx_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfy_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_periodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_periodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_lower_periodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_upper_periodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_nonperiodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_nonperiodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_lower_nonperiodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_upper_nonperiodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_multib_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_multib_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_lower_multib_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_upper_multib_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_2x2v_ser_p1(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *yfield, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfx_no_by_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfy_no_by_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_lower_periodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_upper_periodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_lower_periodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_upper_periodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_lower_nonperiodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_upper_nonperiodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_lower_nonperiodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_upper_nonperiodic_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_lower_multib_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_upper_multib_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_lower_multib_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_upper_multib_2x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_no_by_2x2v_ser_p1(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *yfield, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 

GKYL_CU_DH double gk_collisionless_flux_surfx_2x2v_ser_p2(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfy_2x2v_ser_p2(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_periodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_periodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_lower_periodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_upper_periodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_nonperiodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_nonperiodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_lower_nonperiodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_upper_nonperiodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_multib_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_multib_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_lower_multib_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_upper_multib_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_2x2v_ser_p2(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *yfield, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfx_no_by_2x2v_ser_p2(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfy_no_by_2x2v_ser_p2(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_lower_periodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_upper_periodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_lower_periodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_upper_periodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_lower_nonperiodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_upper_nonperiodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_lower_nonperiodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_upper_nonperiodic_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_lower_multib_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_upper_multib_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_lower_multib_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_upper_multib_2x2v_ser_p2(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_no_by_2x2v_ser_p2(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *yfield, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 

GKYL_CU_DH double gk_collisionless_flux_surfx_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfy_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfz_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_periodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_periodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_lower_periodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_upper_periodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfz_lower_periodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfz_upper_periodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_nonperiodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_nonperiodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_lower_nonperiodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_upper_nonperiodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfz_lower_nonperiodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfz_upper_nonperiodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_multib_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_multib_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_lower_multib_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_upper_multib_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfz_lower_multib_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfz_upper_multib_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_3x2v_ser_p1(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *yfield, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfx_no_by_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfy_no_by_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfz_no_by_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_lower_periodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_upper_periodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_lower_periodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_upper_periodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfz_no_by_lower_periodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfz_no_by_upper_periodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_lower_nonperiodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_upper_nonperiodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_lower_nonperiodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_upper_nonperiodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfz_no_by_lower_nonperiodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfz_no_by_upper_nonperiodic_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_lower_multib_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_no_by_upper_multib_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_lower_multib_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfy_no_by_upper_multib_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfz_no_by_lower_multib_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfz_no_by_upper_multib_3x2v_ser_p1(
                    const double *w, const double *dxv,
                    
                    const double *vmap, const double *vmapSq, const double q_, const double m_,
                    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                    const double *phiL, const double *phiR, const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR, 
                    double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_no_by_3x2v_ser_p1(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *yfield, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 


EXTERN_C_END
