#pragma once

#include <math.h>
#include <gkyl_util.h>
#include <gkyl_dg_geom.h>
#include <gkyl_gk_dg_geom.h>

EXTERN_C_BEG

GKYL_CU_DH double gk_collisionless_passive_flux_surfx_1x1v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_edge_surfx_1x1v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_multib_boundary_surfx_1x1v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_multib_boundary_edge_surfx_1x1v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 

GKYL_CU_DH double gk_collisionless_passive_flux_surfx_1x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_edge_surfx_1x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_multib_boundary_surfx_1x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_multib_boundary_edge_surfx_1x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 

GKYL_CU_DH double gk_collisionless_passive_flux_surfx_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_edge_surfx_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_surfy_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_edge_surfy_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_multib_boundary_surfx_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_multib_boundary_edge_surfx_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_multib_boundary_surfy_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_multib_boundary_edge_surfy_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 

GKYL_CU_DH double gk_collisionless_passive_flux_surfx_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_edge_surfx_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_surfy_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_edge_surfy_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_surfz_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_edge_surfz_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_multib_boundary_surfx_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_multib_boundary_edge_surfx_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_multib_boundary_surfy_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_multib_boundary_edge_surfy_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_multib_boundary_surfz_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_, 
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gk_collisionless_passive_flux_multib_boundary_edge_surfz_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
                const double *speedsL, const double *speedsR, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 


EXTERN_C_END
