#pragma once

#include <math.h>
#include <gkyl_util.h>
#include <gkyl_dg_geom.h>
#include <gkyl_gk_dg_geom.h>

EXTERN_C_BEG

GKYL_CU_DH double dg_gyrokinetic_passive_vol_1x1v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *speeds,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_gyrokinetic_passive_surfx_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_gyrokinetic_passive_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_gyrokinetic_passive_vol_1x2v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *speeds,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_gyrokinetic_passive_surfx_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_gyrokinetic_passive_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_gyrokinetic_passive_vol_2x2v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *speeds,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_gyrokinetic_passive_surfx_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_gyrokinetic_passive_boundary_surfx_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_gyrokinetic_passive_surfy_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_gyrokinetic_passive_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_gyrokinetic_passive_vol_3x2v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *speeds,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_gyrokinetic_passive_surfx_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_gyrokinetic_passive_boundary_surfx_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_gyrokinetic_passive_surfy_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_gyrokinetic_passive_boundary_surfy_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_gyrokinetic_passive_surfz_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_gyrokinetic_passive_boundary_surfz_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 


EXTERN_C_END
