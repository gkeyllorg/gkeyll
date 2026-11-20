#pragma once 

// Generated (gkylcas/pkpm/pkpm/pkpm-lbo-nodal.mac): the Lenard Bernstein operator of the PKPM 
// distributions at one configuration-space node, in the parallel velocity alone (p=2). 
#include <gkyl_util.h> 

EXTERN_C_BEG 

GKYL_CU_DH void lbo_pkpm_nodal_vol_1v_p2(const double *w, const double *dxv, double nu, double nu_u, double nu_vtsq, int stride, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void lbo_pkpm_nodal_surfvpar_1v_p2(const double *w, const double *dxv, double nu, double nu_u, double nu_vtsq, int stride, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void lbo_pkpm_nodal_boundary_surfvpar_1v_p2(const double *w, const double *dxv, double nu, double nu_u, double nu_vtsq, int stride, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

EXTERN_C_END 
