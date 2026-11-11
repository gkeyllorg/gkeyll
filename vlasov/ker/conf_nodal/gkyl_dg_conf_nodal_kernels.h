#pragma once 

// Generated (gkylcas/vlasov/conf_nodal/conf-nodal.mac): conversions between modal expansions and 
// the configuration-nodal, velocity-modal form. The nodes are the 2^cdim interior Gauss-Legendre 
// points of a p=1 configuration-space cell, node q having bit d set when x_d > 0. 
#include <gkyl_util.h> 

EXTERN_C_BEG 

GKYL_CU_DH void dg_conf_nodal_conf_to_nodal_1x_ser_p1(const double *fmodal, double* GKYL_RESTRICT fnodal); 
GKYL_CU_DH void dg_conf_nodal_conf_to_modal_1x_ser_p1(const double *fnodal, double* GKYL_RESTRICT fmodal); 
GKYL_CU_DH void dg_conf_nodal_phase_to_nodal_1x1v_hyb_p1(const double *fmodal, double* GKYL_RESTRICT fnodal); 
GKYL_CU_DH void dg_conf_nodal_phase_to_modal_accumulate_1x1v_hyb_p1(const double *fnodal, double* GKYL_RESTRICT fmodal); 

GKYL_CU_DH void dg_conf_nodal_conf_to_nodal_2x_ser_p1(const double *fmodal, double* GKYL_RESTRICT fnodal); 
GKYL_CU_DH void dg_conf_nodal_conf_to_modal_2x_ser_p1(const double *fnodal, double* GKYL_RESTRICT fmodal); 
GKYL_CU_DH void dg_conf_nodal_phase_to_nodal_2x1v_hyb_p1(const double *fmodal, double* GKYL_RESTRICT fnodal); 
GKYL_CU_DH void dg_conf_nodal_phase_to_modal_accumulate_2x1v_hyb_p1(const double *fnodal, double* GKYL_RESTRICT fmodal); 

GKYL_CU_DH void dg_conf_nodal_conf_to_nodal_3x_ser_p1(const double *fmodal, double* GKYL_RESTRICT fnodal); 
GKYL_CU_DH void dg_conf_nodal_conf_to_modal_3x_ser_p1(const double *fnodal, double* GKYL_RESTRICT fmodal); 
GKYL_CU_DH void dg_conf_nodal_phase_to_nodal_3x1v_hyb_p1(const double *fmodal, double* GKYL_RESTRICT fnodal); 
GKYL_CU_DH void dg_conf_nodal_phase_to_modal_accumulate_3x1v_hyb_p1(const double *fnodal, double* GKYL_RESTRICT fmodal); 

EXTERN_C_END 
