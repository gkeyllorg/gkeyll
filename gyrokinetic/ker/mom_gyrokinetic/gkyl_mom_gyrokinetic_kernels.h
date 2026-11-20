#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void mom_gyrokinetic_M0_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M1_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2PAR_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3PAR_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2M3_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const
double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M1_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2PAR_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3PAR_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2M3_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_gyrokinetic_M0_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M1_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2PAR_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3PAR_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2M3_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const
double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M1_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2PAR_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3PAR_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2M3_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_gyrokinetic_M0_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M1_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2PAR_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3PAR_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2M3_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const
double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2PERP_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3PERP_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M1_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2PAR_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3PAR_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2M3_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2PERP_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3PERP_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_gyrokinetic_M0_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M1_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2PAR_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3PAR_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2M3_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const
double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2PERP_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3PERP_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M1_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2PAR_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3PAR_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2M3_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2PERP_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3PERP_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_gyrokinetic_M0_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M1_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2PAR_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3PAR_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2M3_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const
double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2PERP_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3PERP_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M1_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2PAR_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3PAR_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2M3_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2PERP_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3PERP_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_gyrokinetic_M0_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M1_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2PAR_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3PAR_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2M3_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const
double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M1_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2PAR_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3PAR_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2M3_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_gyrokinetic_M0_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M1_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2PAR_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3PAR_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2M3_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_HAMILTONIAN_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const
double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M2PERP_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M3PERP_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2PARM2PERP_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M1_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2PAR_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3PAR_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2M3_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_HAMILTONIAN_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M2PERP_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M3PERP_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERP_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out); 

EXTERN_C_END 
