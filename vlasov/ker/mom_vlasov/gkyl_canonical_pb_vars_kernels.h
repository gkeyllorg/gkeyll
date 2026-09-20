#pragma once

#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x1v_ser_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x1v_ser_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x1v_ser_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x1v_ser_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x2v_ser_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x2v_ser_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x2v_ser_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x2v_ser_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x3v_ser_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x3v_ser_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x3v_ser_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x3v_ser_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_2x2v_ser_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_2x2v_ser_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_2x2v_ser_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_2x2v_ser_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_2x3v_ser_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_2x3v_ser_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x1v_tensor_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x1v_tensor_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x1v_tensor_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x1v_tensor_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x2v_tensor_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x2v_tensor_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x2v_tensor_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x2v_tensor_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x3v_tensor_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x3v_tensor_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x3v_tensor_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_1x3v_tensor_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_2x2v_tensor_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_2x2v_tensor_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_2x2v_tensor_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_2x2v_tensor_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_2x3v_tensor_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_2x3v_tensor_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_3x3v_tensor_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov);
GKYL_CU_DH void canonical_pb_vars_pressure_3x3v_tensor_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, 
            const double *nv_i, double* GKYL_RESTRICT d_Jv_P);

EXTERN_C_END
