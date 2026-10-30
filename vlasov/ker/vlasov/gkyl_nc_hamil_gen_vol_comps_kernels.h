#pragma once 
// These component kernels evaluate the volume-term Galerkin contraction 
// <grad(phi), alpha f> from a PRECOMPUTED modal phase-space alpha. They are 
// agnostic to the form of the Hamiltonian -- the alpha handed to them has 
// already been projected onto the dense phase basis -- hence the hamil_gen 
// naming. Every Hamiltonian flavor uses this machinery at 5D/6D; below 
// that the sparse (separable) Hamiltonian uses the inline exact-bracket 
// volume kernels instead. 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_alpha_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *vmap, const double *jacob_pos, const double *jacob_vel,
        const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp0_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp1_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp2_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp3_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp4_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp5_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp6_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp7_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp8_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp9_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp10_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp11_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp12_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp13_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp14_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp15_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp16_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp17_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp18_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp19_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp20_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp21_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp22_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp23_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp24_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp25_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp26_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp27_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp28_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp29_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp30_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp31_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp32_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp33_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp34_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp35_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp36_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp37_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp38_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp39_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp40_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp41_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp42_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp43_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp44_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp45_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp46_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp47_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_nc_hamil_phase_vol_alpha_1x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *vmap, const double *jacob_pos, const double *jacob_vel,
        const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_alpha_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *vmap, const double *jacob_pos, const double *jacob_vel,
        const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp0_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp1_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp2_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp3_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp4_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp5_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp6_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp7_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp8_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp9_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp10_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp11_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp12_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp13_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp14_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp15_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp16_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp17_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp18_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp19_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp20_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp21_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp22_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp23_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp24_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp25_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp26_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp27_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp28_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp29_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp30_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp31_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp32_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp33_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp34_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp35_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp36_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp37_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp38_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp39_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp40_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp41_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp42_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp43_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp44_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp45_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp46_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp47_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_nc_hamil_phase_vol_alpha_2x2v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *vmap, const double *jacob_pos, const double *jacob_vel,
        const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_alpha_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *vmap, const double *jacob_pos, const double *jacob_vel,
        const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp0_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp1_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp2_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp3_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp4_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp5_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp6_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp7_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp8_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp9_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp10_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp11_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp12_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp13_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp14_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp15_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp16_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp17_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp18_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp19_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp20_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp21_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp22_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp23_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp24_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp25_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp26_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp27_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp28_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp29_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp30_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp31_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_alpha_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *vmap, const double *jacob_pos, const double *jacob_vel,
        const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_nc_hamil_phase_vol_alpha_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *vmap, const double *jacob_pos, const double *jacob_vel,
        const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_alpha_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *vmap, const double *jacob_pos, const double *jacob_vel,
        const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp0_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp1_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp2_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp3_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp4_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp5_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp6_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp7_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp8_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp9_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp10_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp11_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp12_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp13_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp14_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp15_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp16_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp17_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp18_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp19_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp20_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp21_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp22_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp23_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp24_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp25_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp26_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp27_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp28_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp29_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp30_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp31_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp32_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp33_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp34_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp35_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp36_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp37_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp38_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp39_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp40_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp41_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp42_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp43_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp44_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp45_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp46_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp47_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp48_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp49_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp50_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp51_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp52_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp53_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp54_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp55_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp56_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp57_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp58_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp59_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp60_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp61_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp62_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp63_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp64_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp65_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp66_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp67_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp68_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp69_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp70_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp71_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp72_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp73_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp74_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp75_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp76_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp77_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp78_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp79_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp80_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp81_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp82_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp83_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp84_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp85_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp86_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp87_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp88_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp89_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp90_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp91_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp92_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp93_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp94_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp95_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp96_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp97_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp98_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp99_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp100_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp101_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp102_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp103_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp104_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp105_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp106_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp107_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp108_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp109_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp110_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp111_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_alpha_2x3v_ser_p2(const double *w, const double *dxv, const int dir,
        const double *vmap, const double *jacob_pos, const double *jacob_vel,
        const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out); 


GKYL_CU_DH void vlasov_nc_hamil_vel_dense_vol_alpha_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *vmap, const double *jacob_pos, const double *jacob_vel,
        const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp0_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp1_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp2_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp3_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp4_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp5_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp6_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp7_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp8_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp9_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp10_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp11_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp12_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp13_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp14_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp15_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp16_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp17_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp18_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp19_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp20_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp21_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp22_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp23_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp24_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp25_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp26_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp27_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp28_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp29_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp30_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp31_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp32_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp33_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp34_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp35_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp36_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp37_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp38_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp39_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp40_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp41_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp42_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp43_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp44_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp45_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp46_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp47_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp48_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp49_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp50_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp51_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp52_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp53_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp54_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp55_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp56_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp57_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp58_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp59_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp60_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp61_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp62_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp63_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *jacob_pos, const double *jacob_vel,
        const double *alpha, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_alpha_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *vmap, const double *jacob_pos, const double *jacob_vel,
        const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_nc_hamil_phase_vol_alpha_3x3v_ser_p1(const double *w, const double *dxv, const int dir,
        const double *vmap, const double *jacob_pos, const double *jacob_vel,
        const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out); 

EXTERN_C_END 
