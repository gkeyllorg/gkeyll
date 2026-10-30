#pragma once

#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double canonical_pb_vol_2x_ser_p1(const double *w, const double *dxv, const double *phi, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_vol_2x_ser_p1(const double *w, const double *dxv, const double *phi, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_hasegawa_mima_source_2x_ser_p1(const double *dxv, double alpha, const double *phi, const double *n0, const double *adiabatic_coupling_phi_n, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_hasegawa_wakatani_source_2x_ser_p1(const double *dxv, double alpha, const double *phi, const double *n0, const double *adiabatic_coupling_phi_n, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_subtract_zonal_2x_ser_p1(const double *phi_zonal, const double *n_zonal, double* GKYL_RESTRICT adiabatic_coupling_phi_n); 
GKYL_CU_DH int canonical_pb_alpha_surfx_2x_ser_p1(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x_ser_p1(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_2x_ser_p1(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_surfx_2x_ser_p1(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x_ser_p1(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x_ser_p1(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_2x_ser_p1(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_surfy_2x_ser_p1(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_2x_ser_p2(const double *w, const double *dxv, const double *phi, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_vol_2x_ser_p2(const double *w, const double *dxv, const double *phi, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_hasegawa_mima_source_2x_ser_p2(const double *dxv, double alpha, const double *phi, const double *n0, const double *adiabatic_coupling_phi_n, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_hasegawa_wakatani_source_2x_ser_p2(const double *dxv, double alpha, const double *phi, const double *n0, const double *adiabatic_coupling_phi_n, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_subtract_zonal_2x_ser_p2(const double *phi_zonal, const double *n_zonal, double* GKYL_RESTRICT adiabatic_coupling_phi_n); 
GKYL_CU_DH int canonical_pb_alpha_surfx_2x_ser_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x_ser_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_2x_ser_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_surfx_2x_ser_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x_ser_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x_ser_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_2x_ser_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_surfy_2x_ser_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_2x_tensor_p2(const double *w, const double *dxv, const double *phi, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_vol_2x_tensor_p2(const double *w, const double *dxv, const double *phi, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_hasegawa_mima_source_2x_tensor_p2(const double *dxv, double alpha, const double *phi, const double *n0, const double *adiabatic_coupling_phi_n, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_hasegawa_wakatani_source_2x_tensor_p2(const double *dxv, double alpha, const double *phi, const double *n0, const double *adiabatic_coupling_phi_n, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_fluid_subtract_zonal_2x_tensor_p2(const double *phi_zonal, const double *n_zonal, double* GKYL_RESTRICT adiabatic_coupling_phi_n); 
GKYL_CU_DH int canonical_pb_alpha_surfx_2x_tensor_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x_tensor_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_2x_tensor_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_surfx_2x_tensor_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x_tensor_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x_tensor_p2(const double *w, const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_2x_tensor_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_two_fluid_surfy_2x_tensor_p2(const double *w, const double *dxv, const double *phi, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void canonical_pb_fluid_default_source(const double *dxv, double alpha, const double *phi, const double *n0, const double *f, double* GKYL_RESTRICT out); 

EXTERN_C_END
