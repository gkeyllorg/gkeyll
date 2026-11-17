#pragma once 
#include <math.h> 
#include <stdbool.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

// p=1 serendipity in 1x: values at the 2 Gauss-Legendre nodes, and back. 
GKYL_CU_DH static inline void fluid_em_coupling_1x_ser_p1_nodal(const double *f, double *fq) 
{ 
  fq[0] = 0.7071067811865475*f[0]-0.7071067811865475*f[1]; 
  fq[1] = 0.7071067811865475*f[1]+0.7071067811865475*f[0]; 
} 
GKYL_CU_DH static inline void fluid_em_coupling_1x_ser_p1_modal(const double *fq, double *f) 
{ 
  f[0] = 0.7071067811865475*fq[1]+0.7071067811865475*fq[0]; 
  f[1] = 0.7071067811865475*fq[1]-0.7071067811865475*fq[0]; 
} 
GKYL_CU_DH void fluid_em_coupling_set_1x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES],
  double epsilon0, bool is_field_static, double dt, struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n,
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current,
  const double *fluid[GKYL_MAX_SPECIES], const double *em); 
GKYL_CU_DH void fluid_em_coupling_copy_1x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES],
  bool has_energy[GKYL_MAX_SPECIES], double epsilon0, struct gkyl_nmat *x,
  double* GKYL_RESTRICT fluid[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 

// p=1 serendipity in 2x: values at the 4 Gauss-Legendre nodes, and back. 
GKYL_CU_DH static inline void fluid_em_coupling_2x_ser_p1_nodal(const double *f, double *fq) 
{ 
  fq[0] = 0.5*f[3]-0.5*f[2]-0.5*f[1]+0.5*f[0]; 
  fq[1] = -(0.5*f[3])+0.5*f[2]-0.5*f[1]+0.5*f[0]; 
  fq[2] = -(0.5*f[3])-0.5*f[2]+0.5*f[1]+0.5*f[0]; 
  fq[3] = 0.5*f[3]+0.5*f[2]+0.5*f[1]+0.5*f[0]; 
} 
GKYL_CU_DH static inline void fluid_em_coupling_2x_ser_p1_modal(const double *fq, double *f) 
{ 
  f[0] = 0.5*fq[3]+0.5*fq[2]+0.5*fq[1]+0.5*fq[0]; 
  f[1] = 0.5*fq[3]+0.5*fq[2]-0.5*fq[1]-0.5*fq[0]; 
  f[2] = 0.5*fq[3]-0.5*fq[2]+0.5*fq[1]-0.5*fq[0]; 
  f[3] = 0.5*fq[3]-0.5*fq[2]-0.5*fq[1]+0.5*fq[0]; 
} 
GKYL_CU_DH void fluid_em_coupling_set_2x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES],
  double epsilon0, bool is_field_static, double dt, struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n,
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current,
  const double *fluid[GKYL_MAX_SPECIES], const double *em); 
GKYL_CU_DH void fluid_em_coupling_copy_2x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES],
  bool has_energy[GKYL_MAX_SPECIES], double epsilon0, struct gkyl_nmat *x,
  double* GKYL_RESTRICT fluid[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 

// p=1 serendipity in 3x: values at the 8 Gauss-Legendre nodes, and back. 
GKYL_CU_DH static inline void fluid_em_coupling_3x_ser_p1_nodal(const double *f, double *fq) 
{ 
  fq[0] = -(0.3535533905932737*f[7])+0.3535533905932737*f[6]+0.3535533905932737*f[5]+0.3535533905932737*f[4]-0.3535533905932737*f[3]-0.3535533905932737*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
  fq[1] = 0.3535533905932737*f[7]-0.3535533905932737*f[6]-0.3535533905932737*f[5]+0.3535533905932737*f[4]+0.3535533905932737*f[3]-0.3535533905932737*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
  fq[2] = 0.3535533905932737*f[7]-0.3535533905932737*f[6]+0.3535533905932737*f[5]-0.3535533905932737*f[4]-0.3535533905932737*f[3]+0.3535533905932737*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
  fq[3] = -(0.3535533905932737*f[7])+0.3535533905932737*f[6]-0.3535533905932737*f[5]-0.3535533905932737*f[4]+0.3535533905932737*f[3]+0.3535533905932737*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
  fq[4] = 0.3535533905932737*f[7]+0.3535533905932737*f[6]-0.3535533905932737*f[5]-0.3535533905932737*f[4]-0.3535533905932737*f[3]-0.3535533905932737*f[2]+0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
  fq[5] = -(0.3535533905932737*f[7])-0.3535533905932737*f[6]+0.3535533905932737*f[5]-0.3535533905932737*f[4]+0.3535533905932737*f[3]-0.3535533905932737*f[2]+0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
  fq[6] = -(0.3535533905932737*f[7])-0.3535533905932737*f[6]-0.3535533905932737*f[5]+0.3535533905932737*f[4]-0.3535533905932737*f[3]+0.3535533905932737*f[2]+0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
  fq[7] = 0.3535533905932737*f[7]+0.3535533905932737*f[6]+0.3535533905932737*f[5]+0.3535533905932737*f[4]+0.3535533905932737*f[3]+0.3535533905932737*f[2]+0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
} 
GKYL_CU_DH static inline void fluid_em_coupling_3x_ser_p1_modal(const double *fq, double *f) 
{ 
  f[0] = 0.3535533905932737*fq[7]+0.3535533905932737*fq[6]+0.3535533905932737*fq[5]+0.3535533905932737*fq[4]+0.3535533905932737*fq[3]+0.3535533905932737*fq[2]+0.3535533905932737*fq[1]+0.3535533905932737*fq[0]; 
  f[1] = 0.3535533905932737*fq[7]+0.3535533905932737*fq[6]+0.3535533905932737*fq[5]+0.3535533905932737*fq[4]-0.3535533905932737*fq[3]-0.3535533905932737*fq[2]-0.3535533905932737*fq[1]-0.3535533905932737*fq[0]; 
  f[2] = 0.3535533905932737*fq[7]+0.3535533905932737*fq[6]-0.3535533905932737*fq[5]-0.3535533905932737*fq[4]+0.3535533905932737*fq[3]+0.3535533905932737*fq[2]-0.3535533905932737*fq[1]-0.3535533905932737*fq[0]; 
  f[3] = 0.3535533905932737*fq[7]-0.3535533905932737*fq[6]+0.3535533905932737*fq[5]-0.3535533905932737*fq[4]+0.3535533905932737*fq[3]-0.3535533905932737*fq[2]+0.3535533905932737*fq[1]-0.3535533905932737*fq[0]; 
  f[4] = 0.3535533905932737*fq[7]+0.3535533905932737*fq[6]-0.3535533905932737*fq[5]-0.3535533905932737*fq[4]-0.3535533905932737*fq[3]-0.3535533905932737*fq[2]+0.3535533905932737*fq[1]+0.3535533905932737*fq[0]; 
  f[5] = 0.3535533905932737*fq[7]-0.3535533905932737*fq[6]+0.3535533905932737*fq[5]-0.3535533905932737*fq[4]-0.3535533905932737*fq[3]+0.3535533905932737*fq[2]-0.3535533905932737*fq[1]+0.3535533905932737*fq[0]; 
  f[6] = 0.3535533905932737*fq[7]-0.3535533905932737*fq[6]-0.3535533905932737*fq[5]+0.3535533905932737*fq[4]+0.3535533905932737*fq[3]-0.3535533905932737*fq[2]-0.3535533905932737*fq[1]+0.3535533905932737*fq[0]; 
  f[7] = 0.3535533905932737*fq[7]-0.3535533905932737*fq[6]-0.3535533905932737*fq[5]+0.3535533905932737*fq[4]-0.3535533905932737*fq[3]+0.3535533905932737*fq[2]+0.3535533905932737*fq[1]-0.3535533905932737*fq[0]; 
} 
GKYL_CU_DH void fluid_em_coupling_set_3x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES],
  double epsilon0, bool is_field_static, double dt, struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n,
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current,
  const double *fluid[GKYL_MAX_SPECIES], const double *em); 
GKYL_CU_DH void fluid_em_coupling_copy_3x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES],
  bool has_energy[GKYL_MAX_SPECIES], double epsilon0, struct gkyl_nmat *x,
  double* GKYL_RESTRICT fluid[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 

EXTERN_C_END 
