#pragma once 
#include <math.h> 
#include <stdbool.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

// Resistive layer at x = sigma dt >= 0: g2 = exp(-x), phi = (1 - g2)/x, phi2 = (1 - phi)/x 
// (1, 1 and 1/2 at x = 0); for x < 1 phi2 from its Taylor series, where the closed form cancels. 
GKYL_CU_DH static inline void fluid_em_coupling_layer_coeffs(double x, double *g2, double *phi, double *phi2) 
{ 
  *g2 = exp(-x); 
  if (x < 1.0) { 
    *phi2 = 1.0/2.0 + x*(-1.0/6.0 + x*(1.0/24.0 + x*(-1.0/120.0 + x*(1.0/720.0 + x*(-1.0/5040.0 + x*(1.0/40320.0 + x*(-1.0/362880.0 + x*(1.0/3628800.0 + x*(-1.0/39916800.0 + x*(1.0/479001600.0 + x*(-1.0/6227020800.0 + x*(1.0/87178291200.0 + x*(-1.0/1307674368000.0 + x*(1.0/20922789888000.0 + x*(-1.0/355687428096000.0 + x*(1.0/6402373705728000.0)))))))))))))))); 
    *phi = 1.0 - x*(*phi2); 
  } 
  else { 
    *phi = -expm1(-x)/x; 
    *phi2 = (1.0 - *phi)/x; 
  } 
} 

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
  const double *sigma, const double *fluid[GKYL_MAX_SPECIES], const double *em); 
GKYL_CU_DH void fluid_em_coupling_copy_1x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES],
  bool has_energy[GKYL_MAX_SPECIES], double epsilon0, double dt, const double *sigma,
  struct gkyl_nmat *x, double* GKYL_RESTRICT fluid[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 

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
  const double *sigma, const double *fluid[GKYL_MAX_SPECIES], const double *em); 
GKYL_CU_DH void fluid_em_coupling_copy_2x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES],
  bool has_energy[GKYL_MAX_SPECIES], double epsilon0, double dt, const double *sigma,
  struct gkyl_nmat *x, double* GKYL_RESTRICT fluid[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 

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
  const double *sigma, const double *fluid[GKYL_MAX_SPECIES], const double *em); 
GKYL_CU_DH void fluid_em_coupling_copy_3x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES],
  bool has_energy[GKYL_MAX_SPECIES], double epsilon0, double dt, const double *sigma,
  struct gkyl_nmat *x, double* GKYL_RESTRICT fluid[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em); 

EXTERN_C_END 
