#include <gkyl_pkpm_fluid_kernels.h> 
GKYL_CU_DH double pkpm_fluid_surfy_2x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out) 
{ 
  // flux_l/flux_r: Stored fluxes [F_rho, F_x, F_y, F_z, F_E] on the lower/upper face of the cell. 
  // out:           Incremented increment of [rho, rho ux, rho uy, rho uz, E]. 
  const double dx1 = 2.0/dxv[1]; 
  const double *gl0 = &flux_l[0]; 
  const double *gr0 = &flux_r[0]; 
  double *out0 = &out[0]; 
  const double *gl1 = &flux_l[2]; 
  const double *gr1 = &flux_r[2]; 
  double *out1 = &out[4]; 
  const double *gl2 = &flux_l[4]; 
  const double *gr2 = &flux_r[4]; 
  double *out2 = &out[8]; 
  const double *gl3 = &flux_l[6]; 
  const double *gr3 = &flux_r[6]; 
  double *out3 = &out[12]; 
  const double *gl4 = &flux_l[8]; 
  const double *gr4 = &flux_r[8]; 
  double *out4 = &out[16]; 
  out0[0] += (0.7071067811865475*gl0[0]-0.7071067811865475*gr0[0])*dx1; 
  out0[1] += (0.7071067811865475*gl0[1]-0.7071067811865475*gr0[1])*dx1; 
  out0[2] += (-(1.224744871391589*gr0[0])-1.224744871391589*gl0[0])*dx1; 
  out0[3] += (-(1.224744871391589*gr0[1])-1.224744871391589*gl0[1])*dx1; 
  out1[0] += (0.7071067811865475*gl1[0]-0.7071067811865475*gr1[0])*dx1; 
  out1[1] += (0.7071067811865475*gl1[1]-0.7071067811865475*gr1[1])*dx1; 
  out1[2] += (-(1.224744871391589*gr1[0])-1.224744871391589*gl1[0])*dx1; 
  out1[3] += (-(1.224744871391589*gr1[1])-1.224744871391589*gl1[1])*dx1; 
  out2[0] += (0.7071067811865475*gl2[0]-0.7071067811865475*gr2[0])*dx1; 
  out2[1] += (0.7071067811865475*gl2[1]-0.7071067811865475*gr2[1])*dx1; 
  out2[2] += (-(1.224744871391589*gr2[0])-1.224744871391589*gl2[0])*dx1; 
  out2[3] += (-(1.224744871391589*gr2[1])-1.224744871391589*gl2[1])*dx1; 
  out3[0] += (0.7071067811865475*gl3[0]-0.7071067811865475*gr3[0])*dx1; 
  out3[1] += (0.7071067811865475*gl3[1]-0.7071067811865475*gr3[1])*dx1; 
  out3[2] += (-(1.224744871391589*gr3[0])-1.224744871391589*gl3[0])*dx1; 
  out3[3] += (-(1.224744871391589*gr3[1])-1.224744871391589*gl3[1])*dx1; 
  out4[0] += (0.7071067811865475*gl4[0]-0.7071067811865475*gr4[0])*dx1; 
  out4[1] += (0.7071067811865475*gl4[1]-0.7071067811865475*gr4[1])*dx1; 
  out4[2] += (-(1.224744871391589*gr4[0])-1.224744871391589*gl4[0])*dx1; 
  out4[3] += (-(1.224744871391589*gr4[1])-1.224744871391589*gl4[1])*dx1; 
  return 0.0; 
} 
