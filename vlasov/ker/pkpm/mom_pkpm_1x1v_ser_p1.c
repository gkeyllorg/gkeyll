#include <gkyl_mom_pkpm_kernels.h> 
GKYL_CU_DH void mom_pkpm_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2.0; 
  const double wvpar = w[1], dvpar = dxv[1]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 
  const double wvpar_cu = wvpar*wvpar*wvpar, dvpar_cu = dvpar*dvpar*dvpar; 

  const double *F_0 = &f[0]; 
  const double *G_1 = &f[6]; 
  const double *M0x = &f[12]; 
  const double *M0y = &f[18]; 
  const double *M0z = &f[24]; 
  const double *Nx = &f[30]; 
  const double *Ny = &f[36]; 
  const double *Nz = &f[42]; 
  out[0] += 1.4142135623730951*F_0[0]*mass*volFact; 
  out[1] += 1.4142135623730951*F_0[1]*mass*volFact; 
  out[2] += mass*volFact*(1.4142135623730951*F_0[0]*wvpar_sq+0.8164965809277261*F_0[2]*dvpar*wvpar+0.10540925533894596*F_0[4]*dvpar_sq+0.11785113019775789*F_0[0]*dvpar_sq); 
  out[3] += mass*volFact*(1.4142135623730951*F_0[1]*wvpar_sq+0.8164965809277261*F_0[3]*dvpar*wvpar+0.105409255338946*F_0[5]*dvpar_sq+0.11785113019775789*F_0[1]*dvpar_sq); 
  out[4] += 1.4142135623730951*G_1[0]*mass*volFact; 
  out[5] += 1.4142135623730951*G_1[1]*mass*volFact; 
  out[6] += mass*volFact*(1.4142135623730951*F_0[0]*wvpar+0.408248290463863*F_0[2]*dvpar); 
  out[7] += mass*volFact*(1.4142135623730951*F_0[1]*wvpar+0.408248290463863*F_0[3]*dvpar); 
  out[8] += mass*volFact*(1.224744871391589*F_0[2]*dvpar*wvpar_sq+1.4142135623730951*F_0[0]*wvpar_cu+0.3162277660168379*F_0[4]*dvpar_sq*wvpar+0.3535533905932737*F_0[0]*dvpar_sq*wvpar+0.06123724356957942*F_0[2]*dvpar_cu); 
  out[9] += mass*volFact*(1.224744871391589*F_0[3]*dvpar*wvpar_sq+1.4142135623730951*F_0[1]*wvpar_cu+0.31622776601683794*F_0[5]*dvpar_sq*wvpar+0.3535533905932737*F_0[1]*dvpar_sq*wvpar+0.06123724356957942*F_0[3]*dvpar_cu); 
  out[10] += mass*volFact*(1.4142135623730951*G_1[0]*wvpar+0.408248290463863*G_1[2]*dvpar); 
  out[11] += mass*volFact*(1.4142135623730951*G_1[1]*wvpar+0.408248290463863*G_1[3]*dvpar); 
  out[12] += 1.4142135623730951*M0x[0]*mass*volFact; 
  out[13] += 1.4142135623730951*M0x[1]*mass*volFact; 
  out[14] += 1.4142135623730951*M0y[0]*mass*volFact; 
  out[15] += 1.4142135623730951*M0y[1]*mass*volFact; 
  out[16] += 1.4142135623730951*M0z[0]*mass*volFact; 
  out[17] += 1.4142135623730951*M0z[1]*mass*volFact; 
  out[18] += mass*volFact*(1.4142135623730951*M0x[0]*wvpar+0.408248290463863*M0x[2]*dvpar); 
  out[19] += mass*volFact*(1.4142135623730951*M0x[1]*wvpar+0.408248290463863*M0x[3]*dvpar); 
  out[20] += mass*volFact*(1.4142135623730951*M0y[0]*wvpar+0.408248290463863*M0y[2]*dvpar); 
  out[21] += mass*volFact*(1.4142135623730951*M0y[1]*wvpar+0.408248290463863*M0y[3]*dvpar); 
  out[22] += mass*volFact*(1.4142135623730951*M0z[0]*wvpar+0.408248290463863*M0z[2]*dvpar); 
  out[23] += mass*volFact*(1.4142135623730951*M0z[1]*wvpar+0.408248290463863*M0z[3]*dvpar); 
  out[24] += mass*volFact*(1.4142135623730951*M0x[0]*wvpar_sq+0.8164965809277261*M0x[2]*dvpar*wvpar+0.10540925533894596*M0x[4]*dvpar_sq+0.11785113019775789*M0x[0]*dvpar_sq); 
  out[25] += mass*volFact*(1.4142135623730951*M0x[1]*wvpar_sq+0.8164965809277261*M0x[3]*dvpar*wvpar+0.105409255338946*M0x[5]*dvpar_sq+0.11785113019775789*M0x[1]*dvpar_sq); 
  out[26] += mass*volFact*(1.4142135623730951*M0y[0]*wvpar_sq+0.8164965809277261*M0y[2]*dvpar*wvpar+0.10540925533894596*M0y[4]*dvpar_sq+0.11785113019775789*M0y[0]*dvpar_sq); 
  out[27] += mass*volFact*(1.4142135623730951*M0y[1]*wvpar_sq+0.8164965809277261*M0y[3]*dvpar*wvpar+0.105409255338946*M0y[5]*dvpar_sq+0.11785113019775789*M0y[1]*dvpar_sq); 
  out[28] += mass*volFact*(1.4142135623730951*M0z[0]*wvpar_sq+0.8164965809277261*M0z[2]*dvpar*wvpar+0.10540925533894596*M0z[4]*dvpar_sq+0.11785113019775789*M0z[0]*dvpar_sq); 
  out[29] += mass*volFact*(1.4142135623730951*M0z[1]*wvpar_sq+0.8164965809277261*M0z[3]*dvpar*wvpar+0.105409255338946*M0z[5]*dvpar_sq+0.11785113019775789*M0z[1]*dvpar_sq); 
  out[30] += 1.4142135623730951*Nx[0]*mass*volFact; 
  out[31] += 1.4142135623730951*Nx[1]*mass*volFact; 
  out[32] += 1.4142135623730951*Ny[0]*mass*volFact; 
  out[33] += 1.4142135623730951*Ny[1]*mass*volFact; 
  out[34] += 1.4142135623730951*Nz[0]*mass*volFact; 
  out[35] += 1.4142135623730951*Nz[1]*mass*volFact; 
} 
GKYL_CU_DH void mom_pkpm_diag_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2.0; 
  const double wvpar = w[1], dvpar = dxv[1]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 
  const double wvpar_cu = wvpar*wvpar*wvpar, dvpar_cu = dvpar*dvpar*dvpar; 
  const double wvpar_qu = wvpar*wvpar*wvpar*wvpar, dvpar_qu = dvpar*dvpar*dvpar*dvpar; 

  const double *F_0 = &f[0]; 
  const double *G_1 = &f[6]; 
  out[0] += 1.4142135623730951*F_0[0]*mass*volFact; 
  out[1] += 1.4142135623730951*F_0[1]*mass*volFact; 
  out[2] += mass*volFact*(1.4142135623730951*F_0[0]*wvpar+0.408248290463863*F_0[2]*dvpar); 
  out[3] += mass*volFact*(1.4142135623730951*F_0[1]*wvpar+0.408248290463863*F_0[3]*dvpar); 
  out[4] += mass*volFact*(1.4142135623730951*F_0[0]*wvpar_sq+0.8164965809277261*F_0[2]*dvpar*wvpar+0.10540925533894596*F_0[4]*dvpar_sq+0.11785113019775789*F_0[0]*dvpar_sq); 
  out[5] += mass*volFact*(1.4142135623730951*F_0[1]*wvpar_sq+0.8164965809277261*F_0[3]*dvpar*wvpar+0.105409255338946*F_0[5]*dvpar_sq+0.11785113019775789*F_0[1]*dvpar_sq); 
  out[6] += 1.4142135623730951*G_1[0]*mass*volFact; 
  out[7] += 1.4142135623730951*G_1[1]*mass*volFact; 
  out[8] += mass*volFact*(1.224744871391589*F_0[2]*dvpar*wvpar_sq+1.4142135623730951*F_0[0]*wvpar_cu+0.3162277660168379*F_0[4]*dvpar_sq*wvpar+0.3535533905932737*F_0[0]*dvpar_sq*wvpar+0.06123724356957942*F_0[2]*dvpar_cu); 
  out[9] += mass*volFact*(1.224744871391589*F_0[3]*dvpar*wvpar_sq+1.4142135623730951*F_0[1]*wvpar_cu+0.31622776601683794*F_0[5]*dvpar_sq*wvpar+0.3535533905932737*F_0[1]*dvpar_sq*wvpar+0.06123724356957942*F_0[3]*dvpar_cu); 
  out[10] += mass*volFact*(1.4142135623730951*G_1[0]*wvpar+0.408248290463863*G_1[2]*dvpar); 
  out[11] += mass*volFact*(1.4142135623730951*G_1[1]*wvpar+0.408248290463863*G_1[3]*dvpar); 
  out[12] += mass*volFact*(0.6324555320336759*F_0[4]*dvpar_sq*wvpar_sq+0.7071067811865475*F_0[0]*dvpar_sq*wvpar_sq+1.4142135623730951*F_0[0]*wvpar_qu+1.6329931618554527*F_0[2]*dvpar*wvpar_cu+0.2449489742783178*F_0[2]*dvpar_cu*wvpar+0.022587697572631273*F_0[4]*dvpar_qu+0.01767766952966368*F_0[0]*dvpar_qu); 
  out[13] += mass*volFact*(0.632455532033676*F_0[5]*dvpar_sq*wvpar_sq+0.7071067811865475*F_0[1]*dvpar_sq*wvpar_sq+1.4142135623730951*F_0[1]*wvpar_qu+1.6329931618554527*F_0[3]*dvpar*wvpar_cu+0.2449489742783178*F_0[3]*dvpar_cu*wvpar+0.02258769757263128*F_0[5]*dvpar_qu+0.01767766952966368*F_0[1]*dvpar_qu); 
  out[14] += mass*volFact*(1.4142135623730951*G_1[0]*wvpar_sq+0.8164965809277261*G_1[2]*dvpar*wvpar+0.10540925533894596*G_1[4]*dvpar_sq+0.11785113019775789*G_1[0]*dvpar_sq); 
  out[15] += mass*volFact*(1.4142135623730951*G_1[1]*wvpar_sq+0.8164965809277261*G_1[3]*dvpar*wvpar+0.105409255338946*G_1[5]*dvpar_sq+0.11785113019775789*G_1[1]*dvpar_sq); 
} 
