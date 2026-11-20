#include <gkyl_pkpm_fluid_kernels.h> 
GKYL_CU_DH double pkpm_fluid_vol_1x_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
  const double *vlasov_pkpm_moms, const double *bvar, const double *u, const double *fluid,
  double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:          Cell-center coordinates (computational). 
  // dxv[NDIM]:        Cell spacing (computational). 
  // jacob_pos:        Position-map Jacobian per configuration direction (cdim blocks of 2 equal 
  //                   entries); the divergence in direction d carries 1/J_d. 
  // vlasov_pkpm_moms: Jacobian-weighted kinetic coupling moments J [p_par, p_perp, q_par, q_perp, 
  //                   s (x, y, z), r (x, y, z)]: s the parallel-perpendicular stress, r the 
  //                   perpendicular heat flux. 
  // bvar:             Magnetic field unit vector [bx, by, bz]. 
  // u:                Flow velocity [ux, uy, uz]. 
  // fluid:            Evolved Jacobian-weighted fluid state J [rho, rho ux, rho uy, rho uz, E]. 
  // out:              Incremented increment of the evolved J [rho, rho ux, rho uy, rho uz, E]. 
  const double dx10 = 2.0/(dxv[0]*jacob_pos[0]); 
  const double *ppar = &vlasov_pkpm_moms[0]; 
  const double *pperp = &vlasov_pkpm_moms[2]; 
  const double *qpar = &vlasov_pkpm_moms[4]; 
  const double *qperp = &vlasov_pkpm_moms[6]; 
  const double *sx = &vlasov_pkpm_moms[8]; 
  const double *rx = &vlasov_pkpm_moms[14]; 
  const double *sy = &vlasov_pkpm_moms[10]; 
  const double *ry = &vlasov_pkpm_moms[16]; 
  const double *sz = &vlasov_pkpm_moms[12]; 
  const double *rz = &vlasov_pkpm_moms[18]; 
  const double *E = &fluid[8]; 
  const double *bx = &bvar[0]; 
  const double *ux = &u[0]; 
  const double *Mx = &fluid[2]; 
  double *out1 = &out[2]; 
  const double *by = &bvar[2]; 
  const double *uy = &u[2]; 
  const double *My = &fluid[4]; 
  double *out2 = &out[4]; 
  const double *bz = &bvar[4]; 
  const double *uz = &u[4]; 
  const double *Mz = &fluid[6]; 
  double *out3 = &out[6]; 
  double *out0 = &out[0]; 
  double *out4 = &out[8]; 
 
  out0[1] += 1.7320508075688772*Mx[0]*dx10; 
  out1[1] += (-(0.8660254037844386*pperp[0]*pow(bx[1],2))+0.8660254037844386*ppar[0]*pow(bx[1],2)-0.8660254037844386*pperp[0]*pow(bx[0],2)+0.8660254037844386*ppar[0]*pow(bx[0],2)+1.224744871391589*Mx[1]*ux[1]+2.4494897427831783*bx[1]*sx[1]-1.7320508075688772*bx[0]*bx[1]*pperp[1]+1.7320508075688772*bx[0]*bx[1]*ppar[1]+1.224744871391589*Mx[0]*ux[0]+2.4494897427831783*bx[0]*sx[0]+1.7320508075688772*pperp[0])*dx10; 
  out2[1] += (1.224744871391589*My[1]*ux[1]+1.224744871391589*bx[1]*sy[1]+1.224744871391589*by[1]*sx[1]-0.8660254037844386*bx[0]*by[1]*pperp[1]-0.8660254037844386*by[0]*bx[1]*pperp[1]+0.8660254037844386*bx[0]*by[1]*ppar[1]+0.8660254037844386*by[0]*bx[1]*ppar[1]-0.8660254037844386*pperp[0]*bx[1]*by[1]+0.8660254037844386*ppar[0]*bx[1]*by[1]+1.224744871391589*My[0]*ux[0]+1.224744871391589*bx[0]*sy[0]+1.224744871391589*by[0]*sx[0]-0.8660254037844386*bx[0]*by[0]*pperp[0]+0.8660254037844386*bx[0]*by[0]*ppar[0])*dx10; 
  out3[1] += (1.224744871391589*Mz[1]*ux[1]+1.224744871391589*bx[1]*sz[1]+1.224744871391589*bz[1]*sx[1]-0.8660254037844386*bx[0]*bz[1]*pperp[1]-0.8660254037844386*bz[0]*bx[1]*pperp[1]+0.8660254037844386*bx[0]*bz[1]*ppar[1]+0.8660254037844386*bz[0]*bx[1]*ppar[1]-0.8660254037844386*pperp[0]*bx[1]*bz[1]+0.8660254037844386*ppar[0]*bx[1]*bz[1]+1.224744871391589*Mz[0]*ux[0]+1.224744871391589*bx[0]*sz[0]+1.224744871391589*bz[0]*sx[0]-0.8660254037844386*bx[0]*bz[0]*pperp[0]+0.8660254037844386*bx[0]*bz[0]*ppar[0])*dx10; 
  out4[1] += (-(1.1022703842524295*pperp[1]*ux[1]*pow(bx[1],2))+1.1022703842524295*ppar[1]*ux[1]*pow(bx[1],2)-0.6123724356957944*pperp[0]*ux[0]*pow(bx[1],2)+0.6123724356957944*ppar[0]*ux[0]*pow(bx[1],2)-0.6123724356957944*pperp[1]*ux[1]*pow(bx[0],2)+0.6123724356957944*ppar[1]*ux[1]*pow(bx[0],2)-0.6123724356957944*pperp[0]*ux[0]*pow(bx[0],2)+0.6123724356957944*ppar[0]*ux[0]*pow(bx[0],2)+0.8660254037844386*bx[0]*sz[1]*uz[1]+0.8660254037844386*bz[0]*sx[1]*uz[1]-1.1022703842524295*bx[1]*bz[1]*pperp[1]*uz[1]-0.6123724356957944*bx[0]*bz[0]*pperp[1]*uz[1]+1.1022703842524295*bx[1]*bz[1]*ppar[1]*uz[1]+0.6123724356957944*bx[0]*bz[0]*ppar[1]*uz[1]+0.8660254037844386*sx[0]*bz[1]*uz[1]-0.6123724356957944*bx[0]*pperp[0]*bz[1]*uz[1]+0.6123724356957944*bx[0]*ppar[0]*bz[1]*uz[1]+0.8660254037844386*sz[0]*bx[1]*uz[1]-0.6123724356957944*bz[0]*pperp[0]*bx[1]*uz[1]+0.6123724356957944*bz[0]*ppar[0]*bx[1]*uz[1]+0.8660254037844386*bx[0]*sy[1]*uy[1]+0.8660254037844386*by[0]*sx[1]*uy[1]-1.1022703842524295*bx[1]*by[1]*pperp[1]*uy[1]-0.6123724356957944*bx[0]*by[0]*pperp[1]*uy[1]+1.1022703842524295*bx[1]*by[1]*ppar[1]*uy[1]+0.6123724356957944*bx[0]*by[0]*ppar[1]*uy[1]+0.8660254037844386*sx[0]*by[1]*uy[1]-0.6123724356957944*bx[0]*pperp[0]*by[1]*uy[1]+0.6123724356957944*bx[0]*ppar[0]*by[1]*uy[1]+0.8660254037844386*sy[0]*bx[1]*uy[1]-0.6123724356957944*by[0]*pperp[0]*bx[1]*uy[1]+0.6123724356957944*by[0]*ppar[0]*bx[1]*uy[1]+1.7320508075688772*bx[0]*sx[1]*ux[1]+1.224744871391589*pperp[1]*ux[1]+1.7320508075688772*sx[0]*bx[1]*ux[1]-1.224744871391589*bx[0]*pperp[0]*bx[1]*ux[1]+1.224744871391589*bx[0]*ppar[0]*bx[1]*ux[1]+1.224744871391589*E[1]*ux[1]+0.8660254037844386*uz[0]*bx[1]*sz[1]+0.8660254037844386*uy[0]*bx[1]*sy[1]+0.8660254037844386*uz[0]*bz[1]*sx[1]+0.8660254037844386*uy[0]*by[1]*sx[1]+1.7320508075688772*ux[0]*bx[1]*sx[1]+1.224744871391589*bx[1]*qperp[1]+0.6123724356957944*bx[1]*qpar[1]-0.6123724356957944*bx[0]*uz[0]*bz[1]*pperp[1]-0.6123724356957944*bx[0]*uy[0]*by[1]*pperp[1]-0.6123724356957944*bz[0]*uz[0]*bx[1]*pperp[1]-0.6123724356957944*by[0]*uy[0]*bx[1]*pperp[1]-1.224744871391589*bx[0]*ux[0]*bx[1]*pperp[1]+0.6123724356957944*bx[0]*uz[0]*bz[1]*ppar[1]+0.6123724356957944*bx[0]*uy[0]*by[1]*ppar[1]+0.6123724356957944*bz[0]*uz[0]*bx[1]*ppar[1]+0.6123724356957944*by[0]*uy[0]*bx[1]*ppar[1]+1.224744871391589*bx[0]*ux[0]*bx[1]*ppar[1]-0.6123724356957944*pperp[0]*uz[0]*bx[1]*bz[1]+0.6123724356957944*ppar[0]*uz[0]*bx[1]*bz[1]-0.6123724356957944*pperp[0]*uy[0]*bx[1]*by[1]+0.6123724356957944*ppar[0]*uy[0]*bx[1]*by[1]+0.8660254037844386*bx[0]*sz[0]*uz[0]+0.8660254037844386*bz[0]*sx[0]*uz[0]-0.6123724356957944*bx[0]*bz[0]*pperp[0]*uz[0]+0.6123724356957944*bx[0]*bz[0]*ppar[0]*uz[0]+0.8660254037844386*bx[0]*sy[0]*uy[0]+0.8660254037844386*by[0]*sx[0]*uy[0]-0.6123724356957944*bx[0]*by[0]*pperp[0]*uy[0]+0.6123724356957944*bx[0]*by[0]*ppar[0]*uy[0]+1.7320508075688772*bx[0]*sx[0]*ux[0]+1.224744871391589*pperp[0]*ux[0]+1.224744871391589*E[0]*ux[0]+1.7320508075688772*rx[0]+1.224744871391589*bx[0]*qperp[0]+0.6123724356957944*bx[0]*qpar[0])*dx10; 
 
  double cflFreq_mid = 0.0; 
  cflFreq_mid += 1.5*dx10*fabs(0.7071067811865475*ux[0]); 
  return cflFreq_mid; 
} 
