#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH int pkpm_vars_prim_2x_ser_p1(const double *fluid, double* GKYL_RESTRICT u, double* GKYL_RESTRICT xi) 
{ 
  // fluid: Input fluid state [rho, rho ux, rho uy, rho uz, E]. 
  // u:     Output flow velocity [ux, uy, uz]: the nodal quotient M/rho, the weak division at p=1. 
  // xi:    Output [xi, xi^2, p], xi^2 = (2 E/rho - |u|^2)/3, p = rho xi^2, at the 4 Gauss-Legendre nodes, 
  //        interpolated to the p=1 basis. 
  // Returns the number of nodes where rho <= 0 (all outputs zero there) or xi^2 <= 0 ([xi, xi^2, p] zero there). 
 
  const double *rho = &fluid[0]; 
  const double *Mx = &fluid[4]; 
  const double *My = &fluid[8]; 
  const double *Mz = &fluid[12]; 
  const double *E = &fluid[16]; 
  double *ux = &u[0]; 
  double *uy = &u[4]; 
  double *uz = &u[8]; 
  double *xi_1 = &xi[0]; 
  double *xi_sq = &xi[4]; 
  double *p_sc = &xi[8]; 
 
  double u_n[3][4] = {{0.0}}; 
  double xi_n[4] = {0.0}; 
  double xisq_n[4] = {0.0}; 
  double p_n[4] = {0.0}; 
  int num_null = 0; 
  { 
    const double rho_q = 0.5*rho[3]-0.5*rho[2]-0.5*rho[1]+0.5*rho[0]; 
    const double E_q = 0.5*E[3]-0.5*E[2]-0.5*E[1]+0.5*E[0]; 
    const double M0_q = 0.5*Mx[3]-0.5*Mx[2]-0.5*Mx[1]+0.5*Mx[0]; 
    const double M1_q = 0.5*My[3]-0.5*My[2]-0.5*My[1]+0.5*My[0]; 
    const double M2_q = 0.5*Mz[3]-0.5*Mz[2]-0.5*Mz[1]+0.5*Mz[0]; 
    if (rho_q > 0.0) { 
      const double inv_rho = 1.0/rho_q; 
      const double xisq = (2.0*E_q - (M0_q*M0_q + M1_q*M1_q + M2_q*M2_q)*inv_rho)*inv_rho/3.0; 
      u_n[0][0] = M0_q*inv_rho; u_n[1][0] = M1_q*inv_rho; u_n[2][0] = M2_q*inv_rho; 
      if (xisq > 0.0) { 
        xisq_n[0] = xisq; xi_n[0] = sqrt(xisq); p_n[0] = rho_q*xisq; 
      } else { num_null += 1; } 
    } else { 
      num_null += 1; 
    } 
  } 
  { 
    const double rho_q = -(0.5*rho[3])+0.5*rho[2]-0.5*rho[1]+0.5*rho[0]; 
    const double E_q = -(0.5*E[3])+0.5*E[2]-0.5*E[1]+0.5*E[0]; 
    const double M0_q = -(0.5*Mx[3])+0.5*Mx[2]-0.5*Mx[1]+0.5*Mx[0]; 
    const double M1_q = -(0.5*My[3])+0.5*My[2]-0.5*My[1]+0.5*My[0]; 
    const double M2_q = -(0.5*Mz[3])+0.5*Mz[2]-0.5*Mz[1]+0.5*Mz[0]; 
    if (rho_q > 0.0) { 
      const double inv_rho = 1.0/rho_q; 
      const double xisq = (2.0*E_q - (M0_q*M0_q + M1_q*M1_q + M2_q*M2_q)*inv_rho)*inv_rho/3.0; 
      u_n[0][1] = M0_q*inv_rho; u_n[1][1] = M1_q*inv_rho; u_n[2][1] = M2_q*inv_rho; 
      if (xisq > 0.0) { 
        xisq_n[1] = xisq; xi_n[1] = sqrt(xisq); p_n[1] = rho_q*xisq; 
      } else { num_null += 1; } 
    } else { 
      num_null += 1; 
    } 
  } 
  { 
    const double rho_q = -(0.5*rho[3])-0.5*rho[2]+0.5*rho[1]+0.5*rho[0]; 
    const double E_q = -(0.5*E[3])-0.5*E[2]+0.5*E[1]+0.5*E[0]; 
    const double M0_q = -(0.5*Mx[3])-0.5*Mx[2]+0.5*Mx[1]+0.5*Mx[0]; 
    const double M1_q = -(0.5*My[3])-0.5*My[2]+0.5*My[1]+0.5*My[0]; 
    const double M2_q = -(0.5*Mz[3])-0.5*Mz[2]+0.5*Mz[1]+0.5*Mz[0]; 
    if (rho_q > 0.0) { 
      const double inv_rho = 1.0/rho_q; 
      const double xisq = (2.0*E_q - (M0_q*M0_q + M1_q*M1_q + M2_q*M2_q)*inv_rho)*inv_rho/3.0; 
      u_n[0][2] = M0_q*inv_rho; u_n[1][2] = M1_q*inv_rho; u_n[2][2] = M2_q*inv_rho; 
      if (xisq > 0.0) { 
        xisq_n[2] = xisq; xi_n[2] = sqrt(xisq); p_n[2] = rho_q*xisq; 
      } else { num_null += 1; } 
    } else { 
      num_null += 1; 
    } 
  } 
  { 
    const double rho_q = 0.5*rho[3]+0.5*rho[2]+0.5*rho[1]+0.5*rho[0]; 
    const double E_q = 0.5*E[3]+0.5*E[2]+0.5*E[1]+0.5*E[0]; 
    const double M0_q = 0.5*Mx[3]+0.5*Mx[2]+0.5*Mx[1]+0.5*Mx[0]; 
    const double M1_q = 0.5*My[3]+0.5*My[2]+0.5*My[1]+0.5*My[0]; 
    const double M2_q = 0.5*Mz[3]+0.5*Mz[2]+0.5*Mz[1]+0.5*Mz[0]; 
    if (rho_q > 0.0) { 
      const double inv_rho = 1.0/rho_q; 
      const double xisq = (2.0*E_q - (M0_q*M0_q + M1_q*M1_q + M2_q*M2_q)*inv_rho)*inv_rho/3.0; 
      u_n[0][3] = M0_q*inv_rho; u_n[1][3] = M1_q*inv_rho; u_n[2][3] = M2_q*inv_rho; 
      if (xisq > 0.0) { 
        xisq_n[3] = xisq; xi_n[3] = sqrt(xisq); p_n[3] = rho_q*xisq; 
      } else { num_null += 1; } 
    } else { 
      num_null += 1; 
    } 
  } 
 
  ux[0] = 0.5*u_n[0][3]+0.5*u_n[0][2]+0.5*u_n[0][1]+0.5*u_n[0][0]; 
  ux[1] = 0.5*u_n[0][3]+0.5*u_n[0][2]-0.5*u_n[0][1]-0.5*u_n[0][0]; 
  ux[2] = 0.5*u_n[0][3]-0.5*u_n[0][2]+0.5*u_n[0][1]-0.5*u_n[0][0]; 
  ux[3] = 0.5*u_n[0][3]-0.5*u_n[0][2]-0.5*u_n[0][1]+0.5*u_n[0][0]; 
  uy[0] = 0.5*u_n[1][3]+0.5*u_n[1][2]+0.5*u_n[1][1]+0.5*u_n[1][0]; 
  uy[1] = 0.5*u_n[1][3]+0.5*u_n[1][2]-0.5*u_n[1][1]-0.5*u_n[1][0]; 
  uy[2] = 0.5*u_n[1][3]-0.5*u_n[1][2]+0.5*u_n[1][1]-0.5*u_n[1][0]; 
  uy[3] = 0.5*u_n[1][3]-0.5*u_n[1][2]-0.5*u_n[1][1]+0.5*u_n[1][0]; 
  uz[0] = 0.5*u_n[2][3]+0.5*u_n[2][2]+0.5*u_n[2][1]+0.5*u_n[2][0]; 
  uz[1] = 0.5*u_n[2][3]+0.5*u_n[2][2]-0.5*u_n[2][1]-0.5*u_n[2][0]; 
  uz[2] = 0.5*u_n[2][3]-0.5*u_n[2][2]+0.5*u_n[2][1]-0.5*u_n[2][0]; 
  uz[3] = 0.5*u_n[2][3]-0.5*u_n[2][2]-0.5*u_n[2][1]+0.5*u_n[2][0]; 
  xi_1[0] = 0.5*xi_n[3]+0.5*xi_n[2]+0.5*xi_n[1]+0.5*xi_n[0]; 
  xi_1[1] = 0.5*xi_n[3]+0.5*xi_n[2]-0.5*xi_n[1]-0.5*xi_n[0]; 
  xi_1[2] = 0.5*xi_n[3]-0.5*xi_n[2]+0.5*xi_n[1]-0.5*xi_n[0]; 
  xi_1[3] = 0.5*xi_n[3]-0.5*xi_n[2]-0.5*xi_n[1]+0.5*xi_n[0]; 
  xi_sq[0] = 0.5*xisq_n[3]+0.5*xisq_n[2]+0.5*xisq_n[1]+0.5*xisq_n[0]; 
  xi_sq[1] = 0.5*xisq_n[3]+0.5*xisq_n[2]-0.5*xisq_n[1]-0.5*xisq_n[0]; 
  xi_sq[2] = 0.5*xisq_n[3]-0.5*xisq_n[2]+0.5*xisq_n[1]-0.5*xisq_n[0]; 
  xi_sq[3] = 0.5*xisq_n[3]-0.5*xisq_n[2]-0.5*xisq_n[1]+0.5*xisq_n[0]; 
  p_sc[0] = 0.5*p_n[3]+0.5*p_n[2]+0.5*p_n[1]+0.5*p_n[0]; 
  p_sc[1] = 0.5*p_n[3]+0.5*p_n[2]-0.5*p_n[1]-0.5*p_n[0]; 
  p_sc[2] = 0.5*p_n[3]-0.5*p_n[2]+0.5*p_n[1]-0.5*p_n[0]; 
  p_sc[3] = 0.5*p_n[3]-0.5*p_n[2]-0.5*p_n[1]+0.5*p_n[0]; 
 
  return num_null; 
} 
