#include <gkyl_pkpm_vars_kernels.h> 
#include <math.h> 
GKYL_CU_DH int pkpm_fluid_positivity_2x_ser_p1(double eps, double* GKYL_RESTRICT fluid) 
{ 
  // eps:   Relative floor: the vertex density and internal energy are kept above eps times the cell values. 
  // fluid: In/out fluid state [rho, rho ux, rho uy, rho uz, E]; the slopes are scaled toward the cell averages. 
  // Returns -1 (inadmissible cell average: nothing done), or the sum of 1 (density slopes scaled) 
  // and 2 (all slopes scaled for the internal energy), 0 when untouched. 
  double *rho = &fluid[0]; 
  double *Mx = &fluid[4]; 
  double *My = &fluid[8]; 
  double *Mz = &fluid[12]; 
  double *E = &fluid[16]; 
 
  // Cell averages. 
  const double rho_avg = 0.5*rho[0]; 
  const double Mx_avg = 0.5*Mx[0]; 
  const double My_avg = 0.5*My[0]; 
  const double Mz_avg = 0.5*Mz[0]; 
  const double E_avg = 0.5*E[0]; 
  if (!(rho_avg > 0.0)) return -1; 
  const double e_avg = E_avg - 0.5*(Mx_avg*Mx_avg + My_avg*My_avg + Mz_avg*Mz_avg)/rho_avg; 
  if (!(e_avg > 0.0)) return -1; 
 
  // Vertex values. 
  double rho_v[4]; 
  rho_v[0] = (0.5)*rho[0] + (-0.8660254037844386)*rho[1] + (-0.8660254037844386)*rho[2] + (1.5)*rho[3]; 
  rho_v[1] = (0.5)*rho[0] + (0.8660254037844386)*rho[1] + (-0.8660254037844386)*rho[2] + (-1.5)*rho[3]; 
  rho_v[2] = (0.5)*rho[0] + (-0.8660254037844386)*rho[1] + (0.8660254037844386)*rho[2] + (-1.5)*rho[3]; 
  rho_v[3] = (0.5)*rho[0] + (0.8660254037844386)*rho[1] + (0.8660254037844386)*rho[2] + (1.5)*rho[3]; 
  double Mx_v[4]; 
  Mx_v[0] = (0.5)*Mx[0] + (-0.8660254037844386)*Mx[1] + (-0.8660254037844386)*Mx[2] + (1.5)*Mx[3]; 
  Mx_v[1] = (0.5)*Mx[0] + (0.8660254037844386)*Mx[1] + (-0.8660254037844386)*Mx[2] + (-1.5)*Mx[3]; 
  Mx_v[2] = (0.5)*Mx[0] + (-0.8660254037844386)*Mx[1] + (0.8660254037844386)*Mx[2] + (-1.5)*Mx[3]; 
  Mx_v[3] = (0.5)*Mx[0] + (0.8660254037844386)*Mx[1] + (0.8660254037844386)*Mx[2] + (1.5)*Mx[3]; 
  double My_v[4]; 
  My_v[0] = (0.5)*My[0] + (-0.8660254037844386)*My[1] + (-0.8660254037844386)*My[2] + (1.5)*My[3]; 
  My_v[1] = (0.5)*My[0] + (0.8660254037844386)*My[1] + (-0.8660254037844386)*My[2] + (-1.5)*My[3]; 
  My_v[2] = (0.5)*My[0] + (-0.8660254037844386)*My[1] + (0.8660254037844386)*My[2] + (-1.5)*My[3]; 
  My_v[3] = (0.5)*My[0] + (0.8660254037844386)*My[1] + (0.8660254037844386)*My[2] + (1.5)*My[3]; 
  double Mz_v[4]; 
  Mz_v[0] = (0.5)*Mz[0] + (-0.8660254037844386)*Mz[1] + (-0.8660254037844386)*Mz[2] + (1.5)*Mz[3]; 
  Mz_v[1] = (0.5)*Mz[0] + (0.8660254037844386)*Mz[1] + (-0.8660254037844386)*Mz[2] + (-1.5)*Mz[3]; 
  Mz_v[2] = (0.5)*Mz[0] + (-0.8660254037844386)*Mz[1] + (0.8660254037844386)*Mz[2] + (-1.5)*Mz[3]; 
  Mz_v[3] = (0.5)*Mz[0] + (0.8660254037844386)*Mz[1] + (0.8660254037844386)*Mz[2] + (1.5)*Mz[3]; 
  double E_v[4]; 
  E_v[0] = (0.5)*E[0] + (-0.8660254037844386)*E[1] + (-0.8660254037844386)*E[2] + (1.5)*E[3]; 
  E_v[1] = (0.5)*E[0] + (0.8660254037844386)*E[1] + (-0.8660254037844386)*E[2] + (-1.5)*E[3]; 
  E_v[2] = (0.5)*E[0] + (-0.8660254037844386)*E[1] + (0.8660254037844386)*E[2] + (-1.5)*E[3]; 
  E_v[3] = (0.5)*E[0] + (0.8660254037844386)*E[1] + (0.8660254037844386)*E[2] + (1.5)*E[3]; 
 
  int flag = 0; 
  // (1) Density: scale the density slopes so that rho_v >= eps rho_avg. 
  double theta = 1.0; 
  const double rho_min = eps*rho_avg; 
  for (int k = 0; k < 4; ++k) { 
    if (rho_v[k] < rho_min) theta = fmin(theta, (rho_avg - rho_min)/(rho_avg - rho_v[k])); 
  } 
  if (theta < 1.0) { 
    for (int i = 1; i < 4; ++i) rho[i] *= theta; 
    for (int k = 0; k < 4; ++k) rho_v[k] = rho_avg + theta*(rho_v[k] - rho_avg); 
    flag |= 1; 
  } 
 
  // (2) Internal energy: scale all slopes so that e(q_v) >= eps e(q_avg). 
  const double e_min = eps*e_avg; 
  const double c = rho_avg*(e_avg - e_min); 
  theta = 1.0; 
  for (int k = 0; k < 4; ++k) { 
    const double drho = rho_v[k] - rho_avg, dMx = Mx_v[k] - Mx_avg, dMy = My_v[k] - My_avg, dMz = Mz_v[k] - Mz_avg, dE = E_v[k] - E_avg; 
    // g(t) = rho(t) E(t) - |M(t)|^2/2 - e_min rho(t) = a t^2 + b t + c along q(t) = q_avg + t (q_v - q_avg). 
    const double a = drho*dE - 0.5*(dMx*dMx + dMy*dMy + dMz*dMz); 
    const double b = rho_avg*dE + drho*E_avg - (Mx_avg*dMx + My_avg*dMy + Mz_avg*dMz) - e_min*drho; 
    if (a + b + c < 0.0) { 
      // g(0) = c > 0 > g(1): exactly one root in (0, 1). 
      double t; 
      if (fabs(a) <= 1.0e-14*(fabs(b) + fabs(c))) { 
        t = -c/b; 
      } else { 
        const double sq = sqrt(fmax(b*b - 4.0*a*c, 0.0)); 
        const double qq = -0.5*(b + (b >= 0.0 ? sq : -sq)); 
        const double t1 = qq/a, t2 = qq != 0.0 ? c/qq : t1; 
        t = (t1 > 0.0 && t1 < 1.0) ? t1 : t2; 
      } 
      theta = fmin(theta, fmax(fmin(t, 1.0), 0.0)); 
    } 
  } 
  if (theta < 1.0) { 
    for (int i = 1; i < 4; ++i) rho[i] *= theta; 
    for (int i = 1; i < 4; ++i) Mx[i] *= theta; 
    for (int i = 1; i < 4; ++i) My[i] *= theta; 
    for (int i = 1; i < 4; ++i) Mz[i] *= theta; 
    for (int i = 1; i < 4; ++i) E[i] *= theta; 
    flag |= 2; 
  } 
  return flag; 
} 
