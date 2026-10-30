#include <gkyl_dg_diffusion_gen_kernels.h>

GKYL_CU_DH double
dg_diffusion_gen_vol_3x_ser_p1(const double* w, const double* dx,
  const double* Dij, const double* qIn, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // Dij: Diffusion coefficient in the center cell
  // q: Input field in the left cell
  // out: Incremented output

  const double Jxx = 4/dx[0]/dx[0];
  const double Jxy = 4/dx[0]/dx[1];
  const double Jxz = 4/dx[0]/dx[2];
  const double Jyy = 4/dx[1]/dx[1];
  const double Jyz = 4/dx[1]/dx[2];
  const double Jzz = 4/dx[2]/dx[2];

  const double* Dxx = &Dij[0];
  const double* Dxy = &Dij[8];
  const double* Dxz = &Dij[16];
  const double* Dyy = &Dij[24];
  const double* Dyz = &Dij[32];
  const double* Dzz = &Dij[40];

  const double Dxx_av = 0.3535533905932737*Dxx[0];
  const double Dxy_av = 0.3535533905932737*Dxy[0];
  const double Dxz_av = 0.3535533905932737*Dxz[0];
  const double Dyy_av = 0.3535533905932737*Dyy[0];
  const double Dyz_av = 0.3535533905932737*Dyz[0];
  const double Dzz_av = 0.3535533905932737*Dzz[0];

  return 4.0*(Jxx*Dxx_av + Jxy*fabs(Dxy_av) + Jxz*fabs(Dxz_av) + Jyy*Dyy_av + Jyz*fabs(Dyz_av) + Jzz*Dzz_av);
}
