// Sat Jun 27 17:10:24 2026
#include <gkyl_basis_gkhyb_vel_kernels.h>
GKYL_CU_DH
void
eval_1d_gkhyb_vel_p1(const double *z, double *b )
{
  const double z0 = z[0];
  b[0] = 7.0710678118654757e-01;
  b[1] = 1.2247448713915889e+00*z0;
  b[2] =  2.3717082451262845e+00*(z0*z0)-7.9056941504209488e-01;
}

GKYL_CU_DH
double
eval_expand_1d_gkhyb_vel_p1(const double *z, const double *f )
{
  const double z0 = z[0];
  return  -7.9056941504209488e-01*f[2]+1.2247448713915889e+00*z0*f[1]+2.3717082451262845e+00*(z0*z0)*f[2]+7.0710678118654757e-01*f[0];
}

GKYL_CU_DH
double
eval_grad_expand_1d_gkhyb_vel_p1(int dir, const double *z, const double *f )
{
  const double z0 = z[0];
  if (dir == 0)
    return  1.2247448713915889e+00*f[1]+4.7434164902525691e+00*f[2]*z0;

  return 0.0; // can't happen, suppresses warning

}

GKYL_CU_DH
void
eval_2d_gkhyb_vel_p1(const double *z, double *b )
{
  const double z0 = z[0];
  const double z1 = z[1];
  b[0] = 5.0000000000000000e-01;
  b[1] = 8.6602540378443860e-01*z0;
  b[2] = 8.6602540378443860e-01*z1;
  b[3] = 1.5000000000000000e+00*z1*z0;
  b[4] =  1.6770509831248424e+00*(z0*z0)-5.5901699437494745e-01;
  b[5] =  -9.6824583655185426e-01*z1+2.9047375096555625e+00*z1*(z0*z0);
}

GKYL_CU_DH
double
eval_expand_2d_gkhyb_vel_p1(const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
  return  8.6602540378443860e-01*z0*f[1]+-9.6824583655185426e-01*z1*f[5]+-5.5901699437494745e-01*f[4]+5.0000000000000000e-01*f[0]+8.6602540378443860e-01*z1*f[2]+1.6770509831248424e+00*f[4]*(z0*z0)+1.5000000000000000e+00*z1*z0*f[3]+2.9047375096555625e+00*z1*f[5]*(z0*z0);
}

GKYL_CU_DH
double
eval_grad_expand_2d_gkhyb_vel_p1(int dir, const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
  if (dir == 0)
    return  1.5000000000000000e+00*z1*f[3]+8.6602540378443860e-01*f[1]+3.3541019662496847e+00*z0*f[4]+5.8094750193111251e+00*z1*z0*f[5];

  if (dir == 1)
    return  8.6602540378443860e-01*f[2]+2.9047375096555625e+00*(z0*z0)*f[5]+1.5000000000000000e+00*f[3]*z0+-9.6824583655185426e-01*f[5];

  return 0.0; // can't happen, suppresses warning

}

