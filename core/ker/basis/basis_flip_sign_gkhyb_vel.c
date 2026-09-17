// Sat Jun 27 17:10:24 2026
#include <gkyl_basis_gkhyb_vel_kernels.h>
GKYL_CU_DH
void
flip_odd_sign_1d_gkhyb_vel_p1(int dir, const double *f, double *fout )
{
  if (dir == 0) {
    fout[0] = 1*f[0];
    fout[1] = -1*f[1];
    fout[2] = 1*f[2];
  }
}

GKYL_CU_DH
void
flip_even_sign_1d_gkhyb_vel_p1(int dir, const double *f, double *fout )
{
  if (dir == 0) {
    fout[0] = -1*f[0];
    fout[1] = 1*f[1];
    fout[2] = -1*f[2];
  }
}

GKYL_CU_DH
void
flip_odd_sign_2d_gkhyb_vel_p1(int dir, const double *f, double *fout )
{
  if (dir == 0) {
    fout[0] = 1*f[0];
    fout[1] = -1*f[1];
    fout[2] = 1*f[2];
    fout[3] = -1*f[3];
    fout[4] = 1*f[4];
    fout[5] = 1*f[5];
  }
  if (dir == 1) {
    fout[0] = 1*f[0];
    fout[1] = 1*f[1];
    fout[2] = -1*f[2];
    fout[3] = -1*f[3];
    fout[4] = 1*f[4];
    fout[5] = -1*f[5];
  }
}

GKYL_CU_DH
void
flip_even_sign_2d_gkhyb_vel_p1(int dir, const double *f, double *fout )
{
  if (dir == 0) {
    fout[0] = -1*f[0];
    fout[1] = 1*f[1];
    fout[2] = -1*f[2];
    fout[3] = 1*f[3];
    fout[4] = -1*f[4];
    fout[5] = -1*f[5];
  }
  if (dir == 1) {
    fout[0] = -1*f[0];
    fout[1] = -1*f[1];
    fout[2] = 1*f[2];
    fout[3] = 1*f[3];
    fout[4] = -1*f[4];
    fout[5] = 1*f[5];
  }
}

