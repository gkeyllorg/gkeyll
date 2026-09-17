#include <gkyl_basis_gkhyb_vel_kernels.h>
GKYL_CU_DH
void
node_coords_1d_gkhyb_vel_p1(double *node_coords)
{
  node_coords[0] = -1;
  node_coords[1] = 0;
  node_coords[2] = 1;
}

GKYL_CU_DH
void
node_coords_2d_gkhyb_vel_p1(double *node_coords)
{
  node_coords[0] = -1;
  node_coords[1] = -1;
  node_coords[2] = -1;
  node_coords[3] = 1;
  node_coords[4] = 0;
  node_coords[5] = -1;
  node_coords[6] = 0;
  node_coords[7] = 1;
  node_coords[8] = 1;
  node_coords[9] = -1;
  node_coords[10] = 1;
  node_coords[11] = 1;
}

