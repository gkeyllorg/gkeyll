#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_sr_2v_ser_p2(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // vmap: Velocity-space map (C^1 cubic for tensor bases, C^0 linear for Serendipity bases, stored in the same vdim*4 layout).
  // hamil: Particle Hamiltonian.
  // hamil_inv: Inverse particle Hamiltonian. Only meaningful for relativistic simulations; zeroed for the non-relativistic Hamiltonian (1/H is singular at v = 0 and unused).
 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  double hamil_nodal[8] = {0.0};
  double hamil_inv_nodal[8] = {0.0};
  hamil_nodal[0] = sqrt(1.0 + pow(0.7071067811865475*vmap_vx[0]-1.224744871391589*vmap_vx[1], 2.0) + pow(0.7071067811865475*vmap_vy[0]-1.224744871391589*vmap_vy[1], 2.0));
  hamil_inv_nodal[0] = 1.0/hamil_nodal[0];
  hamil_nodal[1] = sqrt(1.0 + pow(0.7071067811865475*vmap_vx[0], 2.0) + pow(0.7071067811865475*vmap_vy[0]-1.224744871391589*vmap_vy[1], 2.0));
  hamil_inv_nodal[1] = 1.0/hamil_nodal[1];
  hamil_nodal[2] = sqrt(1.0 + pow(1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0) + pow(0.7071067811865475*vmap_vy[0]-1.224744871391589*vmap_vy[1], 2.0));
  hamil_inv_nodal[2] = 1.0/hamil_nodal[2];
  hamil_nodal[3] = sqrt(1.0 + pow(0.7071067811865475*vmap_vx[0]-1.224744871391589*vmap_vx[1], 2.0) + pow(0.7071067811865475*vmap_vy[0], 2.0));
  hamil_inv_nodal[3] = 1.0/hamil_nodal[3];
  hamil_nodal[4] = sqrt(1.0 + pow(1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0) + pow(0.7071067811865475*vmap_vy[0], 2.0));
  hamil_inv_nodal[4] = 1.0/hamil_nodal[4];
  hamil_nodal[5] = sqrt(1.0 + pow(0.7071067811865475*vmap_vx[0]-1.224744871391589*vmap_vx[1], 2.0) + pow(1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 2.0));
  hamil_inv_nodal[5] = 1.0/hamil_nodal[5];
  hamil_nodal[6] = sqrt(1.0 + pow(0.7071067811865475*vmap_vx[0], 2.0) + pow(1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 2.0));
  hamil_inv_nodal[6] = 1.0/hamil_nodal[6];
  hamil_nodal[7] = sqrt(1.0 + pow(1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0) + pow(1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 2.0));
  hamil_inv_nodal[7] = 1.0/hamil_nodal[7];

  hamil[0] = -(0.16666666666666666*hamil_nodal[7])+0.6666666666666666*hamil_nodal[6]-0.16666666666666666*hamil_nodal[5]+0.6666666666666666*hamil_nodal[4]+0.6666666666666666*hamil_nodal[3]-0.16666666666666666*hamil_nodal[2]+0.6666666666666666*hamil_nodal[1]-0.16666666666666666*hamil_nodal[0]; 
  hamil[1] = 0.09622504486493764*hamil_nodal[7]-0.09622504486493764*hamil_nodal[5]+0.3849001794597506*hamil_nodal[4]-0.3849001794597506*hamil_nodal[3]+0.09622504486493764*hamil_nodal[2]-0.09622504486493764*hamil_nodal[0]; 
  hamil[2] = 0.09622504486493764*hamil_nodal[7]+0.3849001794597506*hamil_nodal[6]+0.09622504486493764*hamil_nodal[5]-0.09622504486493764*hamil_nodal[2]-0.3849001794597506*hamil_nodal[1]-0.09622504486493764*hamil_nodal[0]; 
  hamil[3] = 0.16666666666666666*hamil_nodal[7]-0.16666666666666666*hamil_nodal[5]-0.16666666666666666*hamil_nodal[2]+0.16666666666666666*hamil_nodal[0]; 
  hamil[4] = 0.14907119849998596*hamil_nodal[7]-0.2981423969999719*hamil_nodal[6]+0.14907119849998596*hamil_nodal[5]+0.14907119849998596*hamil_nodal[2]-0.2981423969999719*hamil_nodal[1]+0.14907119849998596*hamil_nodal[0]; 
  hamil[5] = 0.14907119849998596*hamil_nodal[7]+0.14907119849998596*hamil_nodal[5]-0.2981423969999719*hamil_nodal[4]-0.2981423969999719*hamil_nodal[3]+0.14907119849998596*hamil_nodal[2]+0.14907119849998596*hamil_nodal[0]; 
  hamil[6] = 0.08606629658238703*hamil_nodal[7]-0.17213259316477406*hamil_nodal[6]+0.08606629658238703*hamil_nodal[5]-0.08606629658238703*hamil_nodal[2]+0.17213259316477406*hamil_nodal[1]-0.08606629658238703*hamil_nodal[0]; 
  hamil[7] = 0.08606629658238703*hamil_nodal[7]-0.08606629658238703*hamil_nodal[5]-0.17213259316477406*hamil_nodal[4]+0.17213259316477406*hamil_nodal[3]+0.08606629658238703*hamil_nodal[2]-0.08606629658238703*hamil_nodal[0]; 

  hamil_inv[0] = -(0.16666666666666666*hamil_inv_nodal[7])+0.6666666666666666*hamil_inv_nodal[6]-0.16666666666666666*hamil_inv_nodal[5]+0.6666666666666666*hamil_inv_nodal[4]+0.6666666666666666*hamil_inv_nodal[3]-0.16666666666666666*hamil_inv_nodal[2]+0.6666666666666666*hamil_inv_nodal[1]-0.16666666666666666*hamil_inv_nodal[0]; 
  hamil_inv[1] = 0.09622504486493764*hamil_inv_nodal[7]-0.09622504486493764*hamil_inv_nodal[5]+0.3849001794597506*hamil_inv_nodal[4]-0.3849001794597506*hamil_inv_nodal[3]+0.09622504486493764*hamil_inv_nodal[2]-0.09622504486493764*hamil_inv_nodal[0]; 
  hamil_inv[2] = 0.09622504486493764*hamil_inv_nodal[7]+0.3849001794597506*hamil_inv_nodal[6]+0.09622504486493764*hamil_inv_nodal[5]-0.09622504486493764*hamil_inv_nodal[2]-0.3849001794597506*hamil_inv_nodal[1]-0.09622504486493764*hamil_inv_nodal[0]; 
  hamil_inv[3] = 0.16666666666666666*hamil_inv_nodal[7]-0.16666666666666666*hamil_inv_nodal[5]-0.16666666666666666*hamil_inv_nodal[2]+0.16666666666666666*hamil_inv_nodal[0]; 
  hamil_inv[4] = 0.14907119849998596*hamil_inv_nodal[7]-0.2981423969999719*hamil_inv_nodal[6]+0.14907119849998596*hamil_inv_nodal[5]+0.14907119849998596*hamil_inv_nodal[2]-0.2981423969999719*hamil_inv_nodal[1]+0.14907119849998596*hamil_inv_nodal[0]; 
  hamil_inv[5] = 0.14907119849998596*hamil_inv_nodal[7]+0.14907119849998596*hamil_inv_nodal[5]-0.2981423969999719*hamil_inv_nodal[4]-0.2981423969999719*hamil_inv_nodal[3]+0.14907119849998596*hamil_inv_nodal[2]+0.14907119849998596*hamil_inv_nodal[0]; 
  hamil_inv[6] = 0.08606629658238703*hamil_inv_nodal[7]-0.17213259316477406*hamil_inv_nodal[6]+0.08606629658238703*hamil_inv_nodal[5]-0.08606629658238703*hamil_inv_nodal[2]+0.17213259316477406*hamil_inv_nodal[1]-0.08606629658238703*hamil_inv_nodal[0]; 
  hamil_inv[7] = 0.08606629658238703*hamil_inv_nodal[7]-0.08606629658238703*hamil_inv_nodal[5]-0.17213259316477406*hamil_inv_nodal[4]+0.17213259316477406*hamil_inv_nodal[3]+0.08606629658238703*hamil_inv_nodal[2]-0.08606629658238703*hamil_inv_nodal[0]; 

} 
