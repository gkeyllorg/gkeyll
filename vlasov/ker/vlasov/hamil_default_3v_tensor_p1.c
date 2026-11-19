#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_default_3v_tensor_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // vmap: Velocity-space map (C^1 cubic for tensor bases, C^0 linear for Serendipity bases, stored in the same vdim*4 layout).
  // hamil: Particle Hamiltonian.
  // hamil_inv: Inverse particle Hamiltonian. Only meaningful for relativistic simulations; zeroed for the non-relativistic Hamiltonian (1/H is singular at v = 0 and unused).
 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  const double *vmap_vz = &vmap[8]; 
  double hamil_vx_nodal[3] = {0.0};
  hamil_vx_nodal[0] = 0.5*pow(-(1.8708286933869707*vmap_vx[3])+1.5811388300841895*vmap_vx[2]-1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0);
  hamil_vx_nodal[1] = 0.5*pow(0.7071067811865475*vmap_vx[0]-0.7905694150420947*vmap_vx[2], 2.0);
  hamil_vx_nodal[2] = 0.5*pow(1.8708286933869707*vmap_vx[3]+1.5811388300841895*vmap_vx[2]+1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0);
  double hamil_vy_nodal[3] = {0.0};
  hamil_vy_nodal[0] = 0.5*pow(-(1.8708286933869707*vmap_vy[3])+1.5811388300841895*vmap_vy[2]-1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 2.0);
  hamil_vy_nodal[1] = 0.5*pow(0.7071067811865475*vmap_vy[0]-0.7905694150420947*vmap_vy[2], 2.0);
  hamil_vy_nodal[2] = 0.5*pow(1.8708286933869707*vmap_vy[3]+1.5811388300841895*vmap_vy[2]+1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 2.0);
  double hamil_vz_nodal[3] = {0.0};
  hamil_vz_nodal[0] = 0.5*pow(-(1.8708286933869707*vmap_vz[3])+1.5811388300841895*vmap_vz[2]-1.224744871391589*vmap_vz[1]+0.7071067811865475*vmap_vz[0], 2.0);
  hamil_vz_nodal[1] = 0.5*pow(0.7071067811865475*vmap_vz[0]-0.7905694150420947*vmap_vz[2], 2.0);
  hamil_vz_nodal[2] = 0.5*pow(1.8708286933869707*vmap_vz[3]+1.5811388300841895*vmap_vz[2]+1.224744871391589*vmap_vz[1]+0.7071067811865475*vmap_vz[0], 2.0);

  hamil[0] = 0.4714045207910317*hamil_vz_nodal[2]+0.4714045207910317*hamil_vy_nodal[2]+0.4714045207910317*hamil_vx_nodal[2]+1.8856180831641274*hamil_vz_nodal[1]+1.8856180831641274*hamil_vy_nodal[1]+1.8856180831641274*hamil_vx_nodal[1]+0.4714045207910317*hamil_vz_nodal[0]+0.4714045207910317*hamil_vy_nodal[0]+0.4714045207910317*hamil_vx_nodal[0]; 
  hamil[1] = 0.8164965809277261*hamil_vx_nodal[2]-0.8164965809277261*hamil_vx_nodal[0]; 
  hamil[2] = 0.8164965809277261*hamil_vy_nodal[2]-0.8164965809277261*hamil_vy_nodal[0]; 
  hamil[3] = 0.8164965809277261*hamil_vz_nodal[2]-0.8164965809277261*hamil_vz_nodal[0]; 
  hamil[4] = 0.0; 
  hamil[5] = 0.0; 
  hamil[6] = 0.0; 
  hamil[7] = 0.421637021355784*hamil_vx_nodal[2]-0.8432740427115681*hamil_vx_nodal[1]+0.421637021355784*hamil_vx_nodal[0]; 
  hamil[8] = 0.421637021355784*hamil_vy_nodal[2]-0.8432740427115681*hamil_vy_nodal[1]+0.421637021355784*hamil_vy_nodal[0]; 
  hamil[9] = 0.421637021355784*hamil_vz_nodal[2]-0.8432740427115681*hamil_vz_nodal[1]+0.421637021355784*hamil_vz_nodal[0]; 
  hamil[10] = 0.0; 
  hamil[11] = 0.0; 
  hamil[12] = 0.0; 
  hamil[13] = 0.0; 
  hamil[14] = 0.0; 
  hamil[15] = 0.0; 
  hamil[16] = 0.0; 
  hamil[17] = 0.0; 
  hamil[18] = 0.0; 
  hamil[19] = 0.0; 
  hamil[20] = 0.0; 
  hamil[21] = 0.0; 
  hamil[22] = 0.0; 
  hamil[23] = 0.0; 
  hamil[24] = 0.0; 
  hamil[25] = 0.0; 
  hamil[26] = 0.0; 

  hamil_inv[0] = 0.0; 
  hamil_inv[1] = 0.0; 
  hamil_inv[2] = 0.0; 
  hamil_inv[3] = 0.0; 
  hamil_inv[4] = 0.0; 
  hamil_inv[5] = 0.0; 
  hamil_inv[6] = 0.0; 
  hamil_inv[7] = 0.0; 
  hamil_inv[8] = 0.0; 
  hamil_inv[9] = 0.0; 
  hamil_inv[10] = 0.0; 
  hamil_inv[11] = 0.0; 
  hamil_inv[12] = 0.0; 
  hamil_inv[13] = 0.0; 
  hamil_inv[14] = 0.0; 
  hamil_inv[15] = 0.0; 
  hamil_inv[16] = 0.0; 
  hamil_inv[17] = 0.0; 
  hamil_inv[18] = 0.0; 
  hamil_inv[19] = 0.0; 
  hamil_inv[20] = 0.0; 
  hamil_inv[21] = 0.0; 
  hamil_inv[22] = 0.0; 
  hamil_inv[23] = 0.0; 
  hamil_inv[24] = 0.0; 
  hamil_inv[25] = 0.0; 
  hamil_inv[26] = 0.0; 

} 
