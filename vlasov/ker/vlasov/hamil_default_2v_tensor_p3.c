#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_default_2v_tensor_p3(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // vmap: Velocity-space map (C^1 cubic for tensor bases, C^0 linear for Serendipity bases, stored in the same vdim*4 layout).
  // hamil: Particle Hamiltonian.
  // hamil_inv: Inverse particle Hamiltonian. Only meaningful for relativistic simulations; zeroed for the non-relativistic Hamiltonian (1/H is singular at v = 0 and unused).
 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  double hamil_vx_nodal[4] = {0.0};
  hamil_vx_nodal[0] = 0.5*pow(-(1.8708286933869707*vmap_vx[3])+1.5811388300841895*vmap_vx[2]-1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0);
  hamil_vx_nodal[1] = 0.5*pow(0.8366600265340755*vmap_vx[3]-0.3162277660168379*vmap_vx[2]-0.5477225575051661*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0);
  hamil_vx_nodal[2] = 0.5*pow(-(0.8366600265340755*vmap_vx[3])-0.3162277660168379*vmap_vx[2]+0.5477225575051661*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0);
  hamil_vx_nodal[3] = 0.5*pow(1.8708286933869707*vmap_vx[3]+1.5811388300841895*vmap_vx[2]+1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0);
  double hamil_vy_nodal[4] = {0.0};
  hamil_vy_nodal[0] = 0.5*pow(-(1.8708286933869707*vmap_vy[3])+1.5811388300841895*vmap_vy[2]-1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 2.0);
  hamil_vy_nodal[1] = 0.5*pow(0.8366600265340755*vmap_vy[3]-0.3162277660168379*vmap_vy[2]-0.5477225575051661*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 2.0);
  hamil_vy_nodal[2] = 0.5*pow(-(0.8366600265340755*vmap_vy[3])-0.3162277660168379*vmap_vy[2]+0.5477225575051661*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 2.0);
  hamil_vy_nodal[3] = 0.5*pow(1.8708286933869707*vmap_vy[3]+1.5811388300841895*vmap_vy[2]+1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 2.0);

  hamil[0] = 0.16666666666666666*hamil_vy_nodal[3]+0.16666666666666666*hamil_vx_nodal[3]+0.8333333333333334*hamil_vy_nodal[2]+0.8333333333333334*hamil_vx_nodal[2]+0.8333333333333334*hamil_vy_nodal[1]+0.8333333333333334*hamil_vx_nodal[1]+0.16666666666666666*hamil_vy_nodal[0]+0.16666666666666666*hamil_vx_nodal[0]; 
  hamil[1] = 0.2886751345948129*hamil_vx_nodal[3]+0.6454972243679029*hamil_vx_nodal[2]-0.6454972243679029*hamil_vx_nodal[1]-0.2886751345948129*hamil_vx_nodal[0]; 
  hamil[2] = 0.2886751345948129*hamil_vy_nodal[3]+0.6454972243679029*hamil_vy_nodal[2]-0.6454972243679029*hamil_vy_nodal[1]-0.2886751345948129*hamil_vy_nodal[0]; 
  hamil[3] = 0.0; 
  hamil[4] = 0.37267799624996495*hamil_vx_nodal[3]-0.37267799624996495*hamil_vx_nodal[2]-0.37267799624996495*hamil_vx_nodal[1]+0.37267799624996495*hamil_vx_nodal[0]; 
  hamil[5] = 0.37267799624996495*hamil_vy_nodal[3]-0.37267799624996495*hamil_vy_nodal[2]-0.37267799624996495*hamil_vy_nodal[1]+0.37267799624996495*hamil_vy_nodal[0]; 
  hamil[6] = 0.0; 
  hamil[7] = 0.0; 
  hamil[8] = 0.1889822365046136*hamil_vx_nodal[3]-0.4225771273642583*hamil_vx_nodal[2]+0.4225771273642583*hamil_vx_nodal[1]-0.1889822365046136*hamil_vx_nodal[0]; 
  hamil[9] = 0.1889822365046136*hamil_vy_nodal[3]-0.4225771273642583*hamil_vy_nodal[2]+0.4225771273642583*hamil_vy_nodal[1]-0.1889822365046136*hamil_vy_nodal[0]; 
  hamil[10] = 0.0; 
  hamil[11] = 0.0; 
  hamil[12] = 0.0; 
  hamil[13] = 0.0; 
  hamil[14] = 0.0; 
  hamil[15] = 0.0; 

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

} 
