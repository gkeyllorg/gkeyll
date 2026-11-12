#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_sr_1v_tensor_p3(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // vmap: Velocity-space map (C^1 cubic for tensor bases, C^0 linear for Serendipity bases, stored in the same vdim*4 layout).
  // hamil: Particle Hamiltonian.
  // hamil_inv: Inverse particle Hamiltonian. Only meaningful for relativistic simulations; zeroed for the non-relativistic Hamiltonian (1/H is singular at v = 0 and unused).
 
  const double *vmap_vx = &vmap[0]; 
  double hamil_nodal[4] = {0.0};
  double hamil_inv_nodal[4] = {0.0};
  hamil_nodal[0] = sqrt(1.0 + pow(-(1.8708286933869707*vmap_vx[3])+1.5811388300841895*vmap_vx[2]-1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0));
  hamil_inv_nodal[0] = 1.0/hamil_nodal[0];
  hamil_nodal[1] = sqrt(1.0 + pow(0.8366600265340755*vmap_vx[3]-0.3162277660168379*vmap_vx[2]-0.5477225575051661*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0));
  hamil_inv_nodal[1] = 1.0/hamil_nodal[1];
  hamil_nodal[2] = sqrt(1.0 + pow(-(0.8366600265340755*vmap_vx[3])-0.3162277660168379*vmap_vx[2]+0.5477225575051661*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0));
  hamil_inv_nodal[2] = 1.0/hamil_nodal[2];
  hamil_nodal[3] = sqrt(1.0 + pow(1.8708286933869707*vmap_vx[3]+1.5811388300841895*vmap_vx[2]+1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0));
  hamil_inv_nodal[3] = 1.0/hamil_nodal[3];

  hamil[0] = 0.11785113019775789*hamil_nodal[3]+0.5892556509887895*hamil_nodal[2]+0.5892556509887895*hamil_nodal[1]+0.11785113019775789*hamil_nodal[0]; 
  hamil[1] = 0.20412414523193148*hamil_nodal[3]+0.4564354645876384*hamil_nodal[2]-0.4564354645876384*hamil_nodal[1]-0.20412414523193148*hamil_nodal[0]; 
  hamil[2] = 0.2635231383473649*hamil_nodal[3]-0.2635231383473649*hamil_nodal[2]-0.2635231383473649*hamil_nodal[1]+0.2635231383473649*hamil_nodal[0]; 
  hamil[3] = 0.13363062095621214*hamil_nodal[3]-0.29880715233359834*hamil_nodal[2]+0.29880715233359834*hamil_nodal[1]-0.13363062095621214*hamil_nodal[0]; 

  hamil_inv[0] = 0.11785113019775789*hamil_inv_nodal[3]+0.5892556509887895*hamil_inv_nodal[2]+0.5892556509887895*hamil_inv_nodal[1]+0.11785113019775789*hamil_inv_nodal[0]; 
  hamil_inv[1] = 0.20412414523193148*hamil_inv_nodal[3]+0.4564354645876384*hamil_inv_nodal[2]-0.4564354645876384*hamil_inv_nodal[1]-0.20412414523193148*hamil_inv_nodal[0]; 
  hamil_inv[2] = 0.2635231383473649*hamil_inv_nodal[3]-0.2635231383473649*hamil_inv_nodal[2]-0.2635231383473649*hamil_inv_nodal[1]+0.2635231383473649*hamil_inv_nodal[0]; 
  hamil_inv[3] = 0.13363062095621214*hamil_inv_nodal[3]-0.29880715233359834*hamil_inv_nodal[2]+0.29880715233359834*hamil_inv_nodal[1]-0.13363062095621214*hamil_inv_nodal[0]; 

} 
