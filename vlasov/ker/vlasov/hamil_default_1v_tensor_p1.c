#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_default_1v_tensor_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // vmap: Velocity-space map (C^1 cubic for tensor bases, C^0 linear for Serendipity bases, stored in the same vdim*4 layout).
  // hamil: Particle Hamiltonian.
  // hamil_inv: Inverse particle Hamiltonian. Only meaningful for relativistic simulations; zeroed for the non-relativistic Hamiltonian (1/H is singular at v = 0 and unused).
 
  const double *vmap_vx = &vmap[0]; 
  double hamil_vx_nodal[3] = {0.0};
  hamil_vx_nodal[0] = 0.5*pow(-(1.8708286933869707*vmap_vx[3])+1.5811388300841895*vmap_vx[2]-1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0);
  hamil_vx_nodal[1] = 0.5*pow(0.7071067811865475*vmap_vx[0]-0.7905694150420947*vmap_vx[2], 2.0);
  hamil_vx_nodal[2] = 0.5*pow(1.8708286933869707*vmap_vx[3]+1.5811388300841895*vmap_vx[2]+1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0);

  hamil[0] = 0.2357022603955158*hamil_vx_nodal[2]+0.9428090415820636*hamil_vx_nodal[1]+0.2357022603955158*hamil_vx_nodal[0]; 
  hamil[1] = 0.408248290463863*hamil_vx_nodal[2]-0.408248290463863*hamil_vx_nodal[0]; 
  hamil[2] = 0.21081851067789195*hamil_vx_nodal[2]-0.421637021355784*hamil_vx_nodal[1]+0.21081851067789195*hamil_vx_nodal[0]; 

  hamil_inv[0] = 0.0; 
  hamil_inv[1] = 0.0; 
  hamil_inv[2] = 0.0; 

} 
