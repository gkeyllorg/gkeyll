#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void rad_compton_2v_tensor_p1(const double *w, const double *dxv, const double *vmap, 
  double t_cool, double p0, double* GKYL_RESTRICT rad) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // vmap: Velocity-space map (C^1 cubic for tensor bases, C^0 linear for Serendipity bases, stored in the same vdim*4 layout).
  // t_cool: Cooling time.
  // p0: Momentum to relax to.
  // rad: Radiation field (-(p - p0)/tcool for Compton, -sign(p)*p^4 = -|p|*p^3 for curvature).
 
  const double *vmap_vx = &vmap[0]; 
  double *rad_vx = &rad[0];
  double rad_nodal_vx[4] = {0.0};
  const double *vmap_vy = &vmap[4]; 
  double *rad_vy = &rad[9];
  double rad_nodal_vy[4] = {0.0};

  rad_nodal_vx[0] = -1.0/t_cool*(-(1.8708286933869707*vmap_vx[3])+1.5811388300841895*vmap_vx[2]-1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0] - p0);
  rad_nodal_vy[0] = -1.0/t_cool*(-(1.8708286933869707*vmap_vy[3])+1.5811388300841895*vmap_vy[2]-1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0] - p0);

  rad_nodal_vx[1] = -1.0/t_cool*(1.8708286933869707*vmap_vx[3]+1.5811388300841895*vmap_vx[2]+1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0] - p0);
  rad_nodal_vy[1] = -1.0/t_cool*(-(1.8708286933869707*vmap_vy[3])+1.5811388300841895*vmap_vy[2]-1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0] - p0);

  rad_nodal_vx[2] = -1.0/t_cool*(-(1.8708286933869707*vmap_vx[3])+1.5811388300841895*vmap_vx[2]-1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0] - p0);
  rad_nodal_vy[2] = -1.0/t_cool*(1.8708286933869707*vmap_vy[3]+1.5811388300841895*vmap_vy[2]+1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0] - p0);

  rad_nodal_vx[3] = -1.0/t_cool*(1.8708286933869707*vmap_vx[3]+1.5811388300841895*vmap_vx[2]+1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0] - p0);
  rad_nodal_vy[3] = -1.0/t_cool*(1.8708286933869707*vmap_vy[3]+1.5811388300841895*vmap_vy[2]+1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0] - p0);

  rad_vx[0] = 0.5*rad_nodal_vx[3]+0.5*rad_nodal_vx[2]+0.5*rad_nodal_vx[1]+0.5*rad_nodal_vx[0]; 
  rad_vx[1] = 0.2886751345948129*rad_nodal_vx[3]-0.2886751345948129*rad_nodal_vx[2]+0.2886751345948129*rad_nodal_vx[1]-0.2886751345948129*rad_nodal_vx[0]; 
  rad_vx[2] = 0.2886751345948129*rad_nodal_vx[3]+0.2886751345948129*rad_nodal_vx[2]-0.2886751345948129*rad_nodal_vx[1]-0.2886751345948129*rad_nodal_vx[0]; 
  rad_vx[3] = 0.16666666666666666*rad_nodal_vx[3]-0.16666666666666666*rad_nodal_vx[2]-0.16666666666666666*rad_nodal_vx[1]+0.16666666666666666*rad_nodal_vx[0]; 

  rad_vy[0] = 0.5*rad_nodal_vy[3]+0.5*rad_nodal_vy[2]+0.5*rad_nodal_vy[1]+0.5*rad_nodal_vy[0]; 
  rad_vy[1] = 0.2886751345948129*rad_nodal_vy[3]-0.2886751345948129*rad_nodal_vy[2]+0.2886751345948129*rad_nodal_vy[1]-0.2886751345948129*rad_nodal_vy[0]; 
  rad_vy[2] = 0.2886751345948129*rad_nodal_vy[3]+0.2886751345948129*rad_nodal_vy[2]-0.2886751345948129*rad_nodal_vy[1]-0.2886751345948129*rad_nodal_vy[0]; 
  rad_vy[3] = 0.16666666666666666*rad_nodal_vy[3]-0.16666666666666666*rad_nodal_vy[2]-0.16666666666666666*rad_nodal_vy[1]+0.16666666666666666*rad_nodal_vy[0]; 

} 
