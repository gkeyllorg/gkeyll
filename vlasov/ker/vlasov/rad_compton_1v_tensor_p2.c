#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void rad_compton_1v_tensor_p2(const double *w, const double *dxv, const double *vmap, 
  double t_cool, double p0, double* GKYL_RESTRICT rad) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // vmap: Velocity-space map for nonuniform meshes.
  // t_cool: Cooling time.
  // p0: Momentum to relax to.
  // rad: Radiation field (-(p - p0)/tcool for Compton, -sign(p)*p^4 = -|p|*p^3 for curvature).
 
  const double *vmap_vx = &vmap[0]; 
  double *rad_vx = &rad[0];
  double rad_nodal_vx[3] = {0.0};

  rad_nodal_vx[0] = -1.0/t_cool*(-(1.8708286933869707*vmap_vx[3])+1.5811388300841895*vmap_vx[2]-1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0] - p0);

  rad_nodal_vx[1] = -1.0/t_cool*(0.7071067811865475*vmap_vx[0]-0.7905694150420947*vmap_vx[2] - p0);

  rad_nodal_vx[2] = -1.0/t_cool*(1.8708286933869707*vmap_vx[3]+1.5811388300841895*vmap_vx[2]+1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0] - p0);

  rad_vx[0] = 0.2357022603955158*rad_nodal_vx[2]+0.9428090415820636*rad_nodal_vx[1]+0.2357022603955158*rad_nodal_vx[0]; 
  rad_vx[1] = 0.408248290463863*rad_nodal_vx[2]-0.408248290463863*rad_nodal_vx[0]; 
  rad_vx[2] = 0.21081851067789195*rad_nodal_vx[2]-0.421637021355784*rad_nodal_vx[1]+0.21081851067789195*rad_nodal_vx[0]; 

} 
