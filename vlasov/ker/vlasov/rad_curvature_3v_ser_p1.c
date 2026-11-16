#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void rad_curvature_3v_ser_p1(const double *w, const double *dxv, const double *vmap, 
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
  double rad_nodal_vx[8] = {0.0};
  const double *vmap_vy = &vmap[4]; 
  double *rad_vy = &rad[8];
  double rad_nodal_vy[8] = {0.0};
  const double *vmap_vz = &vmap[8]; 
  double *rad_vz = &rad[16];
  double rad_nodal_vz[8] = {0.0};

  rad_nodal_vx[0] = -1.0/t_cool*(pow(0.7071067811865475*vmap_vx[0]-1.224744871391589*vmap_vx[1], 3.0)*fabs(0.7071067811865475*vmap_vx[0]-1.224744871391589*vmap_vx[1]));
  rad_nodal_vy[0] = -1.0/t_cool*(pow(0.7071067811865475*vmap_vy[0]-1.224744871391589*vmap_vy[1], 3.0)*fabs(0.7071067811865475*vmap_vy[0]-1.224744871391589*vmap_vy[1]));
  rad_nodal_vz[0] = -1.0/t_cool*(pow(0.7071067811865475*vmap_vz[0]-1.224744871391589*vmap_vz[1], 3.0)*fabs(0.7071067811865475*vmap_vz[0]-1.224744871391589*vmap_vz[1]));

  rad_nodal_vx[1] = -1.0/t_cool*(pow(1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 3.0)*fabs(1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0]));
  rad_nodal_vy[1] = -1.0/t_cool*(pow(0.7071067811865475*vmap_vy[0]-1.224744871391589*vmap_vy[1], 3.0)*fabs(0.7071067811865475*vmap_vy[0]-1.224744871391589*vmap_vy[1]));
  rad_nodal_vz[1] = -1.0/t_cool*(pow(0.7071067811865475*vmap_vz[0]-1.224744871391589*vmap_vz[1], 3.0)*fabs(0.7071067811865475*vmap_vz[0]-1.224744871391589*vmap_vz[1]));

  rad_nodal_vx[2] = -1.0/t_cool*(pow(0.7071067811865475*vmap_vx[0]-1.224744871391589*vmap_vx[1], 3.0)*fabs(0.7071067811865475*vmap_vx[0]-1.224744871391589*vmap_vx[1]));
  rad_nodal_vy[2] = -1.0/t_cool*(pow(1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 3.0)*fabs(1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0]));
  rad_nodal_vz[2] = -1.0/t_cool*(pow(0.7071067811865475*vmap_vz[0]-1.224744871391589*vmap_vz[1], 3.0)*fabs(0.7071067811865475*vmap_vz[0]-1.224744871391589*vmap_vz[1]));

  rad_nodal_vx[3] = -1.0/t_cool*(pow(1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 3.0)*fabs(1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0]));
  rad_nodal_vy[3] = -1.0/t_cool*(pow(1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 3.0)*fabs(1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0]));
  rad_nodal_vz[3] = -1.0/t_cool*(pow(0.7071067811865475*vmap_vz[0]-1.224744871391589*vmap_vz[1], 3.0)*fabs(0.7071067811865475*vmap_vz[0]-1.224744871391589*vmap_vz[1]));

  rad_nodal_vx[4] = -1.0/t_cool*(pow(0.7071067811865475*vmap_vx[0]-1.224744871391589*vmap_vx[1], 3.0)*fabs(0.7071067811865475*vmap_vx[0]-1.224744871391589*vmap_vx[1]));
  rad_nodal_vy[4] = -1.0/t_cool*(pow(0.7071067811865475*vmap_vy[0]-1.224744871391589*vmap_vy[1], 3.0)*fabs(0.7071067811865475*vmap_vy[0]-1.224744871391589*vmap_vy[1]));
  rad_nodal_vz[4] = -1.0/t_cool*(pow(1.224744871391589*vmap_vz[1]+0.7071067811865475*vmap_vz[0], 3.0)*fabs(1.224744871391589*vmap_vz[1]+0.7071067811865475*vmap_vz[0]));

  rad_nodal_vx[5] = -1.0/t_cool*(pow(1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 3.0)*fabs(1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0]));
  rad_nodal_vy[5] = -1.0/t_cool*(pow(0.7071067811865475*vmap_vy[0]-1.224744871391589*vmap_vy[1], 3.0)*fabs(0.7071067811865475*vmap_vy[0]-1.224744871391589*vmap_vy[1]));
  rad_nodal_vz[5] = -1.0/t_cool*(pow(1.224744871391589*vmap_vz[1]+0.7071067811865475*vmap_vz[0], 3.0)*fabs(1.224744871391589*vmap_vz[1]+0.7071067811865475*vmap_vz[0]));

  rad_nodal_vx[6] = -1.0/t_cool*(pow(0.7071067811865475*vmap_vx[0]-1.224744871391589*vmap_vx[1], 3.0)*fabs(0.7071067811865475*vmap_vx[0]-1.224744871391589*vmap_vx[1]));
  rad_nodal_vy[6] = -1.0/t_cool*(pow(1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 3.0)*fabs(1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0]));
  rad_nodal_vz[6] = -1.0/t_cool*(pow(1.224744871391589*vmap_vz[1]+0.7071067811865475*vmap_vz[0], 3.0)*fabs(1.224744871391589*vmap_vz[1]+0.7071067811865475*vmap_vz[0]));

  rad_nodal_vx[7] = -1.0/t_cool*(pow(1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 3.0)*fabs(1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0]));
  rad_nodal_vy[7] = -1.0/t_cool*(pow(1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 3.0)*fabs(1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0]));
  rad_nodal_vz[7] = -1.0/t_cool*(pow(1.224744871391589*vmap_vz[1]+0.7071067811865475*vmap_vz[0], 3.0)*fabs(1.224744871391589*vmap_vz[1]+0.7071067811865475*vmap_vz[0]));

  rad_vx[0] = 0.3535533905932737*rad_nodal_vx[7]+0.3535533905932737*rad_nodal_vx[6]+0.3535533905932737*rad_nodal_vx[5]+0.3535533905932737*rad_nodal_vx[4]+0.3535533905932737*rad_nodal_vx[3]+0.3535533905932737*rad_nodal_vx[2]+0.3535533905932737*rad_nodal_vx[1]+0.3535533905932737*rad_nodal_vx[0]; 
  rad_vx[1] = 0.20412414523193148*rad_nodal_vx[7]-0.20412414523193148*rad_nodal_vx[6]+0.20412414523193148*rad_nodal_vx[5]-0.20412414523193148*rad_nodal_vx[4]+0.20412414523193148*rad_nodal_vx[3]-0.20412414523193148*rad_nodal_vx[2]+0.20412414523193148*rad_nodal_vx[1]-0.20412414523193148*rad_nodal_vx[0]; 
  rad_vx[2] = 0.20412414523193148*rad_nodal_vx[7]+0.20412414523193148*rad_nodal_vx[6]-0.20412414523193148*rad_nodal_vx[5]-0.20412414523193148*rad_nodal_vx[4]+0.20412414523193148*rad_nodal_vx[3]+0.20412414523193148*rad_nodal_vx[2]-0.20412414523193148*rad_nodal_vx[1]-0.20412414523193148*rad_nodal_vx[0]; 
  rad_vx[3] = 0.20412414523193148*rad_nodal_vx[7]+0.20412414523193148*rad_nodal_vx[6]+0.20412414523193148*rad_nodal_vx[5]+0.20412414523193148*rad_nodal_vx[4]-0.20412414523193148*rad_nodal_vx[3]-0.20412414523193148*rad_nodal_vx[2]-0.20412414523193148*rad_nodal_vx[1]-0.20412414523193148*rad_nodal_vx[0]; 
  rad_vx[4] = 0.11785113019775789*rad_nodal_vx[7]-0.11785113019775789*rad_nodal_vx[6]-0.11785113019775789*rad_nodal_vx[5]+0.11785113019775789*rad_nodal_vx[4]+0.11785113019775789*rad_nodal_vx[3]-0.11785113019775789*rad_nodal_vx[2]-0.11785113019775789*rad_nodal_vx[1]+0.11785113019775789*rad_nodal_vx[0]; 
  rad_vx[5] = 0.11785113019775789*rad_nodal_vx[7]-0.11785113019775789*rad_nodal_vx[6]+0.11785113019775789*rad_nodal_vx[5]-0.11785113019775789*rad_nodal_vx[4]-0.11785113019775789*rad_nodal_vx[3]+0.11785113019775789*rad_nodal_vx[2]-0.11785113019775789*rad_nodal_vx[1]+0.11785113019775789*rad_nodal_vx[0]; 
  rad_vx[6] = 0.11785113019775789*rad_nodal_vx[7]+0.11785113019775789*rad_nodal_vx[6]-0.11785113019775789*rad_nodal_vx[5]-0.11785113019775789*rad_nodal_vx[4]-0.11785113019775789*rad_nodal_vx[3]-0.11785113019775789*rad_nodal_vx[2]+0.11785113019775789*rad_nodal_vx[1]+0.11785113019775789*rad_nodal_vx[0]; 
  rad_vx[7] = 0.06804138174397717*rad_nodal_vx[7]-0.06804138174397717*rad_nodal_vx[6]-0.06804138174397717*rad_nodal_vx[5]+0.06804138174397717*rad_nodal_vx[4]-0.06804138174397717*rad_nodal_vx[3]+0.06804138174397717*rad_nodal_vx[2]+0.06804138174397717*rad_nodal_vx[1]-0.06804138174397717*rad_nodal_vx[0]; 

  rad_vy[0] = 0.3535533905932737*rad_nodal_vy[7]+0.3535533905932737*rad_nodal_vy[6]+0.3535533905932737*rad_nodal_vy[5]+0.3535533905932737*rad_nodal_vy[4]+0.3535533905932737*rad_nodal_vy[3]+0.3535533905932737*rad_nodal_vy[2]+0.3535533905932737*rad_nodal_vy[1]+0.3535533905932737*rad_nodal_vy[0]; 
  rad_vy[1] = 0.20412414523193148*rad_nodal_vy[7]-0.20412414523193148*rad_nodal_vy[6]+0.20412414523193148*rad_nodal_vy[5]-0.20412414523193148*rad_nodal_vy[4]+0.20412414523193148*rad_nodal_vy[3]-0.20412414523193148*rad_nodal_vy[2]+0.20412414523193148*rad_nodal_vy[1]-0.20412414523193148*rad_nodal_vy[0]; 
  rad_vy[2] = 0.20412414523193148*rad_nodal_vy[7]+0.20412414523193148*rad_nodal_vy[6]-0.20412414523193148*rad_nodal_vy[5]-0.20412414523193148*rad_nodal_vy[4]+0.20412414523193148*rad_nodal_vy[3]+0.20412414523193148*rad_nodal_vy[2]-0.20412414523193148*rad_nodal_vy[1]-0.20412414523193148*rad_nodal_vy[0]; 
  rad_vy[3] = 0.20412414523193148*rad_nodal_vy[7]+0.20412414523193148*rad_nodal_vy[6]+0.20412414523193148*rad_nodal_vy[5]+0.20412414523193148*rad_nodal_vy[4]-0.20412414523193148*rad_nodal_vy[3]-0.20412414523193148*rad_nodal_vy[2]-0.20412414523193148*rad_nodal_vy[1]-0.20412414523193148*rad_nodal_vy[0]; 
  rad_vy[4] = 0.11785113019775789*rad_nodal_vy[7]-0.11785113019775789*rad_nodal_vy[6]-0.11785113019775789*rad_nodal_vy[5]+0.11785113019775789*rad_nodal_vy[4]+0.11785113019775789*rad_nodal_vy[3]-0.11785113019775789*rad_nodal_vy[2]-0.11785113019775789*rad_nodal_vy[1]+0.11785113019775789*rad_nodal_vy[0]; 
  rad_vy[5] = 0.11785113019775789*rad_nodal_vy[7]-0.11785113019775789*rad_nodal_vy[6]+0.11785113019775789*rad_nodal_vy[5]-0.11785113019775789*rad_nodal_vy[4]-0.11785113019775789*rad_nodal_vy[3]+0.11785113019775789*rad_nodal_vy[2]-0.11785113019775789*rad_nodal_vy[1]+0.11785113019775789*rad_nodal_vy[0]; 
  rad_vy[6] = 0.11785113019775789*rad_nodal_vy[7]+0.11785113019775789*rad_nodal_vy[6]-0.11785113019775789*rad_nodal_vy[5]-0.11785113019775789*rad_nodal_vy[4]-0.11785113019775789*rad_nodal_vy[3]-0.11785113019775789*rad_nodal_vy[2]+0.11785113019775789*rad_nodal_vy[1]+0.11785113019775789*rad_nodal_vy[0]; 
  rad_vy[7] = 0.06804138174397717*rad_nodal_vy[7]-0.06804138174397717*rad_nodal_vy[6]-0.06804138174397717*rad_nodal_vy[5]+0.06804138174397717*rad_nodal_vy[4]-0.06804138174397717*rad_nodal_vy[3]+0.06804138174397717*rad_nodal_vy[2]+0.06804138174397717*rad_nodal_vy[1]-0.06804138174397717*rad_nodal_vy[0]; 

  rad_vz[0] = 0.3535533905932737*rad_nodal_vz[7]+0.3535533905932737*rad_nodal_vz[6]+0.3535533905932737*rad_nodal_vz[5]+0.3535533905932737*rad_nodal_vz[4]+0.3535533905932737*rad_nodal_vz[3]+0.3535533905932737*rad_nodal_vz[2]+0.3535533905932737*rad_nodal_vz[1]+0.3535533905932737*rad_nodal_vz[0]; 
  rad_vz[1] = 0.20412414523193148*rad_nodal_vz[7]-0.20412414523193148*rad_nodal_vz[6]+0.20412414523193148*rad_nodal_vz[5]-0.20412414523193148*rad_nodal_vz[4]+0.20412414523193148*rad_nodal_vz[3]-0.20412414523193148*rad_nodal_vz[2]+0.20412414523193148*rad_nodal_vz[1]-0.20412414523193148*rad_nodal_vz[0]; 
  rad_vz[2] = 0.20412414523193148*rad_nodal_vz[7]+0.20412414523193148*rad_nodal_vz[6]-0.20412414523193148*rad_nodal_vz[5]-0.20412414523193148*rad_nodal_vz[4]+0.20412414523193148*rad_nodal_vz[3]+0.20412414523193148*rad_nodal_vz[2]-0.20412414523193148*rad_nodal_vz[1]-0.20412414523193148*rad_nodal_vz[0]; 
  rad_vz[3] = 0.20412414523193148*rad_nodal_vz[7]+0.20412414523193148*rad_nodal_vz[6]+0.20412414523193148*rad_nodal_vz[5]+0.20412414523193148*rad_nodal_vz[4]-0.20412414523193148*rad_nodal_vz[3]-0.20412414523193148*rad_nodal_vz[2]-0.20412414523193148*rad_nodal_vz[1]-0.20412414523193148*rad_nodal_vz[0]; 
  rad_vz[4] = 0.11785113019775789*rad_nodal_vz[7]-0.11785113019775789*rad_nodal_vz[6]-0.11785113019775789*rad_nodal_vz[5]+0.11785113019775789*rad_nodal_vz[4]+0.11785113019775789*rad_nodal_vz[3]-0.11785113019775789*rad_nodal_vz[2]-0.11785113019775789*rad_nodal_vz[1]+0.11785113019775789*rad_nodal_vz[0]; 
  rad_vz[5] = 0.11785113019775789*rad_nodal_vz[7]-0.11785113019775789*rad_nodal_vz[6]+0.11785113019775789*rad_nodal_vz[5]-0.11785113019775789*rad_nodal_vz[4]-0.11785113019775789*rad_nodal_vz[3]+0.11785113019775789*rad_nodal_vz[2]-0.11785113019775789*rad_nodal_vz[1]+0.11785113019775789*rad_nodal_vz[0]; 
  rad_vz[6] = 0.11785113019775789*rad_nodal_vz[7]+0.11785113019775789*rad_nodal_vz[6]-0.11785113019775789*rad_nodal_vz[5]-0.11785113019775789*rad_nodal_vz[4]-0.11785113019775789*rad_nodal_vz[3]-0.11785113019775789*rad_nodal_vz[2]+0.11785113019775789*rad_nodal_vz[1]+0.11785113019775789*rad_nodal_vz[0]; 
  rad_vz[7] = 0.06804138174397717*rad_nodal_vz[7]-0.06804138174397717*rad_nodal_vz[6]-0.06804138174397717*rad_nodal_vz[5]+0.06804138174397717*rad_nodal_vz[4]-0.06804138174397717*rad_nodal_vz[3]+0.06804138174397717*rad_nodal_vz[2]+0.06804138174397717*rad_nodal_vz[1]-0.06804138174397717*rad_nodal_vz[0]; 

} 
