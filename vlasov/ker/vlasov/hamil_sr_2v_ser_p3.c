#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_sr_2v_ser_p3(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // vmap: Velocity-space map (C^1 cubic for tensor bases, C^0 linear for Serendipity bases, stored in the same vdim*4 layout).
  // hamil: Particle Hamiltonian.
  // hamil_inv: Inverse particle Hamiltonian. Only meaningful for relativistic simulations; zeroed for the non-relativistic Hamiltonian (1/H is singular at v = 0 and unused).
 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  double hamil_nodal[12] = {0.0};
  double hamil_inv_nodal[12] = {0.0};
  hamil_nodal[0] = sqrt(1.0 + pow(0.7071067811865475*vmap_vx[0]-1.224744871391589*vmap_vx[1], 2.0) + pow(0.7071067811865475*vmap_vy[0]-1.224744871391589*vmap_vy[1], 2.0));
  hamil_inv_nodal[0] = 1.0/hamil_nodal[0];
  hamil_nodal[1] = sqrt(1.0 + pow(0.7071067811865475*vmap_vx[0]-0.5477225575051661*vmap_vx[1], 2.0) + pow(0.7071067811865475*vmap_vy[0]-1.224744871391589*vmap_vy[1], 2.0));
  hamil_inv_nodal[1] = 1.0/hamil_nodal[1];
  hamil_nodal[2] = sqrt(1.0 + pow(0.5477225575051661*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0) + pow(0.7071067811865475*vmap_vy[0]-1.224744871391589*vmap_vy[1], 2.0));
  hamil_inv_nodal[2] = 1.0/hamil_nodal[2];
  hamil_nodal[3] = sqrt(1.0 + pow(1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0) + pow(0.7071067811865475*vmap_vy[0]-1.224744871391589*vmap_vy[1], 2.0));
  hamil_inv_nodal[3] = 1.0/hamil_nodal[3];
  hamil_nodal[4] = sqrt(1.0 + pow(0.7071067811865475*vmap_vx[0]-1.224744871391589*vmap_vx[1], 2.0) + pow(0.7071067811865475*vmap_vy[0]-0.5477225575051661*vmap_vy[1], 2.0));
  hamil_inv_nodal[4] = 1.0/hamil_nodal[4];
  hamil_nodal[5] = sqrt(1.0 + pow(1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0) + pow(0.7071067811865475*vmap_vy[0]-0.5477225575051661*vmap_vy[1], 2.0));
  hamil_inv_nodal[5] = 1.0/hamil_nodal[5];
  hamil_nodal[6] = sqrt(1.0 + pow(0.7071067811865475*vmap_vx[0]-1.224744871391589*vmap_vx[1], 2.0) + pow(0.5477225575051661*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 2.0));
  hamil_inv_nodal[6] = 1.0/hamil_nodal[6];
  hamil_nodal[7] = sqrt(1.0 + pow(1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0) + pow(0.5477225575051661*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 2.0));
  hamil_inv_nodal[7] = 1.0/hamil_nodal[7];
  hamil_nodal[8] = sqrt(1.0 + pow(0.7071067811865475*vmap_vx[0]-1.224744871391589*vmap_vx[1], 2.0) + pow(1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 2.0));
  hamil_inv_nodal[8] = 1.0/hamil_nodal[8];
  hamil_nodal[9] = sqrt(1.0 + pow(0.7071067811865475*vmap_vx[0]-0.5477225575051661*vmap_vx[1], 2.0) + pow(1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 2.0));
  hamil_inv_nodal[9] = 1.0/hamil_nodal[9];
  hamil_nodal[10] = sqrt(1.0 + pow(0.5477225575051661*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0) + pow(1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 2.0));
  hamil_inv_nodal[10] = 1.0/hamil_nodal[10];
  hamil_nodal[11] = sqrt(1.0 + pow(1.224744871391589*vmap_vx[1]+0.7071067811865475*vmap_vx[0], 2.0) + pow(1.224744871391589*vmap_vy[1]+0.7071067811865475*vmap_vy[0], 2.0));
  hamil_inv_nodal[11] = 1.0/hamil_nodal[11];

  hamil[0] = -(0.3333333333333333*hamil_nodal[11])+0.4166666666666667*hamil_nodal[10]+0.4166666666666667*hamil_nodal[9]-0.3333333333333333*hamil_nodal[8]+0.4166666666666667*hamil_nodal[7]+0.4166666666666667*hamil_nodal[6]+0.4166666666666667*hamil_nodal[5]+0.4166666666666667*hamil_nodal[4]-0.3333333333333333*hamil_nodal[3]+0.4166666666666667*hamil_nodal[2]+0.4166666666666667*hamil_nodal[1]-0.3333333333333333*hamil_nodal[0]; 
  hamil[1] = -(0.09622504486493764*hamil_nodal[11])+0.32274861218395146*hamil_nodal[10]-0.32274861218395146*hamil_nodal[9]+0.09622504486493764*hamil_nodal[8]+0.24056261216234412*hamil_nodal[7]-0.24056261216234412*hamil_nodal[6]+0.24056261216234412*hamil_nodal[5]-0.24056261216234412*hamil_nodal[4]-0.09622504486493764*hamil_nodal[3]+0.32274861218395146*hamil_nodal[2]-0.32274861218395146*hamil_nodal[1]+0.09622504486493764*hamil_nodal[0]; 
  hamil[2] = -(0.09622504486493764*hamil_nodal[11])+0.24056261216234412*hamil_nodal[10]+0.24056261216234412*hamil_nodal[9]-0.09622504486493764*hamil_nodal[8]+0.32274861218395146*hamil_nodal[7]+0.32274861218395146*hamil_nodal[6]-0.32274861218395146*hamil_nodal[5]-0.32274861218395146*hamil_nodal[4]+0.09622504486493764*hamil_nodal[3]-0.24056261216234412*hamil_nodal[2]-0.24056261216234412*hamil_nodal[1]+0.09622504486493764*hamil_nodal[0]; 
  hamil[3] = 0.18633899812498247*hamil_nodal[10]-0.18633899812498247*hamil_nodal[9]+0.18633899812498247*hamil_nodal[7]-0.18633899812498247*hamil_nodal[6]-0.18633899812498247*hamil_nodal[5]+0.18633899812498247*hamil_nodal[4]-0.18633899812498247*hamil_nodal[2]+0.18633899812498247*hamil_nodal[1]; 
  hamil[4] = 0.18633899812498247*hamil_nodal[11]-0.18633899812498247*hamil_nodal[10]-0.18633899812498247*hamil_nodal[9]+0.18633899812498247*hamil_nodal[8]+0.18633899812498247*hamil_nodal[3]-0.18633899812498247*hamil_nodal[2]-0.18633899812498247*hamil_nodal[1]+0.18633899812498247*hamil_nodal[0]; 
  hamil[5] = 0.18633899812498247*hamil_nodal[11]+0.18633899812498247*hamil_nodal[8]-0.18633899812498247*hamil_nodal[7]-0.18633899812498247*hamil_nodal[6]-0.18633899812498247*hamil_nodal[5]-0.18633899812498247*hamil_nodal[4]+0.18633899812498247*hamil_nodal[3]+0.18633899812498247*hamil_nodal[0]; 
  hamil[6] = 0.1075828707279838*hamil_nodal[11]-0.1075828707279838*hamil_nodal[10]-0.1075828707279838*hamil_nodal[9]+0.1075828707279838*hamil_nodal[8]-0.1075828707279838*hamil_nodal[3]+0.1075828707279838*hamil_nodal[2]+0.1075828707279838*hamil_nodal[1]-0.1075828707279838*hamil_nodal[0]; 
  hamil[7] = 0.1075828707279838*hamil_nodal[11]-0.1075828707279838*hamil_nodal[8]-0.1075828707279838*hamil_nodal[7]+0.1075828707279838*hamil_nodal[6]-0.1075828707279838*hamil_nodal[5]+0.1075828707279838*hamil_nodal[4]+0.1075828707279838*hamil_nodal[3]-0.1075828707279838*hamil_nodal[0]; 
  hamil[8] = 0.0944911182523068*hamil_nodal[11]-0.21128856368212914*hamil_nodal[10]+0.21128856368212914*hamil_nodal[9]-0.0944911182523068*hamil_nodal[8]+0.0944911182523068*hamil_nodal[3]-0.21128856368212914*hamil_nodal[2]+0.21128856368212914*hamil_nodal[1]-0.0944911182523068*hamil_nodal[0]; 
  hamil[9] = 0.0944911182523068*hamil_nodal[11]+0.0944911182523068*hamil_nodal[8]-0.21128856368212914*hamil_nodal[7]-0.21128856368212914*hamil_nodal[6]+0.21128856368212914*hamil_nodal[5]+0.21128856368212914*hamil_nodal[4]-0.0944911182523068*hamil_nodal[3]-0.0944911182523068*hamil_nodal[0]; 
  hamil[10] = 0.0545544725589981*hamil_nodal[11]-0.12198750911856666*hamil_nodal[10]+0.12198750911856666*hamil_nodal[9]-0.0545544725589981*hamil_nodal[8]-0.0545544725589981*hamil_nodal[3]+0.12198750911856666*hamil_nodal[2]-0.12198750911856666*hamil_nodal[1]+0.0545544725589981*hamil_nodal[0]; 
  hamil[11] = 0.0545544725589981*hamil_nodal[11]-0.0545544725589981*hamil_nodal[8]-0.12198750911856666*hamil_nodal[7]+0.12198750911856666*hamil_nodal[6]+0.12198750911856666*hamil_nodal[5]-0.12198750911856666*hamil_nodal[4]-0.0545544725589981*hamil_nodal[3]+0.0545544725589981*hamil_nodal[0]; 

  hamil_inv[0] = -(0.3333333333333333*hamil_inv_nodal[11])+0.4166666666666667*hamil_inv_nodal[10]+0.4166666666666667*hamil_inv_nodal[9]-0.3333333333333333*hamil_inv_nodal[8]+0.4166666666666667*hamil_inv_nodal[7]+0.4166666666666667*hamil_inv_nodal[6]+0.4166666666666667*hamil_inv_nodal[5]+0.4166666666666667*hamil_inv_nodal[4]-0.3333333333333333*hamil_inv_nodal[3]+0.4166666666666667*hamil_inv_nodal[2]+0.4166666666666667*hamil_inv_nodal[1]-0.3333333333333333*hamil_inv_nodal[0]; 
  hamil_inv[1] = -(0.09622504486493764*hamil_inv_nodal[11])+0.32274861218395146*hamil_inv_nodal[10]-0.32274861218395146*hamil_inv_nodal[9]+0.09622504486493764*hamil_inv_nodal[8]+0.24056261216234412*hamil_inv_nodal[7]-0.24056261216234412*hamil_inv_nodal[6]+0.24056261216234412*hamil_inv_nodal[5]-0.24056261216234412*hamil_inv_nodal[4]-0.09622504486493764*hamil_inv_nodal[3]+0.32274861218395146*hamil_inv_nodal[2]-0.32274861218395146*hamil_inv_nodal[1]+0.09622504486493764*hamil_inv_nodal[0]; 
  hamil_inv[2] = -(0.09622504486493764*hamil_inv_nodal[11])+0.24056261216234412*hamil_inv_nodal[10]+0.24056261216234412*hamil_inv_nodal[9]-0.09622504486493764*hamil_inv_nodal[8]+0.32274861218395146*hamil_inv_nodal[7]+0.32274861218395146*hamil_inv_nodal[6]-0.32274861218395146*hamil_inv_nodal[5]-0.32274861218395146*hamil_inv_nodal[4]+0.09622504486493764*hamil_inv_nodal[3]-0.24056261216234412*hamil_inv_nodal[2]-0.24056261216234412*hamil_inv_nodal[1]+0.09622504486493764*hamil_inv_nodal[0]; 
  hamil_inv[3] = 0.18633899812498247*hamil_inv_nodal[10]-0.18633899812498247*hamil_inv_nodal[9]+0.18633899812498247*hamil_inv_nodal[7]-0.18633899812498247*hamil_inv_nodal[6]-0.18633899812498247*hamil_inv_nodal[5]+0.18633899812498247*hamil_inv_nodal[4]-0.18633899812498247*hamil_inv_nodal[2]+0.18633899812498247*hamil_inv_nodal[1]; 
  hamil_inv[4] = 0.18633899812498247*hamil_inv_nodal[11]-0.18633899812498247*hamil_inv_nodal[10]-0.18633899812498247*hamil_inv_nodal[9]+0.18633899812498247*hamil_inv_nodal[8]+0.18633899812498247*hamil_inv_nodal[3]-0.18633899812498247*hamil_inv_nodal[2]-0.18633899812498247*hamil_inv_nodal[1]+0.18633899812498247*hamil_inv_nodal[0]; 
  hamil_inv[5] = 0.18633899812498247*hamil_inv_nodal[11]+0.18633899812498247*hamil_inv_nodal[8]-0.18633899812498247*hamil_inv_nodal[7]-0.18633899812498247*hamil_inv_nodal[6]-0.18633899812498247*hamil_inv_nodal[5]-0.18633899812498247*hamil_inv_nodal[4]+0.18633899812498247*hamil_inv_nodal[3]+0.18633899812498247*hamil_inv_nodal[0]; 
  hamil_inv[6] = 0.1075828707279838*hamil_inv_nodal[11]-0.1075828707279838*hamil_inv_nodal[10]-0.1075828707279838*hamil_inv_nodal[9]+0.1075828707279838*hamil_inv_nodal[8]-0.1075828707279838*hamil_inv_nodal[3]+0.1075828707279838*hamil_inv_nodal[2]+0.1075828707279838*hamil_inv_nodal[1]-0.1075828707279838*hamil_inv_nodal[0]; 
  hamil_inv[7] = 0.1075828707279838*hamil_inv_nodal[11]-0.1075828707279838*hamil_inv_nodal[8]-0.1075828707279838*hamil_inv_nodal[7]+0.1075828707279838*hamil_inv_nodal[6]-0.1075828707279838*hamil_inv_nodal[5]+0.1075828707279838*hamil_inv_nodal[4]+0.1075828707279838*hamil_inv_nodal[3]-0.1075828707279838*hamil_inv_nodal[0]; 
  hamil_inv[8] = 0.0944911182523068*hamil_inv_nodal[11]-0.21128856368212914*hamil_inv_nodal[10]+0.21128856368212914*hamil_inv_nodal[9]-0.0944911182523068*hamil_inv_nodal[8]+0.0944911182523068*hamil_inv_nodal[3]-0.21128856368212914*hamil_inv_nodal[2]+0.21128856368212914*hamil_inv_nodal[1]-0.0944911182523068*hamil_inv_nodal[0]; 
  hamil_inv[9] = 0.0944911182523068*hamil_inv_nodal[11]+0.0944911182523068*hamil_inv_nodal[8]-0.21128856368212914*hamil_inv_nodal[7]-0.21128856368212914*hamil_inv_nodal[6]+0.21128856368212914*hamil_inv_nodal[5]+0.21128856368212914*hamil_inv_nodal[4]-0.0944911182523068*hamil_inv_nodal[3]-0.0944911182523068*hamil_inv_nodal[0]; 
  hamil_inv[10] = 0.0545544725589981*hamil_inv_nodal[11]-0.12198750911856666*hamil_inv_nodal[10]+0.12198750911856666*hamil_inv_nodal[9]-0.0545544725589981*hamil_inv_nodal[8]-0.0545544725589981*hamil_inv_nodal[3]+0.12198750911856666*hamil_inv_nodal[2]-0.12198750911856666*hamil_inv_nodal[1]+0.0545544725589981*hamil_inv_nodal[0]; 
  hamil_inv[11] = 0.0545544725589981*hamil_inv_nodal[11]-0.0545544725589981*hamil_inv_nodal[8]-0.12198750911856666*hamil_inv_nodal[7]+0.12198750911856666*hamil_inv_nodal[6]+0.12198750911856666*hamil_inv_nodal[5]-0.12198750911856666*hamil_inv_nodal[4]-0.0545544725589981*hamil_inv_nodal[3]+0.0545544725589981*hamil_inv_nodal[0]; 

} 
