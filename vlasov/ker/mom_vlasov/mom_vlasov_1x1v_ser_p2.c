#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[4]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  const double *vmap_vx = &vmap[0]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  out[0] += (0.6324555320336759*f[5]*vmap_vx1_sq+0.7071067811865475*f[0]*vmap_vx1_sq+0.7071067811865475*f[0]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*f[2])*volFact; 
  out[1] += (0.632455532033676*f[7]*vmap_vx1_sq+0.7071067811865475*f[1]*vmap_vx1_sq+0.7071067811865475*f[1]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*f[3])*volFact; 
  out[2] += (0.7071067811865475*f[4]*vmap_vx1_sq+0.7071067811865475*f[4]*vmap_vx0_sq+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*f[6])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  const double *vmap_vx = &vmap[0]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx0_cu = vmap_vx[0]*vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vx1_cu = vmap_vx[1]*vmap_vx[1]*vmap_vx[1]; 
  out[0] += (1.3416407864998738*vmap_vx[0]*f[5]*vmap_vx1_sq+1.5*f[0]*vmap_vx[0]*vmap_vx1_sq+0.9*f[2]*vmap_vx1_cu+1.5*vmap_vx[1]*f[2]*vmap_vx0_sq+0.5*f[0]*vmap_vx0_cu)*volFact; 
  out[1] += (1.3416407864998738*vmap_vx[0]*f[7]*vmap_vx1_sq+1.5*vmap_vx[0]*f[1]*vmap_vx1_sq+0.9*f[3]*vmap_vx1_cu+1.5*vmap_vx[1]*f[3]*vmap_vx0_sq+0.5*f[1]*vmap_vx0_cu)*volFact; 
  out[2] += (1.5*vmap_vx[0]*f[4]*vmap_vx1_sq+0.8999999999999998*f[6]*vmap_vx1_cu+1.5*vmap_vx[1]*f[6]*vmap_vx0_sq+0.5*f[4]*vmap_vx0_cu)*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M1i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[0] += (3.872983346207417*f[2]*hamil[2]+1.7320508075688772*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[1] += (3.872983346207417*hamil[2]*f[3]+1.7320508075688772*f[1]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[2] += (3.872983346207417*hamil[2]*f[6]+1.7320508075688772*hamil[1]*f[4])*dv10*jacob_vx_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M2_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  out[0] += (hamil[2]*f[5]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[1] += (1.0000000000000002*hamil[2]*f[7]+hamil[1]*f[3]+hamil[0]*f[1])*volFact; 
  out[2] += (1.0000000000000002*hamil[1]*f[6]+hamil[0]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M3i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx,
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  double dH_dpx[3] = {0.0}; 
  dH_dpx[0] = 1.7320508075688772*hamil[1]*dv10*jacob_vx_inv; 
  dH_dpx[1] = 3.872983346207417*hamil[2]*dv10*jacob_vx_inv; 

  out[0] += (0.7071067811865475*dH_dpx[0]*hamil[2]*f[5]+0.6324555320336759*dH_dpx[1]*hamil[1]*f[5]+0.6324555320336759*dH_dpx[1]*f[2]*hamil[2]+0.7071067811865475*dH_dpx[0]*hamil[1]*f[2]+0.7071067811865475*hamil[0]*dH_dpx[1]*f[2]+0.7071067811865475*f[0]*dH_dpx[1]*hamil[1]+0.7071067811865475*dH_dpx[0]*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865475*dH_dpx[0]*hamil[2]*f[7]+0.632455532033676*dH_dpx[1]*hamil[1]*f[7]+0.6324555320336759*dH_dpx[1]*hamil[2]*f[3]+0.7071067811865475*dH_dpx[0]*hamil[1]*f[3]+0.7071067811865475*hamil[0]*dH_dpx[1]*f[3]+0.7071067811865475*dH_dpx[1]*f[1]*hamil[1]+0.7071067811865475*dH_dpx[0]*hamil[0]*f[1])*volFact; 
  out[2] += (0.632455532033676*dH_dpx[1]*hamil[2]*f[6]+0.7071067811865475*dH_dpx[0]*hamil[1]*f[6]+0.7071067811865475*hamil[0]*dH_dpx[1]*f[6]+0.7071067811865475*dH_dpx[1]*hamil[1]*f[4]+0.7071067811865475*dH_dpx[0]*hamil[0]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_five_moments_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[4]*volFact; 
  out[3] += (3.872983346207417*f[2]*hamil[2]+1.7320508075688772*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[4] += (3.872983346207417*hamil[2]*f[3]+1.7320508075688772*f[1]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[5] += (3.872983346207417*hamil[2]*f[6]+1.7320508075688772*hamil[1]*f[4])*dv10*jacob_vx_inv*volFact; 
  out[6] += (hamil[2]*f[5]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[7] += (1.0000000000000002*hamil[2]*f[7]+hamil[1]*f[3]+hamil[0]*f[1])*volFact; 
  out[8] += (1.0000000000000002*hamil[1]*f[6]+hamil[0]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*0.25; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += (5.477225575051662*f[2]*hamil[2]+2.4494897427831783*f[0]*hamil[1])*dv10*jacob_vx_inv*volFact; 
  out[2] += (1.4142135623730951*hamil[2]*f[5]+1.4142135623730951*hamil[1]*f[2]+1.4142135623730951*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M0_upper_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *vmap, const double *jacob_vel, double v_thresh, double f_thresh, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double *px = &vmap[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double f_nodes[9] = {0.0};
  double fquad = 0.0;
  if ((0.7071067811865475*px[0]-0.9486832980505137*px[1]) > v_thresh) { 
    fquad = -(0.5999999999999995*f[7])-0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])+0.9*f[3]-0.6708203932499369*(f[2]+f[1])+0.5*f[0]; 
    f_nodes[0] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.75*f[6]+0.4472135954999579*f[5]-0.5590169943749475*f[4]-0.6708203932499369*f[2]+0.5*f[0]; 
    f_nodes[3] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.5999999999999995*f[7]-0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])-0.9*f[3]-0.6708203932499369*f[2]+0.6708203932499369*f[1]+0.5*f[0]; 
    f_nodes[6] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
  } 
  if ((0.7071067811865475*px[0]) > v_thresh) { 
    fquad = 0.75*f[7]-0.5590169943749475*f[5]+0.4472135954999579*f[4]-0.6708203932499369*f[1]+0.5*f[0]; 
    f_nodes[1] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.5*f[0]-0.5590169943749475*(f[5]+f[4]); 
    f_nodes[4] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = -(0.75*f[7])-0.5590169943749475*f[5]+0.4472135954999579*f[4]+0.6708203932499369*f[1]+0.5*f[0]; 
    f_nodes[7] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
  } 
  if ((0.9486832980505137*px[1]+0.7071067811865475*px[0]) > v_thresh) { 
    fquad = -(0.5999999999999995*f[7])+0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])-0.9*f[3]+0.6708203932499369*f[2]-0.6708203932499369*f[1]+0.5*f[0]; 
    f_nodes[2] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = -(0.75*f[6])+0.4472135954999579*f[5]-0.5590169943749475*f[4]+0.6708203932499369*f[2]+0.5*f[0]; 
    f_nodes[5] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.5999999999999995*f[7]+0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])+0.9*f[3]+0.6708203932499369*(f[2]+f[1])+0.5*f[0]; 
    f_nodes[8] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
  } 
  out[0] += (0.21824283369955166*f_nodes[8]+0.34918853391928284*f_nodes[7]+0.21824283369955166*f_nodes[6]+0.34918853391928284*f_nodes[5]+0.5587016542708527*f_nodes[4]+0.34918853391928284*f_nodes[3]+0.21824283369955166*f_nodes[2]+0.34918853391928284*f_nodes[1]+0.21824283369955166*f_nodes[0])*volFact; 
  out[1] += (0.2928034870526277*f_nodes[8]+0.46848557928420453*f_nodes[7]+0.2928034870526277*f_nodes[6]-0.2928034870526277*f_nodes[2]-0.46848557928420453*f_nodes[1]-0.2928034870526277*f_nodes[0])*volFact; 
  out[2] += (0.19520232470175186*f_nodes[8]+0.312323719522803*f_nodes[7]+0.19520232470175186*f_nodes[6]-0.3904046494035038*f_nodes[5]-0.624647439045606*f_nodes[4]-0.3904046494035038*f_nodes[3]+0.19520232470175186*f_nodes[2]+0.312323719522803*f_nodes[1]+0.19520232470175186*f_nodes[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M0_lower_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *vmap, const double *jacob_vel, double v_thresh, double f_thresh, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double *px = &vmap[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double f_nodes[9] = {0.0};
  double fquad = 0.0;
  if ((0.7071067811865475*px[0]-0.9486832980505137*px[1]) < -v_thresh) { 
    fquad = -(0.5999999999999995*f[7])-0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])+0.9*f[3]-0.6708203932499369*(f[2]+f[1])+0.5*f[0]; 
    f_nodes[0] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.75*f[6]+0.4472135954999579*f[5]-0.5590169943749475*f[4]-0.6708203932499369*f[2]+0.5*f[0]; 
    f_nodes[3] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.5999999999999995*f[7]-0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])-0.9*f[3]-0.6708203932499369*f[2]+0.6708203932499369*f[1]+0.5*f[0]; 
    f_nodes[6] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
  } 
  if ((0.7071067811865475*px[0]) < -v_thresh) { 
    fquad = 0.75*f[7]-0.5590169943749475*f[5]+0.4472135954999579*f[4]-0.6708203932499369*f[1]+0.5*f[0]; 
    f_nodes[1] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.5*f[0]-0.5590169943749475*(f[5]+f[4]); 
    f_nodes[4] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = -(0.75*f[7])-0.5590169943749475*f[5]+0.4472135954999579*f[4]+0.6708203932499369*f[1]+0.5*f[0]; 
    f_nodes[7] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
  } 
  if ((0.9486832980505137*px[1]+0.7071067811865475*px[0]) < -v_thresh) { 
    fquad = -(0.5999999999999995*f[7])+0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])-0.9*f[3]+0.6708203932499369*f[2]-0.6708203932499369*f[1]+0.5*f[0]; 
    f_nodes[2] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = -(0.75*f[6])+0.4472135954999579*f[5]-0.5590169943749475*f[4]+0.6708203932499369*f[2]+0.5*f[0]; 
    f_nodes[5] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.5999999999999995*f[7]+0.5999999999999999*f[6]+0.4472135954999579*(f[5]+f[4])+0.9*f[3]+0.6708203932499369*(f[2]+f[1])+0.5*f[0]; 
    f_nodes[8] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
  } 
  out[0] += (0.21824283369955166*f_nodes[8]+0.34918853391928284*f_nodes[7]+0.21824283369955166*f_nodes[6]+0.34918853391928284*f_nodes[5]+0.5587016542708527*f_nodes[4]+0.34918853391928284*f_nodes[3]+0.21824283369955166*f_nodes[2]+0.34918853391928284*f_nodes[1]+0.21824283369955166*f_nodes[0])*volFact; 
  out[1] += (0.2928034870526277*f_nodes[8]+0.46848557928420453*f_nodes[7]+0.2928034870526277*f_nodes[6]-0.2928034870526277*f_nodes[2]-0.46848557928420453*f_nodes[1]-0.2928034870526277*f_nodes[0])*volFact; 
  out[2] += (0.19520232470175186*f_nodes[8]+0.312323719522803*f_nodes[7]+0.19520232470175186*f_nodes[6]-0.3904046494035038*f_nodes[5]-0.624647439045606*f_nodes[4]-0.3904046494035038*f_nodes[3]+0.19520232470175186*f_nodes[2]+0.312323719522803*f_nodes[1]+0.19520232470175186*f_nodes[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M1i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[0] += (2.7386127875258306*f[3]*hamil[7]+1.224744871391589*f[4]*hamil[6]+2.7386127875258306*f[2]*hamil[5]+1.224744871391589*f[1]*hamil[3]+1.224744871391589*f[0]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[1] += (2.4494897427831783*f[6]*hamil[7]+2.7386127875258306*f[2]*hamil[7]+1.0954451150103324*f[1]*hamil[6]+2.7386127875258306*f[3]*hamil[5]+1.0954451150103324*hamil[3]*f[4]+1.224744871391589*f[0]*hamil[3]+1.224744871391589*f[1]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[2] += (2.4494897427831783*f[3]*hamil[7]+0.7824607964359517*f[4]*hamil[6]+1.224744871391589*f[0]*hamil[6]+2.7386127875258306*hamil[5]*f[6]+1.224744871391589*hamil[2]*f[4]+1.0954451150103324*f[1]*hamil[3])*dv10*jacob_vx_inv*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M2_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  out[0] += (0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865475*f[5]*hamil[7]+0.7071067811865475*hamil[5]*f[7]+0.632455532033676*f[3]*hamil[6]+0.632455532033676*hamil[3]*f[6]+0.6324555320336759*f[1]*hamil[4]+0.6324555320336759*hamil[1]*f[4]+0.7071067811865475*f[2]*hamil[3]+0.7071067811865475*hamil[2]*f[3]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
  out[2] += (0.6324555320336759*f[7]*hamil[7]+0.45175395145262565*f[6]*hamil[6]+0.7071067811865475*f[2]*hamil[6]+0.7071067811865475*hamil[2]*f[6]+0.45175395145262565*f[4]*hamil[4]+0.7071067811865475*f[0]*hamil[4]+0.7071067811865475*hamil[0]*f[4]+0.6324555320336759*f[3]*hamil[3]+0.6324555320336759*f[1]*hamil[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_five_moments_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[4]*volFact; 
  out[3] += (2.7386127875258306*f[3]*hamil[7]+1.224744871391589*f[4]*hamil[6]+2.7386127875258306*f[2]*hamil[5]+1.224744871391589*f[1]*hamil[3]+1.224744871391589*f[0]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[4] += (2.4494897427831783*f[6]*hamil[7]+2.7386127875258306*f[2]*hamil[7]+1.0954451150103324*f[1]*hamil[6]+2.7386127875258306*f[3]*hamil[5]+1.0954451150103324*hamil[3]*f[4]+1.224744871391589*f[0]*hamil[3]+1.224744871391589*f[1]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[5] += (2.4494897427831783*f[3]*hamil[7]+0.7824607964359517*f[4]*hamil[6]+1.224744871391589*f[0]*hamil[6]+2.7386127875258306*hamil[5]*f[6]+1.224744871391589*hamil[2]*f[4]+1.0954451150103324*f[1]*hamil[3])*dv10*jacob_vx_inv*volFact; 
  out[6] += (0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[7] += (0.7071067811865475*f[5]*hamil[7]+0.7071067811865475*hamil[5]*f[7]+0.632455532033676*f[3]*hamil[6]+0.632455532033676*hamil[3]*f[6]+0.6324555320336759*f[1]*hamil[4]+0.6324555320336759*hamil[1]*f[4]+0.7071067811865475*f[2]*hamil[3]+0.7071067811865475*hamil[2]*f[3]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
  out[8] += (0.6324555320336759*f[7]*hamil[7]+0.45175395145262565*f[6]*hamil[6]+0.7071067811865475*f[2]*hamil[6]+0.7071067811865475*hamil[2]*f[6]+0.45175395145262565*f[4]*hamil[4]+0.7071067811865475*f[0]*hamil[4]+0.7071067811865475*hamil[0]*f[4]+0.6324555320336759*f[3]*hamil[3]+0.6324555320336759*f[1]*hamil[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_int_five_moments_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*0.25; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += (3.872983346207417*f[3]*hamil[7]+1.7320508075688774*f[4]*hamil[6]+3.872983346207417*f[2]*hamil[5]+1.7320508075688772*f[1]*hamil[3]+1.7320508075688772*f[0]*hamil[2])*dv10*jacob_vx_inv*volFact; 
  out[2] += (f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
