#include <gkyl_mat.h> 
#include <gkyl_euler_kernels.h> 
GKYL_CU_DH int fluid_vars_u_set_2x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *fluid) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A:     preallocated LHS matrix. 
  // rhs:   preallocated RHS vector. 
  // fluid: [rho, rho ux, rho uy, rho uz...], Fluid input state vector.
  //        only need rho and momentum to get flow velocity independent of fluid system. 
  //        (isothermal Euler, Euler, Ten moment). 

  struct gkyl_mat A_ux = gkyl_nmat_get(A, count); 
  struct gkyl_mat A_uy = gkyl_nmat_get(A, count+1); 
  struct gkyl_mat A_uz = gkyl_nmat_get(A, count+2); 
  struct gkyl_mat rhs_ux = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_uy = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_uz = gkyl_nmat_get(rhs, count+2); 
  // Clear matrix and rhs for each component of flow velocity being solved for 
  gkyl_mat_clear(&A_ux, 0.0); gkyl_mat_clear(&rhs_ux, 0.0); 
  gkyl_mat_clear(&A_uy, 0.0); gkyl_mat_clear(&rhs_uy, 0.0); 
  gkyl_mat_clear(&A_uz, 0.0); gkyl_mat_clear(&rhs_uz, 0.0); 
  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[9]; 
  const double *rhouy = &fluid[18]; 
  const double *rhouz = &fluid[27]; 

  int cell_avg = 0;
  // Check if rho < 0 at control points. 
  if (2.5*rho[8]-1.936491673103709*rho[7]-1.936491673103709*rho[6]+1.118033988749895*rho[5]+1.118033988749895*rho[4]+1.5*rho[3]-0.8660254037844386*rho[2]-0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if ((-1.25*rho[8])+0.9682458365518543*rho[6]+1.118033988749895*rho[5]-0.5590169943749475*rho[4]-0.8660254037844386*rho[2]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if (2.5*rho[8]+1.936491673103709*rho[7]-1.936491673103709*rho[6]+1.118033988749895*rho[5]+1.118033988749895*rho[4]-1.5*rho[3]-0.8660254037844386*rho[2]+0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if ((-1.25*rho[8])+0.9682458365518543*rho[7]-0.5590169943749475*rho[5]+1.118033988749895*rho[4]-0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if (0.625*rho[8]-0.5590169943749475*rho[5]-0.5590169943749475*rho[4]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if ((-1.25*rho[8])-0.9682458365518543*rho[7]-0.5590169943749475*rho[5]+1.118033988749895*rho[4]+0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if (2.5*rho[8]-1.936491673103709*rho[7]+1.936491673103709*rho[6]+1.118033988749895*rho[5]+1.118033988749895*rho[4]-1.5*rho[3]+0.8660254037844386*rho[2]-0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if ((-1.25*rho[8])-0.9682458365518543*rho[6]+1.118033988749895*rho[5]-0.5590169943749475*rho[4]+0.8660254037844386*rho[2]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if (2.5*rho[8]+1.936491673103709*rho[7]+1.936491673103709*rho[6]+1.118033988749895*rho[5]+1.118033988749895*rho[4]+1.5*rho[3]+0.8660254037844386*rho[2]+0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
 
  gkyl_mat_set(&rhs_ux,0,0,rhoux[0]); 
  gkyl_mat_set(&rhs_uy,0,0,rhouy[0]); 
  gkyl_mat_set(&rhs_uz,0,0,rhouz[0]); 
  gkyl_mat_set(&rhs_ux,1,0,rhoux[1]); 
  gkyl_mat_set(&rhs_uy,1,0,rhouy[1]); 
  gkyl_mat_set(&rhs_uz,1,0,rhouz[1]); 
  gkyl_mat_set(&rhs_ux,2,0,rhoux[2]); 
  gkyl_mat_set(&rhs_uy,2,0,rhouy[2]); 
  gkyl_mat_set(&rhs_uz,2,0,rhouz[2]); 
  gkyl_mat_set(&rhs_ux,3,0,rhoux[3]); 
  gkyl_mat_set(&rhs_uy,3,0,rhouy[3]); 
  gkyl_mat_set(&rhs_uz,3,0,rhouz[3]); 
  gkyl_mat_set(&rhs_ux,4,0,rhoux[4]); 
  gkyl_mat_set(&rhs_uy,4,0,rhouy[4]); 
  gkyl_mat_set(&rhs_uz,4,0,rhouz[4]); 
  gkyl_mat_set(&rhs_ux,5,0,rhoux[5]); 
  gkyl_mat_set(&rhs_uy,5,0,rhouy[5]); 
  gkyl_mat_set(&rhs_uz,5,0,rhouz[5]); 
  gkyl_mat_set(&rhs_ux,6,0,rhoux[6]); 
  gkyl_mat_set(&rhs_uy,6,0,rhouy[6]); 
  gkyl_mat_set(&rhs_uz,6,0,rhouz[6]); 
  gkyl_mat_set(&rhs_ux,7,0,rhoux[7]); 
  gkyl_mat_set(&rhs_uy,7,0,rhouy[7]); 
  gkyl_mat_set(&rhs_uz,7,0,rhouz[7]); 
  gkyl_mat_set(&rhs_ux,8,0,rhoux[8]); 
  gkyl_mat_set(&rhs_uy,8,0,rhouy[8]); 
  gkyl_mat_set(&rhs_uz,8,0,rhouz[8]); 
 
  double temp_rho = 0.0; 
  temp_rho = 0.5*rho[0]; 
  gkyl_mat_set(&A_ux,0,0,temp_rho); 
  gkyl_mat_set(&A_uy,0,0,temp_rho); 
  gkyl_mat_set(&A_uz,0,0,temp_rho); 
 
  temp_rho = 0.5*rho[1]; 
  gkyl_mat_set(&A_ux,0,1,temp_rho); 
  gkyl_mat_set(&A_uy,0,1,temp_rho); 
  gkyl_mat_set(&A_uz,0,1,temp_rho); 
 
  temp_rho = 0.5*rho[2]; 
  gkyl_mat_set(&A_ux,0,2,temp_rho); 
  gkyl_mat_set(&A_uy,0,2,temp_rho); 
  gkyl_mat_set(&A_uz,0,2,temp_rho); 
 
  temp_rho = 0.5*rho[3]; 
  gkyl_mat_set(&A_ux,0,3,temp_rho); 
  gkyl_mat_set(&A_uy,0,3,temp_rho); 
  gkyl_mat_set(&A_uz,0,3,temp_rho); 
 
  temp_rho = 0.5*rho[4]; 
  gkyl_mat_set(&A_ux,0,4,temp_rho); 
  gkyl_mat_set(&A_uy,0,4,temp_rho); 
  gkyl_mat_set(&A_uz,0,4,temp_rho); 
 
  temp_rho = 0.5*rho[5]; 
  gkyl_mat_set(&A_ux,0,5,temp_rho); 
  gkyl_mat_set(&A_uy,0,5,temp_rho); 
  gkyl_mat_set(&A_uz,0,5,temp_rho); 
 
  temp_rho = 0.5*rho[6]; 
  gkyl_mat_set(&A_ux,0,6,temp_rho); 
  gkyl_mat_set(&A_uy,0,6,temp_rho); 
  gkyl_mat_set(&A_uz,0,6,temp_rho); 
 
  temp_rho = 0.5*rho[7]; 
  gkyl_mat_set(&A_ux,0,7,temp_rho); 
  gkyl_mat_set(&A_uy,0,7,temp_rho); 
  gkyl_mat_set(&A_uz,0,7,temp_rho); 
 
  temp_rho = 0.5*rho[8]; 
  gkyl_mat_set(&A_ux,0,8,temp_rho); 
  gkyl_mat_set(&A_uy,0,8,temp_rho); 
  gkyl_mat_set(&A_uz,0,8,temp_rho); 
 
  temp_rho = 0.5*rho[1]; 
  gkyl_mat_set(&A_ux,1,0,temp_rho); 
  gkyl_mat_set(&A_uy,1,0,temp_rho); 
  gkyl_mat_set(&A_uz,1,0,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[4]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,1,1,temp_rho); 
  gkyl_mat_set(&A_uy,1,1,temp_rho); 
  gkyl_mat_set(&A_uz,1,1,temp_rho); 
 
  temp_rho = 0.5*rho[3]; 
  gkyl_mat_set(&A_ux,1,2,temp_rho); 
  gkyl_mat_set(&A_uy,1,2,temp_rho); 
  gkyl_mat_set(&A_uz,1,2,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[6]+0.5*rho[2]; 
  gkyl_mat_set(&A_ux,1,3,temp_rho); 
  gkyl_mat_set(&A_uy,1,3,temp_rho); 
  gkyl_mat_set(&A_uz,1,3,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[1]; 
  gkyl_mat_set(&A_ux,1,4,temp_rho); 
  gkyl_mat_set(&A_uy,1,4,temp_rho); 
  gkyl_mat_set(&A_uz,1,4,temp_rho); 
 
  temp_rho = 0.5000000000000001*rho[7]; 
  gkyl_mat_set(&A_ux,1,5,temp_rho); 
  gkyl_mat_set(&A_uy,1,5,temp_rho); 
  gkyl_mat_set(&A_uz,1,5,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[3]; 
  gkyl_mat_set(&A_ux,1,6,temp_rho); 
  gkyl_mat_set(&A_uy,1,6,temp_rho); 
  gkyl_mat_set(&A_uz,1,6,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[8]+0.5000000000000001*rho[5]; 
  gkyl_mat_set(&A_ux,1,7,temp_rho); 
  gkyl_mat_set(&A_uy,1,7,temp_rho); 
  gkyl_mat_set(&A_uz,1,7,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[7]; 
  gkyl_mat_set(&A_ux,1,8,temp_rho); 
  gkyl_mat_set(&A_uy,1,8,temp_rho); 
  gkyl_mat_set(&A_uz,1,8,temp_rho); 
 
  temp_rho = 0.5*rho[2]; 
  gkyl_mat_set(&A_ux,2,0,temp_rho); 
  gkyl_mat_set(&A_uy,2,0,temp_rho); 
  gkyl_mat_set(&A_uz,2,0,temp_rho); 
 
  temp_rho = 0.5*rho[3]; 
  gkyl_mat_set(&A_ux,2,1,temp_rho); 
  gkyl_mat_set(&A_uy,2,1,temp_rho); 
  gkyl_mat_set(&A_uz,2,1,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[5]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,2,2,temp_rho); 
  gkyl_mat_set(&A_uy,2,2,temp_rho); 
  gkyl_mat_set(&A_uz,2,2,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[7]+0.5*rho[1]; 
  gkyl_mat_set(&A_ux,2,3,temp_rho); 
  gkyl_mat_set(&A_uy,2,3,temp_rho); 
  gkyl_mat_set(&A_uz,2,3,temp_rho); 
 
  temp_rho = 0.5000000000000001*rho[6]; 
  gkyl_mat_set(&A_ux,2,4,temp_rho); 
  gkyl_mat_set(&A_uy,2,4,temp_rho); 
  gkyl_mat_set(&A_uz,2,4,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[2]; 
  gkyl_mat_set(&A_ux,2,5,temp_rho); 
  gkyl_mat_set(&A_uy,2,5,temp_rho); 
  gkyl_mat_set(&A_uz,2,5,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[8]+0.5000000000000001*rho[4]; 
  gkyl_mat_set(&A_ux,2,6,temp_rho); 
  gkyl_mat_set(&A_uy,2,6,temp_rho); 
  gkyl_mat_set(&A_uz,2,6,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[3]; 
  gkyl_mat_set(&A_ux,2,7,temp_rho); 
  gkyl_mat_set(&A_uy,2,7,temp_rho); 
  gkyl_mat_set(&A_uz,2,7,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[6]; 
  gkyl_mat_set(&A_ux,2,8,temp_rho); 
  gkyl_mat_set(&A_uy,2,8,temp_rho); 
  gkyl_mat_set(&A_uz,2,8,temp_rho); 
 
  temp_rho = 0.5*rho[3]; 
  gkyl_mat_set(&A_ux,3,0,temp_rho); 
  gkyl_mat_set(&A_uy,3,0,temp_rho); 
  gkyl_mat_set(&A_uz,3,0,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[6]+0.5*rho[2]; 
  gkyl_mat_set(&A_ux,3,1,temp_rho); 
  gkyl_mat_set(&A_uy,3,1,temp_rho); 
  gkyl_mat_set(&A_uz,3,1,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[7]+0.5*rho[1]; 
  gkyl_mat_set(&A_ux,3,2,temp_rho); 
  gkyl_mat_set(&A_uy,3,2,temp_rho); 
  gkyl_mat_set(&A_uz,3,2,temp_rho); 
 
  temp_rho = 0.4*rho[8]+0.4472135954999579*rho[5]+0.4472135954999579*rho[4]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,3,3,temp_rho); 
  gkyl_mat_set(&A_uy,3,3,temp_rho); 
  gkyl_mat_set(&A_uz,3,3,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[3]; 
  gkyl_mat_set(&A_ux,3,4,temp_rho); 
  gkyl_mat_set(&A_uy,3,4,temp_rho); 
  gkyl_mat_set(&A_uz,3,4,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[3]; 
  gkyl_mat_set(&A_ux,3,5,temp_rho); 
  gkyl_mat_set(&A_uy,3,5,temp_rho); 
  gkyl_mat_set(&A_uz,3,5,temp_rho); 
 
  temp_rho = 0.4*rho[7]+0.447213595499958*rho[1]; 
  gkyl_mat_set(&A_ux,3,6,temp_rho); 
  gkyl_mat_set(&A_uy,3,6,temp_rho); 
  gkyl_mat_set(&A_uz,3,6,temp_rho); 
 
  temp_rho = 0.4*rho[6]+0.447213595499958*rho[2]; 
  gkyl_mat_set(&A_ux,3,7,temp_rho); 
  gkyl_mat_set(&A_uy,3,7,temp_rho); 
  gkyl_mat_set(&A_uz,3,7,temp_rho); 
 
  temp_rho = 0.4*rho[3]; 
  gkyl_mat_set(&A_ux,3,8,temp_rho); 
  gkyl_mat_set(&A_uy,3,8,temp_rho); 
  gkyl_mat_set(&A_uz,3,8,temp_rho); 
 
  temp_rho = 0.5*rho[4]; 
  gkyl_mat_set(&A_ux,4,0,temp_rho); 
  gkyl_mat_set(&A_uy,4,0,temp_rho); 
  gkyl_mat_set(&A_uz,4,0,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[1]; 
  gkyl_mat_set(&A_ux,4,1,temp_rho); 
  gkyl_mat_set(&A_uy,4,1,temp_rho); 
  gkyl_mat_set(&A_uz,4,1,temp_rho); 
 
  temp_rho = 0.5000000000000001*rho[6]; 
  gkyl_mat_set(&A_ux,4,2,temp_rho); 
  gkyl_mat_set(&A_uy,4,2,temp_rho); 
  gkyl_mat_set(&A_uz,4,2,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[3]; 
  gkyl_mat_set(&A_ux,4,3,temp_rho); 
  gkyl_mat_set(&A_uy,4,3,temp_rho); 
  gkyl_mat_set(&A_uz,4,3,temp_rho); 
 
  temp_rho = 0.31943828249997*rho[4]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,4,4,temp_rho); 
  gkyl_mat_set(&A_uy,4,4,temp_rho); 
  gkyl_mat_set(&A_uz,4,4,temp_rho); 
 
  temp_rho = 0.5*rho[8]; 
  gkyl_mat_set(&A_ux,4,5,temp_rho); 
  gkyl_mat_set(&A_uy,4,5,temp_rho); 
  gkyl_mat_set(&A_uz,4,5,temp_rho); 
 
  temp_rho = 0.31943828249997*rho[6]+0.5000000000000001*rho[2]; 
  gkyl_mat_set(&A_ux,4,6,temp_rho); 
  gkyl_mat_set(&A_uy,4,6,temp_rho); 
  gkyl_mat_set(&A_uz,4,6,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[7]; 
  gkyl_mat_set(&A_ux,4,7,temp_rho); 
  gkyl_mat_set(&A_uy,4,7,temp_rho); 
  gkyl_mat_set(&A_uz,4,7,temp_rho); 
 
  temp_rho = 0.31943828249997*rho[8]+0.5*rho[5]; 
  gkyl_mat_set(&A_ux,4,8,temp_rho); 
  gkyl_mat_set(&A_uy,4,8,temp_rho); 
  gkyl_mat_set(&A_uz,4,8,temp_rho); 
 
  temp_rho = 0.5*rho[5]; 
  gkyl_mat_set(&A_ux,5,0,temp_rho); 
  gkyl_mat_set(&A_uy,5,0,temp_rho); 
  gkyl_mat_set(&A_uz,5,0,temp_rho); 
 
  temp_rho = 0.5000000000000001*rho[7]; 
  gkyl_mat_set(&A_ux,5,1,temp_rho); 
  gkyl_mat_set(&A_uy,5,1,temp_rho); 
  gkyl_mat_set(&A_uz,5,1,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[2]; 
  gkyl_mat_set(&A_ux,5,2,temp_rho); 
  gkyl_mat_set(&A_uy,5,2,temp_rho); 
  gkyl_mat_set(&A_uz,5,2,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[3]; 
  gkyl_mat_set(&A_ux,5,3,temp_rho); 
  gkyl_mat_set(&A_uy,5,3,temp_rho); 
  gkyl_mat_set(&A_uz,5,3,temp_rho); 
 
  temp_rho = 0.5*rho[8]; 
  gkyl_mat_set(&A_ux,5,4,temp_rho); 
  gkyl_mat_set(&A_uy,5,4,temp_rho); 
  gkyl_mat_set(&A_uz,5,4,temp_rho); 
 
  temp_rho = 0.31943828249997*rho[5]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,5,5,temp_rho); 
  gkyl_mat_set(&A_uy,5,5,temp_rho); 
  gkyl_mat_set(&A_uz,5,5,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[6]; 
  gkyl_mat_set(&A_ux,5,6,temp_rho); 
  gkyl_mat_set(&A_uy,5,6,temp_rho); 
  gkyl_mat_set(&A_uz,5,6,temp_rho); 
 
  temp_rho = 0.31943828249997*rho[7]+0.5000000000000001*rho[1]; 
  gkyl_mat_set(&A_ux,5,7,temp_rho); 
  gkyl_mat_set(&A_uy,5,7,temp_rho); 
  gkyl_mat_set(&A_uz,5,7,temp_rho); 
 
  temp_rho = 0.31943828249997*rho[8]+0.5*rho[4]; 
  gkyl_mat_set(&A_ux,5,8,temp_rho); 
  gkyl_mat_set(&A_uy,5,8,temp_rho); 
  gkyl_mat_set(&A_uz,5,8,temp_rho); 
 
  temp_rho = 0.5*rho[6]; 
  gkyl_mat_set(&A_ux,6,0,temp_rho); 
  gkyl_mat_set(&A_uy,6,0,temp_rho); 
  gkyl_mat_set(&A_uz,6,0,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[3]; 
  gkyl_mat_set(&A_ux,6,1,temp_rho); 
  gkyl_mat_set(&A_uy,6,1,temp_rho); 
  gkyl_mat_set(&A_uz,6,1,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[8]+0.5000000000000001*rho[4]; 
  gkyl_mat_set(&A_ux,6,2,temp_rho); 
  gkyl_mat_set(&A_uy,6,2,temp_rho); 
  gkyl_mat_set(&A_uz,6,2,temp_rho); 
 
  temp_rho = 0.4*rho[7]+0.447213595499958*rho[1]; 
  gkyl_mat_set(&A_ux,6,3,temp_rho); 
  gkyl_mat_set(&A_uy,6,3,temp_rho); 
  gkyl_mat_set(&A_uz,6,3,temp_rho); 
 
  temp_rho = 0.31943828249997*rho[6]+0.5000000000000001*rho[2]; 
  gkyl_mat_set(&A_ux,6,4,temp_rho); 
  gkyl_mat_set(&A_uy,6,4,temp_rho); 
  gkyl_mat_set(&A_uz,6,4,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[6]; 
  gkyl_mat_set(&A_ux,6,5,temp_rho); 
  gkyl_mat_set(&A_uy,6,5,temp_rho); 
  gkyl_mat_set(&A_uz,6,5,temp_rho); 
 
  temp_rho = 0.2857142857142857*rho[8]+0.4472135954999579*rho[5]+0.31943828249997*rho[4]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,6,6,temp_rho); 
  gkyl_mat_set(&A_uy,6,6,temp_rho); 
  gkyl_mat_set(&A_uz,6,6,temp_rho); 
 
  temp_rho = 0.4*rho[3]; 
  gkyl_mat_set(&A_ux,6,7,temp_rho); 
  gkyl_mat_set(&A_uy,6,7,temp_rho); 
  gkyl_mat_set(&A_uz,6,7,temp_rho); 
 
  temp_rho = 0.2857142857142857*rho[6]+0.447213595499958*rho[2]; 
  gkyl_mat_set(&A_ux,6,8,temp_rho); 
  gkyl_mat_set(&A_uy,6,8,temp_rho); 
  gkyl_mat_set(&A_uz,6,8,temp_rho); 
 
  temp_rho = 0.5*rho[7]; 
  gkyl_mat_set(&A_ux,7,0,temp_rho); 
  gkyl_mat_set(&A_uy,7,0,temp_rho); 
  gkyl_mat_set(&A_uz,7,0,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[8]+0.5000000000000001*rho[5]; 
  gkyl_mat_set(&A_ux,7,1,temp_rho); 
  gkyl_mat_set(&A_uy,7,1,temp_rho); 
  gkyl_mat_set(&A_uz,7,1,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[3]; 
  gkyl_mat_set(&A_ux,7,2,temp_rho); 
  gkyl_mat_set(&A_uy,7,2,temp_rho); 
  gkyl_mat_set(&A_uz,7,2,temp_rho); 
 
  temp_rho = 0.4*rho[6]+0.447213595499958*rho[2]; 
  gkyl_mat_set(&A_ux,7,3,temp_rho); 
  gkyl_mat_set(&A_uy,7,3,temp_rho); 
  gkyl_mat_set(&A_uz,7,3,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[7]; 
  gkyl_mat_set(&A_ux,7,4,temp_rho); 
  gkyl_mat_set(&A_uy,7,4,temp_rho); 
  gkyl_mat_set(&A_uz,7,4,temp_rho); 
 
  temp_rho = 0.31943828249997*rho[7]+0.5000000000000001*rho[1]; 
  gkyl_mat_set(&A_ux,7,5,temp_rho); 
  gkyl_mat_set(&A_uy,7,5,temp_rho); 
  gkyl_mat_set(&A_uz,7,5,temp_rho); 
 
  temp_rho = 0.4*rho[3]; 
  gkyl_mat_set(&A_ux,7,6,temp_rho); 
  gkyl_mat_set(&A_uy,7,6,temp_rho); 
  gkyl_mat_set(&A_uz,7,6,temp_rho); 
 
  temp_rho = 0.2857142857142857*rho[8]+0.31943828249997*rho[5]+0.4472135954999579*rho[4]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,7,7,temp_rho); 
  gkyl_mat_set(&A_uy,7,7,temp_rho); 
  gkyl_mat_set(&A_uz,7,7,temp_rho); 
 
  temp_rho = 0.2857142857142857*rho[7]+0.447213595499958*rho[1]; 
  gkyl_mat_set(&A_ux,7,8,temp_rho); 
  gkyl_mat_set(&A_uy,7,8,temp_rho); 
  gkyl_mat_set(&A_uz,7,8,temp_rho); 
 
  temp_rho = 0.5*rho[8]; 
  gkyl_mat_set(&A_ux,8,0,temp_rho); 
  gkyl_mat_set(&A_uy,8,0,temp_rho); 
  gkyl_mat_set(&A_uz,8,0,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[7]; 
  gkyl_mat_set(&A_ux,8,1,temp_rho); 
  gkyl_mat_set(&A_uy,8,1,temp_rho); 
  gkyl_mat_set(&A_uz,8,1,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[6]; 
  gkyl_mat_set(&A_ux,8,2,temp_rho); 
  gkyl_mat_set(&A_uy,8,2,temp_rho); 
  gkyl_mat_set(&A_uz,8,2,temp_rho); 
 
  temp_rho = 0.4*rho[3]; 
  gkyl_mat_set(&A_ux,8,3,temp_rho); 
  gkyl_mat_set(&A_uy,8,3,temp_rho); 
  gkyl_mat_set(&A_uz,8,3,temp_rho); 
 
  temp_rho = 0.31943828249997*rho[8]+0.5*rho[5]; 
  gkyl_mat_set(&A_ux,8,4,temp_rho); 
  gkyl_mat_set(&A_uy,8,4,temp_rho); 
  gkyl_mat_set(&A_uz,8,4,temp_rho); 
 
  temp_rho = 0.31943828249997*rho[8]+0.5*rho[4]; 
  gkyl_mat_set(&A_ux,8,5,temp_rho); 
  gkyl_mat_set(&A_uy,8,5,temp_rho); 
  gkyl_mat_set(&A_uz,8,5,temp_rho); 
 
  temp_rho = 0.2857142857142857*rho[6]+0.447213595499958*rho[2]; 
  gkyl_mat_set(&A_ux,8,6,temp_rho); 
  gkyl_mat_set(&A_uy,8,6,temp_rho); 
  gkyl_mat_set(&A_uz,8,6,temp_rho); 
 
  temp_rho = 0.2857142857142857*rho[7]+0.447213595499958*rho[1]; 
  gkyl_mat_set(&A_ux,8,7,temp_rho); 
  gkyl_mat_set(&A_uy,8,7,temp_rho); 
  gkyl_mat_set(&A_uz,8,7,temp_rho); 
 
  temp_rho = 0.2040816326530612*rho[8]+0.31943828249997*rho[5]+0.31943828249997*rho[4]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,8,8,temp_rho); 
  gkyl_mat_set(&A_uy,8,8,temp_rho); 
  gkyl_mat_set(&A_uz,8,8,temp_rho); 
 
  return cell_avg;
} 
