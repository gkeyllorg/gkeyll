#include <gkyl_mat.h> 
#include <gkyl_euler_kernels.h> 
GKYL_CU_DH int fluid_vars_u_set_1x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
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
  const double *rhoux = &fluid[4]; 
  const double *rhouy = &fluid[8]; 
  const double *rhouz = &fluid[12]; 

  int cell_avg = 0;
  // Check if rho < 0 at control points. 
  if ((-1.870828693386971*rho[3])+1.58113883008419*rho[2]-1.224744871391589*rho[1]+0.7071067811865475*rho[0] < 0.0) cell_avg = 1; 
  if (0.7621894676761731*rho[3]-0.5270462766947298*rho[2]-0.408248290463863*rho[1]+0.7071067811865475*rho[0] < 0.0) cell_avg = 1; 
  if ((-0.7621894676761731*rho[3])-0.5270462766947298*rho[2]+0.408248290463863*rho[1]+0.7071067811865475*rho[0] < 0.0) cell_avg = 1; 
  if (1.870828693386971*rho[3]+1.58113883008419*rho[2]+1.224744871391589*rho[1]+0.7071067811865475*rho[0] < 0.0) cell_avg = 1; 
 
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
 
  double temp_rho = 0.0; 
  temp_rho = 0.7071067811865475*rho[0]; 
  gkyl_mat_set(&A_ux,0,0,temp_rho); 
  gkyl_mat_set(&A_uy,0,0,temp_rho); 
  gkyl_mat_set(&A_uz,0,0,temp_rho); 
 
  temp_rho = 0.7071067811865475*rho[1]; 
  gkyl_mat_set(&A_ux,0,1,temp_rho); 
  gkyl_mat_set(&A_uy,0,1,temp_rho); 
  gkyl_mat_set(&A_uz,0,1,temp_rho); 
 
  temp_rho = 0.7071067811865475*rho[2]; 
  gkyl_mat_set(&A_ux,0,2,temp_rho); 
  gkyl_mat_set(&A_uy,0,2,temp_rho); 
  gkyl_mat_set(&A_uz,0,2,temp_rho); 
 
  temp_rho = 0.7071067811865475*rho[3]; 
  gkyl_mat_set(&A_ux,0,3,temp_rho); 
  gkyl_mat_set(&A_uy,0,3,temp_rho); 
  gkyl_mat_set(&A_uz,0,3,temp_rho); 
 
  temp_rho = 0.7071067811865475*rho[1]; 
  gkyl_mat_set(&A_ux,1,0,temp_rho); 
  gkyl_mat_set(&A_uy,1,0,temp_rho); 
  gkyl_mat_set(&A_uz,1,0,temp_rho); 
 
  temp_rho = 0.6324555320336759*rho[2]+0.7071067811865475*rho[0]; 
  gkyl_mat_set(&A_ux,1,1,temp_rho); 
  gkyl_mat_set(&A_uy,1,1,temp_rho); 
  gkyl_mat_set(&A_uz,1,1,temp_rho); 
 
  temp_rho = 0.6210590034081186*rho[3]+0.6324555320336759*rho[1]; 
  gkyl_mat_set(&A_ux,1,2,temp_rho); 
  gkyl_mat_set(&A_uy,1,2,temp_rho); 
  gkyl_mat_set(&A_uz,1,2,temp_rho); 
 
  temp_rho = 0.6210590034081186*rho[2]; 
  gkyl_mat_set(&A_ux,1,3,temp_rho); 
  gkyl_mat_set(&A_uy,1,3,temp_rho); 
  gkyl_mat_set(&A_uz,1,3,temp_rho); 
 
  temp_rho = 0.7071067811865475*rho[2]; 
  gkyl_mat_set(&A_ux,2,0,temp_rho); 
  gkyl_mat_set(&A_uy,2,0,temp_rho); 
  gkyl_mat_set(&A_uz,2,0,temp_rho); 
 
  temp_rho = 0.6210590034081186*rho[3]+0.6324555320336759*rho[1]; 
  gkyl_mat_set(&A_ux,2,1,temp_rho); 
  gkyl_mat_set(&A_uy,2,1,temp_rho); 
  gkyl_mat_set(&A_uz,2,1,temp_rho); 
 
  temp_rho = 0.4517539514526256*rho[2]+0.7071067811865475*rho[0]; 
  gkyl_mat_set(&A_ux,2,2,temp_rho); 
  gkyl_mat_set(&A_uy,2,2,temp_rho); 
  gkyl_mat_set(&A_uz,2,2,temp_rho); 
 
  temp_rho = 0.421637021355784*rho[3]+0.6210590034081186*rho[1]; 
  gkyl_mat_set(&A_ux,2,3,temp_rho); 
  gkyl_mat_set(&A_uy,2,3,temp_rho); 
  gkyl_mat_set(&A_uz,2,3,temp_rho); 
 
  temp_rho = 0.7071067811865475*rho[3]; 
  gkyl_mat_set(&A_ux,3,0,temp_rho); 
  gkyl_mat_set(&A_uy,3,0,temp_rho); 
  gkyl_mat_set(&A_uz,3,0,temp_rho); 
 
  temp_rho = 0.6210590034081186*rho[2]; 
  gkyl_mat_set(&A_ux,3,1,temp_rho); 
  gkyl_mat_set(&A_uy,3,1,temp_rho); 
  gkyl_mat_set(&A_uz,3,1,temp_rho); 
 
  temp_rho = 0.421637021355784*rho[3]+0.6210590034081186*rho[1]; 
  gkyl_mat_set(&A_ux,3,2,temp_rho); 
  gkyl_mat_set(&A_uy,3,2,temp_rho); 
  gkyl_mat_set(&A_uz,3,2,temp_rho); 
 
  temp_rho = 0.421637021355784*rho[2]+0.7071067811865475*rho[0]; 
  gkyl_mat_set(&A_ux,3,3,temp_rho); 
  gkyl_mat_set(&A_uy,3,3,temp_rho); 
  gkyl_mat_set(&A_uz,3,3,temp_rho); 
 
  return cell_avg;
} 
