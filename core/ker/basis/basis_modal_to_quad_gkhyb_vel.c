#include <gkyl_basis_gkhyb_vel_kernels.h> 
GKYL_CU_DH 
void 
modal_to_quad_1d_gkhyb_vel_p1(const double* fmodal, double* GKYL_RESTRICT fquad, long linc2) { 
  switch (linc2) { 
    case 0: 
      fquad[0] = 0.7071067811865475*fmodal[0]-0.9486832980505137*fmodal[1]; 
    break; 

    case 1: 
      fquad[1] = 0.7071067811865475*fmodal[0]; 
    break; 

    case 2: 
      fquad[2] = 0.9486832980505137*fmodal[1]+0.7071067811865475*fmodal[0]; 
    break; 

  } 
} 

GKYL_CU_DH 
void 
modal_to_quad_2d_gkhyb_vel_p1(const double* fmodal, double* GKYL_RESTRICT fquad, long linc2) { 
  switch (linc2) { 
    case 0: 
      fquad[0] = 0.6708203932499369*fmodal[3]-0.5*fmodal[2]-0.6708203932499369*fmodal[1]+0.5*fmodal[0]; 
    break; 

    case 1: 
      fquad[1] = 0.5*fmodal[0]-0.5*fmodal[2]; 
    break; 

    case 2: 
      fquad[2] = -(0.6708203932499369*fmodal[3])-0.5*fmodal[2]+0.6708203932499369*fmodal[1]+0.5*fmodal[0]; 
    break; 

    case 3: 
      fquad[3] = -(0.6708203932499369*fmodal[3])+0.5*fmodal[2]-0.6708203932499369*fmodal[1]+0.5*fmodal[0]; 
    break; 

    case 4: 
      fquad[4] = 0.5*fmodal[2]+0.5*fmodal[0]; 
    break; 

    case 5: 
      fquad[5] = 0.6708203932499369*fmodal[3]+0.5*fmodal[2]+0.6708203932499369*fmodal[1]+0.5*fmodal[0]; 
    break; 

  } 
} 

