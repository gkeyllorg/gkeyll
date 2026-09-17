#include <gkyl_basis_gkhyb_vel_kernels.h> 
GKYL_CU_DH 
void 
quad_to_modal_1d_gkhyb_vel_p1(const double* fquad, double* GKYL_RESTRICT fmodal, long linc2) { 
  switch (linc2) { 
    case 0: 
      fmodal[0] = 0.39283710065919303*fquad[2]+0.6285393610547091*fquad[1]+0.39283710065919303*fquad[0]; 
    break; 

    case 1: 
      fmodal[1] = 0.5270462766947298*fquad[2]-0.5270462766947298*fquad[0]; 
    break; 

  } 
} 

GKYL_CU_DH 
void 
quad_to_modal_2d_gkhyb_vel_p1(const double* fquad, double* GKYL_RESTRICT fmodal, long linc2) { 
  switch (linc2) { 
    case 0: 
      fmodal[0] = 0.2777777777777778*fquad[5]+0.4444444444444444*fquad[4]+0.2777777777777778*fquad[3]+0.2777777777777778*fquad[2]+0.4444444444444444*fquad[1]+0.2777777777777778*fquad[0]; 
    break; 

    case 1: 
      fmodal[1] = 0.37267799624996495*fquad[5]-0.37267799624996495*fquad[3]+0.37267799624996495*fquad[2]-0.37267799624996495*fquad[0]; 
    break; 

    case 2: 
      fmodal[2] = 0.2777777777777778*fquad[5]+0.4444444444444444*fquad[4]+0.2777777777777778*fquad[3]-0.2777777777777778*fquad[2]-0.4444444444444444*fquad[1]-0.2777777777777778*fquad[0]; 
    break; 

    case 3: 
      fmodal[3] = 0.37267799624996495*fquad[5]-0.37267799624996495*fquad[3]-0.37267799624996495*fquad[2]+0.37267799624996495*fquad[0]; 
    break; 

  } 
} 

