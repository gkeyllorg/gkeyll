#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vy_ser_p1_lovy(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf) {
  // drag_coeff_surf: Surface projection of drag coefficient at lower boundary.
  // sgn_drag_coeff_surf: sign(drag_coeff_surf) evaluated at quadrature points along lower surface.
  // const_sgn_drag_coeff_surf: 1 if sign(drag_coeff_surf) is constant along lower boundary, 0 otherwise. 

  const double *alpha_surf = &drag_coeff_surf[8]; 
  double *sgn_alpha_surf = &sgn_drag_coeff_surf[8]; 

  sgn_alpha_surf[0] = 0.0;
  sgn_alpha_surf[1] = 0.0;
  sgn_alpha_surf[2] = 0.0;
  sgn_alpha_surf[3] = 0.0;
  sgn_alpha_surf[4] = 0.0;
  sgn_alpha_surf[5] = 0.0;
  sgn_alpha_surf[6] = 0.0;
  sgn_alpha_surf[7] = 0.0;
  int const_sgn_alpha_surf = 0; 
  *const_sgn_drag_coeff_surf = const_sgn_alpha_surf; 
} 
