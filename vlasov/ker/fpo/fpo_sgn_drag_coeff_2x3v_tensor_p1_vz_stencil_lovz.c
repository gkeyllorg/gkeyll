#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_sgn_drag_coeff_2x3v_vz_tensor_p1_lovz(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf) {
  // drag_coeff_surf: Surface projection of drag coefficient at lower boundary.
  // sgn_drag_coeff_surf: sign(drag_coeff_surf) evaluated at quadrature points along lower surface.
  // const_sgn_drag_coeff_surf: 1 if sign(drag_coeff_surf) is constant along lower boundary, 0 otherwise. 

  const double *alpha_surf = &drag_coeff_surf[72]; 
  double *sgn_alpha_surf = &sgn_drag_coeff_surf[72]; 

  sgn_alpha_surf[0] = 0.0;
  sgn_alpha_surf[1] = 0.0;
  sgn_alpha_surf[2] = 0.0;
  sgn_alpha_surf[3] = 0.0;
  sgn_alpha_surf[4] = 0.0;
  sgn_alpha_surf[5] = 0.0;
  sgn_alpha_surf[6] = 0.0;
  sgn_alpha_surf[7] = 0.0;
  sgn_alpha_surf[8] = 0.0;
  sgn_alpha_surf[9] = 0.0;
  sgn_alpha_surf[10] = 0.0;
  sgn_alpha_surf[11] = 0.0;
  sgn_alpha_surf[12] = 0.0;
  sgn_alpha_surf[13] = 0.0;
  sgn_alpha_surf[14] = 0.0;
  sgn_alpha_surf[15] = 0.0;
  sgn_alpha_surf[16] = 0.0;
  sgn_alpha_surf[17] = 0.0;
  sgn_alpha_surf[18] = 0.0;
  sgn_alpha_surf[19] = 0.0;
  sgn_alpha_surf[20] = 0.0;
  sgn_alpha_surf[21] = 0.0;
  sgn_alpha_surf[22] = 0.0;
  sgn_alpha_surf[23] = 0.0;
  sgn_alpha_surf[24] = 0.0;
  sgn_alpha_surf[25] = 0.0;
  sgn_alpha_surf[26] = 0.0;
  sgn_alpha_surf[27] = 0.0;
  sgn_alpha_surf[28] = 0.0;
  sgn_alpha_surf[29] = 0.0;
  sgn_alpha_surf[30] = 0.0;
  sgn_alpha_surf[31] = 0.0;
  sgn_alpha_surf[32] = 0.0;
  sgn_alpha_surf[33] = 0.0;
  sgn_alpha_surf[34] = 0.0;
  sgn_alpha_surf[35] = 0.0;
  int const_sgn_alpha_surf = 0; 
  *const_sgn_drag_coeff_surf = const_sgn_alpha_surf; 
} 
