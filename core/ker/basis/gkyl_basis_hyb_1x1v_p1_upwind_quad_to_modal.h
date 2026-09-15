GKYL_CU_DH static inline void 
hyb_1x1v_p1_xdir_upwind_quad_to_modal(const double* fUpwindQuad, double* GKYL_RESTRICT fUpwind) { 
  fUpwind[0] = 0.39283710065919303*fUpwindQuad[2]+0.6285393610547091*fUpwindQuad[1]+0.39283710065919303*fUpwindQuad[0]; 
  fUpwind[1] = 0.5270462766947298*fUpwindQuad[2]-0.5270462766947298*fUpwindQuad[0]; 

} 

GKYL_CU_DH static inline void 
hyb_1x1v_p1_vdir_upwind_quad_to_modal(const double* fUpwindQuad, double* GKYL_RESTRICT fUpwind) { 
  fUpwind[0] = 0.7071067811865475*fUpwindQuad[1]+0.7071067811865475*fUpwindQuad[0]; 
  fUpwind[1] = 0.7071067811865475*fUpwindQuad[1]-0.7071067811865475*fUpwindQuad[0]; 

} 

