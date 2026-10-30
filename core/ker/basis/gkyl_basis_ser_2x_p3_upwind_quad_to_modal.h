GKYL_CU_DH static inline void 
ser_2x_p3_upwind_quad_to_modal(const double* fUpwindQuad, double* GKYL_RESTRICT fUpwind) { 
  fUpwind[0] = 0.24597051986529*fUpwindQuad[3]+0.46113626132125757*fUpwindQuad[2]+0.46113626132125757*fUpwindQuad[1]+0.24597051986529*fUpwindQuad[0]; 
  fUpwind[1] = 0.3668728630454641*fUpwindQuad[3]+0.271546746793545*fUpwindQuad[2]-0.271546746793545*fUpwindQuad[1]-0.3668728630454641*fUpwindQuad[0]; 
  fUpwind[2] = 0.3367876570272824*fUpwindQuad[3]-0.3367876570272824*fUpwindQuad[2]-0.3367876570272824*fUpwindQuad[1]+0.3367876570272824*fUpwindQuad[0]; 
  fUpwind[3] = 0.19832227542449957*fUpwindQuad[3]-0.5023295150965308*fUpwindQuad[2]+0.5023295150965308*fUpwindQuad[1]-0.19832227542449957*fUpwindQuad[0]; 

} 
