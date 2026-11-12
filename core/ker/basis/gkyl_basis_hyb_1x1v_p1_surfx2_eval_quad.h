GKYL_CU_DH static inline double 
hyb_1x1v_p1_surfx2_eval_quad_node_0_r(const double* GKYL_RESTRICT f) { 
  return -(1.1180339887498951*f[5])+1.118033988749895*f[4]-0.8660254037844386*f[3]+0.8660254037844386*f[2]-0.5*f[1]+0.5*f[0]; 
} 
GKYL_CU_DH static inline double 
hyb_1x1v_p1_surfx2_eval_quad_node_0_l(const double* GKYL_RESTRICT f) { 
    return -(1.1180339887498951*f[5])+1.118033988749895*f[4]+0.8660254037844386*f[3]-0.8660254037844386*f[2]-0.5*f[1]+0.5*f[0]; 
} 
GKYL_CU_DH static inline double 
hyb_1x1v_p1_surfx2_eval_quad_node_1_r(const double* GKYL_RESTRICT f) { 
  return 1.1180339887498951*f[5]+1.118033988749895*f[4]+0.8660254037844386*(f[3]+f[2])+0.5*(f[1]+f[0]); 
} 
GKYL_CU_DH static inline double 
hyb_1x1v_p1_surfx2_eval_quad_node_1_l(const double* GKYL_RESTRICT f) { 
    return 1.1180339887498951*f[5]+1.118033988749895*f[4]-0.8660254037844386*(f[3]+f[2])+0.5*(f[1]+f[0]); 
} 
