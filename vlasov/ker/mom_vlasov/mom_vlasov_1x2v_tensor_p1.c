#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vx2_sq = vmap_vx[2]*vmap_vx[2]; 
  const double vmap_vx3_sq = vmap_vx[3]*vmap_vx[3]; 
  const double vmap_vy0_sq = vmap_vy[0]*vmap_vy[0]; 
  const double vmap_vy1_sq = vmap_vy[1]*vmap_vy[1]; 
  const double vmap_vy2_sq = vmap_vy[2]*vmap_vy[2]; 
  const double vmap_vy3_sq = vmap_vy[3]*vmap_vy[3]; 
  out[0] += (0.5962847939999438*f[7]*vmap_vx3_sq+f[0]*vmap_vx3_sq+0.6388765649999399*f[7]*vmap_vx2_sq+f[0]*vmap_vx2_sq+0.8944271909999159*f[7]*vmap_vx1_sq+f[0]*vmap_vx1_sq+f[0]*vmap_vx0_sq+1.7566201313073593*vmap_vx[1]*vmap_vx[3]*f[7]+2.0*vmap_vx[0]*vmap_vx[2]*f[7]+1.7566201313073593*f[2]*vmap_vx[2]*vmap_vx[3]+1.7888543819998317*vmap_vx[1]*f[2]*vmap_vx[2]+2.0*vmap_vx[0]*vmap_vx[1]*f[2])*volFact; 
  out[1] += (0.5962847939999442*f[10]*vmap_vx3_sq+f[1]*vmap_vx3_sq+0.63887656499994*f[10]*vmap_vx2_sq+f[1]*vmap_vx2_sq+0.8944271909999161*f[10]*vmap_vx1_sq+f[1]*vmap_vx1_sq+f[1]*vmap_vx0_sq+1.7566201313073595*vmap_vx[1]*vmap_vx[3]*f[10]+2.0000000000000004*vmap_vx[0]*vmap_vx[2]*f[10]+1.7566201313073593*vmap_vx[2]*vmap_vx[3]*f[4]+1.7888543819998317*vmap_vx[1]*vmap_vx[2]*f[4]+2.0*vmap_vx[0]*vmap_vx[1]*f[4])*volFact; 
  out[2] += (vmap_vx[2]*vmap_vy[2]*f[16]+1.0000000000000002*vmap_vx[1]*vmap_vy[2]*f[13]+1.0000000000000002*vmap_vy[1]*vmap_vx[2]*f[11]+vmap_vx[0]*vmap_vy[2]*f[8]+vmap_vy[0]*vmap_vx[2]*f[7]+vmap_vx[1]*vmap_vy[1]*f[6]+vmap_vx[0]*vmap_vy[1]*f[3]+vmap_vy[0]*vmap_vx[1]*f[2]+f[0]*vmap_vx[0]*vmap_vy[0])*volFact; 
  out[3] += (vmap_vx[2]*vmap_vy[2]*f[17]+vmap_vx[1]*vmap_vy[2]*f[15]+vmap_vy[1]*vmap_vx[2]*f[14]+1.0000000000000002*vmap_vx[0]*vmap_vy[2]*f[12]+1.0000000000000002*vmap_vy[0]*vmap_vx[2]*f[10]+vmap_vx[1]*vmap_vy[1]*f[9]+vmap_vx[0]*vmap_vy[1]*f[5]+vmap_vy[0]*vmap_vx[1]*f[4]+vmap_vx[0]*vmap_vy[0]*f[1])*volFact; 
  out[4] += (0.5962847939999438*f[8]*vmap_vy3_sq+f[0]*vmap_vy3_sq+0.6388765649999399*f[8]*vmap_vy2_sq+f[0]*vmap_vy2_sq+0.8944271909999159*f[8]*vmap_vy1_sq+f[0]*vmap_vy1_sq+f[0]*vmap_vy0_sq+1.7566201313073593*vmap_vy[1]*vmap_vy[3]*f[8]+2.0*vmap_vy[0]*vmap_vy[2]*f[8]+1.7566201313073593*vmap_vy[2]*f[3]*vmap_vy[3]+1.7888543819998317*vmap_vy[1]*vmap_vy[2]*f[3]+2.0*vmap_vy[0]*vmap_vy[1]*f[3])*volFact; 
  out[5] += (0.5962847939999442*f[12]*vmap_vy3_sq+f[1]*vmap_vy3_sq+0.63887656499994*f[12]*vmap_vy2_sq+f[1]*vmap_vy2_sq+0.8944271909999161*f[12]*vmap_vy1_sq+f[1]*vmap_vy1_sq+f[1]*vmap_vy0_sq+1.7566201313073595*vmap_vy[1]*vmap_vy[3]*f[12]+2.0000000000000004*vmap_vy[0]*vmap_vy[2]*f[12]+1.7566201313073593*vmap_vy[2]*vmap_vy[3]*f[5]+1.7888543819998317*vmap_vy[1]*vmap_vy[2]*f[5]+2.0*vmap_vy[0]*vmap_vy[1]*f[5])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx0_cu = vmap_vx[0]*vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vx1_cu = vmap_vx[1]*vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vx2_sq = vmap_vx[2]*vmap_vx[2]; 
  const double vmap_vx2_cu = vmap_vx[2]*vmap_vx[2]*vmap_vx[2]; 
  const double vmap_vx3_sq = vmap_vx[3]*vmap_vx[3]; 
  const double vmap_vx3_cu = vmap_vx[3]*vmap_vx[3]*vmap_vx[3]; 
  const double vmap_vy0_sq = vmap_vy[0]*vmap_vy[0]; 
  const double vmap_vy0_cu = vmap_vy[0]*vmap_vy[0]*vmap_vy[0]; 
  const double vmap_vy1_sq = vmap_vy[1]*vmap_vy[1]; 
  const double vmap_vy1_cu = vmap_vy[1]*vmap_vy[1]*vmap_vy[1]; 
  const double vmap_vy2_sq = vmap_vy[2]*vmap_vy[2]; 
  const double vmap_vy2_cu = vmap_vy[2]*vmap_vy[2]*vmap_vy[2]; 
  const double vmap_vy3_sq = vmap_vy[3]*vmap_vy[3]; 
  const double vmap_vy3_cu = vmap_vy[3]*vmap_vy[3]*vmap_vy[3]; 
  out[0] += (3.9212285138526726*vmap_vx[2]*f[7]*vmap_vx3_sq+1.264911064067352*vmap_vx[0]*f[7]*vmap_vx3_sq+1.264911064067352*f[0]*vmap_vx[2]*vmap_vx3_sq+3.252691193458118*vmap_vx[1]*f[2]*vmap_vx3_sq+2.1213203435596424*f[0]*vmap_vx[0]*vmap_vx3_sq+0.7069898943717666*f[2]*vmap_vx3_cu+1.355261854357877*vmap_vx[0]*f[7]*vmap_vx2_sq+2.777460299317654*f[2]*vmap_vx[3]*vmap_vx2_sq+3.3335033970222954*vmap_vx[1]*f[2]*vmap_vx2_sq+2.1213203435596424*f[0]*vmap_vx[0]*vmap_vx2_sq+1.515228816828316*f[7]*vmap_vx2_cu+0.45175395145262565*f[0]*vmap_vx2_cu+3.3335033970222954*vmap_vx[2]*f[7]*vmap_vx1_sq+1.8973665961010278*vmap_vx[0]*f[7]*vmap_vx1_sq+1.6664761795905922*f[2]*vmap_vx[3]*vmap_vx1_sq+1.8973665961010278*f[0]*vmap_vx[2]*vmap_vx1_sq+2.1213203435596424*f[0]*vmap_vx[0]*vmap_vx1_sq+1.2727922061357855*f[2]*vmap_vx1_cu+2.1213203435596424*vmap_vx[2]*f[7]*vmap_vx0_sq+2.1213203435596424*vmap_vx[1]*f[2]*vmap_vx0_sq+0.7071067811865475*f[0]*vmap_vx0_cu+5.5549205986353085*vmap_vx[1]*vmap_vx[2]*vmap_vx[3]*f[7]+3.7263540204487113*vmap_vx[0]*vmap_vx[1]*vmap_vx[3]*f[7]+3.7263540204487113*vmap_vx[0]*f[2]*vmap_vx[2]*vmap_vx[3]+3.7263540204487113*f[0]*vmap_vx[1]*vmap_vx[2]*vmap_vx[3]+3.7947331922020555*vmap_vx[0]*vmap_vx[1]*f[2]*vmap_vx[2])*volFact; 
  out[1] += (3.921228513852674*vmap_vx[2]*f[10]*vmap_vx3_sq+1.264911064067352*vmap_vx[0]*f[10]*vmap_vx3_sq+3.252691193458118*vmap_vx[1]*f[4]*vmap_vx3_sq+1.264911064067352*f[1]*vmap_vx[2]*vmap_vx3_sq+2.1213203435596424*vmap_vx[0]*f[1]*vmap_vx3_sq+0.7069898943717666*f[4]*vmap_vx3_cu+1.355261854357877*vmap_vx[0]*f[10]*vmap_vx2_sq+2.777460299317654*vmap_vx[3]*f[4]*vmap_vx2_sq+3.3335033970222954*vmap_vx[1]*f[4]*vmap_vx2_sq+2.1213203435596424*vmap_vx[0]*f[1]*vmap_vx2_sq+1.515228816828316*f[10]*vmap_vx2_cu+0.45175395145262565*f[1]*vmap_vx2_cu+3.3335033970222954*vmap_vx[2]*f[10]*vmap_vx1_sq+1.897366596101028*vmap_vx[0]*f[10]*vmap_vx1_sq+1.6664761795905922*vmap_vx[3]*f[4]*vmap_vx1_sq+1.8973665961010278*f[1]*vmap_vx[2]*vmap_vx1_sq+2.1213203435596424*vmap_vx[0]*f[1]*vmap_vx1_sq+1.2727922061357855*f[4]*vmap_vx1_cu+2.1213203435596424*vmap_vx[2]*f[10]*vmap_vx0_sq+2.1213203435596424*vmap_vx[1]*f[4]*vmap_vx0_sq+0.7071067811865475*f[1]*vmap_vx0_cu+5.554920598635309*vmap_vx[1]*vmap_vx[2]*vmap_vx[3]*f[10]+3.726354020448713*vmap_vx[0]*vmap_vx[1]*vmap_vx[3]*f[10]+3.7263540204487113*vmap_vx[0]*vmap_vx[2]*vmap_vx[3]*f[4]+3.7947331922020555*vmap_vx[0]*vmap_vx[1]*vmap_vx[2]*f[4]+3.7263540204487113*f[1]*vmap_vx[1]*vmap_vx[2]*vmap_vx[3])*volFact; 
  out[2] += (0.421637021355784*vmap_vy[2]*f[16]*vmap_vx3_sq+0.4216370213557841*vmap_vy[1]*f[11]*vmap_vx3_sq+0.7071067811865475*vmap_vy[2]*f[8]*vmap_vx3_sq+0.421637021355784*vmap_vy[0]*f[7]*vmap_vx3_sq+0.7071067811865475*vmap_vy[1]*f[3]*vmap_vx3_sq+0.7071067811865475*f[0]*vmap_vy[0]*vmap_vx3_sq+0.45175395145262565*vmap_vy[2]*f[16]*vmap_vx2_sq+0.4517539514526257*vmap_vy[1]*f[11]*vmap_vx2_sq+0.7071067811865475*vmap_vy[2]*f[8]*vmap_vx2_sq+0.45175395145262565*vmap_vy[0]*f[7]*vmap_vx2_sq+0.7071067811865475*vmap_vy[1]*f[3]*vmap_vx2_sq+0.7071067811865475*f[0]*vmap_vy[0]*vmap_vx2_sq+0.6324555320336759*vmap_vy[2]*f[16]*vmap_vx1_sq+0.632455532033676*vmap_vy[1]*f[11]*vmap_vx1_sq+0.7071067811865475*vmap_vy[2]*f[8]*vmap_vx1_sq+0.6324555320336759*vmap_vy[0]*f[7]*vmap_vx1_sq+0.7071067811865475*vmap_vy[1]*f[3]*vmap_vx1_sq+0.7071067811865475*f[0]*vmap_vy[0]*vmap_vx1_sq+0.7071067811865475*vmap_vy[2]*f[8]*vmap_vx0_sq+0.7071067811865475*vmap_vy[1]*f[3]*vmap_vx0_sq+0.7071067811865475*f[0]*vmap_vy[0]*vmap_vx0_sq+1.2421180068162374*vmap_vx[1]*vmap_vy[2]*vmap_vx[3]*f[16]+1.4142135623730951*vmap_vx[0]*vmap_vx[2]*vmap_vy[2]*f[16]+1.2421180068162376*vmap_vx[2]*vmap_vy[2]*vmap_vx[3]*f[13]+1.264911064067352*vmap_vx[1]*vmap_vx[2]*vmap_vy[2]*f[13]+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*vmap_vy[2]*f[13]+1.2421180068162376*vmap_vx[1]*vmap_vy[1]*vmap_vx[3]*f[11]+1.4142135623730951*vmap_vx[0]*vmap_vy[1]*vmap_vx[2]*f[11]+1.2421180068162374*vmap_vy[0]*vmap_vx[1]*vmap_vx[3]*f[7]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vx[2]*f[7]+1.2421180068162374*vmap_vy[1]*vmap_vx[2]*vmap_vx[3]*f[6]+1.264911064067352*vmap_vx[1]*vmap_vy[1]*vmap_vx[2]*f[6]+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[6]+1.2421180068162374*vmap_vy[0]*f[2]*vmap_vx[2]*vmap_vx[3]+1.264911064067352*vmap_vy[0]*vmap_vx[1]*f[2]*vmap_vx[2]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[2])*volFact; 
  out[3] += (0.421637021355784*vmap_vy[2]*f[17]*vmap_vx3_sq+0.421637021355784*vmap_vy[1]*f[14]*vmap_vx3_sq+0.7071067811865475*vmap_vy[2]*f[12]*vmap_vx3_sq+0.4216370213557841*vmap_vy[0]*f[10]*vmap_vx3_sq+0.7071067811865475*vmap_vy[1]*f[5]*vmap_vx3_sq+0.7071067811865475*vmap_vy[0]*f[1]*vmap_vx3_sq+0.45175395145262565*vmap_vy[2]*f[17]*vmap_vx2_sq+0.45175395145262565*vmap_vy[1]*f[14]*vmap_vx2_sq+0.7071067811865475*vmap_vy[2]*f[12]*vmap_vx2_sq+0.4517539514526257*vmap_vy[0]*f[10]*vmap_vx2_sq+0.7071067811865475*vmap_vy[1]*f[5]*vmap_vx2_sq+0.7071067811865475*vmap_vy[0]*f[1]*vmap_vx2_sq+0.6324555320336759*vmap_vy[2]*f[17]*vmap_vx1_sq+0.6324555320336759*vmap_vy[1]*f[14]*vmap_vx1_sq+0.7071067811865475*vmap_vy[2]*f[12]*vmap_vx1_sq+0.632455532033676*vmap_vy[0]*f[10]*vmap_vx1_sq+0.7071067811865475*vmap_vy[1]*f[5]*vmap_vx1_sq+0.7071067811865475*vmap_vy[0]*f[1]*vmap_vx1_sq+0.7071067811865475*vmap_vy[2]*f[12]*vmap_vx0_sq+0.7071067811865475*vmap_vy[1]*f[5]*vmap_vx0_sq+0.7071067811865475*vmap_vy[0]*f[1]*vmap_vx0_sq+1.2421180068162374*vmap_vx[1]*vmap_vy[2]*vmap_vx[3]*f[17]+1.4142135623730951*vmap_vx[0]*vmap_vx[2]*vmap_vy[2]*f[17]+1.2421180068162374*vmap_vx[2]*vmap_vy[2]*vmap_vx[3]*f[15]+1.264911064067352*vmap_vx[1]*vmap_vx[2]*vmap_vy[2]*f[15]+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*vmap_vy[2]*f[15]+1.2421180068162374*vmap_vx[1]*vmap_vy[1]*vmap_vx[3]*f[14]+1.4142135623730951*vmap_vx[0]*vmap_vy[1]*vmap_vx[2]*f[14]+1.2421180068162376*vmap_vy[0]*vmap_vx[1]*vmap_vx[3]*f[10]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vx[2]*f[10]+1.2421180068162374*vmap_vy[1]*vmap_vx[2]*vmap_vx[3]*f[9]+1.264911064067352*vmap_vx[1]*vmap_vy[1]*vmap_vx[2]*f[9]+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*vmap_vy[1]*f[9]+1.2421180068162374*vmap_vy[0]*vmap_vx[2]*vmap_vx[3]*f[4]+1.264911064067352*vmap_vy[0]*vmap_vx[1]*vmap_vx[2]*f[4]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vx[1]*f[4])*volFact; 
  out[4] += (0.421637021355784*vmap_vx[2]*f[16]*vmap_vy3_sq+0.4216370213557841*vmap_vx[1]*f[13]*vmap_vy3_sq+0.421637021355784*vmap_vx[0]*f[8]*vmap_vy3_sq+0.7071067811865475*vmap_vx[2]*f[7]*vmap_vy3_sq+0.7071067811865475*vmap_vx[1]*f[2]*vmap_vy3_sq+0.7071067811865475*f[0]*vmap_vx[0]*vmap_vy3_sq+0.45175395145262565*vmap_vx[2]*f[16]*vmap_vy2_sq+0.4517539514526257*vmap_vx[1]*f[13]*vmap_vy2_sq+0.45175395145262565*vmap_vx[0]*f[8]*vmap_vy2_sq+0.7071067811865475*vmap_vx[2]*f[7]*vmap_vy2_sq+0.7071067811865475*vmap_vx[1]*f[2]*vmap_vy2_sq+0.7071067811865475*f[0]*vmap_vx[0]*vmap_vy2_sq+0.6324555320336759*vmap_vx[2]*f[16]*vmap_vy1_sq+0.632455532033676*vmap_vx[1]*f[13]*vmap_vy1_sq+0.6324555320336759*vmap_vx[0]*f[8]*vmap_vy1_sq+0.7071067811865475*vmap_vx[2]*f[7]*vmap_vy1_sq+0.7071067811865475*vmap_vx[1]*f[2]*vmap_vy1_sq+0.7071067811865475*f[0]*vmap_vx[0]*vmap_vy1_sq+0.7071067811865475*vmap_vx[2]*f[7]*vmap_vy0_sq+0.7071067811865475*vmap_vx[1]*f[2]*vmap_vy0_sq+0.7071067811865475*f[0]*vmap_vx[0]*vmap_vy0_sq+1.2421180068162374*vmap_vy[1]*vmap_vx[2]*vmap_vy[3]*f[16]+1.4142135623730951*vmap_vy[0]*vmap_vx[2]*vmap_vy[2]*f[16]+1.2421180068162376*vmap_vx[1]*vmap_vy[1]*vmap_vy[3]*f[13]+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*vmap_vy[2]*f[13]+1.2421180068162376*vmap_vx[2]*vmap_vy[2]*vmap_vy[3]*f[11]+1.264911064067352*vmap_vy[1]*vmap_vx[2]*vmap_vy[2]*f[11]+1.4142135623730951*vmap_vy[0]*vmap_vy[1]*vmap_vx[2]*f[11]+1.2421180068162374*vmap_vx[0]*vmap_vy[1]*vmap_vy[3]*f[8]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vy[2]*f[8]+1.2421180068162374*vmap_vx[1]*vmap_vy[2]*vmap_vy[3]*f[6]+1.264911064067352*vmap_vx[1]*vmap_vy[1]*vmap_vy[2]*f[6]+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[6]+1.2421180068162374*vmap_vx[0]*vmap_vy[2]*f[3]*vmap_vy[3]+1.264911064067352*vmap_vx[0]*vmap_vy[1]*vmap_vy[2]*f[3]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[3])*volFact; 
  out[5] += (0.421637021355784*vmap_vx[2]*f[17]*vmap_vy3_sq+0.421637021355784*vmap_vx[1]*f[15]*vmap_vy3_sq+0.4216370213557841*vmap_vx[0]*f[12]*vmap_vy3_sq+0.7071067811865475*vmap_vx[2]*f[10]*vmap_vy3_sq+0.7071067811865475*vmap_vx[1]*f[4]*vmap_vy3_sq+0.7071067811865475*vmap_vx[0]*f[1]*vmap_vy3_sq+0.45175395145262565*vmap_vx[2]*f[17]*vmap_vy2_sq+0.45175395145262565*vmap_vx[1]*f[15]*vmap_vy2_sq+0.4517539514526257*vmap_vx[0]*f[12]*vmap_vy2_sq+0.7071067811865475*vmap_vx[2]*f[10]*vmap_vy2_sq+0.7071067811865475*vmap_vx[1]*f[4]*vmap_vy2_sq+0.7071067811865475*vmap_vx[0]*f[1]*vmap_vy2_sq+0.6324555320336759*vmap_vx[2]*f[17]*vmap_vy1_sq+0.6324555320336759*vmap_vx[1]*f[15]*vmap_vy1_sq+0.632455532033676*vmap_vx[0]*f[12]*vmap_vy1_sq+0.7071067811865475*vmap_vx[2]*f[10]*vmap_vy1_sq+0.7071067811865475*vmap_vx[1]*f[4]*vmap_vy1_sq+0.7071067811865475*vmap_vx[0]*f[1]*vmap_vy1_sq+0.7071067811865475*vmap_vx[2]*f[10]*vmap_vy0_sq+0.7071067811865475*vmap_vx[1]*f[4]*vmap_vy0_sq+0.7071067811865475*vmap_vx[0]*f[1]*vmap_vy0_sq+1.2421180068162374*vmap_vy[1]*vmap_vx[2]*vmap_vy[3]*f[17]+1.4142135623730951*vmap_vy[0]*vmap_vx[2]*vmap_vy[2]*f[17]+1.2421180068162374*vmap_vx[1]*vmap_vy[1]*vmap_vy[3]*f[15]+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*vmap_vy[2]*f[15]+1.2421180068162374*vmap_vx[2]*vmap_vy[2]*vmap_vy[3]*f[14]+1.264911064067352*vmap_vy[1]*vmap_vx[2]*vmap_vy[2]*f[14]+1.4142135623730951*vmap_vy[0]*vmap_vy[1]*vmap_vx[2]*f[14]+1.2421180068162376*vmap_vx[0]*vmap_vy[1]*vmap_vy[3]*f[12]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vy[2]*f[12]+1.2421180068162374*vmap_vx[1]*vmap_vy[2]*vmap_vy[3]*f[9]+1.264911064067352*vmap_vx[1]*vmap_vy[1]*vmap_vy[2]*f[9]+1.4142135623730951*vmap_vy[0]*vmap_vx[1]*vmap_vy[1]*f[9]+1.2421180068162374*vmap_vx[0]*vmap_vy[2]*vmap_vy[3]*f[5]+1.264911064067352*vmap_vx[0]*vmap_vy[1]*vmap_vy[2]*f[5]+1.4142135623730951*vmap_vx[0]*vmap_vy[0]*vmap_vy[1]*f[5])*volFact; 
  out[6] += (3.9212285138526726*vmap_vy[2]*f[8]*vmap_vy3_sq+1.264911064067352*vmap_vy[0]*f[8]*vmap_vy3_sq+3.252691193458118*vmap_vy[1]*f[3]*vmap_vy3_sq+1.264911064067352*f[0]*vmap_vy[2]*vmap_vy3_sq+2.1213203435596424*f[0]*vmap_vy[0]*vmap_vy3_sq+0.7069898943717666*f[3]*vmap_vy3_cu+1.355261854357877*vmap_vy[0]*f[8]*vmap_vy2_sq+2.777460299317654*f[3]*vmap_vy[3]*vmap_vy2_sq+3.3335033970222954*vmap_vy[1]*f[3]*vmap_vy2_sq+2.1213203435596424*f[0]*vmap_vy[0]*vmap_vy2_sq+1.515228816828316*f[8]*vmap_vy2_cu+0.45175395145262565*f[0]*vmap_vy2_cu+3.3335033970222954*vmap_vy[2]*f[8]*vmap_vy1_sq+1.8973665961010278*vmap_vy[0]*f[8]*vmap_vy1_sq+1.6664761795905922*f[3]*vmap_vy[3]*vmap_vy1_sq+1.8973665961010278*f[0]*vmap_vy[2]*vmap_vy1_sq+2.1213203435596424*f[0]*vmap_vy[0]*vmap_vy1_sq+1.2727922061357855*f[3]*vmap_vy1_cu+2.1213203435596424*vmap_vy[2]*f[8]*vmap_vy0_sq+2.1213203435596424*vmap_vy[1]*f[3]*vmap_vy0_sq+0.7071067811865475*f[0]*vmap_vy0_cu+5.5549205986353085*vmap_vy[1]*vmap_vy[2]*vmap_vy[3]*f[8]+3.7263540204487113*vmap_vy[0]*vmap_vy[1]*vmap_vy[3]*f[8]+3.7263540204487113*vmap_vy[0]*vmap_vy[2]*f[3]*vmap_vy[3]+3.7263540204487113*f[0]*vmap_vy[1]*vmap_vy[2]*vmap_vy[3]+3.7947331922020555*vmap_vy[0]*vmap_vy[1]*vmap_vy[2]*f[3])*volFact; 
  out[7] += (3.921228513852674*vmap_vy[2]*f[12]*vmap_vy3_sq+1.264911064067352*vmap_vy[0]*f[12]*vmap_vy3_sq+3.252691193458118*vmap_vy[1]*f[5]*vmap_vy3_sq+1.264911064067352*f[1]*vmap_vy[2]*vmap_vy3_sq+2.1213203435596424*vmap_vy[0]*f[1]*vmap_vy3_sq+0.7069898943717666*f[5]*vmap_vy3_cu+1.355261854357877*vmap_vy[0]*f[12]*vmap_vy2_sq+2.777460299317654*vmap_vy[3]*f[5]*vmap_vy2_sq+3.3335033970222954*vmap_vy[1]*f[5]*vmap_vy2_sq+2.1213203435596424*vmap_vy[0]*f[1]*vmap_vy2_sq+1.515228816828316*f[12]*vmap_vy2_cu+0.45175395145262565*f[1]*vmap_vy2_cu+3.3335033970222954*vmap_vy[2]*f[12]*vmap_vy1_sq+1.897366596101028*vmap_vy[0]*f[12]*vmap_vy1_sq+1.6664761795905922*vmap_vy[3]*f[5]*vmap_vy1_sq+1.8973665961010278*f[1]*vmap_vy[2]*vmap_vy1_sq+2.1213203435596424*vmap_vy[0]*f[1]*vmap_vy1_sq+1.2727922061357855*f[5]*vmap_vy1_cu+2.1213203435596424*vmap_vy[2]*f[12]*vmap_vy0_sq+2.1213203435596424*vmap_vy[1]*f[5]*vmap_vy0_sq+0.7071067811865475*f[1]*vmap_vy0_cu+5.554920598635309*vmap_vy[1]*vmap_vy[2]*vmap_vy[3]*f[12]+3.726354020448713*vmap_vy[0]*vmap_vy[1]*vmap_vy[3]*f[12]+3.7263540204487113*vmap_vy[0]*vmap_vy[2]*vmap_vy[3]*f[5]+3.7947331922020555*vmap_vy[0]*vmap_vy[1]*vmap_vy[2]*f[5]+3.7263540204487113*f[1]*vmap_vy[1]*vmap_vy[2]*vmap_vy[3])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M1i_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  const double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  double bdiv_vy[9]; 
  bdiv_vy[0] = 0.2777777777777778/jacob_vy[2]+0.4444444444444444/jacob_vy[1]+0.2777777777777778/jacob_vy[0]; 
  bdiv_vy[1] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[2] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[3] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[4] = 0.5/jacob_vy[2]+0.5/jacob_vy[0]; 
  bdiv_vy[5] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[6] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[7] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[8] = 0.2222222222222222/jacob_vy[2]+0.5555555555555556/jacob_vy[1]+0.2222222222222222/jacob_vy[0]; 
  out[0] += 3.872983346207417*bdiv_vx[5]*hamil[8]*f[16]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[7]*f[16]*dv10*volFact+3.872983346207417*bdiv_vx[4]*hamil[8]*f[13]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[7]*f[13]*dv10*volFact+3.872983346207417*bdiv_vx[5]*hamil[6]*f[11]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[3]*f[11]*dv10*volFact+3.872983346207417*bdiv_vx[3]*f[8]*hamil[8]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*hamil[7]*f[8]*dv10*volFact+3.872983346207417*hamil[4]*bdiv_vx[5]*f[7]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[7]*dv10*volFact+3.872983346207417*bdiv_vx[4]*f[6]*hamil[6]*dv10*volFact+3.872983346207417*bdiv_vx[3]*f[3]*hamil[6]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[3]*f[6]*dv10*volFact+3.872983346207417*f[2]*bdiv_vx[4]*hamil[4]*dv10*volFact+3.872983346207417*f[0]*bdiv_vx[3]*hamil[4]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[3]*hamil[3]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[2]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[0]*hamil[1]*dv10*volFact; 
  out[2] += 3.872983346207417*bdiv_vy[5]*hamil[8]*f[16]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[6]*f[16]*dv11*volFact+3.872983346207417*bdiv_vy[5]*hamil[7]*f[13]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[3]*f[13]*dv11*volFact+3.872983346207417*bdiv_vy[4]*hamil[8]*f[11]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[6]*f[11]*dv11*volFact+3.872983346207417*bdiv_vy[3]*f[7]*hamil[8]*dv11*volFact+3.872983346207417*bdiv_vy[5]*hamil[5]*f[8]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[2]*f[8]*dv11*volFact+3.872983346207417*bdiv_vy[4]*f[6]*hamil[7]*dv11*volFact+3.872983346207417*f[2]*bdiv_vy[3]*hamil[7]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*hamil[6]*f[7]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[3]*f[6]*dv11*volFact+3.872983346207417*f[3]*bdiv_vy[4]*hamil[5]*dv11*volFact+3.872983346207417*f[0]*bdiv_vy[3]*hamil[5]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[2]*hamil[3]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[2]*f[3]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[0]*hamil[2]*dv11*volFact; 
  out[1] += 3.872983346207417*bdiv_vx[5]*hamil[8]*f[17]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[7]*f[17]*dv10*volFact+3.872983346207417*bdiv_vx[4]*hamil[8]*f[15]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[7]*f[15]*dv10*volFact+3.872983346207417*bdiv_vx[5]*hamil[6]*f[14]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[3]*f[14]*dv10*volFact+3.872983346207417*bdiv_vx[3]*hamil[8]*f[12]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*hamil[7]*f[12]*dv10*volFact+3.872983346207417*hamil[4]*bdiv_vx[5]*f[10]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[10]*dv10*volFact+3.872983346207417*bdiv_vx[4]*hamil[6]*f[9]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[3]*f[9]*dv10*volFact+3.872983346207417*bdiv_vx[3]*f[5]*hamil[6]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*hamil[3]*f[5]*dv10*volFact+3.872983346207417*bdiv_vx[4]*f[4]*hamil[4]*dv10*volFact+3.872983346207417*f[1]*bdiv_vx[3]*hamil[4]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[4]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[1]*hamil[1]*dv10*volFact; 
  out[3] += 3.872983346207417*bdiv_vy[5]*hamil[8]*f[17]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[6]*f[17]*dv11*volFact+3.872983346207417*bdiv_vy[5]*hamil[7]*f[15]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[3]*f[15]*dv11*volFact+3.872983346207417*bdiv_vy[4]*hamil[8]*f[14]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[6]*f[14]*dv11*volFact+3.872983346207417*bdiv_vy[5]*hamil[5]*f[12]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[2]*f[12]*dv11*volFact+3.872983346207417*bdiv_vy[3]*hamil[8]*f[10]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*hamil[6]*f[10]*dv11*volFact+3.872983346207417*bdiv_vy[4]*hamil[7]*f[9]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[3]*f[9]*dv11*volFact+3.872983346207417*bdiv_vy[3]*f[4]*hamil[7]*dv11*volFact+3.872983346207417*bdiv_vy[4]*f[5]*hamil[5]*dv11*volFact+3.872983346207417*f[1]*bdiv_vy[3]*hamil[5]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[2]*f[5]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*hamil[3]*f[4]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[1]*hamil[2]*dv11*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M2_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += (hamil[8]*f[16]+hamil[7]*f[13]+hamil[6]*f[11]+hamil[5]*f[8]+hamil[4]*f[7]+hamil[3]*f[6]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[1] += (hamil[8]*f[17]+1.0000000000000002*hamil[7]*f[15]+1.0000000000000002*hamil[6]*f[14]+1.0000000000000002*hamil[5]*f[12]+1.0000000000000002*hamil[4]*f[10]+hamil[3]*f[9]+hamil[2]*f[5]+hamil[1]*f[4]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  const double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  double bdiv_vy[9]; 
  bdiv_vy[0] = 0.2777777777777778/jacob_vy[2]+0.4444444444444444/jacob_vy[1]+0.2777777777777778/jacob_vy[0]; 
  bdiv_vy[1] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[2] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[3] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[4] = 0.5/jacob_vy[2]+0.5/jacob_vy[0]; 
  bdiv_vy[5] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[6] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[7] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[8] = 0.2222222222222222/jacob_vy[2]+0.5555555555555556/jacob_vy[1]+0.2222222222222222/jacob_vy[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[2] += 3.872983346207417*bdiv_vx[5]*hamil[8]*f[16]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[7]*f[16]*dv10*volFact+3.872983346207417*bdiv_vx[4]*hamil[8]*f[13]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[7]*f[13]*dv10*volFact+3.872983346207417*bdiv_vx[5]*hamil[6]*f[11]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[3]*f[11]*dv10*volFact+3.872983346207417*bdiv_vx[3]*f[8]*hamil[8]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*hamil[7]*f[8]*dv10*volFact+3.872983346207417*hamil[4]*bdiv_vx[5]*f[7]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[7]*dv10*volFact+3.872983346207417*bdiv_vx[4]*f[6]*hamil[6]*dv10*volFact+3.872983346207417*bdiv_vx[3]*f[3]*hamil[6]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[3]*f[6]*dv10*volFact+3.872983346207417*f[2]*bdiv_vx[4]*hamil[4]*dv10*volFact+3.872983346207417*f[0]*bdiv_vx[3]*hamil[4]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[3]*hamil[3]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[2]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[0]*hamil[1]*dv10*volFact; 
  out[4] += 3.872983346207417*bdiv_vy[5]*hamil[8]*f[16]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[6]*f[16]*dv11*volFact+3.872983346207417*bdiv_vy[5]*hamil[7]*f[13]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[3]*f[13]*dv11*volFact+3.872983346207417*bdiv_vy[4]*hamil[8]*f[11]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[6]*f[11]*dv11*volFact+3.872983346207417*bdiv_vy[3]*f[7]*hamil[8]*dv11*volFact+3.872983346207417*bdiv_vy[5]*hamil[5]*f[8]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[2]*f[8]*dv11*volFact+3.872983346207417*bdiv_vy[4]*f[6]*hamil[7]*dv11*volFact+3.872983346207417*f[2]*bdiv_vy[3]*hamil[7]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*hamil[6]*f[7]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[3]*f[6]*dv11*volFact+3.872983346207417*f[3]*bdiv_vy[4]*hamil[5]*dv11*volFact+3.872983346207417*f[0]*bdiv_vy[3]*hamil[5]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[2]*hamil[3]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[2]*f[3]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[0]*hamil[2]*dv11*volFact; 
  out[6] += hamil[8]*f[16]*volFact+hamil[7]*f[13]*volFact+hamil[6]*f[11]*volFact+hamil[5]*f[8]*volFact+hamil[4]*f[7]*volFact+hamil[3]*f[6]*volFact+hamil[2]*f[3]*volFact+hamil[1]*f[2]*volFact+f[0]*hamil[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[3] += 3.872983346207417*bdiv_vx[5]*hamil[8]*f[17]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[7]*f[17]*dv10*volFact+3.872983346207417*bdiv_vx[4]*hamil[8]*f[15]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[7]*f[15]*dv10*volFact+3.872983346207417*bdiv_vx[5]*hamil[6]*f[14]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[3]*f[14]*dv10*volFact+3.872983346207417*bdiv_vx[3]*hamil[8]*f[12]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*hamil[7]*f[12]*dv10*volFact+3.872983346207417*hamil[4]*bdiv_vx[5]*f[10]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[10]*dv10*volFact+3.872983346207417*bdiv_vx[4]*hamil[6]*f[9]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[3]*f[9]*dv10*volFact+3.872983346207417*bdiv_vx[3]*f[5]*hamil[6]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*hamil[3]*f[5]*dv10*volFact+3.872983346207417*bdiv_vx[4]*f[4]*hamil[4]*dv10*volFact+3.872983346207417*f[1]*bdiv_vx[3]*hamil[4]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[4]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[1]*hamil[1]*dv10*volFact; 
  out[5] += 3.872983346207417*bdiv_vy[5]*hamil[8]*f[17]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[6]*f[17]*dv11*volFact+3.872983346207417*bdiv_vy[5]*hamil[7]*f[15]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[3]*f[15]*dv11*volFact+3.872983346207417*bdiv_vy[4]*hamil[8]*f[14]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[6]*f[14]*dv11*volFact+3.872983346207417*bdiv_vy[5]*hamil[5]*f[12]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[2]*f[12]*dv11*volFact+3.872983346207417*bdiv_vy[3]*hamil[8]*f[10]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*hamil[6]*f[10]*dv11*volFact+3.872983346207417*bdiv_vy[4]*hamil[7]*f[9]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[3]*f[9]*dv11*volFact+3.872983346207417*bdiv_vy[3]*f[4]*hamil[7]*dv11*volFact+3.872983346207417*bdiv_vy[4]*f[5]*hamil[5]*dv11*volFact+3.872983346207417*f[1]*bdiv_vy[3]*hamil[5]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[2]*f[5]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*hamil[3]*f[4]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[1]*hamil[2]*dv11*volFact; 
  out[7] += hamil[8]*f[17]*volFact+hamil[7]*f[15]*volFact+hamil[6]*f[14]*volFact+hamil[5]*f[12]*volFact+hamil[4]*f[10]*volFact+hamil[3]*f[9]*volFact+hamil[2]*f[5]*volFact+hamil[1]*f[4]*volFact+hamil[0]*f[1]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_int_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  double bdiv_vy[9]; 
  bdiv_vy[0] = 0.2777777777777778/jacob_vy[2]+0.4444444444444444/jacob_vy[1]+0.2777777777777778/jacob_vy[0]; 
  bdiv_vy[1] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[2] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[3] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[4] = 0.5/jacob_vy[2]+0.5/jacob_vy[0]; 
  bdiv_vy[5] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[6] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[7] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[8] = 0.2222222222222222/jacob_vy[2]+0.5555555555555556/jacob_vy[1]+0.2222222222222222/jacob_vy[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 5.477225575051662*bdiv_vx[5]*hamil[8]*f[16]*dv10*volFact+2.4494897427831783*bdiv_vx[2]*hamil[7]*f[16]*dv10*volFact+5.477225575051662*bdiv_vx[4]*hamil[8]*f[13]*dv10*volFact+2.4494897427831783*bdiv_vx[1]*hamil[7]*f[13]*dv10*volFact+5.477225575051662*bdiv_vx[5]*hamil[6]*f[11]*dv10*volFact+2.4494897427831783*bdiv_vx[2]*hamil[3]*f[11]*dv10*volFact+5.477225575051662*bdiv_vx[3]*f[8]*hamil[8]*dv10*volFact+2.4494897427831783*bdiv_vx[0]*hamil[7]*f[8]*dv10*volFact+5.477225575051662*hamil[4]*bdiv_vx[5]*f[7]*dv10*volFact+2.4494897427831783*hamil[1]*bdiv_vx[2]*f[7]*dv10*volFact+5.477225575051662*bdiv_vx[4]*f[6]*hamil[6]*dv10*volFact+5.477225575051662*bdiv_vx[3]*f[3]*hamil[6]*dv10*volFact+2.4494897427831783*bdiv_vx[1]*hamil[3]*f[6]*dv10*volFact+5.477225575051662*f[2]*bdiv_vx[4]*hamil[4]*dv10*volFact+5.477225575051662*f[0]*bdiv_vx[3]*hamil[4]*dv10*volFact+2.4494897427831783*bdiv_vx[0]*f[3]*hamil[3]*dv10*volFact+2.4494897427831783*bdiv_vx[1]*hamil[1]*f[2]*dv10*volFact+2.4494897427831783*bdiv_vx[0]*f[0]*hamil[1]*dv10*volFact; 
  out[2] += 5.477225575051662*bdiv_vy[5]*hamil[8]*f[16]*dv11*volFact+2.4494897427831783*bdiv_vy[2]*hamil[6]*f[16]*dv11*volFact+5.477225575051662*bdiv_vy[5]*hamil[7]*f[13]*dv11*volFact+2.4494897427831783*bdiv_vy[2]*hamil[3]*f[13]*dv11*volFact+5.477225575051662*bdiv_vy[4]*hamil[8]*f[11]*dv11*volFact+2.4494897427831783*bdiv_vy[1]*hamil[6]*f[11]*dv11*volFact+5.477225575051662*bdiv_vy[3]*f[7]*hamil[8]*dv11*volFact+5.477225575051662*bdiv_vy[5]*hamil[5]*f[8]*dv11*volFact+2.4494897427831783*bdiv_vy[2]*hamil[2]*f[8]*dv11*volFact+5.477225575051662*bdiv_vy[4]*f[6]*hamil[7]*dv11*volFact+5.477225575051662*f[2]*bdiv_vy[3]*hamil[7]*dv11*volFact+2.4494897427831783*bdiv_vy[0]*hamil[6]*f[7]*dv11*volFact+2.4494897427831783*bdiv_vy[1]*hamil[3]*f[6]*dv11*volFact+5.477225575051662*f[3]*bdiv_vy[4]*hamil[5]*dv11*volFact+5.477225575051662*f[0]*bdiv_vy[3]*hamil[5]*dv11*volFact+2.4494897427831783*bdiv_vy[0]*f[2]*hamil[3]*dv11*volFact+2.4494897427831783*bdiv_vy[1]*hamil[2]*f[3]*dv11*volFact+2.4494897427831783*bdiv_vy[0]*f[0]*hamil[2]*dv11*volFact; 
  out[3] += 1.4142135623730951*hamil[8]*f[16]*volFact+1.4142135623730951*hamil[7]*f[13]*volFact+1.4142135623730951*hamil[6]*f[11]*volFact+1.4142135623730951*hamil[5]*f[8]*volFact+1.4142135623730951*hamil[4]*f[7]*volFact+1.4142135623730951*hamil[3]*f[6]*volFact+1.4142135623730951*hamil[2]*f[3]*volFact+1.4142135623730951*hamil[1]*f[2]*volFact+1.4142135623730951*f[0]*hamil[0]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M1i_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  const double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  double bdiv_vy[9]; 
  bdiv_vy[0] = 0.2777777777777778/jacob_vy[2]+0.4444444444444444/jacob_vy[1]+0.2777777777777778/jacob_vy[0]; 
  bdiv_vy[1] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[2] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[3] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[4] = 0.5/jacob_vy[2]+0.5/jacob_vy[0]; 
  bdiv_vy[5] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[6] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[7] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[8] = 0.2222222222222222/jacob_vy[2]+0.5555555555555556/jacob_vy[1]+0.2222222222222222/jacob_vy[0]; 
  out[0] += 3.872983346207417*hamil[4]*bdiv_vx[5]*f[7]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[7]*dv10*volFact+3.872983346207417*f[2]*bdiv_vx[4]*hamil[4]*dv10*volFact+3.872983346207417*f[0]*bdiv_vx[3]*hamil[4]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[2]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[0]*hamil[1]*dv10*volFact; 
  out[2] += 3.872983346207417*bdiv_vy[5]*hamil[5]*f[8]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[2]*f[8]*dv11*volFact+3.872983346207417*f[3]*bdiv_vy[4]*hamil[5]*dv11*volFact+3.872983346207417*f[0]*bdiv_vy[3]*hamil[5]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[2]*f[3]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[0]*hamil[2]*dv11*volFact; 
  out[1] += 3.872983346207417*hamil[4]*bdiv_vx[5]*f[10]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[10]*dv10*volFact+3.872983346207417*bdiv_vx[4]*f[4]*hamil[4]*dv10*volFact+3.872983346207417*f[1]*bdiv_vx[3]*hamil[4]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[4]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[1]*hamil[1]*dv10*volFact; 
  out[3] += 3.872983346207417*bdiv_vy[5]*hamil[5]*f[12]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[2]*f[12]*dv11*volFact+3.872983346207417*bdiv_vy[4]*f[5]*hamil[5]*dv11*volFact+3.872983346207417*f[1]*bdiv_vy[3]*hamil[5]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[2]*f[5]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[1]*hamil[2]*dv11*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_M2_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += (hamil[5]*f[8]+hamil[4]*f[7]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[1] += (1.0000000000000002*hamil[5]*f[12]+1.0000000000000002*hamil[4]*f[10]+hamil[2]*f[5]+hamil[1]*f[4]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  const double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  double bdiv_vy[9]; 
  bdiv_vy[0] = 0.2777777777777778/jacob_vy[2]+0.4444444444444444/jacob_vy[1]+0.2777777777777778/jacob_vy[0]; 
  bdiv_vy[1] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[2] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[3] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[4] = 0.5/jacob_vy[2]+0.5/jacob_vy[0]; 
  bdiv_vy[5] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[6] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[7] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[8] = 0.2222222222222222/jacob_vy[2]+0.5555555555555556/jacob_vy[1]+0.2222222222222222/jacob_vy[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[2] += 3.872983346207417*hamil[4]*bdiv_vx[5]*f[7]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[7]*dv10*volFact+3.872983346207417*f[2]*bdiv_vx[4]*hamil[4]*dv10*volFact+3.872983346207417*f[0]*bdiv_vx[3]*hamil[4]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[2]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[0]*hamil[1]*dv10*volFact; 
  out[4] += 3.872983346207417*bdiv_vy[5]*hamil[5]*f[8]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[2]*f[8]*dv11*volFact+3.872983346207417*f[3]*bdiv_vy[4]*hamil[5]*dv11*volFact+3.872983346207417*f[0]*bdiv_vy[3]*hamil[5]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[2]*f[3]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[0]*hamil[2]*dv11*volFact; 
  out[6] += hamil[5]*f[8]*volFact+hamil[4]*f[7]*volFact+hamil[2]*f[3]*volFact+hamil[1]*f[2]*volFact+f[0]*hamil[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[3] += 3.872983346207417*hamil[4]*bdiv_vx[5]*f[10]*dv10*volFact+1.7320508075688772*hamil[1]*bdiv_vx[2]*f[10]*dv10*volFact+3.872983346207417*bdiv_vx[4]*f[4]*hamil[4]*dv10*volFact+3.872983346207417*f[1]*bdiv_vx[3]*hamil[4]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*hamil[1]*f[4]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[1]*hamil[1]*dv10*volFact; 
  out[5] += 3.872983346207417*bdiv_vy[5]*hamil[5]*f[12]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[2]*f[12]*dv11*volFact+3.872983346207417*bdiv_vy[4]*f[5]*hamil[5]*dv11*volFact+3.872983346207417*f[1]*bdiv_vy[3]*hamil[5]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*hamil[2]*f[5]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[1]*hamil[2]*dv11*volFact; 
  out[7] += hamil[5]*f[12]*volFact+hamil[4]*f[10]*volFact+hamil[2]*f[5]*volFact+hamil[1]*f[4]*volFact+hamil[0]*f[1]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_sparse_int_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  double bdiv_vy[9]; 
  bdiv_vy[0] = 0.2777777777777778/jacob_vy[2]+0.4444444444444444/jacob_vy[1]+0.2777777777777778/jacob_vy[0]; 
  bdiv_vy[1] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[2] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[3] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[4] = 0.5/jacob_vy[2]+0.5/jacob_vy[0]; 
  bdiv_vy[5] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[6] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[7] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[8] = 0.2222222222222222/jacob_vy[2]+0.5555555555555556/jacob_vy[1]+0.2222222222222222/jacob_vy[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 5.477225575051662*hamil[4]*bdiv_vx[5]*f[7]*dv10*volFact+2.4494897427831783*hamil[1]*bdiv_vx[2]*f[7]*dv10*volFact+5.477225575051662*f[2]*bdiv_vx[4]*hamil[4]*dv10*volFact+5.477225575051662*f[0]*bdiv_vx[3]*hamil[4]*dv10*volFact+2.4494897427831783*bdiv_vx[1]*hamil[1]*f[2]*dv10*volFact+2.4494897427831783*bdiv_vx[0]*f[0]*hamil[1]*dv10*volFact; 
  out[2] += 5.477225575051662*bdiv_vy[5]*hamil[5]*f[8]*dv11*volFact+2.4494897427831783*bdiv_vy[2]*hamil[2]*f[8]*dv11*volFact+5.477225575051662*f[3]*bdiv_vy[4]*hamil[5]*dv11*volFact+5.477225575051662*f[0]*bdiv_vy[3]*hamil[5]*dv11*volFact+2.4494897427831783*bdiv_vy[1]*hamil[2]*f[3]*dv11*volFact+2.4494897427831783*bdiv_vy[0]*f[0]*hamil[2]*dv11*volFact; 
  out[3] += 1.4142135623730951*hamil[5]*f[8]*volFact+1.4142135623730951*hamil[4]*f[7]*volFact+1.4142135623730951*hamil[2]*f[3]*volFact+1.4142135623730951*hamil[1]*f[2]*volFact+1.4142135623730951*f[0]*hamil[0]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M1i_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  const double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  double bdiv_vy[9]; 
  bdiv_vy[0] = 0.2777777777777778/jacob_vy[2]+0.4444444444444444/jacob_vy[1]+0.2777777777777778/jacob_vy[0]; 
  bdiv_vy[1] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[2] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[3] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[4] = 0.5/jacob_vy[2]+0.5/jacob_vy[0]; 
  bdiv_vy[5] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[6] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[7] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[8] = 0.2222222222222222/jacob_vy[2]+0.5555555555555556/jacob_vy[1]+0.2222222222222222/jacob_vy[0]; 
  for (int q = 0; q < 2; ++q) { 
  const double s0 = ((q >> 0) & 1) ? 1.0 : -1.0; 
  double f_node[9], h_node[9]; 
  f_node[0] = 0.7071067811865475*f[1]*s0+0.7071067811865475*f[0]; 
  f_node[1] = 0.7071067811865475*f[5]*s0+0.7071067811865475*f[3]; 
  f_node[2] = 0.7071067811865475*f[12]*s0+0.7071067811865475*f[8]; 
  f_node[3] = 0.7071067811865475*f[4]*s0+0.7071067811865475*f[2]; 
  f_node[4] = 0.7071067811865475*f[9]*s0+0.7071067811865475*f[6]; 
  f_node[5] = 0.7071067811865475*f[15]*s0+0.7071067811865475*f[13]; 
  f_node[6] = 0.7071067811865475*f[10]*s0+0.7071067811865475*f[7]; 
  f_node[7] = 0.7071067811865475*f[14]*s0+0.7071067811865475*f[11]; 
  f_node[8] = 0.7071067811865475*f[17]*s0+0.7071067811865475*f[16]; 
  h_node[0] = 0.7071067811865475*hamil[1]*s0+0.7071067811865475*hamil[0]; 
  h_node[1] = 0.7071067811865475*hamil[5]*s0+0.7071067811865475*hamil[3]; 
  h_node[2] = 0.7071067811865475*hamil[12]*s0+0.7071067811865475*hamil[8]; 
  h_node[3] = 0.7071067811865475*hamil[4]*s0+0.7071067811865475*hamil[2]; 
  h_node[4] = 0.7071067811865475*hamil[9]*s0+0.7071067811865475*hamil[6]; 
  h_node[5] = 0.7071067811865475*hamil[15]*s0+0.7071067811865475*hamil[13]; 
  h_node[6] = 0.7071067811865475*hamil[10]*s0+0.7071067811865475*hamil[7]; 
  h_node[7] = 0.7071067811865475*hamil[14]*s0+0.7071067811865475*hamil[11]; 
  h_node[8] = 0.7071067811865475*hamil[17]*s0+0.7071067811865475*hamil[16]; 
  const double m1_x = 3.872983346207417*bdiv_vx[5]*f_node[8]*h_node[8]*dv10+3.872983346207417*bdiv_vx[4]*f_node[5]*h_node[8]*dv10+3.872983346207417*f_node[2]*bdiv_vx[3]*h_node[8]*dv10+1.7320508075688772*bdiv_vx[2]*h_node[5]*f_node[8]*dv10+3.872983346207417*bdiv_vx[5]*f_node[7]*h_node[7]*dv10+3.872983346207417*bdiv_vx[4]*f_node[4]*h_node[7]*dv10+3.872983346207417*f_node[1]*bdiv_vx[3]*h_node[7]*dv10+1.7320508075688772*bdiv_vx[2]*h_node[4]*f_node[7]*dv10+3.872983346207417*bdiv_vx[5]*f_node[6]*h_node[6]*dv10+3.872983346207417*f_node[3]*bdiv_vx[4]*h_node[6]*dv10+3.872983346207417*f_node[0]*bdiv_vx[3]*h_node[6]*dv10+1.7320508075688772*bdiv_vx[2]*h_node[3]*f_node[6]*dv10+1.7320508075688772*bdiv_vx[1]*f_node[5]*h_node[5]*dv10+1.7320508075688772*bdiv_vx[0]*f_node[2]*h_node[5]*dv10+1.7320508075688772*bdiv_vx[1]*f_node[4]*h_node[4]*dv10+1.7320508075688772*bdiv_vx[0]*f_node[1]*h_node[4]*dv10+1.7320508075688772*bdiv_vx[1]*f_node[3]*h_node[3]*dv10+1.7320508075688772*bdiv_vx[0]*f_node[0]*h_node[3]*dv10; 
  const double m1_y = 3.872983346207417*bdiv_vy[5]*f_node[8]*h_node[8]*dv11+3.872983346207417*bdiv_vy[4]*f_node[7]*h_node[8]*dv11+3.872983346207417*bdiv_vy[3]*f_node[6]*h_node[8]*dv11+1.7320508075688772*bdiv_vy[2]*h_node[7]*f_node[8]*dv11+1.7320508075688772*bdiv_vy[1]*f_node[7]*h_node[7]*dv11+1.7320508075688772*bdiv_vy[0]*f_node[6]*h_node[7]*dv11+3.872983346207417*bdiv_vy[5]*f_node[5]*h_node[5]*dv11+3.872983346207417*bdiv_vy[4]*f_node[4]*h_node[5]*dv11+3.872983346207417*bdiv_vy[3]*f_node[3]*h_node[5]*dv11+1.7320508075688772*bdiv_vy[2]*h_node[4]*f_node[5]*dv11+3.872983346207417*f_node[2]*h_node[2]*bdiv_vy[5]*dv11+1.7320508075688772*bdiv_vy[1]*f_node[4]*h_node[4]*dv11+1.7320508075688772*bdiv_vy[0]*f_node[3]*h_node[4]*dv11+3.872983346207417*f_node[1]*h_node[2]*bdiv_vy[4]*dv11+3.872983346207417*f_node[0]*h_node[2]*bdiv_vy[3]*dv11+1.7320508075688772*h_node[1]*bdiv_vy[2]*f_node[2]*dv11+1.7320508075688772*bdiv_vy[1]*f_node[1]*h_node[1]*dv11+1.7320508075688772*bdiv_vy[0]*f_node[0]*h_node[1]*dv11; 
  out[0] += 0.7071067811865475*m1_x*volFact; 
  out[2] += 0.7071067811865475*m1_y*volFact; 
  out[1] += 0.7071067811865475*m1_x*s0*volFact; 
  out[3] += 0.7071067811865475*m1_y*s0*volFact; 
  } 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_M2_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += (0.7071067811865475*f[17]*hamil[17]+0.7071067811865475*f[16]*hamil[16]+0.7071067811865475*f[15]*hamil[15]+0.7071067811865475*f[14]*hamil[14]+0.7071067811865475*f[13]*hamil[13]+0.7071067811865475*f[12]*hamil[12]+0.7071067811865475*f[11]*hamil[11]+0.7071067811865475*f[10]*hamil[10]+0.7071067811865475*f[9]*hamil[9]+0.7071067811865475*f[8]*hamil[8]+0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865475*f[16]*hamil[17]+0.7071067811865475*hamil[16]*f[17]+0.7071067811865475*f[13]*hamil[15]+0.7071067811865475*hamil[13]*f[15]+0.7071067811865475*f[11]*hamil[14]+0.7071067811865475*hamil[11]*f[14]+0.7071067811865475*f[8]*hamil[12]+0.7071067811865475*hamil[8]*f[12]+0.7071067811865475*f[7]*hamil[10]+0.7071067811865475*hamil[7]*f[10]+0.7071067811865475*f[6]*hamil[9]+0.7071067811865475*hamil[6]*f[9]+0.7071067811865475*f[3]*hamil[5]+0.7071067811865475*hamil[3]*f[5]+0.7071067811865475*f[2]*hamil[4]+0.7071067811865475*hamil[2]*f[4]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  const double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  double bdiv_vy[9]; 
  bdiv_vy[0] = 0.2777777777777778/jacob_vy[2]+0.4444444444444444/jacob_vy[1]+0.2777777777777778/jacob_vy[0]; 
  bdiv_vy[1] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[2] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[3] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[4] = 0.5/jacob_vy[2]+0.5/jacob_vy[0]; 
  bdiv_vy[5] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[6] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[7] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[8] = 0.2222222222222222/jacob_vy[2]+0.5555555555555556/jacob_vy[1]+0.2222222222222222/jacob_vy[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  for (int q = 0; q < 2; ++q) { 
  const double s0 = ((q >> 0) & 1) ? 1.0 : -1.0; 
  double f_node[9], h_node[9]; 
  f_node[0] = 0.7071067811865475*f[1]*s0+0.7071067811865475*f[0]; 
  f_node[1] = 0.7071067811865475*f[5]*s0+0.7071067811865475*f[3]; 
  f_node[2] = 0.7071067811865475*f[12]*s0+0.7071067811865475*f[8]; 
  f_node[3] = 0.7071067811865475*f[4]*s0+0.7071067811865475*f[2]; 
  f_node[4] = 0.7071067811865475*f[9]*s0+0.7071067811865475*f[6]; 
  f_node[5] = 0.7071067811865475*f[15]*s0+0.7071067811865475*f[13]; 
  f_node[6] = 0.7071067811865475*f[10]*s0+0.7071067811865475*f[7]; 
  f_node[7] = 0.7071067811865475*f[14]*s0+0.7071067811865475*f[11]; 
  f_node[8] = 0.7071067811865475*f[17]*s0+0.7071067811865475*f[16]; 
  h_node[0] = 0.7071067811865475*hamil[1]*s0+0.7071067811865475*hamil[0]; 
  h_node[1] = 0.7071067811865475*hamil[5]*s0+0.7071067811865475*hamil[3]; 
  h_node[2] = 0.7071067811865475*hamil[12]*s0+0.7071067811865475*hamil[8]; 
  h_node[3] = 0.7071067811865475*hamil[4]*s0+0.7071067811865475*hamil[2]; 
  h_node[4] = 0.7071067811865475*hamil[9]*s0+0.7071067811865475*hamil[6]; 
  h_node[5] = 0.7071067811865475*hamil[15]*s0+0.7071067811865475*hamil[13]; 
  h_node[6] = 0.7071067811865475*hamil[10]*s0+0.7071067811865475*hamil[7]; 
  h_node[7] = 0.7071067811865475*hamil[14]*s0+0.7071067811865475*hamil[11]; 
  h_node[8] = 0.7071067811865475*hamil[17]*s0+0.7071067811865475*hamil[16]; 
  const double m1_x = 3.872983346207417*bdiv_vx[5]*f_node[8]*h_node[8]*dv10+3.872983346207417*bdiv_vx[4]*f_node[5]*h_node[8]*dv10+3.872983346207417*f_node[2]*bdiv_vx[3]*h_node[8]*dv10+1.7320508075688772*bdiv_vx[2]*h_node[5]*f_node[8]*dv10+3.872983346207417*bdiv_vx[5]*f_node[7]*h_node[7]*dv10+3.872983346207417*bdiv_vx[4]*f_node[4]*h_node[7]*dv10+3.872983346207417*f_node[1]*bdiv_vx[3]*h_node[7]*dv10+1.7320508075688772*bdiv_vx[2]*h_node[4]*f_node[7]*dv10+3.872983346207417*bdiv_vx[5]*f_node[6]*h_node[6]*dv10+3.872983346207417*f_node[3]*bdiv_vx[4]*h_node[6]*dv10+3.872983346207417*f_node[0]*bdiv_vx[3]*h_node[6]*dv10+1.7320508075688772*bdiv_vx[2]*h_node[3]*f_node[6]*dv10+1.7320508075688772*bdiv_vx[1]*f_node[5]*h_node[5]*dv10+1.7320508075688772*bdiv_vx[0]*f_node[2]*h_node[5]*dv10+1.7320508075688772*bdiv_vx[1]*f_node[4]*h_node[4]*dv10+1.7320508075688772*bdiv_vx[0]*f_node[1]*h_node[4]*dv10+1.7320508075688772*bdiv_vx[1]*f_node[3]*h_node[3]*dv10+1.7320508075688772*bdiv_vx[0]*f_node[0]*h_node[3]*dv10; 
  const double m1_y = 3.872983346207417*bdiv_vy[5]*f_node[8]*h_node[8]*dv11+3.872983346207417*bdiv_vy[4]*f_node[7]*h_node[8]*dv11+3.872983346207417*bdiv_vy[3]*f_node[6]*h_node[8]*dv11+1.7320508075688772*bdiv_vy[2]*h_node[7]*f_node[8]*dv11+1.7320508075688772*bdiv_vy[1]*f_node[7]*h_node[7]*dv11+1.7320508075688772*bdiv_vy[0]*f_node[6]*h_node[7]*dv11+3.872983346207417*bdiv_vy[5]*f_node[5]*h_node[5]*dv11+3.872983346207417*bdiv_vy[4]*f_node[4]*h_node[5]*dv11+3.872983346207417*bdiv_vy[3]*f_node[3]*h_node[5]*dv11+1.7320508075688772*bdiv_vy[2]*h_node[4]*f_node[5]*dv11+3.872983346207417*f_node[2]*h_node[2]*bdiv_vy[5]*dv11+1.7320508075688772*bdiv_vy[1]*f_node[4]*h_node[4]*dv11+1.7320508075688772*bdiv_vy[0]*f_node[3]*h_node[4]*dv11+3.872983346207417*f_node[1]*h_node[2]*bdiv_vy[4]*dv11+3.872983346207417*f_node[0]*h_node[2]*bdiv_vy[3]*dv11+1.7320508075688772*h_node[1]*bdiv_vy[2]*f_node[2]*dv11+1.7320508075688772*bdiv_vy[1]*f_node[1]*h_node[1]*dv11+1.7320508075688772*bdiv_vy[0]*f_node[0]*h_node[1]*dv11; 
  const double m2q = f_node[8]*h_node[8]+f_node[7]*h_node[7]+f_node[6]*h_node[6]+f_node[5]*h_node[5]+f_node[4]*h_node[4]+f_node[3]*h_node[3]+f_node[2]*h_node[2]+f_node[1]*h_node[1]+f_node[0]*h_node[0]; 
  out[2] += 0.7071067811865475*m1_x*volFact; 
  out[4] += 0.7071067811865475*m1_y*volFact; 
  out[6] += 0.7071067811865475*m2q*volFact; 
  out[3] += 0.7071067811865475*m1_x*s0*volFact; 
  out[5] += 0.7071067811865475*m1_y*s0*volFact; 
  out[7] += 0.7071067811865475*m2q*s0*volFact; 
  } 
} 
GKYL_CU_DH void mom_vlasov_hamil_phase_int_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double bdiv_vx[9]; 
  bdiv_vx[0] = 0.2777777777777778/jacob_vx[2]+0.4444444444444444/jacob_vx[1]+0.2777777777777778/jacob_vx[0]; 
  bdiv_vx[1] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[2] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[3] = 0.37267799624996495/jacob_vx[2]-0.37267799624996495/jacob_vx[0]; 
  bdiv_vx[4] = 0.5/jacob_vx[2]+0.5/jacob_vx[0]; 
  bdiv_vx[5] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[6] = 0.24845199749997662/jacob_vx[2]-0.49690399499995325/jacob_vx[1]+0.24845199749997662/jacob_vx[0]; 
  bdiv_vx[7] = 0.3333333333333333/jacob_vx[2]-0.3333333333333333/jacob_vx[0]; 
  bdiv_vx[8] = 0.2222222222222222/jacob_vx[2]+0.5555555555555556/jacob_vx[1]+0.2222222222222222/jacob_vx[0]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *jacob_vy = &jacob_vel[3]; 
  double bdiv_vy[9]; 
  bdiv_vy[0] = 0.2777777777777778/jacob_vy[2]+0.4444444444444444/jacob_vy[1]+0.2777777777777778/jacob_vy[0]; 
  bdiv_vy[1] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[2] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[3] = 0.37267799624996495/jacob_vy[2]-0.37267799624996495/jacob_vy[0]; 
  bdiv_vy[4] = 0.5/jacob_vy[2]+0.5/jacob_vy[0]; 
  bdiv_vy[5] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[6] = 0.24845199749997662/jacob_vy[2]-0.49690399499995325/jacob_vy[1]+0.24845199749997662/jacob_vy[0]; 
  bdiv_vy[7] = 0.3333333333333333/jacob_vy[2]-0.3333333333333333/jacob_vy[0]; 
  bdiv_vy[8] = 0.2222222222222222/jacob_vy[2]+0.5555555555555556/jacob_vy[1]+0.2222222222222222/jacob_vy[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 3.872983346207417*bdiv_vx[5]*f[17]*hamil[17]*dv10*volFact+3.872983346207417*bdiv_vx[4]*f[15]*hamil[17]*dv10*volFact+3.872983346207417*bdiv_vx[3]*f[12]*hamil[17]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[15]*f[17]*dv10*volFact+3.872983346207417*bdiv_vx[5]*f[16]*hamil[16]*dv10*volFact+3.872983346207417*bdiv_vx[4]*f[13]*hamil[16]*dv10*volFact+3.872983346207417*bdiv_vx[3]*f[8]*hamil[16]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[13]*f[16]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*f[15]*hamil[15]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[12]*hamil[15]*dv10*volFact+3.872983346207417*bdiv_vx[5]*f[14]*hamil[14]*dv10*volFact+3.872983346207417*bdiv_vx[4]*f[9]*hamil[14]*dv10*volFact+3.872983346207417*bdiv_vx[3]*f[5]*hamil[14]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[9]*f[14]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*f[13]*hamil[13]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[8]*hamil[13]*dv10*volFact+3.872983346207417*bdiv_vx[5]*f[11]*hamil[11]*dv10*volFact+3.872983346207417*bdiv_vx[4]*f[6]*hamil[11]*dv10*volFact+3.872983346207417*bdiv_vx[3]*f[3]*hamil[11]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[6]*f[11]*dv10*volFact+3.872983346207417*bdiv_vx[5]*f[10]*hamil[10]*dv10*volFact+3.872983346207417*bdiv_vx[4]*f[4]*hamil[10]*dv10*volFact+3.872983346207417*f[1]*bdiv_vx[3]*hamil[10]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[4]*f[10]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*f[9]*hamil[9]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[5]*hamil[9]*dv10*volFact+3.872983346207417*bdiv_vx[5]*f[7]*hamil[7]*dv10*volFact+3.872983346207417*f[2]*bdiv_vx[4]*hamil[7]*dv10*volFact+3.872983346207417*f[0]*bdiv_vx[3]*hamil[7]*dv10*volFact+1.7320508075688772*bdiv_vx[2]*hamil[2]*f[7]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*f[6]*hamil[6]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[3]*hamil[6]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*f[4]*hamil[4]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[1]*hamil[4]*dv10*volFact+1.7320508075688772*bdiv_vx[1]*f[2]*hamil[2]*dv10*volFact+1.7320508075688772*bdiv_vx[0]*f[0]*hamil[2]*dv10*volFact; 
  out[2] += 3.872983346207417*bdiv_vy[5]*f[17]*hamil[17]*dv11*volFact+3.872983346207417*bdiv_vy[4]*f[14]*hamil[17]*dv11*volFact+3.872983346207417*bdiv_vy[3]*f[10]*hamil[17]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[14]*f[17]*dv11*volFact+3.872983346207417*bdiv_vy[5]*f[16]*hamil[16]*dv11*volFact+3.872983346207417*bdiv_vy[4]*f[11]*hamil[16]*dv11*volFact+3.872983346207417*bdiv_vy[3]*f[7]*hamil[16]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[11]*f[16]*dv11*volFact+3.872983346207417*bdiv_vy[5]*f[15]*hamil[15]*dv11*volFact+3.872983346207417*bdiv_vy[4]*f[9]*hamil[15]*dv11*volFact+3.872983346207417*bdiv_vy[3]*f[4]*hamil[15]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[9]*f[15]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*f[14]*hamil[14]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[10]*hamil[14]*dv11*volFact+3.872983346207417*bdiv_vy[5]*f[13]*hamil[13]*dv11*volFact+3.872983346207417*bdiv_vy[4]*f[6]*hamil[13]*dv11*volFact+3.872983346207417*f[2]*bdiv_vy[3]*hamil[13]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[6]*f[13]*dv11*volFact+3.872983346207417*bdiv_vy[5]*f[12]*hamil[12]*dv11*volFact+3.872983346207417*bdiv_vy[4]*f[5]*hamil[12]*dv11*volFact+3.872983346207417*f[1]*bdiv_vy[3]*hamil[12]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[5]*f[12]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*f[11]*hamil[11]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[7]*hamil[11]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*f[9]*hamil[9]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[4]*hamil[9]*dv11*volFact+3.872983346207417*bdiv_vy[5]*f[8]*hamil[8]*dv11*volFact+3.872983346207417*f[3]*bdiv_vy[4]*hamil[8]*dv11*volFact+3.872983346207417*f[0]*bdiv_vy[3]*hamil[8]*dv11*volFact+1.7320508075688772*bdiv_vy[2]*hamil[3]*f[8]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*f[6]*hamil[6]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[2]*hamil[6]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*f[5]*hamil[5]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[1]*hamil[5]*dv11*volFact+1.7320508075688772*bdiv_vy[1]*f[3]*hamil[3]*dv11*volFact+1.7320508075688772*bdiv_vy[0]*f[0]*hamil[3]*dv11*volFact; 
  out[3] += f[17]*hamil[17]*volFact+f[16]*hamil[16]*volFact+f[15]*hamil[15]*volFact+f[14]*hamil[14]*volFact+f[13]*hamil[13]*volFact+f[12]*hamil[12]*volFact+f[11]*hamil[11]*volFact+f[10]*hamil[10]*volFact+f[9]*hamil[9]*volFact+f[8]*hamil[8]*volFact+f[7]*hamil[7]*volFact+f[6]*hamil[6]*volFact+f[5]*hamil[5]*volFact+f[4]*hamil[4]*volFact+f[3]*hamil[3]*volFact+f[2]*hamil[2]*volFact+f[1]*hamil[1]*volFact+f[0]*hamil[0]*volFact; 
} 
