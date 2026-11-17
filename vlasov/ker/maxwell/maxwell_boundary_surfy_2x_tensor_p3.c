#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH double maxwell_boundary_surfy_2x_tensor_p3(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_edge, const double *jacob_pos_skin, const int edge, const double *qedge, const double *qskin, double* GKYL_RESTRICT out) 
{ 
  // meq:            Speed of light and the cleaning speeds over it.
  // w[NDIM]:        Cell-center coordinates.
  // dx[NDIM]:       Cell spacing.
  // jacob_pos_edge: Position-map Jacobian of the interior neighbor of the boundary cell.
  // jacob_pos_skin: Position-map Jacobian of the boundary cell.
  // edge:           -1 on the lower edge of the domain, +1 on the upper edge.
  // qedge:          Field of the interior neighbor of the boundary cell.
  // qskin:          Field of the boundary cell.
  // out:            Incremented output of the boundary cell: the central flux on its
  //                 interior face and the outgoing waves only on the boundary face.

  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
  const double cchi = c*chi, cgamma = c*gamma; 
 
  const double dx1 = 2.0/dx[1]; 
  const double jacob_edge_inv = 1.0/jacob_pos_edge[4]; 
  const double jacob_skin_inv = 1.0/jacob_pos_skin[4]; 
 
  const double *ex_edge = &qedge[0]; 
  const double *ey_edge = &qedge[16]; 
  const double *ez_edge = &qedge[32]; 
  const double *bx_edge = &qedge[48]; 
  const double *by_edge = &qedge[64]; 
  const double *bz_edge = &qedge[80]; 
  const double *ph_edge = &qedge[96]; 
  const double *ps_edge = &qedge[112]; 

  const double *ex_skin = &qskin[0]; 
  const double *ey_skin = &qskin[16]; 
  const double *ez_skin = &qskin[32]; 
  const double *bx_skin = &qskin[48]; 
  const double *by_skin = &qskin[64]; 
  const double *bz_skin = &qskin[80]; 
  const double *ph_skin = &qskin[96]; 
  const double *ps_skin = &qskin[112]; 

  double *outEx = &out[0]; 
  double *outEy = &out[16]; 
  double *outEz = &out[32]; 
  double *outBx = &out[48]; 
  double *outBy = &out[64]; 
  double *outBz = &out[80]; 
  double *outPh = &out[96]; 
  double *outPs = &out[112]; 

  double incr[16]; 
 
  if (edge == -1) { 

  incr[0] = (0.6614378277661477*ex_skin[9]-0.5590169943749475*ex_skin[5]+0.4330127018922193*ex_skin[2]-0.25*ex_skin[0])*c*jacob_skin_inv+c2*((1.3228756555322954*bz_skin[9]+0.8660254037844386*bz_skin[2])*jacob_skin_inv+(-(0.6614378277661477*bz_edge[9])+0.5590169943749475*bz_edge[5]-0.4330127018922193*bz_edge[2]+0.25*bz_edge[0])*jacob_edge_inv); 
  incr[1] = (0.6614378277661477*ex_skin[12]-0.5590169943749476*ex_skin[7]+0.4330127018922193*ex_skin[3]-0.25*ex_skin[1])*c*jacob_skin_inv+c2*((1.3228756555322954*bz_skin[12]+0.8660254037844386*bz_skin[3])*jacob_skin_inv+(-(0.6614378277661477*bz_edge[12])+0.5590169943749476*bz_edge[7]-0.4330127018922193*bz_edge[3]+0.25*bz_edge[1])*jacob_edge_inv); 
  incr[2] = (-(1.14564392373896*ex_skin[9])+0.9682458365518543*ex_skin[5]-0.75*ex_skin[2]+0.4330127018922193*ex_skin[0])*c*jacob_skin_inv+c2*((1.9364916731037085*bz_skin[5]+0.8660254037844386*bz_skin[0])*jacob_skin_inv+(-(1.14564392373896*bz_edge[9])+0.9682458365518543*bz_edge[5]-0.75*bz_edge[2]+0.4330127018922193*bz_edge[0])*jacob_edge_inv); 
  incr[3] = (-(1.14564392373896*ex_skin[12])+0.9682458365518543*ex_skin[7]-0.75*ex_skin[3]+0.4330127018922193*ex_skin[1])*c*jacob_skin_inv+c2*((1.9364916731037085*bz_skin[7]+0.8660254037844386*bz_skin[1])*jacob_skin_inv+(-(1.14564392373896*bz_edge[12])+0.9682458365518543*bz_edge[7]-0.75*bz_edge[3]+0.4330127018922193*bz_edge[1])*jacob_edge_inv); 
  incr[4] = (0.6614378277661477*ex_skin[14]-0.5590169943749475*ex_skin[10]+0.43301270189221935*ex_skin[6]-0.25*ex_skin[4])*c*jacob_skin_inv+c2*((1.3228756555322954*bz_skin[14]+0.8660254037844387*bz_skin[6])*jacob_skin_inv+(-(0.6614378277661477*bz_edge[14])+0.5590169943749475*bz_edge[10]-0.43301270189221935*bz_edge[6]+0.25*bz_edge[4])*jacob_edge_inv); 
  incr[5] = (1.4790199457749043*ex_skin[9]-1.25*ex_skin[5]+0.9682458365518543*ex_skin[2]-0.5590169943749475*ex_skin[0])*c*jacob_skin_inv+c2*((2.9580398915498085*bz_skin[9]+1.9364916731037085*bz_skin[2])*jacob_skin_inv+(-(1.4790199457749043*bz_edge[9])+1.25*bz_edge[5]-0.9682458365518543*bz_edge[2]+0.5590169943749475*bz_edge[0])*jacob_edge_inv); 
  incr[6] = (-(1.14564392373896*ex_skin[14])+0.9682458365518543*ex_skin[10]-0.75*ex_skin[6]+0.43301270189221935*ex_skin[4])*c*jacob_skin_inv+c2*((1.9364916731037085*bz_skin[10]+0.8660254037844387*bz_skin[4])*jacob_skin_inv+(-(1.14564392373896*bz_edge[14])+0.9682458365518543*bz_edge[10]-0.75*bz_edge[6]+0.43301270189221935*bz_edge[4])*jacob_edge_inv); 
  incr[7] = (1.4790199457749038*ex_skin[12]-1.25*ex_skin[7]+0.9682458365518543*ex_skin[3]-0.5590169943749476*ex_skin[1])*c*jacob_skin_inv+c2*((2.9580398915498076*bz_skin[12]+1.9364916731037085*bz_skin[3])*jacob_skin_inv+(-(1.4790199457749038*bz_edge[12])+1.25*bz_edge[7]-0.9682458365518543*bz_edge[3]+0.5590169943749476*bz_edge[1])*jacob_edge_inv); 
  incr[8] = (0.6614378277661477*ex_skin[15]-0.5590169943749473*ex_skin[13]+0.4330127018922193*ex_skin[11]-0.25*ex_skin[8])*c*jacob_skin_inv+c2*((1.3228756555322954*bz_skin[15]+0.8660254037844386*bz_skin[11])*jacob_skin_inv+(-(0.6614378277661477*bz_edge[15])+0.5590169943749473*bz_edge[13]-0.4330127018922193*bz_edge[11]+0.25*bz_edge[8])*jacob_edge_inv); 
  incr[9] = (-(1.75*ex_skin[9])+1.4790199457749043*ex_skin[5]-1.14564392373896*ex_skin[2]+0.6614378277661477*ex_skin[0])*c*jacob_skin_inv+c2*((2.9580398915498085*bz_skin[5]+1.3228756555322954*bz_skin[0])*jacob_skin_inv+(-(1.75*bz_edge[9])+1.4790199457749043*bz_edge[5]-1.14564392373896*bz_edge[2]+0.6614378277661477*bz_edge[0])*jacob_edge_inv); 
  incr[10] = (1.479019945774904*ex_skin[14]-1.25*ex_skin[10]+0.9682458365518543*ex_skin[6]-0.5590169943749475*ex_skin[4])*c*jacob_skin_inv+c2*((2.958039891549808*bz_skin[14]+1.9364916731037085*bz_skin[6])*jacob_skin_inv+(-(1.479019945774904*bz_edge[14])+1.25*bz_edge[10]-0.9682458365518543*bz_edge[6]+0.5590169943749475*bz_edge[4])*jacob_edge_inv); 
  incr[11] = (-(1.14564392373896*ex_skin[15])+0.9682458365518541*ex_skin[13]-0.75*ex_skin[11]+0.4330127018922193*ex_skin[8])*c*jacob_skin_inv+c2*((1.9364916731037083*bz_skin[13]+0.8660254037844386*bz_skin[8])*jacob_skin_inv+(-(1.14564392373896*bz_edge[15])+0.9682458365518541*bz_edge[13]-0.75*bz_edge[11]+0.4330127018922193*bz_edge[8])*jacob_edge_inv); 
  incr[12] = (-(1.75*ex_skin[12])+1.4790199457749038*ex_skin[7]-1.14564392373896*ex_skin[3]+0.6614378277661477*ex_skin[1])*c*jacob_skin_inv+c2*((2.9580398915498076*bz_skin[7]+1.3228756555322954*bz_skin[1])*jacob_skin_inv+(-(1.75*bz_edge[12])+1.4790199457749038*bz_edge[7]-1.14564392373896*bz_edge[3]+0.6614378277661477*bz_edge[1])*jacob_edge_inv); 
  incr[13] = (1.479019945774904*ex_skin[15]-1.25*ex_skin[13]+0.9682458365518541*ex_skin[11]-0.5590169943749473*ex_skin[8])*c*jacob_skin_inv+c2*((2.958039891549808*bz_skin[15]+1.9364916731037083*bz_skin[11])*jacob_skin_inv+(-(1.479019945774904*bz_edge[15])+1.25*bz_edge[13]-0.9682458365518541*bz_edge[11]+0.5590169943749473*bz_edge[8])*jacob_edge_inv); 
  incr[14] = (-(1.75*ex_skin[14])+1.479019945774904*ex_skin[10]-1.14564392373896*ex_skin[6]+0.6614378277661477*ex_skin[4])*c*jacob_skin_inv+c2*((2.958039891549808*bz_skin[10]+1.3228756555322954*bz_skin[4])*jacob_skin_inv+(-(1.75*bz_edge[14])+1.479019945774904*bz_edge[10]-1.14564392373896*bz_edge[6]+0.6614378277661477*bz_edge[4])*jacob_edge_inv); 
  incr[15] = (-(1.75*ex_skin[15])+1.479019945774904*ex_skin[13]-1.14564392373896*ex_skin[11]+0.6614378277661477*ex_skin[8])*c*jacob_skin_inv+c2*((2.958039891549808*bz_skin[13]+1.3228756555322954*bz_skin[8])*jacob_skin_inv+(-(1.75*bz_edge[15])+1.479019945774904*bz_edge[13]-1.14564392373896*bz_edge[11]+0.6614378277661477*bz_edge[8])*jacob_edge_inv); 

  outEx[0] += incr[0]*dx1; 
  outEx[1] += incr[1]*dx1; 
  outEx[2] += incr[2]*dx1; 
  outEx[3] += incr[3]*dx1; 
  outEx[4] += incr[4]*dx1; 
  outEx[5] += incr[5]*dx1; 
  outEx[6] += incr[6]*dx1; 
  outEx[7] += incr[7]*dx1; 
  outEx[8] += incr[8]*dx1; 
  outEx[9] += incr[9]*dx1; 
  outEx[10] += incr[10]*dx1; 
  outEx[11] += incr[11]*dx1; 
  outEx[12] += incr[12]*dx1; 
  outEx[13] += incr[13]*dx1; 
  outEx[14] += incr[14]*dx1; 
  outEx[15] += incr[15]*dx1; 

  incr[0] = c2chi*((-(1.3228756555322954*ph_skin[9])-0.8660254037844386*ph_skin[2])*jacob_skin_inv+(0.6614378277661477*ph_edge[9]-0.5590169943749475*ph_edge[5]+0.4330127018922193*ph_edge[2]-0.25*ph_edge[0])*jacob_edge_inv)+(0.6614378277661477*ey_skin[9]-0.5590169943749475*ey_skin[5]+0.4330127018922193*ey_skin[2]-0.25*ey_skin[0])*cchi*jacob_skin_inv; 
  incr[1] = c2chi*((-(1.3228756555322954*ph_skin[12])-0.8660254037844386*ph_skin[3])*jacob_skin_inv+(0.6614378277661477*ph_edge[12]-0.5590169943749476*ph_edge[7]+0.4330127018922193*ph_edge[3]-0.25*ph_edge[1])*jacob_edge_inv)+(0.6614378277661477*ey_skin[12]-0.5590169943749476*ey_skin[7]+0.4330127018922193*ey_skin[3]-0.25*ey_skin[1])*cchi*jacob_skin_inv; 
  incr[2] = (-(1.14564392373896*ey_skin[9])+0.9682458365518543*ey_skin[5]-0.75*ey_skin[2]+0.4330127018922193*ey_skin[0])*cchi*jacob_skin_inv+c2chi*((-(1.9364916731037085*ph_skin[5])-0.8660254037844386*ph_skin[0])*jacob_skin_inv+(1.14564392373896*ph_edge[9]-0.9682458365518543*ph_edge[5]+0.75*ph_edge[2]-0.4330127018922193*ph_edge[0])*jacob_edge_inv); 
  incr[3] = (-(1.14564392373896*ey_skin[12])+0.9682458365518543*ey_skin[7]-0.75*ey_skin[3]+0.4330127018922193*ey_skin[1])*cchi*jacob_skin_inv+c2chi*((-(1.9364916731037085*ph_skin[7])-0.8660254037844386*ph_skin[1])*jacob_skin_inv+(1.14564392373896*ph_edge[12]-0.9682458365518543*ph_edge[7]+0.75*ph_edge[3]-0.4330127018922193*ph_edge[1])*jacob_edge_inv); 
  incr[4] = c2chi*((-(1.3228756555322954*ph_skin[14])-0.8660254037844387*ph_skin[6])*jacob_skin_inv+(0.6614378277661477*ph_edge[14]-0.5590169943749475*ph_edge[10]+0.43301270189221935*ph_edge[6]-0.25*ph_edge[4])*jacob_edge_inv)+(0.6614378277661477*ey_skin[14]-0.5590169943749475*ey_skin[10]+0.43301270189221935*ey_skin[6]-0.25*ey_skin[4])*cchi*jacob_skin_inv; 
  incr[5] = c2chi*((-(2.9580398915498085*ph_skin[9])-1.9364916731037085*ph_skin[2])*jacob_skin_inv+(1.4790199457749043*ph_edge[9]-1.25*ph_edge[5]+0.9682458365518543*ph_edge[2]-0.5590169943749475*ph_edge[0])*jacob_edge_inv)+(1.4790199457749043*ey_skin[9]-1.25*ey_skin[5]+0.9682458365518543*ey_skin[2]-0.5590169943749475*ey_skin[0])*cchi*jacob_skin_inv; 
  incr[6] = (-(1.14564392373896*ey_skin[14])+0.9682458365518543*ey_skin[10]-0.75*ey_skin[6]+0.43301270189221935*ey_skin[4])*cchi*jacob_skin_inv+c2chi*((-(1.9364916731037085*ph_skin[10])-0.8660254037844387*ph_skin[4])*jacob_skin_inv+(1.14564392373896*ph_edge[14]-0.9682458365518543*ph_edge[10]+0.75*ph_edge[6]-0.43301270189221935*ph_edge[4])*jacob_edge_inv); 
  incr[7] = c2chi*((-(2.9580398915498076*ph_skin[12])-1.9364916731037085*ph_skin[3])*jacob_skin_inv+(1.4790199457749038*ph_edge[12]-1.25*ph_edge[7]+0.9682458365518543*ph_edge[3]-0.5590169943749476*ph_edge[1])*jacob_edge_inv)+(1.4790199457749038*ey_skin[12]-1.25*ey_skin[7]+0.9682458365518543*ey_skin[3]-0.5590169943749476*ey_skin[1])*cchi*jacob_skin_inv; 
  incr[8] = c2chi*((-(1.3228756555322954*ph_skin[15])-0.8660254037844386*ph_skin[11])*jacob_skin_inv+(0.6614378277661477*ph_edge[15]-0.5590169943749473*ph_edge[13]+0.4330127018922193*ph_edge[11]-0.25*ph_edge[8])*jacob_edge_inv)+(0.6614378277661477*ey_skin[15]-0.5590169943749473*ey_skin[13]+0.4330127018922193*ey_skin[11]-0.25*ey_skin[8])*cchi*jacob_skin_inv; 
  incr[9] = (-(1.75*ey_skin[9])+1.4790199457749043*ey_skin[5]-1.14564392373896*ey_skin[2]+0.6614378277661477*ey_skin[0])*cchi*jacob_skin_inv+c2chi*((-(2.9580398915498085*ph_skin[5])-1.3228756555322954*ph_skin[0])*jacob_skin_inv+(1.75*ph_edge[9]-1.4790199457749043*ph_edge[5]+1.14564392373896*ph_edge[2]-0.6614378277661477*ph_edge[0])*jacob_edge_inv); 
  incr[10] = c2chi*((-(2.958039891549808*ph_skin[14])-1.9364916731037085*ph_skin[6])*jacob_skin_inv+(1.479019945774904*ph_edge[14]-1.25*ph_edge[10]+0.9682458365518543*ph_edge[6]-0.5590169943749475*ph_edge[4])*jacob_edge_inv)+(1.479019945774904*ey_skin[14]-1.25*ey_skin[10]+0.9682458365518543*ey_skin[6]-0.5590169943749475*ey_skin[4])*cchi*jacob_skin_inv; 
  incr[11] = (-(1.14564392373896*ey_skin[15])+0.9682458365518541*ey_skin[13]-0.75*ey_skin[11]+0.4330127018922193*ey_skin[8])*cchi*jacob_skin_inv+c2chi*((-(1.9364916731037083*ph_skin[13])-0.8660254037844386*ph_skin[8])*jacob_skin_inv+(1.14564392373896*ph_edge[15]-0.9682458365518541*ph_edge[13]+0.75*ph_edge[11]-0.4330127018922193*ph_edge[8])*jacob_edge_inv); 
  incr[12] = (-(1.75*ey_skin[12])+1.4790199457749038*ey_skin[7]-1.14564392373896*ey_skin[3]+0.6614378277661477*ey_skin[1])*cchi*jacob_skin_inv+c2chi*((-(2.9580398915498076*ph_skin[7])-1.3228756555322954*ph_skin[1])*jacob_skin_inv+(1.75*ph_edge[12]-1.4790199457749038*ph_edge[7]+1.14564392373896*ph_edge[3]-0.6614378277661477*ph_edge[1])*jacob_edge_inv); 
  incr[13] = c2chi*((-(2.958039891549808*ph_skin[15])-1.9364916731037083*ph_skin[11])*jacob_skin_inv+(1.479019945774904*ph_edge[15]-1.25*ph_edge[13]+0.9682458365518541*ph_edge[11]-0.5590169943749473*ph_edge[8])*jacob_edge_inv)+(1.479019945774904*ey_skin[15]-1.25*ey_skin[13]+0.9682458365518541*ey_skin[11]-0.5590169943749473*ey_skin[8])*cchi*jacob_skin_inv; 
  incr[14] = (-(1.75*ey_skin[14])+1.479019945774904*ey_skin[10]-1.14564392373896*ey_skin[6]+0.6614378277661477*ey_skin[4])*cchi*jacob_skin_inv+c2chi*((-(2.958039891549808*ph_skin[10])-1.3228756555322954*ph_skin[4])*jacob_skin_inv+(1.75*ph_edge[14]-1.479019945774904*ph_edge[10]+1.14564392373896*ph_edge[6]-0.6614378277661477*ph_edge[4])*jacob_edge_inv); 
  incr[15] = (-(1.75*ey_skin[15])+1.479019945774904*ey_skin[13]-1.14564392373896*ey_skin[11]+0.6614378277661477*ey_skin[8])*cchi*jacob_skin_inv+c2chi*((-(2.958039891549808*ph_skin[13])-1.3228756555322954*ph_skin[8])*jacob_skin_inv+(1.75*ph_edge[15]-1.479019945774904*ph_edge[13]+1.14564392373896*ph_edge[11]-0.6614378277661477*ph_edge[8])*jacob_edge_inv); 

  outEy[0] += incr[0]*dx1; 
  outEy[1] += incr[1]*dx1; 
  outEy[2] += incr[2]*dx1; 
  outEy[3] += incr[3]*dx1; 
  outEy[4] += incr[4]*dx1; 
  outEy[5] += incr[5]*dx1; 
  outEy[6] += incr[6]*dx1; 
  outEy[7] += incr[7]*dx1; 
  outEy[8] += incr[8]*dx1; 
  outEy[9] += incr[9]*dx1; 
  outEy[10] += incr[10]*dx1; 
  outEy[11] += incr[11]*dx1; 
  outEy[12] += incr[12]*dx1; 
  outEy[13] += incr[13]*dx1; 
  outEy[14] += incr[14]*dx1; 
  outEy[15] += incr[15]*dx1; 

  incr[0] = (0.6614378277661477*ez_skin[9]-0.5590169943749475*ez_skin[5]+0.4330127018922193*ez_skin[2]-0.25*ez_skin[0])*c*jacob_skin_inv+c2*((-(1.3228756555322954*bx_skin[9])-0.8660254037844386*bx_skin[2])*jacob_skin_inv+(0.6614378277661477*bx_edge[9]-0.5590169943749475*bx_edge[5]+0.4330127018922193*bx_edge[2]-0.25*bx_edge[0])*jacob_edge_inv); 
  incr[1] = (0.6614378277661477*ez_skin[12]-0.5590169943749476*ez_skin[7]+0.4330127018922193*ez_skin[3]-0.25*ez_skin[1])*c*jacob_skin_inv+c2*((-(1.3228756555322954*bx_skin[12])-0.8660254037844386*bx_skin[3])*jacob_skin_inv+(0.6614378277661477*bx_edge[12]-0.5590169943749476*bx_edge[7]+0.4330127018922193*bx_edge[3]-0.25*bx_edge[1])*jacob_edge_inv); 
  incr[2] = (-(1.14564392373896*ez_skin[9])+0.9682458365518543*ez_skin[5]-0.75*ez_skin[2]+0.4330127018922193*ez_skin[0])*c*jacob_skin_inv+c2*((-(1.9364916731037085*bx_skin[5])-0.8660254037844386*bx_skin[0])*jacob_skin_inv+(1.14564392373896*bx_edge[9]-0.9682458365518543*bx_edge[5]+0.75*bx_edge[2]-0.4330127018922193*bx_edge[0])*jacob_edge_inv); 
  incr[3] = (-(1.14564392373896*ez_skin[12])+0.9682458365518543*ez_skin[7]-0.75*ez_skin[3]+0.4330127018922193*ez_skin[1])*c*jacob_skin_inv+c2*((-(1.9364916731037085*bx_skin[7])-0.8660254037844386*bx_skin[1])*jacob_skin_inv+(1.14564392373896*bx_edge[12]-0.9682458365518543*bx_edge[7]+0.75*bx_edge[3]-0.4330127018922193*bx_edge[1])*jacob_edge_inv); 
  incr[4] = (0.6614378277661477*ez_skin[14]-0.5590169943749475*ez_skin[10]+0.43301270189221935*ez_skin[6]-0.25*ez_skin[4])*c*jacob_skin_inv+c2*((-(1.3228756555322954*bx_skin[14])-0.8660254037844387*bx_skin[6])*jacob_skin_inv+(0.6614378277661477*bx_edge[14]-0.5590169943749475*bx_edge[10]+0.43301270189221935*bx_edge[6]-0.25*bx_edge[4])*jacob_edge_inv); 
  incr[5] = (1.4790199457749043*ez_skin[9]-1.25*ez_skin[5]+0.9682458365518543*ez_skin[2]-0.5590169943749475*ez_skin[0])*c*jacob_skin_inv+c2*((-(2.9580398915498085*bx_skin[9])-1.9364916731037085*bx_skin[2])*jacob_skin_inv+(1.4790199457749043*bx_edge[9]-1.25*bx_edge[5]+0.9682458365518543*bx_edge[2]-0.5590169943749475*bx_edge[0])*jacob_edge_inv); 
  incr[6] = (-(1.14564392373896*ez_skin[14])+0.9682458365518543*ez_skin[10]-0.75*ez_skin[6]+0.43301270189221935*ez_skin[4])*c*jacob_skin_inv+c2*((-(1.9364916731037085*bx_skin[10])-0.8660254037844387*bx_skin[4])*jacob_skin_inv+(1.14564392373896*bx_edge[14]-0.9682458365518543*bx_edge[10]+0.75*bx_edge[6]-0.43301270189221935*bx_edge[4])*jacob_edge_inv); 
  incr[7] = (1.4790199457749038*ez_skin[12]-1.25*ez_skin[7]+0.9682458365518543*ez_skin[3]-0.5590169943749476*ez_skin[1])*c*jacob_skin_inv+c2*((-(2.9580398915498076*bx_skin[12])-1.9364916731037085*bx_skin[3])*jacob_skin_inv+(1.4790199457749038*bx_edge[12]-1.25*bx_edge[7]+0.9682458365518543*bx_edge[3]-0.5590169943749476*bx_edge[1])*jacob_edge_inv); 
  incr[8] = (0.6614378277661477*ez_skin[15]-0.5590169943749473*ez_skin[13]+0.4330127018922193*ez_skin[11]-0.25*ez_skin[8])*c*jacob_skin_inv+c2*((-(1.3228756555322954*bx_skin[15])-0.8660254037844386*bx_skin[11])*jacob_skin_inv+(0.6614378277661477*bx_edge[15]-0.5590169943749473*bx_edge[13]+0.4330127018922193*bx_edge[11]-0.25*bx_edge[8])*jacob_edge_inv); 
  incr[9] = (-(1.75*ez_skin[9])+1.4790199457749043*ez_skin[5]-1.14564392373896*ez_skin[2]+0.6614378277661477*ez_skin[0])*c*jacob_skin_inv+c2*((-(2.9580398915498085*bx_skin[5])-1.3228756555322954*bx_skin[0])*jacob_skin_inv+(1.75*bx_edge[9]-1.4790199457749043*bx_edge[5]+1.14564392373896*bx_edge[2]-0.6614378277661477*bx_edge[0])*jacob_edge_inv); 
  incr[10] = (1.479019945774904*ez_skin[14]-1.25*ez_skin[10]+0.9682458365518543*ez_skin[6]-0.5590169943749475*ez_skin[4])*c*jacob_skin_inv+c2*((-(2.958039891549808*bx_skin[14])-1.9364916731037085*bx_skin[6])*jacob_skin_inv+(1.479019945774904*bx_edge[14]-1.25*bx_edge[10]+0.9682458365518543*bx_edge[6]-0.5590169943749475*bx_edge[4])*jacob_edge_inv); 
  incr[11] = (-(1.14564392373896*ez_skin[15])+0.9682458365518541*ez_skin[13]-0.75*ez_skin[11]+0.4330127018922193*ez_skin[8])*c*jacob_skin_inv+c2*((-(1.9364916731037083*bx_skin[13])-0.8660254037844386*bx_skin[8])*jacob_skin_inv+(1.14564392373896*bx_edge[15]-0.9682458365518541*bx_edge[13]+0.75*bx_edge[11]-0.4330127018922193*bx_edge[8])*jacob_edge_inv); 
  incr[12] = (-(1.75*ez_skin[12])+1.4790199457749038*ez_skin[7]-1.14564392373896*ez_skin[3]+0.6614378277661477*ez_skin[1])*c*jacob_skin_inv+c2*((-(2.9580398915498076*bx_skin[7])-1.3228756555322954*bx_skin[1])*jacob_skin_inv+(1.75*bx_edge[12]-1.4790199457749038*bx_edge[7]+1.14564392373896*bx_edge[3]-0.6614378277661477*bx_edge[1])*jacob_edge_inv); 
  incr[13] = (1.479019945774904*ez_skin[15]-1.25*ez_skin[13]+0.9682458365518541*ez_skin[11]-0.5590169943749473*ez_skin[8])*c*jacob_skin_inv+c2*((-(2.958039891549808*bx_skin[15])-1.9364916731037083*bx_skin[11])*jacob_skin_inv+(1.479019945774904*bx_edge[15]-1.25*bx_edge[13]+0.9682458365518541*bx_edge[11]-0.5590169943749473*bx_edge[8])*jacob_edge_inv); 
  incr[14] = (-(1.75*ez_skin[14])+1.479019945774904*ez_skin[10]-1.14564392373896*ez_skin[6]+0.6614378277661477*ez_skin[4])*c*jacob_skin_inv+c2*((-(2.958039891549808*bx_skin[10])-1.3228756555322954*bx_skin[4])*jacob_skin_inv+(1.75*bx_edge[14]-1.479019945774904*bx_edge[10]+1.14564392373896*bx_edge[6]-0.6614378277661477*bx_edge[4])*jacob_edge_inv); 
  incr[15] = (-(1.75*ez_skin[15])+1.479019945774904*ez_skin[13]-1.14564392373896*ez_skin[11]+0.6614378277661477*ez_skin[8])*c*jacob_skin_inv+c2*((-(2.958039891549808*bx_skin[13])-1.3228756555322954*bx_skin[8])*jacob_skin_inv+(1.75*bx_edge[15]-1.479019945774904*bx_edge[13]+1.14564392373896*bx_edge[11]-0.6614378277661477*bx_edge[8])*jacob_edge_inv); 

  outEz[0] += incr[0]*dx1; 
  outEz[1] += incr[1]*dx1; 
  outEz[2] += incr[2]*dx1; 
  outEz[3] += incr[3]*dx1; 
  outEz[4] += incr[4]*dx1; 
  outEz[5] += incr[5]*dx1; 
  outEz[6] += incr[6]*dx1; 
  outEz[7] += incr[7]*dx1; 
  outEz[8] += incr[8]*dx1; 
  outEz[9] += incr[9]*dx1; 
  outEz[10] += incr[10]*dx1; 
  outEz[11] += incr[11]*dx1; 
  outEz[12] += incr[12]*dx1; 
  outEz[13] += incr[13]*dx1; 
  outEz[14] += incr[14]*dx1; 
  outEz[15] += incr[15]*dx1; 

  incr[0] = ((0.6614378277661477*bx_skin[9]-0.5590169943749475*bx_skin[5]+0.4330127018922193*bx_skin[2]-0.25*bx_skin[0])*c-1.3228756555322954*ez_skin[9]-0.8660254037844386*ez_skin[2])*jacob_skin_inv+(0.6614378277661477*ez_edge[9]-0.5590169943749475*ez_edge[5]+0.4330127018922193*ez_edge[2]-0.25*ez_edge[0])*jacob_edge_inv; 
  incr[1] = ((0.6614378277661477*bx_skin[12]-0.5590169943749476*bx_skin[7]+0.4330127018922193*bx_skin[3]-0.25*bx_skin[1])*c-1.3228756555322954*ez_skin[12]-0.8660254037844386*ez_skin[3])*jacob_skin_inv+(0.6614378277661477*ez_edge[12]-0.5590169943749476*ez_edge[7]+0.4330127018922193*ez_edge[3]-0.25*ez_edge[1])*jacob_edge_inv; 
  incr[2] = ((-(1.14564392373896*bx_skin[9])+0.9682458365518543*bx_skin[5]-0.75*bx_skin[2]+0.4330127018922193*bx_skin[0])*c-1.9364916731037085*ez_skin[5]-0.8660254037844386*ez_skin[0])*jacob_skin_inv+(1.14564392373896*ez_edge[9]-0.9682458365518543*ez_edge[5]+0.75*ez_edge[2]-0.4330127018922193*ez_edge[0])*jacob_edge_inv; 
  incr[3] = ((-(1.14564392373896*bx_skin[12])+0.9682458365518543*bx_skin[7]-0.75*bx_skin[3]+0.4330127018922193*bx_skin[1])*c-1.9364916731037085*ez_skin[7]-0.8660254037844386*ez_skin[1])*jacob_skin_inv+(1.14564392373896*ez_edge[12]-0.9682458365518543*ez_edge[7]+0.75*ez_edge[3]-0.4330127018922193*ez_edge[1])*jacob_edge_inv; 
  incr[4] = ((0.6614378277661477*bx_skin[14]-0.5590169943749475*bx_skin[10]+0.43301270189221935*bx_skin[6]-0.25*bx_skin[4])*c-1.3228756555322954*ez_skin[14]-0.8660254037844387*ez_skin[6])*jacob_skin_inv+(0.6614378277661477*ez_edge[14]-0.5590169943749475*ez_edge[10]+0.43301270189221935*ez_edge[6]-0.25*ez_edge[4])*jacob_edge_inv; 
  incr[5] = ((1.4790199457749043*bx_skin[9]-1.25*bx_skin[5]+0.9682458365518543*bx_skin[2]-0.5590169943749475*bx_skin[0])*c-2.9580398915498085*ez_skin[9]-1.9364916731037085*ez_skin[2])*jacob_skin_inv+(1.4790199457749043*ez_edge[9]-1.25*ez_edge[5]+0.9682458365518543*ez_edge[2]-0.5590169943749475*ez_edge[0])*jacob_edge_inv; 
  incr[6] = ((-(1.14564392373896*bx_skin[14])+0.9682458365518543*bx_skin[10]-0.75*bx_skin[6]+0.43301270189221935*bx_skin[4])*c-1.9364916731037085*ez_skin[10]-0.8660254037844387*ez_skin[4])*jacob_skin_inv+(1.14564392373896*ez_edge[14]-0.9682458365518543*ez_edge[10]+0.75*ez_edge[6]-0.43301270189221935*ez_edge[4])*jacob_edge_inv; 
  incr[7] = ((1.4790199457749038*bx_skin[12]-1.25*bx_skin[7]+0.9682458365518543*bx_skin[3]-0.5590169943749476*bx_skin[1])*c-2.9580398915498076*ez_skin[12]-1.9364916731037085*ez_skin[3])*jacob_skin_inv+(1.4790199457749038*ez_edge[12]-1.25*ez_edge[7]+0.9682458365518543*ez_edge[3]-0.5590169943749476*ez_edge[1])*jacob_edge_inv; 
  incr[8] = ((0.6614378277661477*bx_skin[15]-0.5590169943749473*bx_skin[13]+0.4330127018922193*bx_skin[11]-0.25*bx_skin[8])*c-1.3228756555322954*ez_skin[15]-0.8660254037844386*ez_skin[11])*jacob_skin_inv+(0.6614378277661477*ez_edge[15]-0.5590169943749473*ez_edge[13]+0.4330127018922193*ez_edge[11]-0.25*ez_edge[8])*jacob_edge_inv; 
  incr[9] = ((-(1.75*bx_skin[9])+1.4790199457749043*bx_skin[5]-1.14564392373896*bx_skin[2]+0.6614378277661477*bx_skin[0])*c-2.9580398915498085*ez_skin[5]-1.3228756555322954*ez_skin[0])*jacob_skin_inv+(1.75*ez_edge[9]-1.4790199457749043*ez_edge[5]+1.14564392373896*ez_edge[2]-0.6614378277661477*ez_edge[0])*jacob_edge_inv; 
  incr[10] = ((1.479019945774904*bx_skin[14]-1.25*bx_skin[10]+0.9682458365518543*bx_skin[6]-0.5590169943749475*bx_skin[4])*c-2.958039891549808*ez_skin[14]-1.9364916731037085*ez_skin[6])*jacob_skin_inv+(1.479019945774904*ez_edge[14]-1.25*ez_edge[10]+0.9682458365518543*ez_edge[6]-0.5590169943749475*ez_edge[4])*jacob_edge_inv; 
  incr[11] = ((-(1.14564392373896*bx_skin[15])+0.9682458365518541*bx_skin[13]-0.75*bx_skin[11]+0.4330127018922193*bx_skin[8])*c-1.9364916731037083*ez_skin[13]-0.8660254037844386*ez_skin[8])*jacob_skin_inv+(1.14564392373896*ez_edge[15]-0.9682458365518541*ez_edge[13]+0.75*ez_edge[11]-0.4330127018922193*ez_edge[8])*jacob_edge_inv; 
  incr[12] = ((-(1.75*bx_skin[12])+1.4790199457749038*bx_skin[7]-1.14564392373896*bx_skin[3]+0.6614378277661477*bx_skin[1])*c-2.9580398915498076*ez_skin[7]-1.3228756555322954*ez_skin[1])*jacob_skin_inv+(1.75*ez_edge[12]-1.4790199457749038*ez_edge[7]+1.14564392373896*ez_edge[3]-0.6614378277661477*ez_edge[1])*jacob_edge_inv; 
  incr[13] = ((1.479019945774904*bx_skin[15]-1.25*bx_skin[13]+0.9682458365518541*bx_skin[11]-0.5590169943749473*bx_skin[8])*c-2.958039891549808*ez_skin[15]-1.9364916731037083*ez_skin[11])*jacob_skin_inv+(1.479019945774904*ez_edge[15]-1.25*ez_edge[13]+0.9682458365518541*ez_edge[11]-0.5590169943749473*ez_edge[8])*jacob_edge_inv; 
  incr[14] = ((-(1.75*bx_skin[14])+1.479019945774904*bx_skin[10]-1.14564392373896*bx_skin[6]+0.6614378277661477*bx_skin[4])*c-2.958039891549808*ez_skin[10]-1.3228756555322954*ez_skin[4])*jacob_skin_inv+(1.75*ez_edge[14]-1.479019945774904*ez_edge[10]+1.14564392373896*ez_edge[6]-0.6614378277661477*ez_edge[4])*jacob_edge_inv; 
  incr[15] = ((-(1.75*bx_skin[15])+1.479019945774904*bx_skin[13]-1.14564392373896*bx_skin[11]+0.6614378277661477*bx_skin[8])*c-2.958039891549808*ez_skin[13]-1.3228756555322954*ez_skin[8])*jacob_skin_inv+(1.75*ez_edge[15]-1.479019945774904*ez_edge[13]+1.14564392373896*ez_edge[11]-0.6614378277661477*ez_edge[8])*jacob_edge_inv; 

  outBx[0] += incr[0]*dx1; 
  outBx[1] += incr[1]*dx1; 
  outBx[2] += incr[2]*dx1; 
  outBx[3] += incr[3]*dx1; 
  outBx[4] += incr[4]*dx1; 
  outBx[5] += incr[5]*dx1; 
  outBx[6] += incr[6]*dx1; 
  outBx[7] += incr[7]*dx1; 
  outBx[8] += incr[8]*dx1; 
  outBx[9] += incr[9]*dx1; 
  outBx[10] += incr[10]*dx1; 
  outBx[11] += incr[11]*dx1; 
  outBx[12] += incr[12]*dx1; 
  outBx[13] += incr[13]*dx1; 
  outBx[14] += incr[14]*dx1; 
  outBx[15] += incr[15]*dx1; 

  incr[0] = ((-(1.3228756555322954*ps_skin[9])-0.8660254037844386*ps_skin[2])*jacob_skin_inv+(0.6614378277661477*ps_edge[9]-0.5590169943749475*ps_edge[5]+0.4330127018922193*ps_edge[2]-0.25*ps_edge[0])*jacob_edge_inv)*gamma+(0.6614378277661477*by_skin[9]-0.5590169943749475*by_skin[5]+0.4330127018922193*by_skin[2]-0.25*by_skin[0])*cgamma*jacob_skin_inv; 
  incr[1] = ((-(1.3228756555322954*ps_skin[12])-0.8660254037844386*ps_skin[3])*jacob_skin_inv+(0.6614378277661477*ps_edge[12]-0.5590169943749476*ps_edge[7]+0.4330127018922193*ps_edge[3]-0.25*ps_edge[1])*jacob_edge_inv)*gamma+(0.6614378277661477*by_skin[12]-0.5590169943749476*by_skin[7]+0.4330127018922193*by_skin[3]-0.25*by_skin[1])*cgamma*jacob_skin_inv; 
  incr[2] = ((-(1.9364916731037085*ps_skin[5])-0.8660254037844386*ps_skin[0])*jacob_skin_inv+(1.14564392373896*ps_edge[9]-0.9682458365518543*ps_edge[5]+0.75*ps_edge[2]-0.4330127018922193*ps_edge[0])*jacob_edge_inv)*gamma+(-(1.14564392373896*by_skin[9])+0.9682458365518543*by_skin[5]-0.75*by_skin[2]+0.4330127018922193*by_skin[0])*cgamma*jacob_skin_inv; 
  incr[3] = ((-(1.9364916731037085*ps_skin[7])-0.8660254037844386*ps_skin[1])*jacob_skin_inv+(1.14564392373896*ps_edge[12]-0.9682458365518543*ps_edge[7]+0.75*ps_edge[3]-0.4330127018922193*ps_edge[1])*jacob_edge_inv)*gamma+(-(1.14564392373896*by_skin[12])+0.9682458365518543*by_skin[7]-0.75*by_skin[3]+0.4330127018922193*by_skin[1])*cgamma*jacob_skin_inv; 
  incr[4] = ((-(1.3228756555322954*ps_skin[14])-0.8660254037844387*ps_skin[6])*jacob_skin_inv+(0.6614378277661477*ps_edge[14]-0.5590169943749475*ps_edge[10]+0.43301270189221935*ps_edge[6]-0.25*ps_edge[4])*jacob_edge_inv)*gamma+(0.6614378277661477*by_skin[14]-0.5590169943749475*by_skin[10]+0.43301270189221935*by_skin[6]-0.25*by_skin[4])*cgamma*jacob_skin_inv; 
  incr[5] = ((-(2.9580398915498085*ps_skin[9])-1.9364916731037085*ps_skin[2])*jacob_skin_inv+(1.4790199457749043*ps_edge[9]-1.25*ps_edge[5]+0.9682458365518543*ps_edge[2]-0.5590169943749475*ps_edge[0])*jacob_edge_inv)*gamma+(1.4790199457749043*by_skin[9]-1.25*by_skin[5]+0.9682458365518543*by_skin[2]-0.5590169943749475*by_skin[0])*cgamma*jacob_skin_inv; 
  incr[6] = ((-(1.9364916731037085*ps_skin[10])-0.8660254037844387*ps_skin[4])*jacob_skin_inv+(1.14564392373896*ps_edge[14]-0.9682458365518543*ps_edge[10]+0.75*ps_edge[6]-0.43301270189221935*ps_edge[4])*jacob_edge_inv)*gamma+(-(1.14564392373896*by_skin[14])+0.9682458365518543*by_skin[10]-0.75*by_skin[6]+0.43301270189221935*by_skin[4])*cgamma*jacob_skin_inv; 
  incr[7] = ((-(2.9580398915498076*ps_skin[12])-1.9364916731037085*ps_skin[3])*jacob_skin_inv+(1.4790199457749038*ps_edge[12]-1.25*ps_edge[7]+0.9682458365518543*ps_edge[3]-0.5590169943749476*ps_edge[1])*jacob_edge_inv)*gamma+(1.4790199457749038*by_skin[12]-1.25*by_skin[7]+0.9682458365518543*by_skin[3]-0.5590169943749476*by_skin[1])*cgamma*jacob_skin_inv; 
  incr[8] = ((-(1.3228756555322954*ps_skin[15])-0.8660254037844386*ps_skin[11])*jacob_skin_inv+(0.6614378277661477*ps_edge[15]-0.5590169943749473*ps_edge[13]+0.4330127018922193*ps_edge[11]-0.25*ps_edge[8])*jacob_edge_inv)*gamma+(0.6614378277661477*by_skin[15]-0.5590169943749473*by_skin[13]+0.4330127018922193*by_skin[11]-0.25*by_skin[8])*cgamma*jacob_skin_inv; 
  incr[9] = ((-(2.9580398915498085*ps_skin[5])-1.3228756555322954*ps_skin[0])*jacob_skin_inv+(1.75*ps_edge[9]-1.4790199457749043*ps_edge[5]+1.14564392373896*ps_edge[2]-0.6614378277661477*ps_edge[0])*jacob_edge_inv)*gamma+(-(1.75*by_skin[9])+1.4790199457749043*by_skin[5]-1.14564392373896*by_skin[2]+0.6614378277661477*by_skin[0])*cgamma*jacob_skin_inv; 
  incr[10] = ((-(2.958039891549808*ps_skin[14])-1.9364916731037085*ps_skin[6])*jacob_skin_inv+(1.479019945774904*ps_edge[14]-1.25*ps_edge[10]+0.9682458365518543*ps_edge[6]-0.5590169943749475*ps_edge[4])*jacob_edge_inv)*gamma+(1.479019945774904*by_skin[14]-1.25*by_skin[10]+0.9682458365518543*by_skin[6]-0.5590169943749475*by_skin[4])*cgamma*jacob_skin_inv; 
  incr[11] = ((-(1.9364916731037083*ps_skin[13])-0.8660254037844386*ps_skin[8])*jacob_skin_inv+(1.14564392373896*ps_edge[15]-0.9682458365518541*ps_edge[13]+0.75*ps_edge[11]-0.4330127018922193*ps_edge[8])*jacob_edge_inv)*gamma+(-(1.14564392373896*by_skin[15])+0.9682458365518541*by_skin[13]-0.75*by_skin[11]+0.4330127018922193*by_skin[8])*cgamma*jacob_skin_inv; 
  incr[12] = ((-(2.9580398915498076*ps_skin[7])-1.3228756555322954*ps_skin[1])*jacob_skin_inv+(1.75*ps_edge[12]-1.4790199457749038*ps_edge[7]+1.14564392373896*ps_edge[3]-0.6614378277661477*ps_edge[1])*jacob_edge_inv)*gamma+(-(1.75*by_skin[12])+1.4790199457749038*by_skin[7]-1.14564392373896*by_skin[3]+0.6614378277661477*by_skin[1])*cgamma*jacob_skin_inv; 
  incr[13] = ((-(2.958039891549808*ps_skin[15])-1.9364916731037083*ps_skin[11])*jacob_skin_inv+(1.479019945774904*ps_edge[15]-1.25*ps_edge[13]+0.9682458365518541*ps_edge[11]-0.5590169943749473*ps_edge[8])*jacob_edge_inv)*gamma+(1.479019945774904*by_skin[15]-1.25*by_skin[13]+0.9682458365518541*by_skin[11]-0.5590169943749473*by_skin[8])*cgamma*jacob_skin_inv; 
  incr[14] = ((-(2.958039891549808*ps_skin[10])-1.3228756555322954*ps_skin[4])*jacob_skin_inv+(1.75*ps_edge[14]-1.479019945774904*ps_edge[10]+1.14564392373896*ps_edge[6]-0.6614378277661477*ps_edge[4])*jacob_edge_inv)*gamma+(-(1.75*by_skin[14])+1.479019945774904*by_skin[10]-1.14564392373896*by_skin[6]+0.6614378277661477*by_skin[4])*cgamma*jacob_skin_inv; 
  incr[15] = ((-(2.958039891549808*ps_skin[13])-1.3228756555322954*ps_skin[8])*jacob_skin_inv+(1.75*ps_edge[15]-1.479019945774904*ps_edge[13]+1.14564392373896*ps_edge[11]-0.6614378277661477*ps_edge[8])*jacob_edge_inv)*gamma+(-(1.75*by_skin[15])+1.479019945774904*by_skin[13]-1.14564392373896*by_skin[11]+0.6614378277661477*by_skin[8])*cgamma*jacob_skin_inv; 

  outBy[0] += incr[0]*dx1; 
  outBy[1] += incr[1]*dx1; 
  outBy[2] += incr[2]*dx1; 
  outBy[3] += incr[3]*dx1; 
  outBy[4] += incr[4]*dx1; 
  outBy[5] += incr[5]*dx1; 
  outBy[6] += incr[6]*dx1; 
  outBy[7] += incr[7]*dx1; 
  outBy[8] += incr[8]*dx1; 
  outBy[9] += incr[9]*dx1; 
  outBy[10] += incr[10]*dx1; 
  outBy[11] += incr[11]*dx1; 
  outBy[12] += incr[12]*dx1; 
  outBy[13] += incr[13]*dx1; 
  outBy[14] += incr[14]*dx1; 
  outBy[15] += incr[15]*dx1; 

  incr[0] = ((0.6614378277661477*bz_skin[9]-0.5590169943749475*bz_skin[5]+0.4330127018922193*bz_skin[2]-0.25*bz_skin[0])*c+1.3228756555322954*ex_skin[9]+0.8660254037844386*ex_skin[2])*jacob_skin_inv+(-(0.6614378277661477*ex_edge[9])+0.5590169943749475*ex_edge[5]-0.4330127018922193*ex_edge[2]+0.25*ex_edge[0])*jacob_edge_inv; 
  incr[1] = ((0.6614378277661477*bz_skin[12]-0.5590169943749476*bz_skin[7]+0.4330127018922193*bz_skin[3]-0.25*bz_skin[1])*c+1.3228756555322954*ex_skin[12]+0.8660254037844386*ex_skin[3])*jacob_skin_inv+(-(0.6614378277661477*ex_edge[12])+0.5590169943749476*ex_edge[7]-0.4330127018922193*ex_edge[3]+0.25*ex_edge[1])*jacob_edge_inv; 
  incr[2] = ((-(1.14564392373896*bz_skin[9])+0.9682458365518543*bz_skin[5]-0.75*bz_skin[2]+0.4330127018922193*bz_skin[0])*c+1.9364916731037085*ex_skin[5]+0.8660254037844386*ex_skin[0])*jacob_skin_inv+(-(1.14564392373896*ex_edge[9])+0.9682458365518543*ex_edge[5]-0.75*ex_edge[2]+0.4330127018922193*ex_edge[0])*jacob_edge_inv; 
  incr[3] = ((-(1.14564392373896*bz_skin[12])+0.9682458365518543*bz_skin[7]-0.75*bz_skin[3]+0.4330127018922193*bz_skin[1])*c+1.9364916731037085*ex_skin[7]+0.8660254037844386*ex_skin[1])*jacob_skin_inv+(-(1.14564392373896*ex_edge[12])+0.9682458365518543*ex_edge[7]-0.75*ex_edge[3]+0.4330127018922193*ex_edge[1])*jacob_edge_inv; 
  incr[4] = ((0.6614378277661477*bz_skin[14]-0.5590169943749475*bz_skin[10]+0.43301270189221935*bz_skin[6]-0.25*bz_skin[4])*c+1.3228756555322954*ex_skin[14]+0.8660254037844387*ex_skin[6])*jacob_skin_inv+(-(0.6614378277661477*ex_edge[14])+0.5590169943749475*ex_edge[10]-0.43301270189221935*ex_edge[6]+0.25*ex_edge[4])*jacob_edge_inv; 
  incr[5] = ((1.4790199457749043*bz_skin[9]-1.25*bz_skin[5]+0.9682458365518543*bz_skin[2]-0.5590169943749475*bz_skin[0])*c+2.9580398915498085*ex_skin[9]+1.9364916731037085*ex_skin[2])*jacob_skin_inv+(-(1.4790199457749043*ex_edge[9])+1.25*ex_edge[5]-0.9682458365518543*ex_edge[2]+0.5590169943749475*ex_edge[0])*jacob_edge_inv; 
  incr[6] = ((-(1.14564392373896*bz_skin[14])+0.9682458365518543*bz_skin[10]-0.75*bz_skin[6]+0.43301270189221935*bz_skin[4])*c+1.9364916731037085*ex_skin[10]+0.8660254037844387*ex_skin[4])*jacob_skin_inv+(-(1.14564392373896*ex_edge[14])+0.9682458365518543*ex_edge[10]-0.75*ex_edge[6]+0.43301270189221935*ex_edge[4])*jacob_edge_inv; 
  incr[7] = ((1.4790199457749038*bz_skin[12]-1.25*bz_skin[7]+0.9682458365518543*bz_skin[3]-0.5590169943749476*bz_skin[1])*c+2.9580398915498076*ex_skin[12]+1.9364916731037085*ex_skin[3])*jacob_skin_inv+(-(1.4790199457749038*ex_edge[12])+1.25*ex_edge[7]-0.9682458365518543*ex_edge[3]+0.5590169943749476*ex_edge[1])*jacob_edge_inv; 
  incr[8] = ((0.6614378277661477*bz_skin[15]-0.5590169943749473*bz_skin[13]+0.4330127018922193*bz_skin[11]-0.25*bz_skin[8])*c+1.3228756555322954*ex_skin[15]+0.8660254037844386*ex_skin[11])*jacob_skin_inv+(-(0.6614378277661477*ex_edge[15])+0.5590169943749473*ex_edge[13]-0.4330127018922193*ex_edge[11]+0.25*ex_edge[8])*jacob_edge_inv; 
  incr[9] = ((-(1.75*bz_skin[9])+1.4790199457749043*bz_skin[5]-1.14564392373896*bz_skin[2]+0.6614378277661477*bz_skin[0])*c+2.9580398915498085*ex_skin[5]+1.3228756555322954*ex_skin[0])*jacob_skin_inv+(-(1.75*ex_edge[9])+1.4790199457749043*ex_edge[5]-1.14564392373896*ex_edge[2]+0.6614378277661477*ex_edge[0])*jacob_edge_inv; 
  incr[10] = ((1.479019945774904*bz_skin[14]-1.25*bz_skin[10]+0.9682458365518543*bz_skin[6]-0.5590169943749475*bz_skin[4])*c+2.958039891549808*ex_skin[14]+1.9364916731037085*ex_skin[6])*jacob_skin_inv+(-(1.479019945774904*ex_edge[14])+1.25*ex_edge[10]-0.9682458365518543*ex_edge[6]+0.5590169943749475*ex_edge[4])*jacob_edge_inv; 
  incr[11] = ((-(1.14564392373896*bz_skin[15])+0.9682458365518541*bz_skin[13]-0.75*bz_skin[11]+0.4330127018922193*bz_skin[8])*c+1.9364916731037083*ex_skin[13]+0.8660254037844386*ex_skin[8])*jacob_skin_inv+(-(1.14564392373896*ex_edge[15])+0.9682458365518541*ex_edge[13]-0.75*ex_edge[11]+0.4330127018922193*ex_edge[8])*jacob_edge_inv; 
  incr[12] = ((-(1.75*bz_skin[12])+1.4790199457749038*bz_skin[7]-1.14564392373896*bz_skin[3]+0.6614378277661477*bz_skin[1])*c+2.9580398915498076*ex_skin[7]+1.3228756555322954*ex_skin[1])*jacob_skin_inv+(-(1.75*ex_edge[12])+1.4790199457749038*ex_edge[7]-1.14564392373896*ex_edge[3]+0.6614378277661477*ex_edge[1])*jacob_edge_inv; 
  incr[13] = ((1.479019945774904*bz_skin[15]-1.25*bz_skin[13]+0.9682458365518541*bz_skin[11]-0.5590169943749473*bz_skin[8])*c+2.958039891549808*ex_skin[15]+1.9364916731037083*ex_skin[11])*jacob_skin_inv+(-(1.479019945774904*ex_edge[15])+1.25*ex_edge[13]-0.9682458365518541*ex_edge[11]+0.5590169943749473*ex_edge[8])*jacob_edge_inv; 
  incr[14] = ((-(1.75*bz_skin[14])+1.479019945774904*bz_skin[10]-1.14564392373896*bz_skin[6]+0.6614378277661477*bz_skin[4])*c+2.958039891549808*ex_skin[10]+1.3228756555322954*ex_skin[4])*jacob_skin_inv+(-(1.75*ex_edge[14])+1.479019945774904*ex_edge[10]-1.14564392373896*ex_edge[6]+0.6614378277661477*ex_edge[4])*jacob_edge_inv; 
  incr[15] = ((-(1.75*bz_skin[15])+1.479019945774904*bz_skin[13]-1.14564392373896*bz_skin[11]+0.6614378277661477*bz_skin[8])*c+2.958039891549808*ex_skin[13]+1.3228756555322954*ex_skin[8])*jacob_skin_inv+(-(1.75*ex_edge[15])+1.479019945774904*ex_edge[13]-1.14564392373896*ex_edge[11]+0.6614378277661477*ex_edge[8])*jacob_edge_inv; 

  outBz[0] += incr[0]*dx1; 
  outBz[1] += incr[1]*dx1; 
  outBz[2] += incr[2]*dx1; 
  outBz[3] += incr[3]*dx1; 
  outBz[4] += incr[4]*dx1; 
  outBz[5] += incr[5]*dx1; 
  outBz[6] += incr[6]*dx1; 
  outBz[7] += incr[7]*dx1; 
  outBz[8] += incr[8]*dx1; 
  outBz[9] += incr[9]*dx1; 
  outBz[10] += incr[10]*dx1; 
  outBz[11] += incr[11]*dx1; 
  outBz[12] += incr[12]*dx1; 
  outBz[13] += incr[13]*dx1; 
  outBz[14] += incr[14]*dx1; 
  outBz[15] += incr[15]*dx1; 

  incr[0] = (0.6614378277661477*ph_skin[9]-0.5590169943749475*ph_skin[5]+0.4330127018922193*ph_skin[2]-0.25*ph_skin[0])*cchi*jacob_skin_inv+chi*((-(1.3228756555322954*ey_skin[9])-0.8660254037844386*ey_skin[2])*jacob_skin_inv+(0.6614378277661477*ey_edge[9]-0.5590169943749475*ey_edge[5]+0.4330127018922193*ey_edge[2]-0.25*ey_edge[0])*jacob_edge_inv); 
  incr[1] = (0.6614378277661477*ph_skin[12]-0.5590169943749476*ph_skin[7]+0.4330127018922193*ph_skin[3]-0.25*ph_skin[1])*cchi*jacob_skin_inv+chi*((-(1.3228756555322954*ey_skin[12])-0.8660254037844386*ey_skin[3])*jacob_skin_inv+(0.6614378277661477*ey_edge[12]-0.5590169943749476*ey_edge[7]+0.4330127018922193*ey_edge[3]-0.25*ey_edge[1])*jacob_edge_inv); 
  incr[2] = (-(1.14564392373896*ph_skin[9])+0.9682458365518543*ph_skin[5]-0.75*ph_skin[2]+0.4330127018922193*ph_skin[0])*cchi*jacob_skin_inv+chi*((-(1.9364916731037085*ey_skin[5])-0.8660254037844386*ey_skin[0])*jacob_skin_inv+(1.14564392373896*ey_edge[9]-0.9682458365518543*ey_edge[5]+0.75*ey_edge[2]-0.4330127018922193*ey_edge[0])*jacob_edge_inv); 
  incr[3] = (-(1.14564392373896*ph_skin[12])+0.9682458365518543*ph_skin[7]-0.75*ph_skin[3]+0.4330127018922193*ph_skin[1])*cchi*jacob_skin_inv+chi*((-(1.9364916731037085*ey_skin[7])-0.8660254037844386*ey_skin[1])*jacob_skin_inv+(1.14564392373896*ey_edge[12]-0.9682458365518543*ey_edge[7]+0.75*ey_edge[3]-0.4330127018922193*ey_edge[1])*jacob_edge_inv); 
  incr[4] = (0.6614378277661477*ph_skin[14]-0.5590169943749475*ph_skin[10]+0.43301270189221935*ph_skin[6]-0.25*ph_skin[4])*cchi*jacob_skin_inv+chi*((-(1.3228756555322954*ey_skin[14])-0.8660254037844387*ey_skin[6])*jacob_skin_inv+(0.6614378277661477*ey_edge[14]-0.5590169943749475*ey_edge[10]+0.43301270189221935*ey_edge[6]-0.25*ey_edge[4])*jacob_edge_inv); 
  incr[5] = (1.4790199457749043*ph_skin[9]-1.25*ph_skin[5]+0.9682458365518543*ph_skin[2]-0.5590169943749475*ph_skin[0])*cchi*jacob_skin_inv+chi*((-(2.9580398915498085*ey_skin[9])-1.9364916731037085*ey_skin[2])*jacob_skin_inv+(1.4790199457749043*ey_edge[9]-1.25*ey_edge[5]+0.9682458365518543*ey_edge[2]-0.5590169943749475*ey_edge[0])*jacob_edge_inv); 
  incr[6] = (-(1.14564392373896*ph_skin[14])+0.9682458365518543*ph_skin[10]-0.75*ph_skin[6]+0.43301270189221935*ph_skin[4])*cchi*jacob_skin_inv+chi*((-(1.9364916731037085*ey_skin[10])-0.8660254037844387*ey_skin[4])*jacob_skin_inv+(1.14564392373896*ey_edge[14]-0.9682458365518543*ey_edge[10]+0.75*ey_edge[6]-0.43301270189221935*ey_edge[4])*jacob_edge_inv); 
  incr[7] = (1.4790199457749038*ph_skin[12]-1.25*ph_skin[7]+0.9682458365518543*ph_skin[3]-0.5590169943749476*ph_skin[1])*cchi*jacob_skin_inv+chi*((-(2.9580398915498076*ey_skin[12])-1.9364916731037085*ey_skin[3])*jacob_skin_inv+(1.4790199457749038*ey_edge[12]-1.25*ey_edge[7]+0.9682458365518543*ey_edge[3]-0.5590169943749476*ey_edge[1])*jacob_edge_inv); 
  incr[8] = (0.6614378277661477*ph_skin[15]-0.5590169943749473*ph_skin[13]+0.4330127018922193*ph_skin[11]-0.25*ph_skin[8])*cchi*jacob_skin_inv+chi*((-(1.3228756555322954*ey_skin[15])-0.8660254037844386*ey_skin[11])*jacob_skin_inv+(0.6614378277661477*ey_edge[15]-0.5590169943749473*ey_edge[13]+0.4330127018922193*ey_edge[11]-0.25*ey_edge[8])*jacob_edge_inv); 
  incr[9] = (-(1.75*ph_skin[9])+1.4790199457749043*ph_skin[5]-1.14564392373896*ph_skin[2]+0.6614378277661477*ph_skin[0])*cchi*jacob_skin_inv+chi*((-(2.9580398915498085*ey_skin[5])-1.3228756555322954*ey_skin[0])*jacob_skin_inv+(1.75*ey_edge[9]-1.4790199457749043*ey_edge[5]+1.14564392373896*ey_edge[2]-0.6614378277661477*ey_edge[0])*jacob_edge_inv); 
  incr[10] = (1.479019945774904*ph_skin[14]-1.25*ph_skin[10]+0.9682458365518543*ph_skin[6]-0.5590169943749475*ph_skin[4])*cchi*jacob_skin_inv+chi*((-(2.958039891549808*ey_skin[14])-1.9364916731037085*ey_skin[6])*jacob_skin_inv+(1.479019945774904*ey_edge[14]-1.25*ey_edge[10]+0.9682458365518543*ey_edge[6]-0.5590169943749475*ey_edge[4])*jacob_edge_inv); 
  incr[11] = (-(1.14564392373896*ph_skin[15])+0.9682458365518541*ph_skin[13]-0.75*ph_skin[11]+0.4330127018922193*ph_skin[8])*cchi*jacob_skin_inv+chi*((-(1.9364916731037083*ey_skin[13])-0.8660254037844386*ey_skin[8])*jacob_skin_inv+(1.14564392373896*ey_edge[15]-0.9682458365518541*ey_edge[13]+0.75*ey_edge[11]-0.4330127018922193*ey_edge[8])*jacob_edge_inv); 
  incr[12] = (-(1.75*ph_skin[12])+1.4790199457749038*ph_skin[7]-1.14564392373896*ph_skin[3]+0.6614378277661477*ph_skin[1])*cchi*jacob_skin_inv+chi*((-(2.9580398915498076*ey_skin[7])-1.3228756555322954*ey_skin[1])*jacob_skin_inv+(1.75*ey_edge[12]-1.4790199457749038*ey_edge[7]+1.14564392373896*ey_edge[3]-0.6614378277661477*ey_edge[1])*jacob_edge_inv); 
  incr[13] = (1.479019945774904*ph_skin[15]-1.25*ph_skin[13]+0.9682458365518541*ph_skin[11]-0.5590169943749473*ph_skin[8])*cchi*jacob_skin_inv+chi*((-(2.958039891549808*ey_skin[15])-1.9364916731037083*ey_skin[11])*jacob_skin_inv+(1.479019945774904*ey_edge[15]-1.25*ey_edge[13]+0.9682458365518541*ey_edge[11]-0.5590169943749473*ey_edge[8])*jacob_edge_inv); 
  incr[14] = (-(1.75*ph_skin[14])+1.479019945774904*ph_skin[10]-1.14564392373896*ph_skin[6]+0.6614378277661477*ph_skin[4])*cchi*jacob_skin_inv+chi*((-(2.958039891549808*ey_skin[10])-1.3228756555322954*ey_skin[4])*jacob_skin_inv+(1.75*ey_edge[14]-1.479019945774904*ey_edge[10]+1.14564392373896*ey_edge[6]-0.6614378277661477*ey_edge[4])*jacob_edge_inv); 
  incr[15] = (-(1.75*ph_skin[15])+1.479019945774904*ph_skin[13]-1.14564392373896*ph_skin[11]+0.6614378277661477*ph_skin[8])*cchi*jacob_skin_inv+chi*((-(2.958039891549808*ey_skin[13])-1.3228756555322954*ey_skin[8])*jacob_skin_inv+(1.75*ey_edge[15]-1.479019945774904*ey_edge[13]+1.14564392373896*ey_edge[11]-0.6614378277661477*ey_edge[8])*jacob_edge_inv); 

  outPh[0] += incr[0]*dx1; 
  outPh[1] += incr[1]*dx1; 
  outPh[2] += incr[2]*dx1; 
  outPh[3] += incr[3]*dx1; 
  outPh[4] += incr[4]*dx1; 
  outPh[5] += incr[5]*dx1; 
  outPh[6] += incr[6]*dx1; 
  outPh[7] += incr[7]*dx1; 
  outPh[8] += incr[8]*dx1; 
  outPh[9] += incr[9]*dx1; 
  outPh[10] += incr[10]*dx1; 
  outPh[11] += incr[11]*dx1; 
  outPh[12] += incr[12]*dx1; 
  outPh[13] += incr[13]*dx1; 
  outPh[14] += incr[14]*dx1; 
  outPh[15] += incr[15]*dx1; 

  incr[0] = (0.6614378277661477*ps_skin[9]-0.5590169943749475*ps_skin[5]+0.4330127018922193*ps_skin[2]-0.25*ps_skin[0])*cgamma*jacob_skin_inv+c2gamma*((-(1.3228756555322954*by_skin[9])-0.8660254037844386*by_skin[2])*jacob_skin_inv+(0.6614378277661477*by_edge[9]-0.5590169943749475*by_edge[5]+0.4330127018922193*by_edge[2]-0.25*by_edge[0])*jacob_edge_inv); 
  incr[1] = (0.6614378277661477*ps_skin[12]-0.5590169943749476*ps_skin[7]+0.4330127018922193*ps_skin[3]-0.25*ps_skin[1])*cgamma*jacob_skin_inv+c2gamma*((-(1.3228756555322954*by_skin[12])-0.8660254037844386*by_skin[3])*jacob_skin_inv+(0.6614378277661477*by_edge[12]-0.5590169943749476*by_edge[7]+0.4330127018922193*by_edge[3]-0.25*by_edge[1])*jacob_edge_inv); 
  incr[2] = (-(1.14564392373896*ps_skin[9])+0.9682458365518543*ps_skin[5]-0.75*ps_skin[2]+0.4330127018922193*ps_skin[0])*cgamma*jacob_skin_inv+c2gamma*((-(1.9364916731037085*by_skin[5])-0.8660254037844386*by_skin[0])*jacob_skin_inv+(1.14564392373896*by_edge[9]-0.9682458365518543*by_edge[5]+0.75*by_edge[2]-0.4330127018922193*by_edge[0])*jacob_edge_inv); 
  incr[3] = (-(1.14564392373896*ps_skin[12])+0.9682458365518543*ps_skin[7]-0.75*ps_skin[3]+0.4330127018922193*ps_skin[1])*cgamma*jacob_skin_inv+c2gamma*((-(1.9364916731037085*by_skin[7])-0.8660254037844386*by_skin[1])*jacob_skin_inv+(1.14564392373896*by_edge[12]-0.9682458365518543*by_edge[7]+0.75*by_edge[3]-0.4330127018922193*by_edge[1])*jacob_edge_inv); 
  incr[4] = (0.6614378277661477*ps_skin[14]-0.5590169943749475*ps_skin[10]+0.43301270189221935*ps_skin[6]-0.25*ps_skin[4])*cgamma*jacob_skin_inv+c2gamma*((-(1.3228756555322954*by_skin[14])-0.8660254037844387*by_skin[6])*jacob_skin_inv+(0.6614378277661477*by_edge[14]-0.5590169943749475*by_edge[10]+0.43301270189221935*by_edge[6]-0.25*by_edge[4])*jacob_edge_inv); 
  incr[5] = (1.4790199457749043*ps_skin[9]-1.25*ps_skin[5]+0.9682458365518543*ps_skin[2]-0.5590169943749475*ps_skin[0])*cgamma*jacob_skin_inv+c2gamma*((-(2.9580398915498085*by_skin[9])-1.9364916731037085*by_skin[2])*jacob_skin_inv+(1.4790199457749043*by_edge[9]-1.25*by_edge[5]+0.9682458365518543*by_edge[2]-0.5590169943749475*by_edge[0])*jacob_edge_inv); 
  incr[6] = (-(1.14564392373896*ps_skin[14])+0.9682458365518543*ps_skin[10]-0.75*ps_skin[6]+0.43301270189221935*ps_skin[4])*cgamma*jacob_skin_inv+c2gamma*((-(1.9364916731037085*by_skin[10])-0.8660254037844387*by_skin[4])*jacob_skin_inv+(1.14564392373896*by_edge[14]-0.9682458365518543*by_edge[10]+0.75*by_edge[6]-0.43301270189221935*by_edge[4])*jacob_edge_inv); 
  incr[7] = (1.4790199457749038*ps_skin[12]-1.25*ps_skin[7]+0.9682458365518543*ps_skin[3]-0.5590169943749476*ps_skin[1])*cgamma*jacob_skin_inv+c2gamma*((-(2.9580398915498076*by_skin[12])-1.9364916731037085*by_skin[3])*jacob_skin_inv+(1.4790199457749038*by_edge[12]-1.25*by_edge[7]+0.9682458365518543*by_edge[3]-0.5590169943749476*by_edge[1])*jacob_edge_inv); 
  incr[8] = (0.6614378277661477*ps_skin[15]-0.5590169943749473*ps_skin[13]+0.4330127018922193*ps_skin[11]-0.25*ps_skin[8])*cgamma*jacob_skin_inv+c2gamma*((-(1.3228756555322954*by_skin[15])-0.8660254037844386*by_skin[11])*jacob_skin_inv+(0.6614378277661477*by_edge[15]-0.5590169943749473*by_edge[13]+0.4330127018922193*by_edge[11]-0.25*by_edge[8])*jacob_edge_inv); 
  incr[9] = (-(1.75*ps_skin[9])+1.4790199457749043*ps_skin[5]-1.14564392373896*ps_skin[2]+0.6614378277661477*ps_skin[0])*cgamma*jacob_skin_inv+c2gamma*((-(2.9580398915498085*by_skin[5])-1.3228756555322954*by_skin[0])*jacob_skin_inv+(1.75*by_edge[9]-1.4790199457749043*by_edge[5]+1.14564392373896*by_edge[2]-0.6614378277661477*by_edge[0])*jacob_edge_inv); 
  incr[10] = (1.479019945774904*ps_skin[14]-1.25*ps_skin[10]+0.9682458365518543*ps_skin[6]-0.5590169943749475*ps_skin[4])*cgamma*jacob_skin_inv+c2gamma*((-(2.958039891549808*by_skin[14])-1.9364916731037085*by_skin[6])*jacob_skin_inv+(1.479019945774904*by_edge[14]-1.25*by_edge[10]+0.9682458365518543*by_edge[6]-0.5590169943749475*by_edge[4])*jacob_edge_inv); 
  incr[11] = (-(1.14564392373896*ps_skin[15])+0.9682458365518541*ps_skin[13]-0.75*ps_skin[11]+0.4330127018922193*ps_skin[8])*cgamma*jacob_skin_inv+c2gamma*((-(1.9364916731037083*by_skin[13])-0.8660254037844386*by_skin[8])*jacob_skin_inv+(1.14564392373896*by_edge[15]-0.9682458365518541*by_edge[13]+0.75*by_edge[11]-0.4330127018922193*by_edge[8])*jacob_edge_inv); 
  incr[12] = (-(1.75*ps_skin[12])+1.4790199457749038*ps_skin[7]-1.14564392373896*ps_skin[3]+0.6614378277661477*ps_skin[1])*cgamma*jacob_skin_inv+c2gamma*((-(2.9580398915498076*by_skin[7])-1.3228756555322954*by_skin[1])*jacob_skin_inv+(1.75*by_edge[12]-1.4790199457749038*by_edge[7]+1.14564392373896*by_edge[3]-0.6614378277661477*by_edge[1])*jacob_edge_inv); 
  incr[13] = (1.479019945774904*ps_skin[15]-1.25*ps_skin[13]+0.9682458365518541*ps_skin[11]-0.5590169943749473*ps_skin[8])*cgamma*jacob_skin_inv+c2gamma*((-(2.958039891549808*by_skin[15])-1.9364916731037083*by_skin[11])*jacob_skin_inv+(1.479019945774904*by_edge[15]-1.25*by_edge[13]+0.9682458365518541*by_edge[11]-0.5590169943749473*by_edge[8])*jacob_edge_inv); 
  incr[14] = (-(1.75*ps_skin[14])+1.479019945774904*ps_skin[10]-1.14564392373896*ps_skin[6]+0.6614378277661477*ps_skin[4])*cgamma*jacob_skin_inv+c2gamma*((-(2.958039891549808*by_skin[10])-1.3228756555322954*by_skin[4])*jacob_skin_inv+(1.75*by_edge[14]-1.479019945774904*by_edge[10]+1.14564392373896*by_edge[6]-0.6614378277661477*by_edge[4])*jacob_edge_inv); 
  incr[15] = (-(1.75*ps_skin[15])+1.479019945774904*ps_skin[13]-1.14564392373896*ps_skin[11]+0.6614378277661477*ps_skin[8])*cgamma*jacob_skin_inv+c2gamma*((-(2.958039891549808*by_skin[13])-1.3228756555322954*by_skin[8])*jacob_skin_inv+(1.75*by_edge[15]-1.479019945774904*by_edge[13]+1.14564392373896*by_edge[11]-0.6614378277661477*by_edge[8])*jacob_edge_inv); 

  outPs[0] += incr[0]*dx1; 
  outPs[1] += incr[1]*dx1; 
  outPs[2] += incr[2]*dx1; 
  outPs[3] += incr[3]*dx1; 
  outPs[4] += incr[4]*dx1; 
  outPs[5] += incr[5]*dx1; 
  outPs[6] += incr[6]*dx1; 
  outPs[7] += incr[7]*dx1; 
  outPs[8] += incr[8]*dx1; 
  outPs[9] += incr[9]*dx1; 
  outPs[10] += incr[10]*dx1; 
  outPs[11] += incr[11]*dx1; 
  outPs[12] += incr[12]*dx1; 
  outPs[13] += incr[13]*dx1; 
  outPs[14] += incr[14]*dx1; 
  outPs[15] += incr[15]*dx1; 

  } else { 

  incr[0] = (-(0.6614378277661477*ex_skin[9])-0.5590169943749475*ex_skin[5]-0.4330127018922193*ex_skin[2]-0.25*ex_skin[0])*c*jacob_skin_inv+c2*((1.3228756555322954*bz_skin[9]+0.8660254037844386*bz_skin[2])*jacob_skin_inv+(-(0.6614378277661477*bz_edge[9])-0.5590169943749475*bz_edge[5]-0.4330127018922193*bz_edge[2]-0.25*bz_edge[0])*jacob_edge_inv); 
  incr[1] = (-(0.6614378277661477*ex_skin[12])-0.5590169943749476*ex_skin[7]-0.4330127018922193*ex_skin[3]-0.25*ex_skin[1])*c*jacob_skin_inv+c2*((1.3228756555322954*bz_skin[12]+0.8660254037844386*bz_skin[3])*jacob_skin_inv+(-(0.6614378277661477*bz_edge[12])-0.5590169943749476*bz_edge[7]-0.4330127018922193*bz_edge[3]-0.25*bz_edge[1])*jacob_edge_inv); 
  incr[2] = (-(1.14564392373896*ex_skin[9])-0.9682458365518543*ex_skin[5]-0.75*ex_skin[2]-0.4330127018922193*ex_skin[0])*c*jacob_skin_inv+c2*((1.9364916731037085*bz_skin[5]+0.8660254037844386*bz_skin[0])*jacob_skin_inv+(1.14564392373896*bz_edge[9]+0.9682458365518543*bz_edge[5]+0.75*bz_edge[2]+0.4330127018922193*bz_edge[0])*jacob_edge_inv); 
  incr[3] = (-(1.14564392373896*ex_skin[12])-0.9682458365518543*ex_skin[7]-0.75*ex_skin[3]-0.4330127018922193*ex_skin[1])*c*jacob_skin_inv+c2*((1.9364916731037085*bz_skin[7]+0.8660254037844386*bz_skin[1])*jacob_skin_inv+(1.14564392373896*bz_edge[12]+0.9682458365518543*bz_edge[7]+0.75*bz_edge[3]+0.4330127018922193*bz_edge[1])*jacob_edge_inv); 
  incr[4] = (-(0.6614378277661477*ex_skin[14])-0.5590169943749475*ex_skin[10]-0.43301270189221935*ex_skin[6]-0.25*ex_skin[4])*c*jacob_skin_inv+c2*((1.3228756555322954*bz_skin[14]+0.8660254037844387*bz_skin[6])*jacob_skin_inv+(-(0.6614378277661477*bz_edge[14])-0.5590169943749475*bz_edge[10]-0.43301270189221935*bz_edge[6]-0.25*bz_edge[4])*jacob_edge_inv); 
  incr[5] = (-(1.4790199457749043*ex_skin[9])-1.25*ex_skin[5]-0.9682458365518543*ex_skin[2]-0.5590169943749475*ex_skin[0])*c*jacob_skin_inv+c2*((2.9580398915498085*bz_skin[9]+1.9364916731037085*bz_skin[2])*jacob_skin_inv+(-(1.4790199457749043*bz_edge[9])-1.25*bz_edge[5]-0.9682458365518543*bz_edge[2]-0.5590169943749475*bz_edge[0])*jacob_edge_inv); 
  incr[6] = (-(1.14564392373896*ex_skin[14])-0.9682458365518543*ex_skin[10]-0.75*ex_skin[6]-0.43301270189221935*ex_skin[4])*c*jacob_skin_inv+c2*((1.9364916731037085*bz_skin[10]+0.8660254037844387*bz_skin[4])*jacob_skin_inv+(1.14564392373896*bz_edge[14]+0.9682458365518543*bz_edge[10]+0.75*bz_edge[6]+0.43301270189221935*bz_edge[4])*jacob_edge_inv); 
  incr[7] = (-(1.4790199457749038*ex_skin[12])-1.25*ex_skin[7]-0.9682458365518543*ex_skin[3]-0.5590169943749476*ex_skin[1])*c*jacob_skin_inv+c2*((2.9580398915498076*bz_skin[12]+1.9364916731037085*bz_skin[3])*jacob_skin_inv+(-(1.4790199457749038*bz_edge[12])-1.25*bz_edge[7]-0.9682458365518543*bz_edge[3]-0.5590169943749476*bz_edge[1])*jacob_edge_inv); 
  incr[8] = (-(0.6614378277661477*ex_skin[15])-0.5590169943749473*ex_skin[13]-0.4330127018922193*ex_skin[11]-0.25*ex_skin[8])*c*jacob_skin_inv+c2*((1.3228756555322954*bz_skin[15]+0.8660254037844386*bz_skin[11])*jacob_skin_inv+(-(0.6614378277661477*bz_edge[15])-0.5590169943749473*bz_edge[13]-0.4330127018922193*bz_edge[11]-0.25*bz_edge[8])*jacob_edge_inv); 
  incr[9] = (-(1.75*ex_skin[9])-1.4790199457749043*ex_skin[5]-1.14564392373896*ex_skin[2]-0.6614378277661477*ex_skin[0])*c*jacob_skin_inv+c2*((2.9580398915498085*bz_skin[5]+1.3228756555322954*bz_skin[0])*jacob_skin_inv+(1.75*bz_edge[9]+1.4790199457749043*bz_edge[5]+1.14564392373896*bz_edge[2]+0.6614378277661477*bz_edge[0])*jacob_edge_inv); 
  incr[10] = (-(1.479019945774904*ex_skin[14])-1.25*ex_skin[10]-0.9682458365518543*ex_skin[6]-0.5590169943749475*ex_skin[4])*c*jacob_skin_inv+c2*((2.958039891549808*bz_skin[14]+1.9364916731037085*bz_skin[6])*jacob_skin_inv+(-(1.479019945774904*bz_edge[14])-1.25*bz_edge[10]-0.9682458365518543*bz_edge[6]-0.5590169943749475*bz_edge[4])*jacob_edge_inv); 
  incr[11] = (-(1.14564392373896*ex_skin[15])-0.9682458365518541*ex_skin[13]-0.75*ex_skin[11]-0.4330127018922193*ex_skin[8])*c*jacob_skin_inv+c2*((1.9364916731037083*bz_skin[13]+0.8660254037844386*bz_skin[8])*jacob_skin_inv+(1.14564392373896*bz_edge[15]+0.9682458365518541*bz_edge[13]+0.75*bz_edge[11]+0.4330127018922193*bz_edge[8])*jacob_edge_inv); 
  incr[12] = (-(1.75*ex_skin[12])-1.4790199457749038*ex_skin[7]-1.14564392373896*ex_skin[3]-0.6614378277661477*ex_skin[1])*c*jacob_skin_inv+c2*((2.9580398915498076*bz_skin[7]+1.3228756555322954*bz_skin[1])*jacob_skin_inv+(1.75*bz_edge[12]+1.4790199457749038*bz_edge[7]+1.14564392373896*bz_edge[3]+0.6614378277661477*bz_edge[1])*jacob_edge_inv); 
  incr[13] = (-(1.479019945774904*ex_skin[15])-1.25*ex_skin[13]-0.9682458365518541*ex_skin[11]-0.5590169943749473*ex_skin[8])*c*jacob_skin_inv+c2*((2.958039891549808*bz_skin[15]+1.9364916731037083*bz_skin[11])*jacob_skin_inv+(-(1.479019945774904*bz_edge[15])-1.25*bz_edge[13]-0.9682458365518541*bz_edge[11]-0.5590169943749473*bz_edge[8])*jacob_edge_inv); 
  incr[14] = (-(1.75*ex_skin[14])-1.479019945774904*ex_skin[10]-1.14564392373896*ex_skin[6]-0.6614378277661477*ex_skin[4])*c*jacob_skin_inv+c2*((2.958039891549808*bz_skin[10]+1.3228756555322954*bz_skin[4])*jacob_skin_inv+(1.75*bz_edge[14]+1.479019945774904*bz_edge[10]+1.14564392373896*bz_edge[6]+0.6614378277661477*bz_edge[4])*jacob_edge_inv); 
  incr[15] = (-(1.75*ex_skin[15])-1.479019945774904*ex_skin[13]-1.14564392373896*ex_skin[11]-0.6614378277661477*ex_skin[8])*c*jacob_skin_inv+c2*((2.958039891549808*bz_skin[13]+1.3228756555322954*bz_skin[8])*jacob_skin_inv+(1.75*bz_edge[15]+1.479019945774904*bz_edge[13]+1.14564392373896*bz_edge[11]+0.6614378277661477*bz_edge[8])*jacob_edge_inv); 

  outEx[0] += incr[0]*dx1; 
  outEx[1] += incr[1]*dx1; 
  outEx[2] += incr[2]*dx1; 
  outEx[3] += incr[3]*dx1; 
  outEx[4] += incr[4]*dx1; 
  outEx[5] += incr[5]*dx1; 
  outEx[6] += incr[6]*dx1; 
  outEx[7] += incr[7]*dx1; 
  outEx[8] += incr[8]*dx1; 
  outEx[9] += incr[9]*dx1; 
  outEx[10] += incr[10]*dx1; 
  outEx[11] += incr[11]*dx1; 
  outEx[12] += incr[12]*dx1; 
  outEx[13] += incr[13]*dx1; 
  outEx[14] += incr[14]*dx1; 
  outEx[15] += incr[15]*dx1; 

  incr[0] = c2chi*((-(1.3228756555322954*ph_skin[9])-0.8660254037844386*ph_skin[2])*jacob_skin_inv+(0.6614378277661477*ph_edge[9]+0.5590169943749475*ph_edge[5]+0.4330127018922193*ph_edge[2]+0.25*ph_edge[0])*jacob_edge_inv)+(-(0.6614378277661477*ey_skin[9])-0.5590169943749475*ey_skin[5]-0.4330127018922193*ey_skin[2]-0.25*ey_skin[0])*cchi*jacob_skin_inv; 
  incr[1] = c2chi*((-(1.3228756555322954*ph_skin[12])-0.8660254037844386*ph_skin[3])*jacob_skin_inv+(0.6614378277661477*ph_edge[12]+0.5590169943749476*ph_edge[7]+0.4330127018922193*ph_edge[3]+0.25*ph_edge[1])*jacob_edge_inv)+(-(0.6614378277661477*ey_skin[12])-0.5590169943749476*ey_skin[7]-0.4330127018922193*ey_skin[3]-0.25*ey_skin[1])*cchi*jacob_skin_inv; 
  incr[2] = (-(1.14564392373896*ey_skin[9])-0.9682458365518543*ey_skin[5]-0.75*ey_skin[2]-0.4330127018922193*ey_skin[0])*cchi*jacob_skin_inv+c2chi*((-(1.9364916731037085*ph_skin[5])-0.8660254037844386*ph_skin[0])*jacob_skin_inv+(-(1.14564392373896*ph_edge[9])-0.9682458365518543*ph_edge[5]-0.75*ph_edge[2]-0.4330127018922193*ph_edge[0])*jacob_edge_inv); 
  incr[3] = (-(1.14564392373896*ey_skin[12])-0.9682458365518543*ey_skin[7]-0.75*ey_skin[3]-0.4330127018922193*ey_skin[1])*cchi*jacob_skin_inv+c2chi*((-(1.9364916731037085*ph_skin[7])-0.8660254037844386*ph_skin[1])*jacob_skin_inv+(-(1.14564392373896*ph_edge[12])-0.9682458365518543*ph_edge[7]-0.75*ph_edge[3]-0.4330127018922193*ph_edge[1])*jacob_edge_inv); 
  incr[4] = c2chi*((-(1.3228756555322954*ph_skin[14])-0.8660254037844387*ph_skin[6])*jacob_skin_inv+(0.6614378277661477*ph_edge[14]+0.5590169943749475*ph_edge[10]+0.43301270189221935*ph_edge[6]+0.25*ph_edge[4])*jacob_edge_inv)+(-(0.6614378277661477*ey_skin[14])-0.5590169943749475*ey_skin[10]-0.43301270189221935*ey_skin[6]-0.25*ey_skin[4])*cchi*jacob_skin_inv; 
  incr[5] = c2chi*((-(2.9580398915498085*ph_skin[9])-1.9364916731037085*ph_skin[2])*jacob_skin_inv+(1.4790199457749043*ph_edge[9]+1.25*ph_edge[5]+0.9682458365518543*ph_edge[2]+0.5590169943749475*ph_edge[0])*jacob_edge_inv)+(-(1.4790199457749043*ey_skin[9])-1.25*ey_skin[5]-0.9682458365518543*ey_skin[2]-0.5590169943749475*ey_skin[0])*cchi*jacob_skin_inv; 
  incr[6] = (-(1.14564392373896*ey_skin[14])-0.9682458365518543*ey_skin[10]-0.75*ey_skin[6]-0.43301270189221935*ey_skin[4])*cchi*jacob_skin_inv+c2chi*((-(1.9364916731037085*ph_skin[10])-0.8660254037844387*ph_skin[4])*jacob_skin_inv+(-(1.14564392373896*ph_edge[14])-0.9682458365518543*ph_edge[10]-0.75*ph_edge[6]-0.43301270189221935*ph_edge[4])*jacob_edge_inv); 
  incr[7] = c2chi*((-(2.9580398915498076*ph_skin[12])-1.9364916731037085*ph_skin[3])*jacob_skin_inv+(1.4790199457749038*ph_edge[12]+1.25*ph_edge[7]+0.9682458365518543*ph_edge[3]+0.5590169943749476*ph_edge[1])*jacob_edge_inv)+(-(1.4790199457749038*ey_skin[12])-1.25*ey_skin[7]-0.9682458365518543*ey_skin[3]-0.5590169943749476*ey_skin[1])*cchi*jacob_skin_inv; 
  incr[8] = c2chi*((-(1.3228756555322954*ph_skin[15])-0.8660254037844386*ph_skin[11])*jacob_skin_inv+(0.6614378277661477*ph_edge[15]+0.5590169943749473*ph_edge[13]+0.4330127018922193*ph_edge[11]+0.25*ph_edge[8])*jacob_edge_inv)+(-(0.6614378277661477*ey_skin[15])-0.5590169943749473*ey_skin[13]-0.4330127018922193*ey_skin[11]-0.25*ey_skin[8])*cchi*jacob_skin_inv; 
  incr[9] = (-(1.75*ey_skin[9])-1.4790199457749043*ey_skin[5]-1.14564392373896*ey_skin[2]-0.6614378277661477*ey_skin[0])*cchi*jacob_skin_inv+c2chi*((-(2.9580398915498085*ph_skin[5])-1.3228756555322954*ph_skin[0])*jacob_skin_inv+(-(1.75*ph_edge[9])-1.4790199457749043*ph_edge[5]-1.14564392373896*ph_edge[2]-0.6614378277661477*ph_edge[0])*jacob_edge_inv); 
  incr[10] = c2chi*((-(2.958039891549808*ph_skin[14])-1.9364916731037085*ph_skin[6])*jacob_skin_inv+(1.479019945774904*ph_edge[14]+1.25*ph_edge[10]+0.9682458365518543*ph_edge[6]+0.5590169943749475*ph_edge[4])*jacob_edge_inv)+(-(1.479019945774904*ey_skin[14])-1.25*ey_skin[10]-0.9682458365518543*ey_skin[6]-0.5590169943749475*ey_skin[4])*cchi*jacob_skin_inv; 
  incr[11] = (-(1.14564392373896*ey_skin[15])-0.9682458365518541*ey_skin[13]-0.75*ey_skin[11]-0.4330127018922193*ey_skin[8])*cchi*jacob_skin_inv+c2chi*((-(1.9364916731037083*ph_skin[13])-0.8660254037844386*ph_skin[8])*jacob_skin_inv+(-(1.14564392373896*ph_edge[15])-0.9682458365518541*ph_edge[13]-0.75*ph_edge[11]-0.4330127018922193*ph_edge[8])*jacob_edge_inv); 
  incr[12] = (-(1.75*ey_skin[12])-1.4790199457749038*ey_skin[7]-1.14564392373896*ey_skin[3]-0.6614378277661477*ey_skin[1])*cchi*jacob_skin_inv+c2chi*((-(2.9580398915498076*ph_skin[7])-1.3228756555322954*ph_skin[1])*jacob_skin_inv+(-(1.75*ph_edge[12])-1.4790199457749038*ph_edge[7]-1.14564392373896*ph_edge[3]-0.6614378277661477*ph_edge[1])*jacob_edge_inv); 
  incr[13] = c2chi*((-(2.958039891549808*ph_skin[15])-1.9364916731037083*ph_skin[11])*jacob_skin_inv+(1.479019945774904*ph_edge[15]+1.25*ph_edge[13]+0.9682458365518541*ph_edge[11]+0.5590169943749473*ph_edge[8])*jacob_edge_inv)+(-(1.479019945774904*ey_skin[15])-1.25*ey_skin[13]-0.9682458365518541*ey_skin[11]-0.5590169943749473*ey_skin[8])*cchi*jacob_skin_inv; 
  incr[14] = (-(1.75*ey_skin[14])-1.479019945774904*ey_skin[10]-1.14564392373896*ey_skin[6]-0.6614378277661477*ey_skin[4])*cchi*jacob_skin_inv+c2chi*((-(2.958039891549808*ph_skin[10])-1.3228756555322954*ph_skin[4])*jacob_skin_inv+(-(1.75*ph_edge[14])-1.479019945774904*ph_edge[10]-1.14564392373896*ph_edge[6]-0.6614378277661477*ph_edge[4])*jacob_edge_inv); 
  incr[15] = (-(1.75*ey_skin[15])-1.479019945774904*ey_skin[13]-1.14564392373896*ey_skin[11]-0.6614378277661477*ey_skin[8])*cchi*jacob_skin_inv+c2chi*((-(2.958039891549808*ph_skin[13])-1.3228756555322954*ph_skin[8])*jacob_skin_inv+(-(1.75*ph_edge[15])-1.479019945774904*ph_edge[13]-1.14564392373896*ph_edge[11]-0.6614378277661477*ph_edge[8])*jacob_edge_inv); 

  outEy[0] += incr[0]*dx1; 
  outEy[1] += incr[1]*dx1; 
  outEy[2] += incr[2]*dx1; 
  outEy[3] += incr[3]*dx1; 
  outEy[4] += incr[4]*dx1; 
  outEy[5] += incr[5]*dx1; 
  outEy[6] += incr[6]*dx1; 
  outEy[7] += incr[7]*dx1; 
  outEy[8] += incr[8]*dx1; 
  outEy[9] += incr[9]*dx1; 
  outEy[10] += incr[10]*dx1; 
  outEy[11] += incr[11]*dx1; 
  outEy[12] += incr[12]*dx1; 
  outEy[13] += incr[13]*dx1; 
  outEy[14] += incr[14]*dx1; 
  outEy[15] += incr[15]*dx1; 

  incr[0] = (-(0.6614378277661477*ez_skin[9])-0.5590169943749475*ez_skin[5]-0.4330127018922193*ez_skin[2]-0.25*ez_skin[0])*c*jacob_skin_inv+c2*((-(1.3228756555322954*bx_skin[9])-0.8660254037844386*bx_skin[2])*jacob_skin_inv+(0.6614378277661477*bx_edge[9]+0.5590169943749475*bx_edge[5]+0.4330127018922193*bx_edge[2]+0.25*bx_edge[0])*jacob_edge_inv); 
  incr[1] = (-(0.6614378277661477*ez_skin[12])-0.5590169943749476*ez_skin[7]-0.4330127018922193*ez_skin[3]-0.25*ez_skin[1])*c*jacob_skin_inv+c2*((-(1.3228756555322954*bx_skin[12])-0.8660254037844386*bx_skin[3])*jacob_skin_inv+(0.6614378277661477*bx_edge[12]+0.5590169943749476*bx_edge[7]+0.4330127018922193*bx_edge[3]+0.25*bx_edge[1])*jacob_edge_inv); 
  incr[2] = (-(1.14564392373896*ez_skin[9])-0.9682458365518543*ez_skin[5]-0.75*ez_skin[2]-0.4330127018922193*ez_skin[0])*c*jacob_skin_inv+c2*((-(1.9364916731037085*bx_skin[5])-0.8660254037844386*bx_skin[0])*jacob_skin_inv+(-(1.14564392373896*bx_edge[9])-0.9682458365518543*bx_edge[5]-0.75*bx_edge[2]-0.4330127018922193*bx_edge[0])*jacob_edge_inv); 
  incr[3] = (-(1.14564392373896*ez_skin[12])-0.9682458365518543*ez_skin[7]-0.75*ez_skin[3]-0.4330127018922193*ez_skin[1])*c*jacob_skin_inv+c2*((-(1.9364916731037085*bx_skin[7])-0.8660254037844386*bx_skin[1])*jacob_skin_inv+(-(1.14564392373896*bx_edge[12])-0.9682458365518543*bx_edge[7]-0.75*bx_edge[3]-0.4330127018922193*bx_edge[1])*jacob_edge_inv); 
  incr[4] = (-(0.6614378277661477*ez_skin[14])-0.5590169943749475*ez_skin[10]-0.43301270189221935*ez_skin[6]-0.25*ez_skin[4])*c*jacob_skin_inv+c2*((-(1.3228756555322954*bx_skin[14])-0.8660254037844387*bx_skin[6])*jacob_skin_inv+(0.6614378277661477*bx_edge[14]+0.5590169943749475*bx_edge[10]+0.43301270189221935*bx_edge[6]+0.25*bx_edge[4])*jacob_edge_inv); 
  incr[5] = (-(1.4790199457749043*ez_skin[9])-1.25*ez_skin[5]-0.9682458365518543*ez_skin[2]-0.5590169943749475*ez_skin[0])*c*jacob_skin_inv+c2*((-(2.9580398915498085*bx_skin[9])-1.9364916731037085*bx_skin[2])*jacob_skin_inv+(1.4790199457749043*bx_edge[9]+1.25*bx_edge[5]+0.9682458365518543*bx_edge[2]+0.5590169943749475*bx_edge[0])*jacob_edge_inv); 
  incr[6] = (-(1.14564392373896*ez_skin[14])-0.9682458365518543*ez_skin[10]-0.75*ez_skin[6]-0.43301270189221935*ez_skin[4])*c*jacob_skin_inv+c2*((-(1.9364916731037085*bx_skin[10])-0.8660254037844387*bx_skin[4])*jacob_skin_inv+(-(1.14564392373896*bx_edge[14])-0.9682458365518543*bx_edge[10]-0.75*bx_edge[6]-0.43301270189221935*bx_edge[4])*jacob_edge_inv); 
  incr[7] = (-(1.4790199457749038*ez_skin[12])-1.25*ez_skin[7]-0.9682458365518543*ez_skin[3]-0.5590169943749476*ez_skin[1])*c*jacob_skin_inv+c2*((-(2.9580398915498076*bx_skin[12])-1.9364916731037085*bx_skin[3])*jacob_skin_inv+(1.4790199457749038*bx_edge[12]+1.25*bx_edge[7]+0.9682458365518543*bx_edge[3]+0.5590169943749476*bx_edge[1])*jacob_edge_inv); 
  incr[8] = (-(0.6614378277661477*ez_skin[15])-0.5590169943749473*ez_skin[13]-0.4330127018922193*ez_skin[11]-0.25*ez_skin[8])*c*jacob_skin_inv+c2*((-(1.3228756555322954*bx_skin[15])-0.8660254037844386*bx_skin[11])*jacob_skin_inv+(0.6614378277661477*bx_edge[15]+0.5590169943749473*bx_edge[13]+0.4330127018922193*bx_edge[11]+0.25*bx_edge[8])*jacob_edge_inv); 
  incr[9] = (-(1.75*ez_skin[9])-1.4790199457749043*ez_skin[5]-1.14564392373896*ez_skin[2]-0.6614378277661477*ez_skin[0])*c*jacob_skin_inv+c2*((-(2.9580398915498085*bx_skin[5])-1.3228756555322954*bx_skin[0])*jacob_skin_inv+(-(1.75*bx_edge[9])-1.4790199457749043*bx_edge[5]-1.14564392373896*bx_edge[2]-0.6614378277661477*bx_edge[0])*jacob_edge_inv); 
  incr[10] = (-(1.479019945774904*ez_skin[14])-1.25*ez_skin[10]-0.9682458365518543*ez_skin[6]-0.5590169943749475*ez_skin[4])*c*jacob_skin_inv+c2*((-(2.958039891549808*bx_skin[14])-1.9364916731037085*bx_skin[6])*jacob_skin_inv+(1.479019945774904*bx_edge[14]+1.25*bx_edge[10]+0.9682458365518543*bx_edge[6]+0.5590169943749475*bx_edge[4])*jacob_edge_inv); 
  incr[11] = (-(1.14564392373896*ez_skin[15])-0.9682458365518541*ez_skin[13]-0.75*ez_skin[11]-0.4330127018922193*ez_skin[8])*c*jacob_skin_inv+c2*((-(1.9364916731037083*bx_skin[13])-0.8660254037844386*bx_skin[8])*jacob_skin_inv+(-(1.14564392373896*bx_edge[15])-0.9682458365518541*bx_edge[13]-0.75*bx_edge[11]-0.4330127018922193*bx_edge[8])*jacob_edge_inv); 
  incr[12] = (-(1.75*ez_skin[12])-1.4790199457749038*ez_skin[7]-1.14564392373896*ez_skin[3]-0.6614378277661477*ez_skin[1])*c*jacob_skin_inv+c2*((-(2.9580398915498076*bx_skin[7])-1.3228756555322954*bx_skin[1])*jacob_skin_inv+(-(1.75*bx_edge[12])-1.4790199457749038*bx_edge[7]-1.14564392373896*bx_edge[3]-0.6614378277661477*bx_edge[1])*jacob_edge_inv); 
  incr[13] = (-(1.479019945774904*ez_skin[15])-1.25*ez_skin[13]-0.9682458365518541*ez_skin[11]-0.5590169943749473*ez_skin[8])*c*jacob_skin_inv+c2*((-(2.958039891549808*bx_skin[15])-1.9364916731037083*bx_skin[11])*jacob_skin_inv+(1.479019945774904*bx_edge[15]+1.25*bx_edge[13]+0.9682458365518541*bx_edge[11]+0.5590169943749473*bx_edge[8])*jacob_edge_inv); 
  incr[14] = (-(1.75*ez_skin[14])-1.479019945774904*ez_skin[10]-1.14564392373896*ez_skin[6]-0.6614378277661477*ez_skin[4])*c*jacob_skin_inv+c2*((-(2.958039891549808*bx_skin[10])-1.3228756555322954*bx_skin[4])*jacob_skin_inv+(-(1.75*bx_edge[14])-1.479019945774904*bx_edge[10]-1.14564392373896*bx_edge[6]-0.6614378277661477*bx_edge[4])*jacob_edge_inv); 
  incr[15] = (-(1.75*ez_skin[15])-1.479019945774904*ez_skin[13]-1.14564392373896*ez_skin[11]-0.6614378277661477*ez_skin[8])*c*jacob_skin_inv+c2*((-(2.958039891549808*bx_skin[13])-1.3228756555322954*bx_skin[8])*jacob_skin_inv+(-(1.75*bx_edge[15])-1.479019945774904*bx_edge[13]-1.14564392373896*bx_edge[11]-0.6614378277661477*bx_edge[8])*jacob_edge_inv); 

  outEz[0] += incr[0]*dx1; 
  outEz[1] += incr[1]*dx1; 
  outEz[2] += incr[2]*dx1; 
  outEz[3] += incr[3]*dx1; 
  outEz[4] += incr[4]*dx1; 
  outEz[5] += incr[5]*dx1; 
  outEz[6] += incr[6]*dx1; 
  outEz[7] += incr[7]*dx1; 
  outEz[8] += incr[8]*dx1; 
  outEz[9] += incr[9]*dx1; 
  outEz[10] += incr[10]*dx1; 
  outEz[11] += incr[11]*dx1; 
  outEz[12] += incr[12]*dx1; 
  outEz[13] += incr[13]*dx1; 
  outEz[14] += incr[14]*dx1; 
  outEz[15] += incr[15]*dx1; 

  incr[0] = ((-(0.6614378277661477*bx_skin[9])-0.5590169943749475*bx_skin[5]-0.4330127018922193*bx_skin[2]-0.25*bx_skin[0])*c-1.3228756555322954*ez_skin[9]-0.8660254037844386*ez_skin[2])*jacob_skin_inv+(0.6614378277661477*ez_edge[9]+0.5590169943749475*ez_edge[5]+0.4330127018922193*ez_edge[2]+0.25*ez_edge[0])*jacob_edge_inv; 
  incr[1] = ((-(0.6614378277661477*bx_skin[12])-0.5590169943749476*bx_skin[7]-0.4330127018922193*bx_skin[3]-0.25*bx_skin[1])*c-1.3228756555322954*ez_skin[12]-0.8660254037844386*ez_skin[3])*jacob_skin_inv+(0.6614378277661477*ez_edge[12]+0.5590169943749476*ez_edge[7]+0.4330127018922193*ez_edge[3]+0.25*ez_edge[1])*jacob_edge_inv; 
  incr[2] = ((-(1.14564392373896*bx_skin[9])-0.9682458365518543*bx_skin[5]-0.75*bx_skin[2]-0.4330127018922193*bx_skin[0])*c-1.9364916731037085*ez_skin[5]-0.8660254037844386*ez_skin[0])*jacob_skin_inv+(-(1.14564392373896*ez_edge[9])-0.9682458365518543*ez_edge[5]-0.75*ez_edge[2]-0.4330127018922193*ez_edge[0])*jacob_edge_inv; 
  incr[3] = ((-(1.14564392373896*bx_skin[12])-0.9682458365518543*bx_skin[7]-0.75*bx_skin[3]-0.4330127018922193*bx_skin[1])*c-1.9364916731037085*ez_skin[7]-0.8660254037844386*ez_skin[1])*jacob_skin_inv+(-(1.14564392373896*ez_edge[12])-0.9682458365518543*ez_edge[7]-0.75*ez_edge[3]-0.4330127018922193*ez_edge[1])*jacob_edge_inv; 
  incr[4] = ((-(0.6614378277661477*bx_skin[14])-0.5590169943749475*bx_skin[10]-0.43301270189221935*bx_skin[6]-0.25*bx_skin[4])*c-1.3228756555322954*ez_skin[14]-0.8660254037844387*ez_skin[6])*jacob_skin_inv+(0.6614378277661477*ez_edge[14]+0.5590169943749475*ez_edge[10]+0.43301270189221935*ez_edge[6]+0.25*ez_edge[4])*jacob_edge_inv; 
  incr[5] = ((-(1.4790199457749043*bx_skin[9])-1.25*bx_skin[5]-0.9682458365518543*bx_skin[2]-0.5590169943749475*bx_skin[0])*c-2.9580398915498085*ez_skin[9]-1.9364916731037085*ez_skin[2])*jacob_skin_inv+(1.4790199457749043*ez_edge[9]+1.25*ez_edge[5]+0.9682458365518543*ez_edge[2]+0.5590169943749475*ez_edge[0])*jacob_edge_inv; 
  incr[6] = ((-(1.14564392373896*bx_skin[14])-0.9682458365518543*bx_skin[10]-0.75*bx_skin[6]-0.43301270189221935*bx_skin[4])*c-1.9364916731037085*ez_skin[10]-0.8660254037844387*ez_skin[4])*jacob_skin_inv+(-(1.14564392373896*ez_edge[14])-0.9682458365518543*ez_edge[10]-0.75*ez_edge[6]-0.43301270189221935*ez_edge[4])*jacob_edge_inv; 
  incr[7] = ((-(1.4790199457749038*bx_skin[12])-1.25*bx_skin[7]-0.9682458365518543*bx_skin[3]-0.5590169943749476*bx_skin[1])*c-2.9580398915498076*ez_skin[12]-1.9364916731037085*ez_skin[3])*jacob_skin_inv+(1.4790199457749038*ez_edge[12]+1.25*ez_edge[7]+0.9682458365518543*ez_edge[3]+0.5590169943749476*ez_edge[1])*jacob_edge_inv; 
  incr[8] = ((-(0.6614378277661477*bx_skin[15])-0.5590169943749473*bx_skin[13]-0.4330127018922193*bx_skin[11]-0.25*bx_skin[8])*c-1.3228756555322954*ez_skin[15]-0.8660254037844386*ez_skin[11])*jacob_skin_inv+(0.6614378277661477*ez_edge[15]+0.5590169943749473*ez_edge[13]+0.4330127018922193*ez_edge[11]+0.25*ez_edge[8])*jacob_edge_inv; 
  incr[9] = ((-(1.75*bx_skin[9])-1.4790199457749043*bx_skin[5]-1.14564392373896*bx_skin[2]-0.6614378277661477*bx_skin[0])*c-2.9580398915498085*ez_skin[5]-1.3228756555322954*ez_skin[0])*jacob_skin_inv+(-(1.75*ez_edge[9])-1.4790199457749043*ez_edge[5]-1.14564392373896*ez_edge[2]-0.6614378277661477*ez_edge[0])*jacob_edge_inv; 
  incr[10] = ((-(1.479019945774904*bx_skin[14])-1.25*bx_skin[10]-0.9682458365518543*bx_skin[6]-0.5590169943749475*bx_skin[4])*c-2.958039891549808*ez_skin[14]-1.9364916731037085*ez_skin[6])*jacob_skin_inv+(1.479019945774904*ez_edge[14]+1.25*ez_edge[10]+0.9682458365518543*ez_edge[6]+0.5590169943749475*ez_edge[4])*jacob_edge_inv; 
  incr[11] = ((-(1.14564392373896*bx_skin[15])-0.9682458365518541*bx_skin[13]-0.75*bx_skin[11]-0.4330127018922193*bx_skin[8])*c-1.9364916731037083*ez_skin[13]-0.8660254037844386*ez_skin[8])*jacob_skin_inv+(-(1.14564392373896*ez_edge[15])-0.9682458365518541*ez_edge[13]-0.75*ez_edge[11]-0.4330127018922193*ez_edge[8])*jacob_edge_inv; 
  incr[12] = ((-(1.75*bx_skin[12])-1.4790199457749038*bx_skin[7]-1.14564392373896*bx_skin[3]-0.6614378277661477*bx_skin[1])*c-2.9580398915498076*ez_skin[7]-1.3228756555322954*ez_skin[1])*jacob_skin_inv+(-(1.75*ez_edge[12])-1.4790199457749038*ez_edge[7]-1.14564392373896*ez_edge[3]-0.6614378277661477*ez_edge[1])*jacob_edge_inv; 
  incr[13] = ((-(1.479019945774904*bx_skin[15])-1.25*bx_skin[13]-0.9682458365518541*bx_skin[11]-0.5590169943749473*bx_skin[8])*c-2.958039891549808*ez_skin[15]-1.9364916731037083*ez_skin[11])*jacob_skin_inv+(1.479019945774904*ez_edge[15]+1.25*ez_edge[13]+0.9682458365518541*ez_edge[11]+0.5590169943749473*ez_edge[8])*jacob_edge_inv; 
  incr[14] = ((-(1.75*bx_skin[14])-1.479019945774904*bx_skin[10]-1.14564392373896*bx_skin[6]-0.6614378277661477*bx_skin[4])*c-2.958039891549808*ez_skin[10]-1.3228756555322954*ez_skin[4])*jacob_skin_inv+(-(1.75*ez_edge[14])-1.479019945774904*ez_edge[10]-1.14564392373896*ez_edge[6]-0.6614378277661477*ez_edge[4])*jacob_edge_inv; 
  incr[15] = ((-(1.75*bx_skin[15])-1.479019945774904*bx_skin[13]-1.14564392373896*bx_skin[11]-0.6614378277661477*bx_skin[8])*c-2.958039891549808*ez_skin[13]-1.3228756555322954*ez_skin[8])*jacob_skin_inv+(-(1.75*ez_edge[15])-1.479019945774904*ez_edge[13]-1.14564392373896*ez_edge[11]-0.6614378277661477*ez_edge[8])*jacob_edge_inv; 

  outBx[0] += incr[0]*dx1; 
  outBx[1] += incr[1]*dx1; 
  outBx[2] += incr[2]*dx1; 
  outBx[3] += incr[3]*dx1; 
  outBx[4] += incr[4]*dx1; 
  outBx[5] += incr[5]*dx1; 
  outBx[6] += incr[6]*dx1; 
  outBx[7] += incr[7]*dx1; 
  outBx[8] += incr[8]*dx1; 
  outBx[9] += incr[9]*dx1; 
  outBx[10] += incr[10]*dx1; 
  outBx[11] += incr[11]*dx1; 
  outBx[12] += incr[12]*dx1; 
  outBx[13] += incr[13]*dx1; 
  outBx[14] += incr[14]*dx1; 
  outBx[15] += incr[15]*dx1; 

  incr[0] = ((-(1.3228756555322954*ps_skin[9])-0.8660254037844386*ps_skin[2])*jacob_skin_inv+(0.6614378277661477*ps_edge[9]+0.5590169943749475*ps_edge[5]+0.4330127018922193*ps_edge[2]+0.25*ps_edge[0])*jacob_edge_inv)*gamma+(-(0.6614378277661477*by_skin[9])-0.5590169943749475*by_skin[5]-0.4330127018922193*by_skin[2]-0.25*by_skin[0])*cgamma*jacob_skin_inv; 
  incr[1] = ((-(1.3228756555322954*ps_skin[12])-0.8660254037844386*ps_skin[3])*jacob_skin_inv+(0.6614378277661477*ps_edge[12]+0.5590169943749476*ps_edge[7]+0.4330127018922193*ps_edge[3]+0.25*ps_edge[1])*jacob_edge_inv)*gamma+(-(0.6614378277661477*by_skin[12])-0.5590169943749476*by_skin[7]-0.4330127018922193*by_skin[3]-0.25*by_skin[1])*cgamma*jacob_skin_inv; 
  incr[2] = ((-(1.9364916731037085*ps_skin[5])-0.8660254037844386*ps_skin[0])*jacob_skin_inv+(-(1.14564392373896*ps_edge[9])-0.9682458365518543*ps_edge[5]-0.75*ps_edge[2]-0.4330127018922193*ps_edge[0])*jacob_edge_inv)*gamma+(-(1.14564392373896*by_skin[9])-0.9682458365518543*by_skin[5]-0.75*by_skin[2]-0.4330127018922193*by_skin[0])*cgamma*jacob_skin_inv; 
  incr[3] = ((-(1.9364916731037085*ps_skin[7])-0.8660254037844386*ps_skin[1])*jacob_skin_inv+(-(1.14564392373896*ps_edge[12])-0.9682458365518543*ps_edge[7]-0.75*ps_edge[3]-0.4330127018922193*ps_edge[1])*jacob_edge_inv)*gamma+(-(1.14564392373896*by_skin[12])-0.9682458365518543*by_skin[7]-0.75*by_skin[3]-0.4330127018922193*by_skin[1])*cgamma*jacob_skin_inv; 
  incr[4] = ((-(1.3228756555322954*ps_skin[14])-0.8660254037844387*ps_skin[6])*jacob_skin_inv+(0.6614378277661477*ps_edge[14]+0.5590169943749475*ps_edge[10]+0.43301270189221935*ps_edge[6]+0.25*ps_edge[4])*jacob_edge_inv)*gamma+(-(0.6614378277661477*by_skin[14])-0.5590169943749475*by_skin[10]-0.43301270189221935*by_skin[6]-0.25*by_skin[4])*cgamma*jacob_skin_inv; 
  incr[5] = ((-(2.9580398915498085*ps_skin[9])-1.9364916731037085*ps_skin[2])*jacob_skin_inv+(1.4790199457749043*ps_edge[9]+1.25*ps_edge[5]+0.9682458365518543*ps_edge[2]+0.5590169943749475*ps_edge[0])*jacob_edge_inv)*gamma+(-(1.4790199457749043*by_skin[9])-1.25*by_skin[5]-0.9682458365518543*by_skin[2]-0.5590169943749475*by_skin[0])*cgamma*jacob_skin_inv; 
  incr[6] = ((-(1.9364916731037085*ps_skin[10])-0.8660254037844387*ps_skin[4])*jacob_skin_inv+(-(1.14564392373896*ps_edge[14])-0.9682458365518543*ps_edge[10]-0.75*ps_edge[6]-0.43301270189221935*ps_edge[4])*jacob_edge_inv)*gamma+(-(1.14564392373896*by_skin[14])-0.9682458365518543*by_skin[10]-0.75*by_skin[6]-0.43301270189221935*by_skin[4])*cgamma*jacob_skin_inv; 
  incr[7] = ((-(2.9580398915498076*ps_skin[12])-1.9364916731037085*ps_skin[3])*jacob_skin_inv+(1.4790199457749038*ps_edge[12]+1.25*ps_edge[7]+0.9682458365518543*ps_edge[3]+0.5590169943749476*ps_edge[1])*jacob_edge_inv)*gamma+(-(1.4790199457749038*by_skin[12])-1.25*by_skin[7]-0.9682458365518543*by_skin[3]-0.5590169943749476*by_skin[1])*cgamma*jacob_skin_inv; 
  incr[8] = ((-(1.3228756555322954*ps_skin[15])-0.8660254037844386*ps_skin[11])*jacob_skin_inv+(0.6614378277661477*ps_edge[15]+0.5590169943749473*ps_edge[13]+0.4330127018922193*ps_edge[11]+0.25*ps_edge[8])*jacob_edge_inv)*gamma+(-(0.6614378277661477*by_skin[15])-0.5590169943749473*by_skin[13]-0.4330127018922193*by_skin[11]-0.25*by_skin[8])*cgamma*jacob_skin_inv; 
  incr[9] = ((-(2.9580398915498085*ps_skin[5])-1.3228756555322954*ps_skin[0])*jacob_skin_inv+(-(1.75*ps_edge[9])-1.4790199457749043*ps_edge[5]-1.14564392373896*ps_edge[2]-0.6614378277661477*ps_edge[0])*jacob_edge_inv)*gamma+(-(1.75*by_skin[9])-1.4790199457749043*by_skin[5]-1.14564392373896*by_skin[2]-0.6614378277661477*by_skin[0])*cgamma*jacob_skin_inv; 
  incr[10] = ((-(2.958039891549808*ps_skin[14])-1.9364916731037085*ps_skin[6])*jacob_skin_inv+(1.479019945774904*ps_edge[14]+1.25*ps_edge[10]+0.9682458365518543*ps_edge[6]+0.5590169943749475*ps_edge[4])*jacob_edge_inv)*gamma+(-(1.479019945774904*by_skin[14])-1.25*by_skin[10]-0.9682458365518543*by_skin[6]-0.5590169943749475*by_skin[4])*cgamma*jacob_skin_inv; 
  incr[11] = ((-(1.9364916731037083*ps_skin[13])-0.8660254037844386*ps_skin[8])*jacob_skin_inv+(-(1.14564392373896*ps_edge[15])-0.9682458365518541*ps_edge[13]-0.75*ps_edge[11]-0.4330127018922193*ps_edge[8])*jacob_edge_inv)*gamma+(-(1.14564392373896*by_skin[15])-0.9682458365518541*by_skin[13]-0.75*by_skin[11]-0.4330127018922193*by_skin[8])*cgamma*jacob_skin_inv; 
  incr[12] = ((-(2.9580398915498076*ps_skin[7])-1.3228756555322954*ps_skin[1])*jacob_skin_inv+(-(1.75*ps_edge[12])-1.4790199457749038*ps_edge[7]-1.14564392373896*ps_edge[3]-0.6614378277661477*ps_edge[1])*jacob_edge_inv)*gamma+(-(1.75*by_skin[12])-1.4790199457749038*by_skin[7]-1.14564392373896*by_skin[3]-0.6614378277661477*by_skin[1])*cgamma*jacob_skin_inv; 
  incr[13] = ((-(2.958039891549808*ps_skin[15])-1.9364916731037083*ps_skin[11])*jacob_skin_inv+(1.479019945774904*ps_edge[15]+1.25*ps_edge[13]+0.9682458365518541*ps_edge[11]+0.5590169943749473*ps_edge[8])*jacob_edge_inv)*gamma+(-(1.479019945774904*by_skin[15])-1.25*by_skin[13]-0.9682458365518541*by_skin[11]-0.5590169943749473*by_skin[8])*cgamma*jacob_skin_inv; 
  incr[14] = ((-(2.958039891549808*ps_skin[10])-1.3228756555322954*ps_skin[4])*jacob_skin_inv+(-(1.75*ps_edge[14])-1.479019945774904*ps_edge[10]-1.14564392373896*ps_edge[6]-0.6614378277661477*ps_edge[4])*jacob_edge_inv)*gamma+(-(1.75*by_skin[14])-1.479019945774904*by_skin[10]-1.14564392373896*by_skin[6]-0.6614378277661477*by_skin[4])*cgamma*jacob_skin_inv; 
  incr[15] = ((-(2.958039891549808*ps_skin[13])-1.3228756555322954*ps_skin[8])*jacob_skin_inv+(-(1.75*ps_edge[15])-1.479019945774904*ps_edge[13]-1.14564392373896*ps_edge[11]-0.6614378277661477*ps_edge[8])*jacob_edge_inv)*gamma+(-(1.75*by_skin[15])-1.479019945774904*by_skin[13]-1.14564392373896*by_skin[11]-0.6614378277661477*by_skin[8])*cgamma*jacob_skin_inv; 

  outBy[0] += incr[0]*dx1; 
  outBy[1] += incr[1]*dx1; 
  outBy[2] += incr[2]*dx1; 
  outBy[3] += incr[3]*dx1; 
  outBy[4] += incr[4]*dx1; 
  outBy[5] += incr[5]*dx1; 
  outBy[6] += incr[6]*dx1; 
  outBy[7] += incr[7]*dx1; 
  outBy[8] += incr[8]*dx1; 
  outBy[9] += incr[9]*dx1; 
  outBy[10] += incr[10]*dx1; 
  outBy[11] += incr[11]*dx1; 
  outBy[12] += incr[12]*dx1; 
  outBy[13] += incr[13]*dx1; 
  outBy[14] += incr[14]*dx1; 
  outBy[15] += incr[15]*dx1; 

  incr[0] = ((-(0.6614378277661477*bz_skin[9])-0.5590169943749475*bz_skin[5]-0.4330127018922193*bz_skin[2]-0.25*bz_skin[0])*c+1.3228756555322954*ex_skin[9]+0.8660254037844386*ex_skin[2])*jacob_skin_inv+(-(0.6614378277661477*ex_edge[9])-0.5590169943749475*ex_edge[5]-0.4330127018922193*ex_edge[2]-0.25*ex_edge[0])*jacob_edge_inv; 
  incr[1] = ((-(0.6614378277661477*bz_skin[12])-0.5590169943749476*bz_skin[7]-0.4330127018922193*bz_skin[3]-0.25*bz_skin[1])*c+1.3228756555322954*ex_skin[12]+0.8660254037844386*ex_skin[3])*jacob_skin_inv+(-(0.6614378277661477*ex_edge[12])-0.5590169943749476*ex_edge[7]-0.4330127018922193*ex_edge[3]-0.25*ex_edge[1])*jacob_edge_inv; 
  incr[2] = ((-(1.14564392373896*bz_skin[9])-0.9682458365518543*bz_skin[5]-0.75*bz_skin[2]-0.4330127018922193*bz_skin[0])*c+1.9364916731037085*ex_skin[5]+0.8660254037844386*ex_skin[0])*jacob_skin_inv+(1.14564392373896*ex_edge[9]+0.9682458365518543*ex_edge[5]+0.75*ex_edge[2]+0.4330127018922193*ex_edge[0])*jacob_edge_inv; 
  incr[3] = ((-(1.14564392373896*bz_skin[12])-0.9682458365518543*bz_skin[7]-0.75*bz_skin[3]-0.4330127018922193*bz_skin[1])*c+1.9364916731037085*ex_skin[7]+0.8660254037844386*ex_skin[1])*jacob_skin_inv+(1.14564392373896*ex_edge[12]+0.9682458365518543*ex_edge[7]+0.75*ex_edge[3]+0.4330127018922193*ex_edge[1])*jacob_edge_inv; 
  incr[4] = ((-(0.6614378277661477*bz_skin[14])-0.5590169943749475*bz_skin[10]-0.43301270189221935*bz_skin[6]-0.25*bz_skin[4])*c+1.3228756555322954*ex_skin[14]+0.8660254037844387*ex_skin[6])*jacob_skin_inv+(-(0.6614378277661477*ex_edge[14])-0.5590169943749475*ex_edge[10]-0.43301270189221935*ex_edge[6]-0.25*ex_edge[4])*jacob_edge_inv; 
  incr[5] = ((-(1.4790199457749043*bz_skin[9])-1.25*bz_skin[5]-0.9682458365518543*bz_skin[2]-0.5590169943749475*bz_skin[0])*c+2.9580398915498085*ex_skin[9]+1.9364916731037085*ex_skin[2])*jacob_skin_inv+(-(1.4790199457749043*ex_edge[9])-1.25*ex_edge[5]-0.9682458365518543*ex_edge[2]-0.5590169943749475*ex_edge[0])*jacob_edge_inv; 
  incr[6] = ((-(1.14564392373896*bz_skin[14])-0.9682458365518543*bz_skin[10]-0.75*bz_skin[6]-0.43301270189221935*bz_skin[4])*c+1.9364916731037085*ex_skin[10]+0.8660254037844387*ex_skin[4])*jacob_skin_inv+(1.14564392373896*ex_edge[14]+0.9682458365518543*ex_edge[10]+0.75*ex_edge[6]+0.43301270189221935*ex_edge[4])*jacob_edge_inv; 
  incr[7] = ((-(1.4790199457749038*bz_skin[12])-1.25*bz_skin[7]-0.9682458365518543*bz_skin[3]-0.5590169943749476*bz_skin[1])*c+2.9580398915498076*ex_skin[12]+1.9364916731037085*ex_skin[3])*jacob_skin_inv+(-(1.4790199457749038*ex_edge[12])-1.25*ex_edge[7]-0.9682458365518543*ex_edge[3]-0.5590169943749476*ex_edge[1])*jacob_edge_inv; 
  incr[8] = ((-(0.6614378277661477*bz_skin[15])-0.5590169943749473*bz_skin[13]-0.4330127018922193*bz_skin[11]-0.25*bz_skin[8])*c+1.3228756555322954*ex_skin[15]+0.8660254037844386*ex_skin[11])*jacob_skin_inv+(-(0.6614378277661477*ex_edge[15])-0.5590169943749473*ex_edge[13]-0.4330127018922193*ex_edge[11]-0.25*ex_edge[8])*jacob_edge_inv; 
  incr[9] = ((-(1.75*bz_skin[9])-1.4790199457749043*bz_skin[5]-1.14564392373896*bz_skin[2]-0.6614378277661477*bz_skin[0])*c+2.9580398915498085*ex_skin[5]+1.3228756555322954*ex_skin[0])*jacob_skin_inv+(1.75*ex_edge[9]+1.4790199457749043*ex_edge[5]+1.14564392373896*ex_edge[2]+0.6614378277661477*ex_edge[0])*jacob_edge_inv; 
  incr[10] = ((-(1.479019945774904*bz_skin[14])-1.25*bz_skin[10]-0.9682458365518543*bz_skin[6]-0.5590169943749475*bz_skin[4])*c+2.958039891549808*ex_skin[14]+1.9364916731037085*ex_skin[6])*jacob_skin_inv+(-(1.479019945774904*ex_edge[14])-1.25*ex_edge[10]-0.9682458365518543*ex_edge[6]-0.5590169943749475*ex_edge[4])*jacob_edge_inv; 
  incr[11] = ((-(1.14564392373896*bz_skin[15])-0.9682458365518541*bz_skin[13]-0.75*bz_skin[11]-0.4330127018922193*bz_skin[8])*c+1.9364916731037083*ex_skin[13]+0.8660254037844386*ex_skin[8])*jacob_skin_inv+(1.14564392373896*ex_edge[15]+0.9682458365518541*ex_edge[13]+0.75*ex_edge[11]+0.4330127018922193*ex_edge[8])*jacob_edge_inv; 
  incr[12] = ((-(1.75*bz_skin[12])-1.4790199457749038*bz_skin[7]-1.14564392373896*bz_skin[3]-0.6614378277661477*bz_skin[1])*c+2.9580398915498076*ex_skin[7]+1.3228756555322954*ex_skin[1])*jacob_skin_inv+(1.75*ex_edge[12]+1.4790199457749038*ex_edge[7]+1.14564392373896*ex_edge[3]+0.6614378277661477*ex_edge[1])*jacob_edge_inv; 
  incr[13] = ((-(1.479019945774904*bz_skin[15])-1.25*bz_skin[13]-0.9682458365518541*bz_skin[11]-0.5590169943749473*bz_skin[8])*c+2.958039891549808*ex_skin[15]+1.9364916731037083*ex_skin[11])*jacob_skin_inv+(-(1.479019945774904*ex_edge[15])-1.25*ex_edge[13]-0.9682458365518541*ex_edge[11]-0.5590169943749473*ex_edge[8])*jacob_edge_inv; 
  incr[14] = ((-(1.75*bz_skin[14])-1.479019945774904*bz_skin[10]-1.14564392373896*bz_skin[6]-0.6614378277661477*bz_skin[4])*c+2.958039891549808*ex_skin[10]+1.3228756555322954*ex_skin[4])*jacob_skin_inv+(1.75*ex_edge[14]+1.479019945774904*ex_edge[10]+1.14564392373896*ex_edge[6]+0.6614378277661477*ex_edge[4])*jacob_edge_inv; 
  incr[15] = ((-(1.75*bz_skin[15])-1.479019945774904*bz_skin[13]-1.14564392373896*bz_skin[11]-0.6614378277661477*bz_skin[8])*c+2.958039891549808*ex_skin[13]+1.3228756555322954*ex_skin[8])*jacob_skin_inv+(1.75*ex_edge[15]+1.479019945774904*ex_edge[13]+1.14564392373896*ex_edge[11]+0.6614378277661477*ex_edge[8])*jacob_edge_inv; 

  outBz[0] += incr[0]*dx1; 
  outBz[1] += incr[1]*dx1; 
  outBz[2] += incr[2]*dx1; 
  outBz[3] += incr[3]*dx1; 
  outBz[4] += incr[4]*dx1; 
  outBz[5] += incr[5]*dx1; 
  outBz[6] += incr[6]*dx1; 
  outBz[7] += incr[7]*dx1; 
  outBz[8] += incr[8]*dx1; 
  outBz[9] += incr[9]*dx1; 
  outBz[10] += incr[10]*dx1; 
  outBz[11] += incr[11]*dx1; 
  outBz[12] += incr[12]*dx1; 
  outBz[13] += incr[13]*dx1; 
  outBz[14] += incr[14]*dx1; 
  outBz[15] += incr[15]*dx1; 

  incr[0] = (-(0.6614378277661477*ph_skin[9])-0.5590169943749475*ph_skin[5]-0.4330127018922193*ph_skin[2]-0.25*ph_skin[0])*cchi*jacob_skin_inv+chi*((-(1.3228756555322954*ey_skin[9])-0.8660254037844386*ey_skin[2])*jacob_skin_inv+(0.6614378277661477*ey_edge[9]+0.5590169943749475*ey_edge[5]+0.4330127018922193*ey_edge[2]+0.25*ey_edge[0])*jacob_edge_inv); 
  incr[1] = (-(0.6614378277661477*ph_skin[12])-0.5590169943749476*ph_skin[7]-0.4330127018922193*ph_skin[3]-0.25*ph_skin[1])*cchi*jacob_skin_inv+chi*((-(1.3228756555322954*ey_skin[12])-0.8660254037844386*ey_skin[3])*jacob_skin_inv+(0.6614378277661477*ey_edge[12]+0.5590169943749476*ey_edge[7]+0.4330127018922193*ey_edge[3]+0.25*ey_edge[1])*jacob_edge_inv); 
  incr[2] = (-(1.14564392373896*ph_skin[9])-0.9682458365518543*ph_skin[5]-0.75*ph_skin[2]-0.4330127018922193*ph_skin[0])*cchi*jacob_skin_inv+chi*((-(1.9364916731037085*ey_skin[5])-0.8660254037844386*ey_skin[0])*jacob_skin_inv+(-(1.14564392373896*ey_edge[9])-0.9682458365518543*ey_edge[5]-0.75*ey_edge[2]-0.4330127018922193*ey_edge[0])*jacob_edge_inv); 
  incr[3] = (-(1.14564392373896*ph_skin[12])-0.9682458365518543*ph_skin[7]-0.75*ph_skin[3]-0.4330127018922193*ph_skin[1])*cchi*jacob_skin_inv+chi*((-(1.9364916731037085*ey_skin[7])-0.8660254037844386*ey_skin[1])*jacob_skin_inv+(-(1.14564392373896*ey_edge[12])-0.9682458365518543*ey_edge[7]-0.75*ey_edge[3]-0.4330127018922193*ey_edge[1])*jacob_edge_inv); 
  incr[4] = (-(0.6614378277661477*ph_skin[14])-0.5590169943749475*ph_skin[10]-0.43301270189221935*ph_skin[6]-0.25*ph_skin[4])*cchi*jacob_skin_inv+chi*((-(1.3228756555322954*ey_skin[14])-0.8660254037844387*ey_skin[6])*jacob_skin_inv+(0.6614378277661477*ey_edge[14]+0.5590169943749475*ey_edge[10]+0.43301270189221935*ey_edge[6]+0.25*ey_edge[4])*jacob_edge_inv); 
  incr[5] = (-(1.4790199457749043*ph_skin[9])-1.25*ph_skin[5]-0.9682458365518543*ph_skin[2]-0.5590169943749475*ph_skin[0])*cchi*jacob_skin_inv+chi*((-(2.9580398915498085*ey_skin[9])-1.9364916731037085*ey_skin[2])*jacob_skin_inv+(1.4790199457749043*ey_edge[9]+1.25*ey_edge[5]+0.9682458365518543*ey_edge[2]+0.5590169943749475*ey_edge[0])*jacob_edge_inv); 
  incr[6] = (-(1.14564392373896*ph_skin[14])-0.9682458365518543*ph_skin[10]-0.75*ph_skin[6]-0.43301270189221935*ph_skin[4])*cchi*jacob_skin_inv+chi*((-(1.9364916731037085*ey_skin[10])-0.8660254037844387*ey_skin[4])*jacob_skin_inv+(-(1.14564392373896*ey_edge[14])-0.9682458365518543*ey_edge[10]-0.75*ey_edge[6]-0.43301270189221935*ey_edge[4])*jacob_edge_inv); 
  incr[7] = (-(1.4790199457749038*ph_skin[12])-1.25*ph_skin[7]-0.9682458365518543*ph_skin[3]-0.5590169943749476*ph_skin[1])*cchi*jacob_skin_inv+chi*((-(2.9580398915498076*ey_skin[12])-1.9364916731037085*ey_skin[3])*jacob_skin_inv+(1.4790199457749038*ey_edge[12]+1.25*ey_edge[7]+0.9682458365518543*ey_edge[3]+0.5590169943749476*ey_edge[1])*jacob_edge_inv); 
  incr[8] = (-(0.6614378277661477*ph_skin[15])-0.5590169943749473*ph_skin[13]-0.4330127018922193*ph_skin[11]-0.25*ph_skin[8])*cchi*jacob_skin_inv+chi*((-(1.3228756555322954*ey_skin[15])-0.8660254037844386*ey_skin[11])*jacob_skin_inv+(0.6614378277661477*ey_edge[15]+0.5590169943749473*ey_edge[13]+0.4330127018922193*ey_edge[11]+0.25*ey_edge[8])*jacob_edge_inv); 
  incr[9] = (-(1.75*ph_skin[9])-1.4790199457749043*ph_skin[5]-1.14564392373896*ph_skin[2]-0.6614378277661477*ph_skin[0])*cchi*jacob_skin_inv+chi*((-(2.9580398915498085*ey_skin[5])-1.3228756555322954*ey_skin[0])*jacob_skin_inv+(-(1.75*ey_edge[9])-1.4790199457749043*ey_edge[5]-1.14564392373896*ey_edge[2]-0.6614378277661477*ey_edge[0])*jacob_edge_inv); 
  incr[10] = (-(1.479019945774904*ph_skin[14])-1.25*ph_skin[10]-0.9682458365518543*ph_skin[6]-0.5590169943749475*ph_skin[4])*cchi*jacob_skin_inv+chi*((-(2.958039891549808*ey_skin[14])-1.9364916731037085*ey_skin[6])*jacob_skin_inv+(1.479019945774904*ey_edge[14]+1.25*ey_edge[10]+0.9682458365518543*ey_edge[6]+0.5590169943749475*ey_edge[4])*jacob_edge_inv); 
  incr[11] = (-(1.14564392373896*ph_skin[15])-0.9682458365518541*ph_skin[13]-0.75*ph_skin[11]-0.4330127018922193*ph_skin[8])*cchi*jacob_skin_inv+chi*((-(1.9364916731037083*ey_skin[13])-0.8660254037844386*ey_skin[8])*jacob_skin_inv+(-(1.14564392373896*ey_edge[15])-0.9682458365518541*ey_edge[13]-0.75*ey_edge[11]-0.4330127018922193*ey_edge[8])*jacob_edge_inv); 
  incr[12] = (-(1.75*ph_skin[12])-1.4790199457749038*ph_skin[7]-1.14564392373896*ph_skin[3]-0.6614378277661477*ph_skin[1])*cchi*jacob_skin_inv+chi*((-(2.9580398915498076*ey_skin[7])-1.3228756555322954*ey_skin[1])*jacob_skin_inv+(-(1.75*ey_edge[12])-1.4790199457749038*ey_edge[7]-1.14564392373896*ey_edge[3]-0.6614378277661477*ey_edge[1])*jacob_edge_inv); 
  incr[13] = (-(1.479019945774904*ph_skin[15])-1.25*ph_skin[13]-0.9682458365518541*ph_skin[11]-0.5590169943749473*ph_skin[8])*cchi*jacob_skin_inv+chi*((-(2.958039891549808*ey_skin[15])-1.9364916731037083*ey_skin[11])*jacob_skin_inv+(1.479019945774904*ey_edge[15]+1.25*ey_edge[13]+0.9682458365518541*ey_edge[11]+0.5590169943749473*ey_edge[8])*jacob_edge_inv); 
  incr[14] = (-(1.75*ph_skin[14])-1.479019945774904*ph_skin[10]-1.14564392373896*ph_skin[6]-0.6614378277661477*ph_skin[4])*cchi*jacob_skin_inv+chi*((-(2.958039891549808*ey_skin[10])-1.3228756555322954*ey_skin[4])*jacob_skin_inv+(-(1.75*ey_edge[14])-1.479019945774904*ey_edge[10]-1.14564392373896*ey_edge[6]-0.6614378277661477*ey_edge[4])*jacob_edge_inv); 
  incr[15] = (-(1.75*ph_skin[15])-1.479019945774904*ph_skin[13]-1.14564392373896*ph_skin[11]-0.6614378277661477*ph_skin[8])*cchi*jacob_skin_inv+chi*((-(2.958039891549808*ey_skin[13])-1.3228756555322954*ey_skin[8])*jacob_skin_inv+(-(1.75*ey_edge[15])-1.479019945774904*ey_edge[13]-1.14564392373896*ey_edge[11]-0.6614378277661477*ey_edge[8])*jacob_edge_inv); 

  outPh[0] += incr[0]*dx1; 
  outPh[1] += incr[1]*dx1; 
  outPh[2] += incr[2]*dx1; 
  outPh[3] += incr[3]*dx1; 
  outPh[4] += incr[4]*dx1; 
  outPh[5] += incr[5]*dx1; 
  outPh[6] += incr[6]*dx1; 
  outPh[7] += incr[7]*dx1; 
  outPh[8] += incr[8]*dx1; 
  outPh[9] += incr[9]*dx1; 
  outPh[10] += incr[10]*dx1; 
  outPh[11] += incr[11]*dx1; 
  outPh[12] += incr[12]*dx1; 
  outPh[13] += incr[13]*dx1; 
  outPh[14] += incr[14]*dx1; 
  outPh[15] += incr[15]*dx1; 

  incr[0] = (-(0.6614378277661477*ps_skin[9])-0.5590169943749475*ps_skin[5]-0.4330127018922193*ps_skin[2]-0.25*ps_skin[0])*cgamma*jacob_skin_inv+c2gamma*((-(1.3228756555322954*by_skin[9])-0.8660254037844386*by_skin[2])*jacob_skin_inv+(0.6614378277661477*by_edge[9]+0.5590169943749475*by_edge[5]+0.4330127018922193*by_edge[2]+0.25*by_edge[0])*jacob_edge_inv); 
  incr[1] = (-(0.6614378277661477*ps_skin[12])-0.5590169943749476*ps_skin[7]-0.4330127018922193*ps_skin[3]-0.25*ps_skin[1])*cgamma*jacob_skin_inv+c2gamma*((-(1.3228756555322954*by_skin[12])-0.8660254037844386*by_skin[3])*jacob_skin_inv+(0.6614378277661477*by_edge[12]+0.5590169943749476*by_edge[7]+0.4330127018922193*by_edge[3]+0.25*by_edge[1])*jacob_edge_inv); 
  incr[2] = (-(1.14564392373896*ps_skin[9])-0.9682458365518543*ps_skin[5]-0.75*ps_skin[2]-0.4330127018922193*ps_skin[0])*cgamma*jacob_skin_inv+c2gamma*((-(1.9364916731037085*by_skin[5])-0.8660254037844386*by_skin[0])*jacob_skin_inv+(-(1.14564392373896*by_edge[9])-0.9682458365518543*by_edge[5]-0.75*by_edge[2]-0.4330127018922193*by_edge[0])*jacob_edge_inv); 
  incr[3] = (-(1.14564392373896*ps_skin[12])-0.9682458365518543*ps_skin[7]-0.75*ps_skin[3]-0.4330127018922193*ps_skin[1])*cgamma*jacob_skin_inv+c2gamma*((-(1.9364916731037085*by_skin[7])-0.8660254037844386*by_skin[1])*jacob_skin_inv+(-(1.14564392373896*by_edge[12])-0.9682458365518543*by_edge[7]-0.75*by_edge[3]-0.4330127018922193*by_edge[1])*jacob_edge_inv); 
  incr[4] = (-(0.6614378277661477*ps_skin[14])-0.5590169943749475*ps_skin[10]-0.43301270189221935*ps_skin[6]-0.25*ps_skin[4])*cgamma*jacob_skin_inv+c2gamma*((-(1.3228756555322954*by_skin[14])-0.8660254037844387*by_skin[6])*jacob_skin_inv+(0.6614378277661477*by_edge[14]+0.5590169943749475*by_edge[10]+0.43301270189221935*by_edge[6]+0.25*by_edge[4])*jacob_edge_inv); 
  incr[5] = (-(1.4790199457749043*ps_skin[9])-1.25*ps_skin[5]-0.9682458365518543*ps_skin[2]-0.5590169943749475*ps_skin[0])*cgamma*jacob_skin_inv+c2gamma*((-(2.9580398915498085*by_skin[9])-1.9364916731037085*by_skin[2])*jacob_skin_inv+(1.4790199457749043*by_edge[9]+1.25*by_edge[5]+0.9682458365518543*by_edge[2]+0.5590169943749475*by_edge[0])*jacob_edge_inv); 
  incr[6] = (-(1.14564392373896*ps_skin[14])-0.9682458365518543*ps_skin[10]-0.75*ps_skin[6]-0.43301270189221935*ps_skin[4])*cgamma*jacob_skin_inv+c2gamma*((-(1.9364916731037085*by_skin[10])-0.8660254037844387*by_skin[4])*jacob_skin_inv+(-(1.14564392373896*by_edge[14])-0.9682458365518543*by_edge[10]-0.75*by_edge[6]-0.43301270189221935*by_edge[4])*jacob_edge_inv); 
  incr[7] = (-(1.4790199457749038*ps_skin[12])-1.25*ps_skin[7]-0.9682458365518543*ps_skin[3]-0.5590169943749476*ps_skin[1])*cgamma*jacob_skin_inv+c2gamma*((-(2.9580398915498076*by_skin[12])-1.9364916731037085*by_skin[3])*jacob_skin_inv+(1.4790199457749038*by_edge[12]+1.25*by_edge[7]+0.9682458365518543*by_edge[3]+0.5590169943749476*by_edge[1])*jacob_edge_inv); 
  incr[8] = (-(0.6614378277661477*ps_skin[15])-0.5590169943749473*ps_skin[13]-0.4330127018922193*ps_skin[11]-0.25*ps_skin[8])*cgamma*jacob_skin_inv+c2gamma*((-(1.3228756555322954*by_skin[15])-0.8660254037844386*by_skin[11])*jacob_skin_inv+(0.6614378277661477*by_edge[15]+0.5590169943749473*by_edge[13]+0.4330127018922193*by_edge[11]+0.25*by_edge[8])*jacob_edge_inv); 
  incr[9] = (-(1.75*ps_skin[9])-1.4790199457749043*ps_skin[5]-1.14564392373896*ps_skin[2]-0.6614378277661477*ps_skin[0])*cgamma*jacob_skin_inv+c2gamma*((-(2.9580398915498085*by_skin[5])-1.3228756555322954*by_skin[0])*jacob_skin_inv+(-(1.75*by_edge[9])-1.4790199457749043*by_edge[5]-1.14564392373896*by_edge[2]-0.6614378277661477*by_edge[0])*jacob_edge_inv); 
  incr[10] = (-(1.479019945774904*ps_skin[14])-1.25*ps_skin[10]-0.9682458365518543*ps_skin[6]-0.5590169943749475*ps_skin[4])*cgamma*jacob_skin_inv+c2gamma*((-(2.958039891549808*by_skin[14])-1.9364916731037085*by_skin[6])*jacob_skin_inv+(1.479019945774904*by_edge[14]+1.25*by_edge[10]+0.9682458365518543*by_edge[6]+0.5590169943749475*by_edge[4])*jacob_edge_inv); 
  incr[11] = (-(1.14564392373896*ps_skin[15])-0.9682458365518541*ps_skin[13]-0.75*ps_skin[11]-0.4330127018922193*ps_skin[8])*cgamma*jacob_skin_inv+c2gamma*((-(1.9364916731037083*by_skin[13])-0.8660254037844386*by_skin[8])*jacob_skin_inv+(-(1.14564392373896*by_edge[15])-0.9682458365518541*by_edge[13]-0.75*by_edge[11]-0.4330127018922193*by_edge[8])*jacob_edge_inv); 
  incr[12] = (-(1.75*ps_skin[12])-1.4790199457749038*ps_skin[7]-1.14564392373896*ps_skin[3]-0.6614378277661477*ps_skin[1])*cgamma*jacob_skin_inv+c2gamma*((-(2.9580398915498076*by_skin[7])-1.3228756555322954*by_skin[1])*jacob_skin_inv+(-(1.75*by_edge[12])-1.4790199457749038*by_edge[7]-1.14564392373896*by_edge[3]-0.6614378277661477*by_edge[1])*jacob_edge_inv); 
  incr[13] = (-(1.479019945774904*ps_skin[15])-1.25*ps_skin[13]-0.9682458365518541*ps_skin[11]-0.5590169943749473*ps_skin[8])*cgamma*jacob_skin_inv+c2gamma*((-(2.958039891549808*by_skin[15])-1.9364916731037083*by_skin[11])*jacob_skin_inv+(1.479019945774904*by_edge[15]+1.25*by_edge[13]+0.9682458365518541*by_edge[11]+0.5590169943749473*by_edge[8])*jacob_edge_inv); 
  incr[14] = (-(1.75*ps_skin[14])-1.479019945774904*ps_skin[10]-1.14564392373896*ps_skin[6]-0.6614378277661477*ps_skin[4])*cgamma*jacob_skin_inv+c2gamma*((-(2.958039891549808*by_skin[10])-1.3228756555322954*by_skin[4])*jacob_skin_inv+(-(1.75*by_edge[14])-1.479019945774904*by_edge[10]-1.14564392373896*by_edge[6]-0.6614378277661477*by_edge[4])*jacob_edge_inv); 
  incr[15] = (-(1.75*ps_skin[15])-1.479019945774904*ps_skin[13]-1.14564392373896*ps_skin[11]-0.6614378277661477*ps_skin[8])*cgamma*jacob_skin_inv+c2gamma*((-(2.958039891549808*by_skin[13])-1.3228756555322954*by_skin[8])*jacob_skin_inv+(-(1.75*by_edge[15])-1.479019945774904*by_edge[13]-1.14564392373896*by_edge[11]-0.6614378277661477*by_edge[8])*jacob_edge_inv); 

  outPs[0] += incr[0]*dx1; 
  outPs[1] += incr[1]*dx1; 
  outPs[2] += incr[2]*dx1; 
  outPs[3] += incr[3]*dx1; 
  outPs[4] += incr[4]*dx1; 
  outPs[5] += incr[5]*dx1; 
  outPs[6] += incr[6]*dx1; 
  outPs[7] += incr[7]*dx1; 
  outPs[8] += incr[8]*dx1; 
  outPs[9] += incr[9]*dx1; 
  outPs[10] += incr[10]*dx1; 
  outPs[11] += incr[11]*dx1; 
  outPs[12] += incr[12]*dx1; 
  outPs[13] += incr[13]*dx1; 
  outPs[14] += incr[14]*dx1; 
  outPs[15] += incr[15]*dx1; 

  } 

  return 0.;

} 
