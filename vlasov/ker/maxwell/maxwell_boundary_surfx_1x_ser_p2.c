#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH double maxwell_boundary_surfx_1x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos_edge, const double *jacob_pos_skin, const int edge, const double *qedge, const double *qskin, double* GKYL_RESTRICT out) 
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
 
  const double dx1 = 2.0/dx[0]; 
  const double jacob_edge_inv = 1.0/jacob_pos_edge[0]; 
  const double jacob_skin_inv = 1.0/jacob_pos_skin[0]; 
 
  const double *ex_edge = &qedge[0]; 
  const double *ey_edge = &qedge[3]; 
  const double *ez_edge = &qedge[6]; 
  const double *bx_edge = &qedge[9]; 
  const double *by_edge = &qedge[12]; 
  const double *bz_edge = &qedge[15]; 
  const double *ph_edge = &qedge[18]; 
  const double *ps_edge = &qedge[21]; 

  const double *ex_skin = &qskin[0]; 
  const double *ey_skin = &qskin[3]; 
  const double *ez_skin = &qskin[6]; 
  const double *bx_skin = &qskin[9]; 
  const double *by_skin = &qskin[12]; 
  const double *bz_skin = &qskin[15]; 
  const double *ph_skin = &qskin[18]; 
  const double *ps_skin = &qskin[21]; 

  double *outEx = &out[0]; 
  double *outEy = &out[3]; 
  double *outEz = &out[6]; 
  double *outBx = &out[9]; 
  double *outBy = &out[12]; 
  double *outBz = &out[15]; 
  double *outPh = &out[18]; 
  double *outPs = &out[21]; 

  double incr[3]; 
 
  if (edge == -1) { 

  incr[0] = (-(0.5590169943749475*ex_skin[2])+0.4330127018922193*ex_skin[1]-0.25*ex_skin[0])*cchi*jacob_skin_inv+c2chi*((-(0.5590169943749475*ph_edge[2])+0.4330127018922193*ph_edge[1]-0.25*ph_edge[0])*jacob_edge_inv-0.8660254037844386*ph_skin[1]*jacob_skin_inv); 
  incr[1] = c2chi*((-(1.9364916731037085*ph_skin[2])-0.8660254037844386*ph_skin[0])*jacob_skin_inv+(-(0.9682458365518543*ph_edge[2])+0.75*ph_edge[1]-0.4330127018922193*ph_edge[0])*jacob_edge_inv)+(0.9682458365518543*ex_skin[2]-0.75*ex_skin[1]+0.4330127018922193*ex_skin[0])*cchi*jacob_skin_inv; 
  incr[2] = (-(1.25*ex_skin[2])+0.9682458365518543*ex_skin[1]-0.5590169943749475*ex_skin[0])*cchi*jacob_skin_inv+c2chi*((-(1.25*ph_edge[2])+0.9682458365518543*ph_edge[1]-0.5590169943749475*ph_edge[0])*jacob_edge_inv-1.9364916731037085*ph_skin[1]*jacob_skin_inv); 

  outEx[0] += incr[0]*dx1; 
  outEx[1] += incr[1]*dx1; 
  outEx[2] += incr[2]*dx1; 

  incr[0] = (-(0.5590169943749475*ey_skin[2])+0.4330127018922193*ey_skin[1]-0.25*ey_skin[0])*c*jacob_skin_inv+c2*((-(0.5590169943749475*bz_edge[2])+0.4330127018922193*bz_edge[1]-0.25*bz_edge[0])*jacob_edge_inv-0.8660254037844386*bz_skin[1]*jacob_skin_inv); 
  incr[1] = (0.9682458365518543*ey_skin[2]-0.75*ey_skin[1]+0.4330127018922193*ey_skin[0])*c*jacob_skin_inv+c2*((-(1.9364916731037085*bz_skin[2])-0.8660254037844386*bz_skin[0])*jacob_skin_inv+(-(0.9682458365518543*bz_edge[2])+0.75*bz_edge[1]-0.4330127018922193*bz_edge[0])*jacob_edge_inv); 
  incr[2] = (-(1.25*ey_skin[2])+0.9682458365518543*ey_skin[1]-0.5590169943749475*ey_skin[0])*c*jacob_skin_inv+c2*((-(1.25*bz_edge[2])+0.9682458365518543*bz_edge[1]-0.5590169943749475*bz_edge[0])*jacob_edge_inv-1.9364916731037085*bz_skin[1]*jacob_skin_inv); 

  outEy[0] += incr[0]*dx1; 
  outEy[1] += incr[1]*dx1; 
  outEy[2] += incr[2]*dx1; 

  incr[0] = (-(0.5590169943749475*ez_skin[2])+0.4330127018922193*ez_skin[1]-0.25*ez_skin[0])*c*jacob_skin_inv+c2*(0.8660254037844386*by_skin[1]*jacob_skin_inv+(0.5590169943749475*by_edge[2]-0.4330127018922193*by_edge[1]+0.25*by_edge[0])*jacob_edge_inv); 
  incr[1] = (0.9682458365518543*ez_skin[2]-0.75*ez_skin[1]+0.4330127018922193*ez_skin[0])*c*jacob_skin_inv+c2*((1.9364916731037085*by_skin[2]+0.8660254037844386*by_skin[0])*jacob_skin_inv+(0.9682458365518543*by_edge[2]-0.75*by_edge[1]+0.4330127018922193*by_edge[0])*jacob_edge_inv); 
  incr[2] = (-(1.25*ez_skin[2])+0.9682458365518543*ez_skin[1]-0.5590169943749475*ez_skin[0])*c*jacob_skin_inv+c2*(1.9364916731037085*by_skin[1]*jacob_skin_inv+(1.25*by_edge[2]-0.9682458365518543*by_edge[1]+0.5590169943749475*by_edge[0])*jacob_edge_inv); 

  outEz[0] += incr[0]*dx1; 
  outEz[1] += incr[1]*dx1; 
  outEz[2] += incr[2]*dx1; 

  incr[0] = ((-(0.5590169943749475*ps_edge[2])+0.4330127018922193*ps_edge[1]-0.25*ps_edge[0])*jacob_edge_inv-0.8660254037844386*ps_skin[1]*jacob_skin_inv)*gamma+(-(0.5590169943749475*bx_skin[2])+0.4330127018922193*bx_skin[1]-0.25*bx_skin[0])*cgamma*jacob_skin_inv; 
  incr[1] = ((-(1.9364916731037085*ps_skin[2])-0.8660254037844386*ps_skin[0])*jacob_skin_inv+(-(0.9682458365518543*ps_edge[2])+0.75*ps_edge[1]-0.4330127018922193*ps_edge[0])*jacob_edge_inv)*gamma+(0.9682458365518543*bx_skin[2]-0.75*bx_skin[1]+0.4330127018922193*bx_skin[0])*cgamma*jacob_skin_inv; 
  incr[2] = ((-(1.25*ps_edge[2])+0.9682458365518543*ps_edge[1]-0.5590169943749475*ps_edge[0])*jacob_edge_inv-1.9364916731037085*ps_skin[1]*jacob_skin_inv)*gamma+(-(1.25*bx_skin[2])+0.9682458365518543*bx_skin[1]-0.5590169943749475*bx_skin[0])*cgamma*jacob_skin_inv; 

  outBx[0] += incr[0]*dx1; 
  outBx[1] += incr[1]*dx1; 
  outBx[2] += incr[2]*dx1; 

  incr[0] = ((-(0.5590169943749475*by_skin[2])+0.4330127018922193*by_skin[1]-0.25*by_skin[0])*c+0.8660254037844386*ez_skin[1])*jacob_skin_inv+(0.5590169943749475*ez_edge[2]-0.4330127018922193*ez_edge[1]+0.25*ez_edge[0])*jacob_edge_inv; 
  incr[1] = ((0.9682458365518543*by_skin[2]-0.75*by_skin[1]+0.4330127018922193*by_skin[0])*c+1.9364916731037085*ez_skin[2]+0.8660254037844386*ez_skin[0])*jacob_skin_inv+(0.9682458365518543*ez_edge[2]-0.75*ez_edge[1]+0.4330127018922193*ez_edge[0])*jacob_edge_inv; 
  incr[2] = ((-(1.25*by_skin[2])+0.9682458365518543*by_skin[1]-0.5590169943749475*by_skin[0])*c+1.9364916731037085*ez_skin[1])*jacob_skin_inv+(1.25*ez_edge[2]-0.9682458365518543*ez_edge[1]+0.5590169943749475*ez_edge[0])*jacob_edge_inv; 

  outBy[0] += incr[0]*dx1; 
  outBy[1] += incr[1]*dx1; 
  outBy[2] += incr[2]*dx1; 

  incr[0] = ((-(0.5590169943749475*bz_skin[2])+0.4330127018922193*bz_skin[1]-0.25*bz_skin[0])*c-0.8660254037844386*ey_skin[1])*jacob_skin_inv+(-(0.5590169943749475*ey_edge[2])+0.4330127018922193*ey_edge[1]-0.25*ey_edge[0])*jacob_edge_inv; 
  incr[1] = ((0.9682458365518543*bz_skin[2]-0.75*bz_skin[1]+0.4330127018922193*bz_skin[0])*c-1.9364916731037085*ey_skin[2]-0.8660254037844386*ey_skin[0])*jacob_skin_inv+(-(0.9682458365518543*ey_edge[2])+0.75*ey_edge[1]-0.4330127018922193*ey_edge[0])*jacob_edge_inv; 
  incr[2] = ((-(1.25*bz_skin[2])+0.9682458365518543*bz_skin[1]-0.5590169943749475*bz_skin[0])*c-1.9364916731037085*ey_skin[1])*jacob_skin_inv+(-(1.25*ey_edge[2])+0.9682458365518543*ey_edge[1]-0.5590169943749475*ey_edge[0])*jacob_edge_inv; 

  outBz[0] += incr[0]*dx1; 
  outBz[1] += incr[1]*dx1; 
  outBz[2] += incr[2]*dx1; 

  incr[0] = (-(0.5590169943749475*ph_skin[2])+0.4330127018922193*ph_skin[1]-0.25*ph_skin[0])*cchi*jacob_skin_inv+chi*((-(0.5590169943749475*ex_edge[2])+0.4330127018922193*ex_edge[1]-0.25*ex_edge[0])*jacob_edge_inv-0.8660254037844386*ex_skin[1]*jacob_skin_inv); 
  incr[1] = (0.9682458365518543*ph_skin[2]-0.75*ph_skin[1]+0.4330127018922193*ph_skin[0])*cchi*jacob_skin_inv+chi*((-(1.9364916731037085*ex_skin[2])-0.8660254037844386*ex_skin[0])*jacob_skin_inv+(-(0.9682458365518543*ex_edge[2])+0.75*ex_edge[1]-0.4330127018922193*ex_edge[0])*jacob_edge_inv); 
  incr[2] = (-(1.25*ph_skin[2])+0.9682458365518543*ph_skin[1]-0.5590169943749475*ph_skin[0])*cchi*jacob_skin_inv+chi*((-(1.25*ex_edge[2])+0.9682458365518543*ex_edge[1]-0.5590169943749475*ex_edge[0])*jacob_edge_inv-1.9364916731037085*ex_skin[1]*jacob_skin_inv); 

  outPh[0] += incr[0]*dx1; 
  outPh[1] += incr[1]*dx1; 
  outPh[2] += incr[2]*dx1; 

  incr[0] = (-(0.5590169943749475*ps_skin[2])+0.4330127018922193*ps_skin[1]-0.25*ps_skin[0])*cgamma*jacob_skin_inv+c2gamma*((-(0.5590169943749475*bx_edge[2])+0.4330127018922193*bx_edge[1]-0.25*bx_edge[0])*jacob_edge_inv-0.8660254037844386*bx_skin[1]*jacob_skin_inv); 
  incr[1] = (0.9682458365518543*ps_skin[2]-0.75*ps_skin[1]+0.4330127018922193*ps_skin[0])*cgamma*jacob_skin_inv+c2gamma*((-(1.9364916731037085*bx_skin[2])-0.8660254037844386*bx_skin[0])*jacob_skin_inv+(-(0.9682458365518543*bx_edge[2])+0.75*bx_edge[1]-0.4330127018922193*bx_edge[0])*jacob_edge_inv); 
  incr[2] = (-(1.25*ps_skin[2])+0.9682458365518543*ps_skin[1]-0.5590169943749475*ps_skin[0])*cgamma*jacob_skin_inv+c2gamma*((-(1.25*bx_edge[2])+0.9682458365518543*bx_edge[1]-0.5590169943749475*bx_edge[0])*jacob_edge_inv-1.9364916731037085*bx_skin[1]*jacob_skin_inv); 

  outPs[0] += incr[0]*dx1; 
  outPs[1] += incr[1]*dx1; 
  outPs[2] += incr[2]*dx1; 

  } else { 

  incr[0] = (-(0.5590169943749475*ex_skin[2])-0.4330127018922193*ex_skin[1]-0.25*ex_skin[0])*cchi*jacob_skin_inv+c2chi*((0.5590169943749475*ph_edge[2]+0.4330127018922193*ph_edge[1]+0.25*ph_edge[0])*jacob_edge_inv-0.8660254037844386*ph_skin[1]*jacob_skin_inv); 
  incr[1] = c2chi*((-(1.9364916731037085*ph_skin[2])-0.8660254037844386*ph_skin[0])*jacob_skin_inv+(-(0.9682458365518543*ph_edge[2])-0.75*ph_edge[1]-0.4330127018922193*ph_edge[0])*jacob_edge_inv)+(-(0.9682458365518543*ex_skin[2])-0.75*ex_skin[1]-0.4330127018922193*ex_skin[0])*cchi*jacob_skin_inv; 
  incr[2] = (-(1.25*ex_skin[2])-0.9682458365518543*ex_skin[1]-0.5590169943749475*ex_skin[0])*cchi*jacob_skin_inv+c2chi*((1.25*ph_edge[2]+0.9682458365518543*ph_edge[1]+0.5590169943749475*ph_edge[0])*jacob_edge_inv-1.9364916731037085*ph_skin[1]*jacob_skin_inv); 

  outEx[0] += incr[0]*dx1; 
  outEx[1] += incr[1]*dx1; 
  outEx[2] += incr[2]*dx1; 

  incr[0] = (-(0.5590169943749475*ey_skin[2])-0.4330127018922193*ey_skin[1]-0.25*ey_skin[0])*c*jacob_skin_inv+c2*((0.5590169943749475*bz_edge[2]+0.4330127018922193*bz_edge[1]+0.25*bz_edge[0])*jacob_edge_inv-0.8660254037844386*bz_skin[1]*jacob_skin_inv); 
  incr[1] = (-(0.9682458365518543*ey_skin[2])-0.75*ey_skin[1]-0.4330127018922193*ey_skin[0])*c*jacob_skin_inv+c2*((-(1.9364916731037085*bz_skin[2])-0.8660254037844386*bz_skin[0])*jacob_skin_inv+(-(0.9682458365518543*bz_edge[2])-0.75*bz_edge[1]-0.4330127018922193*bz_edge[0])*jacob_edge_inv); 
  incr[2] = (-(1.25*ey_skin[2])-0.9682458365518543*ey_skin[1]-0.5590169943749475*ey_skin[0])*c*jacob_skin_inv+c2*((1.25*bz_edge[2]+0.9682458365518543*bz_edge[1]+0.5590169943749475*bz_edge[0])*jacob_edge_inv-1.9364916731037085*bz_skin[1]*jacob_skin_inv); 

  outEy[0] += incr[0]*dx1; 
  outEy[1] += incr[1]*dx1; 
  outEy[2] += incr[2]*dx1; 

  incr[0] = (-(0.5590169943749475*ez_skin[2])-0.4330127018922193*ez_skin[1]-0.25*ez_skin[0])*c*jacob_skin_inv+c2*(0.8660254037844386*by_skin[1]*jacob_skin_inv+(-(0.5590169943749475*by_edge[2])-0.4330127018922193*by_edge[1]-0.25*by_edge[0])*jacob_edge_inv); 
  incr[1] = (-(0.9682458365518543*ez_skin[2])-0.75*ez_skin[1]-0.4330127018922193*ez_skin[0])*c*jacob_skin_inv+c2*((1.9364916731037085*by_skin[2]+0.8660254037844386*by_skin[0])*jacob_skin_inv+(0.9682458365518543*by_edge[2]+0.75*by_edge[1]+0.4330127018922193*by_edge[0])*jacob_edge_inv); 
  incr[2] = (-(1.25*ez_skin[2])-0.9682458365518543*ez_skin[1]-0.5590169943749475*ez_skin[0])*c*jacob_skin_inv+c2*(1.9364916731037085*by_skin[1]*jacob_skin_inv+(-(1.25*by_edge[2])-0.9682458365518543*by_edge[1]-0.5590169943749475*by_edge[0])*jacob_edge_inv); 

  outEz[0] += incr[0]*dx1; 
  outEz[1] += incr[1]*dx1; 
  outEz[2] += incr[2]*dx1; 

  incr[0] = ((0.5590169943749475*ps_edge[2]+0.4330127018922193*ps_edge[1]+0.25*ps_edge[0])*jacob_edge_inv-0.8660254037844386*ps_skin[1]*jacob_skin_inv)*gamma+(-(0.5590169943749475*bx_skin[2])-0.4330127018922193*bx_skin[1]-0.25*bx_skin[0])*cgamma*jacob_skin_inv; 
  incr[1] = ((-(1.9364916731037085*ps_skin[2])-0.8660254037844386*ps_skin[0])*jacob_skin_inv+(-(0.9682458365518543*ps_edge[2])-0.75*ps_edge[1]-0.4330127018922193*ps_edge[0])*jacob_edge_inv)*gamma+(-(0.9682458365518543*bx_skin[2])-0.75*bx_skin[1]-0.4330127018922193*bx_skin[0])*cgamma*jacob_skin_inv; 
  incr[2] = ((1.25*ps_edge[2]+0.9682458365518543*ps_edge[1]+0.5590169943749475*ps_edge[0])*jacob_edge_inv-1.9364916731037085*ps_skin[1]*jacob_skin_inv)*gamma+(-(1.25*bx_skin[2])-0.9682458365518543*bx_skin[1]-0.5590169943749475*bx_skin[0])*cgamma*jacob_skin_inv; 

  outBx[0] += incr[0]*dx1; 
  outBx[1] += incr[1]*dx1; 
  outBx[2] += incr[2]*dx1; 

  incr[0] = ((-(0.5590169943749475*by_skin[2])-0.4330127018922193*by_skin[1]-0.25*by_skin[0])*c+0.8660254037844386*ez_skin[1])*jacob_skin_inv+(-(0.5590169943749475*ez_edge[2])-0.4330127018922193*ez_edge[1]-0.25*ez_edge[0])*jacob_edge_inv; 
  incr[1] = ((-(0.9682458365518543*by_skin[2])-0.75*by_skin[1]-0.4330127018922193*by_skin[0])*c+1.9364916731037085*ez_skin[2]+0.8660254037844386*ez_skin[0])*jacob_skin_inv+(0.9682458365518543*ez_edge[2]+0.75*ez_edge[1]+0.4330127018922193*ez_edge[0])*jacob_edge_inv; 
  incr[2] = ((-(1.25*by_skin[2])-0.9682458365518543*by_skin[1]-0.5590169943749475*by_skin[0])*c+1.9364916731037085*ez_skin[1])*jacob_skin_inv+(-(1.25*ez_edge[2])-0.9682458365518543*ez_edge[1]-0.5590169943749475*ez_edge[0])*jacob_edge_inv; 

  outBy[0] += incr[0]*dx1; 
  outBy[1] += incr[1]*dx1; 
  outBy[2] += incr[2]*dx1; 

  incr[0] = ((-(0.5590169943749475*bz_skin[2])-0.4330127018922193*bz_skin[1]-0.25*bz_skin[0])*c-0.8660254037844386*ey_skin[1])*jacob_skin_inv+(0.5590169943749475*ey_edge[2]+0.4330127018922193*ey_edge[1]+0.25*ey_edge[0])*jacob_edge_inv; 
  incr[1] = ((-(0.9682458365518543*bz_skin[2])-0.75*bz_skin[1]-0.4330127018922193*bz_skin[0])*c-1.9364916731037085*ey_skin[2]-0.8660254037844386*ey_skin[0])*jacob_skin_inv+(-(0.9682458365518543*ey_edge[2])-0.75*ey_edge[1]-0.4330127018922193*ey_edge[0])*jacob_edge_inv; 
  incr[2] = ((-(1.25*bz_skin[2])-0.9682458365518543*bz_skin[1]-0.5590169943749475*bz_skin[0])*c-1.9364916731037085*ey_skin[1])*jacob_skin_inv+(1.25*ey_edge[2]+0.9682458365518543*ey_edge[1]+0.5590169943749475*ey_edge[0])*jacob_edge_inv; 

  outBz[0] += incr[0]*dx1; 
  outBz[1] += incr[1]*dx1; 
  outBz[2] += incr[2]*dx1; 

  incr[0] = (-(0.5590169943749475*ph_skin[2])-0.4330127018922193*ph_skin[1]-0.25*ph_skin[0])*cchi*jacob_skin_inv+chi*((0.5590169943749475*ex_edge[2]+0.4330127018922193*ex_edge[1]+0.25*ex_edge[0])*jacob_edge_inv-0.8660254037844386*ex_skin[1]*jacob_skin_inv); 
  incr[1] = (-(0.9682458365518543*ph_skin[2])-0.75*ph_skin[1]-0.4330127018922193*ph_skin[0])*cchi*jacob_skin_inv+chi*((-(1.9364916731037085*ex_skin[2])-0.8660254037844386*ex_skin[0])*jacob_skin_inv+(-(0.9682458365518543*ex_edge[2])-0.75*ex_edge[1]-0.4330127018922193*ex_edge[0])*jacob_edge_inv); 
  incr[2] = (-(1.25*ph_skin[2])-0.9682458365518543*ph_skin[1]-0.5590169943749475*ph_skin[0])*cchi*jacob_skin_inv+chi*((1.25*ex_edge[2]+0.9682458365518543*ex_edge[1]+0.5590169943749475*ex_edge[0])*jacob_edge_inv-1.9364916731037085*ex_skin[1]*jacob_skin_inv); 

  outPh[0] += incr[0]*dx1; 
  outPh[1] += incr[1]*dx1; 
  outPh[2] += incr[2]*dx1; 

  incr[0] = (-(0.5590169943749475*ps_skin[2])-0.4330127018922193*ps_skin[1]-0.25*ps_skin[0])*cgamma*jacob_skin_inv+c2gamma*((0.5590169943749475*bx_edge[2]+0.4330127018922193*bx_edge[1]+0.25*bx_edge[0])*jacob_edge_inv-0.8660254037844386*bx_skin[1]*jacob_skin_inv); 
  incr[1] = (-(0.9682458365518543*ps_skin[2])-0.75*ps_skin[1]-0.4330127018922193*ps_skin[0])*cgamma*jacob_skin_inv+c2gamma*((-(1.9364916731037085*bx_skin[2])-0.8660254037844386*bx_skin[0])*jacob_skin_inv+(-(0.9682458365518543*bx_edge[2])-0.75*bx_edge[1]-0.4330127018922193*bx_edge[0])*jacob_edge_inv); 
  incr[2] = (-(1.25*ps_skin[2])-0.9682458365518543*ps_skin[1]-0.5590169943749475*ps_skin[0])*cgamma*jacob_skin_inv+c2gamma*((1.25*bx_edge[2]+0.9682458365518543*bx_edge[1]+0.5590169943749475*bx_edge[0])*jacob_edge_inv-1.9364916731037085*bx_skin[1]*jacob_skin_inv); 

  outPs[0] += incr[0]*dx1; 
  outPs[1] += incr[1]*dx1; 
  outPs[2] += incr[2]*dx1; 

  } 

  return 0.;

} 
