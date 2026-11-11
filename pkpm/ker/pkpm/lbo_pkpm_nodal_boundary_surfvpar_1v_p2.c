#include <gkyl_lbo_pkpm_nodal_kernels.h> 
GKYL_CU_DH void lbo_pkpm_nodal_boundary_surfvpar_1v_p2(const double *w, const double *dxv, double nu, double nu_u, double nu_vtsq, int stride, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[1]:    Velocity cell centre. 
  // dxv[1]:  Velocity cell width. 
  // nu:      Collision frequency at the node. 
  // nu_u:    Collision frequency times the drag velocity at the node. 
  // nu_vtsq: Collision frequency times the thermal speed squared at the node. 
  // stride:  Offset between the blocks [F_0, G_1, M0_perp (3), N_perp (3)] of a cell. 
  // edge:    -1 for the cell at the lower end of the velocity grid, 1 for the cell at the upper end. 
  // fedge/fskin: Velocity expansions of the blocks at the node in the interior neighbour and in the cell (w, dxv are the cell's). 
  // out:     Incremented output in the cell. 
  const double rdv = 2.0/dxv[0]; 
  const double d2 = nu_vtsq*rdv*rdv; 
  if (edge == -1) { 
  const double af = rdv*(nu*(w[0] + 0.5*dxv[0]) - nu_u); 
  const double *fu = af < 0.0 ? fskin : fedge; 
  const double su = af < 0.0 ? 1.224744871391589 : -1.224744871391589; 
  for (int b = 0; b < 8; ++b) { 
    const double *he = &fedge[b*stride]; 
    const double *hs = &fskin[b*stride]; 
    const double *hu = &fu[b*stride]; 
    double *o = &out[b*stride]; 
    const double r0 = d2*(0.3535533905932737*hs[0] + 0.49755260400283263*hs[1] + 0.34587411908091625*hs[2] + 0.3535533905932737*he[0] - 0.49755260400283263*he[1] + 0.34587411908091625*he[2]); 
    const double r1 = af*(0.7071067811865475*hu[0] + su*hu[1] + 1.5811388300841895*hu[2]) 
      + d2*(-1.325825214724776*hs[0] - 1.6840241981634343*hs[1] - 0.9486832980505137*hs[2] + 1.325825214724776*he[0] - 1.6840241981634343*he[1] + 0.9486832980505137*he[2]); 
    const double hw = d2*(0.7071067811865475*hs[0] - 1.224744871391589*hs[1] + 1.5811388300841895*hs[2]); 
    o[0] += 0.7071067811865475*r1; 
    o[1] += 1.224744871391589*r1 - 1.224744871391589*r0 + 1.224744871391589*hw; 
    o[2] += 1.5811388300841895*r1 - 4.743416490252569*r0 - 4.743416490252569*hw; 
  } 
  } else { 
  const double af = rdv*(nu*(w[0] - 0.5*dxv[0]) - nu_u); 
  const double *fu = af < 0.0 ? fedge : fskin; 
  const double su = af < 0.0 ? 1.224744871391589 : -1.224744871391589; 
  for (int b = 0; b < 8; ++b) { 
    const double *he = &fedge[b*stride]; 
    const double *hs = &fskin[b*stride]; 
    const double *hu = &fu[b*stride]; 
    double *o = &out[b*stride]; 
    const double r0 = d2*(0.3535533905932737*he[0] + 0.49755260400283263*he[1] + 0.34587411908091625*he[2] + 0.3535533905932737*hs[0] - 0.49755260400283263*hs[1] + 0.34587411908091625*hs[2]); 
    const double r1 = af*(0.7071067811865475*hu[0] + su*hu[1] + 1.5811388300841895*hu[2]) 
      + d2*(-1.325825214724776*he[0] - 1.6840241981634343*he[1] - 0.9486832980505137*he[2] + 1.325825214724776*hs[0] - 1.6840241981634343*hs[1] + 0.9486832980505137*hs[2]); 
    const double hw = d2*(0.7071067811865475*hs[0] + 1.224744871391589*hs[1] + 1.5811388300841895*hs[2]); 
    o[0] += -0.7071067811865475*r1; 
    o[1] += 1.224744871391589*r1 + 1.224744871391589*r0 - 1.224744871391589*hw; 
    o[2] += -1.5811388300841895*r1 - 4.743416490252569*r0 - 4.743416490252569*hw; 
  } 
  } 
} 
