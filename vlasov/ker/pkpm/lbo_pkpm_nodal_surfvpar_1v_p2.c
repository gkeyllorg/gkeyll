#include <gkyl_lbo_pkpm_nodal_kernels.h> 
GKYL_CU_DH void lbo_pkpm_nodal_surfvpar_1v_p2(const double *w, const double *dxv, double nu, double nu_u, double nu_vtsq, int stride, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[1]:    Velocity cell centre. 
  // dxv[1]:  Velocity cell width. 
  // nu:      Collision frequency at the node. 
  // nu_u:    Collision frequency times the drag velocity at the node. 
  // nu_vtsq: Collision frequency times the thermal speed squared at the node. 
  // stride:  Offset between the blocks [F_0, G_1, M0_perp (3), N_perp (3)] of a cell. 
  // fl/fc/fr: Velocity expansions of the blocks at the node in the cell below, the cell and the cell above. 
  // out:     Incremented output in the cell. 
  const double rdv = 2.0/dxv[0]; 
  const double d2 = nu_vtsq*rdv*rdv; 
  // Drag speed (times 2/dv) at the lower and upper faces; the upwind cell is the one above a face where it is >= 0. 
  const double af_l = rdv*(nu*(w[0] - 0.5*dxv[0]) - nu_u); 
  const double af_r = rdv*(nu*(w[0] + 0.5*dxv[0]) - nu_u); 
  const double *fu_l = af_l < 0.0 ? fl : fc; 
  const double su_l = af_l < 0.0 ? 1.224744871391589 : -1.224744871391589; 
  const double *fu_r = af_r < 0.0 ? fc : fr; 
  const double su_r = af_r < 0.0 ? 1.224744871391589 : -1.224744871391589; 
  for (int b = 0; b < 8; ++b) { 
    const double *hl = &fl[b*stride]; 
    const double *hc = &fc[b*stride]; 
    const double *hr = &fr[b*stride]; 
    const double *hu_l = &fu_l[b*stride]; 
    const double *hu_r = &fu_r[b*stride]; 
    double *o = &out[b*stride]; 
    const double r0_l = d2*(0.3535533905932737*hl[0] + 0.49755260400283263*hl[1] + 0.34587411908091625*hl[2] + 0.3535533905932737*hc[0] - 0.49755260400283263*hc[1] + 0.34587411908091625*hc[2]); 
    const double r1_l = af_l*(0.7071067811865475*hu_l[0] + su_l*hu_l[1] + 1.5811388300841895*hu_l[2]) 
      + d2*(-1.325825214724776*hl[0] - 1.6840241981634343*hl[1] - 0.9486832980505137*hl[2] + 1.325825214724776*hc[0] - 1.6840241981634343*hc[1] + 0.9486832980505137*hc[2]); 
    const double r0_r = d2*(0.3535533905932737*hc[0] + 0.49755260400283263*hc[1] + 0.34587411908091625*hc[2] + 0.3535533905932737*hr[0] - 0.49755260400283263*hr[1] + 0.34587411908091625*hr[2]); 
    const double r1_r = af_r*(0.7071067811865475*hu_r[0] + su_r*hu_r[1] + 1.5811388300841895*hu_r[2]) 
      + d2*(-1.325825214724776*hc[0] - 1.6840241981634343*hc[1] - 0.9486832980505137*hc[2] + 1.325825214724776*hr[0] - 1.6840241981634343*hr[1] + 0.9486832980505137*hr[2]); 
    o[0] += -0.7071067811865475*r1_l + 0.7071067811865475*r1_r; 
    o[1] += 1.224744871391589*r1_l + 1.224744871391589*r0_l + 1.224744871391589*r1_r - 1.224744871391589*r0_r; 
    o[2] += -1.5811388300841895*r1_l - 4.743416490252569*r0_l + 1.5811388300841895*r1_r - 4.743416490252569*r0_r; 
  } 
} 
