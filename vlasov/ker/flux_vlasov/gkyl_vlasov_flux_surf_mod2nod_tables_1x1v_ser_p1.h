#pragma once 
// Generated nodal surface-evaluation (MODAL->NODAL) tables for the 
// alpha_quad producers and lax_flux_nodal kernels (see 
// surf-eval-tables.mac). Header-static so both host and device 
// instantiations of GKYL_CU_DH kernels can read them; unused arrays 
// are elided per TU. 
#include <gkyl_util.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_meta.h> 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_conf_ev[4] = { 
  0.7071067811865475, 
  -0.7071067811865475, 
  0.7071067811865475, 
  0.7071067811865475, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_conf_dx0[4] = { 
  0.0, 
  1.224744871391589, 
  0.0, 
  1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_vel_ev_v0[2] = { 
  0.7071067811865475, 
  -1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_vel_sparse_idx[2] = { 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_v0_Cm[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_v0_CmD[4] = { 
  0.0, 
  2.5, 
  0.0, 
  2.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_v0_CmDx0[4] = { 
  0.0, 
  2.5, 
  0.0, 
  2.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_v0_V[2] = { 
  1.0, 
  2.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_v0_coefl[4] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_v0_coefr[4] = { 
  0.5, 
  0.34641016151377546, 
  -0.34641016151377546, 
  -0.24, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_ph_v0_cmap[4] = { 
  0, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_ph_v0_vlmap[4] = { 
  0, 
  0, 
  1, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_ph_v0_vrmap[4] = { 
  0, 
  0, 
  1, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_v0_Vd0[2] = { 
  0.0, 
  2.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_ph_v0_vld0map[4] = { 
  0, 
  0, 
  1, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_ph_v0_vrd0map[4] = { 
  0, 
  0, 
  1, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_v0_dcoefl0[4] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_v0_dcoefr0[4] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x1v_ser_p1_ph_v0_meta = { .nb = 4, .na = 2, .nmv = 2, .no = 2, .ni = 1, .ns = 2, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_vel_vol_dv0[4] = { 
  0.0, 
  1.224744871391589, 
  0.0, 
  1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_confsurf_x0_ev_l[2] = { 
  0.7071067811865475, 
  1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_confsurf_x0_ev_r[2] = { 
  0.7071067811865475, 
  -1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_x0_Cm[1] = { 
  1.0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_x0_V[8] = { 
  1.0, 
  2.5, 
  -1.4433756729740645, 
  -3.6084391824351614, 
  1.0, 
  2.5, 
  1.4433756729740645, 
  3.6084391824351614, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_x0_coefl[4] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_x0_coefr[4] = { 
  0.5, 
  -0.34641016151377546, 
  0.34641016151377546, 
  -0.24, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_ph_x0_cmap[4] = { 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_ph_x0_vlmap[4] = { 
  0, 
  1, 
  2, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_ph_x0_vrmap[4] = { 
  0, 
  1, 
  2, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_x0_Vd0[6] = { 
  0.0, 
  2.5, 
  6.25, 
  0.0, 
  2.5, 
  6.25, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_ph_x0_vld0map[4] = { 
  0, 
  0, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_ph_x0_vrd0map[4] = { 
  0, 
  0, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_x0_dcoefl0[4] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_ph_x0_dcoefr0[4] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  -0.24, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x1v_ser_p1_ph_x0_meta = { .nb = 4, .na = 1, .nmv = 4, .no = 1, .ni = 2, .ns = 2, .fhat_off = 0 }; 
