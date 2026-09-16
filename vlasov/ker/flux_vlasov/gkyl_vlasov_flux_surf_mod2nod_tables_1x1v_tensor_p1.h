#pragma once 
// Generated nodal surface-evaluation (MODAL->NODAL) tables for the 
// alpha_quad producers and lax_flux_nodal kernels (see 
// surf-eval-tables.mac). Header-static so both host and device 
// instantiations of GKYL_CU_DH kernels can read them; unused arrays 
// are elided per TU. 
#include <gkyl_util.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_meta.h> 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_conf_ev[4] = { 
  0.7071067811865475, 
  -0.7071067811865475, 
  0.7071067811865475, 
  0.7071067811865475, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_conf_dx0[4] = { 
  0.0, 
  1.224744871391589, 
  0.0, 
  1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_vel_ev_v0[3] = { 
  0.7071067811865475, 
  -1.224744871391589, 
  1.5811388300841895, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_vel_sparse_idx[3] = { 
  0, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_Cm[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_CmD[4] = { 
  0.0, 
  2.5, 
  0.0, 
  2.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_CmDx0[4] = { 
  0.0, 
  2.5, 
  0.0, 
  2.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_V[3] = { 
  1.0, 
  2.5, 
  -3.8461538461538463, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_coefl[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_coefr[6] = { 
  0.5, 
  0.34641016151377546, 
  -0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ph_v0_cmap[6] = { 
  0, 
  1, 
  0, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ph_v0_vlmap[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ph_v0_vrmap[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_Vd0[3] = { 
  0.0, 
  2.5, 
  -11.538461538461538, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ph_v0_vld0map[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ph_v0_vrd0map[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_dcoefl0[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_dcoefr0[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  0.2906888370749726, 
  0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x1v_tensor_p1_ph_v0_meta = { .nb = 6, .na = 2, .nmv = 3, .no = 2, .ni = 1, .ns = 3, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_vel_vol_dv0[9] = { 
  0.0, 
  1.224744871391589, 
  -3.6742346141747664, 
  0.0, 
  1.224744871391589, 
  0.0, 
  0.0, 
  1.224744871391589, 
  3.6742346141747664, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_confsurf_x0_ev_l[2] = { 
  0.7071067811865475, 
  1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_confsurf_x0_ev_r[2] = { 
  0.7071067811865475, 
  -1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_x0_Cm[1] = { 
  1.0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_x0_V[18] = { 
  1.0, 
  2.5, 
  -1.9364916731037085, 
  -4.841229182759272, 
  -1.5384615384615385, 
  -3.8461538461538463, 
  1.0, 
  2.5, 
  0.0, 
  0.0, 
  1.9230769230769231, 
  4.8076923076923075, 
  1.0, 
  2.5, 
  1.9364916731037085, 
  4.841229182759272, 
  -1.5384615384615385, 
  -3.8461538461538463, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_x0_coefl[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_x0_coefr[6] = { 
  0.5, 
  -0.34641016151377546, 
  0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ph_x0_cmap[6] = { 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ph_x0_vlmap[6] = { 
  0, 
  1, 
  2, 
  3, 
  4, 
  5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ph_x0_vrmap[6] = { 
  0, 
  1, 
  2, 
  3, 
  4, 
  5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_x0_Vd0[15] = { 
  0.0, 
  2.5, 
  6.25, 
  8.93765387586327, 
  22.344134689658176, 
  0.0, 
  2.5, 
  6.25, 
  0.0, 
  0.0, 
  0.0, 
  2.5, 
  6.25, 
  -8.93765387586327, 
  -22.344134689658176, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ph_x0_vld0map[6] = { 
  0, 
  0, 
  1, 
  2, 
  3, 
  4, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ph_x0_vrd0map[6] = { 
  0, 
  0, 
  1, 
  2, 
  3, 
  4, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_x0_dcoefl0[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_x0_dcoefr0[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x1v_tensor_p1_ph_x0_meta = { .nb = 6, .na = 1, .nmv = 6, .no = 1, .ni = 3, .ns = 3, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_conf_ev[4] = { 
  0.7071067811865475, 
  -0.7071067811865475, 
  0.7071067811865475, 
  0.7071067811865475, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_conf_dx0[4] = { 
  0.0, 
  1.224744871391589, 
  0.0, 
  1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_vel_ev_v0[3] = { 
  0.7071067811865475, 
  -1.224744871391589, 
  1.5811388300841895, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_vel_sparse_idx[3] = { 
  0, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_Cm[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_CmD[4] = { 
  0.0, 
  2.5, 
  0.0, 
  2.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_CmDx0[4] = { 
  0.0, 
  2.5, 
  0.0, 
  2.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_V[3] = { 
  1.0, 
  2.5, 
  -3.8461538461538463, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_coefl[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_coefr[6] = { 
  0.5, 
  0.34641016151377546, 
  -0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_ph_v0_cmap[6] = { 
  0, 
  1, 
  0, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_ph_v0_vlmap[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_ph_v0_vrmap[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_Vd0[3] = { 
  0.0, 
  2.5, 
  -11.538461538461538, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_ph_v0_vld0map[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_ph_v0_vrd0map[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_dcoefl0[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_dcoefr0[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  0.2906888370749726, 
  0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x1v_tensor_p1_ho_ph_v0_meta = { .nb = 6, .na = 2, .nmv = 3, .no = 2, .ni = 1, .ns = 3, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_vel_vol_dv0[12] = { 
  0.0, 
  1.224744871391589, 
  -4.084728180770504, 
  0.0, 
  1.224744871391589, 
  -1.6126716885136845, 
  0.0, 
  1.224744871391589, 
  1.6126716885136845, 
  0.0, 
  1.224744871391589, 
  4.084728180770504, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_confsurf_x0_ev_l[2] = { 
  0.7071067811865475, 
  1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_confsurf_x0_ev_r[2] = { 
  0.7071067811865475, 
  -1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_x0_Cm[1] = { 
  1.0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_x0_V[24] = { 
  1.0, 
  2.5, 
  -2.1528407789851314, 
  -5.382101947462829, 
  -2.355129310456592, 
  -5.887823276141479, 
  1.0, 
  2.5, 
  -0.8499526089621406, 
  -2.1248815224053517, 
  1.2562282115554928, 
  3.140570528888732, 
  1.0, 
  2.5, 
  0.8499526089621406, 
  2.1248815224053517, 
  1.2562282115554928, 
  3.140570528888732, 
  1.0, 
  2.5, 
  2.1528407789851314, 
  5.382101947462829, 
  -2.355129310456592, 
  -5.887823276141479, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_x0_coefl[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_x0_coefr[6] = { 
  0.5, 
  -0.34641016151377546, 
  0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_ph_x0_cmap[6] = { 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_ph_x0_vlmap[6] = { 
  0, 
  1, 
  2, 
  3, 
  4, 
  5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_ph_x0_vrmap[6] = { 
  0, 
  1, 
  2, 
  3, 
  4, 
  5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_x0_Vd0[20] = { 
  0.0, 
  2.5, 
  6.25, 
  9.936188210700607, 
  24.840470526751517, 
  0.0, 
  2.5, 
  6.25, 
  3.92285819520988, 
  9.8071454880247, 
  0.0, 
  2.5, 
  6.25, 
  -3.92285819520988, 
  -9.8071454880247, 
  0.0, 
  2.5, 
  6.25, 
  -9.936188210700607, 
  -24.840470526751517, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_ph_x0_vld0map[6] = { 
  0, 
  0, 
  1, 
  2, 
  3, 
  4, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_ph_x0_vrd0map[6] = { 
  0, 
  0, 
  1, 
  2, 
  3, 
  4, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_x0_dcoefl0[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_x0_dcoefr0[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x1v_tensor_p1_ho_ph_x0_meta = { .nb = 6, .na = 1, .nmv = 6, .no = 1, .ni = 4, .ns = 3, .fhat_off = 0 }; 
