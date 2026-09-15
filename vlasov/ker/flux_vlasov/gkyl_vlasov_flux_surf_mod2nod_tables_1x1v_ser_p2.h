#pragma once 
// Generated nodal surface-evaluation (MODAL->NODAL) tables for the 
// alpha_quad producers and lax_flux_nodal kernels (see 
// surf-eval-tables.mac). Header-static so both host and device 
// instantiations of GKYL_CU_DH kernels can read them; unused arrays 
// are elided per TU. 
#include <gkyl_util.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_meta.h> 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_conf_ev[9] = { 
  0.7071067811865475, 
  -0.9486832980505137, 
  0.6324555320336759, 
  0.7071067811865475, 
  0.0, 
  -0.7905694150420947, 
  0.7071067811865475, 
  0.9486832980505137, 
  0.6324555320336759, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_conf_dx0[9] = { 
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
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_vel_ev_v0[3] = { 
  0.7071067811865475, 
  -1.224744871391589, 
  1.5811388300841895, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_vel_sparse_idx[3] = { 
  0, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_v0_Cm[9] = { 
  1.0, 
  -1.9364916731037085, 
  -1.5384615384615385, 
  1.0, 
  0.0, 
  1.9230769230769231, 
  1.0, 
  1.9364916731037085, 
  -1.5384615384615385, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_v0_CmD[9] = { 
  0.0, 
  2.5, 
  8.93765387586327, 
  0.0, 
  2.5, 
  0.0, 
  0.0, 
  2.5, 
  -8.93765387586327, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_v0_CmDx0[9] = { 
  0.0, 
  2.5, 
  8.93765387586327, 
  0.0, 
  2.5, 
  0.0, 
  0.0, 
  2.5, 
  -8.93765387586327, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_v0_V[3] = { 
  1.0, 
  2.5, 
  -3.8461538461538463, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_v0_coefl[8] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_v0_coefr[8] = { 
  0.5, 
  0.34641016151377546, 
  -0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  0.20139513400278566, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ph_v0_cmap[8] = { 
  0, 
  1, 
  0, 
  1, 
  2, 
  0, 
  2, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ph_v0_vlmap[8] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ph_v0_vrmap[8] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_v0_Vd0[3] = { 
  0.0, 
  2.5, 
  -11.538461538461538, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ph_v0_vld0map[8] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ph_v0_vrd0map[8] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_v0_dcoefl0[8] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_v0_dcoefr0[8] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  0.2906888370749726, 
  -0.20139513400278566, 
  0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x1v_ser_p2_ph_v0_meta = { .nb = 8, .na = 3, .nmv = 3, .no = 3, .ni = 1, .ns = 3, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_vel_vol_dv0[9] = { 
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
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_confsurf_x0_ev_l[3] = { 
  0.7071067811865475, 
  1.224744871391589, 
  1.5811388300841895, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_confsurf_x0_ev_r[3] = { 
  0.7071067811865475, 
  -1.224744871391589, 
  1.5811388300841895, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_x0_Cm[1] = { 
  1.0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_x0_V[24] = { 
  1.0, 
  2.5, 
  -1.9364916731037085, 
  -4.841229182759272, 
  -3.8461538461538463, 
  -1.5384615384615385, 
  7.448044896552726, 
  -3.8461538461538463, 
  1.0, 
  2.5, 
  0.0, 
  0.0, 
  -3.8461538461538463, 
  1.9230769230769231, 
  0.0, 
  4.8076923076923075, 
  1.0, 
  2.5, 
  1.9364916731037085, 
  4.841229182759272, 
  -3.8461538461538463, 
  -1.5384615384615385, 
  -7.448044896552726, 
  -3.8461538461538463, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_x0_coefl[8] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_x0_coefr[8] = { 
  0.5, 
  -0.34641016151377546, 
  0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ph_x0_cmap[8] = { 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ph_x0_vlmap[8] = { 
  0, 
  1, 
  2, 
  3, 
  4, 
  5, 
  6, 
  7, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ph_x0_vrmap[8] = { 
  0, 
  1, 
  2, 
  3, 
  4, 
  5, 
  6, 
  7, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_x0_Vd0[18] = { 
  0.0, 
  2.5, 
  6.25, 
  8.93765387586327, 
  -9.615384615384615, 
  22.344134689658176, 
  0.0, 
  2.5, 
  6.25, 
  0.0, 
  -9.615384615384615, 
  0.0, 
  0.0, 
  2.5, 
  6.25, 
  -8.93765387586327, 
  -9.615384615384615, 
  -22.344134689658176, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ph_x0_vld0map[8] = { 
  0, 
  0, 
  1, 
  2, 
  0, 
  3, 
  4, 
  5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ph_x0_vrd0map[8] = { 
  0, 
  0, 
  1, 
  2, 
  0, 
  3, 
  4, 
  5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_x0_dcoefl0[8] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ph_x0_dcoefr0[8] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x1v_ser_p2_ph_x0_meta = { .nb = 8, .na = 1, .nmv = 8, .no = 1, .ni = 3, .ns = 3, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_conf_ev[12] = { 
  0.7071067811865475, 
  -1.0546722811938851, 
  0.9681844646844028, 
  0.7071067811865475, 
  -0.41639003950091297, 
  -0.5164305132317774, 
  0.7071067811865475, 
  0.41639003950091297, 
  -0.5164305132317774, 
  0.7071067811865475, 
  1.0546722811938851, 
  0.9681844646844028, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_conf_dx0[12] = { 
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
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_vel_ev_v0[3] = { 
  0.7071067811865475, 
  -1.224744871391589, 
  1.5811388300841895, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ho_vel_sparse_idx[3] = { 
  0, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_v0_Cm[12] = { 
  1.0, 
  -2.1528407789851314, 
  -2.355129310456592, 
  1.0, 
  -0.8499526089621406, 
  1.2562282115554928, 
  1.0, 
  0.8499526089621406, 
  1.2562282115554928, 
  1.0, 
  2.1528407789851314, 
  -2.355129310456592, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_v0_CmD[12] = { 
  0.0, 
  2.5, 
  9.936188210700607, 
  0.0, 
  2.5, 
  3.92285819520988, 
  0.0, 
  2.5, 
  -3.92285819520988, 
  0.0, 
  2.5, 
  -9.936188210700607, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_v0_CmDx0[12] = { 
  0.0, 
  2.5, 
  9.936188210700607, 
  0.0, 
  2.5, 
  3.92285819520988, 
  0.0, 
  2.5, 
  -3.92285819520988, 
  0.0, 
  2.5, 
  -9.936188210700607, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_v0_V[3] = { 
  1.0, 
  2.5, 
  -3.8461538461538463, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_v0_coefl[8] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_v0_coefr[8] = { 
  0.5, 
  0.34641016151377546, 
  -0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  0.20139513400278566, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ho_ph_v0_cmap[8] = { 
  0, 
  1, 
  0, 
  1, 
  2, 
  0, 
  2, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ho_ph_v0_vlmap[8] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ho_ph_v0_vrmap[8] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_v0_Vd0[3] = { 
  0.0, 
  2.5, 
  -11.538461538461538, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ho_ph_v0_vld0map[8] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ho_ph_v0_vrd0map[8] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_v0_dcoefl0[8] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_v0_dcoefr0[8] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  0.2906888370749726, 
  -0.20139513400278566, 
  0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x1v_ser_p2_ho_ph_v0_meta = { .nb = 8, .na = 3, .nmv = 3, .no = 4, .ni = 1, .ns = 3, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_vel_vol_dv0[12] = { 
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
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_confsurf_x0_ev_l[3] = { 
  0.7071067811865475, 
  1.224744871391589, 
  1.5811388300841895, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_confsurf_x0_ev_r[3] = { 
  0.7071067811865475, 
  -1.224744871391589, 
  1.5811388300841895, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_x0_Cm[1] = { 
  1.0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_x0_V[32] = { 
  1.0, 
  2.5, 
  -2.1528407789851314, 
  -5.382101947462829, 
  -3.8461538461538463, 
  -2.355129310456592, 
  8.280156842250506, 
  -5.887823276141479, 
  1.0, 
  2.5, 
  -0.8499526089621406, 
  -2.1248815224053517, 
  -3.8461538461538463, 
  1.2562282115554928, 
  3.269048496008233, 
  3.140570528888732, 
  1.0, 
  2.5, 
  0.8499526089621406, 
  2.1248815224053517, 
  -3.8461538461538463, 
  1.2562282115554928, 
  -3.269048496008233, 
  3.140570528888732, 
  1.0, 
  2.5, 
  2.1528407789851314, 
  5.382101947462829, 
  -3.8461538461538463, 
  -2.355129310456592, 
  -8.280156842250506, 
  -5.887823276141479, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_x0_coefl[8] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_x0_coefr[8] = { 
  0.5, 
  -0.34641016151377546, 
  0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ho_ph_x0_cmap[8] = { 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ho_ph_x0_vlmap[8] = { 
  0, 
  1, 
  2, 
  3, 
  4, 
  5, 
  6, 
  7, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ho_ph_x0_vrmap[8] = { 
  0, 
  1, 
  2, 
  3, 
  4, 
  5, 
  6, 
  7, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_x0_Vd0[24] = { 
  0.0, 
  2.5, 
  6.25, 
  9.936188210700607, 
  -9.615384615384615, 
  24.840470526751517, 
  0.0, 
  2.5, 
  6.25, 
  3.92285819520988, 
  -9.615384615384615, 
  9.8071454880247, 
  0.0, 
  2.5, 
  6.25, 
  -3.92285819520988, 
  -9.615384615384615, 
  -9.8071454880247, 
  0.0, 
  2.5, 
  6.25, 
  -9.936188210700607, 
  -9.615384615384615, 
  -24.840470526751517, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ho_ph_x0_vld0map[8] = { 
  0, 
  0, 
  1, 
  2, 
  0, 
  3, 
  4, 
  5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p2_ho_ph_x0_vrd0map[8] = { 
  0, 
  0, 
  1, 
  2, 
  0, 
  3, 
  4, 
  5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_x0_dcoefl0[8] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  -0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p2_ho_ph_x0_dcoefr0[8] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  0.20139513400278566, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x1v_ser_p2_ho_ph_x0_meta = { .nb = 8, .na = 1, .nmv = 8, .no = 1, .ni = 4, .ns = 3, .fhat_off = 0 }; 
