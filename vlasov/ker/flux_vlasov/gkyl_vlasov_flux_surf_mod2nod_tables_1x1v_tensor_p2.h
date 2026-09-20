#pragma once 
// Generated nodal surface-evaluation (MODAL->NODAL) tables for the 
// alpha_quad producers and lax_flux_nodal kernels (see 
// surf-eval-tables.mac). Header-static so both host and device 
// instantiations of GKYL_CU_DH kernels can read them; unused arrays 
// are elided per TU. 
#include <gkyl_util.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_meta.h> 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p2_conf_ev[12] = { 
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
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p2_conf_dx0[12] = { 
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
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p2_vel_ev_v0[3] = { 
  0.7071067811865475, 
  -1.224744871391589, 
  1.5811388300841895, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p2_vel_sparse_idx[3] = { 
  0, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_Cm[12] = { 
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
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_CmD[12] = { 
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
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_CmDx0[12] = { 
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
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_V[3] = { 
  1.0, 
  2.5, 
  -3.8461538461538463, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_coefl[9] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  -0.20139513400278566, 
  0.169, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_coefr[9] = { 
  0.5, 
  0.34641016151377546, 
  -0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  0.20139513400278566, 
  -0.20139513400278566, 
  0.169, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p2_ph_v0_cmap[9] = { 
  0, 
  1, 
  0, 
  1, 
  2, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p2_ph_v0_vlmap[9] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p2_ph_v0_vrmap[9] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_Vd0[3] = { 
  0.0, 
  2.5, 
  -11.538461538461538, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p2_ph_v0_vld0map[9] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p2_ph_v0_vrd0map[9] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_dcoefl0[9] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  -0.20139513400278566, 
  0.169, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_dcoefr0[9] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  0.2906888370749726, 
  -0.20139513400278566, 
  0.20139513400278566, 
  -0.169, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x1v_tensor_p2_ph_v0_meta = { .nb = 9, .na = 3, .nmv = 3, .no = 4, .ni = 1, .ns = 3, .fhat_off = 0 }; 
