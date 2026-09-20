#pragma once 
// Generated nodal surface-evaluation (MODAL->NODAL) tables for the 
// alpha_quad producers and lax_flux_nodal kernels (see 
// surf-eval-tables.mac). Header-static so both host and device 
// instantiations of GKYL_CU_DH kernels can read them; unused arrays 
// are elided per TU. 
#include <gkyl_util.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_meta.h> 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_conf_ev[4] = { 
  0.7071067811865475, 
  -0.7071067811865475, 
  0.7071067811865475, 
  0.7071067811865475, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_conf_dx0[4] = { 
  0.0, 
  1.224744871391589, 
  0.0, 
  1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_vel_ev_v0[8] = { 
  0.5, 
  -0.8660254037844386, 
  -0.5, 
  0.8660254037844386, 
  0.5, 
  -0.8660254037844386, 
  0.5, 
  -0.8660254037844386, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_vel_nodes_v0[2] = { 
  -0.5773502691896258, 
  0.5773502691896258, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_vel_dv1_v0[8] = { 
  0.0, 
  0.0, 
  0.8660254037844386, 
  -1.5, 
  0.0, 
  0.0, 
  0.8660254037844386, 
  -1.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_vel_ev_v1[8] = { 
  0.5, 
  -0.5, 
  -0.8660254037844386, 
  0.8660254037844386, 
  0.5, 
  0.5, 
  -0.8660254037844386, 
  -0.8660254037844386, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_vel_nodes_v1[2] = { 
  -0.5773502691896258, 
  0.5773502691896258, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_vel_dv0_v1[8] = { 
  0.0, 
  0.8660254037844386, 
  0.0, 
  -1.5, 
  0.0, 
  0.8660254037844386, 
  0.0, 
  -1.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_vel_sparse_idx[3] = { 
  0, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v0_Cm[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v0_CmD[4] = { 
  0.0, 
  2.5, 
  0.0, 
  2.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v0_CmDx0[4] = { 
  0.0, 
  2.5, 
  0.0, 
  2.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v0_V[8] = { 
  1.0, 
  2.5, 
  -1.4433756729740645, 
  -3.6084391824351614, 
  1.0, 
  2.5, 
  1.4433756729740645, 
  3.6084391824351614, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v0_coefl[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v0_coefr[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  -0.2449489742783178, 
  0.2449489742783178, 
  -0.16970562748477142, 
  0.16970562748477142, 
  -0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_v0_cmap[8] = { 
  0, 
  1, 
  0, 
  0, 
  1, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_v0_vlmap[8] = { 
  0, 
  0, 
  1, 
  2, 
  1, 
  2, 
  3, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_v0_vrmap[8] = { 
  0, 
  0, 
  1, 
  2, 
  1, 
  2, 
  3, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v0_Vd0[6] = { 
  0.0, 
  2.5, 
  -3.6084391824351614, 
  0.0, 
  2.5, 
  3.6084391824351614, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_v0_vld0map[8] = { 
  0, 
  0, 
  1, 
  0, 
  1, 
  0, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_v0_vrd0map[8] = { 
  0, 
  0, 
  1, 
  0, 
  1, 
  0, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v0_dcoefl0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v0_dcoefr0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v0_Vd1[6] = { 
  0.0, 
  2.5, 
  6.25, 
  0.0, 
  2.5, 
  6.25, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_v0_vld1map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_v0_vrd1map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v0_dcoefl1[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v0_dcoefr1[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  -0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x2v_ser_p1_ph_v0_meta = { .nb = 8, .na = 2, .nmv = 4, .no = 2, .ni = 2, .ns = 3, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v1_Cm[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v1_CmDx0[4] = { 
  0.0, 
  2.5, 
  0.0, 
  2.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v1_V[8] = { 
  1.0, 
  -1.4433756729740645, 
  2.5, 
  -3.6084391824351614, 
  1.0, 
  1.4433756729740645, 
  2.5, 
  3.6084391824351614, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v1_coefl[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v1_coefr[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  -0.2449489742783178, 
  0.16970562748477142, 
  -0.16970562748477142, 
  -0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_v1_cmap[8] = { 
  0, 
  1, 
  0, 
  0, 
  1, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_v1_vlmap[8] = { 
  0, 
  0, 
  1, 
  2, 
  1, 
  2, 
  3, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_v1_vrmap[8] = { 
  0, 
  0, 
  1, 
  2, 
  1, 
  2, 
  3, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v1_Vd0[6] = { 
  0.0, 
  2.5, 
  6.25, 
  0.0, 
  2.5, 
  6.25, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_v1_vld0map[8] = { 
  0, 
  0, 
  1, 
  0, 
  1, 
  0, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_v1_vrd0map[8] = { 
  0, 
  0, 
  1, 
  0, 
  1, 
  0, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v1_dcoefl0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v1_dcoefr0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  -0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v1_Vd1[6] = { 
  0.0, 
  2.5, 
  -3.6084391824351614, 
  0.0, 
  2.5, 
  3.6084391824351614, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_v1_vld1map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_v1_vrd1map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v1_dcoefl1[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_v1_dcoefr1[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x2v_ser_p1_ph_v1_meta = { .nb = 8, .na = 2, .nmv = 4, .no = 2, .ni = 2, .ns = 3, .fhat_off = 4 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_vel_vol_dv0[16] = { 
  0.0, 
  0.8660254037844386, 
  0.0, 
  -0.8660254037844386, 
  0.0, 
  0.8660254037844386, 
  0.0, 
  0.8660254037844386, 
  0.0, 
  0.8660254037844386, 
  0.0, 
  -0.8660254037844386, 
  0.0, 
  0.8660254037844386, 
  0.0, 
  0.8660254037844386, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_vel_vol_dv1[16] = { 
  0.0, 
  0.0, 
  0.8660254037844386, 
  -0.8660254037844386, 
  0.0, 
  0.0, 
  0.8660254037844386, 
  -0.8660254037844386, 
  0.0, 
  0.0, 
  0.8660254037844386, 
  0.8660254037844386, 
  0.0, 
  0.0, 
  0.8660254037844386, 
  0.8660254037844386, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_confsurf_x0_ev_l[2] = { 
  0.7071067811865475, 
  1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_confsurf_x0_ev_r[2] = { 
  0.7071067811865475, 
  -1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_x0_Cm[1] = { 
  1.0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_x0_V[32] = { 
  1.0, 
  2.5, 
  -1.4433756729740645, 
  -1.4433756729740645, 
  -3.6084391824351614, 
  -3.6084391824351614, 
  2.0833333333333335, 
  5.208333333333333, 
  1.0, 
  2.5, 
  -1.4433756729740645, 
  1.4433756729740645, 
  -3.6084391824351614, 
  3.6084391824351614, 
  -2.0833333333333335, 
  -5.208333333333333, 
  1.0, 
  2.5, 
  1.4433756729740645, 
  -1.4433756729740645, 
  3.6084391824351614, 
  -3.6084391824351614, 
  -2.0833333333333335, 
  -5.208333333333333, 
  1.0, 
  2.5, 
  1.4433756729740645, 
  1.4433756729740645, 
  3.6084391824351614, 
  3.6084391824351614, 
  2.0833333333333335, 
  5.208333333333333, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_x0_coefl[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_x0_coefr[8] = { 
  0.3535533905932737, 
  -0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  -0.16970562748477142, 
  -0.16970562748477142, 
  0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_x0_cmap[8] = { 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_x0_vlmap[8] = { 
  0, 
  1, 
  2, 
  3, 
  4, 
  5, 
  6, 
  7, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_x0_vrmap[8] = { 
  0, 
  1, 
  2, 
  3, 
  4, 
  5, 
  6, 
  7, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_x0_Vd0[20] = { 
  0.0, 
  2.5, 
  6.25, 
  -3.6084391824351614, 
  -9.021097956087903, 
  0.0, 
  2.5, 
  6.25, 
  3.6084391824351614, 
  9.021097956087903, 
  0.0, 
  2.5, 
  6.25, 
  -3.6084391824351614, 
  -9.021097956087903, 
  0.0, 
  2.5, 
  6.25, 
  3.6084391824351614, 
  9.021097956087903, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_x0_vld0map[8] = { 
  0, 
  0, 
  1, 
  0, 
  2, 
  0, 
  3, 
  4, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_x0_vrd0map[8] = { 
  0, 
  0, 
  1, 
  0, 
  2, 
  0, 
  3, 
  4, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_x0_dcoefl0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_x0_dcoefr0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  -0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_x0_Vd1[20] = { 
  0.0, 
  2.5, 
  6.25, 
  -3.6084391824351614, 
  -9.021097956087903, 
  0.0, 
  2.5, 
  6.25, 
  -3.6084391824351614, 
  -9.021097956087903, 
  0.0, 
  2.5, 
  6.25, 
  3.6084391824351614, 
  9.021097956087903, 
  0.0, 
  2.5, 
  6.25, 
  3.6084391824351614, 
  9.021097956087903, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_x0_vld1map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  2, 
  3, 
  4, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_ph_x0_vrd1map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  2, 
  3, 
  4, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_x0_dcoefl1[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_ph_x0_dcoefr1[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  -0.16970562748477142, 
  0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x2v_ser_p1_ph_x0_meta = { .nb = 8, .na = 1, .nmv = 8, .no = 1, .ni = 4, .ns = 3, .fhat_off = 0 }; 
