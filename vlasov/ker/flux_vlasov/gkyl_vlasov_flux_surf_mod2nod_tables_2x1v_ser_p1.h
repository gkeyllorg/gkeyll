#pragma once 
// Generated nodal surface-evaluation (MODAL->NODAL) tables for the 
// alpha_quad producers and lax_flux_nodal kernels (see 
// surf-eval-tables.mac). Header-static so both host and device 
// instantiations of GKYL_CU_DH kernels can read them; unused arrays 
// are elided per TU. 
#include <gkyl_util.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_meta.h> 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_conf_ev[16] = { 
  0.5, 
  -0.5, 
  -0.5, 
  0.5, 
  0.5, 
  -0.5, 
  0.5, 
  -0.5, 
  0.5, 
  0.5, 
  -0.5, 
  -0.5, 
  0.5, 
  0.5, 
  0.5, 
  0.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_conf_dx0[16] = { 
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
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_vel_ev_v0[2] = { 
  0.7071067811865475, 
  -1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_vel_sparse_idx[2] = { 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_Cm[16] = { 
  1.0, 
  -1.4433756729740645, 
  -1.4433756729740645, 
  2.0833333333333335, 
  1.0, 
  -1.4433756729740645, 
  1.4433756729740645, 
  -2.0833333333333335, 
  1.0, 
  1.4433756729740645, 
  -1.4433756729740645, 
  -2.0833333333333335, 
  1.0, 
  1.4433756729740645, 
  1.4433756729740645, 
  2.0833333333333335, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_CmD[16] = { 
  0.0, 
  2.5, 
  0.0, 
  -3.6084391824351614, 
  0.0, 
  2.5, 
  0.0, 
  3.6084391824351614, 
  0.0, 
  2.5, 
  0.0, 
  -3.6084391824351614, 
  0.0, 
  2.5, 
  0.0, 
  3.6084391824351614, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_CmDx0[16] = { 
  0.0, 
  2.5, 
  0.0, 
  -3.6084391824351614, 
  0.0, 
  2.5, 
  0.0, 
  3.6084391824351614, 
  0.0, 
  2.5, 
  0.0, 
  -3.6084391824351614, 
  0.0, 
  2.5, 
  0.0, 
  3.6084391824351614, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_CmDx1[16] = { 
  0.0, 
  0.0, 
  2.5, 
  -3.6084391824351614, 
  0.0, 
  0.0, 
  2.5, 
  -3.6084391824351614, 
  0.0, 
  0.0, 
  2.5, 
  3.6084391824351614, 
  0.0, 
  0.0, 
  2.5, 
  3.6084391824351614, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_V[2] = { 
  1.0, 
  2.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_coefl[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_coefr[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  -0.2449489742783178, 
  0.16970562748477142, 
  -0.16970562748477142, 
  -0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_v0_cmap[8] = { 
  0, 
  1, 
  2, 
  0, 
  3, 
  1, 
  2, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_v0_vlmap[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  1, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_v0_vrmap[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  1, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_Vd0[2] = { 
  0.0, 
  2.5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_v0_vld0map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  1, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_v0_vrd0map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  1, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_dcoefl0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_dcoefr0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_2x1v_ser_p1_ph_v0_meta = { .nb = 8, .na = 4, .nmv = 2, .no = 4, .ni = 1, .ns = 2, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_vel_vol_dv0[4] = { 
  0.0, 
  1.224744871391589, 
  0.0, 
  1.224744871391589, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_confsurf_x0_ev_l[8] = { 
  0.5, 
  0.8660254037844386, 
  -0.5, 
  -0.8660254037844386, 
  0.5, 
  0.8660254037844386, 
  0.5, 
  0.8660254037844386, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_confsurf_x0_ev_r[8] = { 
  0.5, 
  -0.8660254037844386, 
  -0.5, 
  0.8660254037844386, 
  0.5, 
  -0.8660254037844386, 
  0.5, 
  -0.8660254037844386, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x0_Cm[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x0_V[8] = { 
  1.0, 
  2.5, 
  -1.4433756729740645, 
  -3.6084391824351614, 
  1.0, 
  2.5, 
  1.4433756729740645, 
  3.6084391824351614, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x0_coefl[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x0_coefr[8] = { 
  0.3535533905932737, 
  -0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  -0.16970562748477142, 
  -0.16970562748477142, 
  0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x0_cmap[8] = { 
  0, 
  0, 
  1, 
  0, 
  1, 
  0, 
  1, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x0_vlmap[8] = { 
  0, 
  1, 
  0, 
  2, 
  1, 
  3, 
  2, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x0_vrmap[8] = { 
  0, 
  1, 
  0, 
  2, 
  1, 
  3, 
  2, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x0_Vd0[6] = { 
  0.0, 
  2.5, 
  6.25, 
  0.0, 
  2.5, 
  6.25, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x0_vld0map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x0_vrd0map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x0_dcoefl0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x0_dcoefr0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  -0.16970562748477142, 
  0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_2x1v_ser_p1_ph_x0_meta = { .nb = 8, .na = 2, .nmv = 4, .no = 2, .ni = 2, .ns = 2, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_confsurf_x1_ev_l[8] = { 
  0.5, 
  -0.5, 
  0.8660254037844386, 
  -0.8660254037844386, 
  0.5, 
  0.5, 
  0.8660254037844386, 
  0.8660254037844386, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_confsurf_x1_ev_r[8] = { 
  0.5, 
  -0.5, 
  -0.8660254037844386, 
  0.8660254037844386, 
  0.5, 
  0.5, 
  -0.8660254037844386, 
  -0.8660254037844386, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x1_Cm[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x1_V[8] = { 
  1.0, 
  2.5, 
  -1.4433756729740645, 
  -3.6084391824351614, 
  1.0, 
  2.5, 
  1.4433756729740645, 
  3.6084391824351614, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x1_coefl[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x1_coefr[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  -0.2449489742783178, 
  0.2449489742783178, 
  -0.16970562748477142, 
  0.16970562748477142, 
  -0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x1_cmap[8] = { 
  0, 
  1, 
  0, 
  0, 
  1, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x1_vlmap[8] = { 
  0, 
  0, 
  1, 
  2, 
  1, 
  2, 
  3, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x1_vrmap[8] = { 
  0, 
  0, 
  1, 
  2, 
  1, 
  2, 
  3, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x1_Vd0[6] = { 
  0.0, 
  2.5, 
  6.25, 
  0.0, 
  2.5, 
  6.25, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x1_vld0map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x1_vrd0map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x1_dcoefl0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x1_dcoefr0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  -0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_2x1v_ser_p1_ph_x1_meta = { .nb = 8, .na = 2, .nmv = 4, .no = 2, .ni = 2, .ns = 2, .fhat_off = 4 }; 
