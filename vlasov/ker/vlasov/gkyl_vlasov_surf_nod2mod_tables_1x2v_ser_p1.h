#pragma once 
// Generated projection (NODAL->MODAL) tables for the from-flux 
// vlasov_(boundary_)(ho_)surf* consumer kernels (see 
// surf-eval-tables.mac). Header-static so both host and device 
// instantiations of GKYL_CU_DH kernels can read them; unused arrays 
// are elided per TU. 
#include <gkyl_util.h> 
#include <gkyl_vlasov_surf_nod2mod_meta.h> 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_v0_Cw[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_v0_Vw[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_v0_kamap[4] = { 
  0, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_v0_kbmap[4] = { 
  0, 
  0, 
  1, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_v0_out_off[5] = { 
  0, 
  2, 
  4, 
  6, 
  8, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_v0_out_mode[8] = { 
  0, 
  2, 
  1, 
  4, 
  3, 
  6, 
  5, 
  7, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_v0_out_cl[8] = { 
  0.3535533905932737, 
  -0.6123724356957944, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.16970562748477142, 
  -0.29393876913398137, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_v0_out_cr[8] = { 
  -0.3535533905932737, 
  -0.6123724356957944, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.16970562748477142, 
  -0.29393876913398137, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_v0_pm_a[8] = { 
  0, 
  1, 
  0, 
  0, 
  1, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_v0_pm_b[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  1, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_v0_pm_cl[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  -0.6123724356957944, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.16970562748477142, 
  -0.42426406871192845, 
  -0.29393876913398137, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_v0_pm_cr[8] = { 
  -0.3535533905932737, 
  -0.2449489742783178, 
  -0.6123724356957944, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.16970562748477142, 
  -0.42426406871192845, 
  -0.29393876913398137, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_1x2v_ser_p1_prj_v0_meta = { .nk = 4, .np = 8, .na = 2, .nmv = 2, .no = 2, .ni = 2, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_v1_Cw[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_v1_Vw[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_v1_kamap[4] = { 
  0, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_v1_kbmap[4] = { 
  0, 
  0, 
  1, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_v1_out_off[5] = { 
  0, 
  2, 
  4, 
  6, 
  8, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_v1_out_mode[8] = { 
  0, 
  3, 
  1, 
  5, 
  2, 
  6, 
  4, 
  7, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_v1_out_cl[8] = { 
  0.3535533905932737, 
  -0.6123724356957944, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.16970562748477142, 
  -0.29393876913398137, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_v1_out_cr[8] = { 
  -0.3535533905932737, 
  -0.6123724356957944, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.16970562748477142, 
  -0.29393876913398137, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_v1_pm_a[8] = { 
  0, 
  1, 
  0, 
  0, 
  1, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_v1_pm_b[8] = { 
  0, 
  0, 
  1, 
  0, 
  1, 
  0, 
  1, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_v1_pm_cl[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  -0.6123724356957944, 
  0.16970562748477142, 
  -0.42426406871192845, 
  -0.42426406871192845, 
  -0.29393876913398137, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_v1_pm_cr[8] = { 
  -0.3535533905932737, 
  -0.2449489742783178, 
  -0.2449489742783178, 
  -0.6123724356957944, 
  -0.16970562748477142, 
  -0.42426406871192845, 
  -0.42426406871192845, 
  -0.29393876913398137, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_1x2v_ser_p1_prj_v1_meta = { .nk = 4, .np = 8, .na = 2, .nmv = 2, .no = 2, .ni = 2, .fhat_off = 4 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_x0_Cw[1] = { 
  1.0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_x0_Vw[16] = { 
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
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_x0_kamap[4] = { 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_x0_kbmap[4] = { 
  0, 
  1, 
  2, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_x0_out_off[5] = { 
  0, 
  2, 
  4, 
  6, 
  8, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_x0_out_mode[8] = { 
  0, 
  1, 
  2, 
  4, 
  3, 
  5, 
  6, 
  7, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_x0_out_cl[8] = { 
  0.3535533905932737, 
  -0.6123724356957944, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.16970562748477142, 
  -0.29393876913398137, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_x0_out_cr[8] = { 
  -0.3535533905932737, 
  -0.6123724356957944, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.16970562748477142, 
  -0.29393876913398137, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_x0_pm_a[8] = { 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x2v_ser_p1_prj_x0_pm_b[8] = { 
  0, 
  0, 
  1, 
  2, 
  1, 
  2, 
  3, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_x0_pm_cl[8] = { 
  0.3535533905932737, 
  -0.6123724356957944, 
  0.2449489742783178, 
  0.2449489742783178, 
  -0.42426406871192845, 
  -0.42426406871192845, 
  0.16970562748477142, 
  -0.29393876913398137, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x2v_ser_p1_prj_x0_pm_cr[8] = { 
  -0.3535533905932737, 
  -0.6123724356957944, 
  -0.2449489742783178, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.42426406871192845, 
  -0.16970562748477142, 
  -0.29393876913398137, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_1x2v_ser_p1_prj_x0_meta = { .nk = 4, .np = 8, .na = 1, .nmv = 4, .no = 1, .ni = 4, .fhat_off = 0 }; 
