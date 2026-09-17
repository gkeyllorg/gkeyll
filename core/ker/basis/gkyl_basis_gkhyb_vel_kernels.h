// Sat Jun 27 17:10:24 2026
#pragma once
#include <gkyl_util.h>
EXTERN_C_BEG
GKYL_CU_DH void eval_1d_gkhyb_vel_p1(const double *z, double *b);
GKYL_CU_DH double eval_expand_1d_gkhyb_vel_p1(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_1d_gkhyb_vel_p1(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_odd_sign_1d_gkhyb_vel_p1(int dir, const double *f, double *fout );
GKYL_CU_DH void flip_even_sign_1d_gkhyb_vel_p1(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_1d_gkhyb_vel_p1(double *node_coords);
GKYL_CU_DH void nodal_to_modal_1d_gkhyb_vel_p1(const double *fnodal, double *fmodal);
GKYL_CU_DH void quad_to_modal_1d_gkhyb_vel_p1(const double *fquad, double *fmodal, long linc2);
GKYL_CU_DH void modal_to_quad_1d_gkhyb_vel_p1(const double *fmodal, double *fquad, long linc2);
GKYL_CU_DH void eval_2d_gkhyb_vel_p1(const double *z, double *b);
GKYL_CU_DH double eval_expand_2d_gkhyb_vel_p1(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_2d_gkhyb_vel_p1(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_odd_sign_2d_gkhyb_vel_p1(int dir, const double *f, double *fout );
GKYL_CU_DH void flip_even_sign_2d_gkhyb_vel_p1(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_2d_gkhyb_vel_p1(double *node_coords);
GKYL_CU_DH void nodal_to_modal_2d_gkhyb_vel_p1(const double *fnodal, double *fmodal);
GKYL_CU_DH void quad_to_modal_2d_gkhyb_vel_p1(const double *fquad, double *fmodal, long linc2);
GKYL_CU_DH void modal_to_quad_2d_gkhyb_vel_p1(const double *fmodal, double *fquad, long linc2);
EXTERN_C_END
