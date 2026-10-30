// Mon Jul 19 01:21:18 2021
#pragma once
#include <gkyl_util.h>
#include <gkyl_mat.h>
EXTERN_C_BEG
GKYL_CU_DH void binop_div_set_1d_ser_p0(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);
GKYL_CU_DH void binop_div_set_1d_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);
GKYL_CU_DH void binop_div_set_1d_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);
GKYL_CU_DH void binop_div_set_1d_ser_p3(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);
GKYL_CU_DH void binop_div_set_2d_ser_p0(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);
GKYL_CU_DH void binop_div_set_2d_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);
GKYL_CU_DH void binop_div_set_2d_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);
GKYL_CU_DH void binop_div_set_2d_ser_p3(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);
GKYL_CU_DH void binop_div_set_3d_ser_p0(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);
GKYL_CU_DH void binop_div_set_3d_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);
GKYL_CU_DH void binop_div_set_3d_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);
GKYL_CU_DH void binop_div_set_3d_ser_p3(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);

GKYL_CU_DH void binop_div_set_2d_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);
GKYL_CU_DH void binop_div_set_2d_tensor_p3(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);
GKYL_CU_DH void binop_div_set_3d_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g);

GKYL_CU_DH void binop_div_copy_sol(const struct gkyl_mat *rhs, double* GKYL_RESTRICT fdivg);
EXTERN_C_END
