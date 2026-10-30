#pragma once 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void binop_cross_mul_1d_2d_tensor_p2(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1d_2d_tensor_p2(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1d_2d_tensor_p2(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_1d_2d_tensor_p3(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1d_2d_tensor_p3(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1d_2d_tensor_p3(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_1d_3d_tensor_p2(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1d_3d_tensor_p2(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1d_3d_tensor_p2(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_1d_4d_tensor_p2(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1d_4d_tensor_p2(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1d_4d_tensor_p2(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_2d_3d_tensor_p2(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_2d_3d_tensor_p2(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_2d_3d_tensor_p2(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_2d_3d_tensor_p3(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_2d_3d_tensor_p3(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_2d_3d_tensor_p3(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_2d_4d_tensor_p2(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_2d_4d_tensor_p2(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_2d_4d_tensor_p2(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_2d_5d_tensor_p2(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_2d_5d_tensor_p2(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_2d_5d_tensor_p2(double a, const double *f, const double *g, double *fg, int linc2); 

EXTERN_C_END 
