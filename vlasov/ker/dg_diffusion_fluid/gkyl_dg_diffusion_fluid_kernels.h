#pragma once
#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_1x_ser_p1_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_1x_ser_p1_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_1x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_1x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_1x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_1x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_ser_p1_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_ser_p1_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_ser_p1_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_2x_ser_p1_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_2x_ser_p1_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_2x_ser_p1_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_2x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_2x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_2x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_2x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfy_2x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfy_2x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfy_2x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfy_2x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_ser_p1_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_ser_p1_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_ser_p1_constcoeff_diffdirsxyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_ser_p1_constcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_ser_p1_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_ser_p1_constcoeff_diffdirsyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_ser_p1_constcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_ser_p1_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_ser_p1_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_ser_p1_constcoeff_diffdirsxyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_ser_p1_constcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_ser_p1_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_ser_p1_constcoeff_diffdirsyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_ser_p1_constcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfy_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfy_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfy_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfy_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfz_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfz_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfz_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfz_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_1x_ser_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_1x_ser_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_1x_ser_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_1x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_1x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_1x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_1x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfx_1x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfx_1x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_ser_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_ser_p2_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_ser_p2_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_2x_ser_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_2x_ser_p2_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_2x_ser_p2_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_2x_ser_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_2x_ser_p2_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_2x_ser_p2_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfx_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfx_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfy_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfy_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfy_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfy_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfy_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfy_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_ser_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_ser_p2_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_ser_p2_constcoeff_diffdirsxyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_ser_p2_constcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_ser_p2_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_ser_p2_constcoeff_diffdirsyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_ser_p2_constcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_ser_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_ser_p2_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_ser_p2_constcoeff_diffdirsxyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_ser_p2_constcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_ser_p2_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_ser_p2_constcoeff_diffdirsyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_ser_p2_constcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_3x_ser_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_3x_ser_p2_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_3x_ser_p2_constcoeff_diffdirsxyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_3x_ser_p2_constcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_3x_ser_p2_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_3x_ser_p2_constcoeff_diffdirsyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_3x_ser_p2_constcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfx_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfx_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfy_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfy_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfy_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfy_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfy_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfy_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfz_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfz_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfz_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfz_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfz_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfz_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_1x_ser_p3_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_1x_ser_p3_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_1x_ser_p3_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_1x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_1x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_1x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_1x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfx_1x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfx_1x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_ser_p3_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_ser_p3_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_ser_p3_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_2x_ser_p3_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_2x_ser_p3_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_2x_ser_p3_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_2x_ser_p3_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_2x_ser_p3_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_2x_ser_p3_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfx_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfx_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfy_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfy_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfy_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfy_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfy_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfy_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_1x_tensor_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_1x_tensor_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_1x_tensor_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_1x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_1x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_1x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_1x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfx_1x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfx_1x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p2_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p2_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_2x_tensor_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_2x_tensor_p2_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_2x_tensor_p2_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_2x_tensor_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_2x_tensor_p2_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_2x_tensor_p2_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfx_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfx_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfy_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfy_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfy_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfy_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfy_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfy_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_tensor_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_tensor_p2_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_tensor_p2_constcoeff_diffdirsxyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_tensor_p2_constcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_tensor_p2_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_tensor_p2_constcoeff_diffdirsyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_3x_tensor_p2_constcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_tensor_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_tensor_p2_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_tensor_p2_constcoeff_diffdirsxyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_tensor_p2_constcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_tensor_p2_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_tensor_p2_constcoeff_diffdirsyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_3x_tensor_p2_constcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_3x_tensor_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_3x_tensor_p2_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_3x_tensor_p2_constcoeff_diffdirsxyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_3x_tensor_p2_constcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_3x_tensor_p2_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_3x_tensor_p2_constcoeff_diffdirsyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_3x_tensor_p2_constcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfx_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfx_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfy_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfy_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfy_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfy_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfy_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfy_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfz_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfz_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfz_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfz_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfz_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfz_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p3_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p3_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_vol_2x_tensor_p3_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_2x_tensor_p3_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_2x_tensor_p3_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_vol_2x_tensor_p3_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_2x_tensor_p3_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_2x_tensor_p3_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_vol_2x_tensor_p3_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_2x_tensor_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_2x_tensor_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfx_2x_tensor_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_2x_tensor_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfx_2x_tensor_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfx_2x_tensor_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_fluid_order2_surfy_2x_tensor_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfy_2x_tensor_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_surfy_2x_tensor_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfy_2x_tensor_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_surfy_2x_tensor_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfy_2x_tensor_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);


EXTERN_C_END
