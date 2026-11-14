// -- Gkyl ---------------------------------------------------------------------
//
// C header for TwistShift kernels.
//
//    _______     ___
// + 6 @ |||| # P ||| +
// -----------------------------------------------------------------------------
#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
EXTERN_C_BEG 
  GKYL_CU_DH void twistshift_xlimdg_2x_ser_p1_yshift_p1(double sFac, const double *xLimLo, const double *xLimUp, double yLimLo, double yLimUp, double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);
  GKYL_CU_DH void twistshift_ylimdg_2x_ser_p1_yshift_p1(double sFac, double xLimLo, double xLimUp, const double *yLimLo, const double *yLimUp, double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);
  GKYL_CU_DH void twistshift_fullcell_2x_ser_p1_yshift_p1(double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);

  GKYL_CU_DH void twistshift_xlimdg_2x_ser_p1_yshift_p2(double sFac, const double *xLimLo, const double *xLimUp, double yLimLo, double yLimUp, double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);
  GKYL_CU_DH void twistshift_ylimdg_2x_ser_p1_yshift_p2(double sFac, double xLimLo, double xLimUp, const double *yLimLo, const double *yLimUp, double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);
  GKYL_CU_DH void twistshift_fullcell_2x_ser_p1_yshift_p2(double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);


  GKYL_CU_DH void twistshift_xlimdg_3x_ser_p1_yshift_p1(double sFac, const double *xLimLo, const double *xLimUp, double yLimLo, double yLimUp, double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);
  GKYL_CU_DH void twistshift_ylimdg_3x_ser_p1_yshift_p1(double sFac, double xLimLo, double xLimUp, const double *yLimLo, const double *yLimUp, double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);
  GKYL_CU_DH void twistshift_fullcell_3x_ser_p1_yshift_p1(double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);

  GKYL_CU_DH void twistshift_xlimdg_3x_ser_p1_yshift_p2(double sFac, const double *xLimLo, const double *xLimUp, double yLimLo, double yLimUp, double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);
  GKYL_CU_DH void twistshift_ylimdg_3x_ser_p1_yshift_p2(double sFac, double xLimLo, double xLimUp, const double *yLimLo, const double *yLimUp, double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);
  GKYL_CU_DH void twistshift_fullcell_3x_ser_p1_yshift_p2(double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);


  GKYL_CU_DH void twistshift_xlimdg_3x2v_ser_p1_yshift_p1(double sFac, const double *xLimLo, const double *xLimUp, double yLimLo, double yLimUp, double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);
  GKYL_CU_DH void twistshift_ylimdg_3x2v_ser_p1_yshift_p1(double sFac, double xLimLo, double xLimUp, const double *yLimLo, const double *yLimUp, double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);
  GKYL_CU_DH void twistshift_fullcell_3x2v_ser_p1_yshift_p1(double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);

  GKYL_CU_DH void twistshift_xlimdg_3x2v_ser_p1_yshift_p2(double sFac, const double *xLimLo, const double *xLimUp, double yLimLo, double yLimUp, double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);
  GKYL_CU_DH void twistshift_ylimdg_3x2v_ser_p1_yshift_p2(double sFac, double xLimLo, double xLimUp, const double *yLimLo, const double *yLimUp, double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);
  GKYL_CU_DH void twistshift_fullcell_3x2v_ser_p1_yshift_p2(double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);


EXTERN_C_END
