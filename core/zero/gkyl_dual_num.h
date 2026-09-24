#pragma once

#include <math.h>
#include <gkyl_util.h>

// Dual or hyperreal number
struct gkyl_dn { double x[2]; };
// Dual or hyperreal 2D number
struct gkyl_dn2 { double x[3]; };

/**
   Use the new0 constructor to create a real number and new1 to create
   a dual-number initialized for taking derivatives.

   Basic operators: neg, add, sadd, sub, ssub, mul, smul, div, sdiv
   
   Functions: sq, cube, npow, sqrt, sin, cos, tan, log
*/

// Construct new dual numbers
GKYL_CU_DH
static inline struct gkyl_dn
gdn_new(double x0, double x1)
{
  return (struct gkyl_dn) { x0, x1 };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_new0(double x0)
{
  return gdn_new(x0, 0.0);
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_new1(double x0)
{
  return gdn_new(x0, 1.0);
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_neg(struct gkyl_dn d1)
{
  return (struct gkyl_dn) { -d1.x[0], -d1.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_add(struct gkyl_dn d1, struct gkyl_dn d2)
{
  return (struct gkyl_dn) { d1.x[0]+d2.x[0], d1.x[1]+d2.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_sadd(double s, struct gkyl_dn d1)
{
  return (struct gkyl_dn) { s+d1.x[0], d1.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_sub(struct gkyl_dn d1, struct gkyl_dn d2)
{
  return (struct gkyl_dn) { d1.x[0]-d2.x[0], d1.x[1]-d2.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_ssub(double s, struct gkyl_dn d1)
{
  return (struct gkyl_dn) { s-d1.x[0], -d1.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_mul(struct gkyl_dn d1, struct gkyl_dn d2)
{
  return (struct gkyl_dn) {
      d1.x[0]*d2.x[0],
      d1.x[0]*d2.x[1] + d1.x[1]*d2.x[0]    
  };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_smul(double s, struct gkyl_dn d1)
{
  return (struct gkyl_dn) { s*d1.x[0], s*d1.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_inv(struct gkyl_dn d1)
{
  return (struct gkyl_dn) {
    1.0/d1.x[0],
    -d1.x[1]/(d1.x[0]*d1.x[0])
  };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_div(struct gkyl_dn d1, struct gkyl_dn d2)
{
  return gdn_mul(d1, gdn_inv(d2));
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_sdiv(double s, struct gkyl_dn d1)
{
  return gdn_smul(s, gdn_inv(d1));
}


GKYL_CU_DH
static inline struct gkyl_dn
gdn_sq(struct gkyl_dn d1)
{
  return (struct gkyl_dn) { d1.x[0]*d1.x[0], 2.0*d1.x[0]*d1.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_cube(struct gkyl_dn d1)
{
  return gdn_mul(d1, gdn_sq(d1));
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_npow(struct gkyl_dn d1, int n)
{
  double pn1 = pow(d1.x[0], n-1);
  return (struct gkyl_dn) { pn1*d1.x[0], n*pn1*d1.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_sqrt(struct gkyl_dn d1)
{
  double sx = sqrt(d1.x[0]);
  return (struct gkyl_dn) { sx, d1.x[1]*0.5/sx  };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_sin(struct gkyl_dn d1)
{
  return (struct gkyl_dn) { sin(d1.x[0]), d1.x[1]*cos(d1.x[0]) };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_cos(struct gkyl_dn d1)
{
  return (struct gkyl_dn) { cos(d1.x[0]), -d1.x[1]*sin(d1.x[0]) };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_tan(struct gkyl_dn d1)
{
  return gdn_div(gdn_sin(d1), gdn_cos(d1));
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_log(struct gkyl_dn d1)
{
  return (struct gkyl_dn) { log(d1.x[0]), d1.x[1]/d1.x[0] };
}

/***************/
/** G[2] Duals */
/***************/

// Construct new dual numbers
GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_new(double x0, double x1, double x2)
{
  return (struct gkyl_dn2) { x0, x1, x2 };
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_new00(double x0)
{
  return gdn2_new(x0, 0.0, 0.0);
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_new01(double x0)
{
  return gdn2_new(x0, 0.0, 1.0);
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_new10(double x0)
{
  return gdn2_new(x0, 1.0, 0.0);
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_new11(double x0)
{
  return (struct gkyl_dn2) { x0, 1.0, 1.0 };
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_neg(struct gkyl_dn2 d1)
{
  return (struct gkyl_dn2) { -d1.x[0], -d1.x[1], -d1.x[2] };
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_add(struct gkyl_dn2 d1, struct gkyl_dn2 d2)
{
  return (struct gkyl_dn2) { d1.x[0]+d2.x[0], d1.x[1]+d2.x[1], d1.x[2]+d2.x[2] };
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_sadd(double s, struct gkyl_dn2 d1)
{
  return (struct gkyl_dn2) { s+d1.x[0], d1.x[1], d1.x[2] };
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_sub(struct gkyl_dn2 d1, struct gkyl_dn2 d2)
{
  return (struct gkyl_dn2) { d1.x[0]-d2.x[0], d1.x[1]-d2.x[1], d1.x[2]-d2.x[2] };
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_ssub(double s, struct gkyl_dn2 d1)
{
  return (struct gkyl_dn2) { s-d1.x[0], -d1.x[1], -d1.x[2] };
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_mul(struct gkyl_dn2 d1, struct gkyl_dn2 d2)
{
  return (struct gkyl_dn2) {
    d1.x[0]*d2.x[0],
    d1.x[0]*d2.x[1] + d1.x[1]*d2.x[0], 
    d1.x[0]*d2.x[2] + d1.x[2]*d2.x[0],
  };
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_smul(double s, struct gkyl_dn2 d1)
{
  return (struct gkyl_dn2) { s*d1.x[0], s*d1.x[1], s*d1.x[2] };
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_inv(struct gkyl_dn2 d1)
{
  return (struct gkyl_dn2) {
    1.0/d1.x[0], -d1.x[1]/(d1.x[0]*d1.x[0]), -d1.x[2]/(d1.x[0]*d1.x[0]),
  };
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_div(struct gkyl_dn2 d1, struct gkyl_dn2 d2)
{
  return gdn2_mul(d1, gdn2_inv(d2));
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_sdiv(double s, struct gkyl_dn2 d1)
{
  return gdn2_smul(s, gdn2_inv(d1));
}


GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_sq(struct gkyl_dn2 d1)
{
  return gdn2_mul(d1, d1);
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_cube(struct gkyl_dn2 d1)
{
  return gdn2_mul(d1, gdn2_sq(d1));
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_npow(struct gkyl_dn2 d1, int n)
{
  double pn1 = pow(d1.x[0], n-1);
  return (struct gkyl_dn2) { pn1*d1.x[0], n*pn1*d1.x[1], n*pn1*d1.x[2] };
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_sqrt(struct gkyl_dn2 d1)
{
  double sx = sqrt(d1.x[0]);
  return (struct gkyl_dn2) { sx, d1.x[1]*0.5/sx, d1.x[2]*0.5/sx  };
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_sin(struct gkyl_dn2 d1)
{
  return (struct gkyl_dn2) { sin(d1.x[0]), d1.x[1]*cos(d1.x[0]), d1.x[2]*cos(d1.x[0]) };
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_cos(struct gkyl_dn2 d1)
{
  return (struct gkyl_dn2) { cos(d1.x[0]), -d1.x[1]*sin(d1.x[0]), -d1.x[2]*sin(d1.x[0]) };
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_tan(struct gkyl_dn2 d1)
{
  return gdn2_div(gdn2_sin(d1), gdn2_cos(d1));
}

GKYL_CU_DH
static inline struct gkyl_dn2
gdn2_log(struct gkyl_dn2 d1)
{
  return (struct gkyl_dn2) { log(d1.x[0]), d1.x[1]/d1.x[0], d1.x[2]/d1.x[0] };
}
