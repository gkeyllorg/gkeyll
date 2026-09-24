#include <gkyl_tok_geo.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_efit.h>
#include <complex.h>

// Function context to pass to root finder
struct arc_length_ctx {
  const struct gkyl_tok_geo *geo;
  double *arc_memo;
  double *arc_memo_left;
  double *arc_memo_right;
  double psi, rclose, zmin, arcL;
  double rleft, rright, zmax;
  double zmin_left, zmin_right; // for single null full SOL only
  double zmax_left, zmax_right; // for PF UP region (or upper SN, not yet implemented)
  double zmin_iwl, zmax_iwl; // for IWL
  double zmin_iwl_plate, zmax_iwl_plate; // for IWL
  double arcL_right; // this is for when we need to switch sides
  double arcL_left; // this is for when we need to switch sides
  double arcL_tot; // total arc length
  double arcL_start; // For core bloks only. arc length between theta=0 and lower turning point
                     // measured counterclockwise from left side
  double phi_right; // this is for when we need to switch sides
  double phi_left; // this is for when we need to switch sides
  double phi_bot; // For new way of trying to do core
  bool right; // on right side of turning points
  bool pre; // on left side starting from theta=0
  double arcL_q1; // For IWL, arc length in quadrant 1
  double arcL_q2; // For IWL, arc length in quadrant 2
  double arcL_q3; // For IWL, arc length in quadrant 3
  double arcL_q4; // For IWL, arc length in quadrant 4
  bool q1, q2, q3, q4; // For IWL
  double zmaxis;
  enum gkyl_tok_geo_type ftype; // type of geometry
};


// Context to pass to endpoint finder
struct plate_ctx{
  const struct gkyl_tok_geo* geo;
  double psi_curr;
  bool lower;
};

// some helper functions
static double
choose_closest(double ref, double* R, double* out, int nr)
{
  //return fabs(R[0]-ref) < fabs(R[1]-ref) ? out[0] : out[1];
  int imin = 0;
  double min = fabs(R[0]-ref);
  for(int i = 1; i< nr; i++){
    if( fabs(R[i] - ref) < min){
      imin = i;
      min = fabs(R[i] - ref);
    }
  }
  return out[imin];
}

static inline double SQ(double x) { return x*x; }
static inline double CUB(double x) { return x*x*x; }
static inline double complex cSQ(double complex x) { return cpow(x,2.0+0.0*I); }
static double complex cCUB(double complex x) { return cpow(x,3.0+0.0*I); }
static double complex croot(double complex x, double n) { double complex power = 1.0/n + 0.0*I; return cpow(x,power); }


static inline int
get_idx(int dir, double x, const struct gkyl_rect_grid *grid, const struct gkyl_range *range)
{
  double xlower = grid->lower[dir], dx = grid->dx[dir];
  int idx = range->lower[dir] + (int) floor((x-xlower)/dx);
  return idx <= range->upper[dir] ? idx : range->upper[dir];
}

// struct for solutions to roots
struct RdRdZ_sol {
  int nsol;
  double R[2], dRdZ[2];
  double dR[2];
  double dZ[2];
};

// Compute roots R(psi,Z) and dR/dZ(psi,Z) in a p=1 DG cell
static inline struct RdRdZ_sol
calc_RdR_p1(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };

  double y = (Z-xc[1])/(dx[1]*0.5);
  
  double rnorm = (-(1.732050807568877*psi[2]*y)/(3.0*psi[3]*y+1.732050807568877*psi[1]))+(2.0*psi0)/(3.0*psi[3]*y+1.732050807568877*psi[1])-(1.0*psi[0])/(3.0*psi[3]*y+1.732050807568877*psi[1]) ;

  if ((-1<=rnorm) && (rnorm < 1)) {
    double drdznorm = -(3.0*(2.0*psi[3]*psi0-1.0*psi[0]*psi[3]+psi[1]*psi[2]))/SQ(3.0*psi[3]*y+1.732050807568877*psi[1]) ;
    
    sol.nsol = 1;
    sol.R[0] = rnorm*dx[0]*0.5 + xc[0];
    sol.dRdZ[0] = drdznorm*dx[0]/dx[1];
  }
  return sol;
}

// Compute roots R(psi,Z) and dR/dZ(psi,Z) in a p=2 DG cell with tensor basis
static inline struct RdRdZ_sol
calc_RdR_p2_tensor(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };
  double y = (Z-xc[1])/(dx[1]*0.5);

  double aq = 0.125*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4]);
  double bq = 0.125*(23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]) ;
  double cq = 0.125*((13.41640786499874*psi[5]-15.0*psi[8])*SQ(y)+(6.928203230275509*psi[2]-7.745966692414834*psi[6])*y+5.0*psi[8]- 4.47213595499958*psi[5]-4.47213595499958*psi[4]+4.0*psi[0] ) - psi0;

  double delta2 = bq*bq - 4*aq*cq;

  if (delta2 > 0) {
    double r1, r2;
    double delta = sqrt(delta2);
    // compute both roots
    if (bq>=0) {
      r1 = (-bq-delta)/(2*aq);
      r2 = 2*cq/(-bq-delta);
    }
    else {
      r1 = 2*cq/(-bq+delta);
      r2 = (-bq+delta)/(2*aq);
    }

    int sidx = 0;
    if ((-1<=r1) && (r1 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r1*dx[0]*0.5 + xc[0];

      double x = r1;
      double C = 0.125*(SQ(x)*(90.0*psi[8]*y+23.2379000772445*psi[6])+x*(46.47580015448901*psi[7]*y+12.0*psi[3])+2* (13.41640786499874*psi[5]-15.0*psi[8])*y-7.745966692414834*psi[6]+6.928203230275509*psi[2]) ;
      double A = 0.125*(2*x*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4])+23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]); 
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
    if ((-1<=r2) && (r2 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r2*dx[0]*0.5 + xc[0];

      double x = r2;
      double C = 0.125*(SQ(x)*(90.0*psi[8]*y+23.2379000772445*psi[6])+x*(46.47580015448901*psi[7]*y+12.0*psi[3])+2* (13.41640786499874*psi[5]-15.0*psi[8])*y-7.745966692414834*psi[6]+6.928203230275509*psi[2]) ;
      double A = 0.125*(2*x*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4])+23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]); 
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
  }
  return sol;
}

// Compute roots R(psi,Z) and dR/dZ(psi,Z) in a p=2 DG cell with tensor basis
// Use more accurate roots from numerical recipes in C 2007 section 5.6
static inline struct RdRdZ_sol
calc_RdR_p2_tensor_nrc(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };
  double y = (Z-xc[1])/(dx[1]*0.5);

  double aq = 0.125*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4]);
  double bq = 0.125*(23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]) ;
  double cq = 0.125*((13.41640786499874*psi[5]-15.0*psi[8])*SQ(y)+(6.928203230275509*psi[2]-7.745966692414834*psi[6])*y+5.0*psi[8]- 4.47213595499958*psi[5]-4.47213595499958*psi[4]+4.0*psi[0] ) - psi0;

  double delta2 = bq*bq - 4*aq*cq;

  if (delta2 > 0) {
    double r1, r2;
    double delta = sqrt(delta2);
    //// compute both roots
    double qq = -0.5*(bq + (bq/fabs(bq)) * delta);
    r1 = qq/aq;
    r2 = cq/qq;

    int sidx = 0;
    if ((-1<=r1) && (r1 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r1*dx[0]*0.5 + xc[0];

      double x = r1;
      double C = 0.125*(SQ(x)*(90.0*psi[8]*y+23.2379000772445*psi[6])+x*(46.47580015448901*psi[7]*y+12.0*psi[3])+2* (13.41640786499874*psi[5]-15.0*psi[8])*y-7.745966692414834*psi[6]+6.928203230275509*psi[2]) ;
      double A = 0.125*(2*x*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4])+23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]); 
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      sol.dR[sidx] = -C*dx[0];
      sol.dZ[sidx] = A*dx[1];
      
      sidx += 1;
    }
    if ((-1<=r2) && (r2 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r2*dx[0]*0.5 + xc[0];

      double x = r2;
      double C = 0.125*(SQ(x)*(90.0*psi[8]*y+23.2379000772445*psi[6])+x*(46.47580015448901*psi[7]*y+12.0*psi[3])+2* (13.41640786499874*psi[5]-15.0*psi[8])*y-7.745966692414834*psi[6]+6.928203230275509*psi[2]) ;
      double A = 0.125*(2*x*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4])+23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]); 
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      sol.dR[sidx] = -C*dx[0];
      sol.dZ[sidx] = A*dx[1];
      
      sidx += 1;
    }
  }
  return sol;
}

// Compute roots R(psi,Z) and dR/dZ(psi,Z) in a p=2 DG cell with tensor basis if delta2 is negative but very small
static inline struct RdRdZ_sol
calc_RdR_p2_tensor_with_tolerance(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };
  double y = (Z-xc[1])/(dx[1]*0.5);

  double aq = 0.125*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4]);
  double bq = 0.125*(23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]) ;
  double cq = 0.125*((13.41640786499874*psi[5]-15.0*psi[8])*SQ(y)+(6.928203230275509*psi[2]-7.745966692414834*psi[6])*y+5.0*psi[8]- 4.47213595499958*psi[5]-4.47213595499958*psi[4]+4.0*psi[0] ) - psi0;

  double delta2 = bq*bq - 4*aq*cq;
  if(delta2 > 0)
    return sol;


  if (fabs(delta2) < 1.0e-20) {
    // x = [-b +/- sqrt(b^2 - 4ac)] / 2a
    // If b^2-4ac = 0 then we have one root x = -b/2a
    double r = -bq/2.0/aq;

    int sidx = 0;
    if ((-1<=r) && (r < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r*dx[0]*0.5 + xc[0];

      double x = r;
      double C = 0.125*(SQ(x)*(90.0*psi[8]*y+23.2379000772445*psi[6])+x*(46.47580015448901*psi[7]*y+12.0*psi[3])+2* (13.41640786499874*psi[5]-15.0*psi[8])*y-7.745966692414834*psi[6]+6.928203230275509*psi[2]) ;
      double A = 0.125*(2*x*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4])+23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]); 
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
  }
  return sol;
}

// Compute roots R(psi,Z) and dR/dZ(psi,Z) in a p=2 DG cell with tensor basis if delta2 is negative but very small
static inline struct RdRdZ_sol
calc_RdR_p3(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };
  double y = (Z-xc[1])/(dx[1]*0.5);

  double coeffs[4];
  // coeffs = [x^0, x^1, x^2, x^3]
  coeffs[3] = 0.125*(175.0*psi[15]*CUB(y)+88.74119674649424*psi[13]*SQ(y)+(45.8257569495584*psi[11]-105.0*psi[15])*y+26.45751311064591*psi[8]-29.58039891549808*psi[13]);
  coeffs[2] = 0.125*(88.74119674649424*psi[14]*CUB(y)+45.0*psi[10]*SQ(y)+(23.2379000772445*psi[6]-53.24471804789655*psi[14])*y+13.41640786499874*psi[4]-15.0*psi[10]);
  coeffs[1] = 0.125*((45.8257569495584*psi[12]-105.0*psi[15])*CUB(y)+(23.2379000772445*psi[7]-53.24471804789655*psi[13])*SQ(y)+(12.0*psi[3]+63.0*psi[15]-27.49545416973504*psi[12]-27.49545416973504*psi[11])*y-15.87450786638754*psi[8]-7.745966692414834*psi[7]+17.74823934929885*psi[13]+6.928203230275509*psi[1]);
  coeffs[0] = 0.125*((26.45751311064591*psi[9]-29.58039891549808*psi[14])*CUB(y)+(13.41640786499874*psi[5]-15.0*psi[10])*SQ(y)+(-15.87450786638754*psi[9]-7.745966692414834*psi[6]+6.928203230275509*psi[2]+17.74823934929885*psi[14])*y-4.47213595499958*psi[5]-4.47213595499958*psi[4]+5.0*psi[10]+4.0*psi[0]) - psi0;

  coeffs[0] = coeffs[0]/coeffs[3];
  coeffs[1] = coeffs[1]/coeffs[3];
  coeffs[2] = coeffs[2]/coeffs[3];
  coeffs[3] = coeffs[3]/coeffs[3];

  struct gkyl_lo_poly_roots rts;
  rts = gkyl_calc_lo_poly_roots(GKYL_LO_POLY_3, coeffs);

  int sidx = 0;
  for(int i =0; i<3; i++){
    if(rts.rpart[i] < 1.0 && rts.rpart[i] > -1.0 && fabs(rts.impart[i])<1e-16){
      sol.nsol += 1;
      sol.R[sidx] = rts.rpart[i]*dx[0]*0.5 + xc[0];

      double x = rts.rpart[i];
      double dpsidx = 6.5625000000000000e+01*(x*x)*(y*y*y)*psi[15]+-9.6824583655185426e-01*psi[7]+-6.6555897559870685e+00*(y*y)*psi[13]+5.7282196186947996e+00*(y*y*y)*psi[12]+2.9047375096555625e+00*psi[7]*(y*y)+3.3277948779935343e+01*(x*x)*(y*y)*psi[13]+5.8094750193111251e+00*psi[6]*x*y+-1.3125000000000000e+01*(y*y*y)*psi[15]+9.9215674164922145e+00*(x*x)*psi[8]+-3.7500000000000000e+00*x*psi[10]+8.6602540378443860e-01*psi[1]+-3.4369317712168801e+00*y*psi[12]+-1.3311179511974137e+01*x*psi[14]*y+-1.9843134832984430e+00*psi[8]+-3.4369317712168801e+00*y*psi[11]+-3.9375000000000000e+01*(x*x)*y*psi[15]+1.7184658856084400e+01*(x*x)*y*psi[11]+2.2185299186623562e+00*psi[13]+2.2185299186623560e+01*x*psi[14]*(y*y*y)+1.5000000000000000e+00*y*psi[3]+-1.1092649593311780e+01*(x*x)*psi[13]+1.1250000000000000e+01*x*(y*y)*psi[10]+7.8750000000000000e+00*y*psi[15]+3.3541019662496847e+00*x*psi[4];
      double dpsidy = -9.6824583655185426e-01*psi[6]+1.5000000000000000e+00*x*psi[3]+2.2185299186623560e+01*(x*x*x)*y*psi[13]+7.8750000000000000e+00*x*psi[15]+3.3277948779935343e+01*(x*x)*psi[14]*(y*y)+2.2185299186623562e+00*psi[14]+-3.4369317712168801e+00*x*psi[12]+9.9215674164922145e+00*(y*y)*psi[9]+-3.4369317712168801e+00*x*psi[11]+-1.3311179511974137e+01*x*y*psi[13]+-1.1092649593311780e+01*psi[14]*(y*y)+-3.9375000000000000e+01*x*(y*y)*psi[15]+-1.3125000000000000e+01*(x*x*x)*psi[15]+-3.7500000000000000e+00*y*psi[10]+2.9047375096555625e+00*psi[6]*(x*x)+5.8094750193111251e+00*psi[7]*x*y+-6.6555897559870685e+00*(x*x)*psi[14]+5.7282196186947996e+00*(x*x*x)*psi[11]+-1.9843134832984430e+00*psi[9]+6.5625000000000000e+01*(x*x*x)*(y*y)*psi[15]+3.3541019662496847e+00*psi[5]*y+1.1250000000000000e+01*(x*x)*y*psi[10]+1.7184658856084400e+01*x*(y*y)*psi[12]+8.6602540378443860e-01*psi[2];

      sol.dRdZ[sidx] = -dpsidy/dpsidx*dx[0]/dx[1];
      sol.dR[sidx] = -dpsidy*dx[0];
      sol.dZ[sidx] = dpsidx*dx[1];
      sidx+=1;
    }
  }

  return sol;
}

// hyperbolic roots version
static inline struct RdRdZ_sol
calc_RdR_p3_hyperbolic(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };
  double y = (Z-xc[1])/(dx[1]*0.5);
  double coeffs[4];
  // coeffs = [x^0, x^1, x^2, x^3]
  coeffs[3] = 0.125*(175.0*psi[15]*CUB(y)+88.74119674649424*psi[13]*SQ(y)+(45.8257569495584*psi[11]-105.0*psi[15])*y+26.45751311064591*psi[8]-29.58039891549808*psi[13]);
  coeffs[2] = 0.125*(88.74119674649424*psi[14]*CUB(y)+45.0*psi[10]*SQ(y)+(23.2379000772445*psi[6]-53.24471804789655*psi[14])*y+13.41640786499874*psi[4]-15.0*psi[10]);
  coeffs[1] = 0.125*((45.8257569495584*psi[12]-105.0*psi[15])*CUB(y)+(23.2379000772445*psi[7]-53.24471804789655*psi[13])*SQ(y)+(12.0*psi[3]+63.0*psi[15]-27.49545416973504*psi[12]-27.49545416973504*psi[11])*y-15.87450786638754*psi[8]-7.745966692414834*psi[7]+17.74823934929885*psi[13]+6.928203230275509*psi[1]);
  coeffs[0] = 0.125*((26.45751311064591*psi[9]-29.58039891549808*psi[14])*CUB(y)+(13.41640786499874*psi[5]-15.0*psi[10])*SQ(y)+(-15.87450786638754*psi[9]-7.745966692414834*psi[6]+6.928203230275509*psi[2]+17.74823934929885*psi[14])*y-4.47213595499958*psi[5]-4.47213595499958*psi[4]+5.0*psi[10]+4.0*psi[0]) - psi0;
  double A = coeffs[3];
  double B = coeffs[2]/3.0;
  double C = coeffs[1];
  double D = coeffs[0];
  if(A<1e-10) { // We have a quadratic, use quadratic solve
    double aq = coeffs[2];
    double bq = coeffs[1];
    double cq = coeffs[0];
    double delta2 = bq*bq - 4*aq*cq;
    if (delta2 > 0) {
      double r1, r2;
      double delta = sqrt(delta2);
      //// compute both roots
      double qq = -0.5*(bq + (bq/fabs(bq)) * delta);
      r1 = qq/aq;
      r2 = cq/qq;
      int sidx = 0;
      if ((-1<=r1) && (r1 < 1)) {
        sol.nsol += 1;
        sol.R[sidx] = r1*dx[0]*0.5 + xc[0];
        double x = r1;
        double dpsidx = 6.5625000000000000e+01*(x*x)*(y*y*y)*psi[15]+-9.6824583655185426e-01*psi[7]+-6.6555897559870685e+00*(y*y)*psi[13]+5.7282196186947996e+00*(y*y*y)*psi[12]+2.9047375096555625e+00*psi[7]*(y*y)+3.3277948779935343e+01*(x*x)*(y*y)*psi[13]+5.8094750193111251e+00*psi[6]*x*y+-1.3125000000000000e+01*(y*y*y)*psi[15]+9.9215674164922145e+00*(x*x)*psi[8]+-3.7500000000000000e+00*x*psi[10]+8.6602540378443860e-01*psi[1]+-3.4369317712168801e+00*y*psi[12]+-1.3311179511974137e+01*x*psi[14]*y+-1.9843134832984430e+00*psi[8]+-3.4369317712168801e+00*y*psi[11]+-3.9375000000000000e+01*(x*x)*y*psi[15]+1.7184658856084400e+01*(x*x)*y*psi[11]+2.2185299186623562e+00*psi[13]+2.2185299186623560e+01*x*psi[14]*(y*y*y)+1.5000000000000000e+00*y*psi[3]+-1.1092649593311780e+01*(x*x)*psi[13]+1.1250000000000000e+01*x*(y*y)*psi[10]+7.8750000000000000e+00*y*psi[15]+3.3541019662496847e+00*x*psi[4];
        double dpsidy = -9.6824583655185426e-01*psi[6]+1.5000000000000000e+00*x*psi[3]+2.2185299186623560e+01*(x*x*x)*y*psi[13]+7.8750000000000000e+00*x*psi[15]+3.3277948779935343e+01*(x*x)*psi[14]*(y*y)+2.2185299186623562e+00*psi[14]+-3.4369317712168801e+00*x*psi[12]+9.9215674164922145e+00*(y*y)*psi[9]+-3.4369317712168801e+00*x*psi[11]+-1.3311179511974137e+01*x*y*psi[13]+-1.1092649593311780e+01*psi[14]*(y*y)+-3.9375000000000000e+01*x*(y*y)*psi[15]+-1.3125000000000000e+01*(x*x*x)*psi[15]+-3.7500000000000000e+00*y*psi[10]+2.9047375096555625e+00*psi[6]*(x*x)+5.8094750193111251e+00*psi[7]*x*y+-6.6555897559870685e+00*(x*x)*psi[14]+5.7282196186947996e+00*(x*x*x)*psi[11]+-1.9843134832984430e+00*psi[9]+6.5625000000000000e+01*(x*x*x)*(y*y)*psi[15]+3.3541019662496847e+00*psi[5]*y+1.1250000000000000e+01*(x*x)*y*psi[10]+1.7184658856084400e+01*x*(y*y)*psi[12]+8.6602540378443860e-01*psi[2];
        sol.dRdZ[sidx] = -dpsidy/dpsidx*dx[0]/dx[1];
        sidx += 1;
      }
      if ((-1<=r2) && (r2 < 1)) {
        sol.nsol += 1;
        sol.R[sidx] = r2*dx[0]*0.5 + xc[0];
        double x = r2;
        double dpsidx = 6.5625000000000000e+01*(x*x)*(y*y*y)*psi[15]+-9.6824583655185426e-01*psi[7]+-6.6555897559870685e+00*(y*y)*psi[13]+5.7282196186947996e+00*(y*y*y)*psi[12]+2.9047375096555625e+00*psi[7]*(y*y)+3.3277948779935343e+01*(x*x)*(y*y)*psi[13]+5.8094750193111251e+00*psi[6]*x*y+-1.3125000000000000e+01*(y*y*y)*psi[15]+9.9215674164922145e+00*(x*x)*psi[8]+-3.7500000000000000e+00*x*psi[10]+8.6602540378443860e-01*psi[1]+-3.4369317712168801e+00*y*psi[12]+-1.3311179511974137e+01*x*psi[14]*y+-1.9843134832984430e+00*psi[8]+-3.4369317712168801e+00*y*psi[11]+-3.9375000000000000e+01*(x*x)*y*psi[15]+1.7184658856084400e+01*(x*x)*y*psi[11]+2.2185299186623562e+00*psi[13]+2.2185299186623560e+01*x*psi[14]*(y*y*y)+1.5000000000000000e+00*y*psi[3]+-1.1092649593311780e+01*(x*x)*psi[13]+1.1250000000000000e+01*x*(y*y)*psi[10]+7.8750000000000000e+00*y*psi[15]+3.3541019662496847e+00*x*psi[4];
        double dpsidy = -9.6824583655185426e-01*psi[6]+1.5000000000000000e+00*x*psi[3]+2.2185299186623560e+01*(x*x*x)*y*psi[13]+7.8750000000000000e+00*x*psi[15]+3.3277948779935343e+01*(x*x)*psi[14]*(y*y)+2.2185299186623562e+00*psi[14]+-3.4369317712168801e+00*x*psi[12]+9.9215674164922145e+00*(y*y)*psi[9]+-3.4369317712168801e+00*x*psi[11]+-1.3311179511974137e+01*x*y*psi[13]+-1.1092649593311780e+01*psi[14]*(y*y)+-3.9375000000000000e+01*x*(y*y)*psi[15]+-1.3125000000000000e+01*(x*x*x)*psi[15]+-3.7500000000000000e+00*y*psi[10]+2.9047375096555625e+00*psi[6]*(x*x)+5.8094750193111251e+00*psi[7]*x*y+-6.6555897559870685e+00*(x*x)*psi[14]+5.7282196186947996e+00*(x*x*x)*psi[11]+-1.9843134832984430e+00*psi[9]+6.5625000000000000e+01*(x*x*x)*(y*y)*psi[15]+3.3541019662496847e+00*psi[5]*y+1.1250000000000000e+01*(x*x)*y*psi[10]+1.7184658856084400e+01*x*(y*y)*psi[12]+8.6602540378443860e-01*psi[2];
        sol.dRdZ[sidx] = -dpsidy/dpsidx*dx[0]/dx[1];
        sidx += 1;
      }
    }
    return sol;
  }
  double gamma = B/A;
  double a = C/3/A - SQ(gamma);
  double b = 2*CUB(gamma) -C*gamma/A + D/A;
  double complex rho = -4.0*a + 0.0*I;
  double complex s = -4.0*b + 0.0*I;
  double complex t = csqrt(cSQ(s) - cCUB(rho));
  double complex alpha = -0.5 + 0.5*sqrt(3.0)*I;
  double complex alpha_0 = 1.0;
  double complex alpha_1 = alpha;
  double complex alpha_2 = cSQ(alpha);
  double complex x_0 = (0.5 + 0.0*I)* (alpha_0*croot(s+t, 3.0) + rho/alpha_0/croot(s+t,3.0) );
  double complex x_1 = (0.5 + 0.0*I)* (alpha_1*croot(s+t, 3.0) + rho/alpha_1/croot(s+t,3.0) );
  double complex x_2 = (0.5 + 0.0*I)* (alpha_2*croot(s+t, 3.0) + rho/alpha_2/croot(s+t,3.0) );
  double complex r0 = x_0 - (gamma + 0.0*I);
  double complex r1 = x_1 - (gamma + 0.0*I);
  double complex r2 = x_2 - (gamma + 0.0*I);
  double complex roots[3] = {r0, r1, r2};
  int sidx = 0;
  for(int i =0; i<3; i++){
    double rpart = creal(roots[i]);
    double impart = cimag(roots[i]);
    if(rpart < 1.0 && rpart > -1.0 && fabs(impart)<1e-10){
      sol.nsol += 1;
      sol.R[sidx] = rpart*dx[0]*0.5 + xc[0];
      double x = rpart;
      double dpsidx = 6.5625000000000000e+01*(x*x)*(y*y*y)*psi[15]+-9.6824583655185426e-01*psi[7]+-6.6555897559870685e+00*(y*y)*psi[13]+5.7282196186947996e+00*(y*y*y)*psi[12]+2.9047375096555625e+00*psi[7]*(y*y)+3.3277948779935343e+01*(x*x)*(y*y)*psi[13]+5.8094750193111251e+00*psi[6]*x*y+-1.3125000000000000e+01*(y*y*y)*psi[15]+9.9215674164922145e+00*(x*x)*psi[8]+-3.7500000000000000e+00*x*psi[10]+8.6602540378443860e-01*psi[1]+-3.4369317712168801e+00*y*psi[12]+-1.3311179511974137e+01*x*psi[14]*y+-1.9843134832984430e+00*psi[8]+-3.4369317712168801e+00*y*psi[11]+-3.9375000000000000e+01*(x*x)*y*psi[15]+1.7184658856084400e+01*(x*x)*y*psi[11]+2.2185299186623562e+00*psi[13]+2.2185299186623560e+01*x*psi[14]*(y*y*y)+1.5000000000000000e+00*y*psi[3]+-1.1092649593311780e+01*(x*x)*psi[13]+1.1250000000000000e+01*x*(y*y)*psi[10]+7.8750000000000000e+00*y*psi[15]+3.3541019662496847e+00*x*psi[4];
      double dpsidy = -9.6824583655185426e-01*psi[6]+1.5000000000000000e+00*x*psi[3]+2.2185299186623560e+01*(x*x*x)*y*psi[13]+7.8750000000000000e+00*x*psi[15]+3.3277948779935343e+01*(x*x)*psi[14]*(y*y)+2.2185299186623562e+00*psi[14]+-3.4369317712168801e+00*x*psi[12]+9.9215674164922145e+00*(y*y)*psi[9]+-3.4369317712168801e+00*x*psi[11]+-1.3311179511974137e+01*x*y*psi[13]+-1.1092649593311780e+01*psi[14]*(y*y)+-3.9375000000000000e+01*x*(y*y)*psi[15]+-1.3125000000000000e+01*(x*x*x)*psi[15]+-3.7500000000000000e+00*y*psi[10]+2.9047375096555625e+00*psi[6]*(x*x)+5.8094750193111251e+00*psi[7]*x*y+-6.6555897559870685e+00*(x*x)*psi[14]+5.7282196186947996e+00*(x*x*x)*psi[11]+-1.9843134832984430e+00*psi[9]+6.5625000000000000e+01*(x*x*x)*(y*y)*psi[15]+3.3541019662496847e+00*psi[5]*y+1.1250000000000000e+01*(x*x)*y*psi[10]+1.7184658856084400e+01*x*(y*y)*psi[12]+8.6602540378443860e-01*psi[2];
      sol.dRdZ[sidx] = -dpsidy/dpsidx*dx[0]/dx[1];
      sol.dR[sidx] = -dpsidy*dx[0];
      sol.dZ[sidx] = dpsidx*dx[1];
      sidx+=1;
    }
  }
  return sol;
}


// Compute R(psi,Z) given a psi and Z. Can return multiple solutions
// or no solutions. The number of roots found is returned and are
// copied in the array R and dR. The calling function must ensure that
// these arrays are big enough to hold all roots required
static int
R_psiZ(const struct gkyl_tok_geo *geo, double psi, double Z, int nmaxroots,
  double *R, double *dRdZ, double* dR, double *dZ)
{
  int zcell = get_idx(1, Z, &geo->rzgrid, &geo->rzlocal);

  int sidx = 0;
  int idx[2] = { 0, zcell };
  double dx[2] = { geo->rzgrid.dx[0], geo->rzgrid.dx[1] };
  
  struct gkyl_range rangeR;
  gkyl_range_deflate(&rangeR, &geo->rzlocal, (int[]) { 0, 1 }, (int[]) { 0, zcell });

  struct gkyl_range_iter riter;
  gkyl_range_iter_init(&riter, &rangeR);
  
  // loop over all R cells to find psi crossing
  while (gkyl_range_iter_next(&riter) && sidx<=nmaxroots) {
    long loc = gkyl_range_idx(&rangeR, riter.idx);
    const double *psih = gkyl_array_cfetch(geo->psiRZ, loc);

    double xc[2];
    idx[0] = riter.idx[0];
    gkyl_rect_grid_cell_center(&geo->rzgrid, idx, xc);

    struct RdRdZ_sol sol = geo->calc_roots(psih, psi, Z, xc, dx);
    
    if (sol.nsol > 0)
      for (int s=0; s<sol.nsol; ++s) {
        if( (sol.R[s] > geo->rmin) && (sol.R[s] < geo->rmax) ) {
          R[sidx] = sol.R[s];
          dRdZ[sidx] = sol.dRdZ[s];
          dR[sidx] = sol.dR[s];
          dZ[sidx] = sol.dZ[s];
          sidx += 1;
        }
      }
  }

  // Try again if we didn't find any
  if (sidx==0 && geo->inexact_roots) {
    gkyl_range_iter_init(&riter, &rangeR);
    while (gkyl_range_iter_next(&riter) && sidx<=nmaxroots) {
      long loc = gkyl_range_idx(&rangeR, riter.idx);
      const double *psih = gkyl_array_cfetch(geo->psiRZ, loc);

      double xc[2];
      idx[0] = riter.idx[0];
      gkyl_rect_grid_cell_center(&geo->rzgrid, idx, xc);

      struct RdRdZ_sol sol = calc_RdR_p2_tensor_with_tolerance(psih, psi, Z, xc, dx);
      
      if (sol.nsol > 0)
        for (int s=0; s<sol.nsol; ++s) {
          if( (sol.R[s] > geo->rmin) && (sol.R[s] < geo->rmax) ) {
            R[sidx] = sol.R[s];
            dRdZ[sidx] = sol.dRdZ[s];
            dR[sidx] = sol.dR[s];
            dZ[sidx] = sol.dZ[s];
            sidx += 1;
          }
        }
    }
  }

  return sidx;
}

// Compute R(psi,Z) given a psi and Z. Can return multiple solutions
// or no solutions. The number of roots found is returned and are
// copied in the array R and dR. The calling function must ensure that
// these arrays are big enough to hold all roots required
static int
R_psiZ_cubic(const struct gkyl_tok_geo *geo, double psi, double Z, int nmaxroots,
  double *R, double *dRdZ, double* dR, double *dZ)
{
  int zcell = get_idx(1, Z, &geo->rzgrid_cubic, &geo->rzlocal_cubic);

  int sidx = 0;
  int idx[2] = { 0, zcell };
  double dx[2] = { geo->rzgrid_cubic.dx[0], geo->rzgrid_cubic.dx[1] };

  struct gkyl_range rangeR;
  gkyl_range_deflate(&rangeR, &geo->rzlocal_cubic, (int[]) { 0, 1 }, (int[]) { 0, zcell });

  struct gkyl_range_iter riter;
  gkyl_range_iter_init(&riter, &rangeR);

  // loop over all R cells to find psi crossing
  while (gkyl_range_iter_next(&riter) && sidx<=nmaxroots) {
    long loc = gkyl_range_idx(&rangeR, riter.idx);
    const double *psih = gkyl_array_cfetch(geo->psiRZ_cubic, loc);

    double xc[2];
    idx[0] = riter.idx[0];
    gkyl_rect_grid_cell_center(&geo->rzgrid_cubic, idx, xc);

    struct RdRdZ_sol sol = geo->calc_roots(psih, psi, Z, xc, dx);
    
    if (sol.nsol > 0)
      for (int s=0; s<sol.nsol; ++s) {
        if( (sol.R[s] > geo->rmin) && (sol.R[s] < geo->rmax) ) {
          dRdZ[sidx] = sol.dRdZ[s];
          dR[sidx] = sol.dR[s];
          dZ[sidx] = sol.dZ[s];
          sidx += 1;
        }
      }
  }

  return sidx;
}

static double
calc_grad_psi_p1(const double *psih, const double eta[2], const double dx[2])
{
  double x = eta[0];
  double y = eta[1];
  double dpsidx = 1.5*psih[3]*y+0.8660254037844386*psih[1];
  double dpsidy = 1.5*psih[3]*x+0.8660254037844386*psih[2];
  dpsidx = dpsidx*2.0/dx[0];
  dpsidy = dpsidy*2.0/dx[1];
  return sqrt(dpsidx*dpsidx + dpsidy*dpsidy);
}

static double
calc_grad_psi_p2_tensor(const double *psih, const double eta[2], const double dx[2])
{
  double x = eta[0];
  double y = eta[1];
  double dpsidx = 5.625*psih[8]*(2.0*x*SQ(y)-0.6666666666666666*x)+2.904737509655563*psih[7]*(SQ(y)-0.3333333333333333)+5.809475019311126*psih[6]*x*y+1.5*psih[3]*y+3.354101966249684*psih[4]*x+0.8660254037844386*psih[1];
  double dpsidy = 5.625*psih[8]*(2.0*SQ(x)*y-0.6666666666666666*y)+5.809475019311126*psih[7]*x*y+3.354101966249684*psih[5]*y+2.904737509655563*psih[6]*(SQ(x)-0.3333333333333333)+1.5*psih[3]*x+0.8660254037844386*psih[2];
  dpsidx = dpsidx*2.0/dx[0];
  dpsidy = dpsidy*2.0/dx[1];
  return sqrt(dpsidx*dpsidx + dpsidy*dpsidy);
}

static double
calc_grad_psi_p3(const double *psih, const double eta[2], const double dx[2])
{
  double x = eta[0];
  double y = eta[1];
  double dpsidx = 6.5625000000000000e+01*(x*x)*(y*y*y)*psih[15]+-9.6824583655185426e-01*psih[7]+-6.6555897559870685e+00*(y*y)*psih[13]+5.7282196186947996e+00*(y*y*y)*psih[12]+2.9047375096555625e+00*psih[7]*(y*y)+3.3277948779935343e+01*(x*x)*(y*y)*psih[13]+5.8094750193111251e+00*psih[6]*x*y+-1.3125000000000000e+01*(y*y*y)*psih[15]+9.9215674164922145e+00*(x*x)*psih[8]+-3.7500000000000000e+00*x*psih[10]+8.6602540378443860e-01*psih[1]+-3.4369317712168801e+00*y*psih[12]+-1.3311179511974137e+01*x*psih[14]*y+-1.9843134832984430e+00*psih[8]+-3.4369317712168801e+00*y*psih[11]+-3.9375000000000000e+01*(x*x)*y*psih[15]+1.7184658856084400e+01*(x*x)*y*psih[11]+2.2185299186623562e+00*psih[13]+2.2185299186623560e+01*x*psih[14]*(y*y*y)+1.5000000000000000e+00*y*psih[3]+-1.1092649593311780e+01*(x*x)*psih[13]+1.1250000000000000e+01*x*(y*y)*psih[10]+7.8750000000000000e+00*y*psih[15]+3.3541019662496847e+00*x*psih[4];
  double dpsidy = -9.6824583655185426e-01*psih[6]+1.5000000000000000e+00*x*psih[3]+2.2185299186623560e+01*(x*x*x)*y*psih[13]+7.8750000000000000e+00*x*psih[15]+3.3277948779935343e+01*(x*x)*psih[14]*(y*y)+2.2185299186623562e+00*psih[14]+-3.4369317712168801e+00*x*psih[12]+9.9215674164922145e+00*(y*y)*psih[9]+-3.4369317712168801e+00*x*psih[11]+-1.3311179511974137e+01*x*y*psih[13]+-1.1092649593311780e+01*psih[14]*(y*y)+-3.9375000000000000e+01*x*(y*y)*psih[15]+-1.3125000000000000e+01*(x*x*x)*psih[15]+-3.7500000000000000e+00*y*psih[10]+2.9047375096555625e+00*psih[6]*(x*x)+5.8094750193111251e+00*psih[7]*x*y+-6.6555897559870685e+00*(x*x)*psih[14]+5.7282196186947996e+00*(x*x*x)*psih[11]+-1.9843134832984430e+00*psih[9]+6.5625000000000000e+01*(x*x*x)*(y*y)*psih[15]+3.3541019662496847e+00*psih[5]*y+1.1250000000000000e+01*(x*x)*y*psih[10]+1.7184658856084400e+01*x*(y*y)*psih[12]+8.6602540378443860e-01*psih[2];
  dpsidx = dpsidx*2.0/dx[0];
  dpsidy = dpsidy*2.0/dx[1];
  return sqrt(dpsidx*dpsidx + dpsidy*dpsidy);
}


// Function context to pass to coutour integration function
struct contour_ctx {
  const struct gkyl_tok_geo *geo;
  double psi, last_R;
  long ncall;
};

// Function to pass to numerical quadrature to integrate along a contour
static inline double
contour_func(double Z, void *ctx)
{
  struct contour_ctx *c = ctx;
  c->ncall += 1;
  double R[4] = { 0 }, dRdZ[4] = { 0 };
  double dR[4] = { 0 }, dZ[4] = { 0 };
  
  int nr = gkyl_tok_geo_R_psiZ(c->geo, c->psi, Z, 4, R, dRdZ, dR, dZ);
  double drdz = nr == 1 ? dRdZ[0] : choose_closest(c->last_R, R, dRdZ,nr);
  
  return nr>0 ? sqrt(1+drdz*drdz) : 0.0;
}

static inline double
phi_contour_func(double Z, void *ctx)
{
  struct contour_ctx *c = ctx;
  c->ncall += 1;
  double R[4] = { 0 }, dRdZ[4] = { 0 };
  double dR[4] = { 0 }, dZ[4] = { 0 };
  
  int nr = gkyl_tok_geo_R_psiZ(c->geo, c->psi, Z, 4, R, dRdZ, dR, dZ);
  double drdz = nr == 1 ? dRdZ[0] : choose_closest(c->last_R, R, dRdZ, nr);
  double r_curr = nr == 1 ? R[0] : choose_closest(c->last_R, R, R, nr);

  if (c->geo->use_cubics) {
    double xn[2] = {r_curr, Z};
    double fout[3];
    c->geo->efit->evf->eval_cubic_wgrad(0.0, xn, fout, c->geo->efit->evf->ctx);
    double dpsidR = fout[1];
    double dpsidZ = fout[2]; 
    double grad_psi_mag = sqrt(dpsidR*dpsidR + dpsidZ*dpsidZ);

    double result  = (1/r_curr/grad_psi_mag) *sqrt(1+drdz*drdz) ;
    return nr>0 ? result : 0.0;
  }
  else {
    int rzidx[2];
    int idxtemp = c->geo->rzlocal.lower[0] + (int) floor((r_curr - c->geo->rzgrid.lower[0])/c->geo->rzgrid.dx[0]);
    idxtemp = GKYL_MIN2(idxtemp, c->geo->rzlocal.upper[0]);
    idxtemp = GKYL_MAX2(idxtemp, c->geo->rzlocal.lower[0]);
    rzidx[0] = idxtemp;
    idxtemp = c->geo->rzlocal.lower[1] + (int) floor((Z - c->geo->rzgrid.lower[1])/c->geo->rzgrid.dx[1]);
    idxtemp = GKYL_MIN2(idxtemp, c->geo->rzlocal.upper[1]);
    idxtemp = GKYL_MAX2(idxtemp, c->geo->rzlocal.lower[1]);
    rzidx[1] = idxtemp;

    long loc = gkyl_range_idx((&c->geo->rzlocal), rzidx);
    const double *psih = gkyl_array_cfetch(c->geo->psiRZ, loc);

    double xc[2];
    gkyl_rect_grid_cell_center((&c->geo->rzgrid), rzidx, xc);
    double x = (r_curr-xc[0])/(c->geo->rzgrid.dx[0]*0.5);
    double y = (Z-xc[1])/(c->geo->rzgrid.dx[1]*0.5);

    double eta[2] = {x,y};
    double grad_psi_mag = c->geo->calc_grad_psi(psih, eta, c->geo->rzgrid.dx);

    double result  = (1/r_curr/grad_psi_mag) *sqrt(1+drdz*drdz) ;
    return nr>0 ? result : 0.0;
  }
}

static inline double
dphidtheta_integrand(double Z, void *ctx)
{
  struct contour_ctx *c = ctx;
  c->ncall += 1;
  double R[4] = { 0 }, dRdZ[4] = { 0 };
  double dR[4] = { 0 }, dZ[4] = { 0 };
  
  int nr = gkyl_tok_geo_R_psiZ(c->geo, c->psi, Z, 4, R, dRdZ, dR, dZ);
  double drdz = nr == 1 ? dRdZ[0] : choose_closest(c->last_R, R, dRdZ, nr);
  double r_curr = nr == 1 ? R[0] : choose_closest(c->last_R, R, R, nr);

  if (c->geo->use_cubics) {
    double xn[2] = {r_curr, Z};
    double fout[3];
    c->geo->efit->evf->eval_cubic_wgrad(0.0, xn, fout, c->geo->efit->evf->ctx);
    double dpsidR = fout[1];
    double dpsidZ = fout[2];
    double grad_psi_mag = sqrt(dpsidR*dpsidR + dpsidZ*dpsidZ);

    double result  = (1/r_curr/grad_psi_mag);
    return nr>0 ? result : 0.0;
  }
  else {
    int rzidx[2];
    int idxtemp = c->geo->rzlocal.lower[0] + (int) floor((r_curr - c->geo->rzgrid.lower[0])/c->geo->rzgrid.dx[0]);
    idxtemp = GKYL_MIN2(idxtemp, c->geo->rzlocal.upper[0]);
    idxtemp = GKYL_MAX2(idxtemp, c->geo->rzlocal.lower[0]);
    rzidx[0] = idxtemp;
    idxtemp = c->geo->rzlocal.lower[1] + (int) floor((Z - c->geo->rzgrid.lower[1])/c->geo->rzgrid.dx[1]);
    idxtemp = GKYL_MIN2(idxtemp, c->geo->rzlocal.upper[1]);
    idxtemp = GKYL_MAX2(idxtemp, c->geo->rzlocal.lower[1]);
    rzidx[1] = idxtemp;

    long loc = gkyl_range_idx((&c->geo->rzlocal), rzidx);
    const double *psih = gkyl_array_cfetch(c->geo->psiRZ, loc);

    double xc[2];
    gkyl_rect_grid_cell_center((&c->geo->rzgrid), rzidx, xc);
    double x = (r_curr-xc[0])/(c->geo->rzgrid.dx[0]*0.5);
    double y = (Z-xc[1])/(c->geo->rzgrid.dx[1]*0.5);

    double eta[2] = {x,y};
    double grad_psi_mag = c->geo->calc_grad_psi(psih, eta, c->geo->rzgrid.dx);
    double result  = (1/r_curr/grad_psi_mag);
    return nr>0 ? result : 0.0;
  }
}

// Integrates along a specified contour, optionally using a "memory"
// of previously computed values, or storing computed values in
// memory. The function basically breaks up the integral into a loop
// over z-cells. This needs to be done as the DG representation is,
// well, discontinuous, and adaptive quadrature struggles with such
// functions.
static double
integrate_psi_contour_memo(const struct gkyl_tok_geo *geo, double psi,
  double zmin, double zmax, double rclose,
  bool use_memo, bool fill_memo, double *memo)
{
  struct contour_ctx ctx = {
    .geo = geo,
    .psi = psi,
    .ncall = 0,
    .last_R = rclose
  };

  int nlevels = geo->quad_param.max_level;
  double eps = geo->quad_param.eps;

  struct gkyl_rect_grid rzgrid ;
  struct gkyl_range rzlocal;
  if(geo->use_cubics) {
    rzgrid = geo->rzgrid_cubic;
    rzlocal = geo->rzlocal_cubic;
  }
  else {
    rzgrid = geo->rzgrid;
    rzlocal = geo->rzlocal;
  }
  
  double dz = rzgrid.dx[1];
  double zlo = rzgrid.lower[1];
  int izlo = rzlocal.lower[1], izup = rzlocal.upper[1];
  
  int ilo = get_idx(1, zmin, &rzgrid, &rzlocal);
  int iup = get_idx(1, zmax, &rzgrid, &rzlocal);

  double res = 0.0;
  for (int i=ilo; i<=iup; ++i) {
    double z1 = gkyl_median(zmin, zlo+(i-izlo)*dz, zlo+(i-izlo+1)*dz);
    double z2 = gkyl_median(zmax, zlo+(i-izlo)*dz, zlo+(i-izlo+1)*dz);
    
    if (z1 < z2) {
      if (use_memo) {
        if (fill_memo) {
          struct gkyl_qr_res res_local =
            gkyl_dbl_exp(contour_func, &ctx, z1, z2, nlevels, eps);
          memo[i-izlo] = res_local.res;
          res += res_local.res;
        }
        else {
          if (z2-z1 == dz) {
            res += memo[i-izlo];
          }
          else {
            struct gkyl_qr_res res_local =
              gkyl_dbl_exp(contour_func, &ctx, z1, z2, nlevels, eps);
            res += res_local.res;
          }
        }
      }
      else {
        struct gkyl_qr_res res_local =
          gkyl_dbl_exp(contour_func, &ctx, z1, z2, nlevels, eps);
        res += res_local.res;
      }
    }
  }

  ((struct gkyl_tok_geo *)geo)->stat.nquad_cont_calls += ctx.ncall;
  return res;
}

static double
integrate_phi_along_psi_contour_memo(const struct gkyl_tok_geo *geo, double psi,
  double zmin, double zmax, double rclose,
  bool use_memo, bool fill_memo, double *memo)
{
  struct contour_ctx ctx = {
    .geo = geo,
    .psi = psi,
    .ncall = 0,
    .last_R = rclose
  };

  int nlevels = geo->quad_param.max_level;
  double eps = geo->quad_param.eps;

  struct gkyl_rect_grid rzgrid ;
  struct gkyl_range rzlocal;
  if(geo->use_cubics) {
    rzgrid = geo->rzgrid_cubic;
    rzlocal = geo->rzlocal_cubic;
  }
  else {
    rzgrid = geo->rzgrid;
    rzlocal = geo->rzlocal;
  }
  
  double dz = rzgrid.dx[1];
  double zlo = rzgrid.lower[1];
  int izlo = rzlocal.lower[1], izup = rzlocal.upper[1];
  
  int ilo = get_idx(1, zmin, &rzgrid, &rzlocal);
  int iup = get_idx(1, zmax, &rzgrid, &rzlocal);

  double res = 0.0;
  for (int i=ilo; i<=iup; ++i) {
    double z1 = gkyl_median(zmin, zlo+(i-izlo)*dz, zlo+(i-izlo+1)*dz);
    double z2 = gkyl_median(zmax, zlo+(i-izlo)*dz, zlo+(i-izlo+1)*dz);
    
    if (z1 < z2) {
      if (use_memo) {
        if (fill_memo) {
          struct gkyl_qr_res res_local =
            gkyl_dbl_exp(phi_contour_func, &ctx, z1, z2, nlevels, eps);
          memo[i-izlo] = res_local.res;
          res += res_local.res;
        }
        else {
          if (z2-z1 == dz) {
            res += memo[i-izlo];
          }
          else {
            struct gkyl_qr_res res_local =
              gkyl_dbl_exp(phi_contour_func, &ctx, z1, z2, nlevels, eps);
            res += res_local.res;
          }
        }
      }
      else {
        struct gkyl_qr_res res_local =
          gkyl_dbl_exp(phi_contour_func, &ctx, z1, z2, nlevels, eps);
        res += res_local.res;
      }
    }
  }

  ((struct gkyl_tok_geo *)geo)->stat.nquad_cont_calls += ctx.ncall;
  return res;
}


double phi_func(double alpha_curr, double Z, void *ctx);
double tok_plate_psi_func(double s, void *ctx);

/*
 * Used to set zmin and zmax and attributes of arc_ctx before looping over arc length
*/
void tok_find_endpoints(struct gkyl_tok_geo_grid_inp* inp, struct gkyl_tok_geo *geo, struct arc_length_ctx* arc_ctx, struct plate_ctx* pctx, double psi_curr, double alpha_curr, double* arc_memo, double* arc_memo_left, double* arc_memo_right);

/*
 * Used to set theta extents when using a global normalization factor
 * */
void tok_geo_set_extent(struct gkyl_tok_geo_grid_inp* inp, struct gkyl_tok_geo *geo, double *theta_lo, double *theta_up);


/*
 * Used to set arc_ctx attributes before using ridders to find z
*/
void tok_set_ridders(struct gkyl_tok_geo_grid_inp* inp, struct arc_length_ctx* arc_ctx, double psi_curr, double arcL_curr, double* rclose, double *ridders_min, double* ridders_max);

