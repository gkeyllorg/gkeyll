#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_math.h>
#include <gkyl_efit.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_basis_ops.h>

#include <gkyl_tok_geo.h>
#include <gkyl_tok_geo_priv.h>

#include <complex.h>
#include <math.h>

// Compute roots R(psi,Z) and dR/dZ(psi,Z) in a p=2 DG cell with tensor basis
// Use more accurate roots from numerical recipes in C 2007 section 5.6
static inline struct RdRdZ_sol
quad_root(const double *psi, double psi0, double Z, double xc[2], double dx[2])
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

static inline double cub(double x) { return x*x*x; }
static inline struct RdRdZ_sol
cub_root(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };
  double y = (Z-xc[1])/(dx[1]*0.5);

  double coeffs[4];
  // coeffs = [x^0, x^1, x^2, x^3]
  coeffs[3] = 0.125*(175.0*psi[15]*cub(y)+88.74119674649424*psi[13]*SQ(y)+(45.8257569495584*psi[11]-105.0*psi[15])*y+26.45751311064591*psi[8]-29.58039891549808*psi[13]);
  coeffs[2] = 0.125*(88.74119674649424*psi[14]*cub(y)+45.0*psi[10]*SQ(y)+(23.2379000772445*psi[6]-53.24471804789655*psi[14])*y+13.41640786499874*psi[4]-15.0*psi[10]);
  coeffs[1] = 0.125*((45.8257569495584*psi[12]-105.0*psi[15])*cub(y)+(23.2379000772445*psi[7]-53.24471804789655*psi[13])*SQ(y)+(12.0*psi[3]+63.0*psi[15]-27.49545416973504*psi[12]-27.49545416973504*psi[11])*y-15.87450786638754*psi[8]-7.745966692414834*psi[7]+17.74823934929885*psi[13]+6.928203230275509*psi[1]);
  coeffs[0] = 0.125*((26.45751311064591*psi[9]-29.58039891549808*psi[14])*cub(y)+(13.41640786499874*psi[5]-15.0*psi[10])*SQ(y)+(-15.87450786638754*psi[9]-7.745966692414834*psi[6]+6.928203230275509*psi[2]+17.74823934929885*psi[14])*y-4.47213595499958*psi[5]-4.47213595499958*psi[4]+5.0*psi[10]+4.0*psi[0]) - psi0;

  coeffs[0] = coeffs[0]/coeffs[3];
  coeffs[1] = coeffs[1]/coeffs[3];
  coeffs[2] = coeffs[2]/coeffs[3];
  coeffs[3] = coeffs[3]/coeffs[3];

  struct gkyl_lo_poly_roots rts;
  rts = gkyl_calc_lo_poly_roots(GKYL_LO_POLY_3, coeffs);

  int sidx = 0;
  for(int i =0; i<3; i++){
    if(rts.rpart[i] < 1.0 && rts.rpart[i] > -1.0 && fabs(rts.impart[i])<1e-10){
      sol.nsol += 1;
      sol.R[sidx] = rts.rpart[i]*dx[0]*0.5 + xc[0];

      double x = rts.rpart[i];
      double dpsidx = 6.5625000000000000e+01*(x*x)*(y*y*y)*psi[15]+-9.6824583655185426e-01*psi[7]+-6.6555897559870685e+00*(y*y)*psi[13]+5.7282196186947996e+00*(y*y*y)*psi[12]+2.9047375096555625e+00*psi[7]*(y*y)+3.3277948779935343e+01*(x*x)*(y*y)*psi[13]+5.8094750193111251e+00*psi[6]*x*y+-1.3125000000000000e+01*(y*y*y)*psi[15]+9.9215674164922145e+00*(x*x)*psi[8]+-3.7500000000000000e+00*x*psi[10]+8.6602540378443860e-01*psi[1]+-3.4369317712168801e+00*y*psi[12]+-1.3311179511974137e+01*x*psi[14]*y+-1.9843134832984430e+00*psi[8]+-3.4369317712168801e+00*y*psi[11]+-3.9375000000000000e+01*(x*x)*y*psi[15]+1.7184658856084400e+01*(x*x)*y*psi[11]+2.2185299186623562e+00*psi[13]+2.2185299186623560e+01*x*psi[14]*(y*y*y)+1.5000000000000000e+00*y*psi[3]+-1.1092649593311780e+01*(x*x)*psi[13]+1.1250000000000000e+01*x*(y*y)*psi[10]+7.8750000000000000e+00*y*psi[15]+3.3541019662496847e+00*x*psi[4];
      double dpsidy = -9.6824583655185426e-01*psi[6]+1.5000000000000000e+00*x*psi[3]+2.2185299186623560e+01*(x*x*x)*y*psi[13]+7.8750000000000000e+00*x*psi[15]+3.3277948779935343e+01*(x*x)*psi[14]*(y*y)+2.2185299186623562e+00*psi[14]+-3.4369317712168801e+00*x*psi[12]+9.9215674164922145e+00*(y*y)*psi[9]+-3.4369317712168801e+00*x*psi[11]+-1.3311179511974137e+01*x*y*psi[13]+-1.1092649593311780e+01*psi[14]*(y*y)+-3.9375000000000000e+01*x*(y*y)*psi[15]+-1.3125000000000000e+01*(x*x*x)*psi[15]+-3.7500000000000000e+00*y*psi[10]+2.9047375096555625e+00*psi[6]*(x*x)+5.8094750193111251e+00*psi[7]*x*y+-6.6555897559870685e+00*(x*x)*psi[14]+5.7282196186947996e+00*(x*x*x)*psi[11]+-1.9843134832984430e+00*psi[9]+6.5625000000000000e+01*(x*x*x)*(y*y)*psi[15]+3.3541019662496847e+00*psi[5]*y+1.1250000000000000e+01*(x*x)*y*psi[10]+1.7184658856084400e+01*x*(y*y)*psi[12]+8.6602540378443860e-01*psi[2];

      sol.dRdZ[sidx] = -dpsidy/dpsidx*dx[0]/dx[1];
      sidx+=1;
    }
  }

  return sol;
}


static int
getR(const struct gkyl_range rzlocal, const struct gkyl_rect_grid rzgrid, struct gkyl_array *psiRZ, double psi, double Z, int nmaxroots,
  double *R, double *dR)
{
  int zcell = get_idx(1, Z, &rzgrid, &rzlocal);

  int sidx = 0;
  int idx[2] = { 0, zcell };
  double dx[2] = {rzgrid.dx[0], rzgrid.dx[1] };
  
  struct gkyl_range rangeR;
  gkyl_range_deflate(&rangeR, &rzlocal, (int[]) { 0, 1 }, (int[]) { 0, zcell });

  struct gkyl_range_iter riter;
  gkyl_range_iter_init(&riter, &rangeR);
  
  // loop over all R cells to find psi crossing
  while (gkyl_range_iter_next(&riter) && sidx<=nmaxroots) {
    long loc = gkyl_range_idx(&rangeR, riter.idx);
    const double *psih = gkyl_array_cfetch(psiRZ, loc);

    double xc[2];
    idx[0] = riter.idx[0];
    gkyl_rect_grid_cell_center(&rzgrid, idx, xc);

    struct RdRdZ_sol sol = quad_root(psih, psi, Z, xc, dx);
    
    if (sol.nsol > 0)
      for (int s=0; s<sol.nsol; ++s) {
          R[sidx] = sol.R[s];
          dR[sidx] = sol.dRdZ[s];
          sidx += 1;
      }
  }
  return sidx;
}

static int
getRcub(const struct gkyl_range rzlocal, const struct gkyl_rect_grid rzgrid, struct gkyl_array *psiRZ, double psi, double Z, int nmaxroots,
  double *R, double *dR)
{
  int zcell = get_idx(1, Z, &rzgrid, &rzlocal);

  int sidx = 0;
  int idx[2] = { 0, zcell };
  double dx[2] = {rzgrid.dx[0], rzgrid.dx[1] };
  
  struct gkyl_range rangeR;
  gkyl_range_deflate(&rangeR, &rzlocal, (int[]) { 0, 1 }, (int[]) { 0, zcell });

  struct gkyl_range_iter riter;
  gkyl_range_iter_init(&riter, &rangeR);
  
  // loop over all R cells to find psi crossing
  while (gkyl_range_iter_next(&riter) && sidx<=nmaxroots) {
    long loc = gkyl_range_idx(&rangeR, riter.idx);
    const double *psih = gkyl_array_cfetch(psiRZ, loc);

    double xc[2];
    idx[0] = riter.idx[0];
    gkyl_rect_grid_cell_center(&rzgrid, idx, xc);

    struct RdRdZ_sol sol = cub_root(psih, psi, Z, xc, dx);
    
    if (sol.nsol > 0)
      for (int s=0; s<sol.nsol; ++s) {
          R[sidx] = sol.R[s];
          dR[sidx] = sol.dRdZ[s];
          sidx += 1;
      }
  }
  return sidx;
}


void
compare_quad_and_cub(void)
{

  clock_t start, end;
  double cpu_time_used;

  struct gkyl_efit_inp inp  = {
    .filepath = "gyrokinetic/data/eqdsk/wham.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect =  true,
  };
  struct gkyl_efit* efit = gkyl_efit_new(&inp);

  // project the cubic on cubic basis: this should result in the same
  // DG expansions
  double lower[2] = {efit->rmin, efit->zmin };
  double upper[2] = {efit->rmax, efit->zmax};
  int cells[2] = {efit->nr-1, efit->nz-1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };  
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);
  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 2, 3);
  gkyl_proj_on_basis *projCub = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .num_ret_vals = 1,
      .ctx = efit->evf->ctx,
      .eval = efit->evf->eval_cubic
    }
  );
  struct gkyl_array *psi_cubic_DG = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis_advance(projCub, 0.0, &local, psi_cubic_DG);
  gkyl_grid_sub_array_write(&grid, &local, 0, psi_cubic_DG, "psi_cubic.gkyl");
  gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "psi_quad.gkyl");

  // Now pick a value of Z. Let's choose Z = 0.0
  // We want to see how long each one takes to find the roots
  double psi0 = 1e-10;
  double Z = 0.0;
  int nmaxroots = 4;
  double R[nmaxroots], dR[nmaxroots];
  int nr = 0;

  start = clock();
  nr = getR(efit->rzlocal, efit->rzgrid, efit->psizr, psi0, Z, nmaxroots, R, dR);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  //printf("Quadratic total time = %g\n", cpu_time_used);
  double quad_time = cpu_time_used;

  //printf("Quadratic: nr = %d, R = ", nr);
  //for(int i=0; i<nr; i++) printf(" %g ", R[i]);
  //printf(" | dRdZ = ");
  //for(int i=0; i<nr; i++) printf(" %g ", dR[i]);
  //printf("\n");

  for(int i = 0; i<nmaxroots; i++){
    R[i] = 0.0;
    dR[i] = 0.0;
  }


  start = clock();
  nr = getRcub(local, grid, psi_cubic_DG, psi0, Z, nmaxroots, R, dR);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  //printf("Cubic total time = %g\n", cpu_time_used);
  double cubic_time = cpu_time_used;

  //printf("Cubic: nr = %d, R = ", nr);
  //for(int i=0; i<nr; i++) printf(" %g ", R[i]);
  //printf(" | dRdZ = ");
  //for(int i=0; i<nr; i++) printf(" %g ", dR[i]);
  //printf("\n");

  double ratio = cubic_time/quad_time;
  //printf("Ratio of cubic to quadratic time = %g\n", ratio);

  gkyl_proj_on_basis_release(projCub);
  gkyl_array_release(psi_cubic_DG);
  gkyl_efit_release(efit);


}

TEST_LIST = {
  { "compare_quad_and_cub", compare_quad_and_cub },
  { NULL, NULL },
};
