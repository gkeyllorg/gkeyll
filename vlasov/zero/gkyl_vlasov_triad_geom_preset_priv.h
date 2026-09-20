#include <gkyl_vlasov_triad_geom.h>
#include <gkyl_vlasov_triad_geom_priv.h>

static void
eval_vierbein_inv_from_vierbein(int vdim, evalf_t eval_vierbein,
  double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double vierbein[GKYL_MAX_DIM*GKYL_MAX_DIM] = { 0.0 };
  vierbein_inv_t compute_vierbein_inv_from_vierbein = choose_vierbein_inv_kern(vdim);

  eval_vierbein(t, xn, vierbein, ctx);
  compute_vierbein_inv_from_vierbein(vierbein, fout);
}

void
eval_flat_vierbein_1v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbeins

  // e_ij = e_i(x) . \sigma_j(x)
  double e_xx = 1.0; // Vierbein Coefficients (x-x coefficient).

  fout[0] = e_xx;
}

void
eval_flat_vierbein_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbeins

  // e_ij = e_i(x) . \sigma_j(x)
  double e_xx = 1.0; // Vierbein Coefficients (x-x coefficient).
  double e_xy = 0.0; // Vierbein Coefficients (x-y coefficient).
  double e_yx = 0.0; // Vierbein Coefficients (y-x coefficient).
  double e_yy = 1.0; // Vierbein Coefficients (y-y coefficient).

  fout[0] = e_xx;
  fout[1] = e_xy;
  fout[2] = e_yx;
  fout[3] = e_yy;
}

void
eval_flat_vierbein_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbeins

  // e_ij = e_i(x) . \sigma_j(x)
  double e_xx = 1.0; // Vierbein Coefficients (x-x coefficient).
  double e_xy = 0.0; // Vierbein Coefficients (x-y coefficient).
  double e_xz = 0.0; // Vierbein Coefficients (x-z coefficient).

  double e_yx = 0.0; // Vierbein Coefficients (y-x coefficient).
  double e_yy = 1.0; // Vierbein Coefficients (y-y coefficient).
  double e_yz = 0.0; // Vierbein Coefficients (y-z coefficient).

  double e_zx = 0.0; // Vierbein Coefficients (z-x coefficient).
  double e_zy = 0.0; // Vierbein Coefficients (z-y coefficient).
  double e_zz = 1.0; // Vierbein Coefficients (z-z coefficient).

  fout[0] = e_xx;
  fout[1] = e_xy;
  fout[2] = e_xz;

  fout[3] = e_yx;
  fout[4] = e_yy;
  fout[5] = e_yz;

  fout[6] = e_zx;
  fout[7] = e_zy;
  fout[8] = e_zz;
}

void
eval_spherical_rtheta_vierbein_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbeins
  double r = xn[0], theta = xn[1];

  // e_ij = e_i(x) . \sigma_j(x)
  double e_rr = 1.0; // Vierbein Coefficients (r-r coefficient).
  double e_rt = 0.0; // Vierbein Coefficients (r-theta coefficient).
  double e_rp = 0.0; // Vierbein Coefficients (r-phi coefficient).

  double e_tr = 0.0; // Vierbein Coefficients (theta-r coefficient).
  double e_tt = r; // Vierbein Coefficients (theta-theta coefficient).
  double e_tp = 0.0; // Vierbein Coefficients (theta-phi coefficient).

  double e_pr = 0.0; // Vierbein Coefficients (phi-r coefficient).
  double e_pt = 0.0; // Vierbein Coefficients (phi-theta coefficient).
  double e_pp = r*sin(theta); // Vierbein Coefficients (phi-phi coefficient).

  fout[0] = e_rr;
  fout[1] = e_rt;
  fout[2] = e_rp;

  fout[3] = e_tr;
  fout[4] = e_tt;
  fout[5] = e_tp;

  fout[6] = e_pr;
  fout[7] = e_pt;
  fout[8] = e_pp;
}

void
eval_flat_vierbein_inv_1v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  eval_vierbein_inv_from_vierbein(1, eval_flat_vierbein_1v, t, xn, fout, ctx);
}

void
eval_flat_vierbein_inv_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  eval_vierbein_inv_from_vierbein(2, eval_flat_vierbein_2v, t, xn, fout, ctx);
}

void
eval_flat_vierbein_inv_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  eval_vierbein_inv_from_vierbein(3, eval_flat_vierbein_3v, t, xn, fout, ctx);
}

void
eval_spherical_rtheta_vierbein_inv_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  eval_vierbein_inv_from_vierbein(3, eval_spherical_rtheta_vierbein_3v, t, xn, fout, ctx);
}

void
eval_flat_vierbein_gradient_1v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = 0.0;
}

void
eval_flat_vierbein_gradient_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  for (int i=0; i<8; ++i) fout[i] = 0.0;
}

void
eval_flat_vierbein_gradient_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  for (int i=0; i<27; ++i) fout[i] = 0.0;
}

void
eval_spherical_rtheta_vierbein_gradient_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbein Gradients
  double r = xn[0], theta = xn[1];

  // d (e_ij) / dr
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;

  fout[3] = 0.0;
  fout[4] = 1.0;
  fout[5] = 0.0;

  fout[6] = 0.0;
  fout[7] = 0.0;
  fout[8] = sin(theta);

  // d (e_ij) / dtheta
  fout[9] = 0.0;
  fout[10] = 0.0;
  fout[11] = 0.0;

  fout[12] = 0.0;
  fout[13] = 0.0;
  fout[14] = 0.0;

  fout[15] = 0.0;
  fout[16] = 0.0;
  fout[17] = r*cos(theta);

  // d (e_ij) / dphi
  fout[18] = 0.0;
  fout[19] = 0.0;
  fout[20] = 0.0;

  fout[21] = 0.0;
  fout[22] = 0.0;
  fout[23] = 0.0;

  fout[24] = 0.0;
  fout[25] = 0.0;
  fout[26] = 0.0;
}

void
eval_annulus_vierbein_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbeins
  double q_r = xn[0];

  // e_ij = e_i(x) . \sigma_j(x)
  double e_rr = 1.0; // Vierbein Coefficients (r-r coefficient).
  double e_rt = 0.0; // Vierbein Coefficients (r-theta coefficient).
  double e_tr = 0.0; // Vierbein Coefficients (theta-r coefficient).
  double e_tt = q_r; // Vierbein Coefficients (theta-theta coefficient).
  
  fout[0] = e_rr;
  fout[1] = e_rt;
  fout[2] = e_tr;
  fout[3] = e_tt;
}

void
eval_annulus_vierbein_inv_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  eval_vierbein_inv_from_vierbein(2, eval_annulus_vierbein_2v, t, xn, fout, ctx);
}

void
eval_annulus_vierbein_gradient_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbein Gradients 
  double q_r = xn[0];

  // d (e_ij) / dr 
  double d_e_rr_dr = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  double d_e_rt_dr = 0.0; // Vierbein Gradient Coefficients (r-theta coefficient).
  double d_e_tr_dr = 0.0; // Vierbein Gradient Coefficients (theta-r coefficient).
  double d_e_tt_dr = 1.0; // Vierbein Gradient Coefficients (theta-theta coefficient).

  // d (e_ij) / dtheta 
  double d_e_rr_dt = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  double d_e_rt_dt = 0.0; // Vierbein Gradient Coefficients (r-theta coefficient).
  double d_e_tr_dt = 0.0; // Vierbein Gradient Coefficients (theta-r coefficient).
  double d_e_tt_dt = 0.0; // Vierbein Gradient Coefficients (theta-theta coefficient).
  
  fout[0] = d_e_rr_dr;
  fout[1] = d_e_rt_dr;
  fout[2] = d_e_tr_dr;
  fout[3] = d_e_tt_dr;

  fout[4] = d_e_rr_dt;
  fout[5] = d_e_rt_dt;
  fout[6] = d_e_tr_dt;
  fout[7] = d_e_tt_dt;

}

static inline double
flat_hamil(const double* GKYL_RESTRICT xn, int cdim, int vdim)
{
  double p_sq = 0.0;
  for (int d=0; d<vdim; ++d) {
    double p_hat = xn[cdim+d];
    p_sq += p_hat*p_hat;
  }
  return sqrt(1.0 + p_sq);
}

void
eval_flat_hamil_1x1v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = flat_hamil(xn, 1, 1);
}

void
eval_flat_hamil_2x1v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = flat_hamil(xn, 2, 1);
}

void
eval_flat_hamil_3x1v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = flat_hamil(xn, 3, 1);
}

void
eval_flat_hamil_1x2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = flat_hamil(xn, 1, 2);
}

void
eval_flat_hamil_2x2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = flat_hamil(xn, 2, 2);
}

void
eval_flat_hamil_3x2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = flat_hamil(xn, 3, 2);
}

void
eval_flat_hamil_1x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = flat_hamil(xn, 1, 3);
}

void
eval_flat_hamil_2x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = flat_hamil(xn, 2, 3);
}

void
eval_flat_hamil_3x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = flat_hamil(xn, 3, 3);
}

void
eval_ks_rphi_hamil_1x2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];

  double pr_hat = xn[1];
  double pphi_hat = xn[2];

  // Build the Hamiltonian for Kerr-Schild in 2D (Assumes theta = pi/2)
  double rho_sq = r * r ;
  double H = (1.0/ sqrt(1 + 2*M*r/rho_sq)) * (
        sqrt(1 + pr_hat * pr_hat + pphi_hat * pphi_hat) 
        - (2*M*r/rho_sq) * pr_hat );

  fout[0] = H;
}

void
eval_ks_rphi_hamil_2x2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];
  double phi = xn[1];

  double pr_hat = xn[2];
  double pphi_hat = xn[3];

  // Build the Hamiltonian for Kerr-Schild in 2D (Assumes theta = pi/2)
  double rho_sq = r * r ;
  double H = (1.0/ sqrt(1 + 2*M*r/rho_sq)) * (
        sqrt(1 + pr_hat * pr_hat + pphi_hat * pphi_hat) 
        - (2*M*r/rho_sq) * pr_hat );

  fout[0] = H;
}

void
eval_ks_r_hamil_1x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];

  double pr_hat = xn[1];
  double ptheta_hat = xn[2];
  double pphi_hat = xn[3];

  // Build the Hamiltonian for Kerr-Schild in 3D
  double rho_sq = r * r;
  double H = (1.0/ sqrt(1 + 2*M*r/rho_sq)) * (
        sqrt(1 + pr_hat * pr_hat + ptheta_hat * ptheta_hat + pphi_hat * pphi_hat) 
        - (2*M*r/rho_sq) * pr_hat );

  fout[0] = H;
}


void
eval_ks_rtheta_hamil_2x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];
  double theta = xn[1];

  double pr_hat = xn[2];
  double ptheta_hat = xn[3];
  double pphi_hat = xn[4];

  // Build the Hamiltonian for Kerr-Schild in 3D
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);
  double H = (1.0/ sqrt(1 + 2*M*r/rho_sq)) * (
        sqrt(1 + pr_hat * pr_hat + ptheta_hat * ptheta_hat + pphi_hat * pphi_hat) 
        - (2*M*r/rho_sq) * pr_hat );

  fout[0] = H;
}

void
eval_ks_rphi_vierbein_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Parameters
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  // Coordinates
  double r = xn[0];
  double phi = xn[1];

  // Intermediate Variables
  double rho_sq = r * r ;

  // Metric spatial covariant components
  double h_rr = ( 1.0 + 2.0 * M * r / rho_sq);
  double h_pp = ( rho_sq + a * a * ( 1.0 + 2.0 * M * r / rho_sq )   );
  double h_thth = rho_sq;
  double h_rp = - a * ( 1.0 + 2.0 * M * r / rho_sq );

  // Vierbein: e_i^a = g_i . sigma^a
  fout[0] = sqrt( h_rr );
  fout[1] = 0.0;
  fout[2] = h_rp / sqrt( h_rr );
  fout[3] = sqrt( h_pp - h_rp * h_rp / h_rr );
}

void
eval_ks_rphi_vierbein_gradient_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  double r = xn[0];
  double phi = xn[1];

  // Intermediate Variables
  double rho_sq = r * r ;

  // Gradient w.r.t. r: d(e_i^a)/dr
  fout[0] = (M * ( - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[1] = 0.0;
  fout[2] = -(M * a * ( - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[3] = r / sqrt(rho_sq);

  // Gradient w.r.t. phi: d(e_i^a)/dphi
  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;

}

void 
eval_ks_r_vierbein_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  // Parameters
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  // Coordinates
  double r = xn[0];

  // Intermediate Variables
  double rho_sq = r * r;

  // Metric spatial covariant components
  double h_rr = ( 1.0 + 2.0 * M * r / rho_sq);
  double h_pp = ( rho_sq + a * a * ( 1.0 + 2.0 * M * r / rho_sq )  );
  double h_thth = rho_sq;
  double h_rp = - a * ( 1.0 + 2.0 * M * r / rho_sq );

  // Vierbein: e_i^a = g_i . sigma^a
  fout[0] = sqrt( h_rr );
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = sqrt( rho_sq );
  fout[5] = 0.0;
  fout[6] = h_rp / sqrt( h_rr );
  fout[7] = 0.0;
  fout[8] = sqrt( h_pp - h_rp * h_rp / h_rr );
}

void 
eval_ks_r_vierbein_gradient_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  double r = xn[0];

  // Intermediate Variables
  double rho_sq = r * r ;

  // Gradient w.r.t. r: d(e_i^a)/dr
  fout[0] = (M * ( - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = r / sqrt(rho_sq);
  fout[5] = 0.0;
  fout[6] = -(M * a * ( - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[7] = 0.0;
  fout[8] = (r) / sqrt(rho_sq);

  // Gradient w.r.t. theta: d(e_i^a)/dtheta
  fout[9]  = 0.0;
  fout[10] = 0.0;
  fout[11] = 0.0;
  fout[12] = 0.0;
  fout[13] = 0.0;
  fout[14] = 0.0;
  fout[15] = 0.0;
  fout[16] = 0.0;
  fout[17] = 0.0;

  // Gradient w.r.t. theta: d(e_i^a)/dphi
  fout[18] = 0.0;
  fout[19] = 0.0;
  fout[20] = 0.0;
  fout[21] = 0.0;
  fout[22] = 0.0;
  fout[23] = 0.0;
  fout[24] = 0.0;
  fout[25] = 0.0;
  fout[26] = 0.0;
}

void 
eval_ks_rtheta_vierbein_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  // Parameters
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  // Coordinates
  double r = xn[0];
  double theta = xn[1];

  // Intermediate Variables
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);

  // Metric spatial covariant components
  double h_rr = ( 1.0 + 2.0 * M * r / rho_sq);
  double h_pp = sin(theta) * sin(theta) * ( rho_sq + a * a * ( 1.0 + 2.0 * M * r / rho_sq ) * sin(theta) * sin(theta)  );
  double h_thth = rho_sq;
  double h_rp = - a * ( 1.0 + 2.0 * M * r / rho_sq ) * sin(theta) * sin(theta);

  // Vierbein: e_i^a = g_i . sigma^a
  fout[0] = sqrt( h_rr );
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = sqrt( rho_sq );
  fout[5] = 0.0;
  fout[6] = h_rp / sqrt( h_rr );
  fout[7] = 0.0;
  fout[8] = sqrt( h_pp - h_rp * h_rp / h_rr );
}

void 
eval_ks_rtheta_vierbein_gradient_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  double r = xn[0];
  double theta = xn[1];

  // Intermediate Variables
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);

  // Gradient w.r.t. r: d(e_i^a)/dr
  fout[0] = (M * (a * a * cos(theta) * cos(theta) - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = r / sqrt(rho_sq);
  fout[5] = 0.0;
  fout[6] = -(M * a * sin(theta) * sin(theta) * (a * a * cos(theta) * cos(theta) - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[7] = 0.0;
  fout[8] = (r * sin(theta)) / sqrt(rho_sq);

  // Gradient w.r.t. theta: d(e_i^a)/dtheta
  fout[9]  = (M * a * a * r * sin(2.0 * theta)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[10] = 0.0;
  fout[11] = 0.0;
  fout[12] = 0.0;
  fout[13] = -(a * a * sin(2.0 * theta)) / (2.0 * sqrt(rho_sq));
  fout[14] = 0.0;
  fout[15] = -(2.0 * a * cos(theta) * sin(theta) * (a * a * a * a * cos(theta) * cos(theta) * cos(theta) * cos(theta) + 2.0 * M * r * r * r + r * r * r * r + 2.0 * a * a * r * r * cos(theta) * cos(theta) + 2.0 * M * a * a * r * cos(theta) * cos(theta) + M * a * a * r * sin(theta) * sin(theta))) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[16] = 0.0;
  fout[17] = (sin(4.0 * theta) * a * a + 2.0 * sin(2.0 * theta) * r * r) / (2.0 * sqrt(1.0 - cos(2.0 * theta)) * sqrt(a * a * cos(2.0 * theta) + a * a + 2.0 * r * r));

  // Gradient w.r.t. theta: d(e_i^a)/dphi
  fout[18] = 0.0;
  fout[19] = 0.0;
  fout[20] = 0.0;
  fout[21] = 0.0;
  fout[22] = 0.0;
  fout[23] = 0.0;
  fout[24] = 0.0;
  fout[25] = 0.0;
  fout[26] = 0.0;
}

void 
eval_ks_vierbein_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  // Parameters
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  // Coordinates
  double r = xn[0];
  double theta = xn[1];
  double phi = xn[2];

  // Intermediate Variables
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);

  // Metric spatial covariant components
  double h_rr = ( 1.0 + 2.0 * M * r / rho_sq);
  double h_pp = sin(theta) * sin(theta) * ( rho_sq + a * a * ( 1.0 + 2.0 * M * r / rho_sq ) * sin(theta) * sin(theta)  );
  double h_thth = rho_sq;
  double h_rp = - a * ( 1.0 + 2.0 * M * r / rho_sq ) * sin(theta) * sin(theta);

  // Vierbein: e_i^a = g_i . sigma^a
  fout[0] = sqrt( h_rr );
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = sqrt( rho_sq );
  fout[5] = 0.0;
  fout[6] = h_rp / sqrt( h_rr );
  fout[7] = 0.0;
  fout[8] = sqrt( h_pp - h_rp * h_rp / h_rr );
}

void 
eval_ks_vierbein_cart_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  // Parameters
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  // Coordinates
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];

  // Intermediate Variables

  // Metric spatial covariant components

  // Vierbein: e_i^a = g_i . sigma^a
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
  fout[8] = 0.0;
}

void 
eval_ks_vierbein_gradient_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  double r = xn[0];
  double theta = xn[1];
  double phi = xn[2];

  // Intermediate Variables
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);

  // Gradient w.r.t. r: d(e_i^a)/dr
  fout[0] = (M * (a * a * cos(theta) * cos(theta) - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = r / sqrt(rho_sq);
  fout[5] = 0.0;
  fout[6] = -(M * a * sin(theta) * sin(theta) * (a * a * cos(theta) * cos(theta) - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[7] = 0.0;
  fout[8] = (r * sin(theta)) / sqrt(rho_sq);

  // Gradient w.r.t. theta: d(e_i^a)/dtheta
  fout[9]  = (M * a * a * r * sin(2.0 * theta)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[10] = 0.0;
  fout[11] = 0.0;
  fout[12] = 0.0;
  fout[13] = -(a * a * sin(2.0 * theta)) / (2.0 * sqrt(rho_sq));
  fout[14] = 0.0;
  fout[15] = -(2.0 * a * cos(theta) * sin(theta) * (a * a * a * a * cos(theta) * cos(theta) * cos(theta) * cos(theta) + 2.0 * M * r * r * r + r * r * r * r + 2.0 * a * a * r * r * cos(theta) * cos(theta) + 2.0 * M * a * a * r * cos(theta) * cos(theta) + M * a * a * r * sin(theta) * sin(theta))) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[16] = 0.0;
  fout[17] = (sin(4.0 * theta) * a * a + 2.0 * sin(2.0 * theta) * r * r) / (2.0 * sqrt(1.0 - cos(2.0 * theta)) * sqrt(a * a * cos(2.0 * theta) + a * a + 2.0 * r * r));

  // Gradient w.r.t. theta: d(e_i^a)/dphi
  fout[18] = 0.0;
  fout[19] = 0.0;
  fout[20] = 0.0;
  fout[21] = 0.0;
  fout[22] = 0.0;
  fout[23] = 0.0;
  fout[24] = 0.0;
  fout[25] = 0.0;
  fout[26] = 0.0;
}


void 
eval_ks_vierbein_gradient_cart_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  double x = xn[0];
  double y = xn[1];
  double z = xn[2];

  // Intermediate Variables

  // Gradient w.r.t. r: d(e_i^a)/dr
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
  fout[8] = 0.0;

  // Gradient w.r.t. theta: d(e_i^a)/dtheta
  fout[9]  = 0.0;
  fout[10] = 0.0;
  fout[11] = 0.0;
  fout[12] = 0.0;
  fout[13] = 0.0;
  fout[14] = 0.0;
  fout[15] = 0.0;
  fout[16] = 0.0;
  fout[17] = 0.0;

  // Gradient w.r.t. theta: d(e_i^a)/dphi
  fout[18] = 0.0;
  fout[19] = 0.0;
  fout[20] = 0.0;
  fout[21] = 0.0;
  fout[22] = 0.0;
  fout[23] = 0.0;
  fout[24] = 0.0;
  fout[25] = 0.0;
  fout[26] = 0.0;
}


void
eval_rz_cylindrical_vierbein_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbeins
  double q_r = xn[0];

  // e_ij = e_i(x) . \sigma_j(x)
  double e_rr = 1.0; // Vierbein Coefficients (r-r coefficient).
  double e_rz = 0.0;
  double e_rt = 0.0; // Vierbein Coefficients (r-theta coefficient).
  
  double e_zr = 0.0;
  double e_zz = 1.0;
  double e_zt = 0.0; 

  double e_tr = 0.0; // Vierbein Coefficients (theta-r coefficient).
  double e_tz = 0.0; 
  double e_tt = q_r; // Vierbein Coefficients (theta-theta coefficient).
  
  fout[0] = e_rr;
  fout[1] = e_rz;
  fout[2] = e_rt;

  fout[3] = e_zr;
  fout[4] = e_zz;
  fout[5] = e_zt;

  fout[6] = e_tr;
  fout[7] = e_tz;
  fout[8] = e_tt;
}

void
eval_rz_cylindrical_vierbein_gradient_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbein Gradients 
  double q_r = xn[0];

  // d (e_ij) / dr 
  double d_e_rr_dr = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  double d_e_rz_dr = 0.0;
  double d_e_rt_dr = 0.0; // Vierbein Gradient Coefficients (r-theta coefficient).
  
  double d_e_zr_dr = 0.0;
  double d_e_zz_dr = 0.0;
  double d_e_zt_dr = 0.0;

  double d_e_tr_dr = 0.0; // Vierbein Gradient Coefficients (theta-r coefficient).
  double d_e_tz_dr = 0.0;
  double d_e_tt_dr = 1.0; // Vierbein Gradient Coefficients (theta-theta coefficient).

  // d(e_ij) / dz
  double d_e_rr_dz = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  double d_e_rz_dz = 0.0;
  double d_e_rt_dz = 0.0; // Vierbein Gradient Coefficients (r-theta coefficient).
  
  double d_e_zr_dz = 0.0;
  double d_e_zz_dz = 0.0;
  double d_e_zt_dz = 0.0;

  double d_e_tr_dz = 0.0; // Vierbein Gradient Coefficients (theta-r coefficient).
  double d_e_tz_dz = 0.0;
  double d_e_tt_dz = 0.0; // Vierbein Gradient Coefficients (theta-theta coefficient).


  // d (e_ij) / dtheta 
  double d_e_rr_dt = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  double d_e_rz_dt = 0.0;
  double d_e_rt_dt = 0.0; // Vierbein Gradient Coefficients (r-theta coefficient).

  double d_e_zr_dt = 0.0;
  double d_e_zz_dt = 0.0;
  double d_e_zt_dt = 0.0;

  double d_e_tr_dt = 0.0; // Vierbein Gradient Coefficients (theta-r coefficient).
  double d_e_tz_dt = 0.0;
  double d_e_tt_dt = 0.0; // Vierbein Gradient Coefficients (theta-theta coefficient).
  
  fout[0] = d_e_rr_dr;
  fout[1] = d_e_rz_dr;
  fout[2] = d_e_rt_dr;
  fout[3] = d_e_zr_dr;
  fout[4] = d_e_zz_dr;
  fout[5] = d_e_zt_dr;
  fout[6] = d_e_tr_dr;
  fout[7] = d_e_tz_dr;
  fout[8] = d_e_tt_dr;

  fout[9+0] = d_e_rr_dz;
  fout[9+1] = d_e_rz_dz;
  fout[9+2] = d_e_rt_dz;
  fout[9+3] = d_e_zr_dz;
  fout[9+4] = d_e_zz_dz;
  fout[9+5] = d_e_zt_dz;
  fout[9+6] = d_e_tr_dz;
  fout[9+7] = d_e_tz_dz;
  fout[9+8] = d_e_tt_dz;

  fout[18+0] = d_e_rr_dt;
  fout[18+1] = d_e_rz_dt;
  fout[18+2] = d_e_rt_dt;
  fout[18+3] = d_e_zr_dt;
  fout[18+4] = d_e_zz_dt;
  fout[18+5] = d_e_zt_dt;
  fout[18+6] = d_e_tr_dt;
  fout[18+7] = d_e_tz_dt;
  fout[18+8] = d_e_tt_dt;

}

void
eval_ks_hamil_3x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];
  double theta = xn[1];
  double phi = xn[2];

  double pr_hat = xn[3];
  double ptheta_hat = xn[4];
  double pphi_hat = xn[5];

  // Build the Hamiltonian for Kerr-Schild in 3D
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);
  double H = (1.0/ sqrt(1 + 2*M*r/rho_sq)) * (
        sqrt(1 + pr_hat * pr_hat + ptheta_hat * ptheta_hat + pphi_hat * pphi_hat) 
        - (2*M*r/rho_sq) * pr_hat );

  fout[0] = H;
}

void
eval_ks_hamil_cart_3x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];

  double px_hat = xn[3];
  double py_hat = xn[4];
  double pz_hat = xn[5];

  // Build the Hamiltonian for Kerr-Schild in 3D
  double H = 0.0;

  fout[0] = H;
}

void
eval_rz_cylindrical_vierbein_inv_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  eval_vierbein_inv_from_vierbein(3, eval_rz_cylindrical_vierbein_3v, t, xn, fout, ctx);
}

void
eval_ks_rphi_vierbein_inv_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  eval_vierbein_inv_from_vierbein(2, eval_ks_rphi_vierbein_2v, t, xn, fout, ctx);
}

void
eval_ks_r_vierbein_inv_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  eval_vierbein_inv_from_vierbein(3, eval_ks_r_vierbein_3v, t, xn, fout, ctx);
}

void
eval_ks_rtheta_vierbein_inv_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  eval_vierbein_inv_from_vierbein(3, eval_ks_rtheta_vierbein_3v, t, xn, fout, ctx);
}

void
eval_ks_vierbein_inv_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  eval_vierbein_inv_from_vierbein(3, eval_ks_vierbein_3v, t, xn, fout, ctx);
}

void
eval_ks_vierbein_inv_cart_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  eval_vierbein_inv_from_vierbein(3, eval_ks_vierbein_cart_3v, t, xn, fout, ctx);
}

// for use in kernel tables
typedef struct { evalf_t kernels[3]; } hamil_kern_list;
typedef struct { evalf_t kernels[3]; } vierbein_kern_list;
typedef struct { evalf_t kernels[3]; } vierbein_gradient_kern_list;

static const hamil_kern_list flat_hamil_list[] = {
  { eval_flat_hamil_1x1v, eval_flat_hamil_2x1v, eval_flat_hamil_3x1v },
  { eval_flat_hamil_1x2v, eval_flat_hamil_2x2v, eval_flat_hamil_3x2v },
  { eval_flat_hamil_1x3v, eval_flat_hamil_2x3v, eval_flat_hamil_3x3v }
};

static const hamil_kern_list ks_rphi_hamil_list[] = {
  { NULL, NULL, NULL },
  { eval_ks_rphi_hamil_1x2v, eval_ks_rphi_hamil_2x2v, NULL },
  { NULL, NULL, NULL }
};

static const hamil_kern_list ks_r_hamil_list[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { eval_ks_r_hamil_1x3v, NULL, NULL }
};

static const hamil_kern_list ks_rtheta_hamil_list[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, eval_ks_rtheta_hamil_2x3v, eval_ks_hamil_3x3v }
};

static const hamil_kern_list ks_hamil_3v_list[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, eval_ks_hamil_3x3v }
};

static const hamil_kern_list ks_hamil_cart_3v_list[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, eval_ks_hamil_cart_3x3v }
};

static const vierbein_kern_list flat_vierbein_list[] = {
  { eval_flat_vierbein_1v },
  { eval_flat_vierbein_2v },
  { eval_flat_vierbein_3v }
};

static const vierbein_kern_list flat_vierbein_inv_list[] = {
  { eval_flat_vierbein_inv_1v },
  { eval_flat_vierbein_inv_2v },
  { eval_flat_vierbein_inv_3v }
};

static const vierbein_gradient_kern_list flat_vierbein_gradient_list[] = {
  { eval_flat_vierbein_gradient_1v },
  { eval_flat_vierbein_gradient_2v },
  { eval_flat_vierbein_gradient_3v }
};

static const vierbein_kern_list annulus_vierbein_list[] = {
  { NULL },
  { eval_annulus_vierbein_2v },
  { NULL }
};

static const vierbein_kern_list annulus_vierbein_inv_list[] = {
  { NULL },
  { eval_annulus_vierbein_inv_2v },
  { NULL }
};

static const vierbein_gradient_kern_list annulus_vierbein_gradient_list[] = {
  { NULL },
  { eval_annulus_vierbein_gradient_2v },
  { NULL }
};

static const vierbein_kern_list ks_rphi_vierbein_list[] = {
  { NULL },
  { eval_ks_rphi_vierbein_2v },
  { NULL }
};

static const vierbein_kern_list ks_rphi_vierbein_inv_list[] = {
  { NULL },
  { eval_ks_rphi_vierbein_inv_2v },
  { NULL }
};

static const vierbein_gradient_kern_list ks_rphi_vierbein_gradient_list[] = {
  { NULL },
  { eval_ks_rphi_vierbein_gradient_2v },
  { NULL }
};

static const vierbein_kern_list ks_r_vierbein_list[] = {
  { NULL },
  { NULL },
  { eval_ks_r_vierbein_3v }
};

static const vierbein_kern_list ks_r_vierbein_inv_list[] = {
  { NULL },
  { NULL },
  { eval_ks_r_vierbein_inv_3v }
};

static const vierbein_gradient_kern_list ks_r_vierbein_gradient_list[] = {
  { NULL },
  { NULL },
  { eval_ks_r_vierbein_gradient_3v }
};

static const vierbein_kern_list ks_rtheta_vierbein_list[] = {
  { NULL },
  { NULL },
  { eval_ks_rtheta_vierbein_3v }
};

static const vierbein_kern_list ks_rtheta_vierbein_inv_list[] = {
  { NULL },
  { NULL },
  { eval_ks_rtheta_vierbein_inv_3v }
};

static const vierbein_gradient_kern_list ks_rtheta_vierbein_gradient_list[] = {
  { NULL },
  { NULL },
  { eval_ks_rtheta_vierbein_gradient_3v }
};

static const vierbein_kern_list ks_vierbein_list[] = {
  { NULL },
  { NULL},
  { eval_ks_vierbein_3v  }
};

static const vierbein_kern_list ks_vierbein_inv_list[] = {
  { NULL },
  { NULL },
  { eval_ks_vierbein_inv_3v }
};

static const vierbein_kern_list ks_vierbein_cart_list[] = {
  { NULL },
  { NULL},
  { eval_ks_vierbein_cart_3v  }
};

static const vierbein_kern_list ks_vierbein_inv_cart_list[] = {
  { NULL },
  { NULL },
  { eval_ks_vierbein_inv_cart_3v }
};

static const vierbein_gradient_kern_list ks_vierbein_gradient_list[] = {
  { NULL },
  { NULL },
  { eval_ks_vierbein_gradient_3v }
};

static const vierbein_gradient_kern_list ks_vierbein_gradient_cart_list[] = {
  { NULL },
  { NULL },
  { eval_ks_vierbein_gradient_cart_3v }
};

static const vierbein_kern_list rz_cylindrical_vierbein_list[] = {
  { NULL },
  { NULL },
  { eval_rz_cylindrical_vierbein_3v }
};

static const vierbein_kern_list rz_cylindrical_vierbein_inv_list[] = {
  { NULL },
  { NULL },
  { eval_rz_cylindrical_vierbein_inv_3v }
};

static const vierbein_gradient_kern_list rz_cylindrical_vierbein_gradient_list[] = {
  { NULL },
  { NULL },
  { eval_rz_cylindrical_vierbein_gradient_3v }
};

static const vierbein_kern_list spherical_rtheta_vierbein_list[] = {
  { NULL },
  { NULL },
  { eval_spherical_rtheta_vierbein_3v }
};

static const vierbein_kern_list spherical_rtheta_vierbein_inv_list[] = {
  { NULL },
  { NULL },
  { eval_spherical_rtheta_vierbein_inv_3v }
};

static const vierbein_gradient_kern_list spherical_rtheta_vierbein_gradient_list[] = {
  { NULL },
  { NULL },
  { eval_spherical_rtheta_vierbein_gradient_3v }
};

static evalf_t
choose_vierbein_kern(enum gkyl_triad_preset_geom_type type, int vdim)
{
  switch(type) {
    case GKYL_TRIAD_FLAT:
      return flat_vierbein_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_ANNULUS:
      return annulus_vierbein_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_CYLINDRICAL_RZ:
      return rz_cylindrical_vierbein_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_SPHERICAL_RTHETA:
      return spherical_rtheta_vierbein_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_RPHI:
      return ks_rphi_vierbein_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_R:
      return ks_r_vierbein_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_RTHETA:
      return ks_rtheta_vierbein_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return ks_vierbein_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_CART_GR_KERR_SCHILD_3V:
      return ks_vierbein_cart_list[vdim-1].kernels[0];
      break;
    default:
      assert(false);
  }
}

static evalf_t
choose_vierbein_inv_preset_kern(enum gkyl_triad_preset_geom_type type, int vdim)
{
  switch(type) {
    case GKYL_TRIAD_FLAT:
      return flat_vierbein_inv_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_ANNULUS:
      return annulus_vierbein_inv_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_CYLINDRICAL_RZ:
      return rz_cylindrical_vierbein_inv_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_SPHERICAL_RTHETA:
      return spherical_rtheta_vierbein_inv_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_RPHI:
      return ks_rphi_vierbein_inv_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_R:
      return ks_r_vierbein_inv_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_RTHETA:
      return ks_rtheta_vierbein_inv_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return ks_vierbein_inv_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_CART_GR_KERR_SCHILD_3V:
      return ks_vierbein_inv_cart_list[vdim-1].kernels[0];
      break;
    default:
      assert(false);
  }
}

static evalf_t
choose_vierbein_gradient_kern(enum gkyl_triad_preset_geom_type type, int vdim)
{
  switch(type) {
    case GKYL_TRIAD_FLAT:
      return flat_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_ANNULUS:
      return annulus_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_CYLINDRICAL_RZ:
      return rz_cylindrical_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_SPHERICAL_RTHETA:
      return spherical_rtheta_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_RPHI:
      return ks_rphi_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_R:
      return ks_r_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_RTHETA:
      return ks_rtheta_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return ks_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_CART_GR_KERR_SCHILD_3V:
      return ks_vierbein_gradient_cart_list[vdim-1].kernels[0];
      break;
    default:
      assert(false);
  }
}

static evalf_t
choose_hamil_kern(enum gkyl_triad_preset_geom_type type, int cdim, int vdim)
{
  switch(type) {
    case GKYL_TRIAD_FLAT:
      return flat_hamil_list[vdim-1].kernels[cdim-1];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_RPHI:
      return ks_rphi_hamil_list[vdim-1].kernels[cdim-1];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_R:
      return ks_r_hamil_list[vdim-1].kernels[cdim-1];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_RTHETA:
      return ks_rtheta_hamil_list[vdim-1].kernels[cdim-1];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return ks_hamil_3v_list[vdim-1].kernels[cdim-1];
      break;
    case GKYL_TRIAD_CART_GR_KERR_SCHILD_3V:
      return ks_hamil_cart_3v_list[vdim-1].kernels[cdim-1];
      break;
    default:
      assert(false);
  }
}
