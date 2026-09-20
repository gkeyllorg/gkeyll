#include <gkyl_dg_gr_maxwell_geom.h>
#include <gkyl_vlasov_triad_geom.h>
#include <assert.h>

void
eval_non_relativistic_lapse(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Flat Geometry, so lapse is 1.
  fout[0] = 1.0;
}

void
eval_non_relativistic_shift(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Flat Geometry, so shift is 0.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
eval_flat_h_ij(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // theta = pi/2
  double r = xn[0];

  // Metric spatial covariant components h_ij
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 1.0;
  fout[4] = 0.0;
  fout[5] = 1.0;
}

void
eval_flat_h_ij_inv(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  
  // theta = pi/2
  double r = xn[0];

  // Metric spatialcontravariant components h_ij_inv
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 1.0;
  fout[4] = 0.0;
  fout[5] = 1.0;
}

void
eval_geom_factor_flat(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // no contributions in cartesian coordinates
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
eval_flat_det_h(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // determinant is simply 1
  fout[0] = 1.0;

}

void
eval_annulus_h_ij(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // theta = pi/2
  double r = xn[0];

  // Metric spatial covariant components  h_ij
  double h_rr = 1.0;
  double h_pp = r * r;
  double h_zz = 1.0;

  fout[0] = h_rr; // h_rr
  fout[1] = 0.0; // h_rp
  fout[2] = 0.0; // h_rt
  fout[3] = h_pp; // h_pp
  fout[4] = 0.0; // h_tp
  fout[5] = h_zz; // h_zz
}


void
eval_annulus_h_ij_inv(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  
  // theta = pi/2
  double r = xn[0];

  // Metric spatialcontravariant components  h_ij_inv
  double h_rr = 1.0;
  double h_pp = r * r;
  double h_zz = 1.0;

  fout[0] = 1/h_rr; // h_rr
  fout[1] = 0.0; // h_rp
  fout[2] = 0.0; // h_rt
  fout[3] = 1/h_pp; // h_pp
  fout[4] = 0.0; // h_tp
  fout[5] = 1/h_zz; // h_zz
}

void
eval_geom_factor_annulus(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  
  // Only radial contributions
  double r = xn[0];

  // no contributions in cartesian coordinates
  fout[0] = 1.0/r;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
eval_annulus_det_h(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // determinant is simply r
  double r = xn[0];

  fout[0] = r;
}

void
eval_spherical_h_ij_rtheta(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // theta = pi/2
  double r = xn[0];
  double theta = xn[1];

  // Metric spatial covariant components  h_ij
  double h_rr = 1.0;
  double h_tt = r * r;
  double h_pp = r * r * sin(theta) * sin(theta);

  fout[0] = h_rr; // h_rr
  fout[1] = 0.0; // h_rt
  fout[2] = 0.0; // h_rp
  fout[3] = h_tt; // h_tt
  fout[4] = 0.0; // h_tp
  fout[5] = h_pp; // h_pp
}

void
eval_spherical_h_ij_inv_rtheta(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  
  // Coordrinates are r, theta which the problem depends on
  double r = xn[0];
  double theta = xn[1];

  // Metric spatialcontravariant components  h_ij_inv
  double h_rr = 1.0;
  double h_tt = r * r;
  double h_pp = r * r * sin(theta) * sin(theta);

  fout[0] = 1/h_rr; // h_rr
  fout[1] = 0.0; // h_rt
  fout[2] = 0.0; // h_rp
  fout[3] = 1/h_tt; // h_tt
  fout[4] = 0.0; // h_tp
  fout[5] = 1/h_pp; // h_pp
}

void
eval_geom_factor_sph_rtheta(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  
  // Only radial contributions
  double r = xn[0];
  double theta = xn[1];

  // no contributions in cartesian coordinates
  fout[0] = 2.0/r;
  fout[1] = cos(theta)/(r*r*sin(theta));
  fout[2] = 0.0;
}

void
eval_spherical_det_h_rtheta(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  
  // determinant is simply r*r sin(theta)
  double r = xn[0];
  double theta = xn[1];

  fout[0] = r * r * sin(theta);
}

void
eval_ks_lapse_r(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Grab the geometry
  const struct gkyl_dg_gr_maxwell_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  // theta = pi/2
  double r = xn[0];

  // Build the lapse
  double rho_sq = r * r;
  double alpha = (1.0/ sqrt( 1.0 + 2.0*M*r/rho_sq));

  fout[0] = alpha;
}

void
eval_ks_lapse_rtheta(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Grab the geometry
  const struct gkyl_dg_gr_maxwell_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];
  double theta = xn[1];

  // Build the lapse
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);
  double alpha = (1.0/ sqrt( 1.0 + 2.0*M*r/rho_sq));

  fout[0] = alpha;
}

void
eval_ks_shift_r(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Grab the geometry
  const struct gkyl_dg_gr_maxwell_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  // theta = pi/2
  double r = xn[0];

  // Build the shift \beta^r
  double rho_sq = r * r;
  double beta_r = (2.0 * M * r) / ( rho_sq + 2.0 * M * r);

  fout[0] = beta_r;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
eval_ks_shift_rtheta(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Grab the geometry
  const struct gkyl_dg_gr_maxwell_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];
  double theta = xn[1];

  // Build the shift \beta^r
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);
  double beta_r = (2.0 * M * r) / ( rho_sq + 2.0 * M * r);

  fout[0] = beta_r;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
eval_ks_h_ij_r(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Grab the geometry
  const struct gkyl_dg_gr_maxwell_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  // theta = pi/2
  double r = xn[0];

  // Metric spatial covariant components  h_ij
  double rho_sq = r * r;
  double h_rr = ( 1.0 + 2.0 * M * r / rho_sq);
  double h_pp = ( rho_sq + a * a * ( 1.0 + 2.0 * M * r / rho_sq )  );
  double h_tt = rho_sq;
  double h_rp = - a * ( 1.0 + 2.0 * M * r / rho_sq );

  fout[0] = h_rr; // h_rr
  fout[1] = 0.0; // h_rt
  fout[2] = h_rp; // h_rp
  fout[3] = h_tt; // h_tt
  fout[4] = 0.0; // h_tp
  fout[5] = h_pp; // h_pp
}

void
eval_ks_h_ij_inv_r(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct gkyl_dg_gr_maxwell_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  // theta = pi/2
  double r = xn[0];

  // Metric spatialcontravariant components  h_ij_inv
  double rho_sq = r * r;
  double h_rr = ( 1.0 + 2.0 * M * r / rho_sq);
  double h_pp = ( rho_sq + a * a * ( 1.0 + 2.0 * M * r / rho_sq )  );
  double h_tt = rho_sq;
  double h_rp = - a * ( 1.0 + 2.0 * M * r / rho_sq );
  double h_rt = 0.0;
  double h_tp = 0.0;

  // Jacobian squared of the spatial metric
  double Jc_sq = h_pp*h_rt*h_rt - 2*h_rt*h_rp*h_tp + h_tt*h_rp*h_rp + h_rr*h_tp*h_tp - h_rr*h_tt*h_pp;

  fout[0] = -(- h_tp*h_tp + h_tt*h_pp)/(Jc_sq); // h_rr_inv
  fout[1] = -(h_rp*h_tp - h_rt*h_pp)/(Jc_sq); // h_rt_inv
  fout[2] = -(h_rt*h_tp - h_rp*h_tt)/(Jc_sq); // h_rp_inv
  fout[3] = -(- h_rp*h_rp + h_rr*h_pp)/(Jc_sq); // h_tt_inv
  fout[4] = -(h_rt*h_rp - h_rr*h_tp)/(Jc_sq); // h_tp_inv
  fout[5] = -(- h_rt*h_rt + h_rr*h_tt)/(Jc_sq); // h_pp_inv
}

void
eval_geom_factor_ks_r(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct gkyl_dg_gr_maxwell_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  // theta = pi/2
  double r = xn[0];

  fout[0] = (-2.0*M*M*a*a - M*a*a*r + 5.0*M*r*r*r + 2.0*r*r*r*r)/((r*r + 2.0*M*r)*(r*r + 2.0*M*r)*r);
  fout[1] = 0.0;
  fout[2] = -(M*a)/((r*r + 2.0*M*r)*r*r);
}



void
eval_ks_h_ij_rtheta(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Grab the geometry
  const struct gkyl_dg_gr_maxwell_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];
  double theta = xn[1];

  // Metric spatial covariant components  h_ij
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);
  double h_rr = ( 1.0 + 2.0 * M * r / rho_sq);
  double h_pp = sin(theta) * sin(theta) * ( rho_sq + a * a * ( 1.0 + 2.0 * M * r / rho_sq ) * sin(theta) * sin(theta)  );
  double h_tt = rho_sq;
  double h_rp = - a * ( 1.0 + 2.0 * M * r / rho_sq ) * sin(theta) * sin(theta);

  fout[0] = h_rr; // h_rr
  fout[1] = 0.0; // h_rt
  fout[2] = h_rp; // h_rp
  fout[3] = h_tt; // h_tt
  fout[4] = 0.0; // h_tp
  fout[5] = h_pp; // h_pp
}

void
eval_ks_h_ij_inv_rtheta(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct gkyl_dg_gr_maxwell_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];
  double theta = xn[1];

  // Metric spatialcontravariant components  h_ij_inv
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);
  double h_rr = ( 1.0 + 2.0 * M * r / rho_sq);
  double h_pp = sin(theta) * sin(theta) * ( rho_sq + a * a * ( 1.0 + 2.0 * M * r / rho_sq ) * sin(theta) * sin(theta)  );
  double h_tt = rho_sq;
  double h_rp = - a * ( 1.0 + 2.0 * M * r / rho_sq ) * sin(theta) * sin(theta);
  double h_rt = 0.0;
  double h_tp = 0.0;

  // Jacobian squared of the spatial metric
  double Jc_sq = h_pp*h_rt*h_rt - 2*h_rt*h_rp*h_tp + h_tt*h_rp*h_rp + h_rr*h_tp*h_tp - h_rr*h_tt*h_pp;

  fout[0] = -(- h_tp*h_tp + h_tt*h_pp)/(Jc_sq); // h_rr_inv
  fout[1] = -(h_rp*h_tp - h_rt*h_pp)/(Jc_sq); // h_rt_inv
  fout[2] = -(h_rt*h_tp - h_rp*h_tt)/(Jc_sq); // h_rp_inv
  fout[3] = -(- h_rp*h_rp + h_rr*h_pp)/(Jc_sq); // h_tt_inv
  fout[4] = -(h_rt*h_rp - h_rr*h_tp)/(Jc_sq); // h_tp_inv
  fout[5] = -(- h_rt*h_rt + h_rr*h_tt)/(Jc_sq); // h_pp_inv
}



void
eval_geom_factor_ks_rtheta(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct gkyl_dg_gr_maxwell_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];
  double theta = xn[1];

  // Some simplfiying terms to make the algebra more readable
  double rho_sq = r*r + a*a*cos(theta)*cos(theta);
  double rho_sq_M  = rho_sq + 2*M*r;
  double term1 = rho_sq*(r*r + a*a) + 2*M*a*a*r*sin(theta)*sin(theta);

  fout[0] = (4*a*a*r + 4*r*r*r + 2*M*a*a*sin(theta)*sin(theta) - 2*a*a*r*sin(theta)*sin(theta))/(rho_sq_M*rho_sq) - (term1*((M + r)/rho_sq_M + r/rho_sq))/(rho_sq_M*rho_sq);
  fout[1] = (cos(theta)*(rho_sq_M*rho_sq + 2*M*a*a*r*sin(theta)*sin(theta)))/(rho_sq_M*rho_sq*rho_sq*sin(theta));
  fout[2] = (M*a*(a*a*cos(theta)*cos(theta) - r*r))/(rho_sq_M*rho_sq*rho_sq);
}

void
eval_ks_det_h_r(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Grab the geometry
  const struct gkyl_dg_gr_maxwell_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  // theta = pi/2
  double r = xn[0];

  // Build the det_h \beta^r
  double rho_sq = r * r;
  double sqrt_det_h = sqrt(rho_sq) * sqrt(2 * M * r + rho_sq);

  fout[0] = sqrt_det_h;
}

void
eval_ks_det_h_rtheta(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Grab the geometry
  const struct gkyl_dg_gr_maxwell_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];
  double theta = xn[1];

  // Build the det_h \beta^r
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);
  double sqrt_det_h = sqrt(rho_sq) * sin(theta) * sqrt(2 * M * r + rho_sq);

  fout[0] = sqrt_det_h;
}

static evalf_t
choose_lapse_kern(enum gkyl_triad_preset_geom_type type)
{
  switch(type) {
    case GKYL_TRIAD_FLAT:
    case GKYL_TRIAD_ANNULUS:
    case GKYL_TRIAD_SPHERICAL_RTHETA:
      return eval_non_relativistic_lapse;
      break;

    case GKYL_TRIAD_GR_KERR_SCHILD_RPHI:
    case GKYL_TRIAD_GR_KERR_SCHILD_R:
      return eval_ks_lapse_r;
      break;

    case GKYL_TRIAD_GR_KERR_SCHILD_RTHETA:
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return eval_ks_lapse_rtheta;
      break;

    default:
      assert(false);
  }
}

static evalf_t
choose_shift_kern(enum gkyl_triad_preset_geom_type type)
{
  switch(type) {
    case GKYL_TRIAD_FLAT:
    case GKYL_TRIAD_ANNULUS:
    case GKYL_TRIAD_SPHERICAL_RTHETA:
      return eval_non_relativistic_shift;
      break;

    case GKYL_TRIAD_GR_KERR_SCHILD_RPHI:
    case GKYL_TRIAD_GR_KERR_SCHILD_R:
      return eval_ks_shift_r;
      break;

    case GKYL_TRIAD_GR_KERR_SCHILD_RTHETA:
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return eval_ks_shift_rtheta;
      break;
      
    default:
      assert(false);
  }
}

static evalf_t
choose_geom_factor_con_kern(enum gkyl_triad_preset_geom_type type)
{
  switch(type) {
    case GKYL_TRIAD_FLAT:
      return eval_geom_factor_flat;
      break;

    case GKYL_TRIAD_ANNULUS:
      return eval_geom_factor_annulus;
      break;

    case GKYL_TRIAD_SPHERICAL_RTHETA:
      return eval_geom_factor_sph_rtheta;
      break;

    case GKYL_TRIAD_GR_KERR_SCHILD_RPHI:
    case GKYL_TRIAD_GR_KERR_SCHILD_R:
      return eval_geom_factor_ks_r;
      break;

    case GKYL_TRIAD_GR_KERR_SCHILD_RTHETA:
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return eval_geom_factor_ks_rtheta;
      break;

    default:
      assert(false);
  }
}

static evalf_t
choose_h_ij_kern(enum gkyl_triad_preset_geom_type type)
{
  switch(type) {
    case GKYL_TRIAD_FLAT:
      return eval_flat_h_ij;
      break;

    case GKYL_TRIAD_ANNULUS:
      return eval_annulus_h_ij;
      break;

    case GKYL_TRIAD_SPHERICAL_RTHETA:
      return eval_spherical_h_ij_rtheta;
      break;

    case GKYL_TRIAD_GR_KERR_SCHILD_RPHI:
    case GKYL_TRIAD_GR_KERR_SCHILD_R:
      return eval_ks_h_ij_r;
      break;

    case GKYL_TRIAD_GR_KERR_SCHILD_RTHETA:
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return eval_ks_h_ij_rtheta;
      break;

    default:
      assert(false);
  }
}


static evalf_t
choose_h_ij_inv_kern(enum gkyl_triad_preset_geom_type type)
{
  switch(type) {
    case GKYL_TRIAD_FLAT:
      return eval_flat_h_ij_inv;
      break;

    case GKYL_TRIAD_ANNULUS:
      return eval_annulus_h_ij_inv;
      break;

    case GKYL_TRIAD_SPHERICAL_RTHETA:
      return eval_spherical_h_ij_inv_rtheta;
      break;

    case GKYL_TRIAD_GR_KERR_SCHILD_RPHI:
    case GKYL_TRIAD_GR_KERR_SCHILD_R:
      return eval_ks_h_ij_inv_r;
      break;

    case GKYL_TRIAD_GR_KERR_SCHILD_RTHETA:
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return eval_ks_h_ij_inv_rtheta;
      break;

    default:
      assert(false);
  }
}

static evalf_t
choose_det_h_kern(enum gkyl_triad_preset_geom_type type)
{
  switch(type) {
    case GKYL_TRIAD_FLAT:
      return eval_flat_det_h;
      break;

    case GKYL_TRIAD_ANNULUS:
      return eval_annulus_det_h;
      break;

    case GKYL_TRIAD_SPHERICAL_RTHETA:
      return eval_spherical_det_h_rtheta;
      break;

    case GKYL_TRIAD_GR_KERR_SCHILD_RPHI:
    case GKYL_TRIAD_GR_KERR_SCHILD_R:
      return eval_ks_det_h_r;
      break;

    case GKYL_TRIAD_GR_KERR_SCHILD_RTHETA:
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return eval_ks_det_h_rtheta;
      break;

    default:
      assert(false);
  }
}
