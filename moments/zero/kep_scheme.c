#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_kep_scheme.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_euler_priv.h>

#include <float.h>
#include <math.h>

#define RHO 0     
#define RHOU d[0]
#define RHOV d[1]
#define RHOW d[2]
#define ER 4
#define PR 4

#define SQ(x) ((x)*(x))

static const int dir_shuffle[][3] = {
  {1, 2, 3},
  {2, 3, 1},
  {3, 1, 2}
};

struct gkyl_kep_scheme {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update

  bool use_hybrid_flux; // should we use shock detector for hybrid flux
  
  double cfl; // CFL number
  
  const struct gkyl_wv_eqn *equation; // equation object
  struct gkyl_wave_geom *geom; // geometry object  
};

gkyl_kep_scheme*
gkyl_kep_scheme_new(const struct gkyl_kep_scheme_inp *inp)
{
  struct gkyl_kep_scheme *up;
  up = gkyl_malloc(sizeof(*up));

  up->grid = *(inp->grid);
  up->ndim = up->grid.ndim;
  
  up->num_up_dirs = inp->num_up_dirs;
  for (int i=0; i<inp->num_up_dirs; ++i)
    up->update_dirs[i] = inp->update_dirs[i];

  up->cfl = inp->cfl;
  up->use_hybrid_flux = inp->use_hybrid_flux;

  up->equation = gkyl_wv_eqn_acquire(inp->equation);
  up->geom = gkyl_wave_geom_acquire(inp->geom);

  return up;
}

// compute kinetic energy (no density factor)
static inline
double euler_ke(const double v[5])
{
  return 0.5*(v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
}

static inline double
euler_flux(int dir, double gas_gamma, const double q[5], double flux[5])
{
  const int *d = dir_shuffle[dir];

  double pr = gkyl_euler_pressure(gas_gamma, q), u = q[RHOU]/q[0];
  flux[0] = q[RHOU]; // rho*u
  flux[RHOU] = q[RHOU]*u + pr; // rho*u*u + pr
  flux[RHOV] = q[RHOV]*u; // rho*v*u
  flux[RHOW] = q[RHOW]*u; // rho*w*u
  flux[4] = (q[4]+pr)*u; // (E+p)*u  

  double u2 = sqrt( SQ(q[1]/q[0]) + SQ(q[2]/q[0]) + SQ(q[3]/q[0]) );
  return fabs(u2) + sqrt(gas_gamma*pr/q[0]);  
}

// Lax fluxes
static inline void
mlax_flux(int dir, double gas_gamma, const double qm[5], const double qp[5], double flux[5])
{
  double fm[5], fp[5];

  double amaxp = euler_flux(dir, gas_gamma, qp, fp);
  double amaxm = euler_flux(dir, gas_gamma, qm, fm);

  for (int i=0; i<5; ++i)
    flux[i] = 0.5*(fp[i]+fm[i]) - 0.5*fmax(amaxm, amaxp)*(qp[i]-qm[i]);
}

// Numerical flux using modified KEP scheme
static inline void
mkep_flux(int dir, double gas_gamma, const double vm[5], const double vp[5], double flux[5])
{
  double vbar[5];
  
  for (int i=0; i<5; ++i) vbar[i] = 0.5*(vm[i]+vp[i]);
  double kebar = 0.5*(euler_ke(vm) + euler_ke(vp));
  
  const int *d = dir_shuffle[dir];
  flux[0] = vbar[RHO]*vbar[RHOU]; // rho*u
  // momentum flux must have this form to ensure KEP property
  flux[RHOU] =  flux[0]*vbar[RHOU] + vbar[PR]; // rho*u*u + pe
  flux[RHOV] =  flux[0]*vbar[RHOV]; // rho*u*v
  flux[RHOW] =  flux[0]*vbar[RHOW]; // rho*u*v
  // following ensure stability of linear perturbations around uniform flow
  flux[ER] = gas_gamma/(gas_gamma-1)*vbar[PR]*vbar[RHOU] + flux[0]*kebar; // (E+p)*u
}

static inline long
get_offset(int dir, int loc, const struct gkyl_range *range)
{
  int idx[GKYL_MAX_CDIM] = { 0, 0, 0 };
  idx[dir] = loc;
  return gkyl_range_offset(range, idx);
}

static void
calc_alpha(const gkyl_kep_scheme *kep, const struct gkyl_range *update_rng,
  const struct gkyl_array *qin,   struct gkyl_array *alpha)
{
  int ndim = update_rng->ndim;
  gkyl_array_clear_range(alpha, 0.0, update_rng);

  double eps = 1e-14; // to prevent div by 0.0
  double nu = 0.1;
  double phi = 0.01; // base diffusion

  enum { IRHO, IRHOU, IRHOV, IRHOW}; // indexing Euler conserved vars

  enum { IC, IL, IR }; // indexing into offsets
  long offsets[GKYL_MAX_CDIM][3];
  
  for (int d=0; d<ndim; ++d) {
    offsets[d][IC] = get_offset(d, 0, update_rng);
    offsets[d][IL] = get_offset(d, -1, update_rng);
    offsets[d][IR] = get_offset(d, 1, update_rng);
  }

  double L = DBL_MAX, dx[GKYL_MAX_CDIM];
  
  for (int d=0; d<ndim; ++d) {
    dx[d] = kep->grid.dx[d];
    L = fmin(L, dx[d]);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_rng);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(update_rng, iter.idx);
    
    double divu = 0.0, curlx = 0.0, curly = 0.0, curlz = 0.0;
    double u2 = 0.0; // u*u

    const double *qc = gkyl_array_cfetch(qin, loc);

    // compute div(u), curl(u) and u*u
    for (int d=0; d<ndim; ++d) {
      
      const double *ql = gkyl_array_cfetch(qin, loc+offsets[d][IL]);
      const double *qr = gkyl_array_cfetch(qin, loc+offsets[d][IR]);

      u2 += SQ(qc[d+1]/qc[0]);
      divu += (qr[d+1]/qr[0] - ql[d+1]/ql[0])/(2.0*dx[d]);
      //divu += (qr[d+1] - ql[d+1])/(2.0*dx[d]*qc[0]);

      if (d == 0) {
        curly += -(qr[IRHOW]/qr[0] - ql[IRHOW]/ql[0])/(2*dx[0]);
        curlz += (qr[IRHOV]/qr[0] - ql[IRHOV]/ql[0])/(2*dx[0]);
      }
      else if (d == 1) {
        curlx += (qr[IRHOW]/qr[0] - ql[IRHOW]/ql[0])/(2*dx[1]);
        curlz += -(qr[IRHOU]/qr[0] - ql[IRHOU]/ql[0])/(2*dx[1]);
      }
      else {
        curlx += -(qr[IRHOV]/qr[0] - ql[IRHOV]/ql[0])/(2*dx[2]);
        curly += (qr[IRHOU]/qr[0] - ql[IRHOU]/ql[0])/(2*dx[2]);
      }
    }

    double div2 = divu*divu;
    double curl2 = curlx*curlx + curly*curly + curlz*curlz;
    double Omega2 = nu*nu*u2/(L*L);
    
    double *al = gkyl_array_fetch(alpha, loc);

    double ducros = fmin(4.0/3.0*div2/(div2+curl2+eps), 1.0);
    double theta = div2/(div2+Omega2+eps);
    al[0] = fmax(ducros-phi, 0)*theta + phi;
  }
}

void
gkyl_kep_scheme_advance(const gkyl_kep_scheme *kep, const struct gkyl_range *update_rng,
  const struct gkyl_array *qin,   struct gkyl_array *alpha,
  struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  int ndim = update_rng->ndim;
  double gas_gamma = gkyl_wv_euler_gas_gamma(kep->equation);

  int idxm[GKYL_MAX_DIM], idxp[GKYL_MAX_DIM];

  double vm[5], vp[5], flux[5], lflux[5];

  gkyl_array_clear_range(rhs, 0.0, update_rng);

  if (kep->use_hybrid_flux)
    calc_alpha(kep, update_rng, qin, alpha);

  for (int d=0; d<kep->num_up_dirs; ++d) {
    int dir = kep->update_dirs[d];
    double dx = kep->grid.dx[dir];

    int loidx = update_rng->lower[dir];
    int upidx = update_rng->upper[dir]+1; // one more edge than cells

    struct gkyl_range perp_range;
    gkyl_range_shorten_from_above(&perp_range, update_rng, dir, 1);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &perp_range);

    while (gkyl_range_iter_next(&iter)) {

      gkyl_copy_int_arr(ndim, iter.idx, idxm);
      gkyl_copy_int_arr(ndim, iter.idx, idxp);
      
      for (int i=loidx; i<=upidx; ++i) { // note upidx is inclusive
        idxm[dir] = i-1; idxp[dir] = i;

        long linm = gkyl_range_idx(update_rng, idxm);
        long linp = gkyl_range_idx(update_rng, idxp);

        const double *qm = gkyl_array_cfetch(qin, linm);
        const double *qp = gkyl_array_cfetch(qin, linp);

        gkyl_euler_prim_vars(gas_gamma, qm, vm);
        gkyl_euler_prim_vars(gas_gamma, qp, vp);

        mkep_flux(dir, gas_gamma, vm, vp, flux);

        if (kep->use_hybrid_flux)
          mlax_flux(dir, gas_gamma, qm, qp, lflux);

        // accumulate contribution of flux to left/right cell
        double *rhsm = gkyl_array_fetch(rhs, linm);
        double *rhsp = gkyl_array_fetch(rhs, linp);

        if (kep->use_hybrid_flux) {
          const double *alpham = gkyl_array_fetch(alpha, linm);
          const double *alphap = gkyl_array_fetch(alpha, linp);
          double am = fmax(alpham[0], alphap[0]);

          for (int m=0; m<5; ++m) {
            rhsp[m] += ((1.0-am)*flux[m] + am*lflux[m])/dx;
            rhsm[m] += -((1.0-am)*flux[m] + am*lflux[m])/dx;
          }
        }
        else {
          for (int m=0; m<5; ++m) {
            rhsp[m] += flux[m]/dx;
            rhsm[m] += -flux[m]/dx;
          }
        }

        double *cflrate_d = gkyl_array_fetch(cflrate, linp);

        // rotate q to local coordinates before computing max speed
        const int *d = dir_shuffle[dir];
        double qp_dir[5];
        qp_dir[0] = qp[0];
        qp_dir[RHOU] = qp[1]; qp_dir[RHOV] = qp[2]; qp_dir[RHOW] = qp[3];
        qp_dir[4] = qp[4];
        
        cflrate_d[0] += gkyl_euler_max_abs_speed(gas_gamma, qp_dir)/dx;
      }
    }
  }
}

double
gkyl_kep_scheme_max_dt(const gkyl_kep_scheme *kep, const struct gkyl_range *update_range,
  const struct gkyl_array *qin)
{
  double max_dt = DBL_MAX;
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {

    for (int d=0; d<kep->num_up_dirs; ++d) {
      int dir = kep->update_dirs[d];
      double dx = kep->grid.dx[dir];

      const double *q = gkyl_array_cfetch(qin, gkyl_range_idx(update_range, iter.idx));
      double maxs = gkyl_wv_eqn_max_speed(kep->equation, q);
      max_dt = fmin(max_dt, kep->cfl*dx/maxs);
    }
    
  }

  return max_dt;  
}

void
gkyl_kep_scheme_release(gkyl_kep_scheme* up)
{
  gkyl_wv_eqn_release(up->equation);
  gkyl_wave_geom_release(up->geom);
  gkyl_free(up);
}
