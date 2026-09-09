#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_gkgeom.h>
#include <gkyl_math.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

#include <math.h>
#include <string.h>

struct gkyl_gkgeom {
  struct gkyl_rect_grid rzgrid; // RZ grid on which psi(R,Z) is defined
  const struct gkyl_array *psiRZ; // psi(R,Z) DG representation
  struct gkyl_range rzlocal; // local range over which psiRZ is defined
  int num_rzbasis; // number of basis functions in RZ

  struct { int max_iter; double eps; } root_param;
  struct { int max_level; double eps; } quad_param;

  // pointer to root finder (depends on polyorder)
  struct RdRdZ_sol (*calc_roots)(const double *psi, double psi0, double Z,
    double xc[2], double dx[2]);

  struct gkyl_gkgeom_stat stat; 
};

// some helper functions
static inline double
choose_closest(double ref, double R[2], double out[2])
{
  return fabs(R[0]-ref) < fabs(R[1]-ref) ? out[0] : out[1];
}

static inline double SQ(double x) { return x*x; }

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

// Compute roots R(psi,Z) and dR/dZ(psi,Z) in a p=2 DG cell
static inline struct RdRdZ_sol
calc_RdR_ser_p2(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };
  double y = (Z-xc[1])/(dx[1]*0.5);

  double aq = 2.904737509655563*psi[6]*y+1.677050983124842*psi[4]; 
  double bq = 2.904737509655563*psi[7]*SQ(y)+1.5*psi[3]*y-0.9682458365518543*psi[7]+0.8660254037844386*psi[1]; 
  double cq = 1.677050983124842*psi[5]*SQ(y)-0.9682458365518543*psi[6]*y+0.8660254037844386*psi[2]*y-1.0*psi0-0.5590169943749475*psi[5]-0.5590169943749475*psi[4]+0.5*psi[0]; 
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
      double C = 5.809475019311126*psi[7]*x*y+3.354101966249685*psi[5]*y+2.904737509655563*psi[6]*SQ(x)+1.5*psi[3]*x-0.9682458365518543*psi[6]+0.8660254037844386*psi[2]; 
      double A = 2.904737509655563*psi[7]*SQ(y)+5.809475019311126*psi[6]*x*y+1.5*psi[3]*y+3.354101966249685*psi[4]*x-0.9682458365518543*psi[7]+0.8660254037844386*psi[1];
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
    if ((-1<=r2) && (r2 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r2*dx[0]*0.5 + xc[0];

      double x = r2;
      double C = 5.809475019311126*psi[7]*x*y+3.354101966249685*psi[5]*y+2.904737509655563*psi[6]*SQ(x)+1.5*psi[3]*x-0.9682458365518543*psi[6]+0.8660254037844386*psi[2]; 
      double A = 2.904737509655563*psi[7]*SQ(y)+5.809475019311126*psi[6]*x*y+1.5*psi[3]*y+3.354101966249685*psi[4]*x-0.9682458365518543*psi[7]+0.8660254037844386*psi[1];
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
  }
  return sol;
}

// Compute roots R(psi,Z) and dR/dZ(psi,Z) in a p=2 DG cell
static inline struct RdRdZ_sol
calc_RdR_ten_p2(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };
  double y = (Z-xc[1])/(dx[1]*0.5);

  double aq = 2.904737509655563*psi[6]*y+1.677050983124842*psi[4]; 
  double bq = 2.904737509655563*psi[7]*SQ(y)+1.5*psi[3]*y-0.9682458365518543*psi[7]+0.8660254037844386*psi[1]; 
  double cq = 1.677050983124842*psi[5]*SQ(y)-0.9682458365518543*psi[6]*y+0.8660254037844386*psi[2]*y-1.0*psi0-0.5590169943749475*psi[5]-0.5590169943749475*psi[4]+0.5*psi[0]; 
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
      double C = 5.809475019311126*psi[7]*x*y+3.354101966249685*psi[5]*y+2.904737509655563*psi[6]*SQ(x)+1.5*psi[3]*x-0.9682458365518543*psi[6]+0.8660254037844386*psi[2]; 
      double A = 2.904737509655563*psi[7]*SQ(y)+5.809475019311126*psi[6]*x*y+1.5*psi[3]*y+3.354101966249685*psi[4]*x-0.9682458365518543*psi[7]+0.8660254037844386*psi[1];
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
    if ((-1<=r2) && (r2 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r2*dx[0]*0.5 + xc[0];

      double x = r2;
      double C = 5.809475019311126*psi[7]*x*y+3.354101966249685*psi[5]*y+2.904737509655563*psi[6]*SQ(x)+1.5*psi[3]*x-0.9682458365518543*psi[6]+0.8660254037844386*psi[2]; 
      double A = 2.904737509655563*psi[7]*SQ(y)+5.809475019311126*psi[6]*x*y+1.5*psi[3]*y+3.354101966249685*psi[4]*x-0.9682458365518543*psi[7]+0.8660254037844386*psi[1];
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
  }
  return sol;
}

// Compute R(psi,Z) given a psi and Z. Can return multiple solutions
// or no solutions. The number of roots found is returned and are
// copied in the array R and dR. The calling function must ensure that
// these arrays are big enough to hold all roots required
static int
R_psiZ(const gkyl_gkgeom *geo, double psi, double Z, int nmaxroots,
  double *R, double *dR)
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
        R[sidx] = sol.R[s];
        dR[sidx] = sol.dRdZ[s];
        sidx += 1;
      }
  }
  return sidx;
}

// Function context to pass to coutour integration function
struct contour_ctx {
  const gkyl_gkgeom *geo;
  double psi, last_R;
  long ncall;
};

// Function to pass to numerical quadrature to integrate along a contour
static inline double
contour_func(double Z, void *ctx)
{
  struct contour_ctx *c = ctx;
  c->ncall += 1;
  double R[2] = { 0 }, dR[2] = { 0 };
  
  int nr = R_psiZ(c->geo, c->psi, Z, 2, R, dR);
  double dRdZ = nr == 1 ? dR[0] : choose_closest(c->last_R, R, dR);
  
  return nr>0 ? sqrt(1+dRdZ*dRdZ) : 0.0;
}

// Integrates along a specified contour, optionally using a "memory"
// of previously computed values, or storing computed values in
// memory. The function basically breaks up the integral into a loop
// over z-cells. This needs to be done as the DG representation is,
// well, discontinuous, and adaptive quadrature struggles with such
// functions.
static double
integrate_psi_contour_memo(const gkyl_gkgeom *geo, double psi,
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
  
  double dz = geo->rzgrid.dx[1];
  double zlo = geo->rzgrid.lower[1];
  int izlo = geo->rzlocal.lower[1], izup = geo->rzlocal.upper[1];
  
  int ilo = get_idx(1, zmin, &geo->rzgrid, &geo->rzlocal);
  int iup = get_idx(1, zmax, &geo->rzgrid, &geo->rzlocal);

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

  ((gkyl_gkgeom *)geo)->stat.nquad_cont_calls += ctx.ncall;
  return res;
}

// Function context to pass to root finder
struct arc_length_ctx {
  const gkyl_gkgeom *geo;
  double *arc_memo;
  double psi, rclose, zmin, arcL;
};


// Function to pass to root-finder to find Z location for given arc-length
static inline double
arc_length_func(double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL;
  double ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose,
    true, false, arc_memo) - arcL;
  return ival;
}

gkyl_gkgeom*
gkyl_gkgeom_new(const struct gkyl_gkgeom_inp *inp)
{
  struct gkyl_gkgeom *geo = gkyl_malloc(sizeof(*geo));

  geo->rzgrid = *inp->rzgrid;
  geo->psiRZ = gkyl_array_acquire(inp->psiRZ);
  geo->num_rzbasis = inp->rzbasis->num_basis;
  memcpy(&geo->rzlocal, inp->rzlocal, sizeof(struct gkyl_range));

  geo->root_param.eps =
    inp->root_param.eps > 0 ? inp->root_param.eps : 1e-10;
  geo->root_param.max_iter =
    inp->root_param.max_iter > 0 ? inp->root_param.max_iter : 100;

  geo->quad_param.max_level =
    inp->quad_param.max_levels > 0 ? inp->quad_param.max_levels : 10;
  geo->quad_param.eps =
    inp->quad_param.eps > 0 ? inp->quad_param.eps : 1e-10;

  if (inp->rzbasis->poly_order == 1) {
    geo->calc_roots = calc_RdR_p1;
  }
  else if (inp->rzbasis->poly_order == 2) {
    if (inp->rzbasis->b_type == GKYL_BASIS_MODAL_SERENDIPITY)
      geo->calc_roots = calc_RdR_ser_p2;
    else
      geo->calc_roots = calc_RdR_ten_p2;
  }

  geo->stat = (struct gkyl_gkgeom_stat) { };
  
  return geo;
}

double
gkyl_gkgeom_integrate_psi_contour(const gkyl_gkgeom *geo, double psi,
  double zmin, double zmax, double rclose)
{
  return integrate_psi_contour_memo(geo, psi, zmin, zmax, rclose,
    false, false, 0);
}

int
gkyl_gkgeom_R_psiZ(const gkyl_gkgeom *geo, double psi, double Z, int nmaxroots,
  double *R, double *dR)
{
  return R_psiZ(geo, psi, Z, nmaxroots, R, dR);
}

// write out nodal coordinates 
static void
write_nodal_coordinates(const char *nm, struct gkyl_range *nrange,
  struct gkyl_array *nodes)
{
  double lower[3] = { 0.0, 0.0, 0.0 };
  double upper[3] = { 1.0, 1.0, 1.0 };
  int cells[3];
  for (int i=0; i<nrange->ndim; ++i)
    cells[i] = gkyl_range_shape(nrange, i);
  
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  gkyl_grid_sub_array_write(&grid, nrange, 0, nodes, nm);
}

void
gkyl_gkgeom_calcgeom(const gkyl_gkgeom *geo,
  const struct gkyl_gkgeom_geo_inp *inp, struct gkyl_array *mapc2p)
{
  int poly_order = inp->cbasis->poly_order;
  int nodes[3] = { 1, 1, 1 };
  if (poly_order == 1)
    for (int d=0; d<inp->cgrid->ndim; ++d)
      nodes[d] = inp->cgrid->cells[d]+1;
  if (poly_order == 2)
    for (int d=0; d<inp->cgrid->ndim; ++d)
      nodes[d] = 2*inp->cgrid->cells[d]+1;

  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, inp->cgrid->ndim, nodes);
  struct gkyl_array *mc2p = gkyl_array_new(GKYL_DOUBLE, inp->cgrid->ndim, nrange.volume);

  enum { TH_IDX, PH_IDX, AL_IDX }; // arrangement of computational coordinates
  enum { R_IDX, Z_IDX }; // arrangement of physical coordinates  
  
  double dtheta = inp->cgrid->dx[TH_IDX],
    dphi = inp->cgrid->dx[PH_IDX],
    dalpha = inp->cgrid->dx[AL_IDX];
  
  double theta_lo = inp->cgrid->lower[TH_IDX],
    phi_lo = inp->cgrid->lower[PH_IDX],
    alpha_lo = inp->cgrid->lower[AL_IDX];

  double dx_fact = poly_order == 1 ? 1 : 0.5;
  dtheta *= dx_fact; dphi *= dx_fact; dalpha *= dx_fact;

  double rclose = inp->rclose;

  int nzcells = geo->rzgrid.cells[1];
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo
  };

  int cidx[2] = { 0 };
  for (int ip=nrange.lower[PH_IDX]; ip<=nrange.upper[PH_IDX]; ++ip) {

    double zmin = inp->zmin, zmax = inp->zmax;

    double psi_curr = phi_lo + ip*dphi;
    double arcL = integrate_psi_contour_memo(geo, psi_curr, zmin, zmax, rclose,
      true, true, arc_memo);

    double delta_arcL = arcL/(poly_order*inp->cgrid->cells[TH_IDX]);

    cidx[PH_IDX] = ip;

    do {
      // set node coordinates of first node
      cidx[TH_IDX] = nrange.lower[TH_IDX];
      double *mc2p_n = gkyl_array_fetch(mc2p, gkyl_range_idx(&nrange, cidx));
      mc2p_n[Z_IDX] = zmin;
      double R[2] = { 0 }, dR[2] = { 0 };    
      int nr = R_psiZ(geo, psi_curr, zmin, 2, R, dR);
      mc2p_n[R_IDX] = choose_closest(rclose, R, R);
    } while(0);

    // set node coordinates of rest of nodes
    double arcL_curr = 0.0;
    for (int it=nrange.lower[TH_IDX]+1; it<nrange.upper[TH_IDX]; ++it) {
      arcL_curr += delta_arcL;

      arc_ctx.psi = psi_curr;
      arc_ctx.rclose = rclose;
      arc_ctx.zmin = zmin;
      arc_ctx.arcL = arcL_curr;

      struct gkyl_qr_res res = gkyl_ridders(arc_length_func, &arc_ctx,
        zmin, zmax, -arcL_curr, arcL-arcL_curr,
        geo->root_param.max_iter, 1e-10);
      double z_curr = res.res;
      ((gkyl_gkgeom *)geo)->stat.nroot_cont_calls += res.nevals;

      double R[2] = { 0 }, dR[2] = { 0 };
      int nr = R_psiZ(geo, psi_curr, z_curr, 2, R, dR);
      double r_curr = choose_closest(rclose, R, R);

      cidx[TH_IDX] = it;
      double *mc2p_n = gkyl_array_fetch(mc2p, gkyl_range_idx(&nrange, cidx));
      mc2p_n[Z_IDX] = z_curr;
      mc2p_n[R_IDX] = r_curr;
    }

    do {
      // set node coordinates of last node
      cidx[TH_IDX] = nrange.upper[TH_IDX];
      double *mc2p_n = gkyl_array_fetch(mc2p, gkyl_range_idx(&nrange, cidx));
      mc2p_n[Z_IDX] = zmax;
      double R[2] = { 0 }, dR[2] = { 0 };    
      int nr = R_psiZ(geo, psi_curr, zmax, 2, R, dR);
      mc2p_n[R_IDX] = choose_closest(rclose, R, R);
    } while (0);
  }

  if (inp->write_node_coord_array)
    write_nodal_coordinates(inp->node_file_nm, &nrange, mc2p);

  gkyl_free(arc_memo);
  gkyl_array_release(mc2p);  
}

struct gkyl_gkgeom_stat
gkyl_gkgeom_get_stat(const gkyl_gkgeom *geo)
{
  return geo->stat;
}

void
gkyl_gkgeom_release(gkyl_gkgeom *geo)
{
  gkyl_array_release(geo->psiRZ);
  gkyl_free(geo);
}
