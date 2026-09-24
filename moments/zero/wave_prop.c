#include <float.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_null_comm.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_prop.h>

#include <gkyl_level_set.h>

struct gkyl_wave_prop {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  enum gkyl_wave_limiter limiter; // limiter to use
  double cfl; // CFL number
  const struct gkyl_wv_eqn *equation; // equation object

  bool force_low_order_flux; // only use Lax flux
  bool check_inv_domain; // flag to indicate if invariant domains are checked

  enum gkyl_wave_split_type split_type; // type of splitting to use

  struct gkyl_wave_geom *geom; // geometry object
  struct gkyl_comm *comm; // communcator
  
  // data for 1D slice update
  struct gkyl_array *waves, *apdq, *amdq, *speeds, *flux2;
  // flags to indicate if fluctuations should be recomputed
  struct gkyl_array *redo_fluct;

  // some stats
  long n_calls; // number of calls to updater
  long n_bad_advance_calls; // number of calls in which positivity had to be fixed
  long n_bad_cells; // number  of cells fixed
  long n_max_bad_cells; // maximum number of cells fixed in a call
};

static inline double
fmax3(double a, double b, double c)
{
  return fmax(fmax(a,b),c);
}

static inline double
fmin3(double a, double b, double c)
{
  return fmin(fmin(a,b),c);
}

// limiter function
static inline double
limiter_function(double r, enum gkyl_wave_limiter limiter)
{
  double theta = 0.0;
  switch (limiter) {
    case GKYL_NO_LIMITER:
      theta = 1.0;
      break;
    
    // ** Fully formally-verified implementation of the minmod flux limiter **
    // ** Proof of symmetry (equivalent action on forward and backward gradients): ../proofs/finite_volume/proof_limiter_minmod_symmetry.rkt **
    // ** Proof of second-order TVD (total variation diminishing): ../proofs/finite_volume/proof_limiter_minmod_tvd.rkt **
    case GKYL_MIN_MOD:
      theta = fmax(0.0, fmin(1.0, r));
      break;

    // ** Partially formally-verified implementation of the superbee flux limiter **
    // ** Proof of symmetry (equivalent action on forward and backward gradients): NOT PROVEN **
    // ** Proof of second-order TVD (total variation diminishing): ../proofs/finite_volume/proof_limiter_superbee_tvd.rkt **
    case GKYL_SUPERBEE:
      theta = fmax3(0.0, fmin((2.0 * r), 1.0), fmin(r, 2.0));
      break;

    // ** Partially formally-verified implementation of the van Leer flux limiter **
    // ** Proof of symmetry (equivalent action on forward and backward gradients): ../proofs/finite_volume/proof_limiter_van_leer_symmetry.rkt **
    // ** Proof of second-order TVD (total variation diminishing): NOT PROVEN **
    case GKYL_VAN_LEER:
      theta = ((r + fabs(r)) / (1.0 + fabs(r)));
      break;

    // ** Fully formally-verified implementation of the monotonized-centered flux limiter **
    // ** Proof of symmetry (equivalent action on forward and backward gradients): ../proofs/finite_volume/proof_limiter_monotonized_centered_symmetry.rkt **
    // ** Proof of second-order TVD (total variation diminishing): ../proofs/finite_volume/proof_limiter_monotonized_centered_tvd.rkt **
    case GKYL_MONOTONIZED_CENTERED:
      theta = fmax(0.0, fmin3((2.0 * r), ((1.0 + r) / 2.0), 2.0));
      break;

    case GKYL_BEAM_WARMING:
      theta = r;
      break;

    case GKYL_ZERO:
      theta = 0;
      break;
  }
  return theta;
}

gkyl_wave_prop*
gkyl_wave_prop_new(const struct gkyl_wave_prop_inp *winp)
{
  gkyl_wave_prop *up = gkyl_malloc(sizeof(*up));

  up->grid = *(winp->grid);
  up->ndim = up->grid.ndim;
  
  up->num_up_dirs = winp->num_up_dirs;
  for (int i=0; i<winp->num_up_dirs; ++i)
    up->update_dirs[i] = winp->update_dirs[i];

  up->limiter = winp->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : winp->limiter;
  up->cfl = winp->cfl;
  up->equation = gkyl_wv_eqn_acquire(winp->equation);

  if (winp->comm)
    up->comm = gkyl_comm_acquire(winp->comm);
  else
    up->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { } );

  up->force_low_order_flux = winp->force_low_order_flux;
  up->check_inv_domain = winp->check_inv_domain;

  up->split_type = winp->split_type;

  int nghost[3] = { 2, 2, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&up->grid, nghost, &ext_range, &range);

  int max_1d = 0;
  for (int d=0; d<ext_range.ndim; ++d) {
    int shape = gkyl_range_shape(&ext_range, d);
    max_1d = max_1d > shape ? max_1d : shape;
  }

  // allocate memory to store 1D slices of waves, speeds and
  // second-order correction flux
  int meqn = winp->equation->num_equations, mwaves = winp->equation->num_waves;
  up->waves = gkyl_array_new(GKYL_DOUBLE, meqn*mwaves, max_1d);
  up->apdq = gkyl_array_new(GKYL_DOUBLE, meqn, max_1d);
  up->amdq = gkyl_array_new(GKYL_DOUBLE, meqn, max_1d);
  up->speeds = gkyl_array_new(GKYL_DOUBLE, mwaves, max_1d);
  up->flux2 = gkyl_array_new(GKYL_DOUBLE, meqn, max_1d);

  up->redo_fluct = gkyl_array_new(GKYL_DOUBLE, meqn, max_1d);

  up->geom = gkyl_wave_geom_acquire(winp->geom);

  up->n_calls = up->n_bad_advance_calls = 0;
  up->n_bad_cells = up->n_max_bad_cells = 0;

  return up;
}

// some helper functions

static inline void
copy_wv_vec(int n, double * GKYL_RESTRICT out, const double * GKYL_RESTRICT inp)
{
  for (int i=0; i<n; ++i) out[i] = inp[i];
}

static inline void
calc_jump(int n, const double *ql, const double *qr, double * GKYL_RESTRICT jump)
{
  for (int d=0; d<n; ++d) jump[d] = qr[d]-ql[d];
}

static inline void
calc_first_order_update(int meqn, double dtdx,
  double * GKYL_RESTRICT q, const double * GKYL_RESTRICT amdq_r, const double * GKYL_RESTRICT apdq_l)
{
  for (int i=0; i<meqn; ++i)
    q[i] = q[i] - dtdx*(apdq_l[i] + amdq_r[i]);
}

static inline double
calc_cfla(int mwaves, double cfla, double dtdx, const double *s)
{
  double c = cfla;
  for (int i=0; i<mwaves; ++i)
    c = fmax(c, dtdx*fabs(s[i]));
  return c;
}

static inline double
wave_dot_prod(int meqn, const double * GKYL_RESTRICT wa, const double * GKYL_RESTRICT wb)
{
  double dot = 0.0;
  for (int i=0; i<meqn; ++i) dot += wa[i]*wb[i];
  return dot;
}

static inline void
wave_rescale(int meqn, double fact, double *w)
{
  for (int i=0; i<meqn; ++i) w[i] *= fact; 
}

static inline void
calc_second_order_qflux(int meqn, double dtdx, double s,
  const double *waves, double * GKYL_RESTRICT flux2)
{
  double sfact = 0.5*fabs(s)*(1-fabs(s)*dtdx);
  for (int i=0; i<meqn; ++i)
    flux2[i] += sfact*waves[i];
}

// this is the sign function for doubles
static inline int sign_double(double val) { return (0.0 < val) - (val < 0.0); }

static inline void
calc_second_order_fflux(int meqn, double dtdx, double s,
  const double *waves, double * GKYL_RESTRICT flux2)
{
  double sfact = 0.5*sign_double(s)*(1-fabs(s)*dtdx);
  for (int i=0; i<meqn; ++i)
    flux2[i] += sfact*waves[i];
}

static inline void
calc_second_order_update(int meqn, double dtdx, double * GKYL_RESTRICT qout,
  const double *fl, const double *fr)
{
  for (int i=0; i<meqn; ++i)
    qout[i] += -dtdx*(fr[i]-fl[i]);
}

static void
limit_waves(const gkyl_wave_prop *wv, int mwaves,
  const struct gkyl_range *slice_range,
  int lower, int upper, struct gkyl_array *waves, const struct gkyl_array *speed)
{
  int meqn = wv->equation->num_equations;

  for (int mw=0; mw<mwaves; ++mw) {
    const double *wl = gkyl_array_cfetch(waves, gkyl_ridx(*slice_range, lower-1));
    const double *wr = gkyl_array_cfetch(waves, gkyl_ridx(*slice_range, lower));

    double dotr = wave_dot_prod(meqn, &wl[mw*meqn], &wr[mw*meqn]);

    for (int i=lower; i<=upper; ++i) {
      double dotl = dotr;
      
      double * GKYL_RESTRICT wi = gkyl_array_fetch(waves, gkyl_ridx(*slice_range, i));
      const double * GKYL_RESTRICT wi1 = gkyl_array_cfetch(waves, gkyl_ridx(*slice_range, i+1));
      
      double wnorm2 = wave_dot_prod(meqn, &wi[mw*meqn], &wi[mw*meqn]);
      dotr = wave_dot_prod(meqn, &wi[mw*meqn], &wi1[mw*meqn]);

      if (wnorm2 > 0) {
        const double *s = gkyl_array_cfetch(speed, gkyl_ridx(*slice_range, i));
        double r = s[mw] > 0 ? dotl/wnorm2 : dotr/wnorm2;
        double theta = limiter_function(r, wv->limiter);
        wave_rescale(meqn, theta, &wi[mw*meqn]);
      }
    }
  }
}

// advance method
struct gkyl_wave_prop_status
gkyl_wave_prop_advance(gkyl_wave_prop *wv,
  double tm, double dt, const struct gkyl_range *update_range,
  struct gkyl_array *phi, const struct gkyl_array *qin, struct gkyl_array *qout)
{
  wv->n_calls += 1;
  
  int ndim = update_range->ndim;
  int meqn = wv->equation->num_equations;
  //  when forced to use Lax fluxes, we only have a single wave
  int mwaves = wv->force_low_order_flux ? 2 :  wv->equation->num_waves;

  double cfla = 0.0, cfl = wv->cfl, cflm = 1.1*cfl;
  double is_cfl_violated = 0.0; // delibrately a double
  
  double ql_local[meqn], qr_local[meqn];
  double fjump_local[meqn];
  double waves_local[meqn*mwaves];
  double amdq_local[meqn], apdq_local[meqn];
  double delta[meqn];
  
  int idxl[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  double max_speed = 0.0;

  // state of the update
  enum update_state {
    WV_FIRST_SWEEP, WV_POSITIVITY_SWEEP, WV_FIN_SWEEP
  } state, next_state;

  for (int d=0; d<wv->num_up_dirs; ++d) {
    int dir = wv->update_dirs[d];

    double dtdx = dt/wv->grid.dx[dir];

    // upper/lower bounds in direction 'd'. These are edge indices
    int loidx = update_range->lower[dir]-1;
    int upidx = update_range->upper[dir]+2;

    // cell indices in 1D slice for interior cells
    int loidx_c = update_range->lower[dir];
    int upidx_c = update_range->upper[dir];

    struct gkyl_range slice_range;
    gkyl_range_init(&slice_range, 1, (int[]) { loidx }, (int[]) { upidx } );

    struct gkyl_range perp_range;
    gkyl_range_shorten_from_above(&perp_range, update_range, dir, 1);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &perp_range);

    // outer loop is over perpendicular directions, inner loop over 1D
    // slice along that direction
    while (gkyl_range_iter_next(&iter)) {
      
      gkyl_copy_int_arr(ndim, iter.idx, idxl);
      gkyl_copy_int_arr(ndim, iter.idx, idxr);

      gkyl_array_clear(wv->redo_fluct, 1.0);
      
      enum gkyl_wv_flux_type ftype = wv->force_low_order_flux ?
        GKYL_WV_LOW_ORDER_FLUX : GKYL_WV_HIGH_ORDER_FLUX;

      state = WV_FIRST_SWEEP;

      // perform 1D sweeps, fixing positivity if required
      while (state != WV_FIN_SWEEP) {

        if (state == WV_POSITIVITY_SWEEP)
          ftype = GKYL_WV_LOW_ORDER_FLUX;

        // copy previous time-step solution
        for (int i=loidx_c; i<=upidx_c; ++i) {
          idxl[dir] = i; // cell index
          long lidx = gkyl_range_idx(update_range, idxl);
          copy_wv_vec(meqn, gkyl_array_fetch(qout, lidx), gkyl_array_cfetch(qin, lidx));
        }

        for (int i=loidx; i<=upidx; ++i) {
          idxl[dir] = i-1; idxr[dir] = i;
          long sidx = gkyl_ridx(slice_range, i);

          const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxr);
          double *s = gkyl_array_fetch(wv->speeds, sidx);
          const double *redo_fluct = gkyl_array_cfetch(wv->redo_fluct, sidx);

          if (redo_fluct[0] > 0.0) {

            // compute fluctuations and waves only if needed (this
            // prevents doing the full 1D sweep with low-order fluxes
            // on positivity violations)
            long lidx = gkyl_range_idx(update_range, idxl);
            long ridx = gkyl_range_idx(update_range, idxr);

            const double *qinl = gkyl_array_cfetch(qin, lidx);
            const double *qinr = gkyl_array_cfetch(qin, ridx);

            const double *phil = gkyl_array_cfetch(phi, lidx);
            const double *phir = gkyl_array_cfetch(phi, ridx);

            gkyl_wv_eqn_rotate_to_local(wv->equation, cg->tau1[dir], cg->tau2[dir], cg->norm[dir], qinl, ql_local);
            gkyl_wv_eqn_rotate_to_local(wv->equation, cg->tau1[dir], cg->tau2[dir], cg->norm[dir], qinr, qr_local);

            if (wv->split_type == GKYL_WAVE_QWAVE)
              calc_jump(meqn, ql_local, qr_local, delta);
            else
              gkyl_wv_eqn_flux_jump(wv->equation, ql_local, qr_local, delta);

            double my_max_speed = gkyl_wv_eqn_waves(wv->equation, ftype, delta,
              ql_local, qr_local, phil[0], phir[0], waves_local, s);
            max_speed = max_speed > my_max_speed ? max_speed : my_max_speed;
            
            double lenr = cg->lenr[dir];
            for (int mw=0; mw<mwaves; ++mw)
              s[mw] *= lenr; // rescale speeds

            // compute fluctuations in local coordinates
            if (wv->split_type == GKYL_WAVE_QWAVE)
              gkyl_wv_eqn_qfluct(wv->equation, ftype, ql_local, qr_local,
                phil[0], phir[0], waves_local, s, amdq_local, apdq_local);
            else
              gkyl_wv_eqn_ffluct(wv->equation, ftype, ql_local, qr_local,
                phil[0], phir[0], waves_local, s, amdq_local, apdq_local);


            double *waves = gkyl_array_fetch(wv->waves, sidx);
            for (int mw=0; mw<mwaves; ++mw)
              // rotate waves back
              gkyl_wv_eqn_rotate_to_global(wv->equation, 
                cg->tau1[dir], cg->tau2[dir], cg->norm[dir], &waves_local[mw*meqn], &waves[mw*meqn]
              );

            // rotate fluctuations
            double *amdq = gkyl_array_fetch(wv->amdq, sidx);
            gkyl_wv_eqn_rotate_to_global(wv->equation, 
              cg->tau1[dir], cg->tau2[dir], cg->norm[dir], amdq_local, amdq);
            
            double *apdq = gkyl_array_fetch(wv->apdq, sidx);
           gkyl_wv_eqn_rotate_to_global(wv->equation, 
              cg->tau1[dir], cg->tau2[dir], cg->norm[dir], apdq_local, apdq);
          }
          
          cfla = calc_cfla(mwaves, cfla, dtdx/cg->kappa, s);
        }

        if (cfla > cflm) // check time-step before any updates are performed
          is_cfl_violated = 1.0;

        if (is_cfl_violated > 0) {
          // we need to use this goto to jump out of this deep loop to
          // avoid potential problems with taking too large a
          // time-step. NOTE: This jump is local to a rank. An
          // all-reduce here can't be done as one may end up with a
          // hang due to missing allreduce from some ranks.
          goto outsideloop;
        }

        // compute first-order update in each cell
        for (int i=loidx_c; i<=upidx_c; ++i) { // loop is over cells
        
          idxl[dir] = i; // cell index and left-edge index
          long lidx = gkyl_range_idx(update_range, idxl);

          const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxl);

          calc_first_order_update(meqn, dtdx/cg->kappa,
            gkyl_array_fetch(qout, lidx), 
            gkyl_array_cfetch(wv->amdq, gkyl_ridx(slice_range, i+1)),
            gkyl_array_cfetch(wv->apdq, gkyl_ridx(slice_range, i))
          );
        }

        if (state == WV_FIRST_SWEEP) {
          // we only compute second-correction if we are in first sweep
          
          // apply limiters to waves for all edges in update range,
          // including edges that are on the range boundary
          limit_waves(wv, mwaves, &slice_range,
            update_range->lower[dir], update_range->upper[dir]+1, wv->waves, wv->speeds);

          // get the kappa in the first ghost cell on left (needed in
          // the second order flux calculation)
          idxl[dir] = update_range->lower[dir]-1;
          const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxl);
          double kappal = cg->kappa;

          gkyl_array_clear(wv->flux2, 0.0);
          // compute second-order correction fluxes at each interface:
          // note that there is one extra edge than cell
          for (int i=loidx_c; i<=upidx_c+1; ++i) {
            long sidx = gkyl_ridx(slice_range, i);

            const double *waves = gkyl_array_cfetch(wv->waves, sidx);
            const double *s = gkyl_array_cfetch(wv->speeds, sidx);
            double *flux2 = gkyl_array_fetch(wv->flux2, sidx);

            idxl[dir] = i;
            const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxl);
            double kappar = cg->kappa;

            if (wv->split_type == GKYL_WAVE_QWAVE) {
              for (int mw=0; mw<mwaves; ++mw)
                calc_second_order_qflux(meqn, dtdx/(0.5*(kappal+kappar)), s[mw], &waves[mw*meqn], flux2);
            }
            else {
              for (int mw=0; mw<mwaves; ++mw)
                calc_second_order_fflux(meqn, dtdx/(0.5*(kappal+kappar)), s[mw], &waves[mw*meqn], flux2);
            }

            kappal = kappar;
          }

          // add second correction flux to solution in each interior cell
          for (int i=loidx_c; i<=upidx_c; ++i) {
            long sidx = gkyl_ridx(slice_range, i);

            idxl[dir] = i;
            const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxl);

            calc_second_order_update(meqn, dtdx/cg->kappa,
              gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl)),
              gkyl_array_cfetch(wv->flux2, gkyl_ridx(slice_range, i)),
              gkyl_array_cfetch(wv->flux2, gkyl_ridx(slice_range, i+1))
            );
          }
        }

        next_state = WV_FIN_SWEEP;
        // check invariant domains if needed
        if ( (state == WV_FIRST_SWEEP) && wv->check_inv_domain) {
          long n_bad_cells = 0;            

          gkyl_array_clear(wv->redo_fluct, 0.0); // by default no edge needs recomputing
          
          // check if invariant domains are violated, flagging edges
          // of each bad cell
          for (int i=loidx_c; i<=upidx_c; ++i) {
            idxl[dir] = i;
            const double *qt = gkyl_array_cfetch(qout, gkyl_range_idx(update_range, idxl));
            if (!gkyl_wv_eqn_check_inv(wv->equation, qt)) {

              double *redo_fluct_l = gkyl_array_fetch(wv->redo_fluct, gkyl_ridx(slice_range, i));
              double *redo_fluct_r = gkyl_array_fetch(wv->redo_fluct, gkyl_ridx(slice_range, i+1));
              // mark left and right edges so fluctuations are redone
              redo_fluct_l[0] = 1.0;
              redo_fluct_r[0] = 1.0;

              n_bad_cells += 1;
            }
          }

          if (n_bad_cells > 0) {
            // we need to resweep the 1D slice again
            next_state = WV_POSITIVITY_SWEEP;
            wv->n_bad_advance_calls += 1;
          }

          wv->n_bad_cells += n_bad_cells;
          wv->n_max_bad_cells = wv->n_max_bad_cells >  n_bad_cells ? wv->n_max_bad_cells : n_bad_cells;
        }

        if (wv->equation->type == GKYL_EQN_EULER_RGFM) {
          euler_rgfm_reinit_level_set(wv, update_range, idxl, loidx_c, upidx_c, qout, dir);
        }
        if (wv->equation->type == GKYL_EQN_GR_MAXWELL) {
          gr_maxwell_impose_gauge(wv, update_range, idxl, loidx_c, upidx_c, qout, dir);
        }
        if (wv->equation->type == GKYL_EQN_GR_MAXWELL_TETRAD) {
          gr_maxwell_tetrad_impose_gauge(wv, update_range, idxl, loidx_c, upidx_c, qout, dir);
        }
        if (wv->equation->type == GKYL_EQN_GR_EULER) {
          gr_euler_impose_gauge(wv, update_range, idxl, loidx_c, upidx_c, qout, dir);
        }
        if (wv->equation->type == GKYL_EQN_GR_EULER_TETRAD) {
          gr_euler_tetrad_impose_gauge(wv, update_range, idxl, loidx_c, upidx_c, qout, dir);
        }
        if (wv->equation->type == GKYL_EQN_GR_ULTRA_REL_EULER) {
          gr_ultra_rel_euler_impose_gauge(wv, update_range, idxl, loidx_c, upidx_c, qout, dir);
        }
        if (wv->equation->type == GKYL_EQN_GR_ULTRA_REL_EULER_TETRAD) {
          gr_ultra_rel_euler_tetrad_impose_gauge(wv, update_range, idxl, loidx_c, upidx_c, qout, dir);
        }
        if (wv->equation->type == GKYL_EQN_GR_TWOFLUID) {
          gr_twofluid_impose_gauge(wv, update_range, idxl, loidx_c, upidx_c, qout, dir);
        }
        if (wv->equation->type == GKYL_EQN_GR_TWOFLUID_TETRAD) {
          gr_twofluid_tetrad_impose_gauge(wv, update_range, idxl, loidx_c, upidx_c, qout, dir);
        }
        if (wv->equation->type == GKYL_EQN_GR_MHD) {
          gr_mhd_impose_gauge(wv, update_range, idxl, loidx_c, upidx_c, qout, dir);
        }
        if (wv->equation->type == GKYL_EQN_GR_MHD_TETRAD) {
          gr_mhd_tetrad_impose_gauge(wv, update_range, idxl, loidx_c, upidx_c, qout, dir);
        }

        state = next_state; // change state for next sweep
        
      } // end loop over sweeps
    } // end loop over perpendicular directions
  } // end loop over directions

  outsideloop:
  ;

  // compute actual CFL, status & max-speed across all domains
  double red_vars[3] = { cfla, is_cfl_violated, max_speed };
  double red_vars_global[3] = { 0.0, 0.0, 0.0 };
  gkyl_comm_allreduce(wv->comm, GKYL_DOUBLE, GKYL_MAX, 3, &red_vars, &red_vars_global);

  cfla = red_vars_global[0];
  is_cfl_violated = red_vars_global[1];
  max_speed = red_vars_global[2];

  double dt_suggested = dt*cfl/fmax(cfla, DBL_MIN);

  if (is_cfl_violated > 0.0)
    // indicate failure, and return smaller stable time-step
    return (struct gkyl_wave_prop_status) {
      .success = 0,
      .dt_suggested = dt_suggested,
      .max_speed = max_speed,
    };
  
  // on success, suggest only bigger time-step; (Only way dt can
  // reduce is if the update fails. If the code comes here the update
  // succeeded and so we should not allow dt to reduce).

  return (struct gkyl_wave_prop_status) {
    .success = is_cfl_violated > 0.0 ? 0 : 1,
    .dt_suggested = dt_suggested > dt ? dt_suggested : dt,
    .max_speed = max_speed,
  };
}

double
gkyl_wave_prop_max_dt(const gkyl_wave_prop *wv, const struct gkyl_range *update_range,
  const struct gkyl_array *qin)
{
  double max_dt = DBL_MAX;
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {

    for (int d=0; d<wv->num_up_dirs; ++d) {
      int dir = wv->update_dirs[d];
      double dx = wv->grid.dx[dir];

      const double *q = gkyl_array_cfetch(qin, gkyl_range_idx(update_range, iter.idx));
      double maxs = gkyl_wv_eqn_max_speed(wv->equation, q);
      max_dt = fmin(max_dt, wv->cfl*dx/maxs);
    }
    
  }

  return max_dt;
}

struct gkyl_wave_prop_stats
gkyl_wave_prop_stats(const gkyl_wave_prop *wv)
{
  return (struct gkyl_wave_prop_stats) {
    .n_calls = wv->n_calls,
    .n_bad_advance_calls = wv->n_bad_advance_calls,
    .n_bad_cells = wv->n_bad_cells,
    .n_max_bad_cells = wv->n_max_bad_cells
  };
}

void
gkyl_wave_prop_release(gkyl_wave_prop* up)
{
  gkyl_wv_eqn_release(up->equation);
  gkyl_array_release(up->waves);
  gkyl_array_release(up->apdq);
  gkyl_array_release(up->amdq);
  gkyl_array_release(up->speeds);
  gkyl_array_release(up->flux2);
  gkyl_array_release(up->redo_fluct);
  gkyl_comm_release(up->comm);
  
  gkyl_wave_geom_release(up->geom);
  
  gkyl_free(up);
}
