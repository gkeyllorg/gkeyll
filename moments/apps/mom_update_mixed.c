#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdint.h>

#include <gkyl_moment_priv.h>
#include <gkyl_gr_einstein_excision_priv.h>
#include <gkyl_sources_explicit_priv.h>
#include <gkyl_wv_gr_euler_priv.h>

// Per-species mixed-scheme time step: 
// - the GR-Euler fluid is advanced by the MP (method-of-lines SSP-RK3) 
// - the vacuum-Einstein spacetime by the wave-propagation
// ,coupled each step.
//

static bool
mixed_env_flag_enabled(const char *name)
{
  const char *v = getenv(name);
  return v != NULL && strcmp(v, "0") != 0 && strcmp(v, "false") != 0 &&
    strcmp(v, "False") != 0 && strcmp(v, "off") != 0 && strcmp(v, "OFF") != 0;
}

static bool
mixed_experimental_enabled(void)
{
  return mixed_env_flag_enabled("GKYL_EXPERIMENTAL");
}

static bool
mixed_experimental_flag_enabled(const char *name)
{
  return mixed_experimental_enabled() && mixed_env_flag_enabled(name);
}

static bool
mixed_force_split_path(void)
{
  return getenv("GKYL_MIXED_FORCE_SPLIT") != NULL
    || getenv("GKYL_WB_SPACETIME") != NULL
    || mixed_experimental_flag_enabled("GKYL_GAMMA_DRIVER")
    || mixed_env_flag_enabled("GKYL_SPACETIME_CELL_SPD_REPAIR");
}

static bool
mixed_dyn_collapse_hold_phase(const gkyl_moment_app *app)
{
  if (getenv("GKYL_DYN_COLLAPSE") == NULL) {
    return false;
  }

  // Stage times straddle the event. The driver's committed WB release, not
  // an RK stage timestamp, determines whether the whole step is held.
  for (int i=0; i<app->num_species; ++i)
    if (app->species[i].has_gr_euler)
      return !gkyl_gr_euler_wb_disabled(app->species[i].equation);
  return false;
}

static bool
mixed_all_species_wave_prop(const gkyl_moment_app *app)
{
  if (app->num_species == 0) {
    return false;
  }
  for (int i = 0; i < app->num_species; ++i) {
    if (app->species[i].scheme_type != GKYL_MOMENT_WAVE_PROP) {
      return false;
    }
  }
  return true;
}

static bool
mixed_no_species_wave_prop(const gkyl_moment_app *app)
{
  for (int i = 0; i < app->num_species; ++i) {
    if (app->species[i].scheme_type == GKYL_MOMENT_WAVE_PROP) {
      return false;
    }
  }
  return true;
}

// Inverse of a symmetric 3x3 matrix. Returns the determinant; fills inv.
static double
sym3_inverse(const double m[3][3], double inv[3][3])
{
  double c00 = m[1][1]*m[2][2] - m[1][2]*m[2][1];
  double c01 = m[1][2]*m[2][0] - m[1][0]*m[2][2];
  double c02 = m[1][0]*m[2][1] - m[1][1]*m[2][0];
  double det = m[0][0]*c00 + m[0][1]*c01 + m[0][2]*c02;
  double idet = 0.0;
  if (det != 0.0) {
    idet = 1.0/det;
  }
  inv[0][0] = c00*idet;
  inv[1][1] = (m[0][0]*m[2][2] - m[0][2]*m[2][0])*idet;
  inv[2][2] = (m[0][0]*m[1][1] - m[0][1]*m[1][0])*idet;
  inv[0][1] = inv[1][0] = (m[0][2]*m[2][1] - m[0][1]*m[2][2])*idet;
  inv[0][2] = inv[2][0] = (m[0][1]*m[1][2] - m[0][2]*m[1][1])*idet;
  inv[1][2] = inv[2][1] = (m[0][2]*m[1][0] - m[0][0]*m[1][2])*idet;
  return det;
}


static bool
mixed_spacetime_in_excision_region(double excision_threshold, bool is_conformal, const double *q)
{
  return gkyl_gr_einstein_is_excised(excision_threshold, q, is_conformal ? 77 : 64);
}

static void
mixed_zero_state(int meqn, double *q)
{
  for (int m = 0; m < meqn; ++m) {
    q[m] = 0.0;
  }
}


// ============================================================================
// Experimental connection-sourced (Gamma-driver) shift. Requires GKYL_EXPERIMENTAL=1 and GKYL_GAMMA_DRIVER.
// Without these settings it is never called. The latter's value is
// the drive strength mu_S (default 0.75); GKYL_GAMMA_DRIVER_ETA is the damping eta
// (default 1.0); GKYL_GAMMA_DRIVER_ALPHA_FLOOR (default = excision threshold) gates
// which cells are driven.
//
// For constant mu_S and eta, integrating beta_tt=mu_S*Gamma_t-eta*beta_t
// with zero initial shift rate gives beta_t=mu_S*(Gamma-Gamma0)-eta*(beta-beta0).
// Store Gamma0 and beta0 at the driver's first activation, including after WB
// release. Omitting this constant instead prescribes a nonzero initial rate.
// The reference is immutable through stages/retries and belongs to this species.
//
// D stores the conformal metric half-derivative. Its Christoffel contraction
// needs no additional finite-difference chi correction.
// Default is the integrated non-advective form (Alcubierre et al., gr-qc/0206072,
// Eqs. 43-46).
// Well-posedness of first-order Bona-Masso + driver is not proven here -- this is an experiment
// (cf. the first-order-reduction + moving-puncture-gauge machinery of FO-CCZ4, arXiv:1707.09910).
static bool
gamma_driver_cell_active(const double *q, double alpha_floor, bool is_conformal)
{
  if (!(isfinite(q[9]) && q[9] >= alpha_floor)) return false;
  if (is_conformal && !(isfinite(q[64]) && q[64] > 0.0)) return false;
  return true;
}

static bool mixed_newrad_enabled(const gkyl_moment_app *app);
static void apply_newrad_bc(gkyl_moment_app *app, struct moment_species *esp,
  const struct gkyl_array *fold, double dt, const struct gkyl_array *equilibrium,
  int first_component, int end_component);

static double
gamma_driver_alpha_floor(const struct moment_species *sp)
{
  const char *value = getenv("GKYL_GAMMA_DRIVER_ALPHA_FLOOR");
  return value ? atof(value) : sp->has_vacuum_einstein_conformal
    ? sp->vacuum_einstein_conformal_excision_threshold : sp->vacuum_einstein_excision_threshold;
}

static void
gamma_driver_advance_derivatives(gkyl_moment_app *app, struct moment_species *esp,
  const struct gkyl_array *fold, double tcurr, double alpha_floor)
{
  struct gkyl_array *sf=esp->f[0];
  bool conformal=esp->has_vacuum_einstein_conformal;
  moment_species_apply_bc(app,tcurr,esp,sf);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter,&app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc=gkyl_range_idx(&app->local_ext,iter.idx);
    double *q=gkyl_array_fetch(sf,loc);
    const double *old=gkyl_array_cfetch(fold,loc);
    if (!gamma_driver_cell_active(q,alpha_floor,conformal)) continue;
    for (int i=0; i<esp->ndim; ++i) {
      int ip[GKYL_MAX_DIM], im[GKYL_MAX_DIM];
      for (int d=0; d<esp->ndim; ++d) ip[d]=im[d]=iter.idx[d];
      ip[i]++; im[i]--;
      const double *qp=gkyl_array_cfetch(sf,gkyl_range_idx(&app->local_ext,ip));
      const double *qm=gkyl_array_cfetch(sf,gkyl_range_idx(&app->local_ext,im));
      const double *oldp=gkyl_array_cfetch(fold,gkyl_range_idx(&app->local_ext,ip));
      const double *oldm=gkyl_array_cfetch(fold,gkyl_range_idx(&app->local_ext,im));
      bool pv=gamma_driver_cell_active(qp,alpha_floor,conformal);
      bool mv=gamma_driver_cell_active(qm,alpha_floor,conformal);
      bool periodic=false;
      for (int d=0; d<app->num_periodic_dir; ++d)
        periodic=periodic || app->periodic_dirs[d]==i;
      int inward=0;
      if (!periodic && iter.idx[i]==app->global.lower[i]) inward=1;
      if (!periodic && iter.idx[i]==app->global.upper[i]) inward=-1;
      const double *q1=NULL, *q2=NULL, *old1=NULL, *old2=NULL;
      if (inward && app->local.upper[i]-app->local.lower[i]>=2) {
        ip[i]=iter.idx[i]+inward; im[i]=iter.idx[i]+2*inward;
        q1=gkyl_array_cfetch(sf,gkyl_range_idx(&app->local_ext,ip));
        q2=gkyl_array_cfetch(sf,gkyl_range_idx(&app->local_ext,im));
        old1=gkyl_array_cfetch(fold,gkyl_range_idx(&app->local_ext,ip));
        old2=gkyl_array_cfetch(fold,gkyl_range_idx(&app->local_ext,im));
        if (!gamma_driver_cell_active(q1,alpha_floor,conformal)
            || !gamma_driver_cell_active(q2,alpha_floor,conformal)) q1=q2=NULL;
      }
      for (int j=0; j<esp->ndim; ++j) {
        double b=q[52+j]-old[52+j];
        double bp=pv && isfinite(qp[52+j]) && isfinite(oldp[52+j]) ? qp[52+j]-oldp[52+j] : b;
        double bm=mv && isfinite(qm[52+j]) && isfinite(oldm[52+j]) ? qm[52+j]-oldm[52+j] : b;
        double derivative=pv && mv ? (bp-bm)/(2*app->grid.dx[i])
          : pv ? (bp-b)/app->grid.dx[i] : mv ? (b-bm)/app->grid.dx[i] : 0;
        // Differentiate the actual NewRad increment using the same one-sided
        // boundary stencil, not a separate ghost extrapolation of beta.
        if (q1 && isfinite(q1[52+j]) && isfinite(q2[52+j])
            && isfinite(old1[52+j]) && isfinite(old2[52+j]))
          derivative=inward*(-3*b+4*(q1[52+j]-old1[52+j])-(q2[52+j]-old2[52+j]))/(2*app->grid.dx[i]);
        // B_t = (1/2) partial_i beta_t. Preserve the initial discrete residual;
        // do not replace an analytic B by a mesh derivative when beta is fixed.
        q[55+3*i+j]=old[55+3*i+j]+0.5*derivative;
      }
    }
  }
  moment_species_apply_bc(app,tcurr,esp,sf);
}

static void
gamma_driver_shift_solve(gkyl_moment_app *app, struct moment_species *esp,
  const struct gkyl_array *fold, double tcurr, double dt)
{
  int ndim = esp->ndim;
  if (ndim < 2 || ndim > 3) return;

  const char *env = getenv("GKYL_GAMMA_DRIVER");
  if (!mixed_experimental_flag_enabled("GKYL_GAMMA_DRIVER")) return;

  // If a mid-run dynamical-collapse trigger is scheduled, keep the shift frozen (driver OFF) during the
  // HOLD phase t < trigger, so the held star's spacetime stays quasi-static; the driver activates only
  // once the star is released at the trigger.
  if (mixed_dyn_collapse_hold_phase(app)) return;

  double mu_S = atof(env);
  if (mu_S == 0.0) mu_S = 0.75;
  double eta = 1.0;
  const char *eta_env = getenv("GKYL_GAMMA_DRIVER_ETA");
  if (eta_env != NULL) eta = atof(eta_env);

  double alpha_floor = gamma_driver_alpha_floor(esp);

  bool is_conformal = esp->has_vacuum_einstein_conformal;
  struct gkyl_array *sf = esp->f[0];

  bool initialize_reference = !esp->gamma_driver_reference;
  if (initialize_reference) {
    esp->gamma_driver_reference = gkyl_array_new(GKYL_DOUBLE, 6, app->local_ext.volume);
    gkyl_array_clear(esp->gamma_driver_reference, 0.0);
    gkyl_moment_app_cout(app, stderr,
      "[GAMMA-REFERENCE] t=%.17g zero initial driver rate; reference captured at driver activation\n", tcurr);
  }


  // Pass 1: reconstruct Gammatilde^i locally and integrate the driver for beta (slots 52-54).
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    double *q = gkyl_array_fetch(sf, loc);
    bool active = gamma_driver_cell_active(q, alpha_floor, is_conformal);
    if (!active && !initialize_reference) continue;

    double gt[3][3], gti[3][3];
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        gt[i][j] = q[3*i + j];
    sym3_inverse(gt, gti);

    double Gam[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < 3; ++i) {
      double acc = 0.0;
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k) {
          double gjk = gti[j][k];
          if (gjk == 0.0) continue;
          for (int l = 0; l < 3; ++l) {
            double Djlk = q[19 + 9*j + 3*l + k];
            double Dklj = q[19 + 9*k + 3*l + j];
            double Dljk = q[19 + 9*l + 3*j + k];
            acc += gjk * gti[i][l] * (Djlk + Dklj - Dljk);
          }
        }
      Gam[i] = isfinite(acc) ? acc : 0.0;
    }


    // Capture even lapse-gated cells now, so later stages/retries cannot seed
    // a reference from a provisional evolved state.
    if (initialize_reference) {
      double *ref = gkyl_array_fetch(esp->gamma_driver_reference, loc);
      for (int i=0; i<3; ++i) {
        ref[i] = Gam[i];
        ref[3+i] = q[52+i];
      }
    }
    if (!active) continue;

    {
      double *ref = gkyl_array_fetch(esp->gamma_driver_reference, loc);
      for (int i = 0; i < ndim; ++i) {
        double b_old = isfinite(q[52+i]) ? q[52+i] : 0.0;
        double b_new = b_old + dt*(mu_S*(Gam[i]-ref[i]) - eta*(b_old-ref[3+i]));
        if (!isfinite(b_new)) b_new = b_old;
        q[52+i] = b_new;
      }
    }
  }

  // Use the actual driver rate for the radiative shift boundary, then derive B.
  // beta/B have zero Einstein flux and source; the completed transport substep
  // preserves them. This gauge subflow owns their evolution.
  if (mixed_newrad_enabled(app))
    apply_newrad_bc(app, esp, fold, dt, NULL, 52, 55);
  gamma_driver_advance_derivatives(app,esp,fold,tcurr,alpha_floor);

}

// SSP-RK2 for the integrated Gamma subflow, including its radiative boundary.
static void
gamma_driver_shift_rk2(gkyl_moment_app *app, struct moment_species *esp,
  double tcurr, double dt)
{
  if (esp->ndim<2 || esp->ndim>3 || !mixed_experimental_flag_enabled("GKYL_GAMMA_DRIVER")
      || mixed_dyn_collapse_hold_phase(app)) return;
  struct gkyl_array *old=gkyl_array_clone(esp->f[0]);
  gamma_driver_shift_solve(app,esp,old,tcurr,dt);
  struct gkyl_array *stage=gkyl_array_clone(esp->f[0]);
  gamma_driver_shift_solve(app,esp,stage,tcurr+dt,dt);
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it,&app->local);
  while (gkyl_range_iter_next(&it)) {
    long loc=gkyl_range_idx(&app->local_ext,it.idx);
    double *q=gkyl_array_fetch(esp->f[0],loc);
    const double *q0=gkyl_array_cfetch(old,loc);
    for (int i=0; i<esp->ndim; ++i) q[52+i]=0.5*(q0[52+i]+q[52+i]);
  }
  gamma_driver_advance_derivatives(app,esp,old,tcurr+dt,gamma_driver_alpha_floor(esp));
  gkyl_array_release(stage);
  gkyl_array_release(old);
}




// One SSP-RK3 forward-Euler stage of the spacetime source (h): gauge source (via the shared plus, if coupled, the frozen-fluid ADM matter source.
static void
spacetime_source_euler_stage(bool is_conformal, double excision_threshold,
  enum gkyl_spacetime_slicing slicing, enum gkyl_spacetime_evolution evolution,
  const double *qfluid, double gas_gamma, double tcurr, double h, int meqn,
  double *in, double *out, const double *wb_src)
{
  // Gauge source: out = in + h * S_gauge(in)
  if (is_conformal) {
    explicit_vacuum_einstein_conformal_source_update_euler(0, excision_threshold, slicing, evolution,
      tcurr, h, in, out);
  }
  else {
    explicit_vacuum_einstein_source_update_euler(0, excision_threshold, slicing, evolution,
      tcurr, h, in, out);
  }

  // Fluid -> spacetime matter source (only in live cells): out += h * M(in, fluid), M touching K_ij and V_k only
  if (qfluid &&
      !mixed_spacetime_in_excision_region(excision_threshold, is_conformal, out)) {
    double qfluid_stage[71];
    for (int m = 0; m < 71; m++) {
      qfluid_stage[m] = qfluid[m];
    }
    sync_fluid_metric_from_einstein(in, qfluid_stage, excision_threshold, is_conformal);

    double mrhs[77] = { 0.0 };
    add_einstein_matter_source(in, qfluid_stage, gas_gamma, excision_threshold, is_conformal, evolution, mrhs);
    for (int m = 0; m < meqn; ++m) {
      out[m] += h * mrhs[m];
    }
  }

  if (wb_src != NULL &&
      !mixed_spacetime_in_excision_region(excision_threshold, is_conformal, in)) {
    for (int m = 0; m < meqn; ++m) {
      out[m] -= h * wb_src[m];
    }
  }

  if (mixed_spacetime_in_excision_region(excision_threshold, is_conformal, out)) {
    mixed_zero_state(meqn, out);
  }
}

// Integrate the spacetime source (gauge + optional frozen-fluid matter) over a half-step h, in place on the spacetime state array `f`, using SSP-RK3 with the per-cell forward-Euler stage above.
// [SPACETIME WB] If wb_src != NULL, integrate (S(q) - S_eq) instead of S(q), where S_eq = wb_src is the
// frozen per-cell discrete source the gauge produces on the analytic KS-TOV equilibrium. The subtraction is
// applied inside each SSP-RK3 Euler stage, so q_eq is an exact fixed point of the source update when S(q_eq)=S_eq.
void
moment_spacetime_source_advance(gkyl_moment_app *app, struct moment_species *einstein_sp,
  struct gkyl_array *f, const struct gkyl_array *fluid_f, double gas_gamma, double tcurr, double h,
  const struct gkyl_array *wb_src)
{
  bool is_conformal = einstein_sp->has_vacuum_einstein_conformal;
  int meqn = einstein_sp->num_equations;

  double excision_threshold = is_conformal ? einstein_sp->vacuum_einstein_conformal_excision_threshold
    : einstein_sp->vacuum_einstein_excision_threshold;
  enum gkyl_spacetime_slicing slicing = is_conformal ? einstein_sp->vacuum_einstein_conformal_spacetime_slicing
    : einstein_sp->vacuum_einstein_spacetime_slicing;
  enum gkyl_spacetime_evolution evolution = is_conformal ? einstein_sp->vacuum_einstein_conformal_spacetime_evolution
    : einstein_sp->vacuum_einstein_spacetime_evolution;
  if (is_conformal) {
    // Grid-wide stages are required: the differentiated chi source contributes
    // to P/H before the next K source is evaluated.
    moment_species_apply_bc(app, tcurr, einstein_sp, f);
    struct gkyl_array *old = gkyl_array_clone(f);
    struct gkyl_array *stage = gkyl_array_clone(f);
    const double times[3] = {tcurr, tcurr+h, tcurr+.5*h};
    const double weights[3] = {1, .25, 2.0/3.0};
    for (int rk=0; rk<3; ++rk) {
      struct gkyl_range_iter cells;
      gkyl_range_iter_init(&cells, &app->local);
      while (gkyl_range_iter_next(&cells)) {
        long loc = gkyl_range_idx(&app->local_ext, cells.idx);
        spacetime_source_euler_stage(true, excision_threshold, slicing, evolution,
          fluid_f ? gkyl_array_cfetch(fluid_f,loc) : NULL, gas_gamma, times[rk], h,
          meqn, gkyl_array_fetch(stage,loc), gkyl_array_fetch(f,loc),
          wb_src ? gkyl_array_cfetch(wb_src,loc) : NULL);
      }
      moment_species_apply_bc(app, times[rk]+h, einstein_sp, f);
      moment_species_conformal_derivatives_advance(app, einstein_sp, stage, f, times[rk]+h);
      gkyl_array_scale_range(f, weights[rk], &app->local);
      gkyl_array_accumulate_range(f, 1-weights[rk], old, &app->local);
      if (rk==2) {
        gkyl_range_iter_init(&cells, &app->local);
        while (gkyl_range_iter_next(&cells)) {
          double *q=gkyl_array_fetch(f,gkyl_range_idx(&app->local_ext,cells.idx));
          if (mixed_spacetime_in_excision_region(excision_threshold,true,q)) mixed_zero_state(meqn,q);
        }
      }
      moment_species_apply_bc(app, rk==0 ? tcurr+h : rk==1 ? tcurr+.5*h : tcurr+h,
        einstein_sp, f);
      if (rk<2) gkyl_array_copy_range(stage, f, &app->local_ext);
    }
    gkyl_array_release(old); gkyl_array_release(stage);
    return;
  }
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    double *q = gkyl_array_fetch(f, loc);
    const double *qfluid = fluid_f ? gkyl_array_cfetch(fluid_f, loc) : NULL;

    double f_old[77], f_stage[77], f_tmp[77];
    for (int m = 0; m < meqn; ++m) {
      f_old[m] = q[m];
    }

    // stage 1
    const double *seq = wb_src ? gkyl_array_cfetch(wb_src, loc) : NULL;
    spacetime_source_euler_stage(is_conformal, excision_threshold, slicing, evolution, qfluid,
      gas_gamma, tcurr, h, meqn, f_old, f_stage, seq);

    // stage 2: 3/4 old + 1/4 euler(stage1)
    spacetime_source_euler_stage(is_conformal, excision_threshold, slicing, evolution, qfluid,
      gas_gamma, tcurr + h, h, meqn, f_stage, f_tmp, seq);
    for (int m = 0; m < meqn; ++m) {
      f_stage[m] = (0.75 * f_old[m]) + (0.25 * f_tmp[m]);
    }

    // stage 3: 1/3 old + 2/3 euler(stage2)
    spacetime_source_euler_stage(is_conformal, excision_threshold, slicing, evolution, qfluid,
      gas_gamma, tcurr + (0.5 * h), h, meqn, f_stage, f_tmp, seq);
    for (int m = 0; m < meqn; ++m) {
      q[m] = ((1.0 / 3.0) * f_old[m]) + ((2.0 / 3.0) * f_tmp[m]);
    }
    if (mixed_spacetime_in_excision_region(excision_threshold, is_conformal, q)) {
      mixed_zero_state(meqn, q);
    }
  }

}

// [SPACETIME WB] Compute S_eq per cell: the discrete forward-Euler gauge source the spacetime scheme
// produces on the analytic KS-TOV equilibrium (esp->f[0] at t=0) with the frozen equilibrium fluid. In the
// continuum this source is exactly cancelled by transport to keep the star stationary; discretely it is the
// O(h) truncation that drives the KS drift. Stored once and subtracted every step (well-balancing). Computed
// over app->local (the range the source half-step touches); excised cells get 0.
static void
compute_spacetime_wb_src(gkyl_moment_app *app, struct moment_species *esp,
  const struct gkyl_array *eq_st, const struct gkyl_array *fluid_f, double gas_gamma, double tcurr,
  struct gkyl_array *wb_src)
{
  bool is_conformal = esp->has_vacuum_einstein_conformal;
  int meqn = esp->num_equations;
  double excision_threshold = is_conformal ? esp->vacuum_einstein_conformal_excision_threshold
    : esp->vacuum_einstein_excision_threshold;
  enum gkyl_spacetime_slicing slicing = is_conformal ? esp->vacuum_einstein_conformal_spacetime_slicing
    : esp->vacuum_einstein_spacetime_slicing;
  enum gkyl_spacetime_evolution evolution = is_conformal ? esp->vacuum_einstein_conformal_spacetime_evolution
    : esp->vacuum_einstein_spacetime_evolution;
  double rate_dt = 1.0e-6;
  const char *rate_dt_env = getenv("GKYL_SPACETIME_WB_SRC_RATE_DT");
  if (rate_dt_env != NULL) {
    double v = atof(rate_dt_env);
    if (v > 0.0) {
      rate_dt = v;
    }
  }

  gkyl_array_clear(wb_src, 0.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *qeq = gkyl_array_cfetch(eq_st, loc);
    double *seq = gkyl_array_fetch(wb_src, loc);
    if (mixed_spacetime_in_excision_region(excision_threshold, is_conformal, qeq)) {
      continue; // excised: leave S_eq = 0 (source is zeroed there anyway)
    }
    const double *qfluid = fluid_f ? gkyl_array_cfetch(fluid_f, loc) : NULL;
    double in[77], out[77];
    for (int m = 0; m < meqn; ++m) in[m] = qeq[m];
    // Extract the source *rate* using a tiny Euler probe. Using h=1 here is not
    // safe because the vacuum-Einstein source updater includes excision cleanup,
    // so a diagnostic rate evaluation can cross a guard and store a clipped
    // finite update instead of S(q_eq).
    spacetime_source_euler_stage(is_conformal, excision_threshold, slicing, evolution, qfluid,
      gas_gamma, tcurr, rate_dt, meqn, in, out, NULL);
    for (int m = 0; m < meqn; ++m) seq[m] = (out[m] - in[m]) / rate_dt;
  }
}

bool
moment_spacetime_set_equilibrium(gkyl_moment_app *app, struct moment_species *esp,
  const struct gkyl_array *fluid, double gas_gamma)
{
  if (!esp->has_vacuum_einstein_conformal || esp->scheme_type != GKYL_MOMENT_WAVE_PROP
      || esp->vacuum_einstein_conformal_excision_threshold != 0
      || esp->spacetime_equilibrium || !fluid) return false;
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local_ext);
  while (gkyl_range_iter_next(&it)) {
    const double *q=gkyl_array_cfetch(esp->f[0],gkyl_range_idx(&app->local_ext,it.idx));
    if (!(q[9]>0 && q[64]>0) || !gkyl_wv_eqn_check_inv(esp->equation,q)) return false;
  }
  esp->spacetime_equilibrium=gkyl_array_clone(esp->f[0]);
  esp->spacetime_equilibrium_source=gkyl_array_clone(esp->f[0]);
  gkyl_array_clear(esp->spacetime_equilibrium_source,0);
  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    long loc=gkyl_range_idx(&app->local_ext,it.idx);
    const double *q=gkyl_array_cfetch(esp->spacetime_equilibrium,loc);
    double *s=gkyl_array_fetch(esp->spacetime_equilibrium_source,loc), f[71];
    memcpy(f,gkyl_array_cfetch(fluid,loc),sizeof f);
    sync_fluid_metric_from_einstein(q,f,0,true);
    gkyl_wv_eqn_source(esp->equation,q,s);
    add_einstein_matter_source(q,f,gas_gamma,0,true,
      esp->vacuum_einstein_conformal_spacetime_evolution,s);
  }
  return true;
}

static void
subtract_spacetime_wb_delta(gkyl_moment_app *app, struct moment_species *esp,
  const struct gkyl_array *delta)
{
  int meqn = esp->num_equations;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    double *q = gkyl_array_fetch(esp->f[0], loc);
    const double *dq = gkyl_array_cfetch(delta, loc);
    for (int m = 0; m < meqn; ++m) {
      q[m] -= dq[m];
    }
  }
}

// beta/B are stationary coefficients of the Einstein transport subflow in both
// layouts. Preserve them instead of retaining HLL/Lax numerical diffusion.
// This is the identity update (zero effective face flux), not B <- grad(beta).
static void
keep_shift_through_transport(gkyl_moment_app *app, struct moment_species *esp)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *old = gkyl_array_cfetch(esp->f[0], loc);
    double *next = gkyl_array_fetch(esp->f[app->ndim], loc);
    for (int m = 52; m < 64; ++m) next[m] = old[m];
  }
}

// Copy the vacuum-Einstein metric into the fluid's metric slots
static void
mixed_sync_metric(gkyl_moment_app *app, int fluid_idx, int einstein_idx)
{
  bool is_conformal = app->species[einstein_idx].has_vacuum_einstein_conformal;
  double excision_threshold = is_conformal ? app->species[einstein_idx].vacuum_einstein_conformal_excision_threshold
    : app->species[einstein_idx].vacuum_einstein_excision_threshold;

  const struct gkyl_array *einstein_f = app->species[einstein_idx].fcurr;
  struct gkyl_array *fluid_f = app->species[fluid_idx].fcurr;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *qe = gkyl_array_cfetch(einstein_f, loc);
    double *qf = gkyl_array_fetch(fluid_f, loc);
    sync_fluid_metric_from_einstein(qe, qf, excision_threshold, is_conformal);
  }
}

static bool
mixed_spacetime_cell_spd_repair_enabled(void)
{
  return mixed_env_flag_enabled("GKYL_SPACETIME_CELL_SPD_REPAIR");
}

static double
mixed_spacetime_cell_spd_repair_max_rel(void)
{
  const char *v = getenv("GKYL_SPACETIME_CELL_SPD_REPAIR_MAX_REL");
  double max_rel = v ? atof(v) : 5.0e-2;
  return max_rel > 0.0 ? max_rel : 5.0e-2;
}

static bool
mixed_metric_symmetrized_spd(const double g[9])
{
  double gs[9];
  for (int m = 0; m < 9; ++m) gs[m] = g[m];
  gs[1] = gs[3] = 0.5 * (g[1] + g[3]);
  gs[2] = gs[6] = 0.5 * (g[2] + g[6]);
  gs[5] = gs[7] = 0.5 * (g[5] + g[7]);
  return gkyl_mp_scheme_gr_metric_is_spd(gs);
}

// [CELL-SPD-REPAIR] Metric cell-average analog of the GR-Euler face SPD limiter.
// If a non-excised cell-average metric leaves the SPD cone by a small amount, blend
// the metric slots toward the previous accepted cell metric. This preserves the
// live lapse/shift/K/matter fields and fixes benign O(h) surface violations without
// turning them into time-step rejection. Gross/nonfinite states are still rejected.
static bool
mixed_spacetime_cell_spd_repair(gkyl_moment_app *app, struct moment_species *esp,
  const struct gkyl_array *anchor)
{
  if (!mixed_spacetime_cell_spd_repair_enabled()) {
    return true;
  }

  bool is_conformal = esp->has_vacuum_einstein_conformal;
  double excision_threshold = is_conformal ? esp->vacuum_einstein_conformal_excision_threshold
    : esp->vacuum_einstein_excision_threshold;
  double max_rel_allowed = mixed_spacetime_cell_spd_repair_max_rel();
  struct gkyl_array *ef = esp->fcurr;

  int rejected_local = 0;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    double *q = gkyl_array_fetch(ef, loc);
    const double *qa = gkyl_array_cfetch(anchor, loc);

    for (int m = 0; m < esp->num_equations; ++m) {
      if (!isfinite(q[m])) {
        rejected_local = 1;
        goto reject_cell;
      }
    }
    if (mixed_spacetime_in_excision_region(excision_threshold, is_conformal, q)) {
      continue;
    }

    double g[9], ga[9];
    bool metric_finite = true, anchor_finite = true;
    for (int m = 0; m < 9; ++m) {
      g[m] = q[m];
      ga[m] = qa[m];
      metric_finite = metric_finite && isfinite(g[m]);
      anchor_finite = anchor_finite && isfinite(ga[m]);
    }
    if (!metric_finite) {
      rejected_local = 1;
      goto reject_cell;
    }
    if (!anchor_finite || !mixed_metric_symmetrized_spd(ga)) {
      rejected_local = 1;
      goto reject_cell;
    }

    double scale = fmax(fmax(fabs(g[0]), fabs(g[4])), fmax(fabs(g[8]), 1.0e-300));
    double asym = fmax(fmax(fabs(g[1] - g[3]), fabs(g[2] - g[6])), fabs(g[5] - g[7]));
    if (asym > max_rel_allowed * scale) {
      rejected_local = 1;
      goto reject_cell;
    }

    g[1] = g[3] = 0.5 * (g[1] + g[3]);
    g[2] = g[6] = 0.5 * (g[2] + g[6]);
    g[5] = g[7] = 0.5 * (g[5] + g[7]);
    ga[1] = ga[3] = 0.5 * (ga[1] + ga[3]);
    ga[2] = ga[6] = 0.5 * (ga[2] + ga[6]);
    ga[5] = ga[7] = 0.5 * (ga[5] + ga[7]);

    if (gkyl_mp_scheme_gr_metric_is_spd(g)) {
      continue;
    }

    double lo = 0.0, hi = 1.0;
    for (int it = 0; it < 50; ++it) {
      double mid = 0.5 * (lo + hi);
      double gm[9];
      for (int m = 0; m < 9; ++m) {
        gm[m] = ga[m] + mid * (g[m] - ga[m]);
      }
      if (gkyl_mp_scheme_gr_metric_is_spd(gm)) {
        lo = mid;
      }
      else {
        hi = mid;
      }
    }

    double theta = (1.0 - 1.0e-12) * lo;
    double grepaired[9];
    double rel = 0.0;
    for (int m = 0; m < 9; ++m) {
      grepaired[m] = ga[m] + theta * (g[m] - ga[m]);
      rel = fmax(rel, fabs(grepaired[m] - q[m]) / scale);
    }
    if (rel > max_rel_allowed || !gkyl_mp_scheme_gr_metric_is_spd(grepaired)) {
      rejected_local = 1;
      goto reject_cell;
    }

    for (int m = 0; m < 9; ++m) {
      q[m] = grepaired[m];
    }
    continue;

reject_cell:
    break;
  }

  int rejected_global = rejected_local;
  gkyl_comm_allreduce(app->comm, GKYL_INT, GKYL_MAX, 1, &rejected_local, &rejected_global);


  return rejected_global == 0;
}



// Hyperbolicity guard: are all non-excised cell-average spacetime metrics valid inputs for the fluid solve?
// A cell is invalid if its spatial metric (slots 0-8) is (a) grossly asymmetric (numerical drift is fine - the SPD test below uses the symmetrized metric, but a large antisymmetric part signals breakdown);
// (b) not SPD (positive leading principal minors of the symmetrized metric, scale-relative floors);
// or, for the conformal system, (c) has a non-positive BSSN conformal factor chi (slot 64).
static bool
mixed_spacetime_metric_valid(gkyl_moment_app *app, const struct moment_species *esp)
{
  bool is_conformal = esp->has_vacuum_einstein_conformal;
  double excision_threshold;
  if (is_conformal) {
    excision_threshold = esp->vacuum_einstein_conformal_excision_threshold;
  }
  else {
    excision_threshold = esp->vacuum_einstein_excision_threshold;
  }
  const struct gkyl_array *ef = esp->fcurr;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  int local_bad = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *e = gkyl_array_cfetch(ef, loc);
    for (int m = 0; m < esp->num_equations; ++m) {
      if (!isfinite(e[m])) {
        local_bad = 1;
        break;
      }
    }
    if (local_bad) {
      break;
    }
    if (mixed_spacetime_in_excision_region(excision_threshold, is_conformal, e)) {
      continue; // excised cell: handled by excision, not the SPD guard
    }

    if (!(isfinite(e[0]) && isfinite(e[1]) && isfinite(e[2]) &&
          isfinite(e[3]) && isfinite(e[4]) && isfinite(e[5]) &&
          isfinite(e[6]) && isfinite(e[7]) && isfinite(e[8]))) {
      local_bad = 1;
      break;
    }

    double scale = fmax(fmax(fabs(e[0]), fabs(e[4])), fmax(fabs(e[8]), 1.0e-300));
    // (a) reject only gross asymmetry
    double asym = fmax(fmax(fabs(e[1] - e[3]), fabs(e[2] - e[6])), fabs(e[5] - e[7]));
    if (asym > 1.0e-2 * scale) {
      local_bad = 1;
      break;
    }
    // (b) SPD via leading principal minors (Sylvester) on the symmetrized metric (Sylvester needs symmetry), scale-relative floors (minors scale as scale^k)
    double g00 = e[0], g11 = e[4], g22 = e[8];
    double g01 = 0.5 * (e[1] + e[3]), g02 = 0.5 * (e[2] + e[6]), g12 = 0.5 * (e[5] + e[7]);
    double ep = 1.0e-12;
    double m1 = g00;
    double m2 = (g00 * g11) - (g01 * g01);
    double m3 = (g00 * ((g11 * g22) - (g12 * g12))) - (g01 * ((g01 * g22) - (g12 * g02))) +
      (g02 * ((g01 * g12) - (g11 * g02)));
    if (!(m1 > ep * scale && m2 > ep * scale * scale && m3 > ep * scale * scale * scale)) {
      local_bad = 1;
      break;
    }
    // (c) conformal: the BSSN conformal factor chi must be positive
    if (is_conformal && !(isfinite(e[64]) && e[64] > 0.0)) {
      local_bad = 1;
      break;
    }
  }


  int global_bad = local_bad;
  gkyl_comm_allreduce(app->comm, GKYL_INT, GKYL_MAX, 1, &local_bad, &global_bad);
  return (global_bad == 0);
}

// Factor E^2-D^2 using E=D+tau: cold atmosphere loses its pressure margin
// when the two rest-mass energy squares are rounded before subtraction.
static double
mixed_gr_euler_momentum_capacity(const double q[5])
{
  return q[4] * (2.0 * q[0] + q[4]);
}

// Invariant-domain guard for the GR-Euler cell averages after the fluid RK step.
// The MP reconstruction PCP only controls face states/fluxes; the geometric
// source update can still push a cell average outside the pressure-positive Valencia cone. In the
// coupled mixed step, catch that before the matter source feeds the bad fluid
// stress into K_ij.
static double mixed_gr_euler_density_marker(const struct moment_species *fsp, const double *q,
  double sqrt_det);

// Atmosphere band: cells whose density marker is at most this multiple of rho_atm.
static bool route1_fx_masked(long l);
static const double mixed_atm_band = 100.0;

static bool
mixed_gr_euler_state_valid(gkyl_moment_app *app, const struct moment_species *fsp)
{
  if (!fsp->has_gr_euler) {
    return true;
  }
  const double rho_atm = container_of(fsp->equation, struct wv_gr_euler, eqn)->rho_atm;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  int local_bad = 0;
  double bad_det = 0.0, bad_de = 0.0, bad_s2 = 0.0, bad_cone = 0.0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *q = gkyl_array_cfetch(fsp->fcurr, loc);
    if (route1_fx_masked(loc)) continue; // [FX] excised fluid cell

    for (int m = 0; m < fsp->num_equations; ++m) {
      if (!isfinite(q[m])) {
        local_bad = 1;
        break;
      }
    }
    if (!local_bad && q[27] >= pow(10.0, -8.0)) {
      if (!(q[0] > 0.0)) {
        local_bad = 1;
      }
      else {
        double g[3][3] = {
          { q[9], q[10], q[11] },
          { q[12], q[13], q[14] },
          { q[15], q[16], q[17] } };
        double gi[3][3];
        bad_det = sym3_inverse(g, gi);
        if (!(isfinite(bad_det) && bad_det > 0.0)) {
          local_bad = 1;
        }
        else {
          double S[3] = { q[1], q[2], q[3] };
          bad_de = q[0] + q[4];
          bad_s2 = 0.0;
          for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
              bad_s2 += gi[i][j] * S[i] * S[j];
          bad_cone = mixed_gr_euler_momentum_capacity(q) - bad_s2;
          // In the atmosphere band only the cone inequality is relaxed: the primitive-recovery floor gives the
          // fluxes and sources a static atmosphere there, and the stored state is not rewritten. Finiteness,
          // D > 0, det > 0 and E > 0 stay required in every cell.
          bool in_band = mixed_gr_euler_density_marker(fsp, q, sqrt(bad_det)) <= mixed_atm_band*rho_atm;
          if (!(bad_cone > 0.0 || in_band)) { // [BARO] (route-1 build, R1_BARO): the recovery's own test (fallback or live atmosphere)
            extern bool gkyl_gr_euler_c2p_baro_handled(double gas_gamma, const double q[71]);
            in_band = gkyl_gr_euler_c2p_baro_handled(fsp->gr_euler_gas_gamma, q);
          }
          if (!(isfinite(bad_de) && isfinite(bad_s2) && isfinite(bad_cone) &&
                bad_de > 0.0 && (bad_cone > 0.0 || in_band))) {
            local_bad = 1;
          }
        }
      }
    }

    if (local_bad) {
      break;
    }
  }


  int global_bad = local_bad;
  gkyl_comm_allreduce(app->comm, GKYL_INT, GKYL_MAX, 1, &local_bad, &global_bad);
  return (global_bad == 0);
}

static double
mixed_gr_euler_density_marker(const struct moment_species *fsp, const double *q,
  double sqrt_det)
{
  if (fsp->num_equations > 71 && isfinite(q[71]) && q[71] >= 0.0 &&
      sqrt_det > 0.0) {
    return q[71] / sqrt_det;
  }
  if (sqrt_det > 0.0) {
    return q[0] / sqrt_det;
  }
  return DBL_MAX;
}

// Flat (Minkowski) asymptotic value of vacuum-Einstein component c: gamma_ij = delta_ij, alpha = 1,
// all else 0. This is f_0 in the NewRad radiative BC (the value the field approaches at infinity).
// For the conformal (BSSN-chi, 77-component) system the asymptotic state is the SAME for the shared
// slots [0-63] -- the conformal metric gamma-tilde_ij -> delta_ij and the lapse alpha -> 1 -- plus the
// BSSN conformal factor chi = q[64] -> 1; its derivatives (first [65-67], second [68-76]) -> 0.
static double
newrad_flat_state(int c, bool is_conformal)
{
  if (c == 0 || c == 4 || c == 8) {
    return 1.0; // gamma_xx, gamma_yy, gamma_zz (conformal: gamma-tilde diagonal)
  }
  if (c == 9) {
    return 1.0; // lapse
  }
  if (is_conformal && c == 64) {
    return 1.0; // BSSN conformal factor chi
  }
  return 0.0;
}

// NweRad Asymptotic characteristic speed of vacuum-Einstein component c, for the radiative BC (Alcubierre gr-qc/0206072)
// Same rule for conformal (BSSN-chi): the gauge sector (lapse [9], lapse-derivative A_k [46-48])
// propagates at the gauge speed sqrt(f); every metric-sector field -- including the conformal factor
// chi [64] and its derivatives -- radiates at the (asymptotic) light speed 1, since gamma-tilde -> delta
// and chi -> 1 make the physical metric flat at the outer boundary.
static double
newrad_speed(int c, double alpha, enum gkyl_spacetime_slicing slicing)
{
  bool gauge = (c == 9 || c == 46 || c == 47 || c == 48);
  if (!gauge) {
    return 1.0;
  }
  double f = 1.0;
  if (slicing == GKYL_1PLUSLOG_SLICING) {
    f = (alpha > 1.0e-12) ? 2.0 / alpha : 2.0;
  }
  // harmonic and geodesic: f = 1.
  return sqrt(f);
}

// 2nd-order radial derivative (x^i/r) d_i f_c at grid index `idx`
static double
newrad_radial_deriv(const struct gkyl_array *arr, const struct gkyl_range *local,
  const struct gkyl_range *local_ext, const struct gkyl_rect_grid *grid,
  const int *idx, const double *xc, double r, int c, int ndim)
{
  double dr = 0.0;
  for (int d = 0; d < ndim; d++) {
    int ip[GKYL_MAX_DIM], im[GKYL_MAX_DIM], ipp[GKYL_MAX_DIM], imm[GKYL_MAX_DIM];
    for (int e = 0; e < ndim; e++) {
      ip[e] = im[e] = ipp[e] = imm[e] = idx[e];
    }
    ip[d] = idx[d] + 1; im[d] = idx[d] - 1;
    ipp[d] = idx[d] + 2; imm[d] = idx[d] - 2;
    double f0 = ((const double *) gkyl_array_cfetch(arr, gkyl_range_idx(local_ext, idx)))[c];
    double di;
    if (idx[d] == local->upper[d]) {
      // inward = -d, 2nd-order backward
      double fm1 = ((const double *) gkyl_array_cfetch(arr, gkyl_range_idx(local_ext, im)))[c];
      double fm2 = ((const double *) gkyl_array_cfetch(arr, gkyl_range_idx(local_ext, imm)))[c];
      di = (3.0 * f0 - 4.0 * fm1 + fm2) / (2.0 * grid->dx[d]);
    }
    else if (idx[d] == local->lower[d]) {
      // inward = +d, 2nd-order forward
      double fp1 = ((const double *) gkyl_array_cfetch(arr, gkyl_range_idx(local_ext, ip)))[c];
      double fp2 = ((const double *) gkyl_array_cfetch(arr, gkyl_range_idx(local_ext, ipp)))[c];
      di = (-3.0 * f0 + 4.0 * fp1 - fp2) / (2.0 * grid->dx[d]);
    }
    else {
      double fp1 = ((const double *) gkyl_array_cfetch(arr, gkyl_range_idx(local_ext, ip)))[c];
      double fm1 = ((const double *) gkyl_array_cfetch(arr, gkyl_range_idx(local_ext, im)))[c];
      di = (fp1 - fm1) / (2.0 * grid->dx[d]);
    }
    dr += (xc[d] / r) * di;
  }
  return dr;
}

static bool
newrad_live_spacetime_state_valid(const double *q, int meqn, bool is_conformal,
  double excision_threshold)
{
  for (int c = 0; c < meqn; ++c) {
    if (!isfinite(q[c])) {
      return false;
    }
  }

  if (mixed_spacetime_in_excision_region(excision_threshold, is_conformal, q)) {
    return false;
  }
  if (is_conformal && !gkyl_gr_conformal_scalars_valid(q)) return false;

  double scale = fmax(fmax(fabs(q[0]), fabs(q[4])), fmax(fabs(q[8]), 1.0e-300));
  double asym = fmax(fmax(fabs(q[1] - q[3]), fabs(q[2] - q[6])), fabs(q[5] - q[7]));
  if (asym > 1.0e-2 * scale) {
    return false;
  }

  double g00 = q[0], g11 = q[4], g22 = q[8];
  double g01 = 0.5 * (q[1] + q[3]);
  double g02 = 0.5 * (q[2] + q[6]);
  double g12 = 0.5 * (q[5] + q[7]);
  double ep = 1.0e-12;
  double m1 = g00;
  double m2 = (g00 * g11) - (g01 * g01);
  double m3 = (g00 * ((g11 * g22) - (g12 * g12))) -
    (g01 * ((g01 * g22) - (g12 * g02))) +
    (g02 * ((g01 * g12) - (g11 * g02)));

  return m1 > ep * scale && m2 > ep * scale * scale &&
    m3 > ep * scale * scale * scale;
}

// NewRad Radiative (Sommerfeld) outer boundary condition, Alcubierre et al. gr-qc/0206072 Eqs. (108)-(113) / ETK NewRad.
// Ansatz f = f_0 + u(r - v t)/r + h(t)/r^n gives the differential BC d_t f = -v (x^i/r) d_i f - v (f - f_0)/r + h'(t)/r^n .
static bool
mixed_newrad_enabled(const gkyl_moment_app *app)
{
  int comm_sz = 1;
  gkyl_comm_get_size(app->comm, &comm_sz);
  // The current stencil addresses local boundary cells, so it is serial-only.
  // Gauge choice must not silently change the outer boundary condition.
  return comm_sz == 1 && getenv("GKYL_NO_NEWRAD") == NULL;
}

static void
apply_newrad_bc(gkyl_moment_app *app, struct moment_species *esp,
  const struct gkyl_array *fold, double dt, const struct gkyl_array *equilibrium,
  int first_component, int end_component)
{
  bool is_conformal = esp->has_vacuum_einstein_conformal;
  int meqn = esp->num_equations;
  enum gkyl_spacetime_slicing slicing = is_conformal ? esp->vacuum_einstein_conformal_spacetime_slicing
    : esp->vacuum_einstein_spacetime_slicing;
  double excision_threshold = is_conformal ? esp->vacuum_einstein_conformal_excision_threshold
    : esp->vacuum_einstein_excision_threshold;
  double nexp = 3.0;
  const char *n_env = getenv("GKYL_NEWRAD_N");
  if (n_env) {
    nexp = atof(n_env);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    bool on_boundary = false;
    for (int d = 0; d < app->ndim; d++) {
      if (iter.idx[d] == app->local.lower[d] || iter.idx[d] == app->local.upper[d]) {
        on_boundary = true;
      }
    }
    if (!on_boundary) {
      continue;
    }

    double x0[3] = { 0.0, 0.0, 0.0 };
    gkyl_rect_grid_cell_center(&app->grid, iter.idx, x0);
    double r0 = sqrt(x0[0]*x0[0] + x0[1]*x0[1] + x0[2]*x0[2]);
    if (r0 < 1.0e-12) {
      continue;
    }

    // p1: one cell inward along every boundary direction (diagonal for edge/corner) -- always interior.
    int p1[GKYL_MAX_DIM];
    for (int d = 0; d < app->ndim; d++) {
      p1[d] = iter.idx[d];
      if (iter.idx[d] == app->local.upper[d]) {
        p1[d] = iter.idx[d] - 1;
      }
      else if (iter.idx[d] == app->local.lower[d]) {
        p1[d] = iter.idx[d] + 1;
      }
    }
    double x1[3] = { 0.0, 0.0, 0.0 };
    gkyl_rect_grid_cell_center(&app->grid, p1, x1);
    double r1 = sqrt(x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2]);

    long loc0 = gkyl_range_idx(&app->local_ext, iter.idx);
    long loc1 = gkyl_range_idx(&app->local_ext, p1);
    const double *f0_old = gkyl_array_cfetch(fold, loc0);
    const double *f1_old = gkyl_array_cfetch(fold, loc1);
    const double *f1_new = gkyl_array_cfetch(esp->f[0], loc1);
    const double *eq0 = equilibrium ? gkyl_array_cfetch(equilibrium, loc0) : NULL;
    const double *eq1 = equilibrium ? gkyl_array_cfetch(equilibrium, loc1) : NULL;
    double alpha0 = f0_old[9];
    double *out = gkyl_array_fetch(esp->f[0], loc0);
    double cand[77];
    memcpy(cand, out, meqn*sizeof(double));

    for (int c = first_component; c < end_component; c++) {
      // The held equilibrium is stationary: radiate its perturbation, not the
      // equilibrium's finite-radius residual. No subtraction when WB is off.
      double finf = equilibrium ? 0.0 : newrad_flat_state(c, is_conformal);
      double f0 = f0_old[c] - (eq0 ? eq0[c] : 0.0);
      double f1 = f1_old[c] - (eq1 ? eq1[c] : 0.0);
      double v = newrad_speed(c, alpha0, slicing);

      // Solve h'(t) from the measured interior rate at p1 (Alcubierre Eq. 112/113).
      double hprime = 0.0;
      if (r1 > 1.0e-12 && dt > 0.0) {
        double dt_f1 = (f1_new[c] - f1_old[c]) / dt;
        double drf1 = newrad_radial_deriv(fold, &app->local, &app->local_ext, &app->grid, p1, x1, r1, c, app->ndim);
        if (equilibrium)
          drf1 -= newrad_radial_deriv(equilibrium, &app->local, &app->local_ext, &app->grid, p1, x1, r1, c, app->ndim);
        double resid1 = dt_f1 + v * drf1 + v * (f1 - finf) / r1;
        hprime = pow(r1, nexp) * resid1;
      }

      double drf0 = newrad_radial_deriv(fold, &app->local, &app->local_ext, &app->grid, iter.idx, x0, r0, c, app->ndim);
      if (equilibrium)
        drf0 -= newrad_radial_deriv(equilibrium, &app->local, &app->local_ext, &app->grid, iter.idx, x0, r0, c, app->ndim);
      double rhs = -v * drf0 - v * (f0 - finf) / r0 + hprime / pow(r0, nexp);
      cand[c] = f0_old[c] + dt * rhs;
    }

    const double *accepted = cand;
    if (!newrad_live_spacetime_state_valid(cand, meqn, is_conformal, excision_threshold)) {
      if (newrad_live_spacetime_state_valid(f1_new, meqn, is_conformal, excision_threshold)) {
        accepted = f1_new;
      }
      else if (newrad_live_spacetime_state_valid(f0_old, meqn, is_conformal, excision_threshold)) {
        accepted = f0_old;
      }
    }

    for (int c = first_component; c < end_component; c++) {
      out[c] = accepted[c];
    }
  }

}

static void
apply_newrad_spacetime_bc(gkyl_moment_app *app, struct moment_species *esp,
  const struct gkyl_array *fold, double tnew, double dt, const struct gkyl_array *equilibrium)
{
  struct gkyl_array *boundary_old = NULL;
  if (esp->has_vacuum_einstein_conformal)
    boundary_old = gkyl_array_clone(esp->f[0]);

  // Active shift solvers own beta/B. P/H follow the radiative chi increment,
  // rather than independent radiation laws for its spatial derivatives.
  apply_newrad_bc(app, esp, fold, dt, equilibrium, 0, 52);
  if (boundary_old)
    apply_newrad_bc(app, esp, fold, dt, equilibrium, 64, 65);
  moment_species_apply_bc(app, tnew, esp, esp->f[0]);
  if (boundary_old) {
    moment_species_conformal_derivatives_advance(app, esp, boundary_old, esp->f[0], tnew);
    moment_species_apply_bc(app, tnew, esp, esp->f[0]);
    gkyl_array_release(boundary_old);
  }
}

static void
mixed_restore_step_backups(gkyl_moment_app *app, int fluid_idx, int einstein_idx,
  const struct gkyl_array *fluid_backup, const struct gkyl_array *spacetime_backup)
{
  if (fluid_backup && fluid_idx >= 0) {
    gkyl_array_copy_range(app->species[fluid_idx].fcurr, fluid_backup, &app->local_ext);
  }
  if (spacetime_backup && einstein_idx >= 0) {
    gkyl_array_copy_range(app->species[einstein_idx].fcurr, spacetime_backup, &app->local_ext);
  }
}

static void
mixed_release_step_backups(struct gkyl_array *fluid_backup,
  struct gkyl_array *spacetime_backup)
{
  if (fluid_backup) {
    gkyl_array_release(fluid_backup);
  }
  if (spacetime_backup) {
    gkyl_array_release(spacetime_backup);
  }
}

// [EXCISION TRANSPORT FREEZE] Excision only gates the SOURCE half-steps (spacetime_source_half_step /
// the GR-Euler matter source, both verified to freeze correctly). Transport (wave_prop's flux-divergence
// update, calc_second_order_update) has NO per-cell excision awareness at all -- it updates every cell
// in range unconditionally, for every species. So a cell that
// is excised going into a step still gets a real flux-divergence update from transport, silently undoing
// the source freeze -- confirmed as the root cause of the post-horizon mass growth and (via the
// lapse-driven fluid excision flag in sync_fluid_metric_from_einstein) the excision-boundary chatter seen
// in scan_tmp/P1_bufferfix_t16_x24. This supersedes the earlier face-level mirror knob
// (GKYL_EXCISION_FACE_BC, removed): that patched only the flux INTO a live cell, per species, behind an
// off-by-default env gate; this freezes the excised cell itself, for both species, unconditionally.
// Verified equivalent in physics to that knob (SumD, M0_in_AH agree to 4+ digits). This restores any cell that
// was ALREADY excised at the START of the step back to its pre-step value, using the backups already
// captured for step-rejection rollback (no extra allocation) -- so excision means frozen for the WHOLE
// step, not just against the source contribution. No-op for a cell that was not excised pre-step.
static void
mixed_restore_excised_spacetime_cells(gkyl_moment_app *app, struct moment_species *esp,
  const struct gkyl_array *spacetime_backup, double excision_threshold, bool is_conformal)
{
  if (!spacetime_backup) {
    return;
  }
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *qb = gkyl_array_cfetch(spacetime_backup, loc);
    if (mixed_spacetime_in_excision_region(excision_threshold, is_conformal, qb)) {
      double *q = gkyl_array_fetch(esp->f[0], loc);
      for (int m = 0; m < esp->num_equations; ++m) {
        q[m] = qb[m];
      }
    }
  }
}

static void
mixed_restore_excised_fluid_cells(gkyl_moment_app *app, struct moment_species *fsp,
  const struct gkyl_array *fluid_backup)
{
  if (!fluid_backup || fsp->num_equations <= 27) {
    return;
  }
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *qb = gkyl_array_cfetch(fluid_backup, loc);
    if (qb[27] < pow(10.0, -8.0)) { // matches explicit_gr_euler_source_update_euler's excision gate
      double *q = gkyl_array_fetch(fsp->fcurr, loc);
      for (int m = 0; m < fsp->num_equations; ++m) {
        q[m] = qb[m];
      }
    }
  }
}


#define ROUTE1_BARO // the route-1 build passed -DROUTE1_BARO (R1_BARO=1)
#include "route1_step.h"

// Take a single mixed-scheme time-step.
struct gkyl_update_status
moment_update_mixed(gkyl_moment_app *app, double dt0)
{
  struct gkyl_update_status st = { .success = true, .dt_actual = dt0, .dt_suggested = DBL_MAX };
  double tcurr = app->tcurr;

  // Identify the coupled GR-Euler (MP) / vacuum-Einstein (wave_prop) pair.
  int fluid_idx = -1, einstein_idx = -1;
  for (int i = 0; i < app->num_species; ++i) {
    if (app->species[i].has_gr_euler && app->species[i].scheme_type == GKYL_MOMENT_MP) {
      fluid_idx = i;
    }
    if ((app->species[i].has_vacuum_einstein || app->species[i].has_vacuum_einstein_conformal)
      && app->species[i].scheme_type == GKYL_MOMENT_WAVE_PROP) {
      einstein_idx = i;
    }
  }
  bool coupled = (fluid_idx >= 0 && einstein_idx >= 0);

  // With no actual mixed coupling/gauge work, use the standard updater for
  // that species family. This keeps vacuum-only and fluid-only control runs
  // comparable to their native GKYL_MOMENT_WAVE_PROP/GKYL_MOMENT_MP drivers,
  // while GKYL_MIXED_FORCE_SPLIT preserves the explicit mixed split for tests.
  if (!coupled && !mixed_force_split_path()) {
    if (mixed_all_species_wave_prop(app)) {
      return moment_update_one_step(app, dt0);
    }
    if (mixed_no_species_wave_prop(app)) {
      return moment_update_ssp_rk3(app, dt0);
    }
  }

  struct gkyl_array *fluid_step_backup = NULL;
  struct gkyl_array *spacetime_step_backup = NULL;
  if (fluid_idx >= 0) {
    struct moment_species *fsp = &app->species[fluid_idx];
    fluid_step_backup = gkyl_array_new(GKYL_DOUBLE, fsp->num_equations, app->local_ext.volume);
    gkyl_array_copy_range(fluid_step_backup, fsp->fcurr, &app->local_ext);
  }
  if (einstein_idx >= 0) {
    struct moment_species *esp = &app->species[einstein_idx];
    spacetime_step_backup = gkyl_array_new(GKYL_DOUBLE, esp->num_equations, app->local_ext.volume);
    gkyl_array_copy_range(spacetime_step_backup, esp->fcurr, &app->local_ext);
  }


  // [SPACETIME WB] Frozen static KS-TOV spacetime/fluid references. These are captured before
  // either member of the coupled pair advances, so the residual subtraction is anchored to the
  // actual t=0 equilibrium rather than a first-step perturbed fluid state.
  static int wb_spacetime = -1;
  static struct gkyl_array *wb_eq_st = NULL;
  static struct gkyl_array *wb_eq_fluid_st = NULL;
  static struct gkyl_array *wb_src_st = NULL;
  static struct gkyl_array *wb_trans_delta_st = NULL;
  static struct gkyl_array *wb_live_pre_st = NULL;
  if (wb_spacetime < 0) wb_spacetime = getenv("GKYL_WB_SPACETIME") ? 1 : 0;
  bool wb_spacetime_active = coupled && wb_spacetime
    && !app->species[einstein_idx].spacetime_equilibrium
    && !gkyl_gr_euler_wb_disabled(app->species[fluid_idx].equation);

  // Cap dt by the spacetime (wave_prop) stable step up front
  double st_max_dt = DBL_MAX;
  if (einstein_idx >= 0) {
    st_max_dt = moment_species_max_dt(app, &app->species[einstein_idx]);
  }
  double dt_use = fmin(dt0, st_max_dt);

  // 1.  Sync the live spacetime metric into the fluid state
  if (coupled) {
    mixed_sync_metric(app, fluid_idx, einstein_idx);
  }


  if (wb_spacetime_active && wb_eq_st == NULL) {
    struct moment_species *esp = &app->species[einstein_idx];
    struct moment_species *fsp = &app->species[fluid_idx];

    wb_eq_st = gkyl_array_new(GKYL_DOUBLE, esp->num_equations, app->local_ext.volume);
    wb_eq_fluid_st = gkyl_array_new(GKYL_DOUBLE, fsp->num_equations, app->local_ext.volume);
    wb_src_st = gkyl_array_new(GKYL_DOUBLE, esp->num_equations, app->local_ext.volume);
    wb_trans_delta_st = gkyl_array_new(GKYL_DOUBLE, esp->num_equations, app->local_ext.volume);
    wb_live_pre_st = gkyl_array_new(GKYL_DOUBLE, esp->num_equations, app->local_ext.volume);

    gkyl_array_copy_range(wb_eq_st, esp->f[0], &app->local_ext);
    gkyl_array_copy_range(wb_eq_fluid_st, fsp->fcurr, &app->local_ext);
    moment_species_apply_bc(app, tcurr, esp, wb_eq_st);
    compute_spacetime_wb_src(app, esp, wb_eq_st, wb_eq_fluid_st,
      fsp->gr_euler_gas_gamma, tcurr, wb_src_st);
  }

  route1_ledger_begin(app, fluid_idx);
  // 2. Advance the fluid (MP SSP-RK3)
  // Symmetric coupling is scoped to the two-sector conformal system. Other
  // species/field configurations retain their existing scheduling and backups.
  bool symmetric_coupling = coupled && app->num_species == 2 && !app->has_field
    && app->species[einstein_idx].has_vacuum_einstein_conformal;
  double fluid_fraction = symmetric_coupling ? 0.5 : 1.0;
  double dt = dt_use;
  if (fluid_idx >= 0) {
    struct gkyl_update_status fluid_st = route1_fluid_rk3(app, fluid_idx, tcurr, fluid_fraction*dt_use);
    route1_ledger_mark(app, fluid_idx, R1LG_FLUX);
    if (!fluid_st.success) {
      mixed_restore_step_backups(app, fluid_idx, einstein_idx,
        fluid_step_backup, spacetime_step_backup);
      mixed_release_step_backups(fluid_step_backup, spacetime_step_backup);
      fluid_st.dt_actual = 0.0;
      fluid_st.dt_suggested /= fluid_fraction;
      return fluid_st;
    }
    dt = fluid_st.dt_actual/fluid_fraction;
    st.dt_suggested = fmin(st.dt_suggested, fluid_st.dt_suggested/fluid_fraction);
    // [EXCISION TRANSPORT FREEZE] moment_update_ssp_rk3's own transport touches every cell,
    // excised or not (see mixed_restore_excised_fluid_cells). Undo that for cells that were
    // already excised before this step, using the pre-step backup.
    mixed_restore_excised_fluid_cells(app, &app->species[fluid_idx], fluid_step_backup);
    route1_ledger_mark(app, fluid_idx, R1LG_EXCISION);
    moment_species_apply_bc(app, tcurr + fluid_fraction*dt,
      &app->species[fluid_idx], app->species[fluid_idx].fcurr);
    if (coupled && !mixed_gr_euler_state_valid(app, &app->species[fluid_idx])) {
      st.success = false;
      st.dt_actual = 0.0;
      st.dt_suggested = 0.5 * dt;
      mixed_restore_step_backups(app, fluid_idx, einstein_idx,
        fluid_step_backup, spacetime_step_backup);
      mixed_release_step_backups(fluid_step_backup, spacetime_step_backup);
      return st;
    }
  }



  // 3. Advance the spacetime
  if (einstein_idx >= 0) {
    struct moment_species *esp = &app->species[einstein_idx];
    int ndim = esp->ndim;
    bool is_conformal = esp->has_vacuum_einstein_conformal;
    double excision_threshold = is_conformal ? esp->vacuum_einstein_conformal_excision_threshold
      : esp->vacuum_einstein_excision_threshold;

    // The symmetric path holds midpoint fluid conserved variables fixed; the
    // source stage reconstructs their primitives with its own live geometry.
    const struct gkyl_array *fluid_f = coupled ? app->species[fluid_idx].fcurr : NULL;
    double gas_gamma = coupled ? app->species[fluid_idx].gr_euler_gas_gamma : 0.0;

    // [SPACETIME WB] Well-balance the gauge/spacetime evolution against the frozen KS-TOV equilibrium: subtract
    // the discrete source S_eq the gauge produces on the analytic equilibrium so the star's spacetime is held
    // (the drift persists with the fluid frozen). S_eq is
    // computed once before the fluid step from the pristine t=0 coupled equilibrium. Env
    // GKYL_WB_SPACETIME; default off => NULL => bit-identical. Respects the collapse-stage WB-off switch
    // (gkyl_gr_euler_wb_disabled).
    const struct gkyl_array *wb_src_pass = NULL;
    const struct gkyl_array *wb_eq_pass = NULL;
    struct gkyl_array *wb_trans_delta_pass = NULL;
    struct gkyl_array *wb_live_pre_pass = NULL;
    if (wb_spacetime_active && wb_eq_st != NULL) {
      wb_src_pass = wb_src_st;
      wb_eq_pass = wb_eq_st;
      wb_trans_delta_pass = wb_trans_delta_st;
      wb_live_pre_pass = wb_live_pre_st;
    }
    if (esp->spacetime_equilibrium) {
      wb_src_pass = esp->spacetime_equilibrium_source;
      wb_eq_pass = esp->spacetime_equilibrium;
    }

    // [ROUTE1] method-of-lines SSP-RK3 spacetime step (route1_step.h); replaces the composition below.
    if (route1_enabled()) {
      if (!route1_step(app, esp, fluid_f, gas_gamma, tcurr, dt, wb_eq_pass, wb_src_pass)) {
        st.success = false;
        st.dt_actual = 0.0;
        st.dt_suggested = 0.5 * dt;
        mixed_restore_step_backups(app, fluid_idx, einstein_idx,
          fluid_step_backup, spacetime_step_backup);
        mixed_release_step_backups(fluid_step_backup, spacetime_step_backup);
        return st;
      }
      goto route1_done;
    }
    // Symmetric Gamma/spacetime composition inside the fluid half-steps.
    gamma_driver_shift_rk2(app,esp,tcurr,0.5*dt);

    // NewRad Snapshot the start-of-step spacetime state (f[0] here is fcurr, ghosts valid from the previous step's apply_bc) so the radiative BC can forward-Euler the boundary layer over the
    // whole step, replacing the hyperbolic result there after 3c.
    // SINGLE-BLOCK ONLY: it clamps to app->local, so it is disabled under domain decomposition (comm size > 1), where per-rank boundaries are not the physical boundary.
    bool use_newrad = mixed_newrad_enabled(app);
    struct gkyl_array *newrad_fold = 0;
    if (use_newrad) {
      newrad_fold = gkyl_array_new(GKYL_DOUBLE, esp->num_equations, app->local_ext.volume);
      gkyl_array_copy_range(newrad_fold, esp->f[0], &app->local_ext);
    }
    struct gkyl_array *cell_spd_anchor = 0;
    if (mixed_spacetime_cell_spd_repair_enabled()) {
      cell_spd_anchor = gkyl_array_new(GKYL_DOUBLE, esp->num_equations, app->local_ext.volume);
      gkyl_array_copy_range(cell_spd_anchor, esp->f[0], &app->local_ext);
    }

    // 3a. Pre-transport half source on f[0] (= fcurr); refresh ghosts for the transport reconstruction.
    moment_spacetime_source_advance(app, esp, esp->f[0], fluid_f, gas_gamma, tcurr, 0.5 * dt, wb_src_pass);
    moment_species_apply_bc(app, tcurr, esp, esp->f[0]);

    // Hyperbolicity guard before the transport solve (so a non-SPD pre-source metric didn't enter wave_prop)
    if (cell_spd_anchor &&
        !mixed_spacetime_cell_spd_repair(app, esp, cell_spd_anchor)) {
      if (newrad_fold) {
        gkyl_array_release(newrad_fold);
      }
      if (cell_spd_anchor) {
        gkyl_array_release(cell_spd_anchor);
      }
      st.success = false;
      st.dt_actual = 0.0;
      st.dt_suggested = 0.5 * dt;
      mixed_restore_step_backups(app, fluid_idx, einstein_idx,
        fluid_step_backup, spacetime_step_backup);
      mixed_release_step_backups(fluid_step_backup, spacetime_step_backup);
      return st;
    }
    if (cell_spd_anchor) {
      moment_species_apply_bc(app, tcurr, esp, esp->f[0]);
    }
    if (!mixed_spacetime_metric_valid(app, esp)) {
      if (newrad_fold) {
        gkyl_array_release(newrad_fold);
      }
      if (cell_spd_anchor) {
        gkyl_array_release(cell_spd_anchor);
      }
      st.success = false;
      st.dt_actual = 0.0; // rejected step: do not advance app->tcurr (moment.c adds dt_actual unconditionally)
      st.dt_suggested = 0.5 * dt;
      mixed_restore_step_backups(app, fluid_idx, einstein_idx,
        fluid_step_backup, spacetime_step_backup);
      mixed_release_step_backups(fluid_step_backup, spacetime_step_backup);
      return st;
    }

    // 3b. Transport: f[0] -> f[ndim] (pure wave_prop, no sources inside).
    if (wb_trans_delta_pass != NULL) {
      gkyl_array_copy_range(wb_live_pre_pass, esp->f[0], &app->local_ext);
      gkyl_array_copy_range(esp->f[0], wb_eq_pass, &app->local_ext);
      moment_species_apply_bc(app, tcurr, esp, esp->f[0]);

      struct gkyl_update_status wb_ts = symmetric_coupling
        ? moment_species_update_symmetric(app, esp, tcurr, dt)
        : moment_species_update(app, esp, tcurr, dt);
      if (!wb_ts.success) {
        gkyl_array_copy_range(esp->f[0], wb_live_pre_pass, &app->local_ext);
        if (newrad_fold) {
          gkyl_array_release(newrad_fold);
        }
        if (cell_spd_anchor) {
          gkyl_array_release(cell_spd_anchor);
        }
        mixed_restore_step_backups(app, fluid_idx, einstein_idx,
          fluid_step_backup, spacetime_step_backup);
        mixed_release_step_backups(fluid_step_backup, spacetime_step_backup);
        return wb_ts;
      }

      gkyl_array_clear(wb_trans_delta_pass, 0.0);
      struct gkyl_range_iter wbtit;
      gkyl_range_iter_init(&wbtit, &app->local);
      while (gkyl_range_iter_next(&wbtit)) {
        long loc = gkyl_range_idx(&app->local_ext, wbtit.idx);
        const double *qtr = gkyl_array_cfetch(esp->f[ndim], loc);
        const double *qeq = gkyl_array_cfetch(wb_eq_pass, loc);
        double *dq = gkyl_array_fetch(wb_trans_delta_pass, loc);
        for (int m = 0; m < esp->num_equations; ++m) {
          // The live transport preserves beta/B, so its WB reference must too.
          dq[m] = (m >= 52 && m < 64) ? 0.0 : qtr[m] - qeq[m];
        }
      }

      gkyl_array_copy_range(esp->f[0], wb_live_pre_pass, &app->local_ext);
    }

    struct gkyl_update_status ts = symmetric_coupling
      ? moment_species_update_symmetric(app, esp, tcurr, dt)
      : moment_species_update(app, esp, tcurr, dt);
    if (!ts.success) {
      if (newrad_fold) {
        gkyl_array_release(newrad_fold);
      }
      if (cell_spd_anchor) {
        gkyl_array_release(cell_spd_anchor);
      }
      mixed_restore_step_backups(app, fluid_idx, einstein_idx,
        fluid_step_backup, spacetime_step_backup);
      mixed_release_step_backups(fluid_step_backup, spacetime_step_backup);
      return ts;
    }
    st.dt_suggested = fmin(st.dt_suggested, ts.dt_suggested);

    // Bring the transported state back into f[0] for the post source.
    keep_shift_through_transport(app, esp);
    gkyl_array_copy_range(esp->f[0], esp->f[ndim], &app->local_ext);
    if (wb_trans_delta_pass != NULL) {
      subtract_spacetime_wb_delta(app, esp, wb_trans_delta_pass);
    }
    // [EXCISION TRANSPORT FREEZE] wave_prop's transport just touched every cell, excised or
    // not (see mixed_restore_excised_spacetime_cells). Undo that for cells already excised
    // before this step, using the pre-step backup -- before anything downstream reads f[0].
    mixed_restore_excised_spacetime_cells(app, esp, spacetime_step_backup,
      excision_threshold, is_conformal);

    // 3c. Post-transport half source on f[0]; reapply BC so ghosts are valid for the next step's sync.
    moment_spacetime_source_advance(app, esp, esp->f[0], fluid_f, gas_gamma, tcurr + 0.5 * dt, 0.5 * dt, wb_src_pass);
    moment_species_apply_bc(app, tcurr + dt, esp, esp->f[0]);

    // NewRad Radiative BC: overwrite the outermost interior layer with the outgoing-wave forward-Euler update from the start-of-step snapshot, then refresh ghosts. Replaces the hyperbolic result there.
    if (use_newrad) {
      apply_newrad_spacetime_bc(app, esp, newrad_fold, tcurr + dt, dt, wb_eq_pass);
      gkyl_array_release(newrad_fold);
      newrad_fold = 0;
    }


    // 3d. Hyperbolicity guard after the full spacetime step
    if (cell_spd_anchor &&
        !mixed_spacetime_cell_spd_repair(app, esp, cell_spd_anchor)) {
      if (cell_spd_anchor) {
        gkyl_array_release(cell_spd_anchor);
      }
      st.success = false;
      st.dt_actual = 0.0;
      st.dt_suggested = 0.5 * dt;
      mixed_restore_step_backups(app, fluid_idx, einstein_idx,
        fluid_step_backup, spacetime_step_backup);
      mixed_release_step_backups(fluid_step_backup, spacetime_step_backup);
      return st;
    }
    if (cell_spd_anchor) {
      moment_species_apply_bc(app, tcurr + dt, esp, esp->f[0]);
    }
    if (!mixed_spacetime_metric_valid(app, esp)) {
      if (cell_spd_anchor) {
        gkyl_array_release(cell_spd_anchor);
      }
      st.success = false;
      st.dt_actual = 0.0;
      st.dt_suggested = 0.5 * dt;
      mixed_restore_step_backups(app, fluid_idx, einstein_idx,
        fluid_step_backup, spacetime_step_backup);
      mixed_release_step_backups(fluid_step_backup, spacetime_step_backup);
      return st;
    }

    // Complete the symmetric Gamma-driver subflow.
    gamma_driver_shift_rk2(app,esp,tcurr+0.5*dt,0.5*dt);
    if (cell_spd_anchor) {
      gkyl_array_release(cell_spd_anchor);
    }
    route1_done:;
  }
  if (coupled) {
    mixed_sync_metric(app, fluid_idx, einstein_idx);
    if (!mixed_gr_euler_state_valid(app, &app->species[fluid_idx])) {
      st.success = false;
      st.dt_actual = 0.0;
      st.dt_suggested = 0.5 * dt;
      mixed_restore_step_backups(app, fluid_idx, einstein_idx,
        fluid_step_backup, spacetime_step_backup);
      mixed_release_step_backups(fluid_step_backup, spacetime_step_backup);
      return st;
    }
  }
  if (symmetric_coupling) {
    struct moment_species *fsp = &app->species[fluid_idx];
    // Include newly excised cells in the transport freeze for this half-step.
    struct gkyl_array *half_backup = gkyl_array_new(GKYL_DOUBLE,
      fsp->num_equations, app->local_ext.volume);
    gkyl_array_copy_range(half_backup, fsp->fcurr, &app->local_ext);
    moment_species_apply_bc(app, tcurr + 0.5*dt, fsp, fsp->fcurr);
    struct gkyl_update_status fluid_st = route1_fluid_rk3(app, fluid_idx, tcurr + 0.5*dt, 0.5*dt);
    route1_ledger_mark(app, fluid_idx, R1LG_FLUX);
    // A shortened last half-step cannot be paired with the completed spacetime
    // step. Reject the entire composition, including the first fluid half-step.
    if (!fluid_st.success || fluid_st.dt_actual != 0.5*dt) {
      st.success = false;
      st.dt_actual = 0.0;
      // A successful but shortened subflow measured a CFL limit, not a bad
      // state. Retry the full composition at that limit instead of halving it.
      st.dt_suggested = fluid_st.success && fluid_st.dt_actual > 0.0
        ? 2.0*fmin(fluid_st.dt_actual, fluid_st.dt_suggested)
        : fmin(0.5*dt, 2.0*fluid_st.dt_suggested);
      if (fluid_st.success && fluid_st.dt_actual > 0.0)
        gkyl_moment_app_cout(app, stdout,
          " [fluid CFL] final half-step shortened from %g to %g; retry full step at %g ...",
          0.5*dt, fluid_st.dt_actual, st.dt_suggested);
      gkyl_array_release(half_backup);
      mixed_restore_step_backups(app, fluid_idx, einstein_idx,
        fluid_step_backup, spacetime_step_backup);
      mixed_release_step_backups(fluid_step_backup, spacetime_step_backup);
      return st;
    }
    st.dt_suggested = fmin(st.dt_suggested, 2.0*fluid_st.dt_suggested);
    mixed_restore_excised_fluid_cells(app, fsp, half_backup);
    route1_ledger_mark(app, fluid_idx, R1LG_EXCISION);
    gkyl_array_release(half_backup);
    mixed_sync_metric(app, fluid_idx, einstein_idx);
    moment_species_apply_bc(app, tcurr + dt, fsp, fsp->fcurr);
    if (!mixed_gr_euler_state_valid(app, fsp)) {
      st.success = false;
      st.dt_actual = 0.0;
      st.dt_suggested = 0.5*dt;
      mixed_restore_step_backups(app, fluid_idx, einstein_idx,
        fluid_step_backup, spacetime_step_backup);
      mixed_release_step_backups(fluid_step_backup, spacetime_step_backup);
      return st;
    }
  }
  mixed_release_step_backups(fluid_step_backup, spacetime_step_backup);

  route1_ledger_end(app, fluid_idx, tcurr, dt);
  st.dt_actual = dt;
  if (st.dt_suggested == DBL_MAX) {
    st.dt_suggested = dt;
  }
  return st;
}
