#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include <gkyl_alloc.h>

#include "tov_solver.h"

struct gkyl_tov {
    int N; // how many radial points

    double K; // polytripic EOS - normalizaiton constant
    double Gamma; // polytripic EOS - adiabatic index

    //radial profiles
    double *r_areal; // areal radius
    double *m; 
    double *P;
    double *rho; // rest mass density
    double *eps; // specific internal energy
    double *e; // total energy density
    double *Phi; // potential (gravitaional potential, metric potential: alpha = exp(Phi))

    // global star properties
    double M_star; // total mass of the star (gravitational mass)
    double R_areal; // areal radius of the star surface
    double rho_c; // density at the center
    double P_c; // pressure at the center
    

    double dr; // radial step

};

// EOS is a single polytrope, so: P = K rho^Gamma

static inline double eos_pressure(double K, double Gamma, double rho)
{
    double P;
    if (rho > 0) {
        P = K * pow(rho, Gamma);
    }
    else P = 0.0;
    return P;
}

// eps = K rho^(gamma-1) / (gamma - 1)
static inline double eos_eps(double K, double Gamma, double rho)
{
    double eps;
    if (rho > 0) {
        eps = K * pow(rho, Gamma - 1.0) / (Gamma - 1.0);
    }
    else eps = 0.0;
    return eps;
}

// e = rho * (1 + eps)
static inline double eos_energy(double K, double Gamma, double rho)
{
    double e;
    if (rho > 0) {
        e = rho + K * pow(rho, Gamma) / (Gamma - 1.0);
    }
    else e = 0.0;
    return e;
}

// rho = (P / K)^(1/Gamma)
static inline double eos_rho_from_pressure(double K, double Gamma, double P)
{
    double rho;
    if (P > 0) {
        rho = pow(P / K, 1.0 / Gamma);
    }
    else rho = 0.0;
    return rho;
}

// TOV integrator
static void tov_rhs(double r, const double y[3], double f[3], double K, double Gamma)
{
    double m_enc = y[0]; // mass enclosed
    double P = fmax(y[1], 0.0);

    if (r <= 0.0) {
        f[0] = f[1] = f[2] = 0.0;
        return;
    }

    double rho = eos_rho_from_pressure(K, Gamma, P);
    double e = eos_energy(K, Gamma, rho);

    if (r * (r - 2.0 * m_enc) <= 0.0) { // denominator from dP/dr
        f[0] = f[1] = f[2] = 0.0;  // if (r * (r - 2.0 * m_enc) <= 0.0) then we're at/inside the horizon - unphysical for a stable star
        return;
    }
    
    f[0] = 4.0 * M_PI * r * r * e; // dm/dr RHS
    f[1] = -(e + P) * (m_enc + 4.0 * M_PI * r * r * r * P) / (r * (r - 2.0 * m_enc)); // dP/dr RHS
    // An RK stage can cross the surface. Vacuum has dm/dr=dP/dr=0,
    // but Phi must still follow the exterior Schwarzschild derivative.
    f[2] = (m_enc + 4.0 * M_PI * r * r * r * P) / (r * (r - 2.0 * m_enc)); // dPhi/dr 

}

static void implicit_tov_source_update(double r, double y[3], double h, double K, double Gamma) // it's an explicit method, not implicit, it's ssp-rk3
{
    double y_old[3], y_stage1[3], y_stage2[3], y_new[3];
    double f[3];

    for (int i = 0; i < 3; i ++) {
        y_old[i] = y[i]; // y_old = y^n
    }

    // stage 1 — Euler step from y_old
    tov_rhs(r, y_old, f, K, Gamma);
    for (int i = 0; i < 3; i++) {
        y_stage1[i] = y_old[i] + h * f[i];
    }

    // stage 2 — RHS at stage1, then mix with y_old
    tov_rhs(r + h, y_stage1, f, K, Gamma);
    for (int i = 0; i < 3; i++) {
        //y_stage2[i] = 0.75*y_old[i] + 0.25*(y_old[i] + h * f[i]);
        y_stage2[i] = 0.75*y_old[i] + 0.25*(y_stage1[i] + h * f[i]);
    }

    // stage 3 — RHS at stage2, final mix
    tov_rhs(r + 0.5*h, y_stage2, f, K, Gamma);
    for (int i = 0; i < 3; i++) {
        y_new[i] = (1.0/3.0)*y_old[i] + (2.0/3.0)*(y_stage2[i] + h * f[i]);
    }

    for (int i = 0; i < 3; i ++) {
        y[i] = y_new[i];
    }

}

// Interpolate on the Boyler Lindquist 
static double interp(const struct gkyl_tov *tov, const double *val, double r)
{
    if (r <= tov->r_areal[0]) {
        return val[0];
    }
    if (r >= tov->r_areal[tov->N - 1]) {
        return val[tov->N - 1];
    }

    int lo = 0, hi = tov->N - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (tov->r_areal[mid] <= r) lo = mid;
        else hi = mid;
    }

    double dr = tov->r_areal[hi] - tov->r_areal[lo];
    double t = 0.0;
    if (dr > 0.0) {
        t = (r - tov->r_areal[lo]) / dr;
    }
    else t = 0.0;

    return (1.0 - t) * val[lo] + t * val[hi];
}

bool
gkyl_tov_eval_mass_derivs(const struct gkyl_tov *tov, double r,
    double *dm, double *d2m)
{
    struct tov_eval_bl bl;
    if (!gkyl_tov_eval_bl(tov, r, &bl) || !dm || !d2m) return false;
    double e_prime = 0.0;
    if (r > 0.0 && bl.rho > 0.0 && bl.P > 0.0) {
        double denom = r * (r - 2.0 * bl.m);
        if (denom <= 0.0) return false;
        double p_prime = -(bl.e + bl.P) *
            (bl.m + 4.0 * M_PI * r*r*r * bl.P) / denom;
        double dp_drho = tov->Gamma * tov->K * pow(bl.rho, tov->Gamma - 1.0);
        e_prime = p_prime * (1.0/dp_drho + 1.0/(tov->Gamma - 1.0));
    }
    *dm = 4.0 * M_PI * r*r * bl.e;
    *d2m = 8.0 * M_PI * r * bl.e + 4.0 * M_PI * r*r * e_prime;
    return true;
}

// eval TOV in BL coords at arbitrary radius
bool gkyl_tov_eval_bl(const struct gkyl_tov *tov, double r, struct tov_eval_bl *out)
{
    if (!tov || !out || r < 0.0) {
        return false;
    }

    out->rho = interp(tov, tov->rho, r);
    out->P = interp(tov, tov->P, r);
    out->eps = interp(tov, tov->eps, r);
    out->e = interp(tov, tov->e, r);
    out->m = interp(tov, tov->m, r);
    out->Phi = interp(tov, tov->Phi, r);

    return true;
}


// eval fluid relative to star center - r is the areal radius shared between Schwarzschild and Kerr-Schild
bool gkyl_tov_eval(const struct gkyl_tov *tov, double x, double y, double z,
    struct tov_eval_fluid *out)
{
    if (!tov || !out) return false;
    double r = sqrt(x*x + y*y + z*z);
    
    if (r >= tov->R_areal) {
        out->rho = 0.0;
        out->P   = 0.0;
        out->eps = 0.0;
        out->e   = 0.0;
        return true;
    }
    
    out->rho = interp(tov, tov->rho, r);
    out->P   = interp(tov, tov->P,   r);
    out->eps = interp(tov, tov->eps, r);
    out->e   = interp(tov, tov->e,   r);
    return true;
}

// helpers to acces star mass and radius
double gkyl_tov_star_mass(const struct gkyl_tov *tov)
{
    if(tov) {
        return tov->M_star;
    }
    else {
        return 0.0;
    }
}

double gkyl_tov_star_radius(const struct gkyl_tov *tov)
{
    if(tov) {
        return tov->R_areal;
    }
    else {
        return 0.0;
    }
}

//main
struct gkyl_tov * gkyl_tov_new(double K, double Gamma, double rho_c, double dr)
{
    struct gkyl_tov *tov = (struct gkyl_tov *) gkyl_malloc(sizeof(struct gkyl_tov));

    if (!tov) {
        return NULL;
    }

    tov->K = K;
    tov->Gamma = Gamma;
    tov->rho_c = rho_c;
    tov->dr = dr;

    tov->r_areal = gkyl_malloc(sizeof(double[TOV_MAX_POINTS]));
    tov->m = gkyl_malloc(sizeof(double[TOV_MAX_POINTS]));
    tov->P  = gkyl_malloc(sizeof(double[TOV_MAX_POINTS]));
    tov->rho = gkyl_malloc(sizeof(double[TOV_MAX_POINTS]));
    tov->eps = gkyl_malloc(sizeof(double[TOV_MAX_POINTS]));
    tov->e = gkyl_malloc(sizeof(double[TOV_MAX_POINTS]));
    tov->Phi = gkyl_malloc(sizeof(double[TOV_MAX_POINTS]));

    if (!tov->r_areal || !tov->m || !tov->P || !tov->rho || !tov->eps || !tov->e || !tov->Phi) {
        gkyl_tov_solution_release(tov);
        return NULL;
    }

    // First we integrate the TOV eqns in the areal coordinates
    double P_c = eos_pressure(K, Gamma, rho_c);
    double e_c = eos_energy(K, Gamma, rho_c);
    tov->P_c = P_c;

    // center
    tov->r_areal[0] = 0.0;
    tov->m[0] = 0.0;
    tov->P[0] = P_c;
    tov->rho[0] = rho_c;
    tov->eps[0] = eos_eps(K, Gamma, rho_c);
    tov->e[0] = e_c;
    tov->Phi[0] = 0.0;  // "provisional"; shifted below

    // Taylor expansion at r = dr  (because the TOV equations are singular-looking at r = 0)
    double r_n = dr;
    double y[3];
    y[0] = (4.0 / 3.0) * M_PI * e_c * r_n * r_n * r_n; //m 
    y[1] = P_c - 2.0 * M_PI * (e_c + P_c) * (e_c / 3.0 + P_c) * r_n * r_n; //P
    y[2] = 2.0 * M_PI * (e_c / 3.0 + P_c) * r_n * r_n; //Phi

    if (y[1] <= 0.0) {
            y[1] = 0.0;
        }

    double rho_n = eos_rho_from_pressure(K, Gamma, y[1]);
    tov->r_areal[1] = r_n;
    tov->m[1]       = y[0];
    tov->P[1]       = y[1];
    tov->rho[1]     = rho_n;
    tov->eps[1]     = eos_eps(K, Gamma, rho_n);
    tov->e[1]       = eos_energy(K, Gamma, rho_n);
    tov->Phi[1]     = y[2];
    int n = 2;     

    // Actual integration call - call to the SSP-RK3 (implicit_tov_source_update)
    double r_max = 1500.0;
    while (r_n < r_max && y[1] > 0.0 && n < TOV_MAX_POINTS) {
        implicit_tov_source_update(r_n, y, dr, K, Gamma);
        r_n += dr;

        if (y[1] <= 0.0) {
            y[1] = 0.0;
        }

        rho_n = eos_rho_from_pressure(K, Gamma, y[1]);

        tov->r_areal[n] = r_n;
        tov->m[n] = y[0];
        tov->P[n]  = y[1];
        tov->rho[n] = rho_n;
        tov->eps[n] = eos_eps(K, Gamma, rho_n);
        tov->e[n] = eos_energy(K, Gamma, rho_n);
        tov->Phi[n] = y[2];
        n++;

        // printf("r=%g m=%g P=%g Phi=%g\n", r_n, y[0], y[1], y[2]);
        if (y[1] <= 0.0) break;
    }

    tov->N = n;
    int i_surf = n - 1; // stellar surface
    tov->M_star  = tov->m[i_surf];
    tov->R_areal = tov->r_areal[i_surf];

    // We fix lapse: Phi(R) = 0.5 ln(1 - 2M/R). Lapse = e^Phi. At the stellar surface, the interior must match exterior Schwarzschild.
    // g_tt = -(1 - 2M/r) = -e^2Phi. => Phi(R) = 0.5 * ln(1 - 2M/R)
    double Phi_exact = 0.5 * log(1.0 - 2.0 * tov->M_star / tov->R_areal);
    double Phi_shift = Phi_exact - tov->Phi[i_surf]; // We integrated with Phi(0) = Phi_c = 0 before (arbitrary), so we need to shift the profile by a constant
    for (int i = 0; i < n; i++) {
        tov->Phi[i] += Phi_shift;
    }   
    // After the Phi shift -- sanity check that the convention is what we expect
    double alpha_surf_int = exp(tov->Phi[i_surf]);
    double alpha_surf_ext = sqrt(1.0 - 2.0*tov->M_star/tov->R_areal);

    assert(fabs(alpha_surf_int - alpha_surf_ext) < 1e-12 && "Lapse doesn't match on the surface"); // this should pass

    // Extend to vacuum exterior out to r_max
    while (r_n < r_max && n < TOV_MAX_POINTS) {
        r_n += dr;
        double V_ext = 2.0 * tov->M_star / r_n;
        tov->r_areal[n] = r_n;
        tov->m[n]       = tov->M_star;
        tov->P[n]       = 0.0;
        tov->rho[n]     = 0.0;
        tov->eps[n]     = 0.0;
        tov->e[n]       = 0.0;

        if (1.0 - V_ext > 0.0) {
            tov->Phi[n] = 0.5 * log(1.0 - V_ext);
        }
        else {
            tov->Phi[n] = tov->Phi[n-1];
        }
        n++;
    }

    tov->N = n;

    return tov;
}

void
gkyl_gr_tov_refresh_geometry_from_state(double gas_gamma, double p_atm, int ncells, double *q)
{
  double e_before = 0.0;
  double m_before = 0.0;
  double phi_before = 0.0;

  double r_before = 0.0;
  double dphi_dr_before = 0.0;
  int surface_idx = ncells - 1;

  for (int i = 0; i < ncells; ++i) {
    double *qi = &q[8*i];

    double r = qi[5];
    double r2 = r * r;
    if (r2 < 1.0e-300) {
      r2 = 1.0e-300;
    }

    double D = qi[0] / r2;
    double tau = qi[1] / r2;   // tau = Etot - D (Valencia energy)
    double mom_r = qi[2] / r2;

    double rho_floor = 1.0e-15;
    double p_floor = 1.0e-300;
    double D_safe = fmax(D, rho_floor);
    double p_min = fmax(p_floor, fabs(mom_r) - tau - D_safe + 1.0e-16);
    double p = fmax((gas_gamma - 1.0) * tau, p_min);
    double vel = 0.0;
    double rho = D_safe;

    for (int iter = 0; iter < 100; iter++) {
      double Q = tau + D_safe + p;
      vel = mom_r / Q;
      if (vel * vel > 1.0 - 1.0e-12) {
        vel = copysign(sqrt(1.0 - 1.0e-12), vel);
      }
      double W2 = 1.0 / (1.0 - vel * vel);
      double W = sqrt(W2);
      rho = D_safe / W;
      double eps = (tau - D_safe * (W - 1.0) - p * (W2 - 1.0)) / (D_safe * W);
      if (eps < 0.0) {
        eps = 0.0;
      }

      double p_eos = (gas_gamma - 1.0) * rho * eps;
      double f_val = p_eos - p;
      double h = 1.0 + eps + p / rho;
      double cs2 = gas_gamma * p / (rho * h);
      if (cs2 < 0.0) {
        cs2 = 0.0;
      }
      if (cs2 > 1.0 - 1.0e-12) {
        cs2 = 1.0 - 1.0e-12;
      }
      double f_prime = vel * vel * cs2 - 1.0;

      double dp = -f_val / f_prime;
      p += dp;
      if (p < p_min) {
        p = p_min;
      }
      if (fabs(dp) < 1.0e-12 * fmax(fabs(p), p_floor)) {
        break;
      }
    }

    double denom = tau + D_safe + p;
    vel = mom_r / denom;
    if (vel * vel > 1.0 - 1.0e-12) {
      vel = copysign(sqrt(1.0 - 1.0e-12), vel);
    }
    double W = 1.0 / sqrt(1.0 - vel * vel);
    rho = D_safe / W;

    bool in_atmosphere = false;
    if (p_atm > 0.0 && p <= (1.0 + 1.0e-12) * p_atm) {
      in_atmosphere = true;
    }

    double p_geom = 0.0;
    if (p_atm > 0.0) {
      p_geom = fmax(0.0, p - p_atm);
    }
    else {
      p_geom = p;
    }

    double rho_geom = 0.0;
    double e_geom = 0.0;
    if (p_geom > 0.0) {
      rho_geom = rho;
      e_geom = rho_geom + p_geom / (gas_gamma - 1.0);
    }

    if (in_atmosphere && surface_idx == ncells - 1) {
      if (i > 0) {
        surface_idx = i - 1;
      }
      else {
        surface_idx = 0;
      }
    }

    // double m = qi[4];
    // for refreshing m start:
    double m = 0.0;
    if (!in_atmosphere) {
      double dm_dr = 4.0 * M_PI * r * r * e_geom;
      if (i == 0) {
        m = (4.0 / 3.0) * M_PI * e_geom * r * r * r;
      }
      else {
        double dm_dr_before = 4.0 * M_PI * r_before * r_before * e_before;
        double dr = r - r_before;
        m = m_before + 0.5 * dr * (dm_dr_before + dm_dr);
      }
    }
    else {
      m = m_before;
    }
    qi[4] = m;
    // for refreshing m stop

    double dphi_dr = 0.0;
    if (!in_atmosphere && r > 0.0 && r - 2.0 * m > 0.0) {
      double radial_stress = mom_r * vel + p_geom; // fluid is moving, no longer just pressure
      dphi_dr = (m + 4.0 * M_PI * r * r * r * radial_stress) / (r * (r - 2.0 * m));
    }

    double phi = 0.0;
    if (i == 0) {
      phi = 0.0; // provisional; shift after sweep
    }
    else {
      double dr = r - r_before;
      phi = phi_before + 0.5*dr*(dphi_dr_before + dphi_dr);
    }

    qi[3] = phi;

    r_before = r;
    e_before = e_geom;
    m_before = m;
    phi_before = phi;
    dphi_dr_before = dphi_dr;
  }

  double *qsurf = &q[8*surface_idx];
  double R = qsurf[5];
  double M = qsurf[4];

  if (R > 0.0 && 1.0 - 2.0*M/R > 0.0) {
    double phi_match = 0.5 * log(1.0 - 2.0 * M / R);
    double phi_shift = phi_match - qsurf[3];

    for (int i = 0; i <= surface_idx; ++i) {
      q[8*i + 3] += phi_shift;
    }
  }

  for (int i = surface_idx + 1; i < ncells; ++i) {
    double *qi = &q[8*i];
    double r = qi[5];
    double m = qi[4];
    if (r > 0.0 && 1.0 - 2.0 * m / r > 0.0) {
      qi[3] = 0.5 * log(1.0 - 2.0 * m / r);
    }
  }
}

void gkyl_tov_solution_release(struct gkyl_tov *tov)
{
    if (!tov) return;
    gkyl_free(tov->r_areal);
    gkyl_free(tov->m);
    gkyl_free(tov->P);
    gkyl_free(tov->rho);
    gkyl_free(tov->eps);
    gkyl_free(tov->e);
    gkyl_free(tov->Phi);
    gkyl_free(tov);
}

