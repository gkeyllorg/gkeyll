// ** Fully formally-verified solvers for the perfectly hyperbolic Maxwell equations **
// ** Lax-Friedrichs Solver: **
// ** Proof of hyperbolicity preservation (Ex and phi components): ../proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_hyperbolicity.rkt **
// ** Proof of hyperbolicity preservation (Ey and Bz components): ../proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_hyperbolicity.rkt **
// ** Proof of hyperbolicity preservation (Ez and By components): ../proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_hyperbolicity.rkt **
// ** Proof of hyperbolicity preservation (Bx and psi components): ../proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Ex and phi components): ../proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_strict_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Ey and Bz components): ../proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_strict_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Ez and By components): ../proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_strict_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Bx and psi components): ../proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_strict_hyperbolicity.rkt **
// ** Proof of CFL stability (Ex and phi components): ../proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_cfl_stability.rkt **
// ** Proof of CFL stability (Ey and Bz components): ../proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_cfl_stability.rkt **
// ** Proof of CFL stability (Ez and By components): ../proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_cfl_stability.rkt **
// ** Proof of CFL stability (Bx and psi components): ../proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_cfl_stability.rkt **
// ** Proof of local Lipschitz continuity of discrete flux function (Ex and phi components): ../proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_local_lipschitz.rkt **
// ** Proof of local Lipschitz continuity of discrete flux function (Ey and Bz components): ../proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_local_lipschitz.rkt **
// ** Proof of local Lipschitz continuity of discrete flux function (Ez and By components): ../proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_local_lipschitz.rkt **
// ** Proof of local Lipschitz continuity of discrete flux function (Bx and psi components): ../proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_local_lipschitz.rkt **
// ** Roe Solver: **
// ** Proof of hyperbolicity preservation (Ex and phi components): ../proofs/finite_volume/proof_maxwell_1d_Ex_phi_roe_hyperbolicity.rkt **
// ** Proof of hyperbolicity preservation (Ey and Bz components): ../proofs/finite_volume/proof_maxwell_1d_Ey_Bz_roe_hyperbolicity.rkt **
// ** Proof of hyperbolicity preservation (Ez and By components): ../proofs/finite_volume/proof_maxwell_1d_Ez_By_roe_hyperbolicity.rkt **
// ** Proof of hyperbolicity preservation (Bx and psi components): ../proofs/finite_volume/proof_maxwell_1d_Bx_psi_roe_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Ex and phi components): ../proofs/finite_volume/proof_maxwell_1d_Ex_phi_roe_strict_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Ey and Bz components): ../proofs/finite_volume/proof_maxwell_1d_Ey_Bz_roe_strict_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Ez and By components): ../proofs/finite_volume/proof_maxwell_1d_Ez_By_roe_strict_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Bx and psi components): ../proofs/finite_volume/proof_maxwell_1d_Bx_psi_roe_strict_hyperbolicity.rkt **
// ** Proof of flux conservation (jump continuity, Ex and phi components): ../proofs/finite_volume/proof_maxwell_1d_Ex_phi_roe_flux_conservation.rkt **
// ** Proof of flux conservation (jump continuity, Ey and Bz components): ../proofs/finite_volume/proof_maxwell_1d_Ey_Bz_roe_flux_conservation.rkt **
// ** Proof of flux conservation (jump continuity, Ez and By components): ../proofs/finite_volume/proof_maxwell_1d_Ez_By_roe_flux_conservation.rkt **
// ** Proof of flux conservation (jump continuity, Bx and psi components): ../proofs/finite_volume/proof_maxwell_1d_Bx_psi_roe_flux_conservation.rkt **

#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_maxwell.h>
#include <gkyl_wv_maxwell_priv.h>

static inline void
maxwell_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  for (int i = 0; i < 8; i++) {
    sout[i] = 0.0;
  }
}

void
gkyl_wv_maxwell_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_maxwell *maxwell = container_of(base->on_dev, struct wv_maxwell, eqn);
    gkyl_cu_free(maxwell);
  }

  struct wv_maxwell *maxwell = container_of(base, struct wv_maxwell, eqn);
  gkyl_free(maxwell);
}

struct gkyl_wv_eqn*
gkyl_wv_maxwell_new(double c, double e_fact, double b_fact, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_wv_maxwell_cu_dev_new(c, e_fact, b_fact);
  } 
#endif

  return gkyl_wv_maxwell_inew(&(struct gkyl_wv_maxwell_inp) {
      .c = c,
      .e_fact = e_fact,
      .b_fact = b_fact,
      .rp_type = WV_MAXWELL_RP_ROE,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_maxwell_inew(const struct gkyl_wv_maxwell_inp* inp)
{
  struct wv_maxwell *maxwell = gkyl_malloc(sizeof(struct wv_maxwell));

  maxwell->eqn.type = GKYL_EQN_MAXWELL;
  maxwell->eqn.num_equations = 8;
  maxwell->eqn.num_diag = 6;

  maxwell->c = inp->c;
  maxwell->e_fact = inp->e_fact;
  maxwell->b_fact = inp->b_fact;

  if (inp->rp_type == WV_MAXWELL_RP_ROE) {
    maxwell->eqn.num_waves = 6;
    maxwell->eqn.waves_func = wave;
    maxwell->eqn.qfluct_func = qfluct;
  }
  else if (inp->rp_type == WV_MAXWELL_RP_LAX) {
    maxwell->eqn.num_waves = 2;
    maxwell->eqn.waves_func = wave_lax_l;
    maxwell->eqn.qfluct_func = qfluct_lax_l;
  }

  maxwell->eqn.flux_jump = flux_jump;
  maxwell->eqn.check_inv_func = check_inv;
  maxwell->eqn.max_speed_func = max_speed;
  maxwell->eqn.rotate_to_local_func = rot_to_local;
  maxwell->eqn.rotate_to_global_func = rot_to_global;
  
  maxwell->eqn.wall_bc_func = maxwell_wall;
  maxwell->eqn.no_slip_bc_func = maxwell_no_slip;

  maxwell->eqn.cons_to_riem = cons_to_riem;
  maxwell->eqn.riem_to_cons = riem_to_cons;

  maxwell->eqn.cons_to_diag = maxwell_cons_to_diag;

  maxwell->eqn.source_func = maxwell_source;

  maxwell->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(maxwell->eqn.flags);
  maxwell->eqn.ref_count = gkyl_ref_count_init(gkyl_wv_maxwell_free);
  maxwell->eqn.on_dev = &maxwell->eqn; // On the CPU, the equation object points to itself.

  maxwell->eqn.embed_geo = inp->embed_geo;
  if (maxwell->eqn.embed_geo) {
    switch (maxwell->eqn.embed_geo->type) {
      case GKYL_EMBED_COPY_B:
        maxwell->eqn.embed_geo->embed_func = wave_embed_copy_B;
        break;

      case GKYL_EMBED_PEC:
        maxwell->eqn.embed_geo->embed_func = wave_embed_pec;
        break;

      case GKYL_EMBED_FUNC:
        break; // already set by gkyl_wv_embed_geo_new

      default:
        assert(false);
        break;
    }
  }

  return &maxwell->eqn;
}
