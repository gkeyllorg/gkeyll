#ifdef GKYL_HAVE_LUA

#include <gkyl_alloc.h>
#include <gkyl_eqn_type.h>
#include <gkyl_lua_utils.h>
#include <gkyl_lw_priv.h>
#include <gkyl_moment.h>
#include <gkyl_moment_lw.h>
#include <gkyl_moment_priv.h>
#include <gkyl_null_comm.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wv_coldfluid.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_euler_mixture.h>
#include <gkyl_wv_euler_rgfm.h>
#include <gkyl_wv_iso_euler.h>
#include <gkyl_wv_iso_euler_mixture.h>
#include <gkyl_wv_mhd.h>
#include <gkyl_wv_reactive_euler.h>
#include <gkyl_wv_sr_euler.h>
#include <gkyl_wv_ten_moment.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_blackhole.h>
#include <gkyl_gr_neutronstar.h>
#include <gkyl_gr_brill_lindquist.h>
#include <gkyl_wv_gr_maxwell.h>
#include <gkyl_wv_gr_maxwell_tetrad.h>
#include <gkyl_wv_gr_ultra_rel_euler.h>
#include <gkyl_wv_gr_ultra_rel_euler_tetrad.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_tetrad.h>
#include <gkyl_wv_gr_medium.h>
#include <gkyl_wv_vacuum_einstein.h>
#include <gkyl_wv_vacuum_einstein_conformal.h>
#include <gkyl_wv_gr_twofluid.h>
#include <gkyl_wv_gr_twofluid_tetrad.h>
#include <gkyl_wv_gr_mhd.h>
#include <gkyl_wv_gr_mhd_tetrad.h>
#include <gkyl_wv_advect.h>
#include <gkyl_wv_burgers.h>
#include <gkyl_zero_lw.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
#include <limits.h>

#include <limits.h>
#include <stdio.h>
#include <string.h>

#include <stc/coption.h>
#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

// Moment scheme type -> enum map.
static const struct gkyl_str_int_pair moment_scheme_type[] = {
  { "WaveProp", GKYL_MOMENT_WAVE_PROP },
  { "MP", GKYL_MOMENT_MP },
  { "KEP", GKYL_MOMENT_KEP },
  { 0, 0 }
};

// Wave limiter -> enum map.
static const struct gkyl_str_int_pair wave_limiter[] = {
  { "NoLimiter", GKYL_NO_LIMITER },
  { "MonotonizedCentered", GKYL_MONOTONIZED_CENTERED },
  { "MinMod", GKYL_MIN_MOD },
  { "SuperBee", GKYL_SUPERBEE },
  { "VanLeer", GKYL_VAN_LEER },
  { "BeamWarming", GKYL_BEAM_WARMING },
  { "Zero", GKYL_ZERO },
  { 0, 0 }
};

// Euler Riemann problem -> enum map.
static const struct gkyl_str_int_pair euler_rp_type[] = {
  { "Roe", WV_EULER_RP_ROE },
  { "HLLC", WV_EULER_RP_HLLC },
  { "Lax", WV_EULER_RP_LAX },
  { "HLL", WV_EULER_RP_HLL },
  { 0, 0 }
};

// MHD Riemann problem -> enum map.
static const struct gkyl_str_int_pair mhd_rp_type[] = {
  { "Roe", WV_MHD_RP_ROE },
  { "HLLD", WV_MHD_RP_HLLD },
  { "Lax", WV_MHD_RP_LAX },
  { 0, 0 }
};

// MHD divergence correction -> enum map.
static const struct gkyl_str_int_pair mhd_divb_type[] = {
  { "None", GKYL_MHD_DIVB_NONE },
  { "GLM", GKYL_MHD_DIVB_GLM },
  { "EightWaves", GKYL_MHD_DIVB_EIGHT_WAVES },
  { 0, 0 }
};

// Braginskii type -> enum map.
static const struct gkyl_str_int_pair braginskii_type[] = {
  { "Mag", GKYL_BRAG_MAG },
  { "Visc", GKYL_BRAG_VISC },
  { "HeatFlux", GKYL_BRAG_HEATFLUX },
  { "UnmagFull", GKYL_BRAG_UNMAG_FULL },
  { "MagFull", GKYL_BRAG_MAG_FULL },
  { 0, 0 }
};

// Spacetime gauge type -> enum map.
static const struct gkyl_str_int_pair spacetime_gauge_type[] = {
  { "Static", GKYL_STATIC_GAUGE },
  { "BlackHoleCollapse", GKYL_BLACKHOLE_COLLAPSE_GAUGE },
  { 0, 0 }
};

// Spacetime slicing type -> enum map.
static const struct gkyl_str_int_pair spacetime_slicing_type[] = {
  { "Geodesic", GKYL_GEODESIC_SLICING },
  { "Harmonic", GKYL_HARMONIC_SLICING },
  { "OnePlusLog", GKYL_1PLUSLOG_SLICING },
  { 0, 0 }
};

// Spacetime evolution type -> enum map.
static const struct gkyl_str_int_pair spacetime_evolution_type[] = {
  { "Ricci", GKYL_RICCI_EVOLUTION },
  { "Einstein", GKYL_EINSTEIN_EVOLUTION },
  { 0, 0 }
};

// Embedded boundary type -> enum map.
static const struct gkyl_str_int_pair embed_geo_type[] = {
  { "Absorb", GKYL_EMBED_ABSORB },
  { "Reflect", GKYL_EMBED_REFLECT },
  { "CopyB", GKYL_EMBED_COPY_B },
  { "PEC", GKYL_EMBED_PEC },
  { "Func", GKYL_EMBED_FUNC },
  { 0, 0 }
};

void
gkyl_register_moment_scheme_types(lua_State *L)
{
  register_types(L, moment_scheme_type, "SchemeType");
}

void
gkyl_register_wave_limiter_types(lua_State *L)
{
  register_types(L, wave_limiter, "WaveLimiter");
}

void
gkyl_register_euler_rp_types(lua_State *L)
{
  register_types(L, euler_rp_type, "EulerRP");
}

void
gkyl_register_mhd_rp_types(lua_State *L)
{
  register_types(L, mhd_rp_type, "MHDRP");
}

void
gkyl_register_mhd_divb_types(lua_State *L)
{
  register_types(L, mhd_divb_type, "DivB");
}

void
gkyl_register_braginskii_types(lua_State *L)
{
  register_types(L, braginskii_type, "Braginskii");
}

void
gkyl_register_spacetime_gauge_types(lua_State *L)
{
  register_types(L, spacetime_gauge_type, "SpacetimeGauge");
}

void
gkyl_register_spacetime_slicing_types(lua_State *L)
{
  register_types(L, spacetime_slicing_type, "SpacetimeSlicing");
}

void
gkyl_register_spacetime_evolution_types(lua_State *L)
{
  register_types(L, spacetime_evolution_type, "SpacetimeEvolution");
}

void
gkyl_register_embed_geo_types(lua_State *L)
{
  register_types(L, embed_geo_type, "EmbedBC");
}

// Magic IDs for use in distinguishing various species and field types.
enum moment_magic_ids {
  MOMENT_SPECIES_DEFAULT = 100, // Fluid species.
  MOMENT_FIELD_DEFAULT, // Maxwell equations.
  MOMENT_EQN_DEFAULT, // Equation object.
  MOMENT_SPACETIME_DEFAULT, // Spacetime object.
};

// Edge-splitting -> enum map.
static const struct gkyl_str_int_pair wave_split_type[] = {
  { "qwave", GKYL_WAVE_QWAVE },
  { "fwave", GKYL_WAVE_FWAVE },
  { 0, 0 }
};

// Reactive Euler Riemann problem -> enum map.
static const struct gkyl_str_int_pair reactive_euler_rp_type[] = {
  { "roe", WV_REACTIVE_EULER_RP_ROE },
  { "lax", WV_REACTIVE_EULER_RP_LAX },
  { 0, 0 }
};

// Euler mixture Riemann problem -> enum map.
static const struct gkyl_str_int_pair euler_mixture_rp_type[] = {
  { "roe", WV_EULER_MIXTURE_RP_ROE },
  { "lax", WV_EULER_MIXTURE_RP_LAX },
  { 0, 0 }
};

// Euler RGFM Riemann problem -> enum map.
static const struct gkyl_str_int_pair euler_rgfm_rp_type[] = {
  { "lax", WV_EULER_RGFM_RP_LAX },
  { 0, 0 }
};

// Isothermal Euler mixture Riemann problem -> enum map.
static const struct gkyl_str_int_pair iso_euler_mixture_rp_type[] = {
  { "roe", WV_ISO_EULER_MIXTURE_RP_ROE },
  { "lax", WV_ISO_EULER_MIXTURE_RP_LAX },
  { 0, 0 }
};

// General relativistic Maxwell Riemann problem -> enum map.
static const struct gkyl_str_int_pair gr_maxwell_rp_type[] = {
  { "roe", WV_GR_MAXWELL_RP_ROE },
  { "lax", WV_GR_MAXWELL_RP_LAX },
  { 0, 0 }
};

// General relativistic Maxwell Riemann problem in the tetrad basis -> enum map.
static const struct gkyl_str_int_pair gr_maxwell_tetrad_rp_type[] = {
  { "roe", WV_GR_MAXWELL_TETRAD_RP_ROE },
  { "lax", WV_GR_MAXWELL_TETRAD_RP_LAX },
  { 0, 0 }
};

// General relativistic Euler Riemann problem (ultra-relativistic equation of state) -> enum map.
static const struct gkyl_str_int_pair gr_ultra_rel_euler_rp_type[] = {
  { "hll", WV_GR_ULTRA_REL_EULER_RP_HLL },
  { "roe", WV_GR_ULTRA_REL_EULER_RP_ROE },
  { "lax", WV_GR_ULTRA_REL_EULER_RP_LAX },
  { 0, 0 }
};

// General relativistic Euler Riemann problem in the tetrad basis (ultra-relativistic equation of state) -> enum map.
static const struct gkyl_str_int_pair gr_ultra_rel_euler_tetrad_rp_type[] = {
  { "hll", WV_GR_ULTRA_REL_EULER_TETRAD_RP_HLL },
  { "roe", WV_GR_ULTRA_REL_EULER_TETRAD_RP_ROE },
  { "lax", WV_GR_ULTRA_REL_EULER_TETRAD_RP_LAX },
  { 0, 0 }
};

// General relativistic Euler Riemann problem (general equation of state) -> enum map.
static const struct gkyl_str_int_pair gr_euler_rp_type[] = {
  { "hll", WV_GR_EULER_RP_HLL },
  { "roe", WV_GR_EULER_RP_ROE },
  { "lax", WV_GR_EULER_RP_LAX },
  { 0, 0 }
};

// General relativistic Euler Riemann problem in the tetrad basis (general equation of state) -> enum map.
static const struct gkyl_str_int_pair gr_euler_tetrad_rp_type[] = {
  { "hll", WV_GR_EULER_TETRAD_RP_HLL },
  { "roe", WV_GR_EULER_TETRAD_RP_ROE },
  { "lax", WV_GR_EULER_TETRAD_RP_LAX },
  { 0, 0 }
};

// Coupled fluid-Einstein Riemann problem (plane-polarized Gowdy spacetimes) -> enum map.
static const struct gkyl_str_int_pair gr_medium_rp_type[] = {
  { "lax", WV_GR_MEDIUM_RP_LAX },
  { 0, 0 }
};

// Vacuum Einstein Riemann problem (Bona-Masso formalism) -> enum map.
static const struct gkyl_str_int_pair vacuum_einstein_rp_type[] = {
  { "hll", WV_VACUUM_EINSTEIN_RP_HLL },
  { "lax", WV_VACUUM_EINSTEIN_RP_LAX },
  { 0, 0 }
};

// Vacuum Einstein Riemann problem (Conformal Bona-Masso formalism) -> enum map.
static const struct gkyl_str_int_pair vacuum_einstein_conformal_rp_type[] = {
  { "hll", WV_VACUUM_EINSTEIN_CONFORMAL_RP_HLL },
  { "lax", WV_VACUUM_EINSTEIN_CONFORMAL_RP_LAX },
  { 0, 0 }
};

// General relativistic two-fluid Riemann problem -> enum map.
static const struct gkyl_str_int_pair gr_twofluid_rp_type[] = {
  { "hll", WV_GR_TWOFLUID_RP_HLL },
  { "lax", WV_GR_TWOFLUID_RP_LAX },
  { 0, 0 }
};

// General relativistic two-fluid Riemann problem in the tetrad basis -> enum map.
static const struct gkyl_str_int_pair gr_twofluid_tetrad_rp_type[] = {
  { "hll", WV_GR_TWOFLUID_TETRAD_RP_HLL },
  { "lax", WV_GR_TWOFLUID_TETRAD_RP_LAX },
  { 0, 0 }
};

// General relativistic magnetohydrodynamics Riemann problem -> enum map.
static const struct gkyl_str_int_pair gr_mhd_rp_type[] = {
  { "hll", WV_GR_MHD_RP_HLL },
  { "lax", WV_GR_MHD_RP_LAX },
  { 0, 0 }
};

// General relativistic magnetohydrodynamics Riemann problem in the tetrad basis -> enum map.
static const struct gkyl_str_int_pair gr_mhd_tetrad_rp_type[] = {
  { "hll", WV_GR_MHD_TETRAD_RP_HLL },
  { "lax", WV_GR_MHD_TETRAD_RP_LAX },
  { 0, 0 }
};

// Metatable name for equation object input struct.
#define MOMENT_WAVE_EQN_METATABLE_NM "GkeyllZero.App.Moments.Eq"

// Methods for manipulating gkyl_wv_eqn objects.

// Lua userdata object for constructing wave equation objects.
struct wv_eqn_lw {
  int magic; // This must be the first element in the struct.
  struct gkyl_wv_eqn *eqn; // Equation object.
  bool has_nn; // Does this equation object have a neural network attached?
  struct gkyl_kann_net *ann; // Neural network.
  bool has_spacetime; // Does this equation object have a spacetime object attached?
  struct gkyl_gr_spacetime *spacetime; // Spacetime.
};

// Clean up memory allocated for equation object.
static int
wv_eqn_lw_gc(lua_State *L)
{
  struct wv_eqn_lw **l_wv_lw = GKYL_CHECK_UDATA(L, MOMENT_WAVE_EQN_METATABLE_NM);
  struct wv_eqn_lw *wv_lw = *l_wv_lw;

  gkyl_wv_eqn_release(wv_lw->eqn);
  if (wv_lw->has_nn) {
    gkyl_kann_net_release(wv_lw->ann);
  }
  if (wv_lw->has_spacetime) {
    gkyl_gr_spacetime_release(wv_lw->spacetime);
  }
  gkyl_free(*l_wv_lw);
  
  return 0;
}

// Acquire equation object.
static struct gkyl_wv_eqn*
wv_eqn_get(lua_State *L)
{
  struct wv_eqn_lw **l_wv_lw = luaL_checkudata(L, -1, MOMENT_WAVE_EQN_METATABLE_NM);
  struct wv_eqn_lw *wv_lw = *l_wv_lw;

  return wv_lw->eqn;
}

/* *************** */
/* Euler Equations */
/* *************** */

// Euler.new { gasGamma = 1.4, rpType = G0.EulerRP.Roe }
// where rpType is one of G0.EulerRP.Roe, G0.EulerRP.Lax, G0.EulerRP.HLL or G0.EulerRP.HLLC.
static int
eqn_euler_lw_new(lua_State *L)
{
  struct wv_eqn_lw *euler_lw = gkyl_malloc(sizeof(*euler_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 1.4);
  enum gkyl_wv_euler_rp rp_type = glua_tbl_get_integer(L, "rpType", WV_EULER_RP_ROE);

  euler_lw->magic = MOMENT_EQN_DEFAULT;
  euler_lw->eqn = gkyl_wv_euler_inew( &(struct gkyl_wv_euler_inp) {
      .gas_gamma = gas_gamma,
      .rp_type = rp_type,
      .use_gpu = false
    }
  );
  euler_lw->has_nn = false;
  euler_lw->ann = 0;
  euler_lw->has_spacetime = false;
  euler_lw->spacetime = 0;

  // Create Lua userdata.
  struct wv_eqn_lw **l_euler_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_euler_lw = euler_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_euler_ctor[] = {
  { "new", eqn_euler_lw_new },
  { 0, 0 }
};

/* ************************** */
/* Isothermal Euler Equations */
/* ************************** */

// IsoEuler.new { vThermal = 1.0 }
static int
eqn_iso_euler_lw_new(lua_State *L)
{
  struct wv_eqn_lw *iso_euler_lw = gkyl_malloc(sizeof(*iso_euler_lw));

  double vt = glua_tbl_get_number(L, "vThermal", -1.0);
  if (vt < 0) {
    return luaL_error(L, "Thermal velocity \"vThermal\" not specified properly!");
  }
  
  iso_euler_lw->magic = MOMENT_EQN_DEFAULT;
  iso_euler_lw->eqn = gkyl_wv_iso_euler_new(vt, false);
  iso_euler_lw->has_nn = false;
  iso_euler_lw->ann = 0;
  iso_euler_lw->has_spacetime = false;
  iso_euler_lw->spacetime = 0;

  // Create Lua userdata.
  struct wv_eqn_lw **l_iso_euler_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_iso_euler_lw = iso_euler_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_iso_euler_ctor[] = {
  { "new", eqn_iso_euler_lw_new },
  { 0, 0 }
};

/* ************************************ */
/* Special Relativistic Euler Equations */
/* ************************************ */

// SrEuler.new { gasgamma = 1.4 }
static int
eqn_sr_euler_lw_new(lua_State *L)
{
  struct wv_eqn_lw *sr_euler_lw = gkyl_malloc(sizeof(*sr_euler_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 5.0 / 3.0);

  sr_euler_lw->magic = MOMENT_EQN_DEFAULT;
  sr_euler_lw->eqn = gkyl_wv_sr_euler_new(gas_gamma);
  sr_euler_lw->has_nn = false;
  sr_euler_lw->ann = 0;
  sr_euler_lw->has_spacetime = false;
  sr_euler_lw->spacetime = 0;

  // Create Lua userdata.
  struct wv_eqn_lw **l_sr_euler_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_sr_euler_lw = sr_euler_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_sr_euler_ctor[] = {
  { "new", eqn_sr_euler_lw_new },
  { 0, 0 }
};

/* ******************** */
/* Cold Fluid Equations */
/* ******************** */

// ColdFluid.new {  }
static int
eqn_coldfluid_lw_new(lua_State *L)
{
  struct wv_eqn_lw *coldfluid_lw = gkyl_malloc(sizeof(*coldfluid_lw));

  coldfluid_lw->magic = MOMENT_EQN_DEFAULT;
  coldfluid_lw->eqn = gkyl_wv_coldfluid_new();
  coldfluid_lw->has_nn = false;
  coldfluid_lw->ann = 0;
  coldfluid_lw->has_spacetime = false;
  coldfluid_lw->spacetime = 0;

  // Create Lua userdata.
  struct wv_eqn_lw **l_coldfluid_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_coldfluid_lw = coldfluid_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_coldfluid_ctor[] = {
  { "new", eqn_coldfluid_lw_new },
  { 0, 0 }
};

/* ******************** */
/* Ten-moment Equations */
/* ******************** */

// TenMoment.new { k0 = 1, hasGradClosure = false }
static int
eqn_tenmoment_lw_new(lua_State *L)
{
  struct wv_eqn_lw *tenm_lw = gkyl_malloc(sizeof(*tenm_lw));

  double k0 = glua_tbl_get_number(L, "k0", 1.0);
  bool has_grad_closure = glua_tbl_get_bool(L, "hasGradClosure", false);
  bool has_nn_closure = glua_tbl_get_bool(L, "hasNNClosure", false);
  int poly_order = glua_tbl_get_integer(L, "polyOrder", 1);
  const char* nn_species_name = glua_tbl_get_string(L, "NNSpeciesName", "");
  const char* nn_closure_file = glua_tbl_get_string(L, "NNClosureFile", "");
  
  struct gkyl_kann_net *ann = 0;
  if (has_nn_closure) {
    const char *fmt = "%s-%s.dat";
    int sz = gkyl_calc_strlen(fmt, nn_closure_file, nn_species_name);
    char fileNm[sz + 1];
    snprintf(fileNm, sizeof fileNm, fmt, nn_closure_file, nn_species_name);
    FILE *file = fopen(fileNm, "r");
    if (file != NULL) {
      ann = gkyl_kann_net_load(fileNm, false);
      fclose(file);
    }
    else {
      ann = 0;
      has_nn_closure = false;
      fprintf(stderr, "Neural network for %s species not found! Disabling NN-based closure.\n", nn_species_name);
    }
  }

  tenm_lw->magic = MOMENT_EQN_DEFAULT;
  tenm_lw->eqn = gkyl_wv_ten_moment_new(k0, has_grad_closure, has_nn_closure, poly_order, ann, false);
  tenm_lw->has_nn = has_nn_closure;
  tenm_lw->ann = ann;
  tenm_lw->has_spacetime = false;
  tenm_lw->spacetime = 0;

  // Create Lua userdata.
  struct wv_eqn_lw **l_tenm_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_tenm_lw = tenm_lw; // Point userdata to the equation object.
  
  // Set metatable
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_tenmoment_ctor[] = {
  { "new", eqn_tenmoment_lw_new },
  { 0, 0 }
};

/* ************* */
/* MHD Equations */
/* ************* */

// Mhd.new { gasgamma = 1.4, rpType = "roe", divB = "glm", glmCh = 0.0, glmAlpha = 0.0 }
// rpType is one of "roe", "hlld", "lax"
// divB is "none", "glm", "eight_waves"
static int
eqn_mhd_lw_new(lua_State *L)
{
  struct wv_eqn_lw *mhd_lw = gkyl_malloc(sizeof(*mhd_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 1.4);
  
  enum gkyl_wv_mhd_rp rp_type = glua_tbl_get_integer(L, "rpType", WV_MHD_RP_ROE);
  enum gkyl_wv_mhd_div_constraint divb = glua_tbl_get_integer(L, "divergenceConstraint", GKYL_MHD_DIVB_NONE);

  double glm_ch = glua_tbl_get_number(L, "glmCh", 1.0);
  double glm_alpha = glua_tbl_get_number(L, "glmAlpha", 0.4);

  mhd_lw->magic = MOMENT_EQN_DEFAULT;
  mhd_lw->eqn = gkyl_wv_mhd_new( &(struct gkyl_wv_mhd_inp) {
      .gas_gamma = gas_gamma,
      .rp_type = rp_type,
      .divergence_constraint = divb,
      .glm_alpha = glm_alpha,
      .glm_ch = glm_ch
    }
  );
  mhd_lw->has_nn = false;
  mhd_lw->ann = 0;
  mhd_lw->has_spacetime = false;
  mhd_lw->spacetime = 0;

  // Create Lua userdata.
  struct wv_eqn_lw **l_mhd_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_mhd_lw = mhd_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_mhd_ctor[] = {
  { "new", eqn_mhd_lw_new },
  { 0, 0 }
};

/* ************************ */
/* Reactive Euler Equations */
/* ************************ */

// ReactiveEuler.new { gasGamma = 1.4, specificHeatCapacity = 2.5, energyOfFormation = 1.0, ignitionTemperature = 0.25, reactionRate = 250.0, rpType = "roe" }
// where rpType is one of "roe" or "lax".
static int
eqn_reactive_euler_lw_new(lua_State *L)
{
  struct wv_eqn_lw *reactive_euler_lw = gkyl_malloc(sizeof(*reactive_euler_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 1.4);
  double specific_heat_capacity = glua_tbl_get_number(L, "specificHeatCapacity", 2.5);
  double energy_of_formation = glua_tbl_get_number(L, "energyOfFormation", 1.0);
  double ignition_temperature = glua_tbl_get_number(L, "ignitionTemperature", 0.25);
  double reaction_rate = glua_tbl_get_number(L, "reactionRate", 250.0);

  const char *rp_str = glua_tbl_get_string(L, "rpType", "lax");
  enum gkyl_wv_reactive_euler_rp rp_type = gkyl_search_str_int_pair_by_str(reactive_euler_rp_type, rp_str, WV_REACTIVE_EULER_RP_LAX);

  reactive_euler_lw->magic = MOMENT_EQN_DEFAULT;
  reactive_euler_lw->eqn = gkyl_wv_reactive_euler_inew( &(struct gkyl_wv_reactive_euler_inp) {
      .gas_gamma = gas_gamma,
      .specific_heat_capacity = specific_heat_capacity,
      .energy_of_formation = energy_of_formation,
      .ignition_temperature = ignition_temperature,
      .reaction_rate = reaction_rate,
      .rp_type = rp_type,
      .use_gpu = false
    }
  );
  reactive_euler_lw->has_nn = false;
  reactive_euler_lw->ann = 0;
  reactive_euler_lw->has_spacetime = false;
  reactive_euler_lw->spacetime = 0;

  // Create Lua userdata.
  struct wv_eqn_lw **l_reactive_euler_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_reactive_euler_lw = reactive_euler_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_reactive_euler_ctor[] = {
  { "new", eqn_reactive_euler_lw_new },
  { 0, 0 }
};

/* *********************** */
/* Euler Mixture Equations */
/* *********************** */

// EulerMixture.new { numComponents = 2, gasGamma = {1.4, 1.4}, rpType = "roe" }
// where rpType is one of "roe" or "lax".
static int
eqn_euler_mixture_lw_new(lua_State *L)
{
  struct wv_eqn_lw *euler_mixture_lw = gkyl_malloc(sizeof(*euler_mixture_lw));

  int num_components = glua_tbl_get_integer(L, "numComponents", 2);

  double *gas_gamma_s = gkyl_malloc(sizeof(double[num_components]));
  with_lua_tbl_tbl(L, "gasGamma") {
    for (int i = 0; i < num_components; i++) {
      gas_gamma_s[i] = glua_tbl_iget_number(L, i + 1, 1.4);
    }
  }
  
  const char *rp_str = glua_tbl_get_string(L, "rpType", "lax");
  enum gkyl_wv_euler_mixture_rp rp_type = gkyl_search_str_int_pair_by_str(euler_mixture_rp_type, rp_str, WV_EULER_MIXTURE_RP_LAX);

  euler_mixture_lw->magic = MOMENT_EQN_DEFAULT;
  euler_mixture_lw->eqn = gkyl_wv_euler_mixture_inew( & (struct gkyl_wv_euler_mixture_inp) {
      .num_species = num_components,
      .gas_gamma_s = gas_gamma_s,
      .rp_type = rp_type,
      .use_gpu = false
    }
  );
  euler_mixture_lw->has_nn = false;
  euler_mixture_lw->ann = 0;
  euler_mixture_lw->has_spacetime = false;
  euler_mixture_lw->spacetime = 0;

  // Create Lua userdata.
  struct wv_eqn_lw **l_euler_mixture_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_euler_mixture_lw = euler_mixture_lw; // Point userdata to the equation object.

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static const luaL_Reg eqn_euler_mixture_ctor[] = {
  { "new", eqn_euler_mixture_lw_new },
  { 0, 0 }
};

/* *********************** */
/* Euler RGFM Equations */
/* *********************** */

// EulerRGFM.new { numComponents = 2, gasGamma = {1.4, 1.4}, reinitFreq = 3, rpType = "lax" }
// where rpType is "lax".
static int
eqn_euler_rgfm_lw_new(lua_State *L)
{
  struct wv_eqn_lw *euler_rgfm_lw = gkyl_malloc(sizeof(*euler_rgfm_lw));

  int num_components = glua_tbl_get_integer(L, "numComponents", 2);

  double *gas_gamma_s = gkyl_malloc(sizeof(double[num_components]));
  with_lua_tbl_tbl(L, "gasGamma") {
    for (int i = 0; i < num_components; i++) {
      gas_gamma_s[i] = glua_tbl_iget_number(L, i + 1, 1.4);
    }
  }
  int reinit_freq = glua_tbl_get_integer(L, "reinitFreq", 3);
  
  const char *rp_str = glua_tbl_get_string(L, "rpType", "lax");
  enum gkyl_wv_euler_rgfm_rp rp_type = gkyl_search_str_int_pair_by_str(euler_rgfm_rp_type, rp_str, WV_EULER_RGFM_RP_LAX);

  euler_rgfm_lw->magic = MOMENT_EQN_DEFAULT;
  euler_rgfm_lw->eqn = gkyl_wv_euler_rgfm_inew( & (struct gkyl_wv_euler_rgfm_inp) {
      .num_species = num_components,
      .gas_gamma_s = gas_gamma_s,
      .reinit_freq = reinit_freq,
      .rp_type = rp_type,
      .use_gpu = false
    }
  );
  euler_rgfm_lw->has_nn = false;
  euler_rgfm_lw->ann = 0;
  euler_rgfm_lw->has_spacetime = false;
  euler_rgfm_lw->spacetime = 0;

  // Create Lua userdata.
  struct wv_eqn_lw **l_euler_rgfm_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_euler_rgfm_lw = euler_rgfm_lw; // Point userdata to the equation object.

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static const luaL_Reg eqn_euler_rgfm_ctor[] = {
  { "new", eqn_euler_rgfm_lw_new },
  { 0, 0 }
};

/* ********************************** */
/* Isothermal Euler Mixture Equations */
/* ********************************** */

// IsoEulerMixture.new { numComponents = 2, vThermal = {1.0, 1.0}, rpType = "roe" }
// where rpType is one of "roe" or "lax".
static int
eqn_iso_euler_mixture_lw_new(lua_State *L)
{
  struct wv_eqn_lw *iso_euler_mixture_lw = gkyl_malloc(sizeof(*iso_euler_mixture_lw));

  int num_components = glua_tbl_get_integer(L, "numComponents", 2);

  double *vt_s = gkyl_malloc(sizeof(double[num_components]));
  with_lua_tbl_tbl(L, "vThermal") {
    for (int i = 0; i < num_components; i++) {
      vt_s[i] = glua_tbl_iget_number(L, i + 1, 1.0);
    }
  }
  
  const char *rp_str = glua_tbl_get_string(L, "rpType", "lax");
  enum gkyl_wv_iso_euler_mixture_rp rp_type = gkyl_search_str_int_pair_by_str(iso_euler_mixture_rp_type, rp_str, WV_ISO_EULER_MIXTURE_RP_LAX);

  iso_euler_mixture_lw->magic = MOMENT_EQN_DEFAULT;
  iso_euler_mixture_lw->eqn = gkyl_wv_iso_euler_mixture_inew( & (struct gkyl_wv_iso_euler_mixture_inp) {
      .num_species = num_components,
      .vt_s = vt_s,
      .rp_type = rp_type,
      .use_gpu = false
    }
  );
  iso_euler_mixture_lw->has_nn = false;
  iso_euler_mixture_lw->ann = 0;
  iso_euler_mixture_lw->has_spacetime = false;
  iso_euler_mixture_lw->spacetime = 0;

  // Create Lua userdata.
  struct wv_eqn_lw **l_iso_euler_mixture_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_iso_euler_mixture_lw = iso_euler_mixture_lw; // Point userdata to the equation object.

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static const luaL_Reg eqn_iso_euler_mixture_ctor[] = {
  { "new", eqn_iso_euler_mixture_lw_new },
  { 0, 0 }
};

/* ************************************** */
/* General Relativistic Maxwell Equations */
/* ************************************** */

// GRMaxwell.new { lightSpeed = 1.0, elcErrorSpeedFactor = 0.0, mgnErrorSpeedFactor = 0.0, rpType = "roe" }
// where rpType is one of "roe" or "lax".
static int
eqn_gr_maxwell_lw_new(lua_State *L)
{
  struct wv_eqn_lw *gr_maxwell_lw = gkyl_malloc(sizeof(*gr_maxwell_lw));

  double light_speed = glua_tbl_get_number(L, "lightSpeed", 1.0);
  double e_fact = glua_tbl_get_number(L, "elcErrorSpeedFactor", 0.0);
  double b_fact = glua_tbl_get_number(L, "mgnErrorSpeedFactor", 0.0);

  const char *rp_str = glua_tbl_get_string(L, "rpType", "lax");
  enum gkyl_wv_gr_maxwell_rp rp_type = gkyl_search_str_int_pair_by_str(gr_maxwell_rp_type, rp_str, WV_GR_MAXWELL_RP_LAX);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  int reinit_freq = glua_tbl_get_integer(L, "reinitFreq", 100);

  with_lua_tbl_tbl(L, "blackHoleParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_blackhole_new(false, mass, spin, pos_x, pos_y, pos_z);
  }

  with_lua_tbl_tbl(L, "neutronStarParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double mass_quadrupole = glua_tbl_get_number(L, "massQuadrupole", 0.0);
    double spin_octupole = glua_tbl_get_number(L, "spinOctupole", 0.0);
    double mass_hexadecapole = glua_tbl_get_number(L, "massHexadecapole", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_neutronstar_new(false, mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z);
  }

  gr_maxwell_lw->magic = MOMENT_EQN_DEFAULT;
  gr_maxwell_lw->eqn = gkyl_wv_gr_maxwell_inew( &(struct gkyl_wv_gr_maxwell_inp) {
      .light_speed = light_speed,
      .e_fact = e_fact,
      .b_fact = b_fact,
      .spacetime = spacetime,
      .reinit_freq = reinit_freq,
      .rp_type = rp_type,
      .use_gpu = false,
    }
  );
  gr_maxwell_lw->has_nn = false;
  gr_maxwell_lw->ann = 0;
  gr_maxwell_lw->has_spacetime = true;
  gr_maxwell_lw->spacetime = spacetime;

  // Create Lua userdata.
  struct wv_eqn_lw **l_gr_maxwell_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_gr_maxwell_lw = gr_maxwell_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_gr_maxwell_ctor[] = {
  { "new", eqn_gr_maxwell_lw_new },
  { 0, 0 }
};

/* ********************************************************** */
/* General Relativistic Maxwell Equations in the Tetrad Basis */
/* ********************************************************** */

// GRMaxwellTetrad.new { lightSpeed = 1.0, elcErrorSpeedFactor = 0.0, mgnErrorSpeedFactor = 0.0, rpType = "roe" }
// where rpType is one of "roe" or "lax".
static int
eqn_gr_maxwell_tetrad_lw_new(lua_State *L)
{
  struct wv_eqn_lw *gr_maxwell_tetrad_lw = gkyl_malloc(sizeof(*gr_maxwell_tetrad_lw));

  double light_speed = glua_tbl_get_number(L, "lightSpeed", 1.0);
  double e_fact = glua_tbl_get_number(L, "elcErrorSpeedFactor", 0.0);
  double b_fact = glua_tbl_get_number(L, "mgnErrorSpeedFactor", 0.0);

  const char *rp_str = glua_tbl_get_string(L, "rpType", "lax");
  enum gkyl_wv_gr_maxwell_tetrad_rp rp_type = gkyl_search_str_int_pair_by_str(gr_maxwell_tetrad_rp_type, rp_str, WV_GR_MAXWELL_TETRAD_RP_LAX);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  int reinit_freq = glua_tbl_get_integer(L, "reinitFreq", 100);

  with_lua_tbl_tbl(L, "blackHoleParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_blackhole_new(false, mass, spin, pos_x, pos_y, pos_z);
  }

  with_lua_tbl_tbl(L, "neutronStarParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double mass_quadrupole = glua_tbl_get_number(L, "massQuadrupole", 0.0);
    double spin_octupole = glua_tbl_get_number(L, "spinOctupole", 0.0);
    double mass_hexadecapole = glua_tbl_get_number(L, "massHexadecapole", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_neutronstar_new(false, mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z);
  }

  gr_maxwell_tetrad_lw->magic = MOMENT_EQN_DEFAULT;
  gr_maxwell_tetrad_lw->eqn = gkyl_wv_gr_maxwell_tetrad_inew( &(struct gkyl_wv_gr_maxwell_tetrad_inp) {
      .light_speed = light_speed,
      .e_fact = e_fact,
      .b_fact = b_fact,
      .spacetime = spacetime,
      .reinit_freq = reinit_freq,
      .rp_type = rp_type,
      .use_gpu = false,
    }
  );
  gr_maxwell_tetrad_lw->has_nn = false;
  gr_maxwell_tetrad_lw->ann = 0;
  gr_maxwell_tetrad_lw->has_spacetime = true;
  gr_maxwell_tetrad_lw->spacetime = spacetime;

  // Create Lua userdata.
  struct wv_eqn_lw **l_gr_maxwell_tetrad_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_gr_maxwell_tetrad_lw = gr_maxwell_tetrad_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_gr_maxwell_tetrad_ctor[] = {
  { "new", eqn_gr_maxwell_tetrad_lw_new },
  { 0, 0 }
};

/* *************************************************************************** */
/* General Relativistic Euler Equations (Ultra-Relativistic Equation of State) */
/* *************************************************************************** */

// GRUltraRelativisticEuler.new { gasGamma = 4.0 / 3.0, rpType = "hll" }
// where rpType is one of "hll", "roe" or "lax".
static int
eqn_gr_ultra_rel_euler_lw_new(lua_State *L)
{
  struct wv_eqn_lw *gr_ultra_rel_euler_lw = gkyl_malloc(sizeof(*gr_ultra_rel_euler_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 4.0 / 3.0);

  const char *rp_str = glua_tbl_get_string(L, "rpType", "hll");
  enum gkyl_wv_gr_ultra_rel_euler_rp rp_type = gkyl_search_str_int_pair_by_str(gr_ultra_rel_euler_rp_type, rp_str, WV_GR_ULTRA_REL_EULER_RP_HLL);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  int reinit_freq = glua_tbl_get_integer(L, "reinitFreq", 100);

  with_lua_tbl_tbl(L, "blackHoleParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_blackhole_new(false, mass, spin, pos_x, pos_y, pos_z);
  }

  with_lua_tbl_tbl(L, "neutronStarParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double mass_quadrupole = glua_tbl_get_number(L, "massQuadrupole", 0.0);
    double spin_octupole = glua_tbl_get_number(L, "spinOctupole", 0.0);
    double mass_hexadecapole = glua_tbl_get_number(L, "massHexadecapole", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_neutronstar_new(false, mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z);
  }

  gr_ultra_rel_euler_lw->magic = MOMENT_EQN_DEFAULT;
  gr_ultra_rel_euler_lw->eqn = gkyl_wv_gr_ultra_rel_euler_inew( &(struct gkyl_wv_gr_ultra_rel_euler_inp) {
      .gas_gamma = gas_gamma,
      .spacetime = spacetime,
      .reinit_freq = reinit_freq,
      .rp_type = rp_type,
      .use_gpu = false,
    }
  );
  gr_ultra_rel_euler_lw->has_nn = false;
  gr_ultra_rel_euler_lw->ann = 0;
  gr_ultra_rel_euler_lw->has_spacetime = true;
  gr_ultra_rel_euler_lw->spacetime = spacetime;

  // Create Lua userdata.
  struct wv_eqn_lw **l_gr_ultra_rel_euler_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_gr_ultra_rel_euler_lw = gr_ultra_rel_euler_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_gr_ultra_rel_euler_ctor[] = {
  { "new", eqn_gr_ultra_rel_euler_lw_new },
  { 0, 0 }
};

/* *********************************************************************************************** */
/* General Relativistic Euler Equations in the Tetrad Basis (Ultra-Relativistic Equation of State) */
/* *********************************************************************************************** */

// GRUltraRelativisticEulerTetrad.new { gasGamma = 4.0 / 3.0, rpType = "hll" }
// where rpType is one of "hll", "roe" or "lax".
static int
eqn_gr_ultra_rel_euler_tetrad_lw_new(lua_State *L)
{
  struct wv_eqn_lw *gr_ultra_rel_euler_tetrad_lw = gkyl_malloc(sizeof(*gr_ultra_rel_euler_tetrad_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 4.0 / 3.0);

  const char *rp_str = glua_tbl_get_string(L, "rpType", "hll");
  enum gkyl_wv_gr_ultra_rel_euler_tetrad_rp rp_type = gkyl_search_str_int_pair_by_str(gr_ultra_rel_euler_tetrad_rp_type, rp_str, WV_GR_ULTRA_REL_EULER_TETRAD_RP_HLL);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  int reinit_freq = glua_tbl_get_integer(L, "reinitFreq", 100);

  with_lua_tbl_tbl(L, "blackHoleParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_blackhole_new(false, mass, spin, pos_x, pos_y, pos_z);
  }

  with_lua_tbl_tbl(L, "neutronStarParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double mass_quadrupole = glua_tbl_get_number(L, "massQuadrupole", 0.0);
    double spin_octupole = glua_tbl_get_number(L, "spinOctupole", 0.0);
    double mass_hexadecapole = glua_tbl_get_number(L, "massHexadecapole", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_neutronstar_new(false, mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z);
  }

  gr_ultra_rel_euler_tetrad_lw->magic = MOMENT_EQN_DEFAULT;
  gr_ultra_rel_euler_tetrad_lw->eqn = gkyl_wv_gr_ultra_rel_euler_tetrad_inew( &(struct gkyl_wv_gr_ultra_rel_euler_tetrad_inp) {
      .gas_gamma = gas_gamma,
      .spacetime = spacetime,
      .reinit_freq = reinit_freq,
      .rp_type = rp_type,
      .use_gpu = false,
    }
  );
  gr_ultra_rel_euler_tetrad_lw->has_nn = false;
  gr_ultra_rel_euler_tetrad_lw->ann = 0;
  gr_ultra_rel_euler_tetrad_lw->has_spacetime = true;
  gr_ultra_rel_euler_tetrad_lw->spacetime = spacetime;

  // Create Lua userdata.
  struct wv_eqn_lw **l_gr_ultra_rel_euler_tetrad_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_gr_ultra_rel_euler_tetrad_lw = gr_ultra_rel_euler_tetrad_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_gr_ultra_rel_euler_tetrad_ctor[] = {
  { "new", eqn_gr_ultra_rel_euler_tetrad_lw_new },
  { 0, 0 }
};

/* **************************************************************** */
/* General Relativistic Euler Equations (General Equation of State) */
/* **************************************************************** */

// GREuler.new { gasGamma = 5.0 / 3.0, rpType = "hll" }
// where rpType is one of "hll", "roe" or "lax".
static int
eqn_gr_euler_lw_new(lua_State *L)
{
  struct wv_eqn_lw *gr_euler_lw = gkyl_malloc(sizeof(*gr_euler_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 5.0 / 3.0);

  const char *rp_str = glua_tbl_get_string(L, "rpType", "hll");
  enum gkyl_wv_gr_euler_rp rp_type = gkyl_search_str_int_pair_by_str(gr_euler_rp_type, rp_str, WV_GR_EULER_RP_HLL);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  int reinit_freq = glua_tbl_get_integer(L, "reinitFreq", 100);

  with_lua_tbl_tbl(L, "blackHoleParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_blackhole_new(false, mass, spin, pos_x, pos_y, pos_z);
  }

  with_lua_tbl_tbl(L, "neutronStarParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double mass_quadrupole = glua_tbl_get_number(L, "massQuadrupole", 0.0);
    double spin_octupole = glua_tbl_get_number(L, "spinOctupole", 0.0);
    double mass_hexadecapole = glua_tbl_get_number(L, "massHexadecapole", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_neutronstar_new(false, mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z);
  }

  gr_euler_lw->magic = MOMENT_EQN_DEFAULT;
  gr_euler_lw->eqn = gkyl_wv_gr_euler_inew( &(struct gkyl_wv_gr_euler_inp) {
      .gas_gamma = gas_gamma,
      .spacetime = spacetime,
      .reinit_freq = reinit_freq,
      .rp_type = rp_type,
      .use_gpu = false,
    }
  );
  gr_euler_lw->has_nn = false;
  gr_euler_lw->ann = 0;
  gr_euler_lw->has_spacetime = true;
  gr_euler_lw->spacetime = spacetime;

  // Create Lua userdata.
  struct wv_eqn_lw **l_gr_euler_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_gr_euler_lw = gr_euler_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_gr_euler_ctor[] = {
  { "new", eqn_gr_euler_lw_new },
  { 0, 0 }
};

/* ************************************************************************************ */
/* General Relativistic Euler Equations in the Tetrad Basis (General Equation of State) */
/* ************************************************************************************ */

// GREulerTetrad.new { gasGamma = 5.0 / 3.0, rpType = "hll" }
// where rpType is one of "hll", "roe" or "lax".
static int
eqn_gr_euler_tetrad_lw_new(lua_State *L)
{
  struct wv_eqn_lw *gr_euler_tetrad_lw = gkyl_malloc(sizeof(*gr_euler_tetrad_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 5.0 / 3.0);

  const char *rp_str = glua_tbl_get_string(L, "rpType", "hll");
  enum gkyl_wv_gr_euler_tetrad_rp rp_type = gkyl_search_str_int_pair_by_str(gr_euler_tetrad_rp_type, rp_str, WV_GR_EULER_TETRAD_RP_HLL);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  int reinit_freq = glua_tbl_get_integer(L, "reinitFreq", 100);

  with_lua_tbl_tbl(L, "blackHoleParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_blackhole_new(false, mass, spin, pos_x, pos_y, pos_z);
  }

  with_lua_tbl_tbl(L, "neutronStarParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double mass_quadrupole = glua_tbl_get_number(L, "massQuadrupole", 0.0);
    double spin_octupole = glua_tbl_get_number(L, "spinOctupole", 0.0);
    double mass_hexadecapole = glua_tbl_get_number(L, "massHexadecapole", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_neutronstar_new(false, mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z);
  }

  gr_euler_tetrad_lw->magic = MOMENT_EQN_DEFAULT;
  gr_euler_tetrad_lw->eqn = gkyl_wv_gr_euler_tetrad_inew( &(struct gkyl_wv_gr_euler_tetrad_inp) {
      .gas_gamma = gas_gamma,
      .spacetime = spacetime,
      .reinit_freq = reinit_freq,
      .rp_type = rp_type,
      .use_gpu = false,
    }
  );
  gr_euler_tetrad_lw->has_nn = false;
  gr_euler_tetrad_lw->ann = 0;
  gr_euler_tetrad_lw->has_spacetime = true;
  gr_euler_tetrad_lw->spacetime = spacetime;

  // Create Lua userdata.
  struct wv_eqn_lw **l_gr_euler_tetrad_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_gr_euler_tetrad_lw = gr_euler_tetrad_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_gr_euler_tetrad_ctor[] = {
  { "new", eqn_gr_euler_tetrad_lw_new },
  { 0, 0 }
};

/* ******************************************************************* */
/* Coupled Fluid-Einstein Equations (Plane-Polarized Gowdy Spacetimes) */
/* ******************************************************************* */

// GRMedium.new { gasGamma = 4.0 / 3.0, kappa = 8.0 * pi, rpType = "lax" }.
static int
eqn_gr_medium_lw_new(lua_State *L)
{
  struct wv_eqn_lw *gr_medium_lw = gkyl_malloc(sizeof(*gr_medium_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 4.0 / 3.0);
  double kappa = glua_tbl_get_number(L, "kappa", 8.0 * M_PI);

  const char *rp_str = glua_tbl_get_string(L, "rpType", "lax");
  enum gkyl_wv_gr_medium_rp rp_type = gkyl_search_str_int_pair_by_str(gr_medium_rp_type, rp_str, WV_GR_MEDIUM_RP_LAX);

  gr_medium_lw->magic = MOMENT_EQN_DEFAULT;
  gr_medium_lw->eqn = gkyl_wv_gr_medium_inew( &(struct gkyl_wv_gr_medium_inp) {
      .gas_gamma = gas_gamma,
      .kappa = kappa,
      .rp_type = rp_type,
      .use_gpu = false,
    }
  );
  gr_medium_lw->has_nn = false;
  gr_medium_lw->ann = 0;
  gr_medium_lw->has_spacetime = false;
  gr_medium_lw->spacetime = 0;

  // Create Lua userdata.
  struct wv_eqn_lw **l_gr_medium_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_gr_medium_lw = gr_medium_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_gr_medium_ctor[] = {
  { "new", eqn_gr_medium_lw_new },
  { 0, 0 }
};

/* ************************************************ */
/* Vacuum Einstein Equations (Bona-Masso Formalism) */
/* ************************************************ */

// VacuumEinstein.new { excisionThreshold = 0.3, rpType = "hll" }.
static int
eqn_vacuum_einstein_lw_new(lua_State *L)
{
  struct wv_eqn_lw *vacuum_einstein_lw = gkyl_malloc(sizeof(*vacuum_einstein_lw));

  double excision_threshold = glua_tbl_get_number(L, "excisionThreshold", 0.3);

  const char *rp_str = glua_tbl_get_string(L, "rpType", "hll");
  enum gkyl_wv_vacuum_einstein_rp rp_type = gkyl_search_str_int_pair_by_str(vacuum_einstein_rp_type, rp_str, WV_VACUUM_EINSTEIN_RP_HLL);
  enum gkyl_spacetime_slicing spacetime_slicing = glua_tbl_get_integer(L, "spacetimeSlicing", GKYL_GEODESIC_SLICING);
  enum gkyl_spacetime_evolution spacetime_evolution = glua_tbl_get_integer(L, "spacetimeEvolution", GKYL_RICCI_EVOLUTION);

  vacuum_einstein_lw->magic = MOMENT_EQN_DEFAULT;
  vacuum_einstein_lw->eqn = gkyl_wv_vacuum_einstein_inew( &(struct gkyl_wv_vacuum_einstein_inp) {
      .excision_threshold = excision_threshold,
      .spacetime_slicing = spacetime_slicing,
      .spacetime_evolution = spacetime_evolution,
      .rp_type = rp_type,
      .use_gpu = false,
    }
  );
  vacuum_einstein_lw->has_nn = false;
  vacuum_einstein_lw->ann = 0;

  // Create Lua userdata.
  struct wv_eqn_lw **l_vacuum_einstein_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_vacuum_einstein_lw = vacuum_einstein_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_vacuum_einstein_ctor[] = {
  { "new", eqn_vacuum_einstein_lw_new },
  { 0, 0 }
};

/* ********************************************************** */
/* Vacuum Einstein Equations (Conformal Bona-Masso Formalism) */
/* ********************************************************** */

// VacuumEinsteinConformal.new { excisionThreshold = 0.3, rpType = "hll" }.
static int
eqn_vacuum_einstein_conformal_lw_new(lua_State *L)
{
  struct wv_eqn_lw *vacuum_einstein_conformal_lw = gkyl_malloc(sizeof(*vacuum_einstein_conformal_lw));

  double excision_threshold = glua_tbl_get_number(L, "excisionThreshold", 0.3);

  const char *rp_str = glua_tbl_get_string(L, "rpType", "hll");
  enum gkyl_wv_vacuum_einstein_conformal_rp rp_type = gkyl_search_str_int_pair_by_str(vacuum_einstein_conformal_rp_type, rp_str, WV_VACUUM_EINSTEIN_CONFORMAL_RP_HLL);
  enum gkyl_spacetime_slicing spacetime_slicing = glua_tbl_get_integer(L, "spacetimeSlicing", GKYL_GEODESIC_SLICING);
  enum gkyl_spacetime_evolution spacetime_evolution = glua_tbl_get_integer(L, "spacetimeEvolution", GKYL_RICCI_EVOLUTION);

  vacuum_einstein_conformal_lw->magic = MOMENT_EQN_DEFAULT;
  vacuum_einstein_conformal_lw->eqn = gkyl_wv_vacuum_einstein_conformal_inew( &(struct gkyl_wv_vacuum_einstein_conformal_inp) {
      .excision_threshold = excision_threshold,
      .spacetime_slicing = spacetime_slicing,
      .spacetime_evolution = spacetime_evolution,
      .rp_type = rp_type,
      .use_gpu = false,
    }
  );
  vacuum_einstein_conformal_lw->has_nn = false;
  vacuum_einstein_conformal_lw->ann = 0;

  // Create Lua userdata.
  struct wv_eqn_lw **l_vacuum_einstein_conformal_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_vacuum_einstein_conformal_lw = vacuum_einstein_conformal_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_vacuum_einstein_conformal_ctor[] = {
  { "new", eqn_vacuum_einstein_conformal_lw_new },
  { 0, 0 }
};

/* ******************************************************************** */
/* General Relativistic Two-Fluid Equations (General Equation of State) */
/* ******************************************************************** */

// GRTwoFluid.new { massElc = 1.0 / 1836.2, massIon = 1.0, chargeElc = -1.0, chargeIon = 1.0, gasGammaElc = 5.0 / 3.0, gasGammaIon = 5.0 / 3.0,
// lightSpeed = 1.0, elcErrorSpeedFactor = 0.0, mgnErrorSpeedFactor = 0.0, rpType = "hll" }
// where rpType is one of "hll" or "lax".
static int
eqn_gr_twofluid_lw_new(lua_State *L)
{
  struct wv_eqn_lw *gr_twofluid_lw = gkyl_malloc(sizeof(*gr_twofluid_lw));

  double mass_elc = glua_tbl_get_number(L, "massElc", 1.0 / 1836.2);
  double mass_ion = glua_tbl_get_number(L, "massIon", 1.0);
  double charge_elc = glua_tbl_get_number(L, "chargeElc", -1.0);
  double charge_ion = glua_tbl_get_number(L, "chargeIon", 1.0);
  double gas_gamma_elc = glua_tbl_get_number(L, "gasGammaElc", 5.0 / 3.0);
  double gas_gamma_ion = glua_tbl_get_number(L, "gasGammaIon", 5.0 / 3.0);

  double light_speed = glua_tbl_get_number(L, "lightSpeed", 1.0);
  double e_fact = glua_tbl_get_number(L, "elcErrorSpeedFactor", 0.0);
  double b_fact = glua_tbl_get_number(L, "mgnErrorSpeedFactor", 0.0);

  const char *rp_str = glua_tbl_get_string(L, "rpType", "hll");
  enum gkyl_wv_gr_twofluid_rp rp_type = gkyl_search_str_int_pair_by_str(gr_twofluid_rp_type, rp_str, WV_GR_TWOFLUID_RP_HLL);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  int reinit_freq = glua_tbl_get_integer(L, "reinitFreq", 100);

  with_lua_tbl_tbl(L, "blackHoleParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_blackhole_new(false, mass, spin, pos_x, pos_y, pos_z);
  }

  with_lua_tbl_tbl(L, "neutronStarParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double mass_quadrupole = glua_tbl_get_number(L, "massQuadrupole", 0.0);
    double spin_octupole = glua_tbl_get_number(L, "spinOctupole", 0.0);
    double mass_hexadecapole = glua_tbl_get_number(L, "massHexadecapole", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_neutronstar_new(false, mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z);
  }

  gr_twofluid_lw->magic = MOMENT_EQN_DEFAULT;
  gr_twofluid_lw->eqn = gkyl_wv_gr_twofluid_inew( &(struct gkyl_wv_gr_twofluid_inp) {
      .mass_elc = mass_elc,
      .mass_ion = mass_ion,
      .charge_elc = charge_elc,
      .charge_ion = charge_ion,
      .gas_gamma_elc = gas_gamma_elc,
      .gas_gamma_ion = gas_gamma_ion,
      .light_speed = light_speed,
      .e_fact = e_fact,
      .b_fact = b_fact,
      .spacetime = spacetime,
      .reinit_freq = reinit_freq,
      .rp_type = rp_type,
      .use_gpu = false,
    }
  );
  gr_twofluid_lw->has_nn = false;
  gr_twofluid_lw->ann = 0;
  gr_twofluid_lw->has_spacetime = true;
  gr_twofluid_lw->spacetime = spacetime;

  // Create Lua userdata.
  struct wv_eqn_lw **l_gr_twofluid_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_gr_twofluid_lw = gr_twofluid_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_gr_twofluid_ctor[] = {
  { "new", eqn_gr_twofluid_lw_new },
  { 0, 0 }
};

/* **************************************************************************************** */
/* General Relativistic Two-Fluid Equations in the Tetrad Basis (General Equation of State) */
/* **************************************************************************************** */

// GRTwoFluidTetrad.new { massElc = 1.0 / 1836.2, massIon = 1.0, chargeElc = -1.0, chargeIon = 1.0, gasGammaElc = 5.0 / 3.0, gasGammaIon = 5.0 / 3.0,
// lightSpeed = 1.0, elcErrorSpeedFactor = 0.0, mgnErrorSpeedFactor = 0.0, rpType = "hll" }
// where rpType is one of "hll" or "lax".
static int
eqn_gr_twofluid_tetrad_lw_new(lua_State *L)
{
  struct wv_eqn_lw *gr_twofluid_tetrad_lw = gkyl_malloc(sizeof(*gr_twofluid_tetrad_lw));

  double mass_elc = glua_tbl_get_number(L, "massElc", 1.0 / 1836.2);
  double mass_ion = glua_tbl_get_number(L, "massIon", 1.0);
  double charge_elc = glua_tbl_get_number(L, "chargeElc", -1.0);
  double charge_ion = glua_tbl_get_number(L, "chargeIon", 1.0);
  double gas_gamma_elc = glua_tbl_get_number(L, "gasGammaElc", 5.0 / 3.0);
  double gas_gamma_ion = glua_tbl_get_number(L, "gasGammaIon", 5.0 / 3.0);

  double light_speed = glua_tbl_get_number(L, "lightSpeed", 1.0);
  double e_fact = glua_tbl_get_number(L, "elcErrorSpeedFactor", 0.0);
  double b_fact = glua_tbl_get_number(L, "mgnErrorSpeedFactor", 0.0);

  const char *rp_str = glua_tbl_get_string(L, "rpType", "hll");
  enum gkyl_wv_gr_twofluid_tetrad_rp rp_type = gkyl_search_str_int_pair_by_str(gr_twofluid_tetrad_rp_type, rp_str, WV_GR_TWOFLUID_TETRAD_RP_HLL);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  int reinit_freq = glua_tbl_get_integer(L, "reinitFreq", 100);

  with_lua_tbl_tbl(L, "blackHoleParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_blackhole_new(false, mass, spin, pos_x, pos_y, pos_z);
  }

  with_lua_tbl_tbl(L, "neutronStarParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double mass_quadrupole = glua_tbl_get_number(L, "massQuadrupole", 0.0);
    double spin_octupole = glua_tbl_get_number(L, "spinOctupole", 0.0);
    double mass_hexadecapole = glua_tbl_get_number(L, "massHexadecapole", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_neutronstar_new(false, mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z);
  }

  gr_twofluid_tetrad_lw->magic = MOMENT_EQN_DEFAULT;
  gr_twofluid_tetrad_lw->eqn = gkyl_wv_gr_twofluid_tetrad_inew( &(struct gkyl_wv_gr_twofluid_tetrad_inp) {
      .mass_elc = mass_elc,
      .mass_ion = mass_ion,
      .charge_elc = charge_elc,
      .charge_ion = charge_ion,
      .gas_gamma_elc = gas_gamma_elc,
      .gas_gamma_ion = gas_gamma_ion,
      .light_speed = light_speed,
      .e_fact = e_fact,
      .b_fact = b_fact,
      .spacetime = spacetime,
      .reinit_freq = reinit_freq,
      .rp_type = rp_type,
      .use_gpu = false,
    }
  );
  gr_twofluid_tetrad_lw->has_nn = false;
  gr_twofluid_tetrad_lw->ann = 0;
  gr_twofluid_tetrad_lw->has_spacetime = true;
  gr_twofluid_tetrad_lw->spacetime = spacetime;

  // Create Lua userdata.
  struct wv_eqn_lw **l_gr_twofluid_tetrad_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_gr_twofluid_tetrad_lw = gr_twofluid_tetrad_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_gr_twofluid_tetrad_ctor[] = {
  { "new", eqn_gr_twofluid_tetrad_lw_new },
  { 0, 0 }
};

/* ********************************************************************************* */
/* General Relativistic Magnetohydrodynamics Equations (Ideal Gas Equation of State) */
/* ********************************************************************************* */

// GRMHD.new { gasGamma = 5.0 / 3.0, lightSpeed = 1.0, mgnErrorSpeedFactor = 0.0, rpType = "hll" }
// where rpType is one of "hll" or "lax".
static int
eqn_gr_mhd_lw_new(lua_State *L)
{
  struct wv_eqn_lw *gr_mhd_lw = gkyl_malloc(sizeof(*gr_mhd_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 5.0 / 3.0);
  double light_speed = glua_tbl_get_number(L, "lightSpeed", 1.0);
  double b_fact = glua_tbl_get_number(L, "mgnErrorSpeedFactor", 0.0);

  const char *rp_str = glua_tbl_get_string(L, "rpType", "hll");
  enum gkyl_wv_gr_mhd_rp rp_type = gkyl_search_str_int_pair_by_str(gr_mhd_rp_type, rp_str, WV_GR_MHD_RP_HLL);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  int reinit_freq = glua_tbl_get_integer(L, "reinitFreq", 100);

  with_lua_tbl_tbl(L, "blackHoleParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_blackhole_new(false, mass, spin, pos_x, pos_y, pos_z);
  }

  with_lua_tbl_tbl(L, "neutronStarParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double mass_quadrupole = glua_tbl_get_number(L, "massQuadrupole", 0.0);
    double spin_octupole = glua_tbl_get_number(L, "spinOctupole", 0.0);
    double mass_hexadecapole = glua_tbl_get_number(L, "massHexadecapole", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_neutronstar_new(false, mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z);
  }

  int spacetime_gauge = glua_tbl_get_integer(L, "spacetimeGauge", GKYL_STATIC_GAUGE);

  gr_mhd_lw->magic = MOMENT_EQN_DEFAULT;
  gr_mhd_lw->eqn = gkyl_wv_gr_mhd_inew( &(struct gkyl_wv_gr_mhd_inp) {
      .gas_gamma = gas_gamma,
      .light_speed = light_speed,
      .b_fact = b_fact,
      .spacetime = spacetime,
      .spacetime_gauge = spacetime_gauge,
      .reinit_freq = reinit_freq,
      .rp_type = rp_type,
      .use_gpu = false,
    }
  );
  gr_mhd_lw->has_nn = false;
  gr_mhd_lw->ann = 0;
  gr_mhd_lw->has_spacetime = true;
  gr_mhd_lw->spacetime = spacetime;

  // Create Lua userdata.
  struct wv_eqn_lw **l_gr_mhd_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_gr_mhd_lw = gr_mhd_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_gr_mhd_ctor[] = {
  { "new", eqn_gr_mhd_lw_new },
  { 0, 0 }
};

/* ***************************************************************************************************** */
/* General Relativistic Magnetohydrodynamics Equations in the Tetrad Basis (Ideal Gas Equation of State) */
/* ***************************************************************************************************** */

// GRMHDTetrad.new { gasGamma = 5.0 / 3.0, lightSpeed = 1.0, mgnErrorSpeedFactor = 0.0, rpType = "hll" }
// where rpType is one of "hll" or "lax".
static int
eqn_gr_mhd_tetrad_lw_new(lua_State *L)
{
  struct wv_eqn_lw *gr_mhd_tetrad_lw = gkyl_malloc(sizeof(*gr_mhd_tetrad_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 5.0 / 3.0);
  double light_speed = glua_tbl_get_number(L, "lightSpeed", 1.0);
  double b_fact = glua_tbl_get_number(L, "mgnErrorSpeedFactor", 0.0);

  const char *rp_str = glua_tbl_get_string(L, "rpType", "hll");
  enum gkyl_wv_gr_mhd_tetrad_rp rp_type = gkyl_search_str_int_pair_by_str(gr_mhd_tetrad_rp_type, rp_str, WV_GR_MHD_TETRAD_RP_HLL);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  int reinit_freq = glua_tbl_get_integer(L, "reinitFreq", 100);

  with_lua_tbl_tbl(L, "blackHoleParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_blackhole_new(false, mass, spin, pos_x, pos_y, pos_z);
  }

  with_lua_tbl_tbl(L, "neutronStarParameters") {
    double mass = glua_tbl_get_number(L, "mass", 1.0);
    double spin = glua_tbl_get_number(L, "spin", 0.0);
    double mass_quadrupole = glua_tbl_get_number(L, "massQuadrupole", 0.0);
    double spin_octupole = glua_tbl_get_number(L, "spinOctupole", 0.0);
    double mass_hexadecapole = glua_tbl_get_number(L, "massHexadecapole", 0.0);
    double pos_x = glua_tbl_get_number(L, "posX", 0.0);
    double pos_y = glua_tbl_get_number(L, "posY", 0.0);
    double pos_z = glua_tbl_get_number(L, "posZ", 0.0);

    spacetime = gkyl_gr_neutronstar_new(false, mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z);
  }

  int spacetime_gauge = glua_tbl_get_integer(L, "spacetimeGauge", GKYL_STATIC_GAUGE);

  gr_mhd_tetrad_lw->magic = MOMENT_EQN_DEFAULT;
  gr_mhd_tetrad_lw->eqn = gkyl_wv_gr_mhd_tetrad_inew( &(struct gkyl_wv_gr_mhd_tetrad_inp) {
      .gas_gamma = gas_gamma,
      .light_speed = light_speed,
      .b_fact = b_fact,
      .spacetime = spacetime,
      .spacetime_gauge = spacetime_gauge,
      .reinit_freq = reinit_freq,
      .rp_type = rp_type,
      .use_gpu = false,
    }
  );
  gr_mhd_tetrad_lw->has_nn = false;
  gr_mhd_tetrad_lw->ann = 0;
  gr_mhd_tetrad_lw->has_spacetime = true;
  gr_mhd_tetrad_lw->spacetime = spacetime;

  // Create Lua userdata.
  struct wv_eqn_lw **l_gr_mhd_tetrad_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_gr_mhd_tetrad_lw = gr_mhd_tetrad_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_gr_mhd_tetrad_ctor[] = {
  { "new", eqn_gr_mhd_tetrad_lw_new },
  { 0, 0 }
};

/* ************************* */
/* Linear Advection Equation */
/* ************************* */

// Advection.new { advectionSpeed = 1.0 }
static int
eqn_advect_lw_new(lua_State *L)
{
  struct wv_eqn_lw *advect_lw = gkyl_malloc(sizeof(*advect_lw));

  double c = glua_tbl_get_number(L, "advectionSpeed", 1.0);

  advect_lw->magic = MOMENT_EQN_DEFAULT;
  advect_lw->eqn = gkyl_wv_advect_new(c, false);
  advect_lw->has_nn = false;
  advect_lw->ann = 0;
  advect_lw->has_spacetime = false;
  advect_lw->spacetime = 0;

  // Create Lua userdata.
  struct wv_eqn_lw **l_advect_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_advect_lw = advect_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_advect_ctor[] = {
  { "new", eqn_advect_lw_new },
  { 0, 0 }
};

/* ************************* */
/* Inviscid Burgers Equation */
/* ************************* */

// Burgers.new { }
static int
eqn_burgers_lw_new(lua_State *L)
{
  struct wv_eqn_lw *burgers_lw = gkyl_malloc(sizeof(*burgers_lw));

  burgers_lw->magic = MOMENT_EQN_DEFAULT;
  burgers_lw->eqn = gkyl_wv_burgers_new(false);
  burgers_lw->has_nn = false;
  burgers_lw->ann = 0;
  burgers_lw->has_spacetime = false;
  burgers_lw->spacetime = 0;

  // Create Lua userdata.
  struct wv_eqn_lw **l_burgers_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_burgers_lw = burgers_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_burgers_ctor[] = {
  { "new", eqn_burgers_lw_new },
  { 0, 0 }
};

// Register and load all wave equation objects.
static void
eqn_openlibs(lua_State *L)
{
  luaL_newmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);

  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, wv_eqn_lw_gc);
  lua_settable(L, -3);

  luaL_register(L, "G0.Moments.Eq.Euler", eqn_euler_ctor);
  luaL_register(L, "G0.Moments.Eq.IsoEuler", eqn_iso_euler_ctor);
  luaL_register(L, "G0.Moments.Eq.SrEuler", eqn_sr_euler_ctor);
  luaL_register(L, "G0.Moments.Eq.ColdFluid", eqn_coldfluid_ctor);
  luaL_register(L, "G0.Moments.Eq.TenMoment", eqn_tenmoment_ctor); 
  luaL_register(L, "G0.Moments.Eq.MHD", eqn_mhd_ctor);
  luaL_register(L, "G0.Moments.Eq.ReactiveEuler", eqn_reactive_euler_ctor);
  luaL_register(L, "G0.Moments.Eq.EulerMixture", eqn_euler_mixture_ctor);
  luaL_register(L, "G0.Moments.Eq.EulerRGFM", eqn_euler_rgfm_ctor);
  luaL_register(L, "G0.Moments.Eq.IsoEulerMixture", eqn_iso_euler_mixture_ctor);
  luaL_register(L, "G0.Moments.Eq.GRMaxwell", eqn_gr_maxwell_ctor);
  luaL_register(L, "G0.Moments.Eq.GRMaxwellTetrad", eqn_gr_maxwell_tetrad_ctor);
  luaL_register(L, "G0.Moments.Eq.GRUltraRelEuler", eqn_gr_ultra_rel_euler_ctor);
  luaL_register(L, "G0.Moments.Eq.GRUltraRelEulerTetrad", eqn_gr_ultra_rel_euler_tetrad_ctor);
  luaL_register(L, "G0.Moments.Eq.GREuler", eqn_gr_euler_ctor);
  luaL_register(L, "G0.Moments.Eq.GREulerTetrad", eqn_gr_euler_tetrad_ctor);
  luaL_register(L, "G0.Moments.Eq.GRMedium", eqn_gr_medium_ctor);
  luaL_register(L, "G0.Moments.Eq.VacuumEinstein", eqn_vacuum_einstein_ctor);
  luaL_register(L, "G0.Moments.Eq.VacuumEinsteinConformal", eqn_vacuum_einstein_conformal_ctor);
  luaL_register(L, "G0.Moments.Eq.GRTwoFluid", eqn_gr_twofluid_ctor);
  luaL_register(L, "G0.Moments.Eq.GRTwoFluidTetrad", eqn_gr_twofluid_tetrad_ctor);
  luaL_register(L, "G0.Moments.Eq.GRMHD", eqn_gr_mhd_ctor);
  luaL_register(L, "G0.Moments.Eq.GRMHDTetrad", eqn_gr_mhd_tetrad_ctor);
  luaL_register(L, "G0.Moments.Eq.LinearAdvection", eqn_advect_ctor);
  luaL_register(L, "G0.Moments.Eq.Burgers", eqn_burgers_ctor);
}

// Metatable name for spacetime object input struct.
#define MOMENT_SPACETIME_METATABLE_NM "GkeyllZero.App.Moments.Spacetime"

// Methods for manipulating gkyl_gr_spacetime objects.

// Lua userdata object for constructing spacetime objects.
struct gr_spacetime_lw {
  int magic; // This must be the first element in the struct.
  struct gkyl_gr_spacetime *spacetime; // Spacetime object.
};

// Clean up memory allocated for spacetime object.
static int
gr_spacetime_lw_gc(lua_State *L)
{
  struct gr_spacetime_lw **l_gr_lw = GKYL_CHECK_UDATA(L, MOMENT_SPACETIME_METATABLE_NM);
  struct gr_spacetime_lw *gr_lw = *l_gr_lw;

  gkyl_gr_spacetime_release(gr_lw->spacetime);
  gkyl_free(*l_gr_lw);

  return 0;
}

// Acquire spacetime object.
static struct gkyl_gr_spacetime*
gr_spacetime_get(lua_State *L)
{
  struct gr_spacetime_lw **l_gr_lw = luaL_checkudata(L, -1, MOMENT_SPACETIME_METATABLE_NM);
  struct gr_spacetime_lw *gr_lw = *l_gr_lw;

  return gr_lw->spacetime;
}

/* ******************* */
/* Minkowski Spacetime */
/* ******************* */

static int
spacetime_minkowski_lw_new(lua_State *L)
{
  struct gr_spacetime_lw *minkowski_lw = gkyl_malloc(sizeof(*minkowski_lw));

  minkowski_lw->magic = MOMENT_SPACETIME_DEFAULT;
  minkowski_lw->spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  // Create Lua userdata.
  struct gr_spacetime_lw **l_minkowski_lw = lua_newuserdata(L, sizeof(struct gr_spacetime_lw*));
  *l_minkowski_lw = minkowski_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_SPACETIME_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

static int
spacetime_minkowski_lw_spatial_metric_tensor(lua_State *L)
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 1);
  double x = luaL_checknumber(L, 2);
  double y = luaL_checknumber(L, 3);
  double z = luaL_checknumber(L, 4);

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->spatial_metric_tensor_func(spacetime, t, x, y, z, &spatial_metric);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, spatial_metric[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
  }
  gkyl_free(spatial_metric);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_minkowski_lw_inv_spatial_metric_tensor(lua_State *L)
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 1);
  double x = luaL_checknumber(L, 2);
  double y = luaL_checknumber(L, 3);
  double z = luaL_checknumber(L, 4);

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->spatial_inv_metric_tensor_func(spacetime, t, x, y, z, &inv_spatial_metric);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, inv_spatial_metric[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(inv_spatial_metric);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_minkowski_lw_spatial_metric_det(lua_State *L)
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 1);
  double x = luaL_checknumber(L, 2);
  double y = luaL_checknumber(L, 3);
  double z = luaL_checknumber(L, 4);

  double spatial_metric_det;
  spacetime->spatial_metric_det_func(spacetime, t, x, y, z, &spatial_metric_det);

  lua_pushnumber(L, spatial_metric_det);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_minkowski_lw_lapse_function(lua_State *L)
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 1);
  double x = luaL_checknumber(L, 2);
  double y = luaL_checknumber(L, 3);
  double z = luaL_checknumber(L, 4);

  double lapse;
  spacetime->lapse_function_func(spacetime, t, x, y, z, &lapse);

  lua_pushnumber(L, lapse);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_minkowski_lw_shift_vector(lua_State *L)
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 1);
  double x = luaL_checknumber(L, 2);
  double y = luaL_checknumber(L, 3);
  double z = luaL_checknumber(L, 4);

  double *shift = gkyl_malloc(sizeof(double[3]));
  spacetime->shift_vector_func(spacetime, t, x, y, z, &shift);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, shift[i]);
    lua_rawset(L, -3);
  }

  gkyl_free(shift);
  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_minkowski_lw_extrinsic_curvature_tensor(lua_State *L)
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 1);
  double x = luaL_checknumber(L, 2);
  double y = luaL_checknumber(L, 3);
  double z = luaL_checknumber(L, 4);

  double dx = luaL_checknumber(L, 5);
  double dy = luaL_checknumber(L, 6);
  double dz = luaL_checknumber(L, 7);

  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->extrinsic_curvature_tensor_func(spacetime, t, x, y, z, dx, dy, dz, &extrinsic_curvature);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, extrinsic_curvature[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(extrinsic_curvature[i]);
  }
  gkyl_free(extrinsic_curvature);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_minkowski_lw_lapse_function_der(lua_State *L)
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 1);
  double x = luaL_checknumber(L, 2);
  double y = luaL_checknumber(L, 3);
  double z = luaL_checknumber(L, 4);

  double dx = luaL_checknumber(L, 5);
  double dy = luaL_checknumber(L, 6);
  double dz = luaL_checknumber(L, 7);

  double *lapse_der = gkyl_malloc(sizeof(double[3]));
  spacetime->lapse_function_der_func(spacetime, t, x, y, z, dx, dy, dz, &lapse_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, lapse_der[i]);
    lua_rawset(L, -3);
  }

  gkyl_free(lapse_der);
  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_minkowski_lw_shift_vector_der(lua_State *L)
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 1);
  double x = luaL_checknumber(L, 2);
  double y = luaL_checknumber(L, 3);
  double z = luaL_checknumber(L, 4);

  double dx = luaL_checknumber(L, 5);
  double dy = luaL_checknumber(L, 6);
  double dz = luaL_checknumber(L, 7);

  double **shift_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    shift_der[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->shift_vector_der_func(spacetime, t, x, y, z, dx, dy, dz, &shift_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, shift_der[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(shift_der[i]);
  }
  gkyl_free(shift_der);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_minkowski_lw_spatial_metric_tensor_der(lua_State *L)
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 1);
  double x = luaL_checknumber(L, 2);
  double y = luaL_checknumber(L, 3);
  double z = luaL_checknumber(L, 4);

  double dx = luaL_checknumber(L, 5);
  double dy = luaL_checknumber(L, 6);
  double dz = luaL_checknumber(L, 7);

  double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

    for (int j = 0; j < 3; j++) {
      spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
    }
  }

  spacetime->spatial_metric_tensor_der_func(spacetime, t, x, y, z, dx, dy, dz, &spatial_metric_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);

      lua_createtable(L, 3, 0);
      for (int k = 0; k < 3; k++) {
        lua_pushinteger(L, k + 1);
        lua_pushnumber(L, spatial_metric_der[i][j][k]);
        lua_rawset(L, -3);
      }

      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      gkyl_free(spatial_metric_der[i][j]);
    }
    gkyl_free(spatial_metric_der[i]);
  }
  gkyl_free(spatial_metric_der);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_minkowski_lw_excision_region(lua_State *L)
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 1);
  double x = luaL_checknumber(L, 2);
  double y = luaL_checknumber(L, 3);
  double z = luaL_checknumber(L, 4);

  bool in_excision_region;
  spacetime->excision_region_func(spacetime, t, x, y, z, &in_excision_region);

  lua_pushboolean(L, in_excision_region);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_minkowski_lw_conformal_factor(lua_State *L)
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 1);
  double x = luaL_checknumber(L, 2);
  double y = luaL_checknumber(L, 3);
  double z = luaL_checknumber(L, 4);

  double conformal_factor;
  spacetime->conformal_factor_func(spacetime, t, x, y, z, &conformal_factor);

  lua_pushnumber(L, conformal_factor);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_minkowski_lw_bssn_conformal_factor(lua_State *L)
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 1);
  double x = luaL_checknumber(L, 2);
  double y = luaL_checknumber(L, 3);
  double z = luaL_checknumber(L, 4);

  double bssn_conformal_factor;
  spacetime->bssn_conformal_factor_func(spacetime, t, x, y, z, &bssn_conformal_factor);

  lua_pushnumber(L, bssn_conformal_factor);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_minkowski_lw_conformal_factor_der(lua_State *L)
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 1);
  double x = luaL_checknumber(L, 2);
  double y = luaL_checknumber(L, 3);
  double z = luaL_checknumber(L, 4);

  double dx = luaL_checknumber(L, 5);
  double dy = luaL_checknumber(L, 6);
  double dz = luaL_checknumber(L, 7);

  double *conformal_factor_der = gkyl_malloc(sizeof(double[3]));
  spacetime->conformal_factor_der_func(spacetime, t, x, y, z, dx, dy, dz, &conformal_factor_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, conformal_factor_der[i]);
    lua_rawset(L, -3);
  }

  gkyl_free(conformal_factor_der);
  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_minkowski_lw_bssn_conformal_factor_der(lua_State *L)
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 1);
  double x = luaL_checknumber(L, 2);
  double y = luaL_checknumber(L, 3);
  double z = luaL_checknumber(L, 4);

  double dx = luaL_checknumber(L, 5);
  double dy = luaL_checknumber(L, 6);
  double dz = luaL_checknumber(L, 7);

  double *bssn_conformal_factor_der = gkyl_malloc(sizeof(double[3]));
  spacetime->bssn_conformal_factor_der_func(spacetime, t, x, y, z, dx, dy, dz, &bssn_conformal_factor_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, bssn_conformal_factor_der[i]);
    lua_rawset(L, -3);
  }

  gkyl_free(bssn_conformal_factor_der);
  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_minkowski_lw_bssn_conformal_factor_der2(lua_State *L)
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_inew( &(struct gkyl_gr_minkowski_inp) {
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 1);
  double x = luaL_checknumber(L, 2);
  double y = luaL_checknumber(L, 3);
  double z = luaL_checknumber(L, 4);

  double dx = luaL_checknumber(L, 5);
  double dy = luaL_checknumber(L, 6);
  double dz = luaL_checknumber(L, 7);

  double **bssn_conformal_factor_der2 = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    bssn_conformal_factor_der2[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->bssn_conformal_factor_der2_func(spacetime, t, x, y, z, dx, dy, dz, &bssn_conformal_factor_der2);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, bssn_conformal_factor_der2[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(bssn_conformal_factor_der2[i]);
  }
  gkyl_free(bssn_conformal_factor_der2);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

// Spacetime constructor.
static struct luaL_Reg spacetime_minkowski_ctor[] = {
  { "new", spacetime_minkowski_lw_new },
  { "spatialMetricTensor", spacetime_minkowski_lw_spatial_metric_tensor },
  { "invSpatialMetricTensor", spacetime_minkowski_lw_inv_spatial_metric_tensor },
  { "spatialMetricDeterminant", spacetime_minkowski_lw_spatial_metric_det },
  { "lapseFunction", spacetime_minkowski_lw_lapse_function },
  { "shiftVector", spacetime_minkowski_lw_shift_vector },
  { "extrinsicCurvatureTensor", spacetime_minkowski_lw_extrinsic_curvature_tensor },
  { "excisionRegion", spacetime_minkowski_lw_excision_region },
  { "lapseFunctionDer", spacetime_minkowski_lw_lapse_function_der },
  { "shiftVectorDer", spacetime_minkowski_lw_shift_vector_der },
  { "spatialMetricTensorDer", spacetime_minkowski_lw_spatial_metric_tensor_der },
  { "conformalFactor", spacetime_minkowski_lw_conformal_factor },
  { "bssnConformalFactor", spacetime_minkowski_lw_bssn_conformal_factor },
  { "conformalFactorDer", spacetime_minkowski_lw_conformal_factor_der },
  { "bssnConformalFactorDer", spacetime_minkowski_lw_bssn_conformal_factor_der },
  { "bssnConformalFactorDer2", spacetime_minkowski_lw_bssn_conformal_factor_der2 },
  { 0, 0 }
};

/* ******************** */
/* Black Hole Spacetime */
/* ******************** */

static int
spacetime_blackhole_lw_new(lua_State *L)
{
  struct gr_spacetime_lw *blackhole_lw = gkyl_malloc(sizeof(*blackhole_lw));

  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  blackhole_lw->magic = MOMENT_SPACETIME_DEFAULT;
  blackhole_lw->spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  // Create Lua userdata.
  struct gr_spacetime_lw **l_blackhole_lw = lua_newuserdata(L, sizeof(struct gr_spacetime_lw*));
  *l_blackhole_lw = blackhole_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_SPACETIME_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

static int
spacetime_blackhole_lw_spatial_metric_tensor(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 6);
  double x = luaL_checknumber(L, 7);
  double y = luaL_checknumber(L, 8);
  double z = luaL_checknumber(L, 9);

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->spatial_metric_tensor_func(spacetime, t, x, y, z, &spatial_metric);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, spatial_metric[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
  }
  gkyl_free(spatial_metric);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_blackhole_lw_inv_spatial_metric_tensor(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 6);
  double x = luaL_checknumber(L, 7);
  double y = luaL_checknumber(L, 8);
  double z = luaL_checknumber(L, 9);

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->spatial_inv_metric_tensor_func(spacetime, t, x, y, z, &inv_spatial_metric);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, inv_spatial_metric[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(inv_spatial_metric);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_blackhole_lw_spatial_metric_det(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 6);
  double x = luaL_checknumber(L, 7);
  double y = luaL_checknumber(L, 8);
  double z = luaL_checknumber(L, 9);

  double spatial_metric_det;
  spacetime->spatial_metric_det_func(spacetime, t, x, y, z, &spatial_metric_det);

  lua_pushnumber(L, spatial_metric_det);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_blackhole_lw_lapse_function(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 6);
  double x = luaL_checknumber(L, 7);
  double y = luaL_checknumber(L, 8);
  double z = luaL_checknumber(L, 9);

  double lapse;
  spacetime->lapse_function_func(spacetime, t, x, y, z, &lapse);

  lua_pushnumber(L, lapse);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_blackhole_lw_shift_vector(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 6);
  double x = luaL_checknumber(L, 7);
  double y = luaL_checknumber(L, 8);
  double z = luaL_checknumber(L, 9);

  double *shift = gkyl_malloc(sizeof(double[3]));
  spacetime->shift_vector_func(spacetime, t, x, y, z, &shift);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, shift[i]);
    lua_rawset(L, -3);
  }

  gkyl_free(shift);
  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_blackhole_lw_extrinsic_curvature_tensor(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 6);
  double x = luaL_checknumber(L, 7);
  double y = luaL_checknumber(L, 8);
  double z = luaL_checknumber(L, 9);

  double dx = luaL_checknumber(L, 10);
  double dy = luaL_checknumber(L, 11);
  double dz = luaL_checknumber(L, 12);

  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->extrinsic_curvature_tensor_func(spacetime, t, x, y, z, dx, dy, dz, &extrinsic_curvature);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, extrinsic_curvature[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(extrinsic_curvature[i]);
  }
  gkyl_free(extrinsic_curvature);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_blackhole_lw_lapse_function_der(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 6);
  double x = luaL_checknumber(L, 7);
  double y = luaL_checknumber(L, 8);
  double z = luaL_checknumber(L, 9);

  double dx = luaL_checknumber(L, 10);
  double dy = luaL_checknumber(L, 11);
  double dz = luaL_checknumber(L, 12);

  double *lapse_der = gkyl_malloc(sizeof(double[3]));
  spacetime->lapse_function_der_func(spacetime, t, x, y, z, dx, dy, dz, &lapse_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, lapse_der[i]);
    lua_rawset(L, -3);
  }

  gkyl_free(lapse_der);
  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_blackhole_lw_shift_vector_der(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 6);
  double x = luaL_checknumber(L, 7);
  double y = luaL_checknumber(L, 8);
  double z = luaL_checknumber(L, 9);

  double dx = luaL_checknumber(L, 10);
  double dy = luaL_checknumber(L, 11);
  double dz = luaL_checknumber(L, 12);

  double **shift_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    shift_der[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->shift_vector_der_func(spacetime, t, x, y, z, dx, dy, dz, &shift_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, shift_der[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(shift_der[i]);
  }
  gkyl_free(shift_der);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_blackhole_lw_spatial_metric_tensor_der(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 6);
  double x = luaL_checknumber(L, 7);
  double y = luaL_checknumber(L, 8);
  double z = luaL_checknumber(L, 9);

  double dx = luaL_checknumber(L, 10);
  double dy = luaL_checknumber(L, 11);
  double dz = luaL_checknumber(L, 12);

  double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

    for (int j = 0; j < 3; j++) {
      spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
    }
  }

  spacetime->spatial_metric_tensor_der_func(spacetime, t, x, y, z, dx, dy, dz, &spatial_metric_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);

      lua_createtable(L, 3, 0);
      for (int k = 0; k < 3; k++) {
        lua_pushinteger(L, k + 1);
        lua_pushnumber(L, spatial_metric_der[i][j][k]);
        lua_rawset(L, -3);
      }

      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      gkyl_free(spatial_metric_der[i][j]);
    }
    gkyl_free(spatial_metric_der[i]);
  }
  gkyl_free(spatial_metric_der);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_blackhole_lw_excision_region(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 6);
  double x = luaL_checknumber(L, 7);
  double y = luaL_checknumber(L, 8);
  double z = luaL_checknumber(L, 9);

  bool in_excision_region;
  spacetime->excision_region_func(spacetime, t, x, y, z, &in_excision_region);

  lua_pushboolean(L, in_excision_region);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_blackhole_lw_conformal_factor(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 6);
  double x = luaL_checknumber(L, 7);
  double y = luaL_checknumber(L, 8);
  double z = luaL_checknumber(L, 9);

  double conformal_factor;
  spacetime->conformal_factor_func(spacetime, t, x, y, z, &conformal_factor);

  lua_pushnumber(L, conformal_factor);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_blackhole_lw_bssn_conformal_factor(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 6);
  double x = luaL_checknumber(L, 7);
  double y = luaL_checknumber(L, 8);
  double z = luaL_checknumber(L, 9);

  double bssn_conformal_factor;
  spacetime->bssn_conformal_factor_func(spacetime, t, x, y, z, &bssn_conformal_factor);

  lua_pushnumber(L, bssn_conformal_factor);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_blackhole_lw_conformal_factor_der(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 6);
  double x = luaL_checknumber(L, 7);
  double y = luaL_checknumber(L, 8);
  double z = luaL_checknumber(L, 9);

  double dx = luaL_checknumber(L, 10);
  double dy = luaL_checknumber(L, 11);
  double dz = luaL_checknumber(L, 12);

  double *conformal_factor_der = gkyl_malloc(sizeof(double[3]));
  spacetime->conformal_factor_der_func(spacetime, t, x, y, z, dx, dy, dz, &conformal_factor_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, conformal_factor_der[i]);
    lua_rawset(L, -3);
  }

  gkyl_free(conformal_factor_der);
  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_blackhole_lw_bssn_conformal_factor_der(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 6);
  double x = luaL_checknumber(L, 7);
  double y = luaL_checknumber(L, 8);
  double z = luaL_checknumber(L, 9);

  double dx = luaL_checknumber(L, 10);
  double dy = luaL_checknumber(L, 11);
  double dz = luaL_checknumber(L, 12);

  double *bssn_conformal_factor_der = gkyl_malloc(sizeof(double[3]));
  spacetime->bssn_conformal_factor_der_func(spacetime, t, x, y, z, dx, dy, dz, &bssn_conformal_factor_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, bssn_conformal_factor_der[i]);
    lua_rawset(L, -3);
  }

  gkyl_free(bssn_conformal_factor_der);
  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_blackhole_lw_bssn_conformal_factor_der2(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double pos_x = luaL_checknumber(L, 3);
  double pos_y = luaL_checknumber(L, 4);
  double pos_z = luaL_checknumber(L, 5);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_inew( &(struct gkyl_gr_blackhole_inp) {
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 6);
  double x = luaL_checknumber(L, 7);
  double y = luaL_checknumber(L, 8);
  double z = luaL_checknumber(L, 9);

  double dx = luaL_checknumber(L, 10);
  double dy = luaL_checknumber(L, 11);
  double dz = luaL_checknumber(L, 12);

  double **bssn_conformal_factor_der2 = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    bssn_conformal_factor_der2[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->bssn_conformal_factor_der2_func(spacetime, t, x, y, z, dx, dy, dz, &bssn_conformal_factor_der2);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, bssn_conformal_factor_der2[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(bssn_conformal_factor_der2[i]);
  }
  gkyl_free(bssn_conformal_factor_der2);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

// Spacetime constructor.
static struct luaL_Reg spacetime_blackhole_ctor[] = {
  { "new", spacetime_blackhole_lw_new },
  { "spatialMetricTensor", spacetime_blackhole_lw_spatial_metric_tensor },
  { "invSpatialMetricTensor", spacetime_blackhole_lw_inv_spatial_metric_tensor },
  { "spatialMetricDeterminant", spacetime_blackhole_lw_spatial_metric_det },
  { "lapseFunction", spacetime_blackhole_lw_lapse_function },
  { "shiftVector", spacetime_blackhole_lw_shift_vector },
  { "extrinsicCurvatureTensor", spacetime_blackhole_lw_extrinsic_curvature_tensor },
  { "lapseFunctionDer", spacetime_blackhole_lw_lapse_function_der },
  { "shiftVectorDer", spacetime_blackhole_lw_shift_vector_der },
  { "spatialMetricTensorDer", spacetime_blackhole_lw_spatial_metric_tensor_der },
  { "excisionRegion", spacetime_blackhole_lw_excision_region },
  { "conformalFactor", spacetime_blackhole_lw_conformal_factor },
  { "bssnConformalFactor", spacetime_blackhole_lw_bssn_conformal_factor },
  { "conformalFactorDer", spacetime_blackhole_lw_conformal_factor_der },
  { "bssnConformalFactorDer", spacetime_blackhole_lw_bssn_conformal_factor_der },
  { "bssnConformalFactorDer2", spacetime_blackhole_lw_bssn_conformal_factor_der2 },
  { 0, 0 }
};

/* ********************** */
/* Neutron Star Spacetime */
/* ********************** */

static int
spacetime_neutronstar_lw_new(lua_State *L)
{
  struct gr_spacetime_lw *neutronstar_lw = gkyl_malloc(sizeof(*neutronstar_lw));

  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double mass_quadrupole = luaL_checknumber(L, 3);
  double spin_octupole = luaL_checknumber(L, 4);
  double mass_hexadecapole = luaL_checknumber(L, 5);
  double pos_x = luaL_checknumber(L, 6);
  double pos_y = luaL_checknumber(L, 7);
  double pos_z = luaL_checknumber(L, 8);

  neutronstar_lw->magic = MOMENT_SPACETIME_DEFAULT;
  neutronstar_lw->spacetime = gkyl_gr_neutronstar_inew( &(struct gkyl_gr_neutronstar_inp) {
      .mass = mass,
      .spin = spin,
      .mass_quadrupole = mass_quadrupole,
      .spin_octupole = spin_octupole,
      .mass_hexadecapole = mass_hexadecapole,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  // Create Lua userdata.
  struct gr_spacetime_lw **l_neutronstar_lw = lua_newuserdata(L, sizeof(struct gr_spacetime_lw*));
  *l_neutronstar_lw = neutronstar_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_SPACETIME_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

static int
spacetime_neutronstar_lw_spatial_metric_tensor(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double mass_quadrupole = luaL_checknumber(L, 3);
  double spin_octupole = luaL_checknumber(L, 4);
  double mass_hexadecapole = luaL_checknumber(L, 5);
  double pos_x = luaL_checknumber(L, 6);
  double pos_y = luaL_checknumber(L, 7);
  double pos_z = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_neutronstar_inew( &(struct gkyl_gr_neutronstar_inp) {
      .mass = mass,
      .spin = spin,
      .mass_quadrupole = mass_quadrupole,
      .spin_octupole = spin_octupole,
      .mass_hexadecapole = mass_hexadecapole,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->spatial_metric_tensor_func(spacetime, t, x, y, z, &spatial_metric);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, spatial_metric[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
  }
  gkyl_free(spatial_metric);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_neutronstar_lw_spatial_metric_det(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double mass_quadrupole = luaL_checknumber(L, 3);
  double spin_octupole = luaL_checknumber(L, 4);
  double mass_hexadecapole = luaL_checknumber(L, 5);
  double pos_x = luaL_checknumber(L, 6);
  double pos_y = luaL_checknumber(L, 7);
  double pos_z = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_neutronstar_inew( &(struct gkyl_gr_neutronstar_inp) {
      .mass = mass,
      .spin = spin,
      .mass_quadrupole = mass_quadrupole,
      .spin_octupole = spin_octupole,
      .mass_hexadecapole = mass_hexadecapole,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double spatial_metric_det;
  spacetime->spatial_metric_det_func(spacetime, t, x, y, z, &spatial_metric_det);

  lua_pushnumber(L, spatial_metric_det);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_neutronstar_lw_lapse_function(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double mass_quadrupole = luaL_checknumber(L, 3);
  double spin_octupole = luaL_checknumber(L, 4);
  double mass_hexadecapole = luaL_checknumber(L, 5);
  double pos_x = luaL_checknumber(L, 6);
  double pos_y = luaL_checknumber(L, 7);
  double pos_z = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_neutronstar_inew( &(struct gkyl_gr_neutronstar_inp) {
      .mass = mass,
      .spin = spin,
      .mass_quadrupole = mass_quadrupole,
      .spin_octupole = spin_octupole,
      .mass_hexadecapole = mass_hexadecapole,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double lapse;
  spacetime->lapse_function_func(spacetime, t, x, y, z, &lapse);

  lua_pushnumber(L, lapse);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_neutronstar_lw_shift_vector(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double mass_quadrupole = luaL_checknumber(L, 3);
  double spin_octupole = luaL_checknumber(L, 4);
  double mass_hexadecapole = luaL_checknumber(L, 5);
  double pos_x = luaL_checknumber(L, 6);
  double pos_y = luaL_checknumber(L, 7);
  double pos_z = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_neutronstar_inew( &(struct gkyl_gr_neutronstar_inp) {
      .mass = mass,
      .spin = spin,
      .mass_quadrupole = mass_quadrupole,
      .spin_octupole = spin_octupole,
      .mass_hexadecapole = mass_hexadecapole,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double *shift = gkyl_malloc(sizeof(double[3]));
  spacetime->shift_vector_func(spacetime, t, x, y, z, &shift);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, shift[i]);
    lua_rawset(L, -3);
  }

  gkyl_free(shift);
  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_neutronstar_lw_extrinsic_curvature_tensor(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double mass_quadrupole = luaL_checknumber(L, 3);
  double spin_octupole = luaL_checknumber(L, 4);
  double mass_hexadecapole = luaL_checknumber(L, 5);
  double pos_x = luaL_checknumber(L, 6);
  double pos_y = luaL_checknumber(L, 7);
  double pos_z = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_neutronstar_inew( &(struct gkyl_gr_neutronstar_inp) {
      .mass = mass,
      .spin = spin,
      .mass_quadrupole = mass_quadrupole,
      .spin_octupole = spin_octupole,
      .mass_hexadecapole = mass_hexadecapole,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double dx = luaL_checknumber(L, 13);
  double dy = luaL_checknumber(L, 14);
  double dz = luaL_checknumber(L, 15);

  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->extrinsic_curvature_tensor_func(spacetime, t, x, y, z, dx, dy, dz, &extrinsic_curvature);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, extrinsic_curvature[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(extrinsic_curvature[i]);
  }
  gkyl_free(extrinsic_curvature);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_neutronstar_lw_lapse_function_der(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double mass_quadrupole = luaL_checknumber(L, 3);
  double spin_octupole = luaL_checknumber(L, 4);
  double mass_hexadecapole = luaL_checknumber(L, 5);
  double pos_x = luaL_checknumber(L, 6);
  double pos_y = luaL_checknumber(L, 7);
  double pos_z = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_neutronstar_inew( &(struct gkyl_gr_neutronstar_inp) {
      .mass = mass,
      .spin = spin,
      .mass_quadrupole = mass_quadrupole,
      .spin_octupole = spin_octupole,
      .mass_hexadecapole = mass_hexadecapole,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double dx = luaL_checknumber(L, 13);
  double dy = luaL_checknumber(L, 14);
  double dz = luaL_checknumber(L, 15);

  double *lapse_der = gkyl_malloc(sizeof(double[3]));
  spacetime->lapse_function_der_func(spacetime, t, x, y, z, dx, dy, dz, &lapse_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, lapse_der[i]);
    lua_rawset(L, -3);
  }

  gkyl_free(lapse_der);
  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_neutronstar_lw_shift_vector_der(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double mass_quadrupole = luaL_checknumber(L, 3);
  double spin_octupole = luaL_checknumber(L, 4);
  double mass_hexadecapole = luaL_checknumber(L, 5);
  double pos_x = luaL_checknumber(L, 6);
  double pos_y = luaL_checknumber(L, 7);
  double pos_z = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_neutronstar_inew( &(struct gkyl_gr_neutronstar_inp) {
      .mass = mass,
      .spin = spin,
      .mass_quadrupole = mass_quadrupole,
      .spin_octupole = spin_octupole,
      .mass_hexadecapole = mass_hexadecapole,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double dx = luaL_checknumber(L, 13);
  double dy = luaL_checknumber(L, 14);
  double dz = luaL_checknumber(L, 15);

  double **shift_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    shift_der[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->shift_vector_der_func(spacetime, t, x, y, z, dx, dy, dz, &shift_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, shift_der[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(shift_der[i]);
  }
  gkyl_free(shift_der);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_neutronstar_lw_spatial_metric_tensor_der(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double mass_quadrupole = luaL_checknumber(L, 3);
  double spin_octupole = luaL_checknumber(L, 4);
  double mass_hexadecapole = luaL_checknumber(L, 5);
  double pos_x = luaL_checknumber(L, 6);
  double pos_y = luaL_checknumber(L, 7);
  double pos_z = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_neutronstar_inew( &(struct gkyl_gr_neutronstar_inp) {
      .mass = mass,
      .spin = spin,
      .mass_quadrupole = mass_quadrupole,
      .spin_octupole = spin_octupole,
      .mass_hexadecapole = mass_hexadecapole,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double dx = luaL_checknumber(L, 13);
  double dy = luaL_checknumber(L, 14);
  double dz = luaL_checknumber(L, 15);

  double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

    for (int j = 0; j < 3; j++) {
      spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
    }
  }

  spacetime->spatial_metric_tensor_der_func(spacetime, t, x, y, z, dx, dy, dz, &spatial_metric_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);

      lua_createtable(L, 3, 0);
      for (int k = 0; k < 3; k++) {
        lua_pushinteger(L, k + 1);
        lua_pushnumber(L, spatial_metric_der[i][j][k]);
        lua_rawset(L, -3);
      }

      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      gkyl_free(spatial_metric_der[i][j]);
    }
    gkyl_free(spatial_metric_der[i]);
  }
  gkyl_free(spatial_metric_der);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_neutronstar_lw_excision_region(lua_State *L)
{
  double mass = luaL_checknumber(L, 1);
  double spin = luaL_checknumber(L, 2);
  double mass_quadrupole = luaL_checknumber(L, 3);
  double spin_octupole = luaL_checknumber(L, 4);
  double mass_hexadecapole = luaL_checknumber(L, 5);
  double pos_x = luaL_checknumber(L, 6);
  double pos_y = luaL_checknumber(L, 7);
  double pos_z = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_neutronstar_inew( &(struct gkyl_gr_neutronstar_inp) {
      .mass = mass,
      .spin = spin,
      .mass_quadrupole = mass_quadrupole,
      .spin_octupole = spin_octupole,
      .mass_hexadecapole = mass_hexadecapole,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  bool in_excision_region;
  spacetime->excision_region_func(spacetime, t, x, y, z, &in_excision_region);

  lua_pushboolean(L, in_excision_region);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

// Spacetime constructor.
static struct luaL_Reg spacetime_neutronstar_ctor[] = {
  { "new", spacetime_neutronstar_lw_new },
  { "spatialMetricTensor", spacetime_neutronstar_lw_spatial_metric_tensor },
  { "spatialMetricDeterminant", spacetime_neutronstar_lw_spatial_metric_det },
  { "lapseFunction", spacetime_neutronstar_lw_lapse_function },
  { "shiftVector", spacetime_neutronstar_lw_shift_vector },
  { "extrinsicCurvatureTensor", spacetime_neutronstar_lw_extrinsic_curvature_tensor },
  { "lapseFunctionDer", spacetime_neutronstar_lw_lapse_function_der },
  { "shiftVectorDer", spacetime_neutronstar_lw_shift_vector_der },
  { "spatialMetricTensorDer", spacetime_neutronstar_lw_spatial_metric_tensor_der },
  { "excisionRegion", spacetime_neutronstar_lw_excision_region },
  { 0, 0 }
};

/* ********************************************* */
/* Brill-Lindquist (Binary Black Hole) Spacetime */
/* ********************************************* */

static int
spacetime_brill_lindquist_lw_new(lua_State *L)
{
  struct gr_spacetime_lw *brill_lindquist_lw = gkyl_malloc(sizeof(*brill_lindquist_lw));

  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  brill_lindquist_lw->magic = MOMENT_SPACETIME_DEFAULT;
  brill_lindquist_lw->spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  // Create Lua userdata.
  struct gr_spacetime_lw **l_brill_lindquist_lw = lua_newuserdata(L, sizeof(struct gr_spacetime_lw*));
  *l_brill_lindquist_lw = brill_lindquist_lw;

  // Set metatable.
  luaL_getmetatable(L, MOMENT_SPACETIME_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

static int
spacetime_brill_lindquist_lw_spatial_metric_tensor(lua_State *L)
{
  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->spatial_metric_tensor_func(spacetime, t, x, y, z, &spatial_metric);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, spatial_metric[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
  }
  gkyl_free(spatial_metric);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_brill_lindquist_lw_inv_spatial_metric_tensor(lua_State *L)
{
  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->spatial_inv_metric_tensor_func(spacetime, t, x, y, z, &inv_spatial_metric);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, inv_spatial_metric[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(inv_spatial_metric);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_brill_lindquist_lw_spatial_metric_det(lua_State *L)
{
  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double spatial_metric_det;
  spacetime->spatial_metric_det_func(spacetime, t, x, y, z, &spatial_metric_det);

  lua_pushnumber(L, spatial_metric_det);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_brill_lindquist_lw_lapse_function(lua_State *L)
{
  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double lapse;
  spacetime->lapse_function_func(spacetime, t, x, y, z, &lapse);

  lua_pushnumber(L, lapse);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_brill_lindquist_lw_shift_vector(lua_State *L)
{
  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double *shift = gkyl_malloc(sizeof(double[3]));
  spacetime->shift_vector_func(spacetime, t, x, y, z, &shift);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, shift[i]);
    lua_rawset(L, -3);
  }

  gkyl_free(shift);
  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_brill_lindquist_lw_extrinsic_curvature_tensor(lua_State *L)
{
  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double dx = luaL_checknumber(L, 13);
  double dy = luaL_checknumber(L, 14);
  double dz = luaL_checknumber(L, 15);

  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->extrinsic_curvature_tensor_func(spacetime, t, x, y, z, dx, dy, dz, &extrinsic_curvature);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, extrinsic_curvature[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(extrinsic_curvature[i]);
  }
  gkyl_free(extrinsic_curvature);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_brill_lindquist_lw_lapse_function_der(lua_State *L)
{
  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double dx = luaL_checknumber(L, 13);
  double dy = luaL_checknumber(L, 14);
  double dz = luaL_checknumber(L, 15);

  double *lapse_der = gkyl_malloc(sizeof(double[3]));
  spacetime->lapse_function_der_func(spacetime, t, x, y, z, dx, dy, dz, &lapse_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, lapse_der[i]);
    lua_rawset(L, -3);
  }

  gkyl_free(lapse_der);
  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_brill_lindquist_lw_shift_vector_der(lua_State *L)
{
  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double dx = luaL_checknumber(L, 13);
  double dy = luaL_checknumber(L, 14);
  double dz = luaL_checknumber(L, 15);

  double **shift_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    shift_der[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->shift_vector_der_func(spacetime, t, x, y, z, dx, dy, dz, &shift_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, shift_der[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(shift_der[i]);
  }
  gkyl_free(shift_der);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_brill_lindquist_lw_spatial_metric_tensor_der(lua_State *L)
{
  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double dx = luaL_checknumber(L, 13);
  double dy = luaL_checknumber(L, 14);
  double dz = luaL_checknumber(L, 15);

  double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

    for (int j = 0; j < 3; j++) {
      spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
    }
  }

  spacetime->spatial_metric_tensor_der_func(spacetime, t, x, y, z, dx, dy, dz, &spatial_metric_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);

      lua_createtable(L, 3, 0);
      for (int k = 0; k < 3; k++) {
        lua_pushinteger(L, k + 1);
        lua_pushnumber(L, spatial_metric_der[i][j][k]);
        lua_rawset(L, -3);
      }

      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      gkyl_free(spatial_metric_der[i][j]);
    }
    gkyl_free(spatial_metric_der[i]);
  }
  gkyl_free(spatial_metric_der);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_brill_lindquist_lw_excision_region(lua_State *L)
{
  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  bool in_excision_region;
  spacetime->excision_region_func(spacetime, t, x, y, z, &in_excision_region);

  lua_pushboolean(L, in_excision_region);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_brill_lindquist_lw_conformal_factor(lua_State *L)
{
  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double conformal_factor;
  spacetime->conformal_factor_func(spacetime, t, x, y, z, &conformal_factor);

  lua_pushnumber(L, conformal_factor);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_brill_lindquist_lw_bssn_conformal_factor(lua_State *L)
{
  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double bssn_conformal_factor;
  spacetime->bssn_conformal_factor_func(spacetime, t, x, y, z, &bssn_conformal_factor);

  lua_pushnumber(L, bssn_conformal_factor);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_brill_lindquist_lw_conformal_factor_der(lua_State *L)
{
  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double dx = luaL_checknumber(L, 13);
  double dy = luaL_checknumber(L, 14);
  double dz = luaL_checknumber(L, 15);

  double *conformal_factor_der = gkyl_malloc(sizeof(double[3]));
  spacetime->conformal_factor_der_func(spacetime, t, x, y, z, dx, dy, dz, &conformal_factor_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, conformal_factor_der[i]);
    lua_rawset(L, -3);
  }

  gkyl_free(conformal_factor_der);
  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_brill_lindquist_lw_bssn_conformal_factor_der(lua_State *L)
{
  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double dx = luaL_checknumber(L, 13);
  double dy = luaL_checknumber(L, 14);
  double dz = luaL_checknumber(L, 15);

  double *bssn_conformal_factor_der = gkyl_malloc(sizeof(double[3]));
  spacetime->bssn_conformal_factor_der_func(spacetime, t, x, y, z, dx, dy, dz, &bssn_conformal_factor_der);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, bssn_conformal_factor_der[i]);
    lua_rawset(L, -3);
  }

  gkyl_free(bssn_conformal_factor_der);
  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

static int
spacetime_brill_lindquist_lw_bssn_conformal_factor_der2(lua_State *L)
{
  double mass1 = luaL_checknumber(L, 1);
  double mass2 = luaL_checknumber(L, 2);
  double pos_x1 = luaL_checknumber(L, 3);
  double pos_y1 = luaL_checknumber(L, 4);
  double pos_z1 = luaL_checknumber(L, 5);
  double pos_x2 = luaL_checknumber(L, 6);
  double pos_y2 = luaL_checknumber(L, 7);
  double pos_z2 = luaL_checknumber(L, 8);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_inew( &(struct gkyl_gr_brill_lindquist_inp) {
      .mass1 = mass1,
      .mass2 = mass2,
      .pos_x1 = pos_x1,
      .pos_y1 = pos_y1,
      .pos_z1 = pos_z1,
      .pos_x2 = pos_x2,
      .pos_y2 = pos_y2,
      .pos_z2 = pos_z2,
      .use_gpu = false
    }
  );

  double t = luaL_checknumber(L, 9);
  double x = luaL_checknumber(L, 10);
  double y = luaL_checknumber(L, 11);
  double z = luaL_checknumber(L, 12);

  double dx = luaL_checknumber(L, 13);
  double dy = luaL_checknumber(L, 14);
  double dz = luaL_checknumber(L, 15);

  double **bssn_conformal_factor_der2 = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    bssn_conformal_factor_der2[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->bssn_conformal_factor_der2_func(spacetime, t, x, y, z, dx, dy, dz, &bssn_conformal_factor_der2);

  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);

    lua_createtable(L, 3, 0);
    for (int j = 0; j < 3; j++) {
      lua_pushinteger(L, j + 1);
      lua_pushnumber(L, bssn_conformal_factor_der2[i][j]);
      lua_rawset(L, -3);
    }

    lua_rawset(L, -3);
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(bssn_conformal_factor_der2[i]);
  }
  gkyl_free(bssn_conformal_factor_der2);

  gkyl_gr_spacetime_release(spacetime);

  return 1;
}

// Spacetime constructor.
static struct luaL_Reg spacetime_brill_lindquist_ctor[] = {
  { "new", spacetime_brill_lindquist_lw_new },
  { "spatialMetricTensor", spacetime_brill_lindquist_lw_spatial_metric_tensor },
  { "invSpatialMetricTensor", spacetime_brill_lindquist_lw_inv_spatial_metric_tensor },
  { "spatialMetricDeterminant", spacetime_brill_lindquist_lw_spatial_metric_det },
  { "lapseFunction", spacetime_brill_lindquist_lw_lapse_function },
  { "shiftVector", spacetime_brill_lindquist_lw_shift_vector },
  { "extrinsicCurvatureTensor", spacetime_brill_lindquist_lw_extrinsic_curvature_tensor },
  { "lapseFunctionDer", spacetime_brill_lindquist_lw_lapse_function_der },
  { "shiftVectorDer", spacetime_brill_lindquist_lw_shift_vector_der },
  { "spatialMetricTensorDer", spacetime_brill_lindquist_lw_spatial_metric_tensor_der },
  { "excisionRegion", spacetime_brill_lindquist_lw_excision_region },
  { "conformalFactor", spacetime_brill_lindquist_lw_conformal_factor },
  { "bssnConformalFactor", spacetime_brill_lindquist_lw_bssn_conformal_factor },
  { "conformalFactorDer", spacetime_brill_lindquist_lw_conformal_factor_der },
  { "bssnConformalFactorDer", spacetime_brill_lindquist_lw_bssn_conformal_factor_der },
  { "bssnConformalFactorDer2", spacetime_brill_lindquist_lw_bssn_conformal_factor_der2 },
  { 0, 0 }
};

// Register and load all GR spacetime objects.
static void
spacetime_openlibs(lua_State *L)
{
  luaL_newmetatable(L, MOMENT_SPACETIME_METATABLE_NM);

  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, gr_spacetime_lw_gc);
  lua_settable(L, -3);

  luaL_register(L, "G0.Moments.Spacetime.Minkowski", spacetime_minkowski_ctor);
  luaL_register(L, "G0.Moments.Spacetime.BlackHole", spacetime_blackhole_ctor);
  luaL_register(L, "G0.Moments.Spacetime.NeutronStar", spacetime_neutronstar_ctor);
  luaL_register(L, "G0.Moments.Spacetime.BrillLindquist", spacetime_brill_lindquist_ctor);
}

/* *************** */
/* Species methods */
/* *************** */

// Metatable name for species input struct.
#define MOMENT_SPECIES_METATABLE_NM "GkeyllZero.App.Moments.Species"

// Lua userdata object for constructing species input.
struct moment_species_lw {
  int magic; // This must be first element in the struct.
  
  struct gkyl_moment_species mom_species; // Input struct to construct species.
  bool evolve; // Is this species evolved?

  struct lua_func_ctx init_ctx; // Lua registry reference to initialization function.
  
  bool has_applied_acceleration_func; // Is there an applied acceleration initialization function?
  struct lua_func_ctx applied_acceleration_func_ref; // Lua registry reference to applied acceleration initialization function.
  bool evolve_applied_acceleration; // Is the applied acceleration evolved?

  bool has_nT_source; // Is there a temperature source function?
  struct lua_func_ctx nT_source_func_ctx; // Lua registry reference to temperature source function.
};

static int
moment_species_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_moment_species mom_species = { };

  mom_species.charge = glua_tbl_get_number(L, "charge", 0.0);
  mom_species.mass = glua_tbl_get_number(L, "mass", 1.0);

  bool has_eqn = false;
  with_lua_tbl_key(L, "equation") {
    mom_species.equation = wv_eqn_get(L);
    has_eqn = true;
  }

  if (has_eqn) {
    // We need to store a reference to equation object in a global
    // table as otherwise it gets garbage-collected while the
    // simulation is being set up.
    lua_getfield(L, -1, "equation");
    int eq_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    
    lua_getglobal(L, "__moment_ref_table");
    int eqtbl_len = glua_objlen(L);
    lua_pushinteger(L, eqtbl_len+1);    
    lua_rawgeti(L, LUA_REGISTRYINDEX, eq_ref);
    lua_rawset(L, -3);

    lua_pop(L, 1);
  }
  else {
    return luaL_error(L, "Species \"equation\" not specfied or incorrect type!");
  }

  mom_species.limiter = glua_tbl_get_integer(L, "limiter", GKYL_MONOTONIZED_CENTERED);

  const char *split_str = glua_tbl_get_string(L, "splitType", "qwave");
  mom_species.split_type = gkyl_search_str_int_pair_by_str(wave_split_type, split_str, GKYL_WAVE_QWAVE);

  bool evolve = glua_tbl_get_bool(L, "evolve", true);
  mom_species.is_static = !evolve; 
  mom_species.force_low_order_flux = glua_tbl_get_bool(L, "forceLowOrderFlux", false);

  int init_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "init")) {
    init_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }
  else {
    return luaL_error(L, "Species must have an \"init\" function for initial conditions!");
  }

  with_lua_tbl_tbl(L, "bcx") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      mom_species.bcx[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcy") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      mom_species.bcy[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcz") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      mom_species.bcz[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  bool has_applied_acceleration_func = false;
  int applied_acceleration_func_ref = LUA_NOREF;
  bool evolve_applied_acceleration = false;

  if (glua_tbl_get_func(L, "appliedAcceleration")) {
    applied_acceleration_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_applied_acceleration_func = true;

    evolve_applied_acceleration = glua_tbl_get_bool(L, "evolveAppliedAcceleration", false);
  }

  bool has_nT_source = false;
  int nT_source_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "nTSource")) {
    has_nT_source = true;
    nT_source_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }
  mom_species.nT_source_set_only_once = glua_tbl_get_bool(L, "nTSourceSetOnlyOnce", false);

  mom_species.has_reactivity = glua_tbl_get_bool(L, "hasReactivity", false);
  if (mom_species.has_reactivity) {
    mom_species.reactivity_gas_gamma = glua_tbl_get_number(L, "reactivityGasGamma", 1.4);
    mom_species.reactivity_specific_heat_capacity = glua_tbl_get_number(L, "reactivitySpecificHeatCapacity", 2.5);
    mom_species.reactivity_energy_of_formation = glua_tbl_get_number(L, "reactivityEnergyOfFormation", 1.0);
    mom_species.reactivity_ignition_temperature = glua_tbl_get_number(L, "reactivityIgnitionTemperature", 0.25);
    mom_species.reactivity_reaction_rate = glua_tbl_get_number(L, "reactivityReactionRate", 250.0);
  }

  mom_species.has_volume_sources = glua_tbl_get_bool(L, "hasVolumeSources", false);
  if (mom_species.has_volume_sources) {
    mom_species.volume_gas_gamma = glua_tbl_get_number(L, "volumeGasGamma", 5.0 / 3.0);
    mom_species.volume_U0 = glua_tbl_get_number(L, "volumeU0", 1.0);
    mom_species.volume_R0 = glua_tbl_get_number(L, "volumeR0", 1.0);
  }

  mom_species.has_einstein_medium = glua_tbl_get_bool(L, "hasEinsteinMedium", false);
  if (mom_species.has_einstein_medium) {
    mom_species.medium_gas_gamma = glua_tbl_get_number(L, "mediumGasGamma", 4.0 / 3.0);
    mom_species.medium_kappa = glua_tbl_get_number(L, "mediumKappa", 8.0 * M_PI);
  }

  mom_species.has_vacuum_einstein = glua_tbl_get_bool(L, "hasVacuumEinstein", false);
  if (mom_species.has_vacuum_einstein) {
    mom_species.vacuum_einstein_excision_threshold = glua_tbl_get_number(L, "vacuumEinsteinExcisionThreshold", 0.3);
    mom_species.vacuum_einstein_spacetime_slicing = glua_tbl_get_integer(L, "vacuumEinsteinSpacetimeSlicing", GKYL_GEODESIC_SLICING);
    mom_species.vacuum_einstein_spacetime_evolution = glua_tbl_get_integer(L, "vacuumEinsteinSpacetimeEvolution", GKYL_EINSTEIN_EVOLUTION);
  }

  mom_species.has_vacuum_einstein_conformal = glua_tbl_get_bool(L, "hasVacuumEinsteinConformal", false);
  if (mom_species.has_vacuum_einstein_conformal) {
    mom_species.vacuum_einstein_conformal_excision_threshold = glua_tbl_get_number(L, "vacuumEinsteinConformalExcisionThreshold", 0.3);
    mom_species.vacuum_einstein_conformal_spacetime_slicing = glua_tbl_get_integer(L, "vacuumEinsteinConformalSpacetimeSlicing", GKYL_GEODESIC_SLICING);
    mom_species.vacuum_einstein_conformal_spacetime_evolution = glua_tbl_get_integer(L, "vacuumEinsteinConformalSpacetimeEvolution", GKYL_EINSTEIN_EVOLUTION);
  }

  mom_species.has_gr_ultra_rel = glua_tbl_get_bool(L, "hasGRUltraRel", false);
  if (mom_species.has_gr_ultra_rel) {
    mom_species.gr_ultra_rel_gas_gamma = glua_tbl_get_number(L, "GRUltraRelGasGamma", 5.0 / 3.0);
  }

  mom_species.has_gr_euler = glua_tbl_get_bool(L, "hasGREuler", false);
  if (mom_species.has_gr_euler) {
    mom_species.gr_euler_gas_gamma = glua_tbl_get_number(L, "GREulerGasGamma", 5.0 / 3.0);
  }

  mom_species.has_gr_twofluid = glua_tbl_get_bool(L, "hasGRTwoFluid", false);
  if (mom_species.has_gr_twofluid) {
    mom_species.gr_twofluid_mass_elc = glua_tbl_get_number(L, "GRTwoFluidMassElc", 1.0 / 1836.2);
    mom_species.gr_twofluid_mass_ion = glua_tbl_get_number(L, "GRTwoFluidMassIon", 1.0);
    mom_species.gr_twofluid_charge_elc = glua_tbl_get_number(L, "GRTwoFluidChargeElc", -1.0);
    mom_species.gr_twofluid_charge_ion = glua_tbl_get_number(L, "GRTwoFluidChargeIon", 1.0);
    mom_species.gr_twofluid_gas_gamma_elc = glua_tbl_get_number(L, "GRTwoFluidGasGammaElc", 5.0 / 3.0);
    mom_species.gr_twofluid_gas_gamma_ion = glua_tbl_get_number(L, "GRTwoFluidGasGammaIon", 5.0 / 3.0);
    mom_species.gr_twofluid_e_fact = glua_tbl_get_number(L, "GRTwoFluidEFact", 0.0);
  }

  mom_species.has_gr_mhd = glua_tbl_get_bool(L, "hasGRMHD", false);
  if (mom_species.has_gr_mhd) {
    mom_species.gr_mhd_gas_gamma = glua_tbl_get_number(L, "GRMHDGasGamma", 5.0 / 3.0);
  }
  
  mom_species.type_brag = glua_tbl_get_integer(L, "braginskiiType", 0);

  struct moment_species_lw *moms_lw = lua_newuserdata(L, sizeof(*moms_lw));
  moms_lw->magic = MOMENT_SPECIES_DEFAULT;
  moms_lw->evolve = evolve;
  moms_lw->mom_species = mom_species;

  moms_lw->init_ctx = (struct lua_func_ctx) {
    .func_ref = init_ref,
    .ndim = 0, // This will be set later.
    .nret = mom_species.equation->num_equations,
    .L = L,
  };

  moms_lw->has_applied_acceleration_func = has_applied_acceleration_func;
  moms_lw->applied_acceleration_func_ref = (struct lua_func_ctx) {
    .func_ref = applied_acceleration_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 3,
    .L = L,
  };
  moms_lw->evolve_applied_acceleration = evolve_applied_acceleration;

  moms_lw->has_nT_source = has_nT_source;
  moms_lw->nT_source_func_ctx = (struct lua_func_ctx) {
    .func_ref = nT_source_ref,
    .ndim = 0, // This will be set later.
    .nret = GKYL_MOM_APP_NUM_NT_SOURCE,
    .L = L,
  };    
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_SPECIES_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor.
static struct luaL_Reg mom_species_ctor[] = {
  { "new", moment_species_lw_new },
  { 0, 0 }
};

/* ************* */
/* Field methods */
/* ************* */

// Metatable name for field input struct.
#define MOMENT_FIELD_METATABLE_NM "GkeyllZero.App.Moments.Field"

// Lua userdata object for constructing field input.
struct moment_field_lw {
  int magic; // This must be first element in the struct.

  bool evolve; // Is this field evolved?  
  struct gkyl_moment_field mom_field; // Input struct to construct field.
  
  struct lua_func_ctx init_ctx; // Lua registry reference to initialization function.

  bool has_external_field_func; // Is there an external field initialization function?
  struct lua_func_ctx external_field_func_ref; // Lua registry reference to external field initialization function.
  bool evolve_external_field; // Is the external field evolved?
  double external_field_ramp_time; // Linear ramp for turning on external field without re-projecting. 

  bool has_applied_current_func; // Is there an applied current initialization function?
  struct lua_func_ctx applied_current_func_ref; // Lua registry reference to applied current initialization function.
  bool evolve_applied_current; // Is the applied current evolved?
  double applied_current_ramp_time; // Linear ramp for turning on applied current without re-projecting. 
};

static int
moment_field_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_moment_field mom_field = { };

  mom_field.epsilon0 = glua_tbl_get_number(L, "epsilon0", 1.0);
  mom_field.mu0 = glua_tbl_get_number(L, "mu0", 1.0);
  mom_field.elc_error_speed_fact = glua_tbl_get_number(L, "elcErrorSpeedFactor", 0.0);
  mom_field.mag_error_speed_fact = glua_tbl_get_number(L, "mgnErrorSpeedFactor", 0.0);

  mom_field.limiter = glua_tbl_get_integer(L, "limiter", GKYL_MONOTONIZED_CENTERED);
  
  bool evolve = glua_tbl_get_bool(L, "evolve", true);
  mom_field.is_static = !evolve; 
  mom_field.use_explicit_em_coupling = glua_tbl_get_bool(L, "useExplicitEmCoupling", false);

  int init_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "init")) {
    init_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }
  else {
    return luaL_error(L, "Field must have an \"init\" function for initial conditions!");
  }
  
  with_lua_tbl_tbl(L, "bcx") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      mom_field.bcx[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcy") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      mom_field.bcy[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcz") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      mom_field.bcz[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  bool has_external_field_func = false;
  int external_field_func_ref = LUA_NOREF;
  bool evolve_external_field = false;
  double external_field_ramp_time = 0.0; 

  if (glua_tbl_get_func(L, "externalFieldInit")) {
    external_field_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_external_field_func = true;

    evolve_external_field = glua_tbl_get_bool(L, "evolveExternalField", false);
    external_field_ramp_time = glua_tbl_get_number(L, "externalFieldRampTime", 0.0);
  }

  bool has_applied_current_func = false;
  int applied_current_func_ref = LUA_NOREF;
  bool evolve_applied_current = false;
  double applied_current_ramp_time = 0.0; 

  if (glua_tbl_get_func(L, "appliedCurrent")) {
    applied_current_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_applied_current_func = true;

    evolve_applied_current = glua_tbl_get_bool(L, "evolveAppliedCurrent", false);
    applied_current_ramp_time = glua_tbl_get_number(L, "appliedCurrentRampTime", 0.0);
  }

  struct moment_field_lw *momf_lw = lua_newuserdata(L, sizeof(*momf_lw));

  momf_lw->magic = MOMENT_FIELD_DEFAULT;
  momf_lw->evolve = evolve;
  momf_lw->mom_field = mom_field;

  momf_lw->init_ctx = (struct lua_func_ctx) {
    .func_ref = init_ref,
    .ndim = 0, // This will be set later.
    .nret = 6,
    .L = L,
  };  

  momf_lw->has_external_field_func = has_external_field_func;
  momf_lw->external_field_func_ref = (struct lua_func_ctx) {
    .func_ref = external_field_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 6,
    .L = L,
  };
  momf_lw->evolve_external_field = evolve_external_field;
  momf_lw->external_field_ramp_time = external_field_ramp_time;

  momf_lw->has_applied_current_func = has_applied_current_func;
  momf_lw->applied_current_func_ref = (struct lua_func_ctx) {
    .func_ref = applied_current_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 3,
    .L = L,
  };
  momf_lw->evolve_applied_current = evolve_applied_current;
  momf_lw->applied_current_ramp_time = applied_current_ramp_time;
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_FIELD_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor.
static struct luaL_Reg mom_field_ctor[] = {
  { "new", moment_field_lw_new },
  { 0, 0 }
};

/* *********** */
/* App methods */
/* *********** */

// Metatable name for top-level Moment App.
#define MOMENT_APP_METATABLE_NM "GkeyllZero.App.Moment"

// Lua userdata object for holding Moment app and simulation parameters.
struct moment_app_lw {
  gkyl_moment_app *app; // Moment app object.
  
  struct lua_func_ctx mapc2p_ctx; // Function context for mapc2p.

  struct lua_func_ctx species_init_ctx[GKYL_MAX_SPECIES]; // Function context for species initial conditions.
  struct lua_func_ctx applied_acceleration_func_ctx[GKYL_MAX_SPECIES]; // Function context for applied acceleration.
  struct lua_func_ctx species_nT_source_func_ctx[GKYL_MAX_SPECIES]; // Function context for temperature sources.

  struct lua_func_ctx field_init_ctx; // Function context for field initial conditions.
  struct lua_func_ctx external_field_func_ctx; // Function context for external field.
  struct lua_func_ctx applied_current_func_ctx; // Function context for applied current.
  
  double t_start, t_end; // Start and end times of simulation.
  int num_frames; // Number of data frames to write.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

// Find index of species in table given its name.
static int
app_find_species(const gkyl_moment_app *app, const char *nm)
{
  for (int i = 0; i < app->num_species; i++) {
    if (strcmp(app->species[i].name, nm) == 0) {
      return i;
    }
  }
  
  return -1;
}

// Gets all species objects from the App table, which must on top of
// the stack. The number of species is returned and the appropriate
// pointers set in the species pointer array.
static int
get_species_inp(lua_State *L, int cdim, struct moment_species_lw *species[GKYL_MAX_SPECIES])
{
  enum { TKEY = -2, TVAL = -1};
  
  int curr = 0;
  lua_pushnil(L); // Initial key is nil.
  while (lua_next(L, TKEY) != 0) {
    // Key at TKEY and value at TVAL.
    if (lua_type(L, TVAL) == LUA_TUSERDATA) {
      struct moment_species_lw *vms = lua_touserdata(L, TVAL);
      if (vms->magic == MOMENT_SPECIES_DEFAULT) {
        
        vms->init_ctx.ndim = cdim;
        vms->applied_acceleration_func_ref.ndim = cdim;
        vms->nT_source_func_ctx.ndim = cdim;
        
        if (lua_type(L,TKEY) == LUA_TSTRING) {
          const char *key = lua_tolstring(L, TKEY, 0);
          strcpy(vms->mom_species.name, key);
        }
        species[curr++] = vms;
      }
    }
    lua_pop(L, 1);
  }
  return curr;
}

// comparison method to sort species array by species name
static int
species_compare_func(const void *a, const void *b)
{
  const struct moment_species_lw *const *spa = a;
  const struct moment_species_lw *const *spb = b;
  return strcmp((*spa)->mom_species.name, (*spb)->mom_species.name);
}

// Create top-level App object.
static int
mom_app_new(lua_State *L)
{
  struct moment_app_lw *app_lw = gkyl_malloc(sizeof(*app_lw));

  // The output prefix to use is stored in the global
  // GKYL_OUT_PREFIX. If this is not found then "g0-moment" is used.
  const char *sim_name = "g0-moment";
  with_lua_global(L, "GKYL_OUT_PREFIX") {
    if (lua_isstring(L, -1)) {
      sim_name = lua_tostring(L, -1);
    }
  }
  
  // Initialize app using table inputs (table is on top of stack).

  app_lw->t_start = glua_tbl_get_number(L, "tStart", 0.0);
  app_lw->t_end = glua_tbl_get_number(L, "tEnd", 1.0);
  app_lw->num_frames = glua_tbl_get_integer(L, "nFrame", 1);
  app_lw->field_energy_calcs = glua_tbl_get_integer(L, "fieldEnergyCalcs", INT_MAX);
  app_lw->integrated_mom_calcs = glua_tbl_get_integer(L, "integratedMomentCalcs", INT_MAX);
  app_lw->dt_failure_tol = glua_tbl_get_number(L, "dtFailureTol", 1.0e-4);
  app_lw->num_failures_max = glua_tbl_get_integer(L, "numFailuresMax", 20);

  struct gkyl_moment mom = { }; // Input table for app.

  strcpy(mom.name, sim_name);
  
  int cdim = 0;
  with_lua_tbl_tbl(L, "cells") {
    mom.ndim = cdim = glua_objlen(L);

    for (int d = 0; d < cdim; d++) {
      mom.cells[d] = glua_tbl_iget_integer(L, d + 1, 0);
    }
  }

  int cuts[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; d++) {
    cuts[d] = 1;
  }
  
  with_lua_tbl_tbl(L, "decompCuts") {
    int ncuts = glua_objlen(L);

    for (int d = 0; d < ncuts; d++) {
      cuts[d] = glua_tbl_iget_integer(L, d + 1, 0);
    }
  }  

  with_lua_tbl_tbl(L, "lower") {
    for (int d = 0; d < cdim; d++) {
      mom.lower[d] = glua_tbl_iget_number(L, d + 1, 0.0);
    }
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d = 0; d < cdim; d++) {
      mom.upper[d] = glua_tbl_iget_number(L, d + 1, 0.0);
    }
  }

  mom.cfl_frac = glua_tbl_get_number(L, "cflFrac", 0.95);

  mom.scheme_type = glua_tbl_get_integer(L, "schemeType", 0);

  mom.num_periodic_dir = 0;
  if (glua_tbl_has_key(L, "periodicDirs")) {
    with_lua_tbl_tbl(L, "periodicDirs") {
      mom.num_periodic_dir = glua_objlen(L);

      for (int d = 0; d < mom.num_periodic_dir; d++) {
        // Indices are off by 1 between Lua and C.
        mom.periodic_dirs[d] = glua_tbl_iget_integer(L, d + 1, 0) - 1;
      }
    }
  }

  mom.has_braginskii = glua_tbl_get_bool(L, "hasBraginskii", false);
  mom.coll_fac = glua_tbl_get_number(L, "collisionFactor", 0.0);

  // mapc2p function.
  mom.c2p_ctx = 0;
  mom.mapc2p = 0;
  bool has_mapc2p = false;
  int mapc2p_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "mapc2p")) {
    has_mapc2p = true;
    mapc2p_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }

  if (has_mapc2p) {
    app_lw->mapc2p_ctx = (struct lua_func_ctx) {
      .func_ref = mapc2p_ref,
      .ndim = cdim,
      .nret = cdim,
      .L = L,
    };
    mom.mapc2p = gkyl_lw_eval_cb;
    mom.c2p_ctx = &app_lw->mapc2p_ctx;
  }

  struct moment_species_lw *species[GKYL_MAX_SPECIES];

  // Set all species input.
  mom.num_species = get_species_inp(L, cdim, species);

  // need to sort the species[] array by name of the species before
  // proceeding as there is no way to ensure that all cores loop over
  // Lua tables in the same order
  qsort(species, mom.num_species, sizeof(struct moment_species_lw *), species_compare_func);
  
  for (int s = 0; s < mom.num_species; s++) {
    mom.species[s] = species[s]->mom_species;
    
    app_lw->species_init_ctx[s] = species[s]->init_ctx;
    mom.species[s].init = gkyl_lw_eval_cb;
    mom.species[s].ctx = &app_lw->species_init_ctx[s];

    app_lw->applied_acceleration_func_ctx[s] = species[s]->applied_acceleration_func_ref;
    if (species[s]->has_applied_acceleration_func) {
      mom.species[s].app_accel = gkyl_lw_eval_cb;
      mom.species[s].app_accel_ctx = &app_lw->applied_acceleration_func_ctx[s];
      mom.species[s].app_accel_evolve = species[s]->evolve_applied_acceleration;
    }

    if (species[s]->has_nT_source) {
      app_lw->species_nT_source_func_ctx[s] = species[s]->nT_source_func_ctx;
      mom.species[s].nT_source_func = gkyl_lw_eval_cb;
      mom.species[s].nT_source_ctx = &app_lw->species_nT_source_func_ctx[s];
    }
  }

  mom.has_collision = glua_tbl_get_bool(L, "hasCollision", false);
  with_lua_tbl_tbl(L, "nuBase") {
    for (int s = 0; s < mom.num_species; s++) {
      if (glua_tbl_iget_tbl(L, s + 1)) {
        for (int s2 = 0; s2 < mom.num_species; s2++) {
          mom.nu_base[s][s2] = glua_tbl_iget_number(L, s2 + 1, 0.0);
        }

        lua_pop(L, 1);
      }
    }
  }

  // Set field input.
  with_lua_tbl_key(L, "field") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct moment_field_lw *momf = lua_touserdata(L, -1);
      
      if (momf->magic == MOMENT_FIELD_DEFAULT) {
        momf->init_ctx.ndim = cdim;
        
        mom.field = momf->mom_field;

        app_lw->field_init_ctx = momf->init_ctx;
        mom.field.init = gkyl_lw_eval_cb;
        mom.field.ctx = &app_lw->field_init_ctx;

        if (momf->has_external_field_func) {
          momf->external_field_func_ref.ndim = cdim;

          app_lw->external_field_func_ctx = momf->external_field_func_ref;
          mom.field.ext_em = gkyl_lw_eval_cb;
          mom.field.ext_em_ctx = &app_lw->external_field_func_ctx;

          mom.field.ext_em_evolve = momf->evolve_external_field;
          mom.field.t_ramp_E = momf->external_field_ramp_time;
        }

        if (momf->has_applied_current_func) {
          momf->applied_current_func_ref.ndim = cdim;

          app_lw->applied_current_func_ctx = momf->applied_current_func_ref;
          mom.field.app_current = gkyl_lw_eval_cb;
          mom.field.app_current_ctx = &app_lw->applied_current_func_ctx;

          mom.field.app_current_evolve = momf->evolve_applied_current;
          mom.field.t_ramp_curr = momf->applied_current_ramp_time;
        }
      }
    }
  }

  // Create parallelism.
  struct gkyl_comm *comm = 0;
  bool has_mpi = false;

  for (int d = 0; d < cdim; d++) {
    mom.parallelism.cuts[d] = cuts[d]; 
  }

#ifdef GKYL_HAVE_MPI
  with_lua_global(L, "GKYL_MPI_COMM") {
    if (lua_islightuserdata(L, -1)) {
      has_mpi = true;
      struct { MPI_Comm comm;} *lw_mpi_comm_world
        = lua_touserdata(L, -1);
      MPI_Comm mpi_comm = lw_mpi_comm_world->comm;
      comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
          .mpi_comm = mpi_comm,
          .sync_corners = true
        }
      );

    }
  }
#endif

  if (!has_mpi) {
    // If there is no proper MPI_Comm specifed, then assume we are a
    // serial sim.
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .sync_corners = true
      }
    );
  }
  mom.parallelism.comm = comm;

  int rank;
  gkyl_comm_get_rank(comm, &rank);

  int comm_sz;
  gkyl_comm_get_size(comm, &comm_sz);

  int tot_cuts = 1; for (int d=0; d<cdim; ++d) tot_cuts *= cuts[d];

  if (tot_cuts != comm_sz) {
    if (0 == rank)
      fprintf(stderr, "tot_cuts = %d (%d)\n", tot_cuts, comm_sz);
    luaL_error(L, "Number of ranks and cuts do not match!");
  }
  
  app_lw->app = gkyl_moment_app_new(&mom); // Create the Moment app.

  gkyl_comm_release(comm);
  
  // Create Lua userdata.
  struct moment_app_lw **l_app_lw = lua_newuserdata(L, sizeof(struct moment_app_lw*));
  *l_app_lw = app_lw; // Point it to the Lua app pointer.

  // Set metatable.
  luaL_getmetatable(L, MOMENT_APP_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Return number of species () -> int.
static int
mom_app_num_species(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  lua_pushnumber(L, app_lw->app->num_species);

  return 1;
}

// Return number of species (int) -> string.
static int
mom_app_species_name(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  int idx = luaL_checkinteger(L, 2) - 1;
  lua_pushstring(L, app_lw->app->species[idx].name);

  return 1;
}

// Compute maximum time-step () -> double.
static int
mom_app_max_dt(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double maxdt = gkyl_moment_app_max_dt(app_lw->app);

  lua_pushnumber(L, maxdt);  
  return 1;
}

// Apply initial conditions. (time) -> bool.
static int
mom_app_apply_ic(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->t_start);
  gkyl_moment_app_apply_ic(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to field. (time) -> bool.
static int
mom_app_apply_ic_field(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->t_start);
  gkyl_moment_app_apply_ic_field(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to species. (name, time) -> bool.
static int
mom_app_apply_ic_species(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  const char *sp_name = luaL_checkstring(L, 2);
  int sidx = app_find_species(app_lw->app, sp_name);
  if (sidx < 0) {
    return luaL_error(L, "Incorrect species name '%s' in apply_ic_species!", sp_name);
  }
  
  double t0 = luaL_optnumber(L, 3, app_lw->t_start);
  gkyl_moment_app_apply_ic_species(app_lw->app, sidx, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// The status table for reads has the following structure:
// {
//   io_status = true or false.
//   io_status_str = "success", "bad-version", "fopen-failed", "fread-failed", "data-mismatch"
//   frame = frame_num (integer)
//   stime = time in file read
// }

static const struct gkyl_str_int_pair rio_status[] = {
  { "success", GKYL_ARRAY_RIO_SUCCESS },
  { "bad-version", GKYL_ARRAY_RIO_BAD_VERSION },
  { "fopen-failed", GKYL_ARRAY_RIO_FOPEN_FAILED },
  { "fread-failed", GKYL_ARRAY_RIO_FREAD_FAILED },
  { "data-mismatch", GKYL_ARRAY_RIO_DATA_MISMATCH },
  { 0, 0 }
};

// Pushes table with status on stack. Table is left on stack.
static void
push_restart_status_table(lua_State *L, struct gkyl_app_restart_status status)
{
  lua_newtable(L);

  lua_pushstring(L, "io_status");
  lua_pushboolean(L, status.io_status == GKYL_ARRAY_RIO_SUCCESS ? true : false);
  lua_rawset(L, -3);

  lua_pushstring(L, "io_status_str");
  lua_pushstring(L, gkyl_search_str_int_pair_by_int(rio_status, status.io_status, "success"));
  lua_rawset(L, -3);

  lua_pushstring(L, "frame");
  lua_pushinteger(L, status.frame);
  lua_rawset(L, -3);

  lua_pushstring(L, "stime");
  lua_pushnumber(L, status.stime);
  lua_rawset(L, -3);  
}

// Read field from file. (file-name) -> status table. See above.
static int
mom_app_from_file_field(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  const char *fname = luaL_checkstring(L, 2);
  struct gkyl_app_restart_status status = gkyl_moment_app_from_file_field(app_lw->app, fname);

  push_restart_status_table(L, status);
  
  return 1;
}

// Read field from file. (file-name, species-name) -> status table. See above.
static int
mom_app_from_file_species(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  const char *fname = luaL_checkstring(L, 2);
  const char *sp_name = luaL_checkstring(L, 3);

  int sidx = app_find_species(app_lw->app, sp_name);
  if (sidx < 0) {
    return luaL_error(L, "Incorrect species name '%s' in from_file_species!", sp_name);
  }
  
  struct gkyl_app_restart_status status = gkyl_moment_app_from_file_species(app_lw->app, sidx, fname);

  push_restart_status_table(L, status);
  
  return 1;
}

// Read field from file. (frame) -> status table. See above.
static int
mom_app_from_frame_field(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  int frame = luaL_checkinteger(L, 2);
  struct gkyl_app_restart_status status = gkyl_moment_app_from_frame_field(app_lw->app, frame);

  push_restart_status_table(L, status);
  
  return 1;
}

// Read field from file. (frame, species-name) -> status table. See above.
static int
mom_app_from_frame_species(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  int frame = luaL_checkinteger(L, 2);
  const char *sp_name = luaL_checkstring(L, 3);

  int sidx = app_find_species(app_lw->app, sp_name);
  if (sidx < 0) {
    return luaL_error(L, "Incorrect species name '%s' in from_frame_species!", sp_name);
  }
  
  struct gkyl_app_restart_status status = gkyl_moment_app_from_frame_species(app_lw->app, sidx, frame);

  push_restart_status_table(L, status);
  
  return 1;
}

// Compute integrated moments. (tm) -> bool.
static int
mom_app_calc_integrated_mom(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_moment_app_calc_integrated_mom(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// get number of field-energy diagnostics stored () -> int.
static int
mom_app_field_energy_ndiag(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  lua_pushinteger(L, gkyl_moment_app_field_energy_ndiag(app_lw->app));
  return 1;
}

// Return the field energy as a table () -> table.
static int
mom_app_get_field_energy(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  int nvals = gkyl_moment_app_field_energy_ndiag(app_lw->app);
  double vals[nvals];
  gkyl_moment_app_get_field_energy(app_lw->app, vals);

  lua_createtable(L, nvals, 0);

  for (int i = 0; i < nvals; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, vals[i]);
    lua_rawset(L, -3);
  }  

  return 1;
}

// Compute integrated field energy (L2 norm of each field
// component). (tm) -> bool.
static int
mom_app_calc_field_energy(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_moment_app_calc_field_energy(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Write solution (field and species) to file (time, frame) -> bool.
static int
mom_app_write(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_moment_app_write(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write field to file (time, frame) -> bool.
static int
mom_app_write_field(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_moment_app_write_field(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write species solution to file (name, time, frame) -> bool.
static int
mom_app_write_species(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  const char *sp_name = luaL_checkstring(L, 2);
  int sidx = app_find_species(app_lw->app, sp_name);
  if (sidx < 0) {
    return luaL_error(L, "Incorrect species name '%s' in write_species!", sp_name);
  }

  double tm = luaL_checknumber(L, 3);
  int frame = luaL_checkinteger(L, 4);
  gkyl_moment_app_write_species(app_lw->app, sidx, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated moments to file () -> bool.
static int
mom_app_write_integrated_mom(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  gkyl_moment_app_write_integrated_mom(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated field energy to file () -> bool.
static int
mom_app_write_field_energy(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  gkyl_moment_app_write_field_energy(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write simulation statistics to JSON. () -> bool.
static int
mom_app_stat_write(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  gkyl_moment_app_stat_write(app_lw->app);

  lua_pushboolean(L, status);
  return 1;  
}

// Write data from simulation to file.
static void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_moment_app_write(app, t_curr, frame);
    gkyl_moment_app_write_field_energy(app);
    gkyl_moment_app_write_integrated_mom(app);
  }
}

// Calculate and append field energy to dynvector.
static void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_moment_app* app, double t_curr, double force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_moment_app_calc_field_energy(app, t_curr);
  }
}

// Calculate and append integrated moments to dynvector.
static void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_moment_app_calc_integrated_mom(app, t_curr);
  }
}

// Update status table has the following structure:
// {
//    success = true or false
//    dt_actual = actual time-step taken
//    dt_suggested = suggested stable time-step
// }

// Update the solution by a suggested time-step. (dt) -> update status
// table. See above. For details see the C API doc for this function.
static int
mom_app_update(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double dt = luaL_checknumber(L, 2);
  struct gkyl_update_status status = gkyl_moment_update(app_lw->app, dt);

  // Return status table on stack.
  lua_newtable(L);
  lua_pushstring(L, "success");
  lua_pushboolean(L, status.success);
  lua_rawset(L, -3);  

  lua_pushstring(L, "dt_actual");
  lua_pushnumber(L, status.dt_actual);
  lua_rawset(L, -3);    

  lua_pushstring(L, "dt_suggested");
  lua_pushnumber(L, status.dt_suggested);
  lua_rawset(L, -3);
  
  return 1;
}

// Step message context.
struct step_message_trigs {
  int log_count; // Number of times logging called.
  int tenth, p1c; 
  struct gkyl_tm_trigger log_trig; // 10% trigger.
  struct gkyl_tm_trigger log_trig_1p; // 1% trigger.
};

// Write log message to console.
static void
write_step_message(const struct gkyl_moment_app *app, struct step_message_trigs *trigs, int step, double t_curr, double dt_next)
{
  if (gkyl_tm_trigger_check_and_bump(&trigs->log_trig, t_curr)) {
    if (trigs->log_count > 0) {
      gkyl_moment_app_cout(app, stdout, " Step %6d at time %#11.8g.  Time-step  %.6e.  Completed %g%s\n", step, t_curr, dt_next, trigs->tenth * 10.0, "%");
    }
    else {
      trigs->log_count += 1;
    }
    
    trigs->tenth += 1;
  }
  if (gkyl_tm_trigger_check_and_bump(&trigs->log_trig_1p, t_curr)) {
    gkyl_moment_app_cout(app, stdout, "%d", trigs->p1c);
    trigs->p1c = (trigs->p1c+1) % 10;
  }
}

static void
show_help(const struct gkyl_moment_app *app)
{
  gkyl_moment_app_cout(app, stdout, "Moment script takes the following arguments:\n");
  gkyl_moment_app_cout(app, stdout, " -h   Print this help message and exit\n");
  gkyl_moment_app_cout(app, stdout, " -V   Show verbose output\n");
  gkyl_moment_app_cout(app, stdout, " -rN  Restart simulation from frame N\n");
  gkyl_moment_app_cout(app, stdout, " -sN  Only run N steps of simulation\n");

  gkyl_moment_app_cout(app, stdout, "\n");
}

static struct gkyl_tool_args *
tool_args_from_argv(int optind, int argc, char *const*argv)
{
  struct gkyl_tool_args *targs = gkyl_malloc(sizeof *targs);
  
  targs->argc = argc-optind;
  targs->argv = 0;

  if (targs->argc > 0) {
    targs->argv = gkyl_malloc(targs->argc*sizeof(char *));
      for (int i = optind, j = 0; i < argc; ++i, ++j) {
        targs->argv[j] = gkyl_malloc(strlen(argv[i])+1);
        strcpy(targs->argv[j], argv[i]);
      }
  }

  return targs;
}

// CLI parser for main script
struct script_cli {
  bool help; // show help
  bool step_mode; // run for fixed number of steps? (for valgrind/cuda-memcheck)
  int num_steps; // number of steps
  bool use_verbose; // Should we use verbose output?
  bool is_restart; // Is this a restarted simulation?
  int restart_frame; // Which frame to restart simulation from.  
  
  struct gkyl_tool_args *rest;
};

static struct script_cli
mom_parse_script_cli(struct gkyl_tool_args *acv)
{
  struct script_cli cli = {
    .help =- false,
    .step_mode = false,
    .num_steps = INT_MAX,
    .use_verbose = false,
    .is_restart = false,
    .restart_frame = 0,
  };
  
  coption_long longopts[] = {
    {0}
  };
  const char* shortopts = "+hVs:r:";

  coption opt = coption_init();
  int c;
  while ((c = coption_get(&opt, acv->argc, acv->argv, shortopts, longopts)) != -1) {
    switch (c) {
      case 'h':
        cli.help = true;
        break;
        
      case 'V':
        cli.use_verbose = true;
        break;

      case 's':
        cli.num_steps = atoi(opt.arg);
        break;
      
      case 'r':
        cli.is_restart = true;
        cli.restart_frame = atoi(opt.arg);
        break;        
        
      case '?':
        break;
    }
  }

  cli.rest = tool_args_from_argv(opt.ind, acv->argc, acv->argv);
  
  return cli;
}
    

// Run simulation. (num_steps) -> bool. num_steps is optional.
static int
mom_app_run(lua_State *L)
{
  bool ret_status = true;

  // Create app object.
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;
  struct gkyl_moment_app *app = app_lw->app;

  // Parse command lines arguments passed to input file.
  struct gkyl_tool_args *args = gkyl_tool_args_new(L);

  struct script_cli script_cli = mom_parse_script_cli(args);
  if (script_cli.help) {
    show_help(app);
    gkyl_tool_args_release(script_cli.rest);
    gkyl_tool_args_release(args);
    goto freeresources;
  }

  gkyl_tool_args_release(script_cli.rest);
  gkyl_tool_args_release(args);

  // Initial and final simulation times.
  double t_curr = app_lw->t_start, t_end = app_lw->t_end;
  long num_steps = script_cli.num_steps;

  gkyl_moment_app_cout(app, stdout, "Initializing Moments Simulation ...\n");

  // Initialize simulation.
  bool is_restart = script_cli.is_restart;
  int restart_frame = script_cli.restart_frame;

  int frame_curr = 0;
  if (is_restart) {
    struct gkyl_app_restart_status status = gkyl_moment_app_read_from_frame(app, restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_moment_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_moment_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_moment_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_moment_app_apply_ic(app, t_curr);
  }

  int num_frames = app_lw->num_frames;
  int field_energy_calcs = app_lw->field_energy_calcs;
  int integrated_mom_calcs = app_lw->integrated_mom_calcs;
  // Triggers for IO and logging.
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames, .tcurr = frame_curr * (t_end / num_frames), .curr = frame_curr };
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  struct step_message_trigs m_trig = {
    .log_count = 0,
    .tenth = t_curr > 0.0 ?  (int) floor(t_curr / t_end * 10.0) : 0.0,
    .p1c = t_curr > 0.0 ?  (int) floor(t_curr / t_end * 100.0) % 10 : 0.0,
    .log_trig = { .dt = t_end / 10.0, .tcurr = t_curr },
    .log_trig_1p = { .dt = t_end / 100.0, .tcurr = t_curr },
  };

  struct timespec tm_ic0 = gkyl_wall_clock();
  // Initialize simulation.
  calc_field_energy(&fe_trig, app, t_curr, false);
  calc_integrated_mom(&im_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);
  
  gkyl_moment_app_cout(app, stdout, "Initialization completed in %g sec\n\n", gkyl_time_diff_now_sec(tm_ic0));

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = app_lw->dt_failure_tol;
  int num_failures = 0, num_failures_max = app_lw->num_failures_max;

  bool use_verbose = script_cli.use_verbose;

  long step = 1;
  while ((t_curr < t_end) && (step <= num_steps)) {
    if (use_verbose) {
      gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    }
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    if (use_verbose) {
      gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    }
    
    if (!status.success) {
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_moment_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_moment_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_moment_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_moment_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_moment_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        write_data(&io_trig, app, t_curr, true);

        break;
      }
    }
    else {
      num_failures = 0;
    }

    if (!use_verbose) {
      write_step_message(app, &m_trig, step, t_curr, status.dt_suggested);
    }

    step += 1;
  }

  calc_field_energy(&fe_trig, app, t_curr, false);
  calc_integrated_mom(&im_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  gkyl_moment_app_cout(app, stdout, "\n\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(app, stdout, "Source updates took %g secs\n", stat.sources_tm);
  gkyl_moment_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);
  gkyl_moment_app_cout(app, stdout, "See log file for full statistics\n");

freeresources:

  lua_pushboolean(L, ret_status);
  return 1;
}

// Return ghost cell as table of 2 elements.
static int
mom_app_nghost(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;
  struct gkyl_moment_app *app = app_lw->app;

  int nghost[3] = { 0 };
  gkyl_moment_app_nghost(app, nghost);
  
  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushinteger(L, nghost[i]);
    lua_rawset(L, -3);
  }
  
  return 1;
}

// Clean up memory allocated for simulation.
static int
mom_app_gc(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  gkyl_moment_app_release(app_lw->app);
  gkyl_free(*l_app_lw);
  
  return 0;
}

// App constructor.
static struct luaL_Reg mom_app_ctor[] = {
  { "new",  mom_app_new },
  { 0, 0 }
};

// App methods
static struct luaL_Reg mom_app_funcs[] = {
  { "num_species", mom_app_num_species },
  { "species_name", mom_app_species_name },

  { "max_dt", mom_app_max_dt },  

  { "apply_ic", mom_app_apply_ic },
  { "apply_ic_field", mom_app_apply_ic_field },
  { "apply_ic_species", mom_app_apply_ic_species },
  
  { "from_file_field", mom_app_from_file_field },
  { "from_file_species", mom_app_from_file_species },

  { "from_frame_field", mom_app_from_frame_field },
  { "from_frame_species", mom_app_from_frame_species },

  { "write", mom_app_write },
  { "write_field", mom_app_write_field },
  { "write_species", mom_app_write_species },
  { "write_field_energy", mom_app_write_field_energy },
  { "write_integrated_mom", mom_app_write_integrated_mom },
  { "stat_write", mom_app_stat_write },

  { "calc_field_energy", mom_app_calc_field_energy },
  { "calc_integrated_mom", mom_app_calc_integrated_mom },

  { "update", mom_app_update },
  { "run", mom_app_run },

  // Some low-level functions, typically not used by ordinary users.
  { "nghost", mom_app_nghost },
  { "field_energy_ndiag", mom_app_field_energy_ndiag },
  { "get_field_energy", mom_app_get_field_energy },
  
  { 0, 0 }
};

static void
app_openlibs(lua_State *L)
{
  // Register top-level App.
  do {
    luaL_newmetatable(L, MOMENT_APP_METATABLE_NM);

    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, mom_app_gc);
    lua_settable(L, -3);

    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, mom_app_funcs);
    
    luaL_register(L, "G0.Moments.App", mom_app_ctor);
    
  }
  while (0);

  // Register Species input struct.
  do {
    luaL_newmetatable(L, MOMENT_SPECIES_METATABLE_NM);
    luaL_register(L, "G0.Moments.Species", mom_species_ctor);
  }
  while (0);

  // Register Field input struct.
  do {
    luaL_newmetatable(L, MOMENT_FIELD_METATABLE_NM);
    luaL_register(L, "G0.Moments.Field", mom_field_ctor);
  }
  while (0);

  // Add globals and other parameters.
  do {
    // Table to store references so the objects do not get garbage-collected.
    lua_newtable(L);
    lua_setglobal(L, "__moment_ref_table");
  }
  while (0);
}

void
gkyl_moment_lw_openlibs(lua_State *L)
{
  // Register types for moment reconstruction and wave limiter schemes.
  gkyl_register_moment_scheme_types(L);
  gkyl_register_wave_limiter_types(L);

  // Register types for Riemann solvers and Braginskii types.
  gkyl_register_euler_rp_types(L);
  gkyl_register_mhd_rp_types(L);
  gkyl_register_mhd_divb_types(L);
  gkyl_register_braginskii_types(L);
  gkyl_register_spacetime_gauge_types(L);
  gkyl_register_spacetime_slicing_types(L);
  gkyl_register_spacetime_evolution_types(L);
  gkyl_register_embed_geo_types(L);
  
  eqn_openlibs(L);
  spacetime_openlibs(L);
  app_openlibs(L);
}

#endif
