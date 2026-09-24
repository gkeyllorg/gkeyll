#ifdef GKYL_HAVE_LUA

#include <gkyl_app.h>
#include <gkyl_lw_priv.h>
#include <gkyl_util.h>

// Define options for moments of a distribution function.
static const struct gkyl_str_int_pair distribution_moms[] = {
  { "M0",                 GKYL_F_MOMENT_M0 }, // Number density.
  { "M1",                 GKYL_F_MOMENT_M1 }, // Momentum density.
  { "M2",                 GKYL_F_MOMENT_M2 }, // Kinetic energy density.
  { "M2par",              GKYL_F_MOMENT_M2PAR }, // Parallel kinetic energy density.
  { "M2perp",             GKYL_F_MOMENT_M2PERP }, // Perpendicular kinetic energy density.
  { "M2ij",               GKYL_F_MOMENT_M2IJ }, // Kinetic energy tensor..
  { "M3",                 GKYL_F_MOMENT_M3 }, // Heat flux.
  { "M3par",              GKYL_F_MOMENT_M3PAR }, // Parallel energy flux.
  { "M3perp",             GKYL_F_MOMENT_M3PERP }, // Perpendicular energy flux.
  { "M3ijk",              GKYL_F_MOMENT_M3IJK }, // Heat flux in lab frame.
  { "MaxwellianMoments",  GKYL_F_MOMENT_MAXWELLIAN }, // M0, drift speed, T/m.
  { "BiMaxwellianMoments",GKYL_F_MOMENT_BIMAXWELLIAN }, // M0, drift speed, Tpar/m, Tperp/m.
  { "LTEMoments",         GKYL_F_MOMENT_LTE }, // Maxwellian or Maxwell-Juttner moments.
  { "M0M1M2",             GKYL_F_MOMENT_M0M1M2 },  // M0, M1, M2.
  { "M0M1M2parM2perp",    GKYL_F_MOMENT_M0M1M2PARM2PERP },  // M0, M1, M2par, M2perp.
  { "HamiltonianMoments", GKYL_F_MOMENT_HAMILTONIAN },  // M0, mass*M1, H moments.
  { "M1_from_H",          GKYL_F_MOMENT_M1_FROM_H }, // dH/dv / m moment.
  { "EnergyMoment",       GKYL_F_MOMENT_ENERGY }, // H moment.
  { "M0EnergyM3",         GKYL_F_MOMENT_M0ENERGYM3 }, // M0, Energy (H) and M3 moments.
  { "Ni",                 GKYL_F_MOMENT_NI }, // M0, M1i for-vector.
  { "Tij",                GKYL_F_MOMENT_TIJ }, // Stress-energy tensor.
  { 0, 0 }
};

// Species boundary conditions -> enum map.
static const struct gkyl_str_int_pair species_bcs[] = {
  { "bcCopy", GKYL_SPECIES_COPY },
  { "bcWall", GKYL_SPECIES_REFLECT },
  { "bcReflect", GKYL_SPECIES_REFLECT },
  { "bcAbsorb", GKYL_SPECIES_ABSORB },
  { "bcNoSlip", GKYL_SPECIES_NO_SLIP },
  { "bcWedge", GKYL_SPECIES_WEDGE },
  { "bcFunc", GKYL_SPECIES_FUNC },
  { "bcFixedFunc", GKYL_SPECIES_FIXED_FUNC },
  { "bcZeroFlux", GKYL_SPECIES_ZERO_FLUX },
  { "bcRecycle", GKYL_SPECIES_RECYCLE },
  { 0, 0 }
};

// Field boundary conditions -> enum map.
static const struct gkyl_str_int_pair field_bcs[] = {
  { "bcCopy", GKYL_FIELD_COPY }, 
  { "bcWall", GKYL_FIELD_PEC_WALL },
  { "bcPECWall", GKYL_FIELD_PEC_WALL },
  { "bcSymWall", GKYL_FIELD_SYM_WALL },
  { "bcReservoir", GKYL_FIELD_RESERVOIR },
  { "bcWedge", GKYL_FIELD_WEDGE },
  { "bcFunc", GKYL_FIELD_FUNC },
  { 0, 0 }
};

void
gkyl_register_distribution_moment_types(lua_State *L)
{
  register_types(L, distribution_moms, "Moment");
}

void
gkyl_register_species_bc_types(lua_State *L)
{
  register_types(L, species_bcs, "SpeciesBc");
}

void
gkyl_register_field_bc_types(lua_State *L)
{
  register_types(L, field_bcs, "FieldBc");
}

void
gkyl_lw_eval_cb(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct lua_func_ctx *fr = ctx;
  lua_State *L = fr->L;

  int ndim = fr->ndim;
  int nret = fr->nret;
  lua_rawgeti(L, LUA_REGISTRYINDEX, fr->func_ref);
  lua_pushnumber(L, t);
  lua_createtable(L, GKYL_MAX_DIM, 0);

  for (int i=0; i<ndim; ++i) {
    lua_pushnumber(L, xn[i]);
    lua_rawseti(L, -2, i+1); 
  }

  if (lua_pcall(L, 2, nret, 0)) {
    const char* ret = lua_tostring(L, -1);
    luaL_error(L, "*** gkyl_lw_eval_cb ERROR: %s\n", ret);
  }

  for (int i=nret-1; i>=0; --i) { // need to fetch in reverse order
    fout[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }  
}

#endif
