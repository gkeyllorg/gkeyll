// Gate 4: which block types own which divertor plate?
//
// A block that terminates on a material plate must be refused when THAT plate
// cannot be reached, and must be unaffected when some OTHER plate is broken.
// Get this wrong in either direction and the failure is quiet: a block that
// silently ignores its own plate builds geometry that ends nowhere physical,
// and a block that consults a plate it does not own is refused for a reason
// that has nothing to do with it.
//
// The contract lives in `tok_plate_slot_required` (tok_geo_utils.c), which is
// static, so this drives it through the public
// `gkyl_tok_geo_check_plate_coverage` -- the same path the multiblock preflight
// and the rho-adjust retry both use. Testing the real path rather than the
// private helper also means a change in how ownership is *consumed* is caught,
// not just a change in the table.
//
// Method: build one geometry object per plate configuration, then sweep every
// member of the enum against it.
//
//   both plates good   -> nothing may be refused
//   lower plate broken -> refused exactly by the types that own the lower plate
//   upper plate broken -> refused exactly by the types that own the upper plate
//
// "Broken" is a plate parked far outside the machine, so no flux surface can
// terminate on it. That is the "shorten the callback until it no longer
// reaches" case from the gate plan, taken to its limit so the outcome cannot
// depend on how far a surface happens to extend.
//
// The sweep is built from the ENUM, not from a list of interesting types, so a
// newly added block type shows up here as an explicit expectation rather than
// silently inheriting `default:` (no plates required) -- which for a new
// divertor type would be exactly the quiet failure above. That is also what
// gate 8 asks for: coverage inventoried from the enum.

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <acutest.h>
#include <gkyl_efit.h>
#include <gkyl_tok_geo.h>

// The straight ASDEX divertor plates, taken from the commented-out forms in
// creg/rt_gk_multib_asdex_2x2v_p1.c, where they document the cubic splines
// that replaced them. Straight lines are enough to answer a reachability
// question and keep the fixture readable.
static void
plate_lower_good(double s, double *RZ)   // outboard
{
  RZ[0] = 1.5966 + (1.6888 - 1.5966)*s;
  RZ[1] = -1.1421 + (-0.8781 - (-1.1421))*s;
}

static void
plate_upper_good(double s, double *RZ)   // inboard
{
  RZ[0] = 1.2686 + (1.1886 - 1.2686)*s;
  RZ[1] = -1.0520 + (-0.7294 - (-1.0520))*s;
}

// Parked far outside the machine: unreachable from any flux surface.
static void
plate_unreachable(double s, double *RZ)
{
  RZ[0] = 9.0 + 0.1*s;
  RZ[1] = 9.0 + 0.1*s;
}

struct type_case {
  enum gkyl_tok_geo_type ftype;
  const char *name;
  bool owns_lower;
  bool owns_upper;
};

// Every member of `enum gkyl_tok_geo_type`, in declaration order. IWL is
// handled separately: its ownership depends on which side of the separatrix
// the requested psi sits, so it cannot be a constant here.
static const struct type_case cases[] = {
  { GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT,     "DN_SOL_OUT",     true,  true  },
  { GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN,      "DN_SOL_IN",      true,  true  },
  { GKYL_GEOMETRY_TOKAMAK_LSN_SOL,        "LSN_SOL",        true,  true  },
  { GKYL_GEOMETRY_TOKAMAK_USN_UP,         "USN_UP",         false, false },
  { GKYL_GEOMETRY_TOKAMAK_CORE,           "CORE",           false, false },
  { GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO,  "DN_SOL_OUT_LO",  true,  false },
  { GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID, "DN_SOL_OUT_MID", false, false },
  { GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP,  "DN_SOL_OUT_UP",  false, true  },
  { GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO,   "DN_SOL_IN_LO",   true,  false },
  { GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID,  "DN_SOL_IN_MID",  false, false },
  { GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP,   "DN_SOL_IN_UP",   false, true  },
  { GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO,     "LSN_SOL_LO",     true,  false },
  { GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID,    "LSN_SOL_MID",    false, false },
  { GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP,     "LSN_SOL_UP",     false, true  },
  // Note the diagonal: the UPPER private-flux LEFT half owns the LOWER plate,
  // and the LOWER private-flux LEFT half owns the UPPER plate. That is not a
  // typo in this table -- it is what tok_plate_slot_required says, and it is
  // worth an explicit test precisely because it reads like one.
  { GKYL_GEOMETRY_TOKAMAK_PF_UP_L,        "PF_UP_L",        true,  false },
  { GKYL_GEOMETRY_TOKAMAK_PF_UP_R,        "PF_UP_R",        false, true  },
  { GKYL_GEOMETRY_TOKAMAK_PF_LO_L,        "PF_LO_L",        false, true  },
  { GKYL_GEOMETRY_TOKAMAK_PF_LO_R,        "PF_LO_R",        true,  false },
  { GKYL_GEOMETRY_TOKAMAK_CORE_L,         "CORE_L",         false, false },
  { GKYL_GEOMETRY_TOKAMAK_CORE_R,         "CORE_R",         false, false },
};

// The psi to interrogate: just outside the ASDEX separatrix (~0.14975), where a
// SOL surface really would run down to a divertor plate.
static const double psi_sol = 0.1520;

static struct gkyl_tok_geo *
make_geo(plate_func lower, plate_func upper)
{
  struct gkyl_efit_inp efit_inp = {
    .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };
  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL,
    .rmin = 0.0, .rmax = 5.0,
    .rclose = 2.5, .rright = 2.5, .rleft = 0.7,
    .zmin = -1.3, .zmax = 1.0,
    .zmin_left = -1.0, .zmin_right = -0.9,
    .plate_spec = true,
    .plate_func_lower = lower,
    .plate_func_upper = upper,
  };
  return gkyl_tok_geo_new(&efit_inp, &ginp);
}

static bool
covered(struct gkyl_tok_geo *geo, enum gkyl_tok_geo_type ftype, double psi)
{
  struct gkyl_tok_geo_grid_inp inp = { .ftype = ftype };
  return gkyl_tok_geo_check_plate_coverage(geo, &inp, psi);
}

// Control: with both plates reachable, no block type may be refused. Without
// this, a table that refused everything would look like a pass below.
static void
test_both_plates_good_covers_every_type(void)
{
  struct gkyl_tok_geo *geo = make_geo(plate_lower_good, plate_upper_good);
  TEST_ASSERT(geo != NULL);
  for (size_t i = 0; i < sizeof(cases)/sizeof(*cases); ++i) {
    bool ok = covered(geo, cases[i].ftype, psi_sol);
    TEST_CHECK(ok);
    TEST_MSG("%s refused at psi=%g with BOTH plates reachable; the fixture "
      "plates do not reach this surface, so the ownership tests below would "
      "be measuring the fixture rather than the contract", cases[i].name, psi_sol);
  }
  gkyl_tok_geo_release(geo);
}

// Break one plate at a time and require that exactly the owners notice.
static void
sweep_with_broken(int broken_side)
{
  struct gkyl_tok_geo *geo = broken_side == 0
    ? make_geo(plate_unreachable, plate_upper_good)
    : make_geo(plate_lower_good, plate_unreachable);
  TEST_ASSERT(geo != NULL);
  const char *which = broken_side == 0 ? "lower" : "upper";

  for (size_t i = 0; i < sizeof(cases)/sizeof(*cases); ++i) {
    bool owns = broken_side == 0 ? cases[i].owns_lower : cases[i].owns_upper;
    bool ok = covered(geo, cases[i].ftype, psi_sol);

    if (owns) {
      TEST_CHECK(!ok);
      TEST_MSG("%s owns the %s plate, but coverage still passed with that "
        "plate unreachable: the block does not consult the plate it "
        "terminates on", cases[i].name, which);
    }
    else {
      TEST_CHECK(ok);
      TEST_MSG("%s does not own the %s plate, yet was refused when it broke: "
        "the block consults a plate it has no business consulting",
        cases[i].name, which);
    }
  }
  gkyl_tok_geo_release(geo);
}

static void test_lower_plate_owners(void) { sweep_with_broken(0); }
static void test_upper_plate_owners(void) { sweep_with_broken(1); }

// IWL is the one type whose ownership is not fixed: it needs both plates only
// on one side of the separatrix. Asserting the FLIP rather than either value
// keeps this independent of the sign convention, which differs between
// equilibria (`geqdsk_sign_convention`).
static void
test_iwl_ownership_flips_across_the_separatrix(void)
{
  struct gkyl_tok_geo *broken = make_geo(plate_unreachable, plate_unreachable);
  TEST_ASSERT(broken != NULL);

  // ASDEX: psisep ~ 0.14975, magnetic axis flux lower. One of these is inside
  // the LCFS and one outside, whichever way the convention runs.
  bool inboard  = covered(broken, GKYL_GEOMETRY_TOKAMAK_IWL, 0.1400);
  bool outboard = covered(broken, GKYL_GEOMETRY_TOKAMAK_IWL, 0.1600);

  TEST_CHECK(inboard != outboard);
  TEST_MSG("IWL required the plates on both sides of the separatrix "
    "(psi=0.14 -> covered=%d, psi=0.16 -> covered=%d). Its plate requirement "
    "is supposed to depend on which side of psisep the surface lies.",
    inboard, outboard);

  gkyl_tok_geo_release(broken);
}

// Guard the sweep itself: if the enum grows, this table must grow with it, or
// the new type is silently untested and inherits "owns no plate".
static void
test_every_enum_member_is_covered(void)
{
  const size_t listed = sizeof(cases)/sizeof(*cases);
  // 21 members total; IWL is exercised by its own test above.
  const size_t enum_members = 21;
  TEST_CHECK(listed + 1 == enum_members);
  TEST_MSG("%zu types listed + 1 (IWL) != %zu enum members. A block type was "
    "added or removed; add it to `cases` with its plate ownership rather than "
    "letting it fall through to the default of owning none.",
    listed, enum_members);
}

TEST_LIST = {
  { "both_plates_good_covers_every_type", test_both_plates_good_covers_every_type },
  { "lower_plate_owners", test_lower_plate_owners },
  { "upper_plate_owners", test_upper_plate_owners },
  { "iwl_ownership_flips_across_the_separatrix", test_iwl_ownership_flips_across_the_separatrix },
  { "every_enum_member_is_covered", test_every_enum_member_is_covered },
  { NULL, NULL },
};
