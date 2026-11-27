#pragma once

#include <gkyl_array.h>

/**
 * State pointers used by gyrokinetic RHS evaluation and time stepping.
 * The app owns the top-level pointer arrays; their entries borrow arrays
 * owned by other modules. RHS evaluation writes instantaneous moment rates to
 * smom_out. Time stepping advances moment integrals with the same arithmetic as
 * distributions; smom_in holds the input integral state. Each step (and retry)
 * starts with zero smom. Accepted smoms remain in step_moms.f; separate
 * step_moms.rate arrays hold their averages over the accepted dt for consumers.
 */
struct gkyl_gyrokinetic_fdot_args {
  int num_species; // Number of charged species (per block for multiblock).
  int num_neut_species; // Number of neutral species (per block for multiblock).
  const struct gkyl_array **fin; // Charged input states.
  const struct gkyl_array **fin_neut; // Neutral input states.
  struct gkyl_array **fout; // Charged rates or updated states.
  struct gkyl_array **fout_neut; // Neutral rates or updated states.
  struct gkyl_array ***smom_in; // Charged stepped-moment lists.
  struct gkyl_array ***smom_in_neut; // Neutral stepped-moment lists.
  struct gkyl_array ***smom_out; // Charged output stepped-moment lists.
  struct gkyl_array ***smom_out_neut; // Neutral output stepped-moment lists.
};

/**
 * Multiblock state pointers with the same ownership as the single-block args.
 * Charged entries use b*num_species+i; neutral entries use b*num_neut_species+i,
 * where b is the local block index. These arrays do not own child-app storage.
 */
struct gkyl_gyrokinetic_multib_fdot_args {
  int nblocks_local; // Number of blocks on this rank.
  int num_species; // Number of charged species (per block for multiblock).
  int num_neut_species; // Number of neutral species (per block for multiblock).
  const struct gkyl_array **fin; // Charged input states.
  const struct gkyl_array **fin_neut; // Neutral input states.
  struct gkyl_array **fout; // Charged rates or updated states.
  struct gkyl_array **fout_neut; // Neutral rates or updated states.
  struct gkyl_array ***smom_in; // Charged stepped-moment lists.
  struct gkyl_array ***smom_in_neut; // Neutral stepped-moment lists.
  struct gkyl_array ***smom_out; // Charged output stepped-moment lists.
  struct gkyl_array ***smom_out_neut; // Neutral output stepped-moment lists.
};
