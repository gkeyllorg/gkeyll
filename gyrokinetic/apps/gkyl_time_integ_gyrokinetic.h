#pragma once

#include <gkyl_array.h>

/**
 * State pointers used by gyrokinetic RHS evaluation and time stepping.
 * The app owns the top-level pointer arrays; their entries borrow arrays
 * owned by other modules.
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
  const struct gkyl_array *aparin; // Input A_parallel (NULL if not evolved).
  struct gkyl_array *aparout; // Output dA_parallel/dt or updated A_parallel.
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
