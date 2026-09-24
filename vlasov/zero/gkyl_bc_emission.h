#pragma once

#include <gkyl_emission_spectrum_model.h>
#include <gkyl_emission_yield_model.h>
#include <gkyl_emission_elastic_model.h>

// context for use in computing emission from boundary
struct gkyl_bc_emission_ctx {
  int num_species;
  double t_bound;
  bool elastic;
  struct gkyl_emission_spectrum_model *spectrum_model[GKYL_MAX_SPECIES];
  struct gkyl_emission_yield_model *yield_model[GKYL_MAX_SPECIES];
  struct gkyl_emission_elastic_model *elastic_model;
  char in_species[GKYL_MAX_SPECIES][128];
};


/**
 * Create the ctx struct required for the emitting wall boundary condition
 *
 * @param num_species Number of impacting species causing emission
 * @param t_bound Time scaling factor for the emission
 * @param elastic Flag for elastic emission
 * @param spectrum_model Model type for emission spectrum
 * @param yield_model Model type for emission yield
 * @param elastic_model Model type for elastic emission
 * @param in_species Table of impacting species names
 * @return New ctx structure
 */
struct gkyl_bc_emission_ctx*
gkyl_bc_emission_new(int num_species, double t_bound, bool elastic,
  struct gkyl_emission_spectrum_model *spectrum_model[],
  struct gkyl_emission_yield_model *yield_model[],
  struct gkyl_emission_elastic_model *elastic_model, char in_species[][128]);

/**
 * Copper preset for secondary electron emission
 *
 * @param num_species Number of impacting species causing emission
 * @param t_bound Time scaling factor for the emission
 * @param in_species Table of impacting species names
 * @param use_gpu bool to determine if on GPU
 * @return New ctx structure
 */
struct gkyl_bc_emission_ctx*
gkyl_bc_emission_secondary_electron_copper_new(int num_species, double t_bound,
  char in_species[][128], bool use_gpu);

/**
 * Oxidized lithium preset for secondary electron emission
 *
 * @param num_species Number of impacting species causing emission
 * @param t_bound Time scaling factor for the emission
 * @param in_species Table of impacting species names
 * @param use_gpu bool to determine if on GPU
 * @return New ctx structure
 */
struct gkyl_bc_emission_ctx*
gkyl_bc_emission_secondary_electron_lithium_oxidized_new(int num_species, double t_bound,
  char in_species[][128], bool use_gpu);

/**
 * Clean (no oxidation) lithium preset for secondary electron emission
 *
 * @param num_species Number of impacting species causing emission
 * @param t_bound Time scaling factor for the emission
 * @param in_species Table of impacting species names
 * @param use_gpu bool to determine if on GPU
 * @return New ctx structure
 */
struct gkyl_bc_emission_ctx*
gkyl_bc_emission_secondary_electron_lithium_clean_new(int num_species, double t_bound,
  char in_species[][128], bool use_gpu);

/**
 * Copper preset for ion-impact secondary electron emission
 *
 * @param num_species Number of impacting species causing emission
 * @param t_bound Time scaling factor for the emission
 * @param in_species Table of impacting species names
 * @param use_gpu bool to determine if on GPU
 * @return New ctx structure
 */
struct gkyl_bc_emission_ctx*
gkyl_bc_emission_ion_impact_copper_new(int num_species, double t_bound,
  char in_species[][128], bool use_gpu);

struct gkyl_bc_emission_ctx* gkyl_bc_emission_secondary_electron_copper_new(int num_species,
  double t_bound, char in_species[][128], bool use_gpu);

struct gkyl_bc_emission_ctx* gkyl_bc_emission_ion_impact_copper_new(int num_species, double t_bound,
  char in_species[][128], bool use_gpu);

/**
 * Free memory associated with bc_emission struct.
 *
 * @param ctx BC ctx.
 */
void
gkyl_bc_emission_release(struct gkyl_bc_emission_ctx *ctx);
