#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

#include <stdbool.h>

// Object type.
typedef struct gkyl_deflated_fem_gyroaverage gkyl_deflated_fem_gyroaverage;

// Inputs to create a new deflated gyroaverage updater.
struct gkyl_deflated_fem_gyroaverage_inp {
  const struct gkyl_rect_grid *grid; // Configuration space grid, (x,z) or (x,y,z).
  const struct gkyl_basis *basis; // Configuration space basis (piecewise linear).
  const struct gkyl_range *local; // Local range (decomposed only along z).
  bool is_periodic[2]; // Whether each perpendicular direction is periodic.
  // Extend the field beyond non-periodic boundaries by reflection instead of linearly about
  // its boundary value (see gkyl_fem_gyroaverage_inp).
  bool reflect_bnd;
  const struct gkyl_array *gxx; // Contravariant metric coefficient g^xx.
  const struct gkyl_array *gxy; // Contravariant metric coefficient g^xy (3x only).
  const struct gkyl_array *gyy; // Contravariant metric coefficient g^yy (3x only).
  const struct gkyl_array *bmag; // Magnetic field amplitude.
  double mass; // Species mass.
  double charge; // Species charge.
  int num_mu; // Number of magnetic moments to gyroaverage at.
  const double *mu; // Magnetic moments (m*vperp^2/(2*B)) to gyroaverage at.
  // Whether to run on the GPU. The geometry arrays may be on the host or on the device (the
  // operator is assembled on the host from a copy), and the fields are then on the device.
  bool use_gpu;
};

/**
 * Create a new updater that gyroaverages a field in the planes perpendicular
 * to the magnetic field (see gkyl_fem_gyroaverage), in a 2x (x,z) or 3x (x,y,z)
 * domain whose last coordinate, z, is along the magnetic field. The field is
 * gyroaveraged on the surfaces of constant z between cells (the z nodes), with
 * the geometry evaluated there, and is linear in z within each cell. The
 * gyroaverage of a field that is continuous along z is thus continuous along
 * z. The geometry should be continuous along z.
 *
 * @param inp Input parameters defined in gkyl_deflated_fem_gyroaverage_inp.
 * @return New updater pointer.
 */
struct gkyl_deflated_fem_gyroaverage *gkyl_deflated_fem_gyroaverage_inew(
  const struct gkyl_deflated_fem_gyroaverage_inp *inp
);

/**
 * Gyroaverage a field that is continuous along z at each magnetic moment.
 *
 * @param up Gyroaverage updater.
 * @param fin Input field (num_basis components).
 * @param fout Gyroaveraged field at each magnetic moment (num_mu*num_basis
 *             components, the coefficients for the n-th magnetic moment
 *             starting at n*num_basis).
 */
void gkyl_deflated_fem_gyroaverage_advance(
  struct gkyl_deflated_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
);

/**
 * Gyroaverage a field that may be discontinuous along z at each magnetic moment, cell by cell
 * along z: in each cell the field is linear in z between its values at the lower and upper z
 * surfaces of the cell, which are gyroaveraged on those surfaces. For a field that is
 * continuous along z this equals gkyl_deflated_fem_gyroaverage_advance; in general the output
 * is discontinuous along z, and gkyl_deflated_fem_gyroaverage_advance_adjoint is the exact
 * adjoint of this operation.
 *
 * @param up Gyroaverage updater.
 * @param fin Input field (num_basis components).
 * @param fout Gyroaveraged field at each magnetic moment (num_mu*num_basis
 *             components, the coefficients for the n-th magnetic moment
 *             starting at n*num_basis).
 */
void gkyl_deflated_fem_gyroaverage_advance_dg(
  struct gkyl_deflated_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
);

/**
 * Apply the adjoint of gkyl_deflated_fem_gyroaverage_advance, i.e. compute
 * fout = sum_n G_n^dagger f_n, G_n being the gyroaverage at the n-th magnetic
 * moment, such that for any field g that is continuous along z
 *
 *   sum_n int f_n G_n(g) dx = int fout g dx
 *
 * holds to machine precision. The input need not be continuous along z, and
 * the output is in general discontinuous along z.
 *
 * @param up Gyroaverage updater.
 * @param fin Input field at each magnetic moment (num_mu*num_basis
 *            components, the coefficients for the n-th magnetic moment
 *            starting at n*num_basis).
 * @param fout Sum of the adjoint gyroaverages (num_basis components).
 */
void gkyl_deflated_fem_gyroaverage_advance_adjoint(
  struct gkyl_deflated_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_deflated_fem_gyroaverage_release(struct gkyl_deflated_fem_gyroaverage *up);
