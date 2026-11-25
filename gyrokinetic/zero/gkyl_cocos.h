#pragma once

// Equilibrium convention parameters from Sauter and Medvedev, Table I.
// Plasma-current and toroidal-field signs belong to the equilibrium, not this descriptor.
struct gkyl_cocos {
  int index; // COCOS identifier.
  int sigma_Bp_cocos; // Sign multiplying grad(phi) cross grad(psi).
  int e_Bp_cocos; // 0 for flux divided by 2*pi, 1 for full poloidal flux.
  int sigma_RvarphiZ_cocos; // Handedness of (R,phi,Z).
  int sigma_rhothetavarphi_cocos; // Handedness of reference (rho,theta,phi).
};

// Fixed internal equilibrium convention; computational theta is field-aligned.
// Changing this descriptor alone does not implement another solver convention.
extern const struct gkyl_cocos gkyl_cocos_tokamak;
