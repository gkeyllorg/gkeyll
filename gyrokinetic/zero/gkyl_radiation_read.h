#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <gkyl_array.h>

#define BUFFER_LEN 100
// Data for a single radiation fit.
struct rad_fit_parameters {
  // A, alpha, beta, gamma, and V0 are fitting parameters.
  // C = 8/sqrt(pi)*(2*charge/mass)^(gamma/2)
  // D = A*(alpha+beta)/C
  // vmag = sqrt(vpar^2 + 2*B*mu/mass)
  // nu(vpar,mu) = D*vmag^(gamma)/(beta*(vmag/V0)^-alpha + alpha*(vmag/V0)^beta)
  double electron_density;  // Electron density at which this fit is at
  double A;
  double alpha;
  double beta;
  double gamma;
  double V0;
  int te_intervals;  // Number of temperature intervals for fit emissivity
  double *te;  // electron temperatures at which the fit emissivity is calculated
  double *Lz;  // fit emissivity when assuming a maxwellian of corresponding temperature
};

// Radiation data for a single charge state (for all electron densities).
struct radiating_state {
  bool state_exists;
  int atomic_number;
  int charge_state;
  int number_of_densities;
  struct gkyl_array *electron_densities; // The different densities for this charge state
  struct rad_fit_parameters *rad_fits; // 1 fit per density for this charge state
};

// All of the radiation fits read.
struct all_radiation_states {
  int max_atomic_number;
  struct radiating_state *all_states; // indecies are z-1,charge_state
};

/**
 * Function to read in all the radiation fit parameters stored in "radiation_fit_params.txt".
 *
 * @return pointer to structure of all_radiation_states
 */
struct all_radiation_states* gkyl_radiation_read_rad_fit_params();

/**
 * Get the number of densities used in the fits.
 *
 * @param rad data: Struct containing radiation fit data.
 * @param atomic_z: Z of element for desired fit information.
 * @param charge_state: charge state of element for desired fit information.
 * @param min_ne: Desired minimum density (closest density is used).
 * @param max_ne: Desired maximum density (closest density is used).
 * @param num_densities: maximum number of densities to return fit parameters for.
 * @return 1 if fit doesn't exist
 */
int gkyl_radiation_read_get_num_densities(const struct all_radiation_states rad_data,
  int atomic_z, int charge_state, double min_ne, double max_ne, int *num_densities);

/**
 * Function to return the fit information for a specfied atomic number, charge state, and ne.
 *
 * @param rad data: Struct containing radiation fit data
 * @param atomic_z: Z of element for desired fit information
 * @param charge_state: charge state of element for desired fit information
 * @param a, alpha, beta, gamma, V0: fit parameters to be returned 
 * @param num_densities: maximum number of densities to return fit parameters for
 * @param electron densities: Array of electron densities
 * @param ref_dens: Reference electron density - choose closest fit density to this
 * @param min_ne: Desired minimum density (closest density is used)
 * @param max_ne: Desired maximum density (closest density is used)
 * @return 1 if fit doesn't exist
 */
int gkyl_radiation_read_get_fit_params(const struct all_radiation_states rad_data, int atomic_z,
  int charge_state, double *a, double *alpha, double *beta, double *gamma, double *V0,
  int *num_densities, double *electron_densities, double ref_dens, double min_ne, double max_ne);

/**
 * Function to return the fit emissivity (Lz) and temperature closest to a given input temperature.
 *
 * @param rad data: Struct containing radiation fit data
 * @param atomic_z: Z of element for desired fit information
 * @param charge_state: charge state of element for desired fit information
 * @param ne: The returned te and lz are for density closest to ne.
 * @param te: returns closest fit temperature to input te
 * @param Lz: returns Lz for closest temperature to input te 
 */
int gkyl_radiation_read_get_fit_lz(const struct all_radiation_states rad_data,
  int atomic_z, int charge_state, double ne, double* te, double* Lz);

/**
 * Free memory of all_radiation_states.
 *
 * @param rad data: Struct containing radiation fit data
 */
void gkyl_radiation_read_release_fit_params(struct all_radiation_states *rad_data);

/**
 * Find the nearest index in an array to given value (assumes 1 component in array)
 * @param arr A sorted array
 * @param target The value to find the closest index to
 * @return The index of the array with the value closest to target
 */
static inline int gkyl_radiation_read_find_nearest_idx(const struct gkyl_array* arr, double target){
  int left = 0;
  int right = arr->size - 1;
  double *data = (double*)arr->data;
  while (left < right) {
    if (fabs(data[left] - target)
	<= fabs(data[right] - target)) {
      right--;
    }
    else {
      left++;
    }
  }
  return left;
}
