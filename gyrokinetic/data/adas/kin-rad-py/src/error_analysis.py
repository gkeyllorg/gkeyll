# created and designed by jaime caballero for my phd at differ and the gkeyll team
# to mimic original https://github.com/jRoeltgen/radiation_operator/tree/main by J. Roeltgen
# code written with some help of ai 
#
# error_analysis.py
# this script mimics the error_analysis script from Roeltgen's MATLAB code
# here: https://github.com/jRoeltgen/radiation_operator/blob/main/error_analysis.m

import numpy as np

def error_analysis(ratio, te, Lz):
    """
    Evaluates if the optimization fit passes strict physical criteria.
    Equivalent to error_analysis.m
    """
    # Normalize the target data to its peak
    Lz_max = np.max(Lz)
    Lz_ratio = Lz / Lz_max

    #print("Lz Ratio (Lz / Lz_max):", Lz_ratio)
    
    # 1. Create Boolean Masks for different zones of the radiation curve
    # We only care about temperatures > 1 eV for these checks
    Lz_ind_1e2 = (Lz_ratio > 1e-2) & (te > 1)  # Top 1% of the peak
    Lz_ind_1e4 = (Lz_ratio > 1e-4) & (te > 1)  # Top 0.01%
    Lz_ind_1e8 = (Lz_ratio > 1e-8) & (te > 1)  # Deep tails

    #print(Lz_ratio)

    Lz_mag = np.abs(np.log10(Lz_ratio))
    
    # 2. Run the Success Checks (Start assuming all True)
    successes = [True] * 6
    
    if np.any(ratio[Lz_ind_1e4] > 2.0): 
        successes[0] = False
    if np.any(ratio[Lz_ind_1e8] > 8.0): 
        successes[1] = False
        
    # Error magnitude scaling check
    cond3 = (np.log10(ratio) > 0.5 * Lz_mag) & (Lz_mag > np.log10(8))
    if np.any(cond3): 
        successes[2] = False

    # condition 3 needs to be checked
    # successes[2] = True
        
    if np.any(ratio[Lz_ind_1e2] > 1.2): 
        successes[3] = False
    if np.any(ratio[Lz_ind_1e4] > 1.4): 
        successes[4] = False
    if np.any(ratio[Lz_ind_1e8] > 2.0): 
        successes[5] = False

    # 3. Calculate Maximum Errors per zone
    # We use 'if np.any()' to prevent NumPy from throwing an error if a mask is completely empty
    max_error = np.zeros(3)
    max_error[0] = np.max(ratio[Lz_ind_1e2]) if np.any(Lz_ind_1e2) else np.inf
    max_error[1] = np.max(ratio[Lz_ind_1e4]) if np.any(Lz_ind_1e4) else np.inf
    max_error[2] = np.max(ratio[Lz_ind_1e8]) if np.any(Lz_ind_1e8) else np.inf
    
    return successes, max_error