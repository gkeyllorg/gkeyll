# created and designed by jaime caballero for my phd at differ and the gkeyll team
# to mimic original https://github.com/jRoeltgen/radiation_operator/tree/main by J. Roeltgen
# code written with some help of ai 
#
# data_parser.py
# script to get a function Li(ne, Te) out of interpolated ADAS data
# this the equivalent of Li in eq. 1 of Roeltgen's paper

import numpy as np
from scipy.interpolate import RectBivariateSpline
import os
import glob

def load_roeltgen_formatted(species, charge_state, log10_ne_cm3=13.0, data_dir="data_formatted"):
    """
    Loads the author's formatted ADAS text files (e.g., plt96_h_formatted.txt).
    
    Columns expected:
    0: Charge State
    1: log10(Te [eV])
    2: log10(Ne [cm^-3])
    3: log10(Lz [W cm^3])
    """

    # search for file containing the species data
    # this data was formatted by the original author J. Roeltgen
    search_pattern = f"plt*_{species.lower()}_formatted.txt"
    full_search_path = os.path.join(data_dir, search_pattern)
    matching_files = glob.glob(full_search_path)
    if not matching_files:
        raise FileNotFoundError(f"No files matching pattern {full_search_path} found.")
    if len(matching_files) > 1:
        raise ValueError(f"Multiple files matching pattern {full_search_path} found: {matching_files}")
    filepath = matching_files[0]
    print(f"Loading data from: {filepath}")

    # Read data, skipping the text header line
    data = np.loadtxt(filepath, skiprows=1)
    
    charges = data[:, 0]
    log_te  = data[:, 1]
    log_ne  = data[:, 2]
    log_lz  = data[:, 3]

    # charges start at 1 in the files, adjust this:
    charge_state = charge_state + 1
    
    # 1. Filter by the requested Charge State
    mask_charge = (charges == charge_state)
    if not np.any(mask_charge):
        raise ValueError(f"Charge state {charge_state} not found for species {species}.")
        
    # 2. Filter by the nearest electron density
    # The paper standard is n_e = 1e19 m^-3, which is 1e13 cm^-3 (log10 = 13.0)
    unique_ne = np.unique(log_ne[mask_charge])
    closest_ne = unique_ne[np.argmin(np.abs(unique_ne - log10_ne_cm3))]
    mask_ne = np.isclose(log_ne, closest_ne)
    
    # Combine the masks
    final_mask = mask_charge & mask_ne
    
    # 3. Extract and convert to true physical values (SI Units)
    te_physical = 10 ** log_te[final_mask]
    
    # Lz is in W*cm^3. To convert to W*m^3, we multiply by 1e-6
    lz_physical = (10 ** log_lz[final_mask]) * 1e-6 
    
    return te_physical, lz_physical

def load_adas_plt_h(filepath):
    """
    Parses a standard ADAS ADF11 (PLT) file for Hydrogen.
    Returns an interpolation function: f(log10_ne, log10_te)
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Parse Grid Sizes (usually the first line)
    # Example: '  9  8' (9 Te points, 8 Ne points)
    header = lines[0].split()
    n_ne = int(header[1])
    n_te = int(header[2])

    # Collect all numbers, skipping metadata
    raw_data = []
    for line in lines[1:]:
        if 'RADIATED POWERS' in line or 'OPEN-ADAS' in line:
            break
        # Clean comments and handle negative numbers properly
        clean_line = line.split('/')[0].replace('-', ' -')
        parts = clean_line.split()
        for p in parts:
            try:
                raw_data.append(float(p.replace('D', 'E')))
            except ValueError:
                continue

    # Slicing ne first, then te
    log10_ne = np.array(raw_data[:n_ne])
    log10_te = np.array(raw_data[n_ne : n_ne + n_te])

    # Extract Matrix: shape (n_ne, n_te)
    data_start = n_ne + n_te
    matrix_vals = np.array(raw_data[data_start : data_start + (n_ne * n_te)])
    lz_matrix = matrix_vals.reshape((n_te, n_ne))
    
    # Create Interpolator: f(log10_ne, log10_te_ev)
    interp = RectBivariateSpline(log10_te, log10_ne, lz_matrix) 
    
    return interp


def get_lz_si(ne, te_ev, interp):
    """
    ne: density in m^-3
    te_ev: temperature in eV
    returns: Lz in W * m^3
    """
    # 1. Evaluate log10(W*cm^3)
    # ADAS ne is usually in cm^-3 (10^13 to 10^20 is typical)
    # If your input ne is m^-3, convert to cm^-3 for the lookup
    ne_cm3 = ne * 1e-6
    log_lz = interp(np.log10(te_ev), np.log10(ne_cm3))[0,0]
    
    # 2. Convert result to W*m^3
    # 10**log_lz gives W*cm^3. Multiply by 1e-6 to get W*m^3
    return (10**log_lz) * 1e-6

