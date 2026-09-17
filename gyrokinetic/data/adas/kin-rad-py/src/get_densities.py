# created and designed by jaime caballero for my phd at differ and the gkeyll team
# to mimic original https://github.com/jRoeltgen/radiation_operator/tree/main by J. Roeltgen
# code written with some help of ai 
#
# get_densities.py
# a quick utility script to parse the legacy roeltgen database and extract
# the exact n_e density grids used for each element and charge state.
# it outputs a perfectly formatted python dictionary to paste into fit_batch.py.

import os

ATOMIC_NUMBERS = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Al':13, 'Ar': 18}
# reverse the map to look up species by their atomic number Z
Z_TO_SPECIES = {v: k for k, v in ATOMIC_NUMBERS.items()}

def extract_density_grid(filepath):
    density_grid = {}
    
    if not os.path.exists(filepath):
        print(f"[error] could not find {filepath}")
        return None
        
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
    current_species = None
    current_charge = None
    intervals_left = 0
    lines_to_skip = 0
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        # skip the raw T_e and L_z arrays
        if lines_to_skip > 0:
            lines_to_skip -= 1
            continue
            
        if line.startswith('***Atomic number='):
            # parse the Z number
            z = int(line.split(';')[0].split('=')[1])
            current_species = Z_TO_SPECIES.get(z)
            if current_species and current_species not in density_grid:
                density_grid[current_species] = {}
                
        elif line.startswith('**Charge state='):
            if not current_species:
                continue
            
            # parse the charge state and how many densities we need to read
            parts = line.split(',')
            c_gkeyll = int(parts[0].split('=')[1])
            intervals_left = int(parts[1].split('=')[1])
            
            # convert 1-based gkeyll indexing to 0-based indexing
            current_charge = c_gkeyll - 1
            
            if current_charge not in density_grid[current_species]:
                density_grid[current_species][current_charge] = []
                
        elif intervals_left > 0 and line[0].isdigit() and 'e' in line:
            # parse the density (it is the first number in the parameter line)
            parts = line.split()
            ne_m3_log = float(parts[0])
            
            # convert from m^-3 to cm^-3 and round cleanly (e.g. 19.0 -> 13.0)
            ne_cm3_log = round(ne_m3_log - 6.0, 1)
            
            density_grid[current_species][current_charge].append(ne_cm3_log)
            
            intervals_left -= 1
            lines_to_skip = 2 # skip the next 2 lines (the T_e and L_z arrays)
            
    return density_grid

if __name__ == "__main__":
    filepath = os.path.join("fits_data", "roeltgen_data", "radiation_fit_parameters.txt")
    
    print("scanning legacy database for density grids...\n")
    grid = extract_density_grid(filepath)
    
    if grid:
        print("DENSITY_GRID = {")
        for species in ATOMIC_NUMBERS.keys():
            if species in grid:
                charges = grid[species]
                # format the charge dictionaries cleanly
                charge_strs = []
                for charge, dens in sorted(charges.items()):
                    charge_strs.append(f"{charge}: {dens}")
                
                dict_str = ", ".join(charge_strs)
                print(f"    '{species}': {{ {dict_str} }},")
        print("}")
        print("\n-> copy and paste the dictionary above into fit_batch.py!")