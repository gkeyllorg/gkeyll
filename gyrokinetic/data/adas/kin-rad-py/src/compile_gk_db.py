# created and designed by jaime caballero for my phd at differ and the gkeyll team
# to mimic original https://github.com/jRoeltgen/radiation_operator/tree/main by J. Roeltgen
# code written with some help of ai 
#
# compile_db.py
# this script acts as the "linker" for the repository. it reads the master csv
# database, fetches the original unscaled adas data arrays, and compiles the 
# strict gkeyll-ready .txt file. 

import os
import sys
import csv
import argparse

# hook into our engine scripts to grab the raw data arrays
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))
from data_parser import load_roeltgen_formatted

# map species to atomic numbers for the gkeyll compiler
ATOMIC_NUMBERS = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 
    'F': 9, 'Ne': 10, 'Ar': 18, 'Kr': 36, 'Xe': 54, 'W': 74
}

def main():
    parser = argparse.ArgumentParser(description="compile gkeyll .txt database from master csv.")
    parser.add_argument("--id", type=str, required=True, help="dataset ID to compile (e.g., 'PROD_RUN')")
    args = parser.parse_args()

    base_dir = os.path.join("fits_data", f"fit-db_{args.id}")
    csv_path = os.path.join(base_dir, f"master-fit_db_{args.id}.txt")
    out_path = os.path.join(base_dir, f"gkeyll-db_fit_{args.id}.txt")

    if not os.path.exists(csv_path):
        print(f"[error] could not find master database at {csv_path}")
        sys.exit(1)

    print(f"=========================================")
    print(f"--- COMPILING GKEYLL DATABASE ---")
    print(f"--- DATASET ID: {args.id} ---")
    print(f"=========================================\n")

    memory_db = {}
    row_count = 0

    print("reading csv and fetching array data...")
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            species = row['element']
            charge = int(row['charge'])
            density = float(row['density'])
            min_te = float(row['min_te'])
            
            # build the nested dictionary structure
            if species not in memory_db: memory_db[species] = {}
            if charge not in memory_db[species]: memory_db[species][charge] = {}
            
            try:
                # fetch the original arrays for the gkeyll test plots
                Te_data, target_data_unscaled = load_roeltgen_formatted(species, charge, density, "formatted_data")
                
                # apply the temperature chop if one was used during the fit
                if min_te > 0.0:
                    mask = Te_data >= min_te
                    Te_data = Te_data[mask]
                    target_data_unscaled = target_data_unscaled[mask]
                    
            except Exception as e:
                print(f"[warning] missing raw data for {species}^{charge}+ @ 10^{density}. skipping entry.")
                continue

            # stash everything exactly how the compiler expects it
            memory_db[species][charge][density] = {
                'params': [float(row['A']), float(row['alpha']), float(row['beta']), float(row['V0']), float(row['gamma'])],
                'te': Te_data.tolist(),
                'lz': target_data_unscaled.tolist(),
                'run_id': row['run_id'],
                'weight': row['optimal_weight'],
                'opt': row['optimiser']
            }
            row_count += 1

    print(f"successfully loaded {row_count} fits into memory.")
    print(f"writing to {out_path}...")

    # --- the actual text compilation logic ---
    with open(out_path, 'w') as f:
        # gkeyll needs to know the highest atomic number present in the file
        max_z = max([ATOMIC_NUMBERS.get(sp, 0) for sp in memory_db.keys()])
        f.write(f"Maximum atomic number in file={max_z}, Number of elements={len(memory_db)}\n")
        
        # sort by atomic number so the file is clean and organized
        for species in sorted(memory_db.keys(), key=lambda x: ATOMIC_NUMBERS.get(x, 0)):
            charges = memory_db[species]
            z = ATOMIC_NUMBERS.get(species, 0)
            f.write(f"***Atomic number={z}; Number of charge states={len(charges)}\n")
            
            for charge in sorted(charges.keys()):
                densities = charges[charge]
                gkeyll_charge = int(charge) + 1
                f.write(f"**Charge state={gkeyll_charge}, Number of density intervals={len(densities)}\n")
                
                for dens_log10 in sorted(densities.keys()):
                    data = densities[dens_log10]
                    ne_m3 = dens_log10 + 6.0
                    p = data['params']
                    te_arr = data['te']
                    lz_arr = data['lz']
                    
                    f.write(f"{ne_m3:.6e} {p[0]:.6e} {p[1]:.6e} {p[2]:.6e} {p[3]:.6e} {p[4]:.6e} {len(te_arr)} ")
                    f.write(f"| run_id={data['run_id']} w={data['weight']} opt={data['opt']}\n")
                    f.write(" ".join([f"{te:.6e}" for te in te_arr]) + "\n")
                    f.write(" ".join([f"{lz:.6e}" for lz in lz_arr]) + "\n")

    print("compilation complete.")
    print("if you want to run simulations with this database")
    print("you put it in gkeyll/gyrokinetic/data/adas/ ")
    print("with the name: radiation_fit_parameters.txt ")
    print("and compile gkeyll.")
if __name__ == "__main__":
    main()

    