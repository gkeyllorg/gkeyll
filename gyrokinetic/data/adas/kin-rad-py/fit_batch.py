# created and designed by jaime caballero for my phd at differ and the gkeyll team
# to mimic original https://github.com/jRoeltgen/radiation_operator/tree/main by J. Roeltgen
# code written with some help of ai 
#
# fit_batch.py
# this is the main production engine for the kin-rad-gkeyll repository.
# it iterates through elements, charge states, and densities, performing non-linear
# optimizations to find the parameters for roeltgen's radiation equation (eq 12).
# it safely stores the results in a master database and compiles a gkeyll-ready 
# text file at the end of the run.

import os
import sys
import csv
import argparse
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from datetime import datetime

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 11,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "legend.fontsize": 10
})

# tqdm for progress bar
try:
    from tqdm import tqdm
except ImportError:
    print("[error] tqdm is required for the batch progress bar. run: pip install tqdm")
    sys.exit(1)

# hook into engine scripts
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))
from data_parser import load_roeltgen_formatted
from optimizer_core import run_single_optimization, safe_integrand
from error_analysis import error_analysis


# -------------------------------------------------------------------------
# physics dictionary: the target density grids
# as discussed, roeltgen used specific density grids (usually 10^13 and 10^14 cm^-3).
# this dictionary dictates exactly what densities the batch script will fit.
# format -> species: { charge_state: [list of log10(n_e) in cm^-3] }
# -------------------------------------------------------------------------
# DENSITY_GRID = {
#     'H':  { 0: [13.0, 14.0] },
#     'He': { 0: [13.0, 14.0], 1: [13.0, 14.0] },
#     'Li': { 0: [13.0, 14.0], 1: [13.0, 14.0], 2: [13.0, 14.0] },
#     'Be': { 0: [13.0, 14.0], 1: [13.0, 14.0], 2: [13.0, 14.0], 3: [13.0, 14.0] },
#     # you can expand this for B, C, N, O, Ar as needed...
# }

# current density grid based on the actual densities used in roeltgen's original fits
# to extract a density grid from radiation_fit_parameters.txt
# you can use the get_densities.py script in this repository
# densities are in log10(cm^-3) 
DENSITY_GRID = {
    'H': { 0: [13.0, 15.3] },
    'He': { 0: [13.0], 1: [13.0] },
    'Li': { 0: [8.0, 8.7, 9.0, 10.0, 10.7, 11.0, 11.7, 12.0, 12.7, 13.0, 13.7, 14.0, 14.7, 15.0], 1: [8.0, 8.7, 9.0, 9.7, 10.0, 10.7, 11.0, 11.7, 12.0, 12.7, 13.0, 13.7], 2: [13.0] },
    'Be': { 0: [13.0], 1: [13.0], 2: [13.0], 3: [13.0] },
    'B': { 0: [13.0], 1: [13.0], 2: [13.0], 3: [13.0], 4: [13.0] },
    'C': { 0: [13.0], 1: [13.0], 2: [13.0], 3: [13.0], 4: [13.0], 5: [10.3, 10.7, 11.0, 11.3, 11.7, 12.0, 12.3, 12.7, 13.0, 13.3, 14.0, 14.3, 14.7, 15.0, 15.3] },
    'N': { 0: [13.0], 1: [13.0], 2: [13.0], 3: [13.0], 4: [13.0], 5: [13.0, 13.3, 14.0, 14.3, 14.7, 15.0, 15.3], 6: [12.3, 13.0] },
    'O': { 0: [13.0], 1: [13.0], 2: [13.0], 3: [13.0], 4: [13.0], 5: [13.0], 6: [7.7, 8.0, 8.3, 8.7, 9.0, 9.3, 9.7, 10.0, 10.3, 10.7, 11.0, 11.7, 12.0, 12.3, 12.7, 13.0, 13.3, 14.3, 14.7, 15.0, 15.3], 7: [11.0, 11.3, 11.7, 12.0, 12.3, 12.7, 13.0, 13.3, 14.0, 14.3, 14.7, 15.0, 15.3] },
    'F': { 0: [13.0], 1: [13.0], 2: [13.0], 3: [13.0], 4: [13.0], 5: [13.0], 6: [13.0], 7: [13.0], 8: [13.0] },
    'Ne': { 0: [13.0], 1: [13.0], 2: [13.0], 3: [13.0], 4: [13.0], 5: [13.0], 6: [13.0], 7: [13.0], 8: [8.0, 8.3, 8.7, 9.0, 9.3, 9.7, 10.0, 10.3, 10.7, 11.0, 11.3, 11.7, 12.0, 12.3, 12.7, 13.0, 13.3, 14.0, 14.3, 14.7, 15.0, 15.3], 9: [13.0] },
    'Al': { 0: [13.0], 1: [13.0], 2: [13.0], 3: [13.0], 4: [13.0], 5: [13.0], 6: [13.0], 7: [13.0], 8: [13.0], 9: [13.0], 10: [13.0] },
    'Ar': { 0: [13.0], 1: [13.0], 2: [13.0], 3: [13.0], 4: [13.0], 5: [13.0], 6: [13.0], 7: [13.0], 8: [13.0], 9: [13.0], 10: [13.0], 11: [13.0], 12: [13.0], 13: [13.0], 14: [13.0], 15: [13.0], 16: [13.0], 17: [13.0] },
}

# map species to atomic numbers for the gkeyll compiler
ATOMIC_NUMBERS = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Ar': 18}

def get_model_emissivity(params, Te_data):
    """
    evaluates the optimized integral over the T_e grid. 
    we split the integral at V_0 to handle the piecewise kink in the integrand.
    """
    A_scaled, alpha, beta, V0, gamma = params
    calcY = np.zeros_like(Te_data)
    
    for j, Te_j in enumerate(Te_data):
        # integrate up to a large velocity instead of infinity for numerical stability
        v_max = np.sqrt(40.0 * Te_j)

        if V0 < v_max:
            val_1, _ = quad(safe_integrand, 0, V0, args=(Te_j, A_scaled, alpha, beta, V0, gamma), epsabs=1e-8, epsrel=1e-8)
            val_2, _ = quad(safe_integrand, V0, v_max, args=(Te_j, A_scaled, alpha, beta, V0, gamma), epsabs=1e-8, epsrel=1e-8)
            calcY[j] = val_1 + val_2
        else:
            val, _ = quad(safe_integrand, 0, v_max, args=(Te_j, A_scaled, alpha, beta, V0, gamma), epsabs=1e-8, epsrel=1e-8)
            calcY[j] = val
            
    return calcY

def fetch_legacy_roeltgen_params(species, charge_state, density_log10_cm3):
    """
    parses the legacy gkeyll radiation_fit_parameters.txt to extract comparison params.
    looks exactly where the new filesystem architecture expects it to be.
    """
    filepath = os.path.join("fits_data", "roeltgen_data", "radiation_fit_parameters.txt")
    
    z = ATOMIC_NUMBERS.get(species.capitalize())
    if z is None: 
        return None
    
    # gkeyll file uses 1-based indexing for charge states
    target_charge = int(charge_state) + 1
    # gkeyll file uses m^-3 for density (e.g. 10^13 cm^-3 -> 10^19 m^-3)
    target_ne_m3 = density_log10_cm3 + 6.0 
    
    if not os.path.exists(filepath):
        # fail silently if the database isn't there, we just won't plot the red line
        print(f"[warning] legacy gkeyll parameters not found at {filepath}. skipping legacy comparison plot.")
        return None
        
    with open(filepath, 'r') as f:
        in_correct_z = False
        in_correct_c = False
        
        for line in f:
            if line.startswith('***Atomic number='):
                current_z = int(line.split(';')[0].split('=')[1])
                in_correct_z = (current_z == z)
            elif line.startswith('**Charge state=') and in_correct_z:
                current_c = int(line.split(',')[0].split('=')[1])
                in_correct_c = (current_c == target_charge)
            elif in_correct_z and in_correct_c and line[0].isdigit():
                parts = line.split()
                ne_val = float(parts[0])
                # check if this is the density block we want
                if np.isclose(ne_val, target_ne_m3, atol=0.1):
                    return [float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4]), float(parts[5])]
                    
    return None

def save_batch_plot(params_fit, Te_data, target_data_scaled, species, charge_state, density_log10, w, run_id, optimizer, min_te, Bs, plot_dir):
    """
    generates and saves the double log-log / linear-log plot. 
    uses the invisible legend trick to dodge the data curves.
    now features the legacy roeltgen fit for visual benchmarking.
    """
    calcY_scaled = get_model_emissivity(params_fit, Te_data)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    x_log = np.log10(Te_data)
    y_target = target_data_scaled / Bs
    y_model = calcY_scaled / Bs

    # format axes
    for ax in (ax1, ax2):
        ax.set_xlabel(r'$\log_{10}(T_e \ [\mathrm{eV}])$')
        ax.grid(True, alpha=0.3)

    ax1.set_ylabel(r'$\log_{10}(L_z \ [\mathrm{W \cdot m^3}])$')
    ax2.set_ylabel(r'$L_z \ [\mathrm{W \cdot m^3}]$')
    ax1.set_title('Log-Log Scale')
    ax2.set_title('Linear-Log Scale')

    # plot base lines (ax1 is log, ax2 is linear)
    ax1.plot(x_log, np.log10(y_target), '-b', label='OpenADAS')
    ax1.plot(x_log, np.log10(y_model), '*-g', label='Caballero fit')
    
    ax2.plot(x_log, y_target, '-b', label='OpenADAS')
    ax2.plot(x_log, y_model, '*-g', label='Caballero fit')

    # tighten the y-axes to remove empty white space at the top
    target_log_max = np.max(np.log10(y_target))
    target_log_min = np.min(np.log10(y_target))
    ax1.set_ylim(target_log_min - 2, target_log_max + 2)
    ax2.set_ylim(0, np.max(y_target) * 1.1)

    # invisible shape to trick the legend collision algorithm
    blank_handle = mpatches.Patch(color='none')

    # --- fetch and plot roeltgen legacy parameters ---
    roeltgen_params = fetch_legacy_roeltgen_params(species, charge_state, density_log10)

    if roeltgen_params is not None:
        # scale A back up for the integration physics, then down for plotting
        r_params_scaled = list(roeltgen_params)
        r_params_scaled[0] *= Bs 
        roeltgen_model = get_model_emissivity(r_params_scaled, Te_data) / Bs

        ax1.plot(x_log, np.log10(roeltgen_model), '*r', label='Roeltgen fit')
        ax2.plot(x_log, roeltgen_model, '*r', label='Roeltgen fit')

        r_text = (
            "Roeltgen Parameters:\n"
            rf'$A = {roeltgen_params[0]:.4e}$' + '\n'
            rf'$\alpha = {roeltgen_params[1]:.4f}$' + '\n'
            rf'$\beta = {roeltgen_params[2]:.4f}$' + '\n'
            rf'$V_0 = {roeltgen_params[3]:.4f}$' + '\n'
            rf'$\gamma = {roeltgen_params[4]:.4f}$'
        )
        
        # lock in lines, then deploy the auto-dodging text box on ax1
        leg1_main = ax1.legend(loc='lower right', fontsize=9)
        ax1.add_artist(leg1_main)
        ax1.legend([blank_handle], [r_text], loc='best', handlelength=0, handletextpad=0, 
                   fontsize=9, facecolor='white', edgecolor='gray', framealpha=0.8)
    else:
        # just a normal legend if no legacy data was found
        ax1.legend(loc='lower right', fontsize=9)

    # --- plot your parameters ---
    A_phys = params_fit[0] / Bs
    param_text = (
        f"Fitted Parameters:\n"
        rf'$A = {A_phys:.4e}$' + '\n'
        rf'$\alpha = {params_fit[1]:.4f}$' + '\n'
        rf'$\beta = {params_fit[2]:.4f}$' + '\n'
        rf'$V_0 = {params_fit[3]:.4f}$' + '\n'
        rf'$\gamma = {params_fit[4]:.4f}$'
    )
    
    # lock in lines, then deploy the auto-dodging text box on ax2
    leg2_main = ax2.legend(loc='lower right', fontsize=9)
    ax2.add_artist(leg2_main)
    ax2.legend([blank_handle], [param_text], loc='best', handlelength=0, handletextpad=0, 
               fontsize=9, facecolor='white', edgecolor='gray', framealpha=0.8)

    sup_title = rf'{species}$^{{{charge_state}+}}$ | $n_e = 10^{{{density_log10}}}\ \mathrm{{cm^{{-3}}}}$ | $w = {w:.2f}$ | $\min(T_e) = {min_te}\ \mathrm{{eV}}$ | Opt:{optimizer.upper()} | ID: {run_id}'
    plt.suptitle(sup_title, fontsize=14)

    # save compressed to keep the repo lightweight
    filename = os.path.join(plot_dir, f"{species}_{charge_state}_{density_log10}_{run_id}.png")
    plt.savefig(filename, dpi=100, bbox_inches='tight')
    plt.close()

def compile_gkeyll_database(memory_db, filepath):
    """
    compiles the safely-stored memory dictionary into the strict gkeyll .txt format.
    we do this at the very end to prevent file corruption if the batch crashes.
    """
    with open(filepath, 'w') as f:
        f.write(f"Maximum atomic number in file=18, Number of elements={len(memory_db)}\n")
        
        for species, charges in memory_db.items():
            z = ATOMIC_NUMBERS.get(species, 0)
            f.write(f"***Atomic number={z}; Number of charge states={len(charges)}\n")
            
            for charge, densities in charges.items():
                # gkeyll indexes charge states starting at 1
                gkeyll_charge = int(charge) + 1
                f.write(f"**Charge state={gkeyll_charge}, Number of density intervals={len(densities)}\n")
                
                for dens_log10, data in densities.items():
                    # convert density from cm^-3 to m^-3
                    ne_m3 = dens_log10 + 6.0
                    p = data['params']
                    te_arr = data['te']
                    lz_arr = data['lz']
                    
                    # main parameter line + custom metadata
                    f.write(f"{ne_m3:.6e} {p[0]:.6e} {p[1]:.6e} {p[2]:.6e} {p[3]:.6e} {p[4]:.6e} {len(te_arr)} ")
                    f.write(f"| run_id={data['run_id']} w={data['weight']} opt={data['opt']}\n")
                    
                    # original T_e and L_z arrays for gkeyll validation plotting
                    f.write(" ".join([f"{te:.6e}" for te in te_arr]) + "\n")
                    f.write(" ".join([f"{lz:.6e}" for lz in lz_arr]) + "\n")

# the function above are they associate variables (gkeyll_memory = {}) are no longer being used, 
# instead the function below
# compile_database_from_csv reads the master csv and at the end compiles the database
def compile_database_from_csv(csv_path, out_path):
    """
    reads the master csv (master db inside a folder fit, which actually is a .txt as well), 
    fetches original adas data arrays, and compiles
    the strict gkeyll-ready .txt file. this avoids the append trap by
    treating the csv as the single source of truth.
    """
    if not os.path.exists(csv_path):
        print(f"[error] could not find master database at {csv_path} for compilation.")
        return

    print("compiling gkeyll database from master csv...")
    memory_db = {}
    
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            species = row['element']
            charge = int(row['charge'])
            density = float(row['density'])
            min_te = float(row['min_te'])
            
            if species not in memory_db: memory_db[species] = {}
            if charge not in memory_db[species]: memory_db[species][charge] = {}
            
            try:
                Te_data, target_data_unscaled = load_roeltgen_formatted(species, charge, density, "formatted_data")
                if min_te > 0.0:
                    mask = Te_data >= min_te
                    Te_data = Te_data[mask]
                    target_data_unscaled = target_data_unscaled[mask]
            except Exception:
                continue

            memory_db[species][charge][density] = {
                'params': [float(row['A']), float(row['alpha']), float(row['beta']), float(row['V0']), float(row['gamma'])],
                'te': Te_data.tolist(),
                'lz': target_data_unscaled.tolist(),
                'run_id': row['run_id'],
                'weight': row['optimal_weight'],
                'opt': row['optimiser']
            }

    # text compilation logic
    with open(out_path, 'w') as f:
        max_z = max([ATOMIC_NUMBERS.get(sp, 0) for sp in memory_db.keys()])
        f.write(f"Maximum atomic number in file={max_z}, Number of elements={len(memory_db)}\n")
        
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

    print(f"compilation complete! gkeyll database updated at: {out_path}")

def main():
    parser = argparse.ArgumentParser(description="batch fit adas radiation parameters.")
    parser.add_argument("--id", type=str, required=True, help="unique name for this batch dataset (e.g., 'run1')")
    parser.add_argument("--elements", type=str, default="all", help="comma-separated elements (e.g., 'H,He,Li') or 'all'")
    parser.add_argument("--optimizer", type=str, default="slsqp", choices=['slsqp', 'ipopt', 'fmincon', 'trust-constr'])
    parser.add_argument("--min-te", type=float, default=1.5, help="chop off low-temp tails below this eV")
    parser.add_argument("--overwrite", action="store_true", help="force overwrite of already completed fits")
    args = parser.parse_args()

    # set up the architecture
    base_dir = os.path.join("fits_data", f"fit-db_{args.id}")
    plot_dir = os.path.join(base_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    
    master_csv_path = os.path.join(base_dir, f"master-fit_db_{args.id}.txt")
    gkeyll_db_path = os.path.join(base_dir, f"gkeyll-db_fit_{args.id}.txt")

    print("=========================================")
    print(f"--- BATCH RADIATION FITTER ---")
    print(f"--- DATASET ID: {args.id} ---")
    print(f"--- OPTIMIZER:  {args.optimizer.upper()} ---")
    print("=========================================\n")

    # read existing master db to figure out what to skip
    completed_fits = set()
    file_exists = os.path.isfile(master_csv_path)
    if file_exists and not args.overwrite:
        with open(master_csv_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # store signature: (species, charge, density)
                completed_fits.add((row['element'], int(row['charge']), float(row['density'])))

    # boot matlab engine if requested
    eng = None
    if args.optimizer == 'fmincon':
        print("booting matlab engine (~15s)...")
        import matlab.engine
        eng = matlab.engine.start_matlab()
        print("matlab connected.\n")

    # figure out target elements
    # figure out target elements
    if args.elements.lower() == "all":
        target_elements = list(DENSITY_GRID.keys())
    else:
        # split by comma, strip accidental spaces, and capitalize (e.g. 'h, he ' -> ['H', 'He'])
        raw_list = args.elements.split(',')
        target_elements = [el.strip().capitalize() for el in raw_list]
        
        # quick safety check to make sure the user didn't request an element not in our grid
        invalid_els = [el for el in target_elements if el not in DENSITY_GRID]
        if invalid_els:
            print(f"[error] the following elements are not in DENSITY_GRID: {invalid_els}")
            sys.exit(1)

    # calculate total iterations for the progress bar
    total_fits = sum(len(densities) for sp in target_elements for ch, densities in DENSITY_GRID.get(sp, {}).items())
    
    # memory dictionary to hold successful fits for the gkeyll compiler at the end
    gkeyll_memory = {}
    Bs = 1e30

    try:
        # progress bar setup
        with tqdm(total=total_fits, desc="fitting progress", unit="fit") as pbar:
            
            # --- THE BATCH LOOP ---
            for species in target_elements:

                # --- SAFETY CHECK ---
                # dynamically check if we have formatted data for this element
                import glob
                formatted_files = glob.glob(f"formatted_data/*_{species.lower()}_formatted.txt")
                if not formatted_files:
                    print(f"\n[warning] no formatted data found for {species}. did you run download_data.py? skipping element...")
                    continue
                # ------------------------

                gkeyll_memory[species] = {}
                
                for charge, densities in DENSITY_GRID[species].items():
                    gkeyll_memory[species][charge] = {}
                    
                    for density in densities:
                        signature = (species, charge, density)
                        
                        if signature in completed_fits and not args.overwrite:
                            pbar.set_postfix_str(f"skipped {species}^{charge}+ @ 10^{density}")
                            pbar.update(1)
                            continue
                            
                        pbar.set_postfix_str(f"loading {species}^{charge}+ @ 10^{density}...")
                        
                        try:
                            # load from the formatted data folder we made with download_data.py
                            Te_data, target_data_unscaled = load_roeltgen_formatted(
                                species=species, 
                                charge_state=charge, 
                                log10_ne_cm3=density,
                                data_dir="formatted_data"
                            )
                            
                            # chop physics tail if requested
                            if args.min_te > 0.0:
                                mask = Te_data >= args.min_te
                                Te_data = Te_data[mask]
                                target_data_unscaled = target_data_unscaled[mask]
                                
                        except Exception as e:
                            pbar.set_postfix_str(f"data err: {species}^{charge}+")
                            pbar.update(1)
                            continue
                            
                        target_data_scaled = target_data_unscaled * Bs
                        
                        # optimizer bounds and setup
                        # VERY IMPORTANT VARIABLES TO SET
                        initial_guess = [0.02, 8e3, 0.8, 1.5, -4.0] 
                        weight_powers = np.arange(0.1, 0.30, 0.05) 
                        
                        global_min_error = np.inf
                        best_fit_params = None
                        best_weight = None
                        best_successes = None
                        best_max_error = None
                        
                        # run weight sweep
                        for w in weight_powers:
                            current_guess = list(initial_guess)
                            passed_all_tests = False
                            
                            # the v_0 dodge loop
                            while current_guess[3] < 200 and not passed_all_tests:
                                pbar.set_postfix_str(f"fitting {species}^{charge}+ | w={w:.2f} | v0={current_guess[3]:.1f}")
                                
                                result = run_single_optimization(current_guess, Te_data, target_data_scaled, w, optimizer_choice=args.optimizer, eng=eng)
                                
                                if result.success:
                                    calcY_scaled = get_model_emissivity(result.x, Te_data)
                                    ratio = np.maximum(calcY_scaled / target_data_scaled, target_data_scaled / calcY_scaled)
                                    successes, max_error = error_analysis(ratio, Te_data, target_data_scaled)

                                    # 1. Base fractional errors (Perfect fit = 0.0)
                                    err_z1 = (max_error[0] - 1.0) / 0.2
                                    err_z2 = (max_error[1] - 1.0) / 0.4
                                    err_z3 = (max_error[2] - 1.0) / 1.0 if np.isfinite(max_error[2]) else 100.0
                                    
                                    # 2. Physical Importance Dials
                                    W1 = 10.0   # The Peak (Zone 1) is massively important
                                    W2 = 2.0    # The Slopes (Zone 2) are moderately important
                                    W3 = 0.5    # The Tails (Zone 3) are just a safety net
                                    
                                    # # # 3. Weighted Euclidean Distance
                                    fit_score = np.sqrt(W1 * (err_z1**2) + W2 * (err_z2**2) + W3 * (err_z3**2))

                                    # physics-weighted L2 metric (prioritizes the macro peak)
                                    # relative_error = np.abs(ratio - 1.0)
                                    # physical_weight = target_data_scaled / np.max(target_data_scaled)
                                    # fit_score = np.sum(relative_error * physical_weight)

                                    if fit_score < global_min_error:
                                        global_min_error = fit_score
                                        best_fit_params = result.x
                                        best_weight = w
                                        best_successes = successes
                                        best_max_error = max_error
                                        
                                    if all(successes):
                                        passed_all_tests = True
                                        initial_guess = result.x # chain guess
                                    else:
                                        current_guess[3] *= 1.25 # kick v_0
                                else:
                                    current_guess[3] *= 1.25

                        # --- logging and saving ---
                        if best_fit_params is not None:
                            run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
                            A_phys = best_fit_params[0] / Bs
                            pass_all = all(best_successes)
                            
                            # 1. save the plot securely
                            save_batch_plot(best_fit_params, Te_data, target_data_scaled, species, charge, density, best_weight, run_id, args.optimizer, args.min_te, Bs, plot_dir)
                            
                            # 2. append to master csv
                            file_exists_now = os.path.isfile(master_csv_path)
                            with open(master_csv_path, "a") as f:
                                if not file_exists_now:
                                    f.write("run_id,element,charge,density,min_te,optimiser,optimal_weight,A,alpha,beta,V0,gamma,max_err_z1,max_err_z2,max_err_z3,passed_all,is_manual_override\n")
                                
                                row = [
                                    run_id, species, str(charge), str(density), str(args.min_te), args.optimizer, 
                                    f"{best_weight:.2f}", f"{A_phys:e}", f"{best_fit_params[1]:.4f}", 
                                    f"{best_fit_params[2]:.4f}", f"{best_fit_params[3]:.4f}", f"{best_fit_params[4]:.4f}",
                                    f"{best_max_error[0]:.4f}", f"{best_max_error[1]:.4f}", f"{best_max_error[2]:.4f}", 
                                    str(pass_all), "False"
                                ]
                                f.write(",".join(row) + "\n")
                                
                            # 3. stash arrays in memory for the gkeyll compiler
                            gkeyll_memory[species][charge][density] = {
                                'params': [A_phys, best_fit_params[1], best_fit_params[2], best_fit_params[3], best_fit_params[4]],
                                'te': Te_data.tolist(),
                                'lz': (target_data_scaled / Bs).tolist(), # store unscaled Lz
                                'run_id': run_id,
                                'weight': f"{best_weight:.2f}",
                                'opt': args.optimizer
                            }
                        
                        pbar.update(1)

        # build the final gkeyll text file using all successful fits
        # print("\ncompiling gkeyll database...")
        # compile_gkeyll_database(gkeyll_memory, gkeyll_db_path)
        # print(f"batch complete. data stored in {base_dir}")

        # build the final gkeyll text file using the csv truth
        print(f"\nbatch loop complete. triggering auto-compiler...")
        compile_database_from_csv(master_csv_path, gkeyll_db_path)



    finally:
        # in case matlab optimizer is being used (no longer needed)
        if eng is not None:
            eng.quit()

if __name__ == "__main__":
    main()
