# created and designed by jaime caballero for my phd at differ and the gkeyll team
# to mimic original https://github.com/jRoeltgen/radiation_operator/tree/main by J. Roeltgen
# code written with some help of ai 
#
# replot_latex.py
# a temporary script to read a master database csv, rapidly recalculate the 
# integration curves, and generate publication-quality plots using latex fonts.

import os
import sys
import csv
import argparse
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# point to src for our data tools
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))
from data_parser import load_roeltgen_formatted
from optimizer_core import safe_integrand

# --- MATPLOTLIB LATEX INJECTION ---
# this forces matplotlib to use the 'computer modern' font (the standard latex font)
# and enables math formatting for all text elements.
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 11,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "legend.fontsize": 10
})

def get_model_emissivity(params, Te_data):
    """re-evaluates the integral using the parameters from the csv."""
    A_scaled, alpha, beta, V0, gamma = params
    calcY = np.zeros_like(Te_data)
    for j, Te_j in enumerate(Te_data):
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
    """grabs the legacy data for the red baseline curve."""
    filepath = os.path.join("fits_data", "roeltgen_data", "radiation_fit_parameters.txt")
    atomic_numbers = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'Ne': 10,  'Ar': 18}
    z = atomic_numbers.get(species.capitalize())
    if not z or not os.path.exists(filepath): return None
    
    target_charge = int(charge_state) + 1
    target_ne_m3 = density_log10_cm3 + 6.0 
    
    with open(filepath, 'r') as f:
        in_z, in_c = False, False
        for line in f:
            if line.startswith('***Atomic number='):
                in_z = (int(line.split(';')[0].split('=')[1]) == z)
            elif line.startswith('**Charge state=') and in_z:
                in_c = (int(line.split(',')[0].split('=')[1]) == target_charge)
            elif in_z and in_c and line[0].isdigit():
                parts = line.split()
                if np.isclose(float(parts[0]), target_ne_m3, atol=0.1):
                    return [float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4]), float(parts[5])]
    return None

def main():
    # SET YOUR DATASET ID HERE
    dataset_id = "prod_rs_server_1"  # e.g., "PROD_2026_05_26"
    parser = argparse.ArgumentParser(description="re-plot the fitted curves with latex fonts for publication.")
    parser.add_argument("--id", type=str, default=dataset_id, help="dataset ID to plot (e.g., 'PROD_RUN')")
    args = parser.parse_args()
    dataset_id = args.id
    
    base_dir = os.path.join("fits_data", f"fit-db_{dataset_id}")
    csv_path = os.path.join(base_dir, f"master-fit_db_{dataset_id}.txt")
    
    # create a new folder for the beautiful plots
    out_dir = os.path.join(base_dir, "plots")
    os.makedirs(out_dir, exist_ok=True)
    
    if not os.path.exists(csv_path):
        print(f"[error] could not find {csv_path}")
        return

    print(f"reading master database: {csv_path}")
    print(f"saving latex plots to: {out_dir}\n")

    Bs = 1e30

    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            species = row['element']
            charge = int(row['charge'])
            density = float(row['density'])
            min_te = float(row['min_te'])
            run_id = row['run_id']
            w = float(row['optimal_weight'])
            
            print(f"re-plotting {species}^{charge}+ @ 10^{density}...")

            # rebuild params array
            params = [float(row['A']) * Bs, float(row['alpha']), float(row['beta']), float(row['V0']), float(row['gamma'])]
            
            # load raw data
            try:
                Te_data, target_data_unscaled = load_roeltgen_formatted(species, charge, density, "formatted_data")
            except Exception:
                print(f"  -> missing raw data for {species}, skipping.")
                continue
                
            if min_te > 0.0:
                mask = Te_data >= min_te
                Te_data = Te_data[mask]
                target_data_unscaled = target_data_unscaled[mask]
                
            target_data_scaled = target_data_unscaled * Bs
            calcY_scaled = get_model_emissivity(params, Te_data)
            
            # --- PLOTTING ---
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
            x_log = np.log10(Te_data)
            y_target = target_data_scaled / Bs
            y_model = calcY_scaled / Bs

            for ax in (ax1, ax2):
                # use latex math rendering for labels
                ax.set_xlabel(r'$\log_{10}(T_e \ [\mathrm{eV}])$')
                ax.grid(True, alpha=0.3)
                
            ax1.set_ylabel(r'$\log_{10}(L_z \ [\mathrm{W \cdot m^3}])$')
            ax2.set_ylabel(r'$L_z \ [\mathrm{W \cdot m^3}]$')
            
            # titles
            ax1.set_title('Log-Log Scale')
            ax2.set_title('Linear-Log Scale')

            ax1.plot(x_log, np.log10(y_target), '-b', label=r'OpenADAS')
            ax1.plot(x_log, np.log10(y_model), '*-g', label=r'Caballero Fit')
            ax2.plot(x_log, y_target, '-b', label=r'OpenADAS')
            ax2.plot(x_log, y_model, '*-g', label=r'Caballero Fit')

            target_log_max = np.max(np.log10(y_target))
            target_log_min = np.min(np.log10(y_target))
            ax1.set_ylim(target_log_min - 2, target_log_max + 2)
            ax2.set_ylim(0, np.max(y_target) * 1.1)

            blank_handle = mpatches.Patch(color='none')
            
            # roeltgen baseline
            r_params = fetch_legacy_roeltgen_params(species, charge, density)
            if r_params is not None:
                r_params_scaled = list(r_params)
                r_params_scaled[0] *= Bs 
                r_model = get_model_emissivity(r_params_scaled, Te_data) / Bs

                ax1.plot(x_log, np.log10(r_model), '*r', label=r'Roeltgen fit')
                ax2.plot(x_log, r_model, '*r', label=r'Roeltgen fit')

                r_text = (
                    'Roeltgen Parameters:' + '\n'
                    rf'$A = {r_params[0]:.4e}$' + '\n'
                    rf'$\alpha = {r_params[1]:.4f}$' + '\n'
                    rf'$\beta = {r_params[2]:.4f}$' + '\n'
                    rf'$V_0 = {r_params[3]:.4f}$' + '\n'
                    rf'$\gamma = {r_params[4]:.4f}$'
                )
                leg1_main = ax1.legend(loc='lower right', fontsize=9)
                ax1.add_artist(leg1_main)
                ax1.legend([blank_handle], [r_text], loc='best', handlelength=0, handletextpad=0, 
                           fontsize=10, facecolor='white', edgecolor='gray', framealpha=0.8)
            else:
                ax1.legend(loc='lower right', fontsize=9)

            # my fit parameters
            my_text = (
                'Fitted Parameters:' + '\n'
                rf'$A = {params[0]/Bs:.4e}$' + '\n'
                rf'$\alpha = {params[1]:.4f}$' + '\n'
                rf'$\beta = {params[2]:.4f}$' + '\n'
                rf'$V_0 = {params[3]:.4f}$' + '\n'
                rf'$\gamma = {params[4]:.4f}$'
            )
            
            leg2_main = ax2.legend(loc='lower right', fontsize=9)
            ax2.add_artist(leg2_main)
            ax2.legend([blank_handle], [my_text], loc='best', handlelength=0, handletextpad=0, 
                       fontsize=10, facecolor='white', edgecolor='gray', framealpha=0.8)

            # majestic main title
            sup_title = rf'{species}$^{{{charge}+}}$ | $n_e = 10^{{{density}}}\ \mathrm{{cm^{{-3}}}}$ | $w = {w:.2f}$ | $\min(T_e) = {min_te}\ \mathrm{{eV}}$ | ID: {run_id}'
            plt.suptitle(sup_title, fontsize=14)
            
            filename = os.path.join(out_dir, f"{species}_{charge}_{density}_{run_id}.png")
            plt.savefig(filename, dpi=120, bbox_inches='tight')
            plt.close()

    print("done. your plots are ready for the paper.")

if __name__ == "__main__":
    main()