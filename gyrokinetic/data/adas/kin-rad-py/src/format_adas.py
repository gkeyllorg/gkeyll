# created and designed by jaime caballero for my phd at differ and the gkeyll team
# to mimic original https://github.com/jRoeltgen/radiation_operator/tree/main by J. Roeltgen
# code written with some help of ai 
#
# format_adas.py
# this script contains the logic to parse the raw ADAS ADF11 .dat files and
# convert them into the strict 4-column format defined by Roeltgen. 

import os
import glob
import argparse

def format_adas_to_roeltgen(input_filepath, output_dir="data_formatted"):
    """
    Parses a raw ADAS ADF11 (.dat) file and converts it into 
    Roeltgen's strict 4-column formatted text file.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    filename = os.path.basename(input_filepath)
    name_no_ext = os.path.splitext(filename)[0]
    out_filename = f"{name_no_ext}_formatted.txt"
    out_filepath = os.path.join(output_dir, out_filename)
    
    with open(input_filepath, 'r') as f:
        lines = f.readlines()
        
    # 1. Read the header dimensions
    header = lines[0].split()
    n_ion = int(header[0])
    n_ne = int(header[1])
    n_T = int(header[2])
    
    # Skip the "---" separator line
    current_line = 2
    
    # Check for metastable resolved files (starts with numbers instead of floats)
    line = lines[current_line]
    if all(a.isdigit() for a in line.split()):
        current_line += 2 # Skip metastable info
        
    # 2. Extract log10(Density) array
    logNe = []
    while len(logNe) < n_ne:
        logNe.extend([float(x) for x in lines[current_line].split()])
        current_line += 1
        
    # 3. Extract log10(Temperature) array
    logT = []
    while len(logT) < n_T:
        logT.extend([float(x) for x in lines[current_line].split()])
        current_line += 1
        
    # 4. Prepare output buffer with exact header
    out_lines = [" Charge State log10(Te (eV)) log10(Electron Density (cm-3)) log10(Coefficient) \n"]
    
    # 5. Loop through each charge state block
    charge_state = 1
    while charge_state <= n_ion:
        # Skip the subheader line (e.g., --------/ IGRD= 1 / IPRT= 0 /--------/ Z1= 1...)
        current_line += 1
        
        # Read the entire matrix of coefficients for this charge state
        block_data = []
        while len(block_data) < n_ne * n_T:
            block_data.extend([float(x) for x in lines[current_line].split()])
            current_line += 1
            
        # Write combinations (Te outer loop, Ne inner loop)
        idx = 0
        for t in logT:
            for n in logNe:
                val = block_data[idx]
                # Format to match Roeltgen's precise scientific notation spacing
                out_lines.append(f"{charge_state:12d} {t:21.17f} {n:25.17f} {val:25.17f}\n")
                idx += 1
        
        charge_state += 1
        
    # 6. Save the perfectly formatted file
    with open(out_filepath, 'w') as f_out:
        f_out.writelines(out_lines)
        
    print(f"  -> Converted {filename} to {out_filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format raw ADAS .dat files into Roeltgen style .txt files.")
    parser.add_argument("--element", type=str, default="all", 
                        help="Element symbol (e.g., H, He, Li) or 'all' to format everything.")
    parser.add_argument("--input_dir", type=str, default="data", 
                        help="Directory containing the raw .dat files")
    parser.add_argument("--output_dir", type=str, default="data_formatted", 
                        help="Directory to save the formatted .txt files")
    
    args = parser.parse_args()
    
    print("=========================================")
    print("--- ADAS DATA FORMATTER ---")
    print("=========================================\n")
    
    # Decide what to search for based on the element argument
    if args.element.lower() == "all":
        search_pattern = os.path.join(args.input_dir, "*.dat")
    else:
        # e.g., looks for anything ending in _he.dat
        element_lower = args.element.lower()
        search_pattern = os.path.join(args.input_dir, f"*_{element_lower}.dat")
        
    raw_files = glob.glob(search_pattern)
    
    if not raw_files:
        print(f"No matching files found for '{args.element}' in {args.input_dir}/.")
        print("Please check your spelling or run fetch_adas_plt.py first.")
    else:
        for file in raw_files:
            format_adas_to_roeltgen(file, args.output_dir)
            
    print("\nFormatting complete.")