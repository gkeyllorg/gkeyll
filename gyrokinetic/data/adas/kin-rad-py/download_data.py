# created and designed by jaime caballero for my phd at differ and the gkeyll team
# to mimic original https://github.com/jRoeltgen/radiation_operator/tree/main by J. Roeltgen
# code written with some help of ai 
#
# download_data.py
# this script is the entry point for the adas radiation data ingestion pipeline.
# it reaches out to openadas, grabs the raw bremsstrahlung and line radiation 
# files (.dat), and translates them into the strict 4-column layout defined 
# by roeltgen. this gives us the clean L_z(T_e, n_e) to do non-linear fitting.


import os
import sys
import glob
import argparse

# point python to the src directory so we can import our engine scripts
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

# import the tools we built previously
from fetch_adas_plt import fetch_plt_file, ADAS_PLT_FILES
from format_adas import format_adas_to_roeltgen

def main():
    parser = argparse.ArgumentParser(description="fetch and format adas radiation data.")
    parser.add_argument("--elements", type=str, default="all", 
                        help="element symbol (e.g., H, He, Li) or 'all' to ingest everything.")
    
    args = parser.parse_args()

    # define our pipeline directories
    raw_dir = "raw_data"
    fmt_dir = "formatted_data"

    print("=========================================")
    print("--- ADAS DATA INGESTION PIPELINE ---")
    print("=========================================\n")

    # figure out which elements we are processing
    if args.elements.lower() == "all":
        elements_to_process = list(ADAS_PLT_FILES.keys())
    else:
        element_cap = args.elements.capitalize()
        if element_cap not in ADAS_PLT_FILES:
            print(f"[error] element '{args.elements}' is not supported.")
            print(f"supported elements: {list(ADAS_PLT_FILES.keys())}")
            sys.exit(1)
        elements_to_process = [element_cap]

    # phase 1: fetch the raw data from openadas
    print(f"\n--- phase 1: fetching raw data ---")
    downloaded_files = []
    for el in elements_to_process:
        # fetch_plt_file handles checking if it already exists to save bandwidth
        filepath = fetch_plt_file(el, data_dir=raw_dir)
        if filepath:
            downloaded_files.append(filepath)

    # phase 2: format the raw data for the optimizer
    print(f"\n--- phase 2: formatting data to roeltgen standard ---")
    if not downloaded_files:
        print("no files were downloaded or found. aborting format phase.")
        sys.exit(1)

    for raw_file in downloaded_files:
        # format_adas_to_roeltgen parses the adf11 grid and outputs the T_e, n_e, and L_z arrays
        format_adas_to_roeltgen(raw_file, output_dir=fmt_dir)

    print("\n=========================================")
    print("ingestion complete. ready for batch fitting.")
    print("=========================================")

if __name__ == "__main__":
    main()