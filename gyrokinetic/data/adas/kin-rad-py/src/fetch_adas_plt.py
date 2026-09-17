# created and designed by jaime caballero for my phd at differ and the gkeyll team
# to mimic original https://github.com/jRoeltgen/radiation_operator/tree/main by J. Roeltgen
# code written with some help of ai 
#
# fetch_adas_plt.py
# this script is responsible for fetching the raw ADAS plt (radiation) data files.

import os
import urllib.request
import urllib.error
import argparse

# Dictionary mapped exactly from Gkeyll's process_adas.py, but strictly for 'plt' (radiation) files
ADAS_PLT_FILES = {
    "H":  "plt12_h.dat",
    "He": "plt96_he.dat",
    "Li": "plt96_li.dat",
    "Be": "plt96_be.dat",
    "B":  "plt89_b.dat",
    "C":  "plt96_c.dat",
    "N":  "plt96_n.dat",
    "O":  "plt96_o.dat",
    "Ar": "plt89_ar.dat",
    "Ne": "plt96_ne.dat",
    "W":  "plt50_w.dat",
    "F":  "plt89_f.dat",
}

def fetch_plt_file(element, data_dir="data_raw"):
    """
    Downloads the ADAS plt file for a given element, mimicking Gkeyll's fetch logic.
    """
    # Sanitize input (e.g., 'he' -> 'He')
    element = element.capitalize()
    
    if element not in ADAS_PLT_FILES:
        print(f"Error: Element '{element}' is not supported.")
        print(f"Supported elements: {list(ADAS_PLT_FILES.keys())}")
        return None

    filename = ADAS_PLT_FILES[element]
    
    # ADAS URL structure: base_url / type_folder / filename
    # E.g., .../adf11/plt96/plt96_he.dat
    folder_name = filename.split('_')[0] 
    base_url = "https://open.adas.ac.uk/download/adf11"
    full_url = f"{base_url}/{folder_name}/{filename}"

    # Ensure the data directory exists
    os.makedirs(data_dir, exist_ok=True)
    loc = os.path.join(data_dir, filename)

    # Check if we already have it to avoid spamming the ADAS servers
    if os.path.exists(loc):
        print(f"  -> File {filename} already exists in '{data_dir}/'. Skipping download.")
        return loc

    print(f"  -> Downloading {filename} from {full_url} ...")
    try:
        with urllib.request.urlopen(full_url) as response:
            if response.status != 200:
                raise urllib.error.HTTPError(full_url, response.status, "Failed to download", response.headers, None)
            content_bytes = response.read()

        # Sanity check: ADAS sometimes returns an empty or tiny error page instead of a 404
        if len(content_bytes) < 1000:
            raise ValueError(f"Could not fetch a valid file for {filename} from ADAS! Response was too short.")

        with open(loc, "wb") as f:
            f.write(content_bytes)
            
        print(f"  -> Successfully saved to {loc}")
        return loc
        
    except Exception as e:
        print(f"  -> Error downloading {filename}: {e}")
        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch ADAS plt (radiation) data files.")
    parser.add_argument("--element", type=str, default="all", 
                        help="Element symbol (e.g., H, He, Li) or 'all' to download everything.")
    parser.add_argument("--dir", type=str, default="data", 
                        help="Target directory to save the files (default is 'data/')")
    
    args = parser.parse_args()

    print("=========================================")
    print("--- ADAS RADIATION (PLT) DATA FETCHER ---")
    print("=========================================\n")

    if args.element.lower() == "all":
        for el in ADAS_PLT_FILES.keys():
            fetch_plt_file(el, args.dir)
    else:
        fetch_plt_file(args.element, args.dir)
        
    print("\nFetch operation complete.")