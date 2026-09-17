# kin-rad-gkeyll
**automated kinetic radiation parameter fitting for gkeyll.**

This repository houses the complete data ingestion and non-linear optimization pipeline for calculating plasma radiation parameters according to the [original repository](https://github.com/jRoeltgen/radiation_operator) written on matlab following the methodology outlined by [J. Roeltgen et al. NF 65 2025](https://iopscience.iop.org/article/10.1088/1741-4326/adff28): "A kinetic line-driven radiation operator and its application to Gyrokinetics". It takes raw bremsstrahlung and line radiation data from openadas, fits it to the analytical equations from the mentioned paper, and compiles the results into a lightweight database ready for the gkeyll kinetic simulation framework.

---

### Filesystem architecture

The repository is separated data into ingestion, batch processing, testing and fine-tuning environments:

```text
kin-rad-gkeyll/
├── raw_data/                  # raw openadas .dat files (git-ignored)
├── formatted_data/            # roeltgen-formatted 4-column .txt files (git-ignored)
├── fits_data/                 # the master output directory
│   ├── roeltgen_data/         # legacy database from matlab implementation 
│   ├── fit-db_<ID>/           # output folder for a specific batch run
│   │   ├── master-fit_db_<ID>.txt   # csv database of all fits and metadata
│   │   ├── gkeyll-db_fit_<ID>.txt   # the compiled gkeyll-ready text file
│   │   └── plots/             # dual log-log/linear-log pngs of every fit
│   └── individual_fits/       # directory for manually tuned fits (TBD)
├── src/                       # core python modules
│   ├── compile_gk_db.py       # produces gkeyll-ready .txt from a database id 
│   ├── error_analysis.py      # handles zone-based fractional error physics
│   ├── optimizer_core.py      # handles the scipy/matlab integration wrapper
│   ├── data_parser.py         # gets a function Li(ne, Te) from adas data
│   ├── fetch_adas_plt.py      # helper script to fetch raw ADAS plt files from openadas
│   ├── format_adas.py         # formats raw .dat adas to formatted .txt
│   ├── get_densities.py       # helper scrip to extract ne points use in original roeltgen fits
│   ├── opt.py                 # helper scripts to run fmincon matlab (no longer needed)
│   ├── optimize.m             # helper scripts to run fmincon matlab (no longer needed)
│   └── update-plots.py        # re-produce the plots for a given database id 
├── test/                      
│   └── check_fit.py           # checks fit according to ID (TBD)
├── download_data.py           # downloads and formats raw ADAS data 
├── fit_batch.py               # fits a batch of radiation parameters for given ID 
│                                and elements/charge/densities
├── fit_single.py              # fine-tune a single fit given an ID
├── swap_fits.py               # tool to overwrite batch fits with manual ones (TBD)
├── .gitignore
├── requirements.txt
├── environment.yml
└── README.md
```

### Usage

Producing fit parameters for some elements is a straightforward two-step process: fetch the data, then run the batch optimizer.

#### Data ingestion (`download_data.py`)
Automatically fetches raw OpenAdas radiation data and translates it into our strict formatted `.txt` grids.

```bash
# download and format all supported elements
python download_data.py
# or just target a specific element
python download_data.py --elements He
```

#### Batch optimization 
fit_batch.py it sweeps through the $N_e$ grids, fits the radiation curves to roeltgen's analytical models, and dynamically saves the best parameters and plots to in the fits_data/ directory.

```bash
# fit all available elements and assign a unique dataset id saved in fits_data/fit-db_<ID>/
python fit_batch.py --id run_1 --elements all --min-te 1.5
```

```bash
# add elements to an existing database 
python fit_batch.py --id run_1 --elements Al, He 
```
> *note: in case you want to overwrite the fits of an existing database, check the flag `--overwrite`*

> *note: the current database included in this repository `prod_rs_server_1` is the current radiation parameters used in gkeyll.*

#### Compile gkeyll database and how to use it at runtime
After the batch optimization is done, you can compile the results into a gkeyll-ready text file with:
```bash
python src/compile_gk_db.py --id prod_rs_server_1
```
> *note: This is already done during `fit_batch.py`, but is just a helper script to re generate the database if needed. IMPORTANT: to actually use the generate data in gkeyll, remember to copy the generated `gkeyll-db_fit_<ID>.txt` file into `gkeyll/gyrokinetic/data/adas/` with the specific name: `radiation_fit_parameters.txt` and re-compile gkeyll. This will copy that file to `${PREFIX:=$HOME/gkylsoft}/gkeyll/share/adas/radiation_fit_parameters.txt` and this will be the file used at runtime.*

#### Manual tuning  (TBD)

If ever the automated fits need to be fine-tuned, the next step of this repository is to include a manual fitting environment in which a specific fit can be dialed "by hand" and store in `fits_data/individual_fits/` for later inserting into a database (e.g. `prod_rs_server_1`). The scripts for this environment are still in development, but the general idea is to have a set of tools that allow you to:

- `fit_single.py`: manually dial in weights and $V_0$ bounds for a highly specific $N_e$ / $T_e$ slice.
- `swap_fits.py`: replace a bad batch fit with your handcrafted fit inside the master database.
- `test/check_fit.py`: visual plots for any fit ID directly from the terminal.

> *note: created and designed by jaime caballero for my phd at DIFFER and the Gkeyll team to mimic original https://github.com/jRoeltgen/radiation_operator/tree/main by J. Roeltgen; the scripts in this repository were created with the help of ai.*

