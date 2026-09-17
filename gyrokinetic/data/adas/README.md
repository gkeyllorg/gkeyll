giADAS data
---------

ADAS data is not included in gkylzero and must be downloaded by user.

To download ADAS data and convert to format accessible by gkylzero,
run machines/mkdeps.<machine_name>.sh with flag --download-adas=yes

When running a simulation, gkeyll reads radiation data from `radiation_fit_parameters.txt`. The radiation fits were originally are produced with [code maintained by J. Roeltgen](https://github.com/jRoeltgen/radiation_operator.git), requiring Matlab. The radiation fits from this implementation are in `radiation_fit_parameters_legacy.txt`

Now, we ported the matlab implemenation to python and the fits can be produced natively. The sub-repository to produce and adjust fits is present in the folder `kin-rad-py`, as well as intructions on how to use it. The current `radiation_fit_parameters.txt` file is produced with this new implementation. If fits are adjusted, remember to have the new `radiation_fit_parameters.txt` in this location, then compile gkeyll. This will copy the data to `${PREFIX:=$HOME/gkylsoft}/gkeyll/share/adas/radiation_fit_parameters.txt` and this effectevely the file use at runtime. 


