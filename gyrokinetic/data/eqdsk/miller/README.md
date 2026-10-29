# Miller GEQDSK generation

## Miller generator

[gkeyll_gen_miller_eqdsk.py](gkeyll_gen_miller_eqdsk.py) is both an importable module and a command-line
generator. It requires NumPy and SciPy. Matplotlib is optional, for `--plot`.

From the repository root:

```bash
python3 gyrokinetic/data/eqdsk/miller/gkeyll_gen_miller_eqdsk.py \
  gyrokinetic/data/eqdsk/miller/miller_inp_tcv.py --output tcv_miller.geqdsk \
  --nr 257 --nz 257 --plot tcv_miller.png
```

The Python input defines `geometry={...}` and a scalar `q(R_omp)` function. All
physical parameters belong in that file; the CLI never evaluates a function
expression. `R_omp` is outboard midplane MAJOR radius in meters. Required geometry
keys are `Rmajor_axis`, `bmag_axis` (toroidal field), and `Rmajor_LCFS` (outboard LCFS major radius).
Defaults: `Zvert_axis=0`, `kappa=1`, `delta=0`, `a_shift=0`. Lengths are meters and fields
tesla. `bmag_axis` is the signed toroidal field: negative means negative cylindrical
B_phi, positive means positive B_phi; zero and nonfinite values are rejected.
The supplied `q(R_omp)` is a positive finite magnitude. Flux quadrature uses
`abs(Rmajor_axis*bmag_axis)`, so reversing the field preserves outward-increasing psi,
the poloidal field, and radial conversions. The GEQDSK central field, F profile,
and q profile carry the toroidal sign; the input q and diagnostic q(R_omp) plot
remain magnitudes. Existing positive-field inputs retain their interpretation;
the supplied examples explicitly choose negative toroidal field.

Optional `Rmajor_omp_bounds=(lower,upper)` requests coverage of complete flux surfaces
for a core annulus, SOL-only simulation, or LCFS-crossing simulation. These are
NOT rectangular R bounds. The equilibrium always includes the axis and LCFS;
q must be valid all the way from the axis through the physical guard surface.
Without coverage bounds the LCFS is the outer certified surface. `--rz-bounds RMIN RMAX ZMIN ZMAX` explicitly sets the rectangular grid but does not extend
the certified physical surfaces. Rectangular corners use the documented smooth
flux continuation beyond the guard, identified in the JSON summary and plots.
The default 5-percent guard needs at least three R-Z cells of clearance; coarse
grids may require `--guard-fraction 0.25`. The node counts must be odd, 5..999.

The output includes `<output>.json`, with parameters, input provenance, flux
values, physical coverage, a conversion table, and numerical residuals. The
GEQDSK has constant `F=Rmajor_axis*bmag_axis` and zero pressure/current placeholders. It is
a synthetic magnetic geometry, not a Grad-Shafranov solution or a wall model.
No limiter is written.

Examples [miller_inp_d3d.py](miller_inp_d3d.py) and [miller_inp_tcv.py](miller_inp_tcv.py) provide
DIII-D and TCV shaping and q profiles. TCV includes a vertical offset and uses
the shifted `R_omp(r)` mapping to specify its coverage endpoints.
[miller_inp_ltx.py](miller_inp_ltx.py) replaces the old rational q profile (pole at
`R_omp` about 0.688328 m) with a positive exponential profile preserving LCFS
q and shear.

### Python API

[gen_miller_tcv.py](gen_miller_tcv.py) is a flat standalone script that defines
its own TCV parameters and q function, builds the equilibrium through the
Python API, and writes `tcv_miller.geqdsk` and
`tcv_miller.geqdsk.json` in the current directory. It also prints the requested
radial bounds in psi. Uncomment its plotting lines to save a diagnostic figure.
Run it from the repository root with:

```bash
python3 gyrokinetic/data/eqdsk/miller/gen_miller_tcv.py
```

A minimal circular example is:

```python
from gkeyll_gen_miller_eqdsk import (
  MillerGeometry, MillerEquilibrium, sample_eqdsk, write_eqdsk)
geometry = MillerGeometry(Rmajor_axis=2, bmag_axis=-1.5, Rmajor_LCFS=2.5)
equilibrium = MillerEquilibrium(geometry, lambda R_omp: 2)
data = sample_eqdsk(equilibrium)
write_eqdsk("circle.geqdsk", equilibrium, data)
print(equilibrium.fs_coords_at_R_omp(2.4))  # R_omp, Miller r, psi, rho_psi, q.
```

### Verification

```bash
python3 -m unittest discover -s gyrokinetic/data/eqdsk/miller -v
```

The Python tests check analytic flux, coordinate conversions, GEQDSK formatting,
field reconstruction, mesh refinement, domain coverage, and invalid inputs.
