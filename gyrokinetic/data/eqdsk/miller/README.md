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
The supplied `q(R_omp)` is the **signed COCOS 3 safety factor**, finite,
nonzero, and of constant sign from the axis through the guard. Define
`F=Rmajor_axis*bmag_axis` and `s_psi=sign(F/q)`. The positive internal flux
magnitude is integrated using `abs(F)/abs(q)`; physical flux and its derivatives
carry `s_psi`. QPSI is written directly from the signed callback.

| bmag_axis | q | Toroidal field, from above | Poloidal field, R right and Z up | Outward psi |
|---|---|---|---|---|
| Positive | Positive | CCW | CCW | Increasing |
| Positive | Negative | CCW | CW | Decreasing |
| Negative | Positive | CW | CW | Decreasing |
| Negative | Negative | CW | CCW | Increasing |

The C and Python generators share this convention. To preserve an equilibrium
from the former magnitude-only interface, use
`q_new=sign(bmag_axis)*q_old_magnitude`. The supplied device examples retain
negative toroidal field and counterclockwise poloidal field by using negative q.
Reversing only q reverses the poloidal field; reversing only bmag_axis at fixed q
reverses the whole field; reversing both inputs reverses only the toroidal field.
Field magnitudes and geometry remain unchanged at fixed input magnitudes.

The application's binormal conversion uses positive `Cy=r0/abs(q0)` and keeps
`q0` signed. The tokamak helper instead uses `(R_omp-R_axis)/abs(q0)`.
Physical flux bounds must be sorted numerically; negative outward flux puts the
outer physical surface at the lower computational radial boundary. The Miller
angle stays geometrically CCW, while the computational parallel angle follows
B_pol and reverses when `s_psi<0`.

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
The equilibrium follows COCOS 3: `B_R=-psi_Z/R`, `B_Z=psi_R/R`, and `B_phi=F/R`;
psi is poloidal flux divided by `2*pi`. The header and JSON summary identify this convention.
No limiter is written.

Examples [miller_inp_d3d.py](miller_inp_d3d.py) and [miller_inp_tcv.py](miller_inp_tcv.py) provide
DIII-D and TCV shaping and q profiles. TCV includes a vertical offset and uses
the shifted `R_omp(r)` mapping to specify its coverage endpoints.
[miller_inp_ltx.py](miller_inp_ltx.py) replaces the old rational q profile (pole at
`R_omp` about 0.688328 m) with a negative signed exponential profile preserving LCFS
|q| and magnetic shear.

### Python API

[gen_miller_tcv.py](gen_miller_tcv.py) is a flat standalone script that defines
its own TCV parameters and q function, builds the equilibrium through the
Python API, and writes `tcv_miller.geqdsk` and
`tcv_miller.geqdsk.json` in the current directory. It also prints the requested
radial bounds in psi. Its plotting lines also save a diagnostic figure.
Run it from the repository root with:

```bash
python3 gyrokinetic/data/eqdsk/miller/gen_miller_tcv.py
```

A minimal circular example with both fields clockwise is:

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
signed vector-field reconstruction in all four orientations, mesh refinement,
domain coverage, and invalid or sign-changing inputs.
