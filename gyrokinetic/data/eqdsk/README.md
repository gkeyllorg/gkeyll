# Sample GEQDSK files

See [Miller GEQDSK generation](miller/README.md) for the Miller generator,
example inputs, numerical documentation, and verification commands.

## Sample equilibria

Put the GEQDSK file you want in this folder and name it `input.geqdsk` to run g0 unit tests with it.

### File descriptions

- `asdex.geqdsk`: Similar to what was used in D. Michels, et al. Phys. Plasmas
  29, 032307 (2022).
- `cerfon.geqdsk`: Generated from analytical equilibrium in Antoine J. Cerfon and Jeffrey P. Freidberg, PHYSICS OF PLASMAS 17, 032502 (2010)
- `elliptical.geqdsk`: Fake equilibriuum generated using elliptical flux surfaces
- `nstxu_DN.geqdsk`: Double-null equilibrium (with lowest drsep value) from NSTX-U test run in March 2016 . Original file name on PPPL Flux cluster: g203655.00304 (shot number = 203655, time-stamp = 304 ms).
- `nstxu_SN.geqdsk`: Single-null equilibrium from NSTX-U test run in April 2016. Original file name on PPPL Flux cluster: g204700.00238 (shot number = 204700, time-stamp = 238 ms).
- `solovev.geqdsk`: Generated from analytical equilibrium described in Noah Mandell's thesis Sec. 5.3
- `step.geqdsk`: From STEP. Proposed magnetic equilibrium.
- `wham.geqdsk`: case with beams and ~0.8 T central field.
- `tcv.geqdsk`:  TCV equilibrium from Oliveira and Body et. al., Nuclear Fusion, 2021. https://github.com/SPCData/TCV-X21/blob/main/1.experimental_data/reference_scenario/65402_t1.eqdsk
- `tcv_upper_SN.geqdsk`: TCV upper single null equilibrium g162940.02944_670 from Gabriele Merlo and David Hatch
- `mast.geadsk`: generated with the tutorial from freegs
