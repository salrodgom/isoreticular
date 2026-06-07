# G5 OPES run summary

## Files

| Path | Description |
|---|---|
| input/plumed.dat   | PLUMED collective-variable and OPES_METAD_EXPLORE input |
| input/in.lmp       | LAMMPS NPT input (core-shell SLC potential, T=298.15 K) |
| logs/plumed.log.*  | PLUMED log for each MPI walker (4 walkers) |
| fes/fes_*.dat      | Reweighted free-energy surface F(volume, delta) per walker |

## Parameters

- Force field: Sanders-Leslie-Catlow shell-model for SiO2 (`catlow.lib`)
- Bias: OPES_METAD_EXPLORE (multithermal-multibaric expansion + umbrellas on delta)
- Replicas: 4 MPI walkers, MULTI-STATE shared
- Temperature window: 250-350 K (target 298.15 K)
- Pressure window: 1 bar to 18 kbar
- Canonical replica used here: p = 3000 bar (~0.30 GPa)
- Time step: 0.25 fs, t_run >= 5 ns per walker

## Notes

- The 30 MB `all_Colvar.data` and the per-walker `COLVARS` (~5.4 MB each) and the
  raw LAMMPS trajectories (`MDNPT.*.xtc`, ~hundreds of MB) are NOT in this repository.
- The `DeltaFs.data` (12 MB) used to drive OPES expansion is also not included.
- Those large files will be archived in the companion Zenodo deposit (see ../docs/zenodo_link.md).
- Reweighting templates used to convert biased MD into FES are in ../analysis/.
