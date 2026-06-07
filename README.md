# Isoreticular RHO zeolite family — reproducibility repository

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20583813.svg)](https://doi.org/10.5281/zenodo.20583813)
[![License: MIT](https://img.shields.io/badge/License%20(code)-MIT-blue.svg)](LICENSE-CODE)
[![License: CC BY 4.0](https://img.shields.io/badge/License%20(data)-CC%20BY%204.0-lightgrey.svg)](LICENSE-DATA)

Companion data and code for the manuscript:

> **A topology-tuned pressure valve across the isoreticular RHO zeolite family**
> S. R. G. Balestra, A. Rivas-Blanco, S. Hamad, A. R. Ruíz-Salvador (2026)

This repository archives the inputs, force fields, optimisation logs, PLUMED
logs, free-energy surfaces and figure-generation scripts used in the paper.
Large raw trajectories and full COLVAR streams are deposited separately in
Zenodo (see `docs/zenodo_link.md`).

**Citable DOI of this snapshot:** [10.5281/zenodo.20583813](https://doi.org/10.5281/zenodo.20583813)

---

## Layout

```
isoreticular/
├── structures/
│   ├── initial/                    Input CIF files (Im-3m P1, one per G_k)
│   └── optimized/
│       ├── slc/                    SLC core-shell, p = 0 bar
│       └── dft_r2scan_rvv10/G2_scan/   DFT r2SCAN+rVV10 pressure scan (G2)
├── optimizations/
│   └── slc_gulp/
│       ├── G1/.../G5/              GULP outputs at p=0 (bfgs + RFO when used)
│       ├── symmetric_centric/      Symmetry-preserved Im-3m references
│       └── catlow.lib              Sanders-Leslie-Catlow shell parameters
├── opes/
│   ├── G1/, G2/, G3/, G4/, G5/     OPES_METAD_EXPLORE runs (per G_k)
│   │   ├── input/   plumed.dat, in.lmp
│   │   ├── logs/    plumed.log.0-3 (4 MPI walkers)
│   │   ├── fes/     fes_*.dat reweighted free-energy surfaces
│   │   └── README.md
│   └── G1/reweight_input/          Reweighting PLUMED scripts
├── analysis/
│   ├── fes_extraction/             Barrier extraction from FES
│   ├── reweighting/                Multistate reweighting helpers
│   ├── structural/                 8MR/4MR/6MR distortion analysis
│   ├── tools/                      cif2lammps and trajectory converters
│   ├── joint_pc_fit.py             Fit of p_c(k) decay
│   └── Phase_diagram.py
├── figures/
│   ├── scripts/                    gnuplot .gp + Python prep_*.py
│   ├── data/                       Tabulated .dat used by the scripts
│   ├── Makefile
│   └── README_manuscript_figures.md
└── docs/
    ├── reproduce.md
    ├── methods.md
    └── zenodo_link.md
```

---

## How to navigate

Each subtree is self-contained:

- **structures/** — pure CIF files for the five members of the isoreticular RHO
  family (G_1 through G_5). The `initial/` snapshots are the unit cells used as
  seed for all simulations; `optimized/slc/G_k_opt_slc_0bar.cif` is the final
  zero-pressure equilibrium under the shell-model SLC potential.

- **optimizations/** — full GULP outputs (`*.gout`) for the zero-pressure
  energy minimisation. Two variants are kept whenever they were run:
  the BFGS default (`bfgs.gout`) and the rational-function (`rfo.gout`)
  switch invoked near saddle points. The `symmetric_centric/` subfolder
  contains the parent Im-3m reference optimisations for G_1, G_6, G_7.

- **opes/** — one folder per member with the PLUMED input that drives
  OPES_METAD_EXPLORE on the (volume, delta) collective variables, the LAMMPS
  NPT input file, the 4-walker PLUMED logs and the reweighted free-energy
  surfaces. Run parameters are summarised in each `README.md`. Large COLVARS
  and trajectories are NOT in this repository; see `docs/zenodo_link.md`.

- **analysis/** — Python helpers used in the paper: barrier extraction with
  Gaussian smoothing, multistate reweighting templates, 8/4/6-MR distortion
  histograms, and the joint Cowley-Levanyuk fit of p_c(k).

- **figures/** — every gnuplot script and every `.dat` file needed to redraw
  the main-text and supplementary figures of the paper. The included Makefile
  rebuilds the EPS panels.

---

## Reproducing the published results

See `docs/reproduce.md`. In short:

1. Optimise each `structures/initial/G_k_init.cif` with GULP and `catlow.lib`
   (or skip and use the provided `structures/optimized/slc/G_k_opt_slc_0bar.cif`).
2. Convert to LAMMPS data with `analysis/tools/cif2lammps.f90`.
3. Launch the LAMMPS+PLUMED OPES run from `opes/G_k/input/`.
4. Reweight at T = 298.15 K with the PLUMED driver scripts in
   `analysis/reweighting/` and `opes/G_k/reweight_input/`.
5. Extract barriers and produce figures with the scripts in `analysis/` and
   `figures/`.

---

## Software versions

| Tool | Version used | Role |
|---|---|---|
| GULP | 5.1.1 | Symmetry-preserved and broken-symmetry optimisations |
| LAMMPS | 23 Jun 2022 | NPT MD with shell model |
| PLUMED | 2.8.0 with OPES module | Enhanced sampling and reweighting |
| Python | 3.10+ | Analysis (numpy, scipy, matplotlib, ase, pymatgen) |
| gnuplot | 5.4+ | Figure rendering |

---

## Licence

Code (`*.py`, `*.gp`, `*.sh`, `*.f90`, Makefile): MIT, see `LICENSE-CODE`.
Data (CIF files, logs, FES, .dat): CC-BY-4.0, see `LICENSE-DATA`.
The original `LICENSE` (Unlicense) is preserved for historical reference but
the dual-licensed terms above govern subsequent versions.

---

## Citation

If you reuse any of this material, please cite the paper and this repository.
See `CITATION.cff` for machine-readable metadata.
