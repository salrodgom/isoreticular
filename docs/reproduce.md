# Reproducing the published results

## 1. Zero-pressure energy minimisation (SLC shell model)

```bash
cd optimizations/slc_gulp
# Inputs at structures/initial/G_k_init.cif use Im-3m P1 setting.
# Force-field parameters: catlow.lib (Sanders-Leslie-Catlow).
# The GULP outputs already shipped here used:
#   keywords: opti conp prop
#   minimiser: BFGS, with automatic switch to RFO near saddles.
# Reproduction:
gulp < G_k_input.gin > G_k.gout
```

The expected output should match `bfgs.gout` / `rfo.gout` in each `G_k/` to
within numerical roundoff of GULP 5.1.1.

## 2. LAMMPS NPT setup at T = 298.15 K

```bash
cd opes/G_k/input
# 1. Convert CIF -> LAMMPS data:
analysis/tools/cif2lammps   structures/optimized/slc/G_k_opt_slc_0bar.cif
# 2. Run LAMMPS with PLUMED (4 MPI walkers):
mpirun -np 4  lmp -partition 4x1  -in in.lmp \
       -plumed plumed.dat
```

The OPES_METAD_EXPLORE bias targets the (volume, delta) collective
variables, with a multithermal-multibaric expansion across 1 bar to 18 kbar
and 250 K to 350 K (see `plumed.dat` for the exact ranges). Each walker
writes its own `COLVARS.*` (excluded from the repo) and shares the bias via
`MULTI-STATE`.

## 3. Reweighting at the target thermodynamic state

```bash
cd opes/G_k/reweight_input
# Use the PLUMED driver to reweight the biased samples at p, T of interest:
plumed driver --plumed plumed_REWEIGHT_temperature_pressure.dat \
              --noatoms
```

The output free-energy surfaces are `fes_*.dat` (one per walker). They are
provided in `opes/G_k/fes/`.

## 4. Barrier extraction and figure generation

```bash
cd analysis/fes_extraction
python extract_barriers_isoRHO.py
# Then in figures/:
make
```

`extract_barriers_isoRHO.py` outputs `barriers_extracted.txt`, which
populates Table S in the supplementary information. The `Makefile` in
`figures/` rebuilds the EPS panels from the scripts in `figures/scripts/`
and the data in `figures/data/`.

## 5. Joint fit of p_c(k)

```bash
cd analysis
python joint_pc_fit.py
```

This produces the Cowley-Levanyuk fit p_c(k) = p_0 * exp(-lambda * k)
that appears in the main text.

---

## Software versions (verified)

- GULP 5.1.1 (Gale, Curtin University)
- LAMMPS stable_2Aug2023 with USER-PLUMED
- PLUMED 2.8.0 + OPES module
- Python 3.10, numpy 1.24, scipy 1.10, ase 3.22, pymatgen 2024.6
- gnuplot 5.4.4

If your versions differ, results should be quantitatively identical for
GULP and LAMMPS (deterministic, given the same seed).
