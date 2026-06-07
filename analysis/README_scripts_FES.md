# FES barrier extraction pipeline

End-to-end PLUMED reweighting + grid-based barrier extraction for the SI table
`tab:SI-fes-barriers`. Replaces the visual estimates currently in the
manuscript with grid-derived values at the exact pressures shown in Fig 6.

## Files

| file                              | role                                                          |
|-----------------------------------|---------------------------------------------------------------|
| `template_plumed_REWEIGHT.tpl`    | PLUMED reweighting template (placeholders for walker, T, P)   |
| `run_reweight.sh <G_k>`           | Runs PLUMED driver at each Fig 6 pressure, writes `fes_<P_bar>.dat` |
| `extract_barriers.py`             | Reads all `fes_<P>.dat`, computes ΔG* at each P, outputs LaTeX |
| `barriers_extracted.txt`          | Output snapshot of extract_barriers.py (current state)        |

## Workflow

For each system 1-4, from the project root:

```bash
cd dir_RHO_isoreticular_G1_222_SG_P1/MultiBaric
bash ../../scripts_FES/run_reweight.sh 1
cd -

cd dir_RHO_isoreticular_G2_SG_P1/Multibaric
bash ../../scripts_FES/run_reweight.sh 2
cd -

cd dir_RHO_isoreticular_G3_SG_P1/Multibaric
bash ../../scripts_FES/run_reweight.sh 3
cd -

cd dir_RHO_isoreticular_G4_SG_P1/Multibaric
bash ../../scripts_FES/run_reweight.sh 4
cd -

python3 scripts_FES/extract_barriers.py > scripts_FES/barriers_extracted.txt
```

The extract_barriers.py output ends with a LaTeX block ready to paste into the
SI table; one row per G_k at the pressure that maximises the barrier.

## Pressure ranges per system

|  G_k  | OPES MULTIBARIC range [bar] | Fig 6 pressures [GPa] | Reweighting feasible |
|-------|----------------------------|-----------------------|----------------------|
| G_1   | 1-18000                    | 1.0, 1.2, 1.225, 1.3  | yes (all)            |
| G_2   | 2500-10000                 | 0.5, 0.75, 0.82, 0.9  | yes (all)            |
| G_3   | 1000-8000                  | 0.1, 0.5, 0.7, (1.0)  | yes for 0.1-0.7      |
| G_4   | 500-2000                   | (0.2), 0.5, 0.6, 0.8  | only 0.2; rest out-of-range |
| G_5   | (no Multibaric data)       | 0.3, 0.4, 0.5, 0.7    | needs new MD run     |

Reweighting outside the OPES MULTIBARIC range is mathematically possible but
the statistical uncertainty grows exponentially with the gap; for G_4 figure
pressures above 2000 bar a new OPES Expanded simulation with extended
MULTIBARIC range (e.g. MAX_PRESSURE = 10000 bar) is recommended. G_5
requires generating MDNPT trajectories from scratch.

## Sanity checks before running

The `run_reweight.sh` script ships with hard-coded reference energy and
volume per G_k, estimated from the first 1000 lines of `ALLCOLVARS.0`:

| G_k | E_ref [kJ/mol] | V_ref [nm^3] | P_sim [bar] |
|-----|---------------:|-------------:|------------:|
| G_1 | -4754832       | 21.8         | 10000       |
| G_2 | -2971300       | 14.3         | 7500        |
| G_3 | -8319800       | 39.4         | 3500        |
| G_4 | -17830235      | 86.0         | 1000        |

These shifts only affect numerical precision in REWEIGHT_TEMP_PRESS; the
reweighted distribution is invariant to the choice. Verify with
`head -3 ALLCOLVARS.0` that the values are reasonable for your run.

## Output naming

`run_reweight.sh` writes `fes_<P_bar>.dat` (e.g. `fes_12250.dat` for
1.225 GPa). `extract_barriers.py` parses the integer suffix as P in bar
(divides by 10000 to get GPa). Legacy files named `fes_11400.dat`,
`fes_11500.dat`, `fes_11600.dat` are still parsed by the regex but will be
tagged as 1.14, 1.15, 1.16 GPa regardless of which system they belong to;
delete or rename them before running `extract_barriers.py` for the final
table:

```bash
for d in dir_RHO_isoreticular_G{2,3,4}_SG_P1/Multibaric ; do
  rm -f $d/fes_11400.dat $d/fes_11500.dat $d/fes_11600.dat
done
```
