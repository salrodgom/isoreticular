# Methods quick reference

This is a one-page recap. For the full text and equations, see the manuscript.

## Force field

The classical SLC (Sanders-Leslie-Catlow, 1984) core-shell pair potential
for silica:

- O core-shell harmonic spring k_cs.
- Si-O Buckingham + harmonic Si-O-Si angle.
- Charges: Si +4, O core +0.8482, O shell -2.8482 (formal -2 net).

Stored in `optimizations/slc_gulp/catlow.lib`. The same parameters are
used inside the LAMMPS NPT runs via the `hybrid buck/coul/long/cs +
lj/cut/coul/long` setup (see `opes/G_k/input/in.lmp`).

## Collective variables

For each D8R unit (double 8-membered ring) in the supercell:

- `delta_D8R^t` = half of the largest absolute difference between the two
  diagonals across the ring at time t (see Eq. 7 of the manuscript).
- `Delta` = average of `delta_D8R^t` over all D8R units in the supercell
  (CV biased by PLUMED).
- `volume` = NPT cell volume (also biased).

## Enhanced sampling

OPES_METAD_EXPLORE (PLUMED 2.8.x), with:

- ECV_MULTITHERMAL_MULTIBARIC over (energy, volume), spanning
  T = 250-350 K and p = 1-18000 bar.
- ECV_UMBRELLAS_LINE on delta in [0.0, 3.0] Angstrom, sigma = 0.1 A.
- 4 MPI walkers, MULTI-STATE sharing of the bias.
- Pace 500 steps, sigma_e = 1.0 kJ/mol, kBT_e = 2.479 kJ/mol (298.15 K).

The "enhance" parameter and pressure/temperature mesh follow Invernizzi
& Parrinello, JCTC (2020). The same `plumed.dat` template is used for
all five members; only the canonical pressure shifts.

## Reweighting

Multistate Bennett acceptance ratio (MBAR) on (energy, volume, delta) via
`plumed driver`. Output: free-energy surface F(volume, delta) at fixed
T = 298.15 K and at the chosen pressure.

## Barrier extraction

Gaussian-smoothed (sigma = 2 px) interpolation of F(volume, delta) on a
regular grid, then a minimum-energy path traced between the two basins.
Implementation in `analysis/fes_extraction/extract_barriers_isoRHO.py`.

## p_c(k) fit

Joint nonlinear least squares of

  p_c(k) = p_0 * exp(-lambda * k)

across G_1 through G_5, with optional power-law alternatives compared via
information criterion. See `analysis/joint_pc_fit.py`.
