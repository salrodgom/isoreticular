# Companion Zenodo deposit

The bulky raw data accompanying this repository (LAMMPS NPT trajectories,
full PLUMED COLVARS streams, DeltaFs.data history of the OPES expansion,
GULP phonon dispersions) are not on GitHub. They are archived on Zenodo with
a citable DOI.

> **DOI:** [10.5281/zenodo.20583814](https://doi.org/10.5281/zenodo.20583814)
> URL: https://doi.org/10.5281/zenodo.20583814

## What lives on Zenodo, not here

| Bundle | Approx. size | Per-G_k | Description |
|---|---|---|---|
| MDNPT_xtc.tar | 5-15 GB | each | LAMMPS-PLUMED NPT trajectories (4 walkers per G_k) |
| COLVARS_full.tar | 20-40 MB | each | Full per-step PLUMED COLVARS.0-3 streams |
| DeltaFs_history.tar | 50-200 MB | each | OPES expansion biases history |
| GULP_phonons.tar | 50-100 MB | each | restart_noGamma outputs with full dispersion |

## How to link them together

1. Tag this repository at `v1.0-arxiv`.
2. Upload the bundles to Zenodo and create a release with the same tag,
   so GitHub auto-deposits the source snapshot and Zenodo mints the DOI.
3. Update the DOI placeholder above and in `CITATION.cff`.
4. Add the DOI badge to the top of `README.md`.

## Checksums

Each bundle ships with an accompanying `.sha256` file. To verify after
download:

```bash
sha256sum -c MDNPT_xtc.tar.sha256
```
