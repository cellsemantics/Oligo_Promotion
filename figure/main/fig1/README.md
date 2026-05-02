# Figure 1 Notebooks

This directory contains the notebooks used to build Figure 1 panels.

## Files in this folder

- `fig1b.ipynb`
- `fig1cde.ipynb`
- `fig1f.ipynb`
- `random_fcgr_latest.pdf` (saved by `fig1b.ipynb`)

## What these notebooks do

- `fig1b.ipynb`
  - tests how FCGR frequency differences vary with k-mer Hamming distance,
  - compares structured k-mer pairs vs random pairs,
  - applies Mann-Whitney testing with BH/FDR correction,
  - produces the Figure 1b distribution plot.

- `fig1cde.ipynb`
  - renders FCGR heatmaps for representative species,
  - generates Figure 1c, 1d, and 1e panels.

- `fig1f.ipynb`
  - builds a shuffled-genome FCGR control (nucleotide counts preserved),
  - generates Figure 1f control heatmap.

## Data inputs

Expected under `data/`:

- `GCF_000017985.1_ASM1798v1_genomic.fna` (used in `fig1b.ipynb`)
- `all speceies FCGR in -log(x) scale truncated to four decimel.txt` (used in `fig1cde.ipynb`)
- `562.5708.fna` (used in `fig1f.ipynb`)

## Key parameters

- `fig1b.ipynb`
  - `kmer_length = 10`
  - `no_of_sample = 100000` (per Hamming distance bin)
  - BH/FDR adjustment via `fdr_bh`
  - random seed set to `42`

- `fig1cde.ipynb`
  - matrix shape is `1024 x 1024` per species
  - colormap typically `RdBu_r`

- `fig1f.ipynb`
  - `kmer_length = 10`
  - shuffled-sequence control with nucleotide counts preserved
  - note: shuffle seed is not fixed in current notebook

## Panel map

- **Fig. 1b**
  - grouped distribution plot of FCGR frequency absolute differences
  - compares Hamming-distance categories and pair type.

- **Fig. 1c-e**
  - species FCGR heatmaps (human, yeast, E. coli).

- **Fig. 1f**
  - shuffled-genome FCGR heatmap control.

## Outputs

- Saved by default in current state:
  - `random_fcgr_latest.pdf` (from `fig1b.ipynb`)
- Other figures are primarily rendered inline.

## Dependencies

- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `statsmodels`
- `kaos`
- local helper: `utility.py`

## Notes

- Project root is discovered dynamically by walking up to `utility.py`.
- `fig1f.ipynb` currently does not lock a shuffle seed, so exact control visualization can vary slightly between runs.
