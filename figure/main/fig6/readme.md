
# Figure 6 Notebook

Notebook file: `fig6bc.ipynb`

This notebook builds Figure 6 analyses and panels (6b/6c workflow).

## Files in this folder

- `fig6bc.ipynb`
- `kGain vs population with synonymus mutation type.pdf`
- `kGain vs population AF wise lab data.pdf`
- `kGain vs population AF wise lab data unique.pdf`

## What the notebook does

- Loads annotated VCF-derived mutation records and reference FASTA.
- Computes evolved and WT `kGain` values.
- Builds population labels (`D`, `R1`, `R3`) and mutation metadata.
- Creates category comparisons for:
  - synonymous vs non-synonymous style grouping (as used in the notebook),
  - `AF = 1` vs `AF < 1`.
- Produces publication panels and per-population significance labels.

## Data inputs

Expected under `data/`:

- `GCF_000005845.2_ASM584v2_genomic.fna`
- `vcf/` directory with annotated VCF files used in the notebook

## Key parameters

- `kmer_length = 10`
- read-depth filter in VCF preprocessing: `DP >= 10`
- one-sided testing in panel stats: `alternative='greater'`
- BH/FDR correction: `method='fdr_bh'`
- bootstrap CI settings in utility cell:
  - `n_boot = 2000`
  - `random_seed = 42`
- AF grouping rule:
  - `AF < 1`
  - `AF = 1`
  - NaN AF values are kept as NaN (not forced into `AF = 1`)

## Panel map

- **Fig. 6b**
  - Box/violin-style comparison of `kGain` across populations (`D`, `R1`, `R3`) with mutation-type grouping.

- **Fig. 6c**
  - Box/violin-style comparison of `kGain` for `AF = 1` vs `AF < 1` across populations.

- **Unique-mutation variant panel**
  - Same AF comparison logic on the unique-mutation subset used in the notebook.

## Outputs

Saved by default in current notebook state:

- `kGain vs population with synonymus mutation type.pdf`
- `kGain vs population AF wise lab data.pdf`
- `kGain vs population AF wise lab data unique.pdf`

## Dependencies

- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `scipy`
- `statsmodels`
- `Bio` (`SeqIO`)
- `tqdm`
- `fcgr`
- `kaos`
- local helper: `utility.py`

## Notes

- Population labels are derived from VCF source names and then collapsed to `D`, `R1`, and `R3`.
