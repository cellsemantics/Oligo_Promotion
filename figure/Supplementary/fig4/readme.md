
# Supplementary Figure 4 Notebook

Notebook file: `fig4.ipynb`

This notebook builds Supplementary Figure 4 COG-function and generation-binned `kGain` analyses.

## Files in this folder

- `fig4.ipynb`
- `median_AG_per_COG_function.pdf`

## What the notebook does

- Loads LTEE mutation and allele-frequency tables.
- Computes evolved and WT `kGain` values.
- Merges annotation tables and mutation-level labels.
- Assigns mutator vs non-mutator classes.
- Maps mutations to operon/COG categories.
- Builds generation-binned COG-function summaries of median `kGain`.

## Data inputs

Expected under `data/`:

- `LTEE_mutational_data.csv`
- `MetaData_ecoli_final.xlsx`
- `concat_pop_annotation.csv`
- `gene type.xlsx`
- `GCF_000017985.1_ASM1798v1_genomic.fna`
- `door_operon_list_60K_paper.txt`

## Key parameters

- `kmer_length = 10`
- `target_generation = 57500`
- fixation settings:
  - `freq_threshold = 0.95`
  - `min_last_points = 2`
- logit clipping:
  - `eps = 1e-4`
- mutator grouping:
  - `mutator_list = ['m1', 'm2', 'm3', 'm4', 'p3', 'p6']`
  - `non_mutator_list = ['p1', 'p2', 'p4', 'p5', 'm5', 'm6']`
- generation binning in COG summary block:
  - `5000`-generation bins

## Panel map

- **Supplementary Fig. 4 workflow**
  - COG-function mapping from operon overlap.
  - generation-binned median `kGain` summaries by mutator class.
  - multi-panel COG-function trend plotting across bins.

## Outputs

Saved by notebook cells in current workflow:

- `median_AG_per_COG_function.pdf`

Note:
- The notebook defines a helper function for independent-scale mutator heatmaps with optional `savepath`, but explicit saved filenames for that helper are not currently set in notebook calls.

## Dependencies

- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `scipy`
- `statsmodels`
- `sklearn`
- `kaos`
- `fcgr`
- `pingouin`
- `pysam`
- `PIL`
- `pdf2image`
- local helpers:
  - `utility.py`
  - `r_plot_utils.py`
