
# Supplementary Figure 3 Notebooks

Notebook files: `fig3abc.ipynb`, `fig3de.ipynb`, `fig3fghijk.ipynb`

This folder uses multiple notebooks to build Supplementary Figure 3 panels.

## Files in this folder

- `fig3abc.ipynb`
- `fig3de.ipynb`
- `fig3fghijk.ipynb`
- `Mutator type wise LLR with gene essentiality.pdf`
- `bootstrap_means_fixed_vs_non_fixed_mutator.pdf`
- `bootstrap_means_fixed_vs_non_fixed_non_mutator.pdf`
- `kGain vs coding type for mutator.pdf`
- `kGain vs coding type for non mutator.pdf`
- `Ecoli kGain vs coding type with gene essentiality as hue for mutator.pdf`
- `Evolved median kGain vs generation with mutator type and essentiality as hue.pdf`

## What these notebooks do

- `fig3abc.ipynb`
  - Loads LTEE gain, allele-count, and LLR score tables.
  - Merges mutator/essentiality annotations.
  - Builds LLR comparison and generation-wise LLR trend panels.

- `fig3de.ipynb`
  - Computes evolved and WT `kGain` values.
  - Merges mutation and annotation tables.
  - Assigns fixation status per mutation.
  - Runs bootstrap-style fixed vs non-fixed comparisons by mutator class.

- `fig3fghijk.ipynb`
  - Continues evolved/WT `kGain` analysis with coding-type and essentiality groupings.
  - Builds generation-wise median `kGain` trends by mutator type and gene class.

## Data inputs

Expected under `data/`:

- `LTEE_mutational_data.csv`
- `MetaData_ecoli_final.xlsx`
- `concat_pop_annotation.csv`
- `gene type.xlsx`
- `LLR.xlsx`
- `GCF_000017985.1_ASM1798v1_genomic.fna`

## Key parameters

- shared `kGain` setting:
  - `kmer_length = 10`

- evolved-sequence generation:
  - `target_generation = 57500`

- fixation settings (`fig3de.ipynb`):
  - `freq_threshold = 0.95`
  - `min_last_points = 2`
  - `eps = 1e-4` (logit clipping)

- mutator grouping:
  - `mutator_list = ['m1', 'm2', 'm3', 'm4', 'p3', 'p6']`
  - `non_mutator_list = ['p1', 'p2', 'p4', 'p5', 'm5', 'm6']`

- statistical reporting includes BH/FDR-corrected p-values in relevant comparison cells.

## Panel map

- **Supplementary Fig. 3a-c workflow**
  - LLR vs mutator/essentiality comparisons and generation-wise median LLR trends (`fig3abc.ipynb`).

- **Supplementary Fig. 3d-e workflow**
  - fixed vs non-fixed `kGain` bootstrap comparison panels by mutator class (`fig3de.ipynb`).

- **Supplementary Fig. 3f-k workflow**
  - coding-type and essentiality-stratified `kGain` comparisons plus generation-wise trend panels (`fig3fghijk.ipynb`).

## Outputs

Saved by notebook cells in current workflow include:

- `Mutator type wise LLR with gene essentiality.pdf`
- `bootstrap_means_fixed_vs_non_fixed_mutator.pdf`
- `bootstrap_means_fixed_vs_non_fixed_non_mutator.pdf`
- `kGain vs coding type for mutator.pdf`
- `Evolved median kGain vs generation with mutator type and essentiality as hue.pdf`


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

