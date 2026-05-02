# Figure 2 Notebook

Notebook file: `fig2.ipynb`

This notebook builds Figure 2 analyses and panels.

## Files in this folder

- `fig2.ipynb`
- `number_of_unique_mutation_count_vs_fixation_type.pdf`
- `number_of_unique_mutation_count_vs_essentiality_type.pdf`
- `number_of_unique_mutation_count_vs_mutator_type.pdf`

## What the notebook does

- Loads LTEE mutation and allele-frequency data.
- Computes WT and evolved `kGain`.
- Adds mutation annotations:
  - `mutator_type`
  - `essentiality_status`
  - `fixation_status`
  - `classification`
- Runs mutation-level regression on allele-frequency trajectories.
- Applies BH/FDR correction for selection classification.
- Builds Figure 2 panels (`2e` to `2m`).
- Fits logistic regression for beneficial vs non-beneficial labeling.

## Data inputs

Expected under `data/`:

- `LTEE_mutational_data.csv`
- `MetaData_ecoli_final.xlsx` (`Mastersheet`)
- `concat_pop_annotation.csv`
- `gene type.xlsx`
- `GCF_000017985.1_ASM1798v1_genomic.fna`

## Key parameters

- `kmer_length = 10`
- `target_generation = 57500`
- `freq_threshold = 0.95` (fixation)
- `min_last_points = 2` (fixation)
- `eps = 1e-4` (logit clipping)
- minimum generation points for regression: `10`
- heatmap generation bin size: `2500`
- effect-size CI seed: `42`
- effect-size CI bootstrap iterations: `2000` (default utility setting)

Mutator grouping:

- `mutator_list = ['m1', 'm2', 'm3', 'm4', 'p3', 'p6']`
- `non_mutator_list = ['p1', 'p2', 'p4', 'p5', 'm5', 'm6']`

## Panel map

- **Fig. 2e**
  - mutation `classification` count summary.

- **Fig. 2f**
  - unique mutation counts by `fixation_status`.

- **Fig. 2g**
  - unique mutation counts by `essentiality_status`.

- **Fig. 2h**
  - unique mutation counts by `mutator_type`.

- **Fig. 2i**
  - grouped `evolved_kGain` boxplot by `mutator_type` and `classification`.
  - uses one-sided Mann-Whitney tests with BH/FDR within mutator family.

- **Fig. 2j**
  - heatmap of median `evolved_kGain` by population and generation bin.

- **Fig. 2k-l**
  - stacked mutation-class bars (`AT->GC` vs `Other`) by population and mutator status.

- **Fig. 2m**
  - logistic regression with predictors:
    - `evolved_kGain`
    - `is_AT_to_GC`

## Outputs

Saved by default in current notebook state:

- `number_of_unique_mutation_count_vs_fixation_type.pdf`
- `number_of_unique_mutation_count_vs_essentiality_type.pdf`
- `number_of_unique_mutation_count_vs_mutator_type.pdf`

Other figure cells are rendered inline; several `savefig(...)` lines are commented.

## Dependencies

- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `scipy`
- `statsmodels`
- `statannotations`
- `kaos`
- local helpers: `utility.py`, `r_plot_utils.py`
- optional imports in notebook:
  - `pdf2image`
  - `pingouin`
  - `PIL`
  - `sklearn`

## Notes

- Project root is discovered dynamically by walking up to `utility.py`.
- Fig. 2f/2g/2h rely on `r_plot_utils.plot_custom_bar_r(...)` and require R support.
