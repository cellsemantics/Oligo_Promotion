# Figure 3 Notebook

Notebook file: `fig3.ipynb`

This notebook builds Figure 3 analyses and panels.

## Files in this folder

- `fig3.ipynb`
- `fig3_generation_median_with_ma3.pdf`
- `bootstrap_means_essential_vs_non_essential_mutator.pdf`
- `bootstrap_means_essential_vs_non_essential_non_mutator.pdf`

## What the notebook does

- Loads LTEE mutation and allele-frequency data.
- Computes WT and evolved `kGain`.
- Adds mutation annotations:
  - `mutator_type`
  - `essentiality_status`
  - `fixation_status`
  - `classification`
- Merges external score tables:
  - `esm_score`
  - `alt_evo_score`
- Computes generation-wise trend tests using `utility.dba_stat`.
- Produces Figure 3 panels (`3a` to `3m`).

## Data inputs

Expected under `data/`:

- `LTEE_mutational_data.csv` (mutation catalog)
- `MetaData_ecoli_final.xlsx` (`Mastersheet`, allele counts by generation)
- `concat_pop_annotation.csv` (allele-frequency annotations)
- `gene type.xlsx` (essential vs non-essential labels)
- `GCF_000017985.1_ASM1798v1_genomic.fna` (reference FASTA)
- `LLR.xlsx` (ESM scores)
- `evo_score.csv` (EVO scores)
- `NIHMS908078-supplement-Supplementary_Table_3.xlsx` (parallel-gene list for `3l`, `3m`)

## Key parameters

- `kmer_length = 10`
- `target_generation = 57500`
- `freq_threshold = 0.95` (fixation)
- `min_last_points = 2` (fixation)
- `eps = 1e-4` (logit clipping)
- minimum generation points for regression: `10`
- `utility.dba_stat` permutations: `10000`
- `utility.dba_stat` random seed: `42`
- bootstrap iterations in Fig. `3f-g`: `10000`
- generation bin step in Fig. `3l`, `3m`: `5000`

Mutator grouping:

- `mutator_list = ['m1', 'm2', 'm3', 'm4', 'p3', 'p6']`
- `non_mutator_list = ['p1', 'p2', 'p4', 'p5', 'm5', 'm6']`

## Panel map

- **Fig. 3a-c**
  - Generation-wise trajectories of median `evolved_kGain`, `esm_score`, and `alt_evo_score`.
  - Includes MA(3) overlay and trend testing.

- **Fig. 3d**
  - `mutator` boxplot: `evolved_kGain` by `essentiality_status`.

- **Fig. 3e**
  - `non_mutator` boxplot: `evolved_kGain` by `essentiality_status`.

- **Fig. 3f-g**
  - Bootstrap mean comparison (essential vs non-essential), stratified by mutator type.

- **Fig. 3h-i**
  - Two-panel generation-wise median `evolved_kGain` by essentiality:
    - `non_mutator`
    - `mutator`

- **Fig. 3j**
  - `mutator` boxplot: `evolved_kGain` by `fixation_status`.

- **Fig. 3k**
  - `non_mutator` boxplot: `evolved_kGain` by `fixation_status`.

- **Fig. 3l**
  - Parallel-gene subset:
    - generation-binned mean `allele_count` (with SEM)
    - grouped by `mutator_type`

- **Fig. 3m**
  - Parallel-gene subset:
    - generation-binned mean `evolved_kGain` (with SEM)
    - grouped by `mutator_type`

## Outputs

Saved by default in current notebook state:

- `fig3_generation_median_with_ma3.pdf`
- `bootstrap_means_essential_vs_non_essential_mutator.pdf`
- `bootstrap_means_essential_vs_non_essential_non_mutator.pdf`

Other figures are displayed inline. Some `savefig(...)` lines are currently commented.

## Dependencies

- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `scipy`
- `statsmodels`
- `kaos`
- `fcgr`
- local helpers: `utility.py`, `r_plot_utils.py`
- optional imports in notebook:
  - `pingouin`
  - `pysam`
  - `pdf2image`
  - `PIL`
  - `sklearn`

## Notes

- Project root is discovered dynamically by walking up to `utility.py`.
- Boxplot significance and effect annotation is handled by:
  - `utility.return_box_with_p_effect_size(...)`
