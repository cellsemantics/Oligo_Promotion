# Figure 3 (`fig3.ipynb`)

This directory contains the notebook used to build Figure 3 analyses and panels.

## Contents

- `fig3.ipynb` — full analysis and plotting workflow
- `fig3_generation_median_with_ma3.pdf` — saved panel for Fig. 3a-c
- `bootstrap_means_essential_vs_non_essential_mutator.pdf` — saved bootstrap panel
- `bootstrap_means_essential_vs_non_essential_non_mutator.pdf` — saved bootstrap panel

## What this notebook does

`fig3.ipynb` performs the following:

1. Loads LTEE mutation and allele-frequency data.
2. Computes WT and evolved `kGain`.
3. Adds annotations (`mutator_type`, `essentiality_status`, `fixation_status`, `classification`).
4. Merges external model scores (`esm_score`, `alt_evo_score`).
5. Computes generation-wise trend statistics with permutation-based Kendall tau (`utility.dba_stat`).
6. Produces panels Fig. 3a through Fig. 3m.

## Data inputs

All expected under `data/`:

| File | Purpose |
| --- | --- |
| `LTEE_mutational_data.csv` | Mutation catalog |
| `MetaData_ecoli_final.xlsx` (`Mastersheet`) | Allele counts by generation |
| `concat_pop_annotation.csv` | Allele-frequency annotations |
| `gene type.xlsx` | Essential vs non-essential labels |
| `GCF_000017985.1_ASM1798v1_genomic.fna` | Reference FASTA for `kGain` |
| `LLR.xlsx` | ESM score table |
| `evo_score.csv` | EVO score table |
| `NIHMS908078-supplement-Supplementary_Table_3.xlsx` | Parallel-gene list (used in Fig. 3l, 3m) |

## Key parameters

| Parameter | Value |
| --- | --- |
| `kmer_length` | `10` |
| `target_generation` | `57500` |
| `freq_threshold` (fixation) | `0.95` |
| `min_last_points` (fixation) | `2` |
| `eps` (logit clipping) | `1e-4` |
| Min generation points for regression | `10` |
| `utility.dba_stat` permutations | `10000` |
| `utility.dba_stat` seed | `42` |
| Bootstrap iterations (Fig. 3f-g) | `10000` |
| Generation-bin step (Fig. 3l, 3m) | `5000` |

Mutator grouping:

- `mutator_list = ['m1', 'm2', 'm3', 'm4', 'p3', 'p6']`
- `non_mutator_list = ['p1', 'p2', 'p4', 'p5', 'm5', 'm6']`

## Panel map

| Panel | Description |
| --- | --- |
| Fig. 3a-c | Generation-wise trajectories of median `evolved_kGain`, `esm_score`, `alt_evo_score`; includes MA(3) overlay and trend testing |
| Fig. 3d | `mutator` boxplot: `evolved_kGain` by `essentiality_status` |
| Fig. 3e | `non_mutator` boxplot: `evolved_kGain` by `essentiality_status` |
| Fig. 3f-g | Bootstrap mean comparison (essential vs non-essential), stratified by mutator type |
| Fig. 3h-i | Two-panel generation-wise median `evolved_kGain` by essentiality (`non_mutator`, `mutator`) |
| Fig. 3j | `mutator` boxplot: `evolved_kGain` by `fixation_status` |
| Fig. 3k | `non_mutator` boxplot: `evolved_kGain` by `fixation_status` |
| Fig. 3l | Parallel-gene subset: generation-binned mean `allele_count` (with SEM), grouped by `mutator_type` |
| Fig. 3m | Parallel-gene subset: generation-binned mean `evolved_kGain` (with SEM), grouped by `mutator_type` |

## Outputs

Saved by default in current notebook state:

- `fig3_generation_median_with_ma3.pdf`
- `bootstrap_means_essential_vs_non_essential_mutator.pdf`
- `bootstrap_means_essential_vs_non_essential_non_mutator.pdf`

Other plots are rendered inline. Some `savefig(...)` lines are present but commented.

## Dependencies

- Core: `numpy`, `pandas`, `matplotlib`, `seaborn`
- Stats: `scipy`, `statsmodels`
- Domain/feature helpers: `kaos`, `fcgr`
- Local helpers: `utility.py`, `r_plot_utils.py`
- Optional imports used in notebook: `pingouin`, `pysam`, `pdf2image`, `PIL`, `sklearn`

## Notes

- The notebook discovers project root dynamically by walking up to `utility.py`.
- Boxplot significance/effect annotations are generated via `utility.return_box_with_p_effect_size(...)`.
