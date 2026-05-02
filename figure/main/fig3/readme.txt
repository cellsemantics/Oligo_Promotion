
# Figure 3 Notebook (`fig3.ipynb`)

This folder contains one notebook, `fig3.ipynb`, that builds the Figure 3 analyses and panels.

The documentation below is aligned to the current notebook code and focuses on methods, panel construction, and files.

## Table of contents

- [Overview](#overview)
- [Data inputs](#data-inputs)
- [Key parameters](#key-parameters)
- [Mutator groups](#mutator-groups)
- [Pipeline summary](#pipeline-summary)
- [Figure panels and outputs](#figure-panels-and-outputs)
  - [Fig. 3a-c](#fig-3a-c)
  - [Fig. 3d](#fig-3d)
  - [Fig. 3e](#fig-3e)
  - [Fig. 3f-g](#fig-3f-g)
  - [Fig. 3h-i](#fig-3h-i)
  - [Fig. 3j](#fig-3j)
  - [Fig. 3k](#fig-3k)
  - [Fig. 3l](#fig-3l)
  - [Fig. 3m](#fig-3m)
- [Saved files](#saved-files)
- [Dependencies](#dependencies)
- [Notes](#notes)

## Overview

`fig3.ipynb`:

- loads LTEE mutation and allele-frequency data,
- computes evolved and WT `kGain` (same framework as Fig2),
- derives mutation-level labels (`mutator_type`, `essentiality_status`, `fixation_status`, `classification`),
- merges ESM and EVO scores,
- computes trend statistics using permutation-based Kendall tau (`utility.dba_stat`),
- generates panels `3a` through `3m`.

## Data inputs

Loaded from `data/`:

| File | Usage |
| --- | --- |
| `LTEE_mutational_data.csv` | Base mutation catalog |
| `MetaData_ecoli_final.xlsx` (`Mastersheet`) | Allele counts by generation |
| `concat_pop_annotation.csv` | Allele-frequency annotations |
| `gene type.xlsx` | Essential vs non-essential labels |
| `GCF_000017985.1_ASM1798v1_genomic.fna` | Reference FASTA for `kGain` |
| `LLR.xlsx` | ESM score table (`esm_score`) |
| `evo_score.csv` | EVO score table (`alt_evo_score`) |
| `NIHMS908078-supplement-Supplementary_Table_3.xlsx` | Parallel-gene list used in Fig. 3l and 3m |

## Key parameters

| Parameter | Value |
| --- | --- |
| `kmer_length` | `10` |
| `target_generation` | `57500` |
| `freq_threshold` (fixation) | `0.95` |
| `min_last_points` (fixation) | `2` |
| `eps` (logit clipping) | `1e-4` |
| Min generations for regression | `10` |
| `utility.dba_stat` permutations | `10000` |
| `utility.dba_stat` seed | `42` |
| Bootstrap iterations (Fig. 3f-g) | `10000` |
| Bootstrap sample fraction (Fig. 3f-g) | `0.9 * min(group sizes)` |
| Fig. 3l / 3m generation bin step | `5000` |

## Mutator groups

- `mutator_list = ['m1','m2','m3','m4','p3','p6']`
- `non_mutator_list = ['p1','p2','p4','p5','m5','m6']`

## Pipeline summary

1. **kGain tables**
   - computes evolved `kGain` on mutated population FASTA at generation 57,500,
   - computes WT `kGain` on reference FASTA.

2. **Annotation and labels**
   - adds mutator class and essentiality class,
   - merges allele-frequency trajectories,
   - derives fixation labels from last two allele-frequency points (`>= 0.95`),
   - performs mutation-wise regression:
     - `logit(allele_freq) ~ generation_number`,
     - requires at least 10 generation points,
     - BH/FDR adjustment on raw p-values,
     - classifies each mutation as `positive`, `negative`, or `neutral`.

3. **External score merges**
   - merges ESM (`LLR.xlsx`) and EVO (`evo_score.csv`) onto mutation records.

4. **Trend stats**
   - runs `utility.dba_stat` for generation-wise trends of:
     - `evolved_kGain`,
     - `esm_score`,
     - `alt_evo_score`,
   - applies BH/FDR across those three tests.

## Figure panels and outputs

### Fig. 3a-c

Three generation-wise trajectories:

- median `evolved_kGain`,
- median `esm_score`,
- median `alt_evo_score`.

The publication-style version overlays:

- raw median with 95% CI (Seaborn bootstrap, `n_boot=2000`),
- MA(3) smoothed median.

The notebook also computes permutation-based Kendall trend statistics with BH/FDR correction for the three trajectories.

### Fig. 3d

Boxplot in `mutator` populations:

- x: `essentiality_status` (`essential` vs `non-essential`)
- y: `evolved_kGain`

Statistical annotation (p-value, effect size, CI) is added directly from utility helpers.

### Fig. 3e

Boxplot in `non_mutator` populations:

- x: `essentiality_status`
- y: `evolved_kGain`

Statistical annotation (p-value, effect size, CI) is added directly from utility helpers.

### Fig. 3f-g

Bootstrap comparison of mean `evolved_kGain` (essential vs non-essential), stratified by mutator type.

Each group is repeatedly resampled (without replacement) and summarized as a bootstrap distribution with 95% confidence intervals.

### Fig. 3h-i

Two-panel line plots of generation-wise median `evolved_kGain` by essentiality status:

- panel 1: `non_mutator`
- panel 2: `mutator`

### Fig. 3j

Boxplot in `mutator` populations:

- x: `fixation_status` (`fixed` vs `not_fixed`)
- y: `evolved_kGain`

Statistical annotation (p-value, effect size, CI) is added directly from utility helpers.

### Fig. 3k

Boxplot in `non_mutator` populations:

- x: `fixation_status`
- y: `evolved_kGain`

Statistical annotation (p-value, effect size, CI) is added directly from utility helpers.

### Fig. 3l

Using only genes in `NIHMS908078-supplement-Supplementary_Table_3.xlsx` (labeled `parallel`):

- generation-binned (`5000`) trend of mean `allele_count`,
- grouped by `mutator_type`,
- with SEM error bars.

### Fig. 3m

Same parallel-gene subset and binning as Fig. 3l, plotting:

- mean `evolved_kGain` by generation bin and `mutator_type`,
- with SEM error bars.

## Saved files

The notebook currently saves:

- `fig3_generation_median_with_ma3.pdf` (Fig. 3a-c publication-style panel),
- `bootstrap_means_essential_vs_non_essential_mutator.pdf` (Fig. 3f-g mutator),
- `bootstrap_means_essential_vs_non_essential_non_mutator.pdf` (Fig. 3f-g non-mutator).

Other figure cells are currently rendered inline; several `savefig(...)` lines are commented out.

## Dependencies

- `numpy`, `pandas`
- `matplotlib`, `seaborn`
- `scipy` (`linregress`, `mannwhitneyu`, others imported)
- `statsmodels` (`multipletests`)
- `sklearn` (imported; confusion matrix utilities)
- `kaos`, `fcgr`
- local helpers: `utility.py`, `r_plot_utils.py`
- optional: `pingouin`, `pysam`, `pdf2image`, `PIL`

## Notes

- Project root is discovered dynamically by walking up to `utility.py`; paths are not hardcoded to notebook launch directory.
- Significance/effect-size annotation in boxplot panels uses `utility.return_box_with_p_effect_size(...)` (Mann-Whitney + median-based effect size + bootstrap CI).
