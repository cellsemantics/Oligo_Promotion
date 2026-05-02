# Figure 2 Notebook (`fig2.ipynb`)

This folder contains one notebook, `fig2.ipynb`, that generates the Figure 2 analyses and panels.

The sections below are aligned to the current notebook logic and outputs.

## Table of contents

- [Overview](#overview)
- [Data inputs](#data-inputs)
- [Key parameters](#key-parameters)
- [Mutator groups](#mutator-groups)
- [Pipeline summary](#pipeline-summary)
- [Figure panels and outputs](#figure-panels-and-outputs)
  - [Fig. 2e](#fig-2e)
  - [Fig. 2f](#fig-2f)
  - [Fig. 2g](#fig-2g)
  - [Fig. 2h](#fig-2h)
  - [Fig. 2i](#fig-2i)
  - [Fig. 2j](#fig-2j)
  - [Fig. 2k-l](#fig-2k-l)
  - [Fig. 2m](#fig-2m)
- [Dependencies](#dependencies)
- [Notes](#notes)

## Overview

`fig2.ipynb`:

- loads LTEE mutation and allele-frequency data,
- computes WT and evolved `kGain`,
- annotates each mutation with mutator type, essentiality, fixation status, and selection class,
- generates figure panels `2e` to `2m`,
- fits the logistic regression used in Fig. 2m.

## Data inputs

All loaded from `data/`:

| File | Usage |
| --- | --- |
| `LTEE_mutational_data.csv` | Base mutation catalog (`mutational_data_all_population`) |
| `MetaData_ecoli_final.xlsx` (`Mastersheet`) | Long-format allele counts by generation |
| `concat_pop_annotation.csv` | Allele-frequency annotations merged into mutation records |
| `gene type.xlsx` | Essential vs non-essential gene labels |
| `GCF_000017985.1_ASM1798v1_genomic.fna` | Reference FASTA for `kGain` computation |

## Key parameters

| Parameter | Value |
| --- | --- |
| `kmer_length` | `10` |
| `target_generation` | `57500` |
| `freq_threshold` (fixation) | `0.95` |
| `min_last_points` (fixation) | `2` |
| `eps` (logit clipping) | `1e-4` |
| Min generations required for regression | `10` |
| Heatmap `bin_size` | `2500` |
| Bootstrap CI seed (effect size) | `42` |
| Bootstrap iterations for effect-size CI | `2000` (`utility.median_based_distance_ci` default) |

## Mutator groups

- `mutator_list = ['m1','m2','m3','m4','p3','p6']`
- `non_mutator_list = ['p1','p2','p4','p5','m5','m6']`

## Pipeline summary

1. **kGain computation**
   - Evolved `kGain` on mutated population FASTA at generation `57,500`.
   - WT `kGain` on reference FASTA.

2. **Annotation merges**
   - Adds mutator/non-mutator labels.
   - Adds essentiality labels.
   - Adds allele-frequency trajectories.

3. **Fixation and selection labels**
   - Fixation rule: `fixed` if the final `min_last_points` values of `allele_freq` are all `>= 0.95`; otherwise `not_fixed`.
   - Fits per-mutation slope with `linregress(logit(allele_freq) ~ generation_number)` if at least `10` distinct generation points exist.
   - Applies BH/FDR to mutation-wise raw p-values.
   - Final class:
     - `positive` if `p_adj < 0.05` and slope `> 0`
     - `negative` if `p_adj < 0.05` and slope `< 0`
     - `neutral` otherwise (including insufficient-data rows)

## Figure panels and outputs

### Fig. 2e

- Displays counts of `classification`.
- Current notebook output:
  - `neutral = 23086`
  - `positive = 9424`
  - `negative = 4330`

### Fig. 2f

- Unique mutation counts by fixation type.
- Output file: `number_of_unique_mutation_count_vs_fixation_type.pdf`

### Fig. 2g

- Unique mutation counts by essentiality type.
- Output file: `number_of_unique_mutation_count_vs_essentiality_type.pdf`

### Fig. 2h

- Unique mutation counts by mutator type.
- Output file: `number_of_unique_mutation_count_vs_mutator_type.pdf`

### Fig. 2i

- Grouped boxplot of `evolved_kGain` by `mutator_type` and `classification`.

Pairwise test setup in the notebook:

- one-sided Mann-Whitney (`utility.man_whiteney(..., alternative="greater")`),
- BH/FDR correction within each mutator family (`mutator` and `non_mutator`, 3 tests each),
- effect size and 95% CI from `utility.median_based_distance_ci`.

Current printed statistics:

- `(('non_mutator','positive'),('non_mutator','neutral'))`: raw `2.972e-17`, BH `8.917e-17`, effect `0.757`, 95% CI `[0.482, 0.952]`
- `(('non_mutator','positive'),('non_mutator','negative'))`: raw `3.028e-02`, BH `4.542e-02`, effect `0.375`, 95% CI `[-0.028, 0.597]`
- `(('non_mutator','neutral'),('non_mutator','negative'))`: raw `1.000e+00`, BH `1.000e+00`, effect `-0.370`, 95% CI `[-0.705, -0.185]`
- `(('mutator','positive'),('mutator','neutral'))`: raw `0.000e+00`, BH `0.000e+00`, effect `0.797`, 95% CI `[0.749, 0.838]`
- `(('mutator','positive'),('mutator','negative'))`: raw `1.630e-134`, BH `2.445e-134`, effect `0.741`, 95% CI `[0.674, 0.834]`
- `(('mutator','neutral'),('mutator','negative'))`: raw `9.947e-01`, BH `9.947e-01`, effect `-0.049`, 95% CI `[-0.099, 0.038]`

### Fig. 2j

- Heatmap of median `evolved_kGain` by population and generation bin.
- Uses:
  - `bin_size = 2500`
  - colormap `BrBG_r`
  - fixed population order:
    - `["m2", "m4", "p3", "m3", "p6", "m1", "m5", "m6", "p2", "p4", "p5", "p1"]`
- Save call is currently commented out.

### Fig. 2k-l

- Stacked bars (`AT->GC` vs `Other`) per population, split by mutator status.
- Save call is currently commented out.

### Fig. 2m

- Logistic regression:
  - target: beneficial (`positive`) vs neutral/negative,
  - features:
    - `evolved_kGain`
    - `is_AT_to_GC` (`1` if REF in `{A,T}` and ALT in `{G,C}`, else `0`)

Current OR output:

- `evolved_kGain: OR = 1.10, 95% CI [1.09, 1.10]`
- `is_AT_to_GC: OR = 1.23, 95% CI [1.17, 1.29]`

## Dependencies

- `numpy`, `pandas`
- `matplotlib`, `seaborn`
- `scipy` (`linregress`)
- `statsmodels` (`Logit`, BH/FDR)
- `statannotations`
- `kaos`
- local helpers: `utility.py`, `r_plot_utils.py`
- optional display helper: `pdf2image`

## Notes

- Project-root discovery is dynamic (walks up to find `utility.py`), so execution does not depend on where Jupyter is launched.
- Fig. 2f, 2g, and 2h use `r_plot_utils.plot_custom_bar_r(...)`, so R support is required for those cells.
