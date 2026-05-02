
# Supplementary Figure 2 Notebooks

Notebook files: `fig2a.ipynb`, `fig2bcdefgh.ipynb`

This folder uses multiple notebooks to build Supplementary Figure 2 panels.

## Files in this folder

- `fig2a.ipynb`
- `fig2bcdefgh.ipynb`
- `nucleotide_transformer_scoring.py`
- `nt_scores.csv`
- `nt_llr_generation_wise.pdf`
- `nt_llr_mutator_vs_nonmutator.pdf`
- `nt_llr_vs_evo_score.pdf`
- `regression OLS nucelo.pdf`
- `regression OLS nucelo.png`

## What these notebooks do

- `fig2a.ipynb`
  - Runs nucleotide-transformer masked-marginal LLR scoring for LTEE SNPs.
  - Uses `nucleotide_transformer_scoring.py` as the scoring backend.
  - Writes NT scores table and generates generation-wise trajectory analysis.
  - Includes NT score correlation checks with available external score tables.

- `fig2bcdefgh.ipynb`
  - Loads LTEE gain and allele-count sheets.
  - Builds mono/di/tri-nucleotide flank features.
  - Fits regression workflow to predict `kGain` from flank features.
  - Produces scatter/regression and coefficient heatmap-style panels.

## Data inputs

Expected under `data/`:

- `kgain_all_population_wt.csv`
- `concat_pop_annotation.csv`
- `MetaData_ecoli_final.xlsx` (`Gain score` and `Mastersheet`)
- `m1_annotated_timecourse.xlsx`
- `m2_annotated_timecourse.xlsx`
- `m3_annotated_timecourse.xlsx`
- `m4_annotated_timecourse.xlsx`
- `m5_annotated_timecourse.xlsx`
- `m6_annotated_timecourse.xlsx`
- `p1_annotated_timecourse.xlsx`
- `p2_annotated_timecourse.xlsx`
- `p3_annotated_timecourse.xlsx`
- `p4_annotated_timecourse.xlsx`
- `p5_annotated_timecourse.xlsx`
- `p6_annotated_timecourse.xlsx`
- `GCF_000017985.1_ASM1798v1_genomic.fna`

Optional external table used in `fig2a.ipynb` when present:

- `thermodynamics/sift4g_scores.csv`

## Key parameters

- shared sequence-feature setting:
  - flank k-mer levels: `k = 1, 2, 3`

- `fig2a.ipynb`:
  - output table: `OUT_CSV = 'nt_scores.csv'`
  - model name:
    - `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species`
  - pinned model revision:
    - `06615c1660c892fc199840c18123f8385b3542a8`
  - train/test split in downstream generation-wise block:
    - `random_state = 42`

- `fig2bcdefgh.ipynb`:
  - train/test split:
    - `test_size = 0.2`
    - `random_state = 42`
  - mutator grouping:
    - `mutator_list = ['m1', 'm2', 'm3', 'm4', 'p3', 'p6']`
    - `non_mutator_list = ['p1', 'p2', 'p4', 'p5', 'm5', 'm6']`

## Panel map

- **Supplementary Fig. 2a workflow**
  - Nucleotide-transformer SNP scoring, NT LLR summaries, and generation-wise trend panel (`fig2a.ipynb`).

- **Supplementary Fig. 2b-h workflow**
  - k-mer flank regression and coefficient-heatmap panel set (`fig2bcdefgh.ipynb`).

## Outputs

Saved by notebook cells in current workflow:

- `nt_scores.csv`
- `nt_llr_generation_wise.pdf`
- `regression OLS nucelo.pdf`
- `regression OLS nucelo.png`

Additional generated files currently present in this folder:

- `nt_llr_mutator_vs_nonmutator.pdf`
- `nt_llr_vs_evo_score.pdf`

Note:
- `coefficient.pdf` and `coefficient.png` save lines are present but commented in the current notebook state.

## Dependencies

- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `scipy`
- `statsmodels`
- `sklearn`
- `torch`
- `transformers`
- `kaos`
- local helper:
  - `utility.py`

