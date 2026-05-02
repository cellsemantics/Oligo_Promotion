
# Figure 8 Notebook

Notebook file: `fig8.ipynb`

This notebook builds Figure 8 mutation-tolerance and `kGain` analyses.

## Files in this folder

- `fig8.ipynb`

## What the notebook does

- Loads per-mutation `kGain` results for multiple genes.
- Merges external tolerance-score annotations.
- Aggregates and pivots data for heatmap-style visualization.
- Creates tolerance-category comparisons of accumulated `kGain`.
- Computes one-sided group tests with BH/FDR-corrected p-values.

## Data inputs

Expected under `data/`:

- `fabZ_complete_result_kgain.csv`
- `lpxC_complete_result_kgain.csv`
- `murA_complete_result_kgain.csv`
- `41467_2023_35940_MOESM8_ESM.csv`

## Key parameters

- tolerance-category cutoff:
  - `ToleranceScore >= 0.8` -> `Tolerant`
  - otherwise -> `Intolerant`
- one-sided test direction: `alternative='greater'`
- BH/FDR correction: `method='fdr_bh'`
- amino-acid plotting order is explicitly defined in the notebook.

## Panel map

- **Heatmap workflow panel(s)**
  - Aggregated `kGain` heatmap by amino-acid mutation dimensions.
  - Tolerance-score heatmap by matching mutation dimensions.

- **Tolerance-category comparison panel(s)**
  - `Tolerant` vs `Intolerant` comparison of accumulated `kGain`, including corrected statistics table.

## Outputs

- In current notebook state, figures are rendered inline.
- No explicit `savefig(...)` output file is defined in this notebook.

## Dependencies

- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `statsmodels`
- `sklearn` (`GaussianMixture` import)
- `pingouin`
- local helper: `utility.py`

## Notes

- The notebook currently references the tolerance file with a leading slash in code in one cell (`/41467_2023_35940_MOESM8_ESM.csv`); ensure your runtime path resolves this correctly.
