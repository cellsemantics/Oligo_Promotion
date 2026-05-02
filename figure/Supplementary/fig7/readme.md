
# Supplementary Figure 7 Notebook

Notebook file: `fig7.ipynb`

This notebook builds Supplementary Figure 7 mutation-level `kGain`, tolerance, and RSA-category analyses.

## Files in this folder

- `fig7.ipynb`
- `supp_fig7_lpxc_heatmaps_50x50mm.pdf`
- `supp_fig7_fabz_heatmaps_50x50mm.pdf`
- `supp_fig7_tolerance_boxplots_50x50mm.pdf`
- `supp_fig7_rsa_boxplot_50x50mm.pdf`

## What the notebook does

- Loads per-mutation `kGain` result tables for selected genes.
- Loads mutation tolerance/RSA annotations.
- Merges `kGain` and tolerance metadata by amino-acid mutation key.
- Builds median heatmap summaries for `kGain` and tolerance score.
- Binarizes tolerance categories and compares accumulated `kGain` across groups.
- Compares accumulated `kGain` across RSA categories within genes.

## Data inputs

Expected under `data/`:

- `fabZ_complete_result_kgain.csv`
- `lpxC_complete_result_kgain.csv`
- `murA_complete_result_kgain.csv`
- `41467_2023_35940_MOESM8_ESM.csv`

## Key parameters

- tolerance binarization cutoff:
  - `ToleranceScore >= 1` -> `Tolerant`
  - otherwise -> `Intolerant`
- one-sided group testing direction in boxplot annotation helper:
  - `alternative = 'greater'`
- multiple-testing correction in boxplot annotation helper:
  - `p_adjust_method = 'fdr_bh'`
- amino-acid mutation axes/order are explicitly configured in plotting cells.

## Panel map

- **Supplementary Fig. 7a**
  - lpxC heatmap workflow: median `kGain` and matched tolerance-score heatmaps.

- **Supplementary Fig. 7b**
  - fabZ heatmap workflow: median `kGain` and matched tolerance-score heatmaps.

- **Supplementary Fig. 7c-e**
  - accumulated `kGain` comparisons across tolerance categories.

- **Supplementary Fig. 7f**
  - gene-wise accumulated `kGain` comparison with `RSA_Category` as hue.

## Outputs

Saved by notebook cells in current workflow:

- `supp_fig7_lpxc_heatmaps_50x50mm.pdf`
- `supp_fig7_fabz_heatmaps_50x50mm.pdf`
- `supp_fig7_tolerance_boxplots_50x50mm.pdf`
- `supp_fig7_rsa_boxplot_50x50mm.pdf`

## Dependencies

- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `statsmodels`
- `sklearn` (`GaussianMixture` imported)
- `pingouin`
- local helper:
  - `utility.py`

## Notes


- The notebook uses project-root discovery via a `utility.py` anchor and reads inputs through `DATA_DIR = PROJECT_ROOT / "data"`.
