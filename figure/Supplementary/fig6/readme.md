
# Supplementary Figure 6 Notebook

Notebook file: `fig6a.ipynb`

This notebook builds Supplementary Figure 6 yeast fitness vs evolved `kGain` analyses.

## Files in this folder

- `fig6a.ipynb`

## What the notebook does

- Loads yeast evolved mutation/gain data.
- Loads population-level fitness metadata.
- Creates allele-switch and generation-expanded mutation tables.
- Computes generation-wise median summaries.
- Fits regression between generation-wise median fitness and median evolved `kGain`.
- Builds scatter-style panel with correlation/regression annotation.

## Data inputs

Expected under `data/`:

- `yeast_evolved_genome_metadata.csv`
- `yeast_metadata.xlsx` (`fitness` sheet)

## Key parameters

- generation columns are converted to numeric generation labels (`70`, `1410`, `2640`, `5150`, `7530`, `10150`) in the notebook workflow.
- regression panel includes:
  - Pearson correlation reporting
  - OLS fit
  - `R^2` reporting (`sklearn.metrics.r2_score`)

## Panel map

- **Supplementary Fig. 6a/6b workflow**
  - generation-wise median fitness vs generation-wise median evolved `kGain` scatter/regression panel (`fig6a.ipynb`).

## Outputs

Saved by notebook cell in current workflow:

- `Yeast median fitness vs median ag scatter plot.pdf`

## Dependencies

- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `scipy`
- `statsmodels`
- `sklearn`
- `PIL`
- `pdf2image`
- local helper:
  - `utility.py`

## Notes

- The notebook uses project-root discovery with `Path.cwd()` and a `utility.py` anchor, then reads files via `DATA_DIR = PROJECT_ROOT / "data"` (relative workflow, not hardcoded absolute paths).
- In this folder, the notebook is named `fig6a.ipynb` (there is no `fig6.ipynb`).
