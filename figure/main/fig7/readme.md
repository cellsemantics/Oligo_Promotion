
# Figure 7 Notebooks

Notebook files: `fig7a.ipynb`, `fig7b.ipynb`, `fig7c.ipynb`, `fig7efg.ipynb`

This folder uses multiple notebooks to build Figure 7 panels.

## Files in this folder

- `fig7a.ipynb`
- `fig7b.ipynb`
- `fig7c.ipynb`
- `fig7efg.ipynb`
- `yeast_evolved_kGain_vs_population_with_generation_as_hue_heatmap.pdf`
- `yeast_essential_generation_boxplot_python.pdf`
- `fig7c_kgain_fixed_status_boxplot.pdf`
- `fig7b_stats_summary.csv`
- `fig7c_stats_summary.csv`

## What these notebooks do

- `fig7a.ipynb`
  - Builds generation/population yeast `kGain` heatmap-style panel.
  - Integrates yeast fitness metadata.

- `fig7b.ipynb`
  - Compares `kGain` by generation with essentiality grouping.
  - Exports summary stats table.

- `fig7c.ipynb`
  - Compares `kGain` with fixation-status grouping.
  - Exports summary stats table.

- `fig7efg.ipynb`
  - Performs grouped comparisons on COVID variant metadata:
    - `VOC/VOI` vs `Control`,
    - structural / non-structural / accessory groupings,
    - `kGain` sign categories (`< 0` vs `>= 0`).
  - Exports a grouped stats table.

## Data inputs

Expected under `data/`:

- `yeast_evolved_genome_metadata.csv`
- `yeast_metadata.xlsx` (`fitness` sheet)
- `fixed_mutation.xlsx`
- `COVID_final_mastersheet.xlsx`

## Key parameters

- For `fig7b.ipynb` and `fig7c.ipynb`:
  - one-sided comparison: `alternative='greater'`
  - BH/FDR correction: `p_adjust_method='fdr_bh'`

- For `fig7efg.ipynb`:
  - `Classification` collapsed to:
    - `VOC/VOI`
    - `Control`
  - gene-group order:
    - `Structural`
    - `Non-structural`
    - `Accessory`
  - sign grouping:
    - `< 0`
    - `>= 0`
  - mixed use of corrected and uncorrected reporting depending on the specific test block in the notebook.

## Panel map

- **Fig. 7a**
  - Yeast `kGain` by population/generation heatmap workflow.

- **Fig. 7b**
  - Yeast `kGain` comparison with essentiality/generation grouping.

- **Fig. 7c**
  - Yeast `kGain` comparison with fixation grouping.

- **Fig. 7e / 7f / 7g workflow**
  - COVID variant-group and gene-group comparison panels plus `kGain` sign-category comparison.

## Outputs

Saved by default in current notebook state:

- `yeast_evolved_kGain_vs_population_with_generation_as_hue_heatmap.pdf`
- `yeast_essential_generation_boxplot_python.pdf`
- `fig7c_kgain_fixed_status_boxplot.pdf`
- `fig7b_stats_summary.csv`
- `fig7c_stats_summary.csv`
- `gene_group_stats.csv` (written by `fig7efg.ipynb`)

## Dependencies

- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `scipy`
- `sklearn` (utility metrics import in yeast notebooks)
- `PIL`
- `pdf2image`
- local helper: `utility.py`

## Notes

- `fig7a.ipynb` markdown contains a legacy label (`Fig. 6a`) but the file is part of the Figure 7 workflow.
