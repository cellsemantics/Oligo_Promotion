
# Figure 4 Notebooks

Notebook files: `fig4ab.ipynb`, `fig4cd.ipynb`

These notebooks build Figure 4 analyses and panels.

## Files in this folder

- `fig4ab.ipynb`
- `fig4cd.ipynb`

## What these notebooks do

- `fig4ab.ipynb`
  - Loads LTEE mutation and allele-frequency data.
  - Computes WT and evolved `kGain`.
  - Adds mutation annotations:
    - `mutator_type`
    - `essentiality_status`
    - `fixation_status`
    - `classification`
  - Merges the parallel-gene set.
  - Produces Fig. `4a` and Fig. `4b`.

- `fig4cd.ipynb`
  - Reuses the same mutation/annotation framework.
  - Merges operon and COG-function labels.
  - Produces Fig. `4c` and Fig. `4d` heatmaps.

## Data inputs

Expected under `data/`:

- `LTEE_mutational_data.csv`
- `MetaData_ecoli_final.xlsx` (`Mastersheet`)
- `concat_pop_annotation.csv`
- `gene type.xlsx`
- `GCF_000017985.1_ASM1798v1_genomic.fna`
- `NIHMS908078-supplement-Supplementary_Table_3.xlsx` (parallel-gene list)
- `door_operon_list_60K_paper.txt` (operon/COG mapping input for `fig4cd.ipynb`)

## Key parameters

- `kmer_length = 10`
- `target_generation = 57500`
- `freq_threshold = 0.95` (fixation)
- `min_last_points = 2` (fixation)
- `eps = 1e-4` (logit clipping)
- minimum generation points for regression: `10`
- generation bin step in Fig. `4a-d`: `5000`

Mutator grouping:

- `mutator_list = ['m1', 'm2', 'm3', 'm4', 'p3', 'p6']`
- `non_mutator_list = ['p1', 'p2', 'p4', 'p5', 'm5', 'm6']`

## Panel map

- **Fig. 4a**
  - Generation-binned median `allele_freq` trends.
  - Compares `parallel` vs `non_parallel` gene sets, stratified by mutator class.

- **Fig. 4b**
  - Generation-binned median `evolved_kGain` trends.
  - Compares `parallel` vs `non_parallel` gene sets, stratified by mutator class.

- **Fig. 4c**
  - Mutator-only heatmap of median `evolved_kGain` by `COG_function` and generation bin.

- **Fig. 4d**
  - Non-mutator-only heatmap of median `evolved_kGain` by `COG_function` and generation bin.

## Outputs

- `fig4cd.ipynb` saves:
  - `mutator_AG_heatmap_independent.pdf`
  - `non_mutator_AG_heatmap_independent.pdf`


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
- optional imports in notebooks:
  - `pingouin`
  - `pysam`
  - `pdf2image`
  - `PIL`
  - `sklearn`

## Notes

- Project root is discovered dynamically by walking up to `utility.py`.
- In `fig4ab.ipynb`, per-bin significance markers are based on one-sided Mann-Whitney testing with BH/FDR correction within mutator and non-mutator families.
