
# Supplementary Figure 1 Notebooks

Notebook files: `fig1_1a_1b_1c_1d.ipynb`, `fig1_1e_1f.ipynb`, `fig1h.ipynb`, `fig1gi.ipynb`

This folder uses multiple notebooks to build Supplementary Figure 1 panels.

## Files in this folder

- `fig1_1a_1b_1c_1d.ipynb`
- `fig1_1e_1f.ipynb`
- `fig1h.ipynb`
- `fig1gi.ipynb`
- `population wise LLR.pdf`

## What these notebooks do

- `fig1_1a_1b_1c_1d.ipynb`
  - Loads species FCGR tables.
  - Builds FCGR heatmap panels.
  - Computes species-level correlation matrix and correlation heatmap.

- `fig1_1e_1f.ipynb`
  - Runs k-mer/FCGR profile analysis on reference genomes.
  - Includes E. coli and yeast-focused preparation blocks.

- `fig1h.ipynb`
  - Computes and plots population-wise LLR comparison panel.

- `fig1gi.ipynb`
  - Builds evolved `kGain` workflow and population-level comparisons.
  - Includes mutation annotation merge and statistical comparison reporting.

## Data inputs

Expected under `data/`:

- `all speceies FCGR in -log(x) scale truncated to four decimel.txt`
- `GCF_000005845.2_ASM584v2_genomic.fna`
- `GCF_009858895.2_ASM985889v3_genomic.fna`
- `w303_vlte.fasta`
- `GCF_000017985.1_ASM1798v1_genomic.fna`
- `LTEE_mutational_data.csv`
- `MetaData_ecoli_final.xlsx`
- `concat_pop_annotation.csv`
- `gene type.xlsx`
- `LLR.xlsx`

## Key parameters

- shared k-mer setting in this folder: `kmer_length = 10`
- `fig1gi.ipynb`:
  - `target_generation = 57500`
  - fixation settings:
    - `freq_threshold = 0.95`
    - `min_last_points = 2`
  - logit clipping:
    - `eps = 1e-4`
  - one-sided tests in comparison cells:
    - `alternative='greater'`
  - BH/FDR correction:
    - `p_adjust_method='fdr_bh'`

Mutator grouping used in `fig1gi.ipynb` / `fig1h.ipynb`:

- `mutator_list = ['m1', 'm2', 'm3', 'm4', 'p3', 'p6']`
- `non_mutator_list = ['p1', 'p2', 'p4', 'p5', 'm5', 'm6']`

## Panel map

- **Supplementary Fig. 1a-d**
  - FCGR heatmaps and species-correlation heatmap workflow (`fig1_1a_1b_1c_1d.ipynb`).

- **Supplementary Fig. 1e-f**
  - k-mer/FCGR profile analysis workflow for E. coli and yeast references (`fig1_1e_1f.ipynb`).

- **Supplementary Fig. 1h**
  - population-wise LLR comparison (`fig1h.ipynb`).

- **Supplementary Fig. 1g / 1i**
  - evolved `kGain` population comparison workflow (`fig1gi.ipynb`).

## Outputs

Saved by notebook cells in current workflow:

- `all heatmap trial.png`
- `correlation heatmap.pdf`
- `correlation heatmap.png`
- `kmer.pdf`
- `population wise LLR.pdf`
- `Evolved kGain vs population.pdf`
- `fig1gi_stats_summary.csv`

## Dependencies

- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `scipy`
- `sklearn`
- `Bio` (`SeqIO`)
- `kaos`
- `fcgr`
- `ruptures`
- local helpers:
  - `utility.py`
  - `r_plot_utils.py`

