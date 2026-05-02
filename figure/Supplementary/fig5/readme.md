
# Supplementary Figure 5 Notebook

Notebook file: `fig5bc.ipynb`

This notebook builds Supplementary Figure 5b-c lab-evolution mutation summaries with WT/evolved `kGain` preprocessing.

## Files in this folder

- `fig5bc.ipynb`
- `kGain vs population AF wise lab data.pdf`
- `Lab experiment fix vs non fix.pdf`

## What the notebook does

- Loads annotated VCF records from the lab-evolution input directory.
- Applies mutation-level quality/depth and SNV-only filters.
- Extracts reference/alternate flanks around each mutation.
- Maps mutation `SOURCE` labels to population (`D`, `R1`, `R3`) and generation indices.
- Computes mutation recurrence counts and assigns fixation-style labels (`fixed` vs `not_fixed`).
- Extracts allele frequency (`AF`) from INFO and assigns AF categories.
- Computes WT and evolved (generation 5) mutation-level `kGain` summaries.
- Builds stacked bar summaries used for Supplementary Fig. 5b and 5c.

## Data inputs

Expected under `data/`:

- `GCF_000005845.2_ASM584v2_genomic.fna`
- `vcf/*_annotated.vcf` (glob pattern)

## Key parameters

- `kmer_length = 10`
- VCF filtering:
  - `QUAL >= 20`
  - `DP >= 10` (parsed from INFO)
  - single-nucleotide substitutions only (`len(REF) == 1` and `len(ALT) == 1`)
- recurrence-based mutation label:
  - `fixed` if mutation count `>= 4`
  - otherwise `not_fixed`
- AF category:
  - `AF < 1`
  - `AF = 1`
- evolved kGain reference point:
  - mutated-population FASTA at generation `5`

Stat-annotation helper block defined in notebook:

- one-sided Welch t-test (`alternative='greater'`)
- optional BH/FDR correction (`p_adjust_method='fdr_bh'`)
- median/MAD effect size with bootstrap CI (`n_boot=2000`, `ci=95`, `random_state=42`)

## Panel map

- **Supplementary Fig. 5b**
  - stacked population-wise variant counts by `AF_category`.

- **Supplementary Fig. 5c**
  - stacked population-wise variant counts by `mutation_type` (`fixed` vs `not_fixed`).

## Outputs

Saved by notebook cells in current workflow:

- `Lab experiment fix vs non fix.pdf`

Additional file currently present in this folder:

- `kGain vs population AF wise lab data.pdf`

## Dependencies

- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `scipy`
- `statsmodels`
- `Bio` (`SeqIO`)
- `tqdm`
- `glob`
- `kaos`
- `fcgr`
- local helper:
  - `utility.py`

## Notes

- The notebook discovers project root dynamically by walking up to `utility.py` and resolves data with `DATA_DIR = PROJECT_ROOT / "data"`.
- Inputs are read from relative paths (`data/` and `data/vcf/`), not hardcoded absolute paths.
