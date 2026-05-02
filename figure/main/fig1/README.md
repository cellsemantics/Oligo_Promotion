# Figure 1 — Notebooks

This folder contains three notebooks that together produce all panels of Figure 1.
All notebooks discover the project root dynamically (by walking up from the current
working directory until `utility.py` is found), so they run correctly regardless of
where Jupyter is launched from.

---

## Notebooks

### `fig1b.ipynb` → Figure 1b

**Goal:** Show that FCGR frequency differences between k-mer pairs scale with their
Hamming distance — i.e., sequence-similar k-mers occupy similar regions of frequency
space, while random k-mer pairs do not.

**Approach:**
1. Loads the *E. coli* B str. REL606 genome and builds a 10-mer FCGR frequency dictionary.
2. For each Hamming distance 1–10, samples 100,000 random k-mers and generates a
   paired k-mer that differs by exactly that many positions ("original" type), plus
   100,000 independent random k-mer pairs ("random" type).
3. Filters to HD 1–3 (all rows) and a balanced sample from HD 4–9 for the ">3"
   category (HD 10 is excluded from the plot — note: data is generated for HD 10
   but not included in `df_filtered`).
4. Computes the absolute FCGR frequency difference for each pair.
5. Plots a grouped boxplot (log y-scale) of `freq_absolute_difference` by Hamming
   distance category and pair type.
6. Reports category-wise one-sided Mann-Whitney U tests (original < random) with
   Benjamini-Hochberg FDR correction across the four categories.

**Data input:**

| File | Role |
|------|------|
| `data/GCF_000017985.1_ASM1798v1_genomic.fna` | *E. coli* B str. REL606 genome (NCBI NC_012967.1, FASTA) used to build the 10-mer FCGR frequency dictionary |

**Output:** `random_fcgr_latest.pdf` (saved in this folder)

**Key parameters:**

| Parameter | Value | Notes |
|-----------|-------|-------|
| `kmer_length` | 10 | k-mer size for FCGR |
| `no_of_sample` | 100,000 | pairs sampled per Hamming distance |
| `SEED` | 42 | NumPy and Python random seeds |
| `adjust_method` | `fdr_bh` | Benjamini-Hochberg FDR |

**Dependencies:** `numpy`, `pandas`, `matplotlib`, `seaborn`, `statsmodels`, `kaos`, `utility`

---

### `fig1cde.ipynb` → Figure 1c, 1d, 1e

**Goal:** Visualise the FCGR heatmaps for three representative species — Human,
*S. cerevisiae* (yeast), and *E. coli* — to illustrate that each genome has a
distinctive, species-specific FCGR pattern.

**Approach:**
1. Loads pre-computed FCGR values for all species from a single comma-separated file.
2. Reshapes each species' column from a flat vector into a 1024 × 1024 matrix
   (corresponding to a 10-mer FCGR at 2^10 = 1024 resolution).
3. Renders the three heatmaps side-by-side using the `RdBu_r` colormap at 600 dpi.

**Data input:**

| File | Role |
|------|------|
| `data/all speceies FCGR in -log(x) scale truncated to four decimel.txt` | Pre-computed −log(x) FCGR values for all species, stored as comma-separated columns (`human`, `yeast`, `ecoli`, …). Each column has 1,048,576 rows (1024 × 1024). |

**Output:** Displayed inline (save line is commented out).

**Key parameters:**

| Parameter | Value | Notes |
|-----------|-------|-------|
| `selected_cmap` | `RdBu_r` | Colormap for all heatmaps |
| matrix shape | 1024 × 1024 | Matches a 10-mer CGR space |
| `SEED` | 42 | Set for reproducibility (not stochastic in this notebook) |

**Dependencies:** `numpy`, `pandas`, `matplotlib`, `seaborn`, `utility`

---

### `fig1f.ipynb` → Figure 1f

**Goal:** Produce a negative control FCGR heatmap by shuffling the *E. coli* genome
sequence while preserving its exact nucleotide composition. The resulting heatmap
should appear uniform/random, confirming that FCGR structure is not an artefact of
base composition alone.

**Approach:**
1. Reads the *E. coli* genome FASTA.
2. Shuffles the full concatenated sequence with `random_permutation_with_count_preserved`,
   which preserves per-nucleotide counts (A/C/G/T) while randomising their order.
3. Computes the FCGR matrix on the shuffled sequence using `kaos.chaos_frequency_matrix`
   with pseudo-counts.
4. Plots the −log-normalised heatmap at 600 dpi.

**Data input:**

| File | Role |
|------|------|
| `data/562.5708.fna` | *E. coli* CVM N34086PS genome (FASTA, strain 562.5708) used as the source sequence for shuffling |

**Output:** Displayed inline (save line is commented out).

**Key parameters:**

| Parameter | Value | Notes |
|-----------|-------|-------|
| `kmer_length` | 10 | k-mer size for FCGR |
| Nucleotide counts | A: 1,239,016 · C: 1,275,261 · G: 1,268,680 · T: 1,242,049 | Preserved exactly in the shuffled sequence |
| SEED | **not set** | `random.shuffle` is called without a fixed seed — the heatmap will differ slightly between runs |

**Dependencies:** `numpy`, `pandas`, `matplotlib`, `seaborn`, `kaos`

---

## Shared conventions

- **Project root discovery:** All notebooks walk up the directory tree to find the
  folder containing `utility.py` and add it to `sys.path`. No hardcoded paths.
- **SEED = 42:** Used in fig1b and fig1cde for NumPy and Python `random` to ensure
  reproducible sampling. **fig1f does not set a seed** — the nucleotide shuffle
  produces a slightly different matrix on each run (cosmetically equivalent for
  a control panel, but worth noting).
- **`utility.py`:** Provides `custom_figure_axis`, `man_whiteney`, and
  `median_based_distance_ci` used by fig1b; `custom_figure_axis` used by fig1cde.
- **`kaos` package:** Used by fig1b (frequency dictionary) and fig1f (frequency
  matrix). Must be installed in the active environment.
- **PDF font embedding:** fig1b sets `pdf.fonttype = 42` and `ps.fonttype = 42` so
  text is editable in Illustrator/Inkscape.

## Data files at a glance

| File (under `data/`) | Used by | Description |
|----------------------|---------|-------------|
| `GCF_000017985.1_ASM1798v1_genomic.fna` | fig1b | *E. coli* B str. REL606 genome (NC_012967.1) |
| `all speceies FCGR in -log(x) scale truncated to four decimel.txt` | fig1cde | Pre-computed −log FCGR values, all species |
| `562.5708.fna` | fig1f | *E. coli* CVM N34086PS genome for shuffled control |
