
# Figure 5 Notebooks

Notebook files: `transformer_clean_final_updated.ipynb`, `fig5a.ipynb`, `fig5bcd.ipynb`, `fig5ef.ipynb`

This folder trains an attention-based regression model that predicts kGain from the one-hot encoded mutation context window, and produces all Figure 5 panels.

## Files in this folder

- `transformer_clean_final_updated.ipynb` — model training
- `fig5a.ipynb` — encoding schematic
- `fig5bcd.ipynb` — training loss curve panels
- `fig5ef.ipynb` — model evaluation panels
- `attn_regressor_focal.weights.h5` — saved model weights (output of training)
- `history.pkl` — epoch-wise training history dict (output of training)
- `epoch_wise_loss.csv` — tabular training/validation loss per epoch
- `custom_error_kde.pdf` — prediction-error KDE plot

## What these notebooks do

### `transformer_clean_final_updated.ipynb`
- Loads `kgain_all_population_wt.csv` (one row per LTEE mutation with flanking sequence and wild-type kGain).
- Encodes each mutation as a 19 × 4 one-hot difference matrix (ref channel = `+1`, alt channel = `−1` at mutation position).
- Splits into 80 % train / 20 % test (random seed 42).
- Builds and trains the attention regressor (architecture below).
- Uses a custom focal regression loss (γ = 2) to down-weight easy predictions.
- Also evaluates the trained model on an independent lab-evolution dataset (`Lab_ltte_with_af1.csv`).
- Saves weights to `attn_regressor_focal.weights.h5` and training history to `history.pkl`.

### `fig5a.ipynb`
- Constructs a worked one-hot encoding example for a representative 19-mer.
- Visualises the reference sequence and the ref − alt difference vector using `logomaker`.
- Produces **Fig. 5a** (encoding schematic / sequence-logo panel).

### `fig5bcd.ipynb`
- Loads `history.pkl` written by the training notebook.
- Builds three training-curve plots via R (rpy2 + ggplot2) and saves them as PDFs.
- **Fig. 5b** → focal loss curve → `Transformer_loss_curve_custom.pdf`
- **Fig. 5c** → MAE curve → `Transformer_MAE_loss_curve.pdf`
- **Fig. 5d** → MSE curve → `Transformer_MSE_loss_curve.pdf`

### `fig5ef.ipynb`
- Reloads `attn_regressor_focal.weights.h5` into the attention model built with `return_attention=True`.
- Re-encodes and re-splits the LTEE dataset (same seed 42) to recover the held-out test set.
- Computes R², MAE, MSE, and Pearson correlation on the test set.
- **Fig. 5e** → actual vs predicted kGain scatter (hex-bin, via R/ggplot2) → `LTEE_actual_vs_predicted.pdf`
- **Fig. 5f** → prediction-error KDE → `LTEE_custom_error_plot.pdf`

## Data inputs

Expected under `data/` (two levels up from this folder):

- `kgain_all_population_wt.csv` — LTEE mutations with 19-mer flanking sequence and wild-type kGain values
- `Lab_ltte_with_af1.csv` — independent lab-evolution dataset for out-of-sample evaluation

## Key parameters

| Category | Parameter | Value |
|---|---|---|
| **Input Representation** | Sequence length | 19 nucleotides (2k−1, k=10) |
| | Input encoding | Dual encoded one-hot (A, C, G, T) |
| | Input dimension | 4 channels |
| **Model Architecture** | Embedding dimension | 32 |
| | Positional encoding | Learned (trainable embedding) |
| | Transformer layers | 2 |
| | Attention heads | 4 |
| | Feed-forward hidden dimension | 32 |
| | Feed-forward activation | ReLU |
| | Pooling | Global average pooling |
| | Output activation | Linear (single neuron) |
| | Dropout rate | 0.1 |
| | Weight initialization | Glorot uniform |
| **Training** | Optimizer | Adam |
| | Learning rate | 0.01 |
| | Loss function | Focal regression loss (γ = 2) |
| | Batch size | 8,192 |
| | Epochs | 500 |
| | Train / Validation split | 80% / 20% |
| | Random seed | 19 |
| **Hyperparameter Selection** | Method | Iterative manual adjustment |
| | Criterion | Validation loss minimization |
| | Test set | Held-out (n = 6,879), reserved for final evaluation only |

Weights saved in `attn_regressor_focal.weights.h5`. Model config in `model_config.json`.

## Panel map

- **Fig. 5a** — one-hot encoding schematic of a 19-mer mutation window (sequence logo + heatmap)
- **Fig. 5b** — epoch-wise focal regression loss (train vs validation)
- **Fig. 5c** — epoch-wise MAE (train vs validation)
- **Fig. 5d** — epoch-wise MSE (train vs validation)
- **Fig. 5e** — actual vs predicted kGain on held-out LTEE test set (R², MAE, Pearson r annotated)
- **Fig. 5f** — prediction-error KDE (Gaussian KDE of residuals, mean line overlaid)

## Outputs

| File | Produced by |
|---|---|
| `attn_regressor_focal.weights.h5` | `transformer_clean_final_updated.ipynb` |
| `history.pkl` | `transformer_clean_final_updated.ipynb` |
| `epoch_wise_loss.csv` | `transformer_clean_final_updated.ipynb` |
| `Transformer_loss_curve_custom.pdf` | `fig5bcd.ipynb` |
| `Transformer_MAE_loss_curve.pdf` | `fig5bcd.ipynb` |
| `Transformer_MSE_loss_curve.pdf` | `fig5bcd.ipynb` |
| `LTEE_actual_vs_predicted.pdf` | `fig5ef.ipynb` |
| `LTEE_custom_error_plot.pdf` | `fig5ef.ipynb` |
| `custom_error_kde.pdf` | `fig5ef.ipynb` |

## Recommended execution order

```
1. transformer_clean_final_updated.ipynb   # train model, save weights + history
2. fig5a.ipynb                             # encoding schematic (no model needed)
3. fig5bcd.ipynb                           # loss curves (requires history.pkl)
4. fig5ef.ipynb                            # evaluation panels (requires weights + data)
```

## Dependencies

- `numpy`, `pandas`, `scipy`, `scikit-learn`
- `matplotlib`, `seaborn`
- `tensorflow` / `keras`
- `logomaker` (Fig. 5a only)
- `rpy2` with R packages `ggplot2`, `hexbin` (Fig. 5b–f plotting)
- `Pillow`, `pdf2image` (optional, for inline PDF display)
- local helper: `utility.py`

## Notes

- `fig5ef.ipynb` rebuilds the model graph with `return_attention=True` before loading weights; it does **not** re-train.
- The training notebook caps GPU memory to 40 000 MB by default; lower this value if running on a smaller GPU.
- Project root is discovered dynamically by walking up to `utility.py`, so notebooks run correctly from any subdirectory.
