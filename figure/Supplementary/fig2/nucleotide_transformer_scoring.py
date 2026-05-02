"""
nucleotide_transformer_scoring.py

Score point mutations using Nucleotide Transformer with a masked-marginal
log-likelihood ratio at the mutation position.

For each SNV, we compute:

    LLR(alt) = log P(alt | masked context) - log P(ref | masked context)

where the 6-mer token(s) overlapping the mutation site are masked, the model
predicts a distribution over its full 6-mer vocabulary at each masked
position, and we sum the probabilities of all 6-mers compatible with the
candidate base at the mutation position. This gives a base-resolution LLR
that depends only on the mutation site, unlike a whole-sequence
pseudo-likelihood (which dilutes the per-position signal across all
sequence-identical tokens).

The implementation also retains a fallback to the whole-sequence
pseudo-likelihood approach for tokenizers without 6-mer-style vocabularies
(e.g., single-nucleotide-resolution variants), but issues a warning so the
caller knows the score is being computed under the fallback.

Sign convention: positive LLR -> alt more probable than ref under the model.
This matches kGain's convention (positive kGain -> alt context more common
in the within-genome k-mer background).
"""
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from Bio import SeqIO
from tqdm import tqdm
from transformers import AutoModelForMaskedLM, AutoTokenizer

# --------------------------------------------------------------------------
# Globals (model cache)
# --------------------------------------------------------------------------
_model = None
_tokenizer = None
_device = None
_token_to_base_indicator = None  # cached per-tokenizer index for fast LLR


# --------------------------------------------------------------------------
# Remote-code patching (handles broken cached modeling_esm.py from InstaDeep)
# --------------------------------------------------------------------------
def _patch_cached_remote_esm():
    """Patch known syntax / import bugs in cached InstaDeep modeling_esm.py."""
    base = Path.home() / ".cache/huggingface/modules/transformers_modules/InstaDeepAI"
    if not base.exists():
        return False
    patched = False
    for f in base.rglob("modeling_esm.py"):
        src = f.read_text()
        orig = src

        # Fix 1: unescaped double quotes inside double-quoted string (~line 1125)
        src = src.replace(
            '`getattr(config, "is_decoder", False)=False`',
            "`getattr(config, 'is_decoder', False)=False`",
        )

        # Fix 2: find_pruneable_heads_and_indices / prune_linear_layer moved
        # from transformers.modeling_utils -> transformers.pytorch_utils in >=4.34
        if "find_pruneable_heads_and_indices" in src and \
           "from transformers.modeling_utils import (" in src:
            def _strip_block(m):
                block = m.group(0)
                for sym in ("find_pruneable_heads_and_indices", "prune_linear_layer"):
                    block = re.sub(rf"^\s*{sym},?\s*\n", "", block, flags=re.MULTILINE)
                return block
            src = re.sub(
                r"from transformers\.modeling_utils import \([^)]*\)",
                _strip_block,
                src,
                count=1,
            )
            if "from transformers.pytorch_utils import" not in src:
                src = (
                    "from transformers.pytorch_utils import "
                    "find_pruneable_heads_and_indices, prune_linear_layer\n"
                    + src
                )

        if src != orig:
            f.write_text(src)
            patched = True
    return patched


# --------------------------------------------------------------------------
# Model loading
# --------------------------------------------------------------------------
def _try_load(model_name, revision, hf_token, trust_remote_code):
    kwargs = dict(trust_remote_code=trust_remote_code)
    if revision:
        kwargs["revision"] = revision
    if hf_token:
        kwargs["token"] = hf_token
    tok = AutoTokenizer.from_pretrained(model_name, **kwargs)
    mdl = AutoModelForMaskedLM.from_pretrained(model_name, **kwargs)
    return tok, mdl


def _load_model(model_name, revision=None, hf_token=None):
    """Load Nucleotide Transformer with multi-stage fallback."""
    global _model, _tokenizer, _device, _token_to_base_indicator

    if _model is not None and _tokenizer is not None:
        return _model, _tokenizer, _device

    _device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Loading {model_name} on {_device} ...")
    if revision:
        print(f"[Info] Using pinned revision: {revision}")

    last_err = None

    # Stage 1: built-in ESM (no remote code)
    try:
        _tokenizer, _model = _try_load(model_name, revision, hf_token, trust_remote_code=False)
        _model.to(_device).eval()
        print("[Info] Loaded with trust_remote_code=False (built-in ESM).")
    except Exception as e:
        last_err = e
        print(f"[Stage 1 failed] {type(e).__name__}: {e}")

        # Stage 2: remote code, raw
        try:
            _tokenizer, _model = _try_load(model_name, revision, hf_token, trust_remote_code=True)
            _model.to(_device).eval()
            print("[Info] Loaded with trust_remote_code=True.")
        except (SyntaxError, ImportError) as e:
            last_err = e
            print(f"[Stage 2 failed] {type(e).__name__}: {e}")
            print("Patching cached remote code ...")
            _patch_cached_remote_esm()

            # Stage 3: remote code after patch
            try:
                _tokenizer, _model = _try_load(model_name, revision, hf_token, trust_remote_code=True)
                _model.to(_device).eval()
                print("[Info] Loaded with trust_remote_code=True after patching.")
            except Exception as e2:
                last_err = e2
                print(f"[Stage 3 failed] {type(e2).__name__}: {e2}")
                msg = (
                    f"Failed to load NucleotideTransformer model.\n"
                    f"Model id attempted: {model_name}\n"
                    f"Original error: {type(last_err).__name__}: {last_err}"
                )
                raise RuntimeError(msg) from last_err
        except Exception as e:
            last_err = e
            msg = (
                f"Failed to load NucleotideTransformer model.\n"
                f"Model id attempted: {model_name}\n"
                f"Original error: {type(last_err).__name__}: {last_err}"
            )
            raise RuntimeError(msg) from last_err

    # Pre-compute token <-> nucleotide indicator once per tokenizer
    _token_to_base_indicator = _build_token_base_indicator(_tokenizer)

    return _model, _tokenizer, _device


# --------------------------------------------------------------------------
# Build a (vocab_size, 4, k) indicator tensor:
#     indicator[token_id, base_idx, pos_in_kmer] = 1
#     iff the k-mer string represented by token_id has nucleotide base_idx at pos_in_kmer
# This lets us, given a masked position and a candidate base, sum the
# log-probabilities of all 6-mer tokens whose `pos_in_kmer`-th character
# equals that base. That sum is exactly the marginal P(base | masked context).
# --------------------------------------------------------------------------
_BASES = ("A", "C", "G", "T")
_BASE_TO_IDX = {b: i for i, b in enumerate(_BASES)}


def _build_token_base_indicator(tokenizer):
    """
    Return a tensor of shape (vocab_size, 4, k_max) on CPU, where
        indicator[t, b, j] = 1.0  iff token t is a length-(>=j+1) DNA k-mer
                                  whose j-th character equals base b.
    For non-DNA tokens (special tokens, padding, etc.), indicator is all zero.
    The k-mer length is inferred per token from its string form; k_max is the
    longest pure-DNA token in the vocabulary.
    """
    vocab = tokenizer.get_vocab()

    # Determine k-mer length(s) present in vocabulary
    dna_tokens = {}
    for tok, idx in vocab.items():
        # Strip any leading meta-character (e.g., "â" prefix some tokenizers add)
        clean = tok.lstrip("Ġ▁_ ")
        if len(clean) >= 1 and all(c in "ACGT" for c in clean):
            dna_tokens[idx] = clean

    if not dna_tokens:
        warnings.warn(
            "Tokenizer vocabulary contains no pure ACGT k-mer tokens; "
            "masked-marginal LLR will fall back to whole-sequence pseudo-likelihood.",
            stacklevel=2,
        )
        return None

    k_max = max(len(s) for s in dna_tokens.values())
    vocab_size = len(vocab)

    indicator = torch.zeros(vocab_size, 4, k_max, dtype=torch.float32)
    for idx, kmer in dna_tokens.items():
        for j, ch in enumerate(kmer):
            indicator[idx, _BASE_TO_IDX[ch], j] = 1.0

    return indicator


# --------------------------------------------------------------------------
# FASTA / sequence helpers
# --------------------------------------------------------------------------
def _load_fasta(fasta_path):
    """Load first record from FASTA. Returns uppercase sequence string."""
    rec = next(SeqIO.parse(str(fasta_path), "fasta"))
    seq = str(rec.seq).upper()
    print(f"Loaded FASTA: {len(seq):,} bp")
    return seq


def _extract_context(seq, pos0, flank):
    """
    Extract a window of length 2*flank + 1 centered at 0-based position pos0.
    Pads with 'N' if near edges.
    """
    L = len(seq)
    start = pos0 - flank
    end = pos0 + flank + 1  # exclusive
    left_pad = max(0, -start)
    right_pad = max(0, end - L)
    s = seq[max(0, start):min(L, end)]
    if left_pad or right_pad:
        s = ("N" * left_pad) + s + ("N" * right_pad)
    return s


# --------------------------------------------------------------------------
# Masked-marginal LLR (preferred)
# --------------------------------------------------------------------------
@torch.no_grad()
def _score_batch_masked_marginal(
    contexts, ref_bases, alt_bases, mask_idx_in_seq,
    model, tokenizer, device, indicator,
):
    """
    Proper masked-marginal LLR.

    For each (context, ref, alt):
      1. Tokenize the unmodified context.
      2. Identify the token(s) whose underlying nucleotide span covers
         `mask_idx_in_seq` (the central / mutation position).
      3. For each such token, replace it with [MASK] and run the model.
         The model produces a probability distribution over the full vocab
         at the masked position.
      4. Aggregate vocab probabilities into per-base probabilities by summing
         tokens whose j-th character equals each base, where j is the offset
         of the mutation site inside that token.
      5. LLR = log P(alt | masked) - log P(ref | masked), summed over all
         tokens that overlap the mutation site (typically one token for
         single-nucleotide-resolution tokenizers, up to k tokens for
         k-mer overlapping tokenizers).

    Returns one float per input mutation.
    """
    if indicator is None:
        # Fallback for tokenizers without DNA k-mer tokens
        return _score_batch_pseudo_likelihood(
            contexts, ref_bases, alt_bases, mask_idx_in_seq,
            model, tokenizer, device,
        )

    indicator = indicator.to(device)
    k_max = indicator.shape[2]
    mask_id = tokenizer.mask_token_id
    if mask_id is None:
        raise RuntimeError("Tokenizer has no mask token; cannot do masked-marginal LLR.")

    out = np.full(len(contexts), np.nan, dtype=np.float64)

    for i, (ctx, ref_b, alt_b) in enumerate(zip(contexts, ref_bases, alt_bases)):
        ref_idx = _BASE_TO_IDX.get(ref_b)
        alt_idx = _BASE_TO_IDX.get(alt_b)
        if ref_idx is None or alt_idx is None:
            continue

        # Tokenize raw context
        enc = tokenizer(ctx, return_tensors="pt", add_special_tokens=True)
        input_ids = enc["input_ids"].to(device)
        token_strs = tokenizer.convert_ids_to_tokens(input_ids[0].tolist())

        # Walk along token strings to map each token -> the slice of
        # input characters it covers (excluding special tokens which contribute 0 chars)
        char_pos = 0
        token_spans = []  # list of (token_index, char_start_inclusive, char_end_exclusive)
        for ti, tok_str in enumerate(token_strs):
            clean = tok_str.lstrip("Ġ▁_ ")
            if len(clean) >= 1 and all(c in "ACGT" for c in clean):
                token_spans.append((ti, char_pos, char_pos + len(clean)))
                char_pos += len(clean)
            # else: special token, contributes 0 chars

        # Find tokens overlapping mask_idx_in_seq
        overlapping = [
            (ti, cs, ce) for (ti, cs, ce) in token_spans
            if cs <= mask_idx_in_seq < ce
        ]
        if not overlapping:
            continue  # mask position outside any DNA token; skip

        llr = 0.0
        valid = True
        for (ti, cs, ce) in overlapping:
            j = mask_idx_in_seq - cs  # offset of mutation site within this k-mer token
            if j < 0 or j >= k_max:
                valid = False
                break

            # Mask this token, run model
            ids_masked = input_ids.clone()
            ids_masked[0, ti] = mask_id
            logits = model(ids_masked).logits  # (1, T, V)
            pos_logits = logits[0, ti]  # (V,)

            # log P(base | masked) = logsumexp over tokens whose j-th char == base,
            # minus log partition (logsumexp over all tokens).
            # Using boolean indexing to exclude non-matching tokens cleanly.
            log_Z = torch.logsumexp(pos_logits, dim=0)  # log partition
            ind_j = indicator[:, :, j].to(device)  # (V, 4)
            log_base_probs = torch.full((4,), float('-inf'), device=device)
            for b in range(4):
                mask = ind_j[:, b].bool()
                if mask.any():
                    log_base_probs[b] = torch.logsumexp(pos_logits[mask], dim=0) - log_Z

            llr += float(log_base_probs[alt_idx] - log_base_probs[ref_idx])

        if valid:
            out[i] = llr

    return out


# --------------------------------------------------------------------------
# Whole-sequence pseudo-likelihood fallback
# (kept for tokenizers without DNA k-mer vocab; clearly labeled)
# --------------------------------------------------------------------------
@torch.no_grad()
def _score_batch_pseudo_likelihood(
    contexts, ref_bases, alt_bases, mask_idx_in_seq,
    model, tokenizer, device,
):
    """
    Fallback: whole-sequence pseudo-likelihood difference.

    NOTE: This score is NOT a per-position LLR. It is the difference of mean
    per-token pseudo-log-likelihoods between the alt-substituted and
    ref-substituted full-window sequences. Its magnitude is diluted relative
    to a true masked-marginal LLR, especially for long context windows.
    Use only when the tokenizer does not expose a DNA k-mer vocabulary.
    """
    warnings.warn(
        "Falling back to whole-sequence pseudo-likelihood scoring. "
        "Score is NOT a per-position masked-marginal LLR.",
        stacklevel=2,
    )

    ref_seqs, alt_seqs = [], []
    for ctx, r, a in zip(contexts, ref_bases, alt_bases):
        chars = list(ctx)
        ref_chars = chars.copy()
        alt_chars = chars.copy()
        ref_chars[mask_idx_in_seq] = r
        alt_chars[mask_idx_in_seq] = a
        ref_seqs.append("".join(ref_chars))
        alt_seqs.append("".join(alt_chars))

    def _seq_logprob(seqs):
        enc = tokenizer(seqs, return_tensors="pt", padding=True, truncation=True)
        enc = {k: v.to(device) for k, v in enc.items()}
        out = model(**enc)
        logits = out.logits
        log_probs = torch.log_softmax(logits, dim=-1)
        ids = enc["input_ids"]
        attn = enc["attention_mask"]
        gathered = log_probs.gather(2, ids.unsqueeze(-1)).squeeze(-1)
        seq_lp = (gathered * attn).sum(dim=1) / attn.sum(dim=1).clamp(min=1)
        return seq_lp.cpu().numpy()

    ref_lp = _seq_logprob(ref_seqs)
    alt_lp = _seq_logprob(alt_seqs)
    return alt_lp - ref_lp


# --------------------------------------------------------------------------
# Public pipeline
# --------------------------------------------------------------------------
def run_nt_pipeline(
    fasta_path,
    mutations_df,
    out_csv=None,
    flank_size=100,
    model_name="InstaDeepAI/nucleotide-transformer-v2-500m-multi-species",
    revision=None,
    hf_token=None,
    batch_size=8,
    pos_col="position",
    ref_col="ref",
    alt_col="alt",
    scoring="masked_marginal",
):
    """
    Score each mutation in mutations_df with NT log-likelihood ratio.

    Parameters
    ----------
    fasta_path : str | Path
        Path to genome FASTA (first record is used).
    mutations_df : pd.DataFrame
        Must contain columns: position (1-based), ref, alt.
    out_csv : str | Path | None
        If provided, write resulting DataFrame to CSV.
    flank_size : int
        Number of bases on each side of the mutation. Total context length is
        2*flank_size + 1. Use a window comfortably larger than the tokenizer's
        k-mer length (e.g., 50-200 for NT v2 with 6-mers).
    model_name : str
    revision : str | None
    hf_token : str | None
    batch_size : int
        Note: with masked-marginal scoring, each mutation requires its own
        forward pass (one per overlapping token), so batch_size has limited
        effect. Use the pseudo-likelihood fallback for true batched scoring.
    pos_col, ref_col, alt_col : str
        Column names in mutations_df.
    scoring : {"masked_marginal", "pseudo_likelihood"}
        - "masked_marginal" (default, recommended): per-position LLR by
          masking the token(s) overlapping the mutation site.
        - "pseudo_likelihood": legacy whole-sequence score; faster but
          mathematically NOT a per-position LLR.

    Returns
    -------
    pd.DataFrame
        Copy of mutations_df with an added 'nt_llr' column.
    """
    seq = _load_fasta(fasta_path)
    model, tokenizer, device = _load_model(model_name, revision=revision, hf_token=hf_token)
    indicator = _token_to_base_indicator

    if scoring not in {"masked_marginal", "pseudo_likelihood"}:
        raise ValueError(f"Unknown scoring='{scoring}'.")

    df = mutations_df.copy().reset_index(drop=True)

    # Validate
    for c in (pos_col, ref_col, alt_col):
        if c not in df.columns:
            raise ValueError(f"mutations_df missing column: {c}")

    llrs = np.full(len(df), np.nan, dtype=np.float64)
    n = len(df)

    for start in tqdm(range(0, n, batch_size), desc=f"NT scoring ({scoring})"):
        end = min(start + batch_size, n)
        rows = df.iloc[start:end]

        contexts, refs, alts, valid_idx = [], [], [], []
        for i, row in rows.iterrows():
            pos1 = int(row[pos_col])
            r = str(row[ref_col]).upper()
            a = str(row[alt_col]).upper()
            if len(r) != 1 or len(a) != 1 or r not in "ACGT" or a not in "ACGT":
                continue  # skip indels / non-SNVs
            pos0 = pos1 - 1
            if pos0 < 0 or pos0 >= len(seq):
                continue
            ctx = _extract_context(seq, pos0, flank_size)
            contexts.append(ctx)
            refs.append(r)
            alts.append(a)
            valid_idx.append(i)

        if not contexts:
            continue

        if scoring == "masked_marginal":
            scores = _score_batch_masked_marginal(
                contexts, refs, alts,
                mask_idx_in_seq=flank_size,
                model=model, tokenizer=tokenizer, device=device,
                indicator=indicator,
            )
        else:
            scores = _score_batch_pseudo_likelihood(
                contexts, refs, alts,
                mask_idx_in_seq=flank_size,
                model=model, tokenizer=tokenizer, device=device,
            )
        for j, idx in enumerate(valid_idx):
            llrs[idx] = scores[j]

    df["nt_llr"] = llrs
    df["nt_scoring_method"] = scoring  # provenance for reproducibility

    if out_csv is not None:
        out_csv = Path(out_csv)
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)
        print(f"Wrote: {out_csv}")

    return df
