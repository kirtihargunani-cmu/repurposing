# Canonical formulas

Every score the pipeline computes, in one place, with the implementation that owns it. **If this
file and the code disagree, the code is right and this file is stale — fix it.**

Written 2026-11-22 because the definitions had accumulated across spec amendments and code
comments, and two of them (`rev_gene`'s sign, the `cos_term` normaliser) had already been stated
wrongly once in prose before being corrected against the implementation.

Notation, fixed throughout:

| symbol | meaning |
|---|---|
| `i` | a regulator (master regulator / TF) |
| `j` | a target gene of regulator `i` |
| `m_j`, `l_j` | that edge's `tfmode` (signed, −1…1) and `likelihood` (0…1), from ARACNe |
| `U` | the scored gene universe (viper `eset.filter`): genes that are a target of some scored regulator |
| `nes_i` | the disease master-regulator NES (tumour vs normal), `group_mr.tsv` |
| `s` | the vector of `nes_i` over the cohort's tested regulators |
| `d_i` | the **disease** gene-contribution vector of regulator `i` (`disease_contrib_nes`) |
| `C_i` | the **drug** gene-contribution vector of regulator `i` (`activity_contrib`) |
| `act_i` | the drug's regulator activity, `= Σ_j C_ij` |

---

## 1. The aREA contribution split (both sides)

The quantity everything else is built on. For a profile, rank its genes over `U`:

```
r_j  = rank_j / (|U| + 1)
t2_j = Φ⁻¹(r_j)                                    signed, the directional term
t1_j = Φ⁻¹(|2 r_j − 1| + (1 − max_j|2 r_j − 1|)/2) unsigned, the magnitude term
```

Per edge `(i, j)`:

```
a_j = m_j · l_j · t2_j
b_j = (1 − |m_j|) · l_j · t1_j
S1_i = Σ_{j∈i} a_j      S2_i = Σ_{j∈i} b_j      σ_i = sign(S1_i)   (0 → +1)

c_ij = ( a_j + σ_i · b_j · [S2_i > 0] ) · k_i
k_i  = w_i / Σ_{j∈i} l_j        w_i = sqrt( Σ_{j∈i} (l_j / max_{j∈i} l_j)² )
```

and `Σ_j c_ij` is exactly the regulator's enrichment score.

- **Drug side** (`79e3_window_contributions.py:Net.contrib`): `C_ij = c_ij`, and the sum **is**
  the condition-level VIPER activity — no null model is involved, so the split is **exact**
  (verified against the stored `.f32`, max |diff| 4.8e-07).
- **Disease side** (`50b_regulon_contributions.R`): msviper applies a permutation null that cannot
  be replayed, so the ES is rescaled to the stored NES **proportionally**:

```
d_ij = c_ij · NES_i / ES_i          (contrib_nes)
```

  This is an allocation, not a re-derivation — stated wherever `d` is used.

## 2. Reversal — the production score

```
rev = −(D · s) / ‖s‖ = −Σ_i nes_i · act_i / ‖s‖
rev_contrib_i = −nes_i · act_i / ‖s‖            (sums to rev)
```

`D` is the drug's regulator-activity vector. Implemented in `79e_drug_reversal.py:285`;
window-resolved in `systems_validation/integrate/window_reversal.py` over conditions and in-key
lines only.

## 3. `rev_gene` — reversal that requires the genes to match (D203)

```
rev            = −Σ_i ( Σ_j d_ij )( Σ_j C_ij ) / ‖s‖      a PRODUCT OF SUMS
rev_gene       = −Σ_i   Σ_j d_ij · C_ij        / ‖s‖      a SUM OF PRODUCTS
rev_gene_contrib_i = −⟨d_i , C_i⟩ / ‖s‖                   (sums to rev_gene)
```

`rev` multiplies the totals, so a regulator reaches a large contribution however it got there.
`rev_gene` pairs the genes, so a regulator whose activity comes from genes orthogonal to the
disease pattern has `⟨d_i, C_i⟩ ≈ 0` even when both totals are large — a down-weighting by
geometry, with **no threshold and no tuned parameter**.

**Sign note, recorded because it was wrong once (D203 amendment 1):** `d_i` already carries the
disease direction, since `Σ_j d_ij = nes_i`. Multiplying by `−nes_i` as well double-counts it and
flips the sign. There is no `nes` factor in `rev_gene`.

`‖s‖` is a per-cohort constant and does not affect within-cohort ranking; it is kept only so the
two scores sit on one scale.

Reported but **not** part of the score:

```
act_gene_i = ⟨d_i , C_i⟩ / ‖d_i‖     the drug's gene vector projected onto the disease pattern,
                                      drug magnitude preserved
```

## 4. Contribution cosine — the opposition view (D194/D195)

```
cos_contrib_i = ⟨d_i , C_i⟩ / (‖d_i‖ ‖C_i‖)            per regulon
cos_term_ij   = d_ij · C_ij / (‖d_i‖ ‖C_i‖)            per gene;  Σ_j cos_term_ij = cos_contrib_i
```

Negative = the drug moves the genes that account for that regulator's dysregulation, in the
opposite direction. Exact to 1.6e-06 over 61,000 regulons.

**This divides by `‖C_i‖` as well, and that is the difference from `rev_gene`.** D202 attributed
the entire ranking loss of the cosine score to that per-regulon, per-drug normalisation:
`rev / ‖D‖` cost **+0.004** AUC, adding the per-regulon normalisation cost **−0.035**. Scale-free
per regulon means a regulon the drug never moved still contributes on the full ±1 scale.

## 5. Programs (D198/D200)

A program is a module of regulators from average-linkage clustering on `1 − |r|`, cut at
`|r| ≥ 0.7` (`58_regulator_modules.R`), over **patients** (disease side) or **drugs** (drug side,
per window).

```
disease_score_P = mean_{i∈P} nes_i
drug_score_P    = mean_{i∈P} act_i
opposition_P    = − disease_score_P · drug_score_P / |disease_score_P|
```

## 6. Windows (D180)

```
early : pert_time ≤ 6 h        late : pert_time ≥ 24 h
other and unrecorded (sentinel −1) are DROPPED, never pooled
```

The unrecorded sentinel must never read as `early`; that latent bug was fixed in
`window_reversal.py` before it could fire.

## 7. Reliability ceiling (D202)

For a drug with ≥ 2 conditions in a window, split them in half, build each half's gene vector, and
take the per-regulon cosine:

```
R = median_i cosine( C_i(half A) , C_i(half B) )
```

`R` bounds every gene-vector comparison — nothing can match a vector better than it matches
itself. Measured **0.531** overall, **0.376 early / 0.590 late**, and strongly sample-size driven
(0.326 at 2 conditions → 0.793 at 10+, Spearman +0.528). Any similarity is reported as a fraction
of `R`, never in absolute terms.

## 8. Approval recovery — reported, never optimised against

AUC of a score against `label_phase ≥ 3`, per (cohort, window), on identical rows. **It is weak
validation**: D192 measured it bimodal (9 cohorts significantly above chance, 9 significantly
below, UCEC-EEC z = −7.7) and D196 showed the bare proliferation axis beats the mechanism score in
COAD/PRAD/LUSC, whose approved sets are cytotoxic. `label_phase = 0` is **not** a negative.
