# What part of "proliferation" should the adjustment remove?

**D157, 2026-11-20.** Prompted by the user's question: *"shouldn't we subtract the averaged out
proliferation signal so that we don't subtract the full proliferation that is specific"*.

The answer is yes, and the code was doing the opposite of what its own header claimed.

## The three layers

A drug's effect on proliferation genes (HALLMARK_E2F_TARGETS + HALLMARK_G2M_CHECKPOINT, 303 genes
on the LINCS profiles) decomposes into three nested pieces:

| layer | what it is | should it be removed? |
|---|---|---|
| **generic** | what *any* drug does to this cell line at this timepoint — assay background, handling, medium | **yes** — it is not about the drug |
| **drug-universal** | this drug arrests growth in *every* line it was tested in — a general cytostatic | **yes** — this is caveat C4, the reason growth-signalling inhibitors rank top in every disease |
| **line-selective** | this drug arrests growth in *this* lineage and not others | **no** — this is a candidate mechanism, and under the discovery goal it is a finding |

The adjustment `A_adj = A − β·π` was written against layer 1. It is applied to layers 2 **and** 3.

## Why layer 1 was already gone

Production runs `PROFILE=specific`, and `79i:203` reads `π` off the **specific** profiles — the ones
`79c2` builds by projecting out the generic drug response per `(ModelID, pert_time)` across a
≥30-drug mean. The generic component is removed at the **gene level**, before `π` is computed.

Measured over the 26,130 condition-level profiles:

| | early (≤6 h) | late (≥24 h) |
|---|---|---|
| mean `pi_raw` | −0.279 (sd 0.355) | −0.868 (sd 0.763) |
| mean `pi_spec` | **−0.018** (sd 0.261) | **−0.014** (sd 0.534) |

`corr(pi_spec, group mean of pi_raw) = **0.026**`. `pi_spec` is mean-zero and carries none of the
generic signal — it *is* the drug-specific deviation.

Variance of `pi_raw` decomposed over `(ModelID, pert_time)`:

```
total          0.500
between-group  0.112   22.4%   <- generic, already removed by 79c2
within-group   0.388   77.6%   <- drug-specific, what 79i has been subtracting
```

Mean |pi_raw| 0.701, |group mean| 0.660, |pi_spec| 0.322.

## The split that matters

Within each window, split `pi_spec` by the drug's cross-line mean:

| window | conditions | drugs | total var | drug-universal | line-selective |
|---|---|---|---|---|---|
| early | 9,225 | 3,062 | 0.0683 | 0.0353 (**51.6%**) | 0.0330 (**48.4%**) |
| late | 16,905 | 2,609 | 0.2849 | 0.1153 (**40.5%**) | 0.1696 (**59.5%**) |

About half of what the adjustment discards is a drug arresting growth *selectively*.

## The change

`PROLIF_MODE` in `scripts/79i_split_reversal.py`:

- `full` — subtract `β·π`. The pre-D157 behaviour, still the default. Nothing moved.
- `universal` — subtract `β · mean(π over that drug's lines, within the window)`.
- `none` — no adjustment.

`β` is fitted on the **full** `π` in every mode, so the per-unit slope keeps its statistical power;
only the component it is applied to changes. This is a partial residualization, not a reduced fit.

Measured on BRCA (815 regulators × 26,130 conditions, per-window β per D154):

```
mean |A_adj - A|   full       0.261
                   universal  0.175      (removes 67.1% of what full removes)
mean |universal - full|       0.180      (max 33.7)
```

## Caveats

1. **`universal` degenerates to `full` for single-line drugs**, where the drug-mean is the value
   itself. Median lines per drug is **1 early / 2 late**. Because multi-line drugs carry more
   conditions the affected share is **20.8% of early and 8.7% of late** conditions. `79i` prints
   this count on every `universal` run.
2. `corr(pi_spec, pi_raw − group mean) = 0.50`, not 1.0. `79c2` projects out a full generic **gene
   vector**; this check centres the scalar `π`. Related operations, not the same one — so the
   22.4% / 77.6% split describes `pi_raw`, not exactly what `79c2` removed.
3. "Universal" is a share of variance, not a clean class. Only **2 drugs early and 8 late** (of those
   in ≥4 lines) have mean `π < −0.3` in *every* line they were profiled in.
4. The risk in the other direction is real and unresolved: removing less proliferation reinstates
   part of caveat C4. The three modes exist so the recovery rate can be measured across all of them
   rather than argued.

## Relationship to D154

D154 fixed **when** `β` is fitted (per window — a pooled β changed sign for 25.2% of regulators).
D157 fixes **what it is applied to**. They are independent and compose: `WINDOWS=1 PROLIF_MODE=universal`.
