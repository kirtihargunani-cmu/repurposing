# Storyboard: animated walk-through of the repurposing pipeline

Fill-in template. One row per **scene** (one Manim class, one slide group). Each scene has 2-4 **beats**: one idea each,
ending in a pause the speaker clicks past. Pre-filled from the docs on 2026-11-02; numbers marked `DATA.x` come from
`scene_data.json` (run `export_scene_data.py`), others from `main:docs/PIPELINE_V2.md` and
`docs/PIPELINE_DRUG_SIDE_V2.md`.

**Recurring visual language:** everything stays consistent from scene to scene.
- **A regulator** is a small labeled circle. Its activity is a vertical bar: warm (up in the tumor) or cool (down).
- **A patient** is a dot. **A cell line** is a hollow square. **A drug** is a yellow pill shape.
- **The "healthy" baseline** is a horizontal line at 0. "Reversal" always means bars moving toward that line.
- **One running example** throughout: lung adenocarcinoma (TCGA-LUAD) + erlotinib. Luminal breast + anti-estrogens is
  the honest failure case.

| # | Scene (class) | Question it answers | Visual metaphor | Beats | Data | Takeaway caption |
|---|---|---|---|---|---|---|
| 0 | `Opening` | What is this about? | Title over a faint network | 1 title · 2 "41 cancers · 14 subgroups · 3,058 drugs" | `pipeline_facts` | Find drugs that push a cancer back toward healthy |
| 1 | `CohortsAndGenes` | What data do we start from? | Dots (patients) pour into 49 bins | 1 10,211 patients into 49 groups · 2 the 17,784-gene filter: hard-to-map genes fade out | `pipeline_facts` | Tumors plus healthy tissue (GTEx), same genes for all |
| 2 | `TumorVsHealthy` | What changed in the tumor? | Two clouds (tumor, healthy), a gene's bar grows between them | 1 a single gene's difference · 2 repeat for all genes: a signature · 3 corrections (purity, collection batch) shrink fake differences | fill from PIPELINE_V2 step 4 | A signature = how every gene differs from healthy |
| 3 | `RegulatorNetworks` | Who controls those genes? | A regulator circle grows edges to its target genes | 1 ARACNe learns targets from co-variation · 2 29 networks, TF only · 3 cross-mappable edges cut | `pipeline_facts.aracne_networks` | Each regulator has a set of target genes |
| 4 | `RegulatorActivity` | Which regulators are driving the tumor? | Targets light up; their regulator's bar rises (VIPER) | 1 targets up together → regulator active · 2 a ranked list of regulators per cancer (NES) | worked example `disease_nes` | Master regulators = the cancer's control switches |
| 5 | `Subgroups` | Are patients within a cancer different? | Patient dots cluster, split into colored subgroups | 1 clustering · 2 keep rule (prediction strength >= 0.75) · 3 14 subgroups survive | `pipeline_facts.subgroups` | 14 subgroups get their own regulators |
| 6 | `CellLineMatching` | Which lab cell lines look like these tumors? | Tumors (dots) and lines (squares) are aligned; nearest lines glow, weight bars appear | 1 before/after Celligner alignment · 2 same-tissue rule · 3 LUAD top lines with weights | `celligner_luad_top_lines` | A drug's lab data counts only if the line resembles the tumor |
| 7 | `DrugProfiles` | What does a drug do to a cell? | Pill drops onto a line; gene bars change; the line's "average drug response" subtracted | 1 LINCS signature · 2 remove the generic "every drug slows growth" response · 3 borrowing from other tissues | `pipeline_facts.lincs_drugs`, `low_transfer_share_pct` | Remove the response every drug causes, keep what's specific |
| 8 | `DrugRegulators` | Which switches does the drug flip? | Same regulator circles as scene 4, now yellow arrows | 1 VIPER on the drug profile · 2 two maps (ARACNe, CollecTRI) averaged · 3 why: ESR1 misread example | `controls` (anti-estrogens) | Read the drug through the same regulators as the disease |
| 9 | `ReversalScore` *(template done)* | Does the drug undo the disease? | Disease bars + opposite yellow arrows → running total | 1 disease bars · 2 drug arrows · 3 contributions sum to rev 22.7 (647 of 1,061) | `worked_examples["TCGA-LUAD\|erlotinib"]` | No single gene carries the score |
| 10 | `Specificity` | Is it special to this cancer? | The drug's rev across 41 cancers as a strip; this cancer pulled out | 1 rev in every cancer · 2 subtract the average → spec · 3 32 drugs top >= 20 cancers by rev; 0 by spec_mech | `controls`, `generic_drugs_by_rev` | A drug that "fixes" everything is just toxic |
| 11 | `ProliferationSplit` | Is it more than stopping division? | Score bar splits into a "division" part and a "mechanism" part | 1 split for palbociclib (mostly division) vs erlotinib | `worked_examples["BRCA-LumB\|palbociclib"]` | We lead with the mechanism part (spec_mech) |
| 12 | `TargetAndFlags` | Can we trust this drug-cancer pair? | Checklist ticks on a drug card: target expressed, lines match, transfer OK | 1 target gate · 2 unmatched_lines (basal breast anti-estrogens) · 3 low_transfer | `controls` (BRCA-Basal rows) | Flagged drugs stay listed but never reach the top tiers |
| 13 | `EvidenceStatus` | How good is the data per cancer? | 41 tiles fill green (matched) or grey (borrowed only); then the stricter cell-type check filters them | 1 14 matched / 27 borrowed · 2 data-gap funnel D1 → D2 → D1b → D2b → 14 pass | `evidence_status_counts`, `data_gap_counts` | Say "no good-quality data" instead of a list |
| 14 | `Results` | What comes out? | Cards per cancer with top 3 FDA drugs; known drugs highlighted | 1 known-drug checks (erlotinib 96th, flutamide 85th, irinotecan 98.5th) · 2 the misses · 3 subgroup example (MSI colon → MEK inhibitors) | `controls`, `fda_top3`, `fda_top3_subgroups` | Hypotheses ranked, not yet validated |
| 15 | `WhatsNext` | What's left? | Checklist | 1 PRISM / GDSC validation (on hold) · 2 cell-type-only scoring · 3 review recurring oddities | — | Next: test against measured drug killing |

## Per-scene template (copy for each scene while drafting)

```
### Scene <#>: <ClassName>
Question:            <one question, in the audience's words>
Opening frame:       <what is on screen when the scene starts; reuse objects from the previous scene>
Beat 1:              <what moves> → narration: "<one sentence>"
Beat 2:              ...
Beat 3:              ...
Closing frame:       <what stays for the next scene>
Data keys:           <scene_data.json keys>
Takeaway caption:    <≤ 12 words>
Honesty note:        <caveat to say aloud, if any (C-number)>
```
