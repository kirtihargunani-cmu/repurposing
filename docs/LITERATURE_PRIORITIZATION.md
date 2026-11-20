# Prioritizing repurposing candidates for in vivo testing: literature check

Written 2026-11-20, prompted by the user's question: once a causal link from drug target to dysregulated
gene is identified (e.g. by IPA), how do published studies decide which compounds go to animal work?

Abstracts, PMC and publisher pages were read; no benchmark tables were re-derived. **There is no
consensus threshold.** What follows is the set of criteria that recur, and how each maps onto what
this project already computes.

## 1. The only study that tested an IPA link directly

A rare-disease repurposing study ran IPA upstream regulator analysis on differentially expressed genes
from a neuron drug screen, flagged three drug-linked networks — **levothyroxine/DIO3,
hydroxyurea/FOXM1, dexamethasone/PPARD** — and validated by qPCR in vitro and in vivo. Genes tied to
an IPA network validated at **5× and 10× the rate** of individual genes not tied to one.
([PubMed 34208365](https://pubmed.ncbi.nlm.nih.gov/34208365/), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC8231191))

This is the closest published support for the idea that network membership is itself a useful
prioritization filter — i.e. the thing D159 was weighing IPA for.

**The counterweight, and it is substantial.** A 2023 benchmark over 16 datasets found IPA's upstream
regulator analysis ranked the correct causal regulator first in **5 cases and failed to find it at all
in 5 of 16** ([Sci Rep 2023](https://link.springer.com/10.1038/s41598-023-41374-0)). An IPA link is
therefore not self-validating; roughly a third of the time the true regulator is absent from the
output entirely. No study was found that validated IPA-prioritized repurposing candidates in animal
*disease* models or clinical outcomes.

## 2. The recurring criteria, in the order they are usually applied

| # | criterion | example | our equivalent |
|---|---|---|---|
| 1 | reversal / connectivity strength | most studies | `rev`, `rev_pct` |
| 2 | consistency across independent disease signatures | Alzheimer's, aspirin as lead ([PMC13193131](https://pmc.ncbi.nlm.nih.gov/articles/PMC13193131/)) | `spec`, `spec_mech` (contrast vs reference diseases) |
| 3 | reversal confirmed in a **disease-relevant** model, not any line | AML: RNA-seq confirmation before xenograft; reversal held in one line and not the other ([MDPI 2025](https://www.mdpi.com/1424-8247/19/11/1849)) | celligner line matching, `transfer_confidence`, `inactive_matched_line` |
| 4 | minimal off-target transcriptome perturbation | reverse the disease signature while changing the rest as little as possible ([bioRxiv 2020](https://www.biorxiv.org/content/10.1101/2020.05.13.093468v1.full.pdf)) | `spec_mech`, proliferation split (`79i`) |
| 5 | tissue access | GBM: blood–brain-barrier penetration probability ([Transl Med 2024](https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-024-06046-1)) | **not modelled** — route-of-administration flags only (`local_route`) |
| 6 | not already tried in that indication | GBM study filtered drugs with prior GBM trials | the discovery goal; `79o`/`79o2` literature and trial columns |
| 7 | orthogonal benchmarking before committing | signature-reversion study checked predictions against clinical trials and PRISM ([bioRxiv 2023](https://www.biorxiv.org/content/10.1101/2023.03.10.532074v1.full)) | known-use AUC (C28), D127 biomarker pairs; **A5/PRISM deliberately not read** |

## 3. Two places this project is ahead of published practice

**Achievable exposure.** No established method was found for filtering connectivity hits by achievable
plasma Cmax. Searches surfaced signal-strength and TAS filters but nothing on exposure. Our D109b
`low_exposure` / D111 Cmax work and D109a `plate_confined` implement something the literature names as
desirable but rarely does.

**Target-centric rather than reversal-ranked.** The field's central critique of signature reversion is
that it may mostly capture **general anti-proliferative activity** rather than subtype-specific effect
([PMC9851289](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9851289/)). That is our caveat C4 verbatim.
Our lead ranking (`layered_pct`, D122) is built from dependency / alteration / over-activity and
binding rather than from reversal, which sidesteps the critique rather than inheriting it — and D142's
finding that the two are uncorrelated (Spearman −0.0002, top-50 overlap 0/50) means this is a real
difference in what is being ranked, not a relabelling.

## 4. One incidental calibration point

A 2025 study filtered LINCS signatures at a **TAS threshold between 0.2 and 0.5**. D156's adopted 6 h
floor of **0.25** therefore sits inside published practice rather than outside it — useful, since D156
recorded that the original 0.4 had no derivation at all.

## 5. What is missing from the literature, for our purposes

- No consensus numeric threshold for advancing a candidate to in vivo.
- No study validating IPA-prioritized repurposing candidates in animal disease models.
- Reproducibility of connectivity scoring itself is weak: querying CMap 2 with CMap 1-derived
  signatures succeeded 17% of the time, and profiles for the same compound, dose and cell line
  correlate poorly ([PMC8413422](https://pmc.ncbi.nlm.nih.gov/articles/PMC8413422)). Any triage built
  on connectivity rank inherits that.
