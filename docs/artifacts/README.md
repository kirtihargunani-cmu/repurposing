# Published pages and decks (copies in git)

Snapshots of the claude.ai artifacts used to present this work, saved 2026-11-23 so they can be read without a
claude.ai account. The live versions are private to the project owner; the links are listed for the record.

| What | Folder | Live version | Updated |
|---|---|---|---|
| **FDA Gene-Level Rankings** (production ranking, D211/D212): per cancer, D/A/O target evidence, early and late gene-level reversal, target annotations | [`../../results/lincs/gene_rankings_page/`](../../results/lincs/gene_rankings_page/) (built by `scripts/79as_gene_ranking_page.py`) | https://claude.ai/artifact/86j6Ns4BozZ6TrjZcHgrxm | 2026-11-23 |
| Drug Repurposing: Validation & Confounders (slides) | [`validation-and-confounders/deck.html`](validation-and-confounders/deck.html) | https://claude.ai/artifact/DuanJ5BH4pY1yjifM89nJn | 2026-11-23 |
| Three Pipeline Validations (slides) | [`three-pipeline-validations/deck.html`](three-pipeline-validations/deck.html) | https://claude.ai/artifact/8ATMZPWoJS1GVkjnxgE53v | 2026-11-23 |
| Gene-Level Reversal Scoring (slides) | [`gene-level-reversal-scoring-deck/deck.html`](gene-level-reversal-scoring-deck/deck.html) | https://claude.ai/artifact/5SjKC2fFxHAsx2FfrvpErh | 2026-11-23 |
| Gene-Level Reversal Scoring (page version) | [`gene-level-reversal-scoring-page/index.html`](gene-level-reversal-scoring-page/index.html) | https://claude.ai/artifact/6dtV3UGv8iqSS5Q41HLqCQ | 2026-11-23 |
| Reversal Pipeline Walkthrough: data in, the math, data out | [`reversal-pipeline-walkthrough/index.html`](reversal-pipeline-walkthrough/index.html) (built by `scripts/79x3_pipeline_explainer.py`) | https://claude.ai/artifact/8NzPDNqdWuj2mZNvvzDoqC | 2026-11-23 |

**Opening them.** The decks and the two pages open directly in a browser. Each deck folder keeps the deck's own
files (`project/deck.json`, `deck.css`, `slides/*.html`); `deck.html` is assembled from them by
`scripts/deck_to_html.py` and shows each slide with its speaker notes underneath. The FDA rankings page loads one
data file per cancer, which browsers refuse from a local file, so serve the folder:

```
cd results/lincs/gene_rankings_page && python -m http.server 8000   # then open http://localhost:8000
```

**Read with the dates in mind.** These are snapshots. The "Validation & Confounders" deck quotes D127 as 78% correct
(49 instances); the current figure is 36/49 per compound and 16/20 per target unit on the lead score (D208, D209),
and the D211 gene-level ranking does not pass it (C30 in `docs/notes/drug_side_decisions_and_caveats.md`).

Not included: the two Claude Docs specifications (disease side Steps 1-7, drug side), which sit in a separate
claude.ai organization that this export could not reach. Their content was mirrored into
[`../PIPELINE_V2.md`](../PIPELINE_V2.md), [`../PIPELINE_DRUG_SIDE_V2.md`](../PIPELINE_DRUG_SIDE_V2.md) and the two
`*_PLAIN_LANGUAGE.md` files, which are current as of 2026-11-23.
