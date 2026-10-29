# repurposing

Drug-repurposing project, fresh start 2026-10-01. Plan: `docs/PROJECT_PLAN_steps1-4.md`; handoff and decisions: `docs/PROJECT_SETUP_NOTES.md`.

**This is the `celligner` branch (drug side). Start here (v2, 2026-10-29):** [`docs/PIPELINE_DRUG_SIDE_V2.md`](docs/PIPELINE_DRUG_SIDE_V2.md) explains every drug-side step: what it does and why, the command, inputs, outputs, decisions and a column guide. Run steps with `scripts/run_v2_drug.sh` after the disease side (`main`: `docs/PIPELINE_V2.md`, `scripts/run_v2.sh`). Decisions and caveats: [`docs/notes/drug_side_decisions_and_caveats.md`](docs/notes/drug_side_decisions_and_caveats.md). Browsable FDA rankings: `results/lincs/fda_rankings_v2.html`.

- `config/datasets.yaml`: pinned sources (release, URL, access, checksum)
- `data/raw`: untouched downloads · `data/interim` · `data/processed`
- `scripts/`, `results/`, `logs/`

Gotchas: Toil = GENCODE v23, recount3 = v26; join on unversioned Ensembl IDs. GTEx whole blood is not a marrow reference.
