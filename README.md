# repurposing

Drug-repurposing project, fresh start 2026-10-01. Plan: `docs/PROJECT_PLAN_steps1-4.md`; handoff and decisions: `docs/PROJECT_SETUP_NOTES.md`.

**Current status:** [`docs/PIPELINE_SUMMARY.md`](docs/PIPELINE_SUMMARY.md) (flow chart and per-group numbers) · progress log [`docs/PROGRESS.md`](docs/PROGRESS.md).

- `config/datasets.yaml`: pinned sources (release, URL, access, checksum)
- `data/raw`: untouched downloads · `data/interim` · `data/processed`
- `scripts/`, `results/`, `logs/`

Gotchas: Toil = GENCODE v23, recount3 = v26; join on unversioned Ensembl IDs. GTEx whole blood is not a marrow reference.
