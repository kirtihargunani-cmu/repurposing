# repurposing

Drug-repurposing project, fresh start 2026-10-01. Plan: `docs/PROJECT_PLAN_steps1-4.md`; handoff and decisions: `docs/PROJECT_SETUP_NOTES.md`.

**Start here (v2, 2026-10-29):** [`docs/PIPELINE_V2.md`](docs/PIPELINE_V2.md) explains every step: what it does and why, the command, inputs, outputs and decisions. Run steps with `scripts/run_v2.sh` (`list` / `show` / step names); check outputs with `Rscript scripts/00_status.R`. Drug side: branch `celligner`, `docs/PIPELINE_DRUG_SIDE_V2.md` and `scripts/run_v2_drug.sh`.

Also: decisions [`docs/DECISIONS.md`](docs/DECISIONS.md) · progress log [`docs/PROGRESS.md`](docs/PROGRESS.md) · v2 rebuild spec [`docs/SPEC_rebuild_v2.md`](docs/SPEC_rebuild_v2.md) · reporting requirements [`docs/LLM_COORDINATION_NOTES.md`](docs/LLM_COORDINATION_NOTES.md) · v2 glossary [`docs/GLOSSARY_V2.md`](docs/GLOSSARY_V2.md) · v1 snapshot [`docs/PIPELINE_SUMMARY.md`](docs/PIPELINE_SUMMARY.md).

- `config/datasets.yaml`: pinned sources (release, URL, access, checksum)
- `data/raw`: untouched downloads · `data/interim` · `data/processed`
- `scripts/`, `results/`, `logs/`

Gotchas: Toil = GENCODE v23, recount3 = v26; join on unversioned Ensembl IDs. GTEx whole blood is not a marrow reference.
