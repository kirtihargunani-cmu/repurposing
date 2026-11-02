# repurposing

Drug repurposing from regulator activity: find the transcription factors each cancer turns on or off relative to
healthy tissue (disease side), then drugs whose effect in matching cell lines pushes them back (drug side). v2,
2026-11-02.

**Start here**

| | Technical (commands, inputs, outputs, decisions) | Plain language, with results |
|---|---|---|
| Disease side (Steps 1-10) | [`docs/PIPELINE_V2.md`](docs/PIPELINE_V2.md) | [`docs/DISEASE_SIDE_PLAIN_LANGUAGE.md`](docs/DISEASE_SIDE_PLAIN_LANGUAGE.md) |
| Drug side | [`docs/PIPELINE_DRUG_SIDE_V2.md`](docs/PIPELINE_DRUG_SIDE_V2.md) | [`docs/DRUG_SIDE_PLAIN_LANGUAGE.md`](docs/DRUG_SIDE_PLAIN_LANGUAGE.md) |

Run: `scripts/run_v2.sh` (disease side), then `scripts/run_v2_drug.sh` (drug side); each takes `list`, `show <step>`
or step names. Check outputs with `Rscript scripts/00_status.R`.

Also: decisions and open items [`docs/DECISIONS.md`](docs/DECISIONS.md) · run log [`docs/PROGRESS.md`](docs/PROGRESS.md) ·
drug-side caveats [`docs/notes/drug_side_decisions_and_caveats.md`](docs/notes/drug_side_decisions_and_caveats.md) ·
v2 terms [`docs/GLOSSARY_V2.md`](docs/GLOSSARY_V2.md) · reporting requirements
[`docs/LLM_COORDINATION_NOTES.md`](docs/LLM_COORDINATION_NOTES.md) · v2 rebuild spec
[`docs/SPEC_rebuild_v2.md`](docs/SPEC_rebuild_v2.md) · animated-deck kit [`presentation/`](presentation/) · v1 and
finished docs [`docs/archive/`](docs/archive/).

Layout: `config/` pinned sources and settings (`datasets.yaml`, `builds.yaml`) · `data/raw` untouched downloads,
`data/interim`, `data/processed` (not in git) · `scripts/`, `R/`, `results/`, `logs/`. Docs are kept on `main`; the
`celligner` worktree (`~/projects/repurposing-celligner`) can run drug-side code.

Gotchas: Toil = GENCODE v23, recount3 = v26; join on unversioned Ensembl IDs. GTEx whole blood is not a marrow
reference.
