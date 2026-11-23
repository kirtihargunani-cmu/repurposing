# repurposing

Drug repurposing from regulator activity. For each cancer (TCGA / TARGET, with GTEx as the healthy reference) we find
the regulators the tumor turns on or off and the genes behind them (**disease side**), then rank FDA drugs whose effect
in cancer-matched cell lines (LINCS L1000) pushes those genes back toward healthy, provided the tumor actually leans on
the drug's target (**drug side**). Every decision is numbered (D1-D212) and recorded; tests are pre-registered before
they are run.

This README says what the `main` branch contains, what is deliberately not in git, and where to start reading.
Updated 2026-11-23.

## Where things stand (2026-11-23)

| | |
|---|---|
| **Production ranking** (D211) | A drug is eligible for a cancer when its target evidence `T` (CRISPR dependency, alteration, overactivity of its target in that cancer) is >= 50; eligible drugs are ranked by **gene-level reversal** (`rev_gene`), separately at **early (<= 6 h)** and **late (>= 24 h)** timepoints. Chemotherapy / DNA-acting drugs are ranked separately. |
| **Regulator networks** | ARACNe networks built per cancer from its own tumors, read with VIPER. ARACNe only (D160, D208); the curated sets are archived. |
| **Primary validation** | D127 biomarker pairs (does a targeted drug rank higher in the cancer carrying its biomarker?), counted one vote per target (D209). The previous lead score passes, 16/20 (p = 0.006). The new gene-level ranking does not (early 2/8, late 1/8): adopted by decision and **under review**. |
| **Reproducibility** | `scripts/run_pipeline.sh` rebuilds every result after the network builds; two full runs from the same start match on content (D210). |
| **Results to look at** | FDA rankings page: `results/lincs/gene_rankings_page/` (see [docs/artifacts](docs/artifacts/README.md) for how to open it), plus slide decks and pages in [`docs/artifacts/`](docs/artifacts/). |

## Start here

| To understand | Read |
|---|---|
| The method, step by step (technical) | [`docs/PIPELINE_V2.md`](docs/PIPELINE_V2.md) (disease side) · [`docs/PIPELINE_DRUG_SIDE_V2.md`](docs/PIPELINE_DRUG_SIDE_V2.md) (drug side) |
| The same, in plain language with results | [`docs/DISEASE_SIDE_PLAIN_LANGUAGE.md`](docs/DISEASE_SIDE_PLAIN_LANGUAGE.md) · [`docs/DRUG_SIDE_PLAIN_LANGUAGE.md`](docs/DRUG_SIDE_PLAIN_LANGUAGE.md) |
| Every decision and why (D1-D212) | [`docs/DECISIONS.md`](docs/DECISIONS.md); drug-side short list and caveats C1-C33: [`docs/notes/drug_side_decisions_and_caveats.md`](docs/notes/drug_side_decisions_and_caveats.md) |
| The exact formulas | [`docs/FORMULAS.md`](docs/FORMULAS.md) |
| What runs, in what order | [`docs/PIPELINE_RUN_ORDER.md`](docs/PIPELINE_RUN_ORDER.md) and [`scripts/run_pipeline.sh`](scripts/run_pipeline.sh) |
| Run history with all numbers | [`docs/PROGRESS.md`](docs/PROGRESS.md) |
| Pre-registered tests | `docs/SPEC_*.md` (each committed before its run; the commit is the registration) |
| Presentations | [`docs/artifacts/`](docs/artifacts/README.md) |

## What `main` contains

| Folder | What is in it |
|---|---|
| `scripts/` | All pipeline code (353 files). Numbered by stage: `10-21` disease groups and covariates, `30` signatures, `40-46` regulons / ARACNe, `50-61` VIPER and patient corrections, `70-79*` drug side (Celligner, LINCS, drug VIPER, scoring, flags, ranking, pages), `80-96` clustering and subgroups. `dNNN/` folders hold the scripts of individual decisions and tests. `run_pipeline.sh` is the production runner; `runlog.sh` wraps every build with a provenance record; `check_builds.py` / `check_params.py` are the pre-run gates. |
| `config/` | Pinned data sources (`datasets.yaml`, with versions and URLs), build registry (`builds.yaml`, `build_outputs.tsv`), adopted thresholds and which decision set them (`adopted_parameters.yaml`), drug gate rules, network and tissue maps. |
| `results/` | Every table, summary and page the pipeline writes that is small enough for git (1,422 files): `step1-step10_v2` disease side, `celligner/`, `lincs/` drug side (rankings, flags, labels, tests, D127), `systems_validation/`. Each build folder carries a `.build_manifest.json` saying which run wrote it. |
| `docs/` | Pipeline, decisions, formulas, specs, progress, plain-language docs, `artifacts/` (presentation copies), `archive/` (v1 and finished docs). |
| `R/`, `tests/` | Shared R paths; tests for the D128 module. |
| `archive/collectri/results/` | The retired CollecTRI / DoRothEA results, kept for the record (D208). |
| `presentation/` | Animated-deck kit. |
| `environment*.yml` | Conda environments with exact versions (`environment.lock.yml` disease side, R; `environment.celligner.lock.yml` drug side, Python). |

Branches: **`main` is the source of truth** for code, docs and results. `celligner` is a worktree branch kept
identical to `main` (merged after every change). `feature/77b-hierarchical-transfer-20261113` holds 4 unmerged
commits from an earlier cell-line transfer benchmark.

## What is not in git, and why

Large inputs and intermediate files stay on the lab server (`MAIL-10`, `~/projects/repurposing`); git keeps the
checksums, versions and code needed to rebuild them.

| Not in git | Size | What it is | How to get it back |
|---|---|---|---|
| `data/raw/` | 62 GB | Downloads: Toil TCGA-TARGET-GTEx, TCGA PanCanAtlas, LINCS L1000 (2020), DepMap 24Q4, ChEMBL 37, Celligner, openFDA labels | Sources, versions and URLs in `config/datasets.yaml`; checksums of 85 files in `data/raw/checksums.sha256` (tracked) |
| `data/processed/` | 146 GB | ARACNe networks and regulons, VIPER activity (disease patients, drug x cell line, drug x condition), LINCS profile matrices | Rebuilt by the numbered scripts; ARACNe took ~2.5 days on 100 bootstraps; regulon checksums in `results/step5/regulon_checksums.sha256` |
| `data/interim/`, `data/backup_*`, `data/pre_chain_snapshot_*` | ~1.3 GB | Caches (API pulls, label text) and pre-rebuild safety copies | Caches refill on rerun; backups are not needed |
| `*.gz` under `results/` | 2.9 GB, 91 files | Gene-level contribution tables (the input to `rev_gene`), timepoint reversal tables, pair labels, disease-side contributions | Rebuilt by `scripts/run_pipeline.sh` from `data/processed` |
| `logs/` | 36 MB | Run logs and the build ledger (`logs/BUILD_LOG.tsv`) | Written by every run |
| `archive/collectri/data/` | 117 MB | Retired CollecTRI / DoRothEA processed data (D208) | Not needed for any current result |

So a clone of `main` can read every final ranking, decision and test result, and all of the code; re-running the
pipeline needs the server's `data/` folder.

## Running it

```
# everything after the ARACNe / VIPER builds, in order, with the threshold gate first
scripts/run_pipeline.sh                    # FROM=core|windows|labels|tests to resume; PAGES=1 to rebuild pages
# the staleness gate currently fails on 15 old bookkeeping entries (D210); until they are cleared:
SKIP_GATE=1 scripts/run_pipeline.sh
```

Upstream of the runner (rarely rerun): disease side `scripts/run_v2.sh` (Steps 1-10), ARACNe
(`scripts/43_run_aracne.sh`), disease VIPER (`scripts/50_viper.R`), drug VIPER (`scripts/79d_drug_viper.R`). Use the
interpreters in `environment*.yml`; scripts state theirs in their headers. Thread budget: 6 cores by default
(`N_CORES`).

Gotchas: Toil = GENCODE v23, recount3 = v26, so join on unversioned Ensembl IDs. GTEx whole blood is not a marrow
reference. A `.gz` written by Python carries a timestamp, so compare decompressed content, not file hashes (D210).
