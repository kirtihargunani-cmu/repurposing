#!/usr/bin/env python
"""F-1: frozen TF-activity inputs for CARNIVAL, one file per cohort.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
PLAN: docs/PLAN_systems_validation_implementation.md, Phase 1 step F-1.

Exports the disease-side TF activity per cohort in the form CARNIVAL's `measurements` argument
takes, intersected with the frozen PKN's nodes BEFORE anything is run.

Why the intersection happens here: CARNIVAL 2.20.0 aborts the whole call if ANY measured node is
absent from the prior-knowledge network, and reports it as "None of your measurements are in prior
knowledge networks (PKN)" even when most of them are (AUDIT addendum 3, C3 - verified that the
matching node alone runs fine). Intersecting up front, and recording how many TFs were dropped,
turns a misleading crash into a logged coverage number.

NO SELECTION RULE IS APPLIED HERE. Every tested regulator is exported with its continuous NES, its
p/FDR and its regulon size. Which TFs become measurements, and how a continuous NES is turned into
the sign CARNIVAL's ILP constrains on, is a FROZEN CHOICE made in F-2 before any biological run -
the spec is explicit that a continuous VIPER NES must never become a binary sign without a
reproducible rule.

Outputs (results/systems_validation/fault_isolation/):
  tf_activity_<cohort>.tsv   regulator, nes, p, fdr, size, in_pkn, has_outgoing_edge
  f1_tf_activity_manifest.json   per-cohort counts, dropped TFs, sha256 of every file
"""
import hashlib, json, os, subprocess, sys
from pathlib import Path
import pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
OUT = ROOT / "results/systems_validation/fault_isolation"
PKN = ROOT / "results/systems_validation/pkn/signed_pkn.tsv"
MR = ROOT / "results/step6/aracne_project_v2_tf/group_mr.tsv"
MIN_TFS = 20          # plan F-1 exit criterion
SPEC_PILOTS = {"TCGA-SKCM": "BRAF-V600 melanoma", "BRCA-Her2": "ERBB2-amplified breast",
               "BRCA-LumA": "ESR1 luminal breast (SEEN development case)",
               "TCGA-LUAD": "EGFR-driver lung"}
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "4"


def sha256(p, chunk=1 << 20):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pkn = pd.read_csv(PKN, sep="\t")
    nodes, sources = set(pkn.source) | set(pkn.target), set(pkn.source)
    mr = pd.read_csv(MR, sep="\t")
    print(f"PKN {len(nodes):,} nodes ({len(sources):,} with outgoing edges); "
          f"group_mr {len(mr):,} rows over {mr.group_id.nunique()} cohorts\n")

    rows, files = [], {}
    for g, d in mr.groupby("group_id"):
        d = d.copy()
        d["in_pkn"] = d.regulator.isin(nodes)
        d["has_outgoing_edge"] = d.regulator.isin(sources)
        kept = d[d.in_pkn]
        f = OUT / f"tf_activity_{g}.tsv"
        cols = [c for c in ["regulator", "nes", "p", "fdr", "size", "in_pkn",
                            "has_outgoing_edge"] if c in d.columns]
        d.sort_values("nes", key=abs, ascending=False)[cols].to_csv(f, sep="\t", index=False)
        files[f] = None
        rows.append(dict(
            group_id=g, tested_regulators=len(d), in_pkn=int(d.in_pkn.sum()),
            dropped_not_in_pkn=int((~d.in_pkn).sum()),
            pct_in_pkn=round(100 * d.in_pkn.mean(), 1),
            with_outgoing_edge=int(d.has_outgoing_edge.sum()),
            sig_fdr05_in_pkn=int((kept.fdr < 0.05).sum()) if "fdr" in kept else None,
            max_abs_nes=round(float(d.nes.abs().max()), 3),
            meets_min_tfs=bool(d.in_pkn.sum() >= MIN_TFS),
            spec_pilot=SPEC_PILOTS.get(g, "")))
    cov = pd.DataFrame(rows).sort_values("in_pkn", ascending=False)
    cov.to_csv(OUT / "f1_tf_coverage.tsv", sep="\t", index=False)

    man = dict(step="F-1 frozen TF-activity inputs",
               plan="docs/PLAN_systems_validation_implementation.md",
               git_sha=subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                      capture_output=True, text=True).stdout.strip(),
               pkn=str(PKN.relative_to(ROOT)), pkn_sha256=sha256(PKN),
               source=str(MR.relative_to(ROOT)), source_sha256=sha256(MR),
               min_tfs_exit_criterion=MIN_TFS,
               cohorts=int(len(cov)),
               cohorts_meeting_exit=int(cov.meets_min_tfs.sum()),
               no_selection_rule_applied=True,
               nes_to_sign_rule="NOT SET HERE - frozen in F-2",
               files={str(k.relative_to(ROOT)): sha256(k) for k in files},
               coverage=cov.to_dict("records"))
    (OUT / "f1_tf_activity_manifest.json").write_text(json.dumps(man, indent=1, default=str))

    print("SPEC-named pilot cohorts:")
    print(cov[cov.spec_pilot != ""][["group_id", "tested_regulators", "in_pkn", "pct_in_pkn",
                                     "dropped_not_in_pkn", "sig_fdr05_in_pkn", "spec_pilot"]]
          .to_string(index=False))
    print(f"\nall {len(cov)} cohorts: TFs in PKN median {cov.in_pkn.median():.0f} "
          f"(min {cov.in_pkn.min()}, max {cov.in_pkn.max()}); "
          f"median {cov.pct_in_pkn.median():.1f}% of tested regulators retained")
    print(f"cohorts meeting the >= {MIN_TFS} exit criterion: {int(cov.meets_min_tfs.sum())} of {len(cov)}")
    print("\nwrote", OUT)


if __name__ == "__main__":
    sys.exit(main())
