#!/usr/bin/env python
"""Phase 2 runner: Steps 5-7 (regulons -> ARACNe -> VIPER -> purity correction, plus diagnostics), by running
the existing step scripts in order. Needs the Steps 1-4 outputs (scripts/run_steps1-4.py).

The analysis code stays in the step scripts (Python, R, Java via bash); this file only orders them, checks
inputs, skips steps whose outputs already exist, logs each step and stops at the first failure, so results are
byte-for-byte what the step scripts produce on their own.

Steps (run in this order; [x] = switch):
  download_regulon_sources  03_download_regulons.sh  Lambert TFs, HGNC, MC3 mutations, GISTIC copy number
  id_map                    40_id_map.py             GENCODE v23 Ensembl -> HGNC map, TF universe
  regulons                  41_regulons.R            CollecTRI + DoRothEA A-C (OmniPath), aracne.networks
  aracne_inputs             42_aracne_inputs.py      per-cohort expression + regulator lists for ARACNe
  aracne_build              43_run_aracne.sh         ARACNe-AP, 100 bootstraps x 28 cohorts, ~2.5 days,
                                                     resumable  [--build-aracne]
  aracne_regulons           44_aracne_regulons.R     networks -> VIPER regulons
  edge_support              45_edge_support.R        bootstrap support per edge          [edge diagnostics]
  edge_confounding          46_edge_confounding.R    edge strength after purity/plate/subtype adjustment
                                                                                         [edge diagnostics]
  viper_<regulon>           50_viper.R               group + patient VIPER, one step per --regulons entry
  compare_<regulon>         51_compare_mr.R          cross-group summary + known-biology controls
  download_replication      04_download_replication.sh  METABRIC, GSE39582, CPTAC-BRCA/COAD (Step 9)
  purity_aracne_project     56_purity_correct.R      ESTIMATE regression of ARACNe patient scores
  quality_correct           61_quality_correct.R     gene-level purity + RNA-quality correction, patient VIPER
  driver_classes            60_driver_classes.py     kinase vs TF driver alterations     [exploratory]
  subtype_test              52_subtype_test.R        subtype-split experiment (CollecTRI) [exploratory]
  aracne_precomputed_test   53_aracne_precomputed_test.R  aracne.networks VIPER test     [exploratory]
  patient_diagnostics       54_patient_nes_diagnostics.R  plate / ESTIMATE / redundancy   [exploratory]
  brca_heterogeneity        55_brca_heterogeneity.R  BRCA patient heterogeneity         [exploratory]

Regulons (--regulons, default all three): collectri, dorothea_ABC, aracne_project. The exploratory steps need
collectri (and 54 also aracne_project) results.

Downloads: HGNC is fetched as the current release, and CollecTRI/DoRothEA are pulled live from OmniPath unless a
saved pull exists in data/raw/regulons/ (41_regulons.R reuses the newest one). The download step therefore never
re-runs when its files exist, unless --redownload is given.

Usage:
  scripts/run_steps5-7.py                      run every enabled step whose outputs are missing
  scripts/run_steps5-7.py --list
  scripts/run_steps5-7.py --build-aracne       include the ARACNe-AP build (run under nohup/tmux: days)
  scripts/run_steps5-7.py --regulons aracne_project --from viper_aracne_project
  scripts/run_steps5-7.py --no-exploratory --no-edge-diagnostics
Logs: logs/pipeline/steps5-7_<step>.log per step; logs/pipeline/steps5-7_runs.tsv (one row per step run).
"""
import argparse
import csv
import datetime as dt
import glob
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PHASE = "steps5-7"
LOGDIR = ROOT / "logs/pipeline"
REGULONS = ["collectri", "dorothea_ABC", "aracne_project"]
NETDIR = "data/processed/step5/aracne"


@dataclass
class Step:
    name: str
    cmd: list                      # [interpreter, script, args...]; interpreter: "python" | "Rscript" | "bash"
    outputs: list                  # paths/globs relative to ROOT (or callables returning what is missing);
    inputs: list = field(default_factory=list)  # step is done when nothing is missing
    env: dict = field(default_factory=dict)
    desc: str = ""
    download: bool = False         # never re-run when outputs exist unless --redownload
    optin: str = ""                # switch that must be on for the step to run
    gate_hint: str = ""            # how to turn that switch on/off (shown when skipped)


def networks():
    """ARACNe cohorts from 42_aracne_inputs.py (empty until that step has run)."""
    f = ROOT / "results/step5/aracne_cohorts.tsv"
    if not f.exists():
        return []
    with open(f) as fh:
        return [r["net"] for r in csv.DictReader(fh, delimiter="\t")]


def per_network(fname):
    def check():
        nets = networks()
        if not nets:
            return ["results/step5/aracne_cohorts.tsv"]
        return [p for p in (f"{NETDIR}/{n}/{fname}" for n in nets) if not (ROOT / p).exists()]
    return check


def regulons_fresh():
    """aracne_project regulons must exist and be newer than every finished network."""
    rds = ROOT / "data/processed/step5/regulons_aracne_project.rds"
    if not rds.exists():
        return [str(rds.relative_to(ROOT))]
    nets = [ROOT / NETDIR / n / "network.txt" for n in networks()]
    stale = [str(n.relative_to(ROOT)) for n in nets if n.exists() and n.stat().st_mtime > rds.stat().st_mtime]
    return [f"{rds.relative_to(ROOT)} (older than {len(stale)} network(s))"] if stale else []


def regulon_source(reg):
    return {"aracne_project": "data/processed/step5/regulons_aracne_project.rds"}.get(reg, "results/step5/regulons_core.rds")


def make_steps(args):
    edge = dict(optin="edge_diagnostics", gate_hint="on by default; off with --no-edge-diagnostics")
    expl = dict(optin="exploratory", gate_hint="on by default; off with --no-exploratory")
    sig = ["results/step4/z_tumor_vs_gtex.tsv", "data/processed/step4/counts_protein_coding.rds", "results/step3/covariates_toil.tsv"]
    steps = [
        Step("download_regulon_sources", ["bash", "scripts/03_download_regulons.sh"],
             ["data/raw/tf_lambert/DatabaseExtract_v_1.01.csv", "data/raw/hgnc/hgnc_complete_set_*.txt",
              "data/raw/tcga_pancan/mc3.v0.2.8.PUBLIC.maf.gz",
              "data/raw/tcga_pancan/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"],
             desc="TF list, HGNC, mutations, copy number", download=True),
        Step("download_replication", ["bash", "scripts/04_download_replication.sh"],
             ["data/raw/replication/metabric/data_mrna_illumina_microarray.txt",
              "data/raw/replication/gse39582/GSE39582_series_matrix.txt.gz",
              "data/raw/replication/cptac_brca/gdc_manifest.tsv", "data/raw/replication/cptac_coad/gdc_manifest.tsv",
              "data/raw/replication/cptac_coad/cbioportal/data_clinical_sample.txt"],
             desc="Step 9 replication cohorts (METABRIC, GSE39582, CPTAC-BRCA/COAD)", download=True),
        Step("id_map", ["python", "scripts/40_id_map.py"], ["results/step5/id_map.tsv", "results/step5/tf_universe.tsv"],
             inputs=["data/raw/hgnc/hgnc_complete_set_*.txt", "data/raw/tf_lambert/DatabaseExtract_v_1.01.csv"],
             desc="Step 5a ID map"),
        Step("regulons", ["Rscript", "scripts/41_regulons.R"],
             ["results/step5/regulons_core.rds", "data/processed/step5/regulons_aracne.rds", "results/step5/regulon_coverage.tsv"],
             inputs=["results/step5/id_map.tsv", "results/step5/tf_universe.tsv"], desc="Step 5b CollecTRI, DoRothEA, aracne.networks"),
        Step("aracne_inputs", ["python", "scripts/42_aracne_inputs.py"],
             ["results/step5/aracne_cohorts.tsv", "results/step5/aracne_regulators.tsv"],
             inputs=["results/step4/signature_summary.tsv", "data/processed/step5/regulons_aracne_*.tsv"],
             desc="Step 5c ARACNe inputs"),
        Step("aracne_build", ["bash", "scripts/43_run_aracne.sh"], [per_network("network.txt")],
             inputs=["results/step5/aracne_cohorts.tsv"], desc="Step 5d ARACNe-AP networks (~2.5 days)",
             optin="build_aracne", gate_hint="--build-aracne"),
        Step("aracne_regulons", ["Rscript", "scripts/44_aracne_regulons.R"], [regulons_fresh],
             inputs=[per_network("network.txt")], desc="Step 5e networks -> regulons"),
        Step("edge_support", ["Rscript", "scripts/45_edge_support.R"],
             [per_network("edge_support.tsv.gz"), "results/step5/aracne_edge_support.tsv"],
             inputs=[per_network("network.txt")], desc="Step 5f edge bootstrap support", **edge),
        Step("edge_confounding", ["Rscript", "scripts/46_edge_confounding.R"],
             [per_network("edge_confounding.tsv.gz"), "results/step5/aracne_edge_confounding.tsv"],
             inputs=[per_network("edge_support.tsv.gz")], desc="edge confounding diagnostic", **edge),
    ]
    for reg in args.regulons:
        steps.append(Step(f"viper_{reg}", ["Rscript", "scripts/50_viper.R"],
                          [f"results/step6/{reg}/group_summary.tsv", f"results/step6/{reg}/group_mr.tsv"],
                          inputs=[regulon_source(reg)] + sig, env={"REGULON": reg, "N_CORES": str(args.cores)},
                          desc=f"Step 6 VIPER ({reg})"))
        steps.append(Step(f"compare_{reg}", ["Rscript", "scripts/51_compare_mr.R"],
                          [f"results/step6/{reg}/fdr10_group_summary.tsv", f"results/step6/{reg}/controls.tsv"],
                          inputs=[f"results/step6/{reg}/group_mr.tsv"], env={"REGULON": reg},
                          desc=f"Step 6b cross-group summary ({reg})"))
    if "aracne_project" in args.regulons:
        steps.append(Step("purity_aracne_project", ["Rscript", "scripts/56_purity_correct.R"],
                          ["results/step6/aracne_project/purity_correction/summary.tsv"],
                          inputs=["results/step6/aracne_project/group_summary.tsv"], env={"REGULON": "aracne_project"},
                          desc="Step 6c purity correction (ARACNe patient scores)"))
        steps.append(Step("quality_correct", ["Rscript", "scripts/61_quality_correct.R"],
                          ["results/step6/aracne_project/quality_correction/checks.tsv"],
                          inputs=["results/step6/aracne_project/regulator_filter.tsv",
                                  "results/step6/aracne_project/purity_correction/summary.tsv"],
                          env={"N_CORES": str(args.cores)},
                          desc="Step 6e gene-level purity + RNA-quality correction -> patient VIPER (Step 8 input)"))
    steps += [
        Step("driver_classes", ["python", "scripts/60_driver_classes.py"], ["results/step7/driver_classes.tsv"],
             inputs=["results/step6/collectri/fdr10_group_summary.tsv", "data/raw/tcga_pancan/mc3.v0.2.8.PUBLIC.maf.gz"],
             desc="Step 7 driver classes", **expl),
        Step("subtype_test", ["Rscript", "scripts/52_subtype_test.R"], ["results/step6/collectri/subtype_test.tsv"],
             inputs=["results/step6/collectri/regulator_filter.tsv"], desc="subtype-split experiment", **expl),
        Step("aracne_precomputed_test", ["Rscript", "scripts/53_aracne_precomputed_test.R"],
             ["results/step6/aracne_precomputed/summary.tsv"],
             inputs=["data/processed/step5/regulons_aracne.rds", "results/step6/collectri/group_mr.tsv"],
             desc="aracne.networks VIPER test", **expl),
        Step("patient_diagnostics", ["Rscript", "scripts/54_patient_nes_diagnostics.R"],
             ["results/step6/aracne_project/diagnostics/pc_plate.tsv", "results/step6/aracne_project/diagnostics/redundancy.tsv"],
             inputs=["data/processed/step6/aracne_project/BRCA-LumA.patient_nes.tsv.gz",
                     "data/processed/step6/collectri/BRCA-LumA.patient_nes.tsv.gz"],
             desc="patient-score plate / ESTIMATE / redundancy diagnostics", **expl),
        Step("brca_heterogeneity", ["Rscript", "scripts/55_brca_heterogeneity.R"],
             ["results/step6/aracne_project/diagnostics/brca_heterogeneity_pca.tsv"],
             inputs=["data/processed/step6/aracne_project/BRCA-LumA.patient_nes.tsv.gz"],
             desc="BRCA patient heterogeneity", **expl),
    ]
    return steps


# ---------------------------------------------------------------- shared runner code (same in run_steps5-7.py)
def exists(pattern):
    return bool(glob.glob(str(ROOT / pattern)))


def missing(paths):
    out = []
    for p in paths:
        out += p() if callable(p) else ([] if exists(p) else [p])
    return out


def interpreter(kind, env_prefix):
    if kind == "bash":
        return shutil.which("bash")
    exe = Path(env_prefix) / "bin" / kind
    if exe.exists():
        return str(exe)
    found = shutil.which(kind)
    if not found:
        sys.exit(f"cannot find {kind}: not in {env_prefix}/bin and not on PATH (set --env-prefix)")
    return found


def run_step(step, args, log_tsv):
    log = LOGDIR / f"{PHASE}_{step.name}.log"
    cmd = [interpreter(step.cmd[0], args.env_prefix)] + step.cmd[1:]
    env = {**os.environ, **step.env}
    shown = " ".join(f"{k}={v}" for k, v in step.env.items()) + (" " if step.env else "") + " ".join(step.cmd)
    print(f"[{dt.datetime.now():%H:%M:%S}] RUN  {step.name:24s} {shown}  (log: {log.relative_to(ROOT)})", flush=True)
    if args.dry_run:
        return True
    t0 = time.time()
    with open(log, "w") as fh:
        fh.write(f"# {dt.datetime.now().isoformat(timespec='seconds')}  {shown}\n")
        fh.flush()
        p = subprocess.Popen(cmd, cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in p.stdout:
            fh.write(line)
            if args.verbose:
                print("    " + line, end="")
        rc = p.wait()
    secs = time.time() - t0
    lost = missing(step.outputs) if rc == 0 else []
    status = "ok" if rc == 0 and not lost else ("missing_outputs" if rc == 0 else f"exit_{rc}")
    with open(log_tsv, "a") as fh:
        fh.write(f"{dt.datetime.now().isoformat(timespec='seconds')}\t{step.name}\t{status}\t{secs:.0f}\t{shown}\n")
    print(f"[{dt.datetime.now():%H:%M:%S}] {'DONE' if status == 'ok' else 'FAIL'} {step.name:24s} "
          f"{status} after {secs / 60:.1f} min", flush=True)
    if lost:
        print(f"    expected outputs not written: {', '.join(lost)}")
    if status != "ok":
        print("    last lines of the log:")
        print("".join("    " + l for l in log.read_text().splitlines(True)[-15:]))
    return status == "ok"


def select(args):
    names = [s.name for s in STEPS]
    for n in (args.only or []) + [x for x in (args.from_, args.to) if x]:
        if n not in names:
            sys.exit(f"unknown step '{n}'. Steps: {', '.join(names)}")
    lo = names.index(args.from_) if args.from_ else 0
    hi = names.index(args.to) if args.to else len(names) - 1
    chosen = [s for i, s in enumerate(STEPS) if lo <= i <= hi and (not args.only or s.name in args.only)]
    return chosen


def main_loop(args, chosen, optin_flags):
    LOGDIR.mkdir(parents=True, exist_ok=True)
    log_tsv = LOGDIR / f"{PHASE}_runs.tsv"
    if not log_tsv.exists():
        log_tsv.write_text("time\tstep\tstatus\tseconds\tcommand\n")
    force_all = args.force or bool(args.from_)
    if args.list:
        for s in STEPS:
            state = "done" if not missing(s.outputs) else "todo"
            gate = f" [{s.gate_hint or '--' + s.optin.replace('_', '-')}]" if s.optin else ""
            print(f"  {s.name:24s} {state:5s} {s.desc}{gate}")
        return 0
    for s in chosen:
        explicitly = bool(args.only) and s.name in args.only
        if s.optin and not optin_flags.get(s.optin) and not explicitly:
            print(f"           SKIP {s.name:24s} not enabled ({s.gate_hint or '--' + s.optin.replace('_', '-')})")
            continue
        done = not missing(s.outputs)
        if done and s.download and not args.redownload:
            print(f"           SKIP {s.name:24s} files present (download steps only re-run with --redownload)")
            continue
        if done and not force_all:
            print(f"           SKIP {s.name:24s} outputs exist (use --force or --from to re-run)")
            continue
        need = missing(s.inputs)
        if need and not args.dry_run:
            print(f"           FAIL {s.name:24s} missing inputs: {', '.join(need)}")
            return 1
        if not run_step(s, args, log_tsv):
            print(f"Stopped at '{s.name}'. Fix the error and re-run; finished steps are skipped automatically.")
            return 1
    print("ALL REQUESTED STEPS DONE")
    return 0


def base_parser(doc):
    ap = argparse.ArgumentParser(description=doc.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=doc.split("\n\n", 1)[1])
    ap.add_argument("--list", action="store_true", help="show steps and whether each is done, then exit")
    ap.add_argument("--only", nargs="+", metavar="STEP", help="run only these steps")
    ap.add_argument("--from", dest="from_", metavar="STEP", help="start at this step (re-runs it and later steps)")
    ap.add_argument("--to", metavar="STEP", help="stop after this step")
    ap.add_argument("--force", action="store_true", help="re-run selected steps even if outputs exist")
    ap.add_argument("--redownload", action="store_true", help="allow download steps to fetch again (may change pinned inputs)")
    ap.add_argument("--dry-run", action="store_true", help="print what would run")
    ap.add_argument("--verbose", action="store_true", help="also stream each step's output to the terminal")
    ap.add_argument("--env-prefix", default=os.environ.get("REPURPOSING_ENV", str(Path.home() / "miniforge3/envs/repurposing")),
                    help="conda env holding Rscript and python (default: ~/miniforge3/envs/repurposing)")
    return ap
# ---------------------------------------------------------------- end shared runner code


STEPS = []


def main():
    global STEPS
    ap = base_parser(__doc__)
    ap.add_argument("--regulons", nargs="+", default=REGULONS, choices=REGULONS, help="regulon sets to run VIPER on")
    ap.add_argument("--cores", type=int, default=20, help="cores for VIPER (50_viper.R N_CORES; default 20)")
    ap.add_argument("--build-aracne", action="store_true", help="include the ARACNe-AP build (~2.5 days)")
    ap.add_argument("--no-edge-diagnostics", action="store_true", help="skip edge support / confounding (45, 46)")
    ap.add_argument("--no-exploratory", action="store_true", help="skip exploratory diagnostics (52-55, 60)")
    args = ap.parse_args()
    STEPS = make_steps(args)
    flags = {"build_aracne": args.build_aracne, "edge_diagnostics": not args.no_edge_diagnostics,
             "exploratory": not args.no_exploratory}
    return main_loop(args, select(args), flags)


if __name__ == "__main__":
    sys.exit(main())
