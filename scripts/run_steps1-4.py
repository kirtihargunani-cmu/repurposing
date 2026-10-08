#!/usr/bin/env python
"""Phase 1 runner: Steps 1-4 (downloads -> disease grounding -> subtype splits -> covariates/ESTIMATE ->
limma-voom signatures), by running the existing step scripts in order.

The analysis code stays in the step scripts (Python, R, bash); this file only orders them, checks inputs,
skips steps whose outputs already exist, logs each step and stops at the first failure, so results are
byte-for-byte what the step scripts produce on their own.

Steps (run in this order):
  setup_env              00_setup_env.sh        create conda env "repurposing" (opt-in: --setup-env)
  download               01_download.sh         Toil, BeatAML, MONDO, UBERON, GENCODE, PanCanAtlas, OncoTree
  download_covariates    02_download_covariates.sh  GTEx annotations, TCGA plates/ABSOLUTE/LAML, TARGET clinical
  ground                 10_ground_toil.py      Step 1: disease grounding
  split                  11_split_subtypes.py   Step 1b: subtype splits, sample -> group
  estimate               21_estimate.R          Step 3: ESTIMATE stromal/immune scores (48 cores)
  covariates             20_covariates.py       Step 3: covariate tables (merges ESTIMATE)
  signatures             30_signatures.R        Step 4: protein-coding signatures (drug matching)
  signatures_extended    30_signatures.R        Step 4: + lincRNA/antisense (GENESET=extended; mechanism only)

Downloads: some sources are fetched as "latest" (OncoTree, GDC aliquot list), so re-downloading can change
pinned inputs. Download steps therefore never re-run when their files exist, even with --force, unless
--redownload is given.

Usage:
  scripts/run_steps1-4.py                     run every step whose outputs are missing
  scripts/run_steps1-4.py --list              show steps and whether each is done
  scripts/run_steps1-4.py --from covariates   (re)run from a step onward (implies --force for those steps)
  scripts/run_steps1-4.py --only signatures --force
  scripts/run_steps1-4.py --dry-run
Logs: logs/pipeline/steps1-4_<step>.log per step; logs/pipeline/steps1-4_runs.tsv (one row per step run).
"""
import argparse
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
PHASE = "steps1-4"
LOGDIR = ROOT / "logs/pipeline"


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


STEPS = [
    Step("setup_env", ["bash", "scripts/00_setup_env.sh"], ["environment.lock.yml"],
         desc="create conda env 'repurposing' + ESTIMATE", optin="setup_env"),
    Step("download", ["bash", "scripts/01_download.sh"],
         ["data/raw/solid_toil/TcgaTargetGtex_gene_expected_count.gz", "data/raw/solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz",
          "data/raw/solid_toil/TcgaTargetGTEX_phenotype.txt.gz", "data/raw/mondo/mondo.obo", "data/raw/uberon/uberon-basic.obo",
          "data/raw/gencode/gencode.v23.annotation.gtf.gz", "data/raw/tcga_pancan/TCGASubtype.20170308.tsv.gz",
          "data/raw/oncotree/tumorTypes.json", "data/raw/checksums.sha256"],
         desc="raw expression + ontologies", download=True),
    Step("download_covariates", ["bash", "scripts/02_download_covariates.sh"],
         ["data/raw/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
          "data/raw/tcga_pancan/TCGA_mastercalls.abs_tables_JSedit.fixed.txt", "data/raw/tcga_pancan/LAML_clinicalMatrix",
          "data/raw/tcga_gdc/rnaseq_aliquots_*.tsv", "data/raw/target_clinical/TARGET_WT_ClinicalData_Validation_20230322.xlsx"],
         desc="covariate sources", download=True),
    Step("ground", ["python", "scripts/10_ground_toil.py"],
         ["results/step1/disease_grounding_toil.tsv", "results/step1/mondo_descendants.tsv"],
         inputs=["config/toil_disease_curation.tsv", "data/raw/mondo/mondo.obo"], desc="Step 1 disease grounding"),
    Step("split", ["python", "scripts/11_split_subtypes.py"],
         ["results/step1/disease_groups_toil.tsv", "results/step1/sample_groups_toil.tsv"],
         inputs=["config/subtype_splits.tsv", "results/step1/disease_grounding_toil.tsv"], desc="Step 1b subtype splits"),
    Step("estimate", ["Rscript", "scripts/21_estimate.R"], ["results/step3/estimate_toil.tsv"],
         inputs=["results/step1/sample_groups_toil.tsv", "data/raw/solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz"],
         desc="Step 3 ESTIMATE scores"),
    Step("covariates", ["python", "scripts/20_covariates.py"],
         ["results/step3/covariates_toil.tsv", "results/step3/covariates_beataml.tsv",
          "results/step3/genes_protein_coding_v23.tsv", "results/step3/genes_extended_v23.tsv"],
         inputs=["results/step3/estimate_toil.tsv", "results/step1/sample_groups_toil.tsv"], desc="Step 3 covariates"),
    Step("signatures", ["Rscript", "scripts/30_signatures.R"],
         ["results/step4/signature_summary.tsv", "results/step4/z_tumor_vs_gtex.tsv", "results/step4/z_tumor_vs_adjacent.tsv",
          "data/processed/step4/counts_protein_coding.rds"],
         inputs=["results/step3/covariates_toil.tsv", "results/step3/genes_protein_coding_v23.tsv"],
         env={"GENESET": "protein_coding"}, desc="Step 4 signatures (protein-coding)"),
    Step("signatures_extended", ["Rscript", "scripts/30_signatures.R"],
         ["results/step4_extended/signature_summary.tsv", "results/step4_extended/z_tumor_vs_gtex.tsv"],
         inputs=["results/step3/covariates_toil.tsv", "results/step3/genes_extended_v23.tsv"],
         env={"GENESET": "extended"}, desc="Step 4 signatures (extended gene set)"),
]


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


def main():
    ap = base_parser(__doc__)
    ap.add_argument("--setup-env", action="store_true", help="include 00_setup_env.sh (creates the conda env)")
    args = ap.parse_args()
    return main_loop(args, select(args), {"setup_env": args.setup_env})


if __name__ == "__main__":
    sys.exit(main())
