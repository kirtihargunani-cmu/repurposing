#!/usr/bin/env python
"""Rebuild the D198 programs (79ai) and D199 flags (79aj) for Experiment E, without touching production.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
Amendment 1e. Runs the production scripts' own code, imported unchanged from the input checkout, with their path
constants redirected: every input is read from `input_root` (the 2026-11-23 run, read-only) and every output goes to
results/lincs/mechanistic_incremental_v1/E_program_access/programs/. Records the module file's checksum.
"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import yaml

WT = Path(__file__).resolve().parents[2]
CFG = yaml.safe_load((WT / "config/mechanistic_incremental_v1.yaml").read_text())
IN = Path(CFG["input_root"])
OUT = WT / CFG["output_root"] / "E_program_access" / "programs"
OUT.mkdir(parents=True, exist_ok=True)
MODS = IN / "results/step6/aracne_project_v2_tf/diagnostics/regulator_modules_v2_tf.tsv"
ACT = IN / "results/lincs/regulon_contributions_v2_tf/regulator_contributions_by_window_all.tsv.gz"


def load(name):
    spec = importlib.util.spec_from_file_location(name, IN / "scripts" / f"{name}.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


ai = load("79ai_programs")
ai.ROOT, ai.RES, ai.OUT, ai.MODS, ai.ACT = IN, IN / "results/lincs", OUT, MODS, ACT
ai.EDGES = IN / "data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz"
ai.GMT = IN / "data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt"
rc = ai.main()
aj = load("79aj_program_filter")
aj.ROOT, aj.RES, aj.OUT, aj.MODS, aj.ACT = IN, IN / "results/lincs", OUT, MODS, ACT
rc2 = aj.main()
(OUT / "inputs.json").write_text(json.dumps({
    "modules": {"path": str(MODS), "sha256": hashlib.sha256(MODS.read_bytes()).hexdigest()},
    "activity": str(ACT), "scripts": ["79ai_programs.py", "79aj_program_filter.py"], "rc": [rc, rc2]}, indent=1))
sys.exit(max(rc or 0, rc2 or 0))
