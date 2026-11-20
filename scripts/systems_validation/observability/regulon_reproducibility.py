#!/usr/bin/env python
"""D149: separate ESTIMATION NOISE from CONTEXT-SPECIFIC REWIRING in ARACNe regulons.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python  (R part: build_aracne_pkn-style read)

D138.1 reported that the same target's footprint shares ~1% of its genes between ARACNe networks,
and I read that as ARACNe failing to reproduce itself. The user challenged that reading on 2026-11-19:
context-dependent targets are expected biology, because a regulator's accessible sites and co-factors
differ by cell type, so the same drug SHOULD engage different downstream genes in different diseases.

That is testable. `regulons_aracne_project_v2_tf_brca_controls.rds` holds BRCA_SUB1/2/3 - three
networks built from 769 RANDOM BRCA tumours each (D81). Same tissue, same disease mix, no biological
difference. Any regulon disagreement between them is estimation noise with the biology held fixed.

Output: results/systems_validation/observability/regulon_reproducibility.tsv
"""
import itertools, json, subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
OUT = ROOT / "results/systems_validation/observability"
RS = "/home/kirti/miniforge3/envs/repurposing/bin/Rscript"
GENES = ["ESR1", "MYC", "FOXA1", "GATA3", "SOX10", "JUN", "ETS1", "E2F1"]
BREAST = ["BRCA", "BRCA_LUM", "BRCA_LUMA", "BRCA_LUMB", "BRCA_BAS"]


def jac(a, b):
    u = a | b
    return len(a & b) / len(u) if u else np.nan


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # control networks live only in the .rds; dump the membership we need
    r = subprocess.run([RS, "-e", f'''
suppressPackageStartupMessages(library(data.table))
C <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf_brca_controls.rds")
g <- c({",".join(repr(g) for g in GENES)})
out <- rbindlist(lapply(names(C), function(n) rbindlist(lapply(g, function(x)
  if (x %in% names(C[[n]])) data.table(network=n, regulator=x, gene_id=names(C[[n]][[x]]$tfmode))
  else NULL))))
fwrite(out, "/dev/stdout", sep="\\t")'''], cwd=ROOT, capture_output=True, text=True)
    ctl = pd.read_csv(pd.io.common.StringIO(r.stdout), sep="\t")
    E = pd.read_csv(ROOT / "data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz", sep="\t")
    E = E[E.in_lincs_gene_space]

    rows = []
    for g in GENES:
        prod = {n: set(d.gene_id) for n, d in E[E.regulator == g].groupby("network")}
        c = {n: set(d.gene_id) for n, d in ctl[ctl.regulator == g].groupby("network")}
        def med(sets):
            v = [jac(a, b) for a, b in itertools.combinations(sets, 2)]
            return (round(float(np.median(v)), 4), len(sets)) if len(v) else (np.nan, len(sets))
        noise, n_ctl = med(list(c.values()))
        br, n_br = med([prod[n] for n in BREAST if n in prod])
        allj, n_all = med(list(prod.values()))
        rows.append(dict(regulator=g, noise_floor_same_cohort=noise, n_control_nets=n_ctl,
                         breast_nested=br, n_breast_nets=n_br,
                         all_networks=allj, n_networks=n_all,
                         context_drop_breast=None if (noise != noise or br != br) else round(noise - br, 4),
                         context_drop_all=None if (noise != noise or allj != allj) else round(noise - allj, 4)))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "regulon_reproducibility.tsv", sep="\t", index=False)
    print(df.to_string(index=False))
    ok = df.dropna(subset=["noise_floor_same_cohort"])
    print(f"\nnoise floor (same cohort, biology fixed): median {ok.noise_floor_same_cohort.median():.3f}")
    print(f"across all tissues:                        median {ok.all_networks.median():.3f}")
    (OUT / "regulon_reproducibility.json").write_text(json.dumps(dict(
        step="D149 estimation noise vs context-specific rewiring",
        control="BRCA_SUB1/2/3, 769 random BRCA tumours each (D81)",
        noise_floor_median=round(float(ok.noise_floor_same_cohort.median()), 4),
        all_networks_median=round(float(ok.all_networks.median()), 4),
        rows=rows), indent=1, default=str))


if __name__ == "__main__":
    sys.exit(main())
