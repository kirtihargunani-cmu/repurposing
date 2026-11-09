#!/usr/bin/env python
"""D117: final tiers and evidence status, recomputed from scratch after every gate and flag (run after 79f3 STAGE=exposure).

Split tables (79i; SFX = _v2_consensus, _v2_consensus_d100early, _v2_consensus_d100late):
  tier_mech (rev_mech_pct), tier_spec_mech (spec_mech_pct), subgroups also tier_spec_sib_mech (spec_sib_mech_pct):
  hit = rev_mech > 0 and the current percentile (after the 79i2 / 79i3 re-ranking) >= 100 - TOP_PCT (79i);
  eligible = fda_rows.eligible (D115) without its rev > 0 term (the hit rule uses rev_mech): passes_gate, not
  unmatched_lines / low_transfer / target_absent / receptor_mismatch / local_route (D107) / D109 weak flags;
  tier 1 = hit and eligible and approved, 2 = hit and eligible, 3 otherwise. Idempotent: never reads the old tiers.
D122: tier_lead (layered_pct >= 90, D122 gates; approved -> 1) and, in the summaries, evidence_status = candidates /
  no adequate compound from tier_lead with layer1_coverage (the reversal-era status kept as evidence_status_reversal).
Summaries (drug_reversal_summary / stratum_reversal_summary): n_t12_matched / n_t12_borrowed / evidence_status /
  evidence_note from the final lead-score tier (groups tier_spec_mech, subgroups tier_spec_sib_mech) with the 79e
  definitions; the 79e values are kept once as *_79e.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
sys.path.insert(0, str(ROOT / "scripts"))
from fda_rows import truthy  # noqa: E402
from route_flag import apply_route  # noqa: E402
from evidence_flags import apply_evidence  # noqa: E402

TOP_PCT = 10            # 79i
LINE_MATCH_MIN = 0.1    # 79e
SFXS = os.environ.get("SFXS", "_v2_consensus _v2_consensus_d100early _v2_consensus_d100late").split()
KEEP = ["n_t12_matched", "n_t12_borrowed", "evidence_status", "evidence_note"]


def note(r):   # 79e wording
    lm = (f"No compound was tested in a line resembling this disease (best line match {r.best_line_match:.4f}, "
          f"{r.best_line if isinstance(r.best_line, str) else 'no line'}; cutoff {LINE_MATCH_MIN})") if r.best_line_match < LINE_MATCH_MIN else \
         f"Compounds tested in matched lines (best {r.best_line}, line match {r.best_line_match:.2f})"
    if r.evidence_status == "matched":
        return f"{lm}; {r.n_t12_matched} tier 1-2 drug(s) rest on matched-line data."
    if r.evidence_status == "borrowed only":
        return (f"{lm}. The {r.n_t12_borrowed} candidate(s) are borrowed from other tissues and rest on the drug acting "
                f"consistently across cell lines (transfer confidence >= 0.5), not on data from this disease's lines.")
    return f"{lm}, and no borrowed drug passes; an adequate compound cannot be reported."


def tiers(df, unit, cols):
    x = apply_evidence(apply_route(df, unit), unit)
    t = lambda c: x[c].map(truthy) if c in x else pd.Series(False, index=x.index)
    rs = x["receptor_status"] if "receptor_status" in x else pd.Series("", index=x.index)
    ok = (x.passes_gate.map(truthy) & ~t("unmatched_lines") & ~t("low_transfer") & ~t("target_absent")
          & ~rs.eq("receptor_mismatch") & ~t("local_route") & ~t("weak_d109"))
    appr = x.approved.map(truthy)
    for tier, pct in cols:
        hit = (x.rev_mech > 0) & (x[pct] >= 100 - TOP_PCT)
        df[tier] = np.where(hit & ok & appr, 1, np.where(hit & ok, 2, 3))
    if "layered_pct" in x:   # D122: lead tier from the layered score and the D122 gates
        ok2 = x.passes_gate.map(truthy) & ~t("target_absent") & ~t("local_route") & ~x.low_exposure.astype(bool)
        hit2 = x.layered_pct >= 100 - TOP_PCT
        df["tier_lead"] = np.where(hit2 & ok2 & appr, 1, np.where(hit2 & ok2, 2, 3))
    return df


def summary(f, df, unit, lead):
    sm = pd.read_csv(f, sep="\t")
    for c in KEEP:
        if c + "_79e" not in sm:
            sm[c + "_79e"] = sm[c]
    t12 = df[lead] <= 2
    nm = df[t12 & df.source.eq("direct")].groupby(unit).size()
    nb = df[t12 & df.source.eq("borrowed")].groupby(unit).size()
    sm["n_t12_matched"] = sm.group_id.map(nm).fillna(0).astype(int)
    sm["n_t12_borrowed"] = sm.group_id.map(nb).fillna(0).astype(int)
    sm["evidence_status"] = np.where(sm.n_t12_matched > 0, "matched", np.where(sm.n_t12_borrowed > 0, "borrowed only", "no adequate compound"))
    sm["evidence_note"] = sm.apply(note, axis=1)
    if "tier_lead" in df:   # D122: status from the layered lead; the reversal-era status kept
        sm["evidence_status_reversal"], sm["evidence_note_reversal"] = sm.evidence_status, sm.evidence_note
        n12 = df[df.tier_lead <= 2].groupby(unit).size()
        cov = df.groupby(unit).apply(lambda d: " ".join(k for k in ("D", "A", "O") if d[f"layered_{k}"].notna().any()) or "none")
        sm["n_t12_lead"] = sm.group_id.map(n12).fillna(0).astype(int)
        sm["layer1_coverage"] = sm.group_id.map(cov).fillna("none")
        sm["evidence_status"] = np.where(sm.n_t12_lead > 0, "candidates", "no adequate compound")
        sm["evidence_note"] = [f"{n} candidate(s) in the top {TOP_PCT}% of the layered score that pass the D122 gates; layer 1 data: {c}"
                               for n, c in zip(sm.n_t12_lead, sm.layer1_coverage)]
    sm.to_csv(f, sep="\t", index=False)
    ch = sm[sm.evidence_status != sm.evidence_status_79e]
    return len(sm), ch[["group_id", "evidence_status_79e", "evidence_status"]]


def main():
    for sfx in SFXS:
        for split, summ, unit, cols, lead in (
                ("drug_reversal_split", "drug_reversal_summary", "group_id",
                 [("tier_mech", "rev_mech_pct"), ("tier_spec_mech", "spec_mech_pct")], "tier_spec_mech"),
                ("stratum_reversal_split", "stratum_reversal_summary", "stratum",
                 [("tier_mech", "rev_mech_pct"), ("tier_spec_mech", "spec_mech_pct"), ("tier_spec_sib_mech", "spec_sib_mech_pct")],
                 "tier_spec_sib_mech")):
            f = RES / f"{split}{sfx}.tsv"
            if not f.exists():
                continue
            df = pd.read_csv(f, sep="\t", low_memory=False)
            old = df[[c for c, _ in cols]].copy()
            old_lead = df["tier_lead"].copy() if "tier_lead" in df else None
            df = tiers(df, unit, cols)
            if "tier_lead" in df:
                cols = cols + [("tier_lead", "layered_pct")]
            df.to_csv(f, sep="\t", index=False)
            if old_lead is not None:
                old["tier_lead"] = old_lead
            msg = "; ".join(f"{c}: t1-2 {int((old[c] <= 2).sum()) if c in old else 0} -> {int((df[c] <= 2).sum())}" for c, _ in cols)
            print(f"{f.name}: {msg}")
            fs = RES / f"{summ}{sfx}.tsv"
            if fs.exists():
                n, ch = summary(fs, df, unit, lead)   # D122: evidence_status from tier_lead when present (inside summary)
                print(f"  {fs.name}: evidence_status changed in {len(ch)} of {n}" + ("".join(f"\n    {r.group_id}: {r.evidence_status_79e} -> {r.evidence_status}" for r in ch.itertuples())))


if __name__ == "__main__":
    main()
