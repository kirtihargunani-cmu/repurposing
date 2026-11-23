#!/usr/bin/env python
"""A self-contained HTML walkthrough: real data in, the math, real data out, per stage.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

User, 2026-11-22: "create an html file that contains snapshots of the data going into each stage
and the math (specifically show the mathematical formulas/transformations)".

Every table on the page is READ FROM THE ACTUAL FILE, not typed in - so the page cannot drift from
what the pipeline produced, and a stage whose output is missing says so instead of showing a
plausible-looking example. Formulas are the ones in docs/FORMULAS.md, which is checked against the
implementations.

Self-contained: no CDN, no external fonts, math rendered as HTML. Offline-openable.

Output: results/lincs/pipeline_explainer.html
"""
import html
import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
CDIR = RES / "regulon_contributions_v2_tf"
OUT = RES / "pipeline_explainer.html"
N = 6


def snap(path, cols=None, n=N, where=None):
    """First n rows of a real file, or an explicit 'missing' marker."""
    p = ROOT / path if not str(path).startswith("/") else Path(path)
    if not p.exists():
        return None, f"not built yet — {path}"
    try:
        kw = dict(sep="\t", low_memory=False)
        if str(p).endswith(".gz") or p.stat().st_size > 50_000_000:
            kw["nrows"] = 400_000 if where else 4000
        # Degrade gracefully: a column this page asks for may not exist yet (a stage rebuilt
        # since, or still running). Showing the columns that DO exist and naming the absent ones
        # is more useful than failing the whole table, and makes the gap visible rather than
        # silently dropping a stage.
        have = list(pd.read_csv(p, sep="\t", nrows=0).columns)
        miss = [c for c in (cols or []) if c not in have]
        use = [c for c in (cols or have) if c in have] or None
        d = pd.read_csv(p, usecols=use, **kw)
        if where:
            d = d.query(where)
        note = str(path)
        if miss:
            note += f"   [not yet written: {', '.join(miss)}]"
        return d.head(n), note
    except Exception as e:                                    # noqa: BLE001
        return None, f"unreadable — {path} ({type(e).__name__})"


def table(df, note):
    if df is None:
        return f'<p class="missing">{html.escape(note)}</p>'
    h = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    rows = ""
    for _, r in df.iterrows():
        cells = "".join(
            f"<td>{html.escape(f'{v:.4g}' if isinstance(v, float) else str(v)[:40])}</td>"
            for v in r)
        rows += f"<tr>{cells}</tr>"
    return (f'<div class="src">{html.escape(note)}</div>'
            f'<div class="scroll"><table><thead><tr>{h}</tr></thead><tbody>{rows}</tbody></table></div>')


# ---------------------------------------------------------------- the stages
def stages():
    S = []

    S.append(dict(
        id="viper", n="1", title="Disease master regulators (upstream)",
        script="50_viper.R → group_mr.tsv",
        math=r"""<p>msviper over the cohort's ARACNe network, tumour vs its GTEx-matched normals.
        The output is one <code>nes<sub>i</sub></code> per regulator — a single number summarising
        how dysregulated that regulator's whole programme is.</p>
        <p class="key">This is the number the rest of the pipeline tries to reverse, and the number
        D202/C1 showed is <em>lossy</em>: within a direction of effect it leaves ~45% of the
        reproducible gene pattern undetermined.</p>""",
        ins=[], outs=[("results/step6/aracne_project_v2_tf/group_mr.tsv", None, None)]))

    S.append(dict(
        id="split", n="2", title="The aREA split — one NES into its genes",
        script="50b_regulon_contributions.R (disease) · 79e3 Net.contrib (drug)",
        math=r"""<p>Rank the profile's genes over the scored universe <code>U</code>:</p>
        <div class="eq">r<sub>j</sub> = rank<sub>j</sub> / (|U| + 1)<br>
        t2<sub>j</sub> = &Phi;<sup>&minus;1</sup>(r<sub>j</sub>) &nbsp;&nbsp;<span class="c">directional</span><br>
        t1<sub>j</sub> = &Phi;<sup>&minus;1</sup>(|2r<sub>j</sub> &minus; 1| + (1 &minus; max|2r &minus; 1|)/2)
        &nbsp;&nbsp;<span class="c">magnitude</span></div>
        <p>Per edge <code>(i,j)</code> with mode <code>m<sub>j</sub></code> and likelihood
        <code>l<sub>j</sub></code>:</p>
        <div class="eq">a<sub>j</sub> = m<sub>j</sub> l<sub>j</sub> t2<sub>j</sub> &nbsp;&nbsp;&nbsp;
        b<sub>j</sub> = (1 &minus; |m<sub>j</sub>|) l<sub>j</sub> t1<sub>j</sub><br>
        S1<sub>i</sub> = &Sigma;<sub>j</sub> a<sub>j</sub> &nbsp;&nbsp; S2<sub>i</sub> = &Sigma;<sub>j</sub> b<sub>j</sub>
        &nbsp;&nbsp; &sigma;<sub>i</sub> = sign(S1<sub>i</sub>)<br><br>
        <b>c<sub>ij</sub> = ( a<sub>j</sub> + &sigma;<sub>i</sub> b<sub>j</sub> [S2<sub>i</sub> &gt; 0] ) &middot; k<sub>i</sub></b>
        &nbsp;&nbsp; k<sub>i</sub> = w<sub>i</sub> / &Sigma;<sub>j</sub> l<sub>j</sub>,
        w<sub>i</sub> = &radic;(&Sigma;<sub>j</sub> (l<sub>j</sub>/l<sub>max</sub>)&sup2;)</div>
        <p><b>Drug side: exact</b> — condition-level VIPER has no null model, so
        &Sigma;<sub>j</sub> c<sub>ij</sub> <em>is</em> the activity (verified, max |diff| 4.8e&minus;07).
        <b>Disease side: an allocation</b> — msviper's permutation null cannot be replayed, so
        <code>d<sub>ij</sub> = c<sub>ij</sub> &middot; NES<sub>i</sub>/ES<sub>i</sub></code>.</p>""",
        ins=[("data/processed/step5/regulons_aracne_project_v2_tf_edges.tsv.gz",
              ["net", "regulator", "target", "tfmode", "likelihood"], None)],
        outs=[("results/step6/aracne_project_v2_tf/regulon_contributions.tsv.gz", None, None)]))

    S.append(dict(
        id="window", n="3", title="Timepoint windows — never pooled",
        script="systems_validation/integrate/window_reversal.py",
        math=r"""<div class="eq">early: pert_time &le; 6 h &nbsp;&nbsp;&nbsp; late: pert_time &ge; 24 h<br>
        <span class="c">`other` and the unrecorded sentinel (&minus;1) are DROPPED, never pooled</span></div>
        <p class="key">D183 measured ~80% of a window's reliable variance as absent from the other,
        and 37.1% of regulator terms flip sign between windows. Pooling averages a correct reading
        with an inverted one — for anti-estrogens in MCF7 the inverted window outvotes the correct
        one roughly 10:1 (D169).</p>""",
        ins=[("data/processed/lincs/drug_condition_profiles_index.tsv",
              ["pert_id", "ModelID", "pert_dose", "pert_time", "n_sig", "median_tas"], None)],
        outs=[("results/systems_validation/integrated/window_reversal.tsv.gz", None, None)]))

    S.append(dict(
        id="rev", n="4", title="Reversal — the production score",
        script="79e_drug_reversal.py:285",
        math=r"""<div class="eq"><b>rev = &minus;(D &middot; s) / &#8214;s&#8214;
        = &minus;&Sigma;<sub>i</sub> nes<sub>i</sub> &middot; act<sub>i</sub> / &#8214;s&#8214;</b></div>
        <p><code>D</code> = the drug's regulator activity, <code>s</code> = the disease NES vector.
        Each regulator's share is <code>rev_contrib<sub>i</sub></code>, and they sum to
        <code>rev</code> exactly.</p>
        <p class="key">Note what this multiplies: the <em>total</em> disease dysregulation by the
        <em>total</em> drug activity. It never asks whether the same genes are involved.</p>""",
        ins=[], outs=[(str(CDIR / "regulator_contributions_by_window.tsv.gz"),
                       ["cohort", "cmap_name", "window", "regulator", "disease_nes",
                        "drug_activity", "rev_contrib"], None)]))

    S.append(dict(
        id="revgene", n="5", title="rev_gene — the genes must match",
        script="79e3_window_contributions.py (COSINE=1)",
        math=r"""<div class="eq">
        rev &nbsp;&nbsp;&nbsp;&nbsp;= &minus;&Sigma;<sub>i</sub>
        (&Sigma;<sub>j</sub> d<sub>ij</sub>)(&Sigma;<sub>j</sub> C<sub>ij</sub>) / &#8214;s&#8214;
        &nbsp;&nbsp;<span class="c">a PRODUCT OF SUMS</span><br>
        <b>rev_gene = &minus;&Sigma;<sub>i</sub> &Sigma;<sub>j</sub>
        d<sub>ij</sub> &middot; C<sub>ij</sub> / &#8214;s&#8214;</b>
        &nbsp;&nbsp;<span class="c">a SUM OF PRODUCTS</span></div>
        <p>A regulator whose activity comes from genes orthogonal to the disease pattern has
        <code>&Sigma;<sub>j</sub> d<sub>ij</sub>C<sub>ij</sub> &asymp; 0</code> even when both
        totals are large — down-weighted <b>by geometry, with no threshold and no tuned
        parameter</b>.</p>
        <p class="warn"><b>Sign note (recorded because it was wrong once).</b>
        <code>d<sub>i</sub></code> already carries the disease direction, since
        &Sigma;<sub>j</sub> d<sub>ij</sub> = nes<sub>i</sub>. Multiplying by &minus;nes<sub>i</sub>
        as well double-counts it and flips the sign. There is <em>no</em> nes factor in
        rev_gene.</p>""",
        ins=[], outs=[(str(CDIR / "regulator_contributions_by_window.tsv.gz"),
                       ["cohort", "cmap_name", "window", "regulator", "disease_nes",
                        "drug_activity", "rev_contrib", "act_gene", "rev_gene_contrib",
                        "cos_contrib"], None)]))

    S.append(dict(
        id="cos", n="6", title="Opposition — per regulon, then per gene",
        script="79e3 (COSINE=1), D194/D195",
        math=r"""<div class="eq">
        cos_contrib<sub>i</sub> = &#10216;d<sub>i</sub>, C<sub>i</sub>&#10217; /
        (&#8214;d<sub>i</sub>&#8214; &#8214;C<sub>i</sub>&#8214;)<br>
        cos_term<sub>ij</sub> = d<sub>ij</sub> C<sub>ij</sub> /
        (&#8214;d<sub>i</sub>&#8214; &#8214;C<sub>i</sub>&#8214;)
        &nbsp;&nbsp;&nbsp; <b>&Sigma;<sub>j</sub> cos_term<sub>ij</sub> = cos_contrib<sub>i</sub></b>
        <span class="c">(exact to 1.6e&minus;06 over 61,000 regulons)</span></div>
        <p>Negative = the drug moves the very genes that account for that regulator's
        dysregulation, in the opposite direction.</p>
        <p class="warn">This divides by <code>&#8214;C<sub>i</sub>&#8214;</code> too, and that is
        the difference from rev_gene. D202 attributed the cosine score's whole ranking loss to it:
        <code>rev/&#8214;D&#8214;</code> cost <b>+0.004</b> AUC, adding the per-regulon
        normalisation cost <b>&minus;0.035</b>. Scale-free per regulon means a regulon the drug
        never moved still contributes on the full &plusmn;1 scale.</p>""",
        ins=[], outs=[(str(CDIR / "target_contributions_by_window.tsv.gz"),
                       ["cohort", "cmap_name", "window", "regulator", "target_symbol", "tfmode",
                        "likelihood", "disease_contrib_nes", "activity_contrib", "cos_term"], None)]))

    S.append(dict(
        id="prog", n="7", title="Programs — regulators that move together",
        script="58_regulator_modules.R → 79ai_programs.py / 79ak_drug_programs.py",
        math=r"""<div class="eq">average-linkage clustering on 1 &minus; |r|, cut at |r| &ge; 0.7<br>
        disease_score<sub>P</sub> = mean<sub>i&isin;P</sub> nes<sub>i</sub> &nbsp;&nbsp;
        drug_score<sub>P</sub> = mean<sub>i&isin;P</sub> act<sub>i</sub><br>
        opposition<sub>P</sub> = &minus; disease_score<sub>P</sub> &middot; drug_score<sub>P</sub> /
        |disease_score<sub>P</sub>|</div>
        <p>Disease programs cluster regulators over <b>patients</b>; drug programs over
        <b>drugs</b>, per window. 461 and 1,035 respectively.</p>
        <p class="warn">Programs span only <b>27.4%</b> of a cohort's regulators (clean ones
        14.1%), which is why collapsing to them <em>lost</em> ranking power (D201).</p>""",
        ins=[], outs=[("results/lincs/programs/programs.tsv",
                       ["cohort", "program", "n_regulators", "hub", "disease_score", "hallmark",
                        "is_proliferation", "is_technical"], None)]))

    S.append(dict(
        id="ceiling", n="8", title="The reliability ceiling — what bounds every comparison",
        script="79e3 (SPLITHALF=1), D202",
        math=r"""<div class="eq">R = median<sub>i</sub> cosine( C<sub>i</sub>(half A),
        C<sub>i</sub>(half B) )</div>
        <p>Split one drug's conditions in half, build each half's gene vector independently, and
        compare them. A drug against <em>itself</em>.</p>
        <p class="key">R = <b>0.531</b> overall (<b>0.376</b> early, <b>0.590</b> late), and it is
        sample-size driven: 0.326 at 2 conditions &rarr; <b>0.793 at 10+</b> (Spearman +0.528).
        Nothing can match a vector better than it matches itself, so every similarity is reported
        as a fraction of R, never in absolute terms.</p>""",
        ins=[], outs=[(str(CDIR / "splithalf_reliability.tsv.gz"),
                       ["cohort", "pert_id", "window", "n_conditions", "regulator", "cos_half"], None)]))

    return S


CSS = """
:root{--bg:#fbfbfa;--fg:#1b1b1a;--mut:#6b6b68;--line:#e3e3e0;--card:#fff;--acc:#1f6f6b;
--warn:#8a4b1f;--warnbg:#fdf6ef;--keybg:#f1f7f6}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#17181a;--fg:#e8e8e6;
--mut:#9b9b97;--line:#2e3034;--card:#1e2023;--acc:#5fbdb6;--warn:#d79a62;--warnbg:#241c14;
--keybg:#17221f}}
:root[data-theme=dark]{--bg:#17181a;--fg:#e8e8e6;--mut:#9b9b97;--line:#2e3034;--card:#1e2023;
--acc:#5fbdb6;--warn:#d79a62;--warnbg:#241c14;--keybg:#17221f}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.6 ui-sans-serif,system-ui,-apple-system,Segoe UI,Roboto,sans-serif}
.wrap{max-width:1100px;margin:0 auto;padding-left:max(16px,env(safe-area-inset-left,0px));padding-right:max(16px,env(safe-area-inset-right,0px))}
header{border-bottom:1px solid var(--line);padding:32px 0 20px;margin-bottom:24px}
h1{margin:0 0 6px;font-size:1.6rem;letter-spacing:-.01em}
.sub{color:var(--mut);margin:0;max-width:70ch}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:18px 20px;margin:0 0 20px}
.card h2{margin:0 0 2px;font-size:1.12rem}
.num{display:inline-block;min-width:1.7em;color:var(--acc);font-variant-numeric:tabular-nums;font-weight:600}
.script{color:var(--mut);font-size:.84rem;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;margin:0 0 12px}
.eq{background:var(--keybg);border-left:3px solid var(--acc);padding:12px 14px;margin:10px 0;
font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.92rem;line-height:1.9;overflow-x:auto}
.eq b{color:var(--acc)}
.c{color:var(--mut);font-style:italic}
.key{background:var(--keybg);border-radius:6px;padding:10px 12px;margin:10px 0}
.warn{background:var(--warnbg);border-left:3px solid var(--warn);border-radius:0 6px 6px 0;padding:10px 12px;margin:10px 0}
h3{font-size:.78rem;text-transform:uppercase;letter-spacing:.07em;color:var(--mut);margin:18px 0 6px}
.src{font-family:ui-monospace,Menlo,monospace;font-size:.76rem;color:var(--mut);margin-bottom:4px;word-break:break-all}
.scroll{overflow-x:auto;border:1px solid var(--line);border-radius:6px}
table{border-collapse:collapse;width:100%;font-size:.8rem;font-variant-numeric:tabular-nums}
th{background:var(--keybg);text-align:left;padding:6px 9px;white-space:nowrap;font-weight:600;border-bottom:1px solid var(--line)}
td{padding:5px 9px;white-space:nowrap;border-bottom:1px solid var(--line);color:var(--fg)}
tr:last-child td{border-bottom:0}
.missing{color:var(--warn);font-size:.85rem;font-style:italic}
footer{color:var(--mut);font-size:.82rem;border-top:1px solid var(--line);padding:18px 0 40px;margin-top:10px}
@media(max-width:560px){h1{font-size:1.3rem}.card{padding:14px 15px}}
"""


def main() -> int:
    parts = []
    for s in stages():
        body = [f'<section class="card" id="{s["id"]}">',
                f'<h2><span class="num">{s["n"]}</span>{html.escape(s["title"])}</h2>',
                f'<p class="script">{html.escape(s["script"])}</p>']
        if s["ins"]:
            body.append("<h3>Data in</h3>")
            for path, cols, where in s["ins"]:
                body.append(table(*snap(path, cols, where=where)))
        body.append("<h3>The transformation</h3>")
        body.append(s["math"])
        if s["outs"]:
            body.append("<h3>Data out</h3>")
            for path, cols, where in s["outs"]:
                body.append(table(*snap(path, cols, where=where)))
        body.append("</section>")
        parts.append("".join(body))

    # ARTIFACT=1 drops the document skeleton: the Artifact platform supplies doctype/head/body,
    # and a page that ships its own renders in quirks mode. One generator, two wrappers, so the
    # published page cannot drift from the repo one.
    art = os.environ.get("ARTIFACT") == "1"
    head = ('<title>Reversal Pipeline Walkthrough</title>\n'
            '<meta name="claude:viewport-fit-sides" content="cover">\n'
            f'<style>{CSS}</style>\n') if art else (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width,initial-scale=1">\n'
        f'<title>Reversal Pipeline Walkthrough</title><style>{CSS}</style></head><body>\n')
    tail = "" if art else "\n</body></html>"
    page = f"""{head}<header><div class="wrap"><h1>Pipeline walkthrough — data in, the math, data out</h1>
<p class="sub">Every table below is read from the actual file the pipeline wrote, so this page
cannot drift from what was produced; a stage that has not been built says so rather than showing a
plausible example. Formulas are those in <code>docs/FORMULAS.md</code>, which is checked against
the implementations.</p></div></header>
<main class="wrap">{''.join(parts)}</main>
<footer><div class="wrap">Generated by <code>scripts/79x3_pipeline_explainer.py</code>.
Stages upstream of this walkthrough (ARACNe networks, VIPER, drug VIPER) are taken as given —
see <code>docs/PIPELINE_RUN_ORDER.md</code> for the boundary and the full run order.</div></footer>{tail}"""
    out = (RES / "pipeline_explainer_artifact.html") if art else OUT
    out.write_text(page)
    print(f"wrote {out}  ({len(page)/1e6:.2f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
