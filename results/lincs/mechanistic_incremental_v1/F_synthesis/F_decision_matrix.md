# Decision matrix — mechanistic_incremental_v1

Resolved = the frozen rule (95% interval excludes the matched null). Every figure is from the stage's committed outputs
on the 2026-11-23 run.

| Stage | Question | Verdict | Key numbers | Next action (spec section 13, as amended) |
|---|---|---|---|---|
| 0 Audit | Are the inputs the latest run, joins safe, decompositions exact? | **Pass** | 10/10 run ids; regulator sums exact to 5.7e-5; `rev_gene` to 2.1e-6 | none |
| A Resolvable signal | Does each processing step carry drug-specific, mechanism-consistent signal? | **Resolved** at every step, both windows | same- vs different-drug cosine 0.35-0.62; `rev_gene` reliability 0.61 early / 0.74 late; PROGENy on-target 0.76-0.77 vs 0.64-0.65; shared target +0.28 to +0.44 | keep the production score under study; signal is there to build on |
| B Temporal | Do early and late differ beyond noise, and by mechanism? | **Resolved**: delayed responses, class time profiles; reversals rare but real | calibration 95.3%; late-only 15.8-26.7% vs 11.3%; PI3K late (OR 4.7), CDK persistent (5.0), HDAC early (3.4), HSP early (5.8), FDR 0.023 | keep windows separate; a prospective time-course for the four classes |
| C Target recovery | Can targets be recovered from the regulator response? | **Resolved** for primary targets at 24 h and 2-hop at 6 h; **24 h inverts network signs** | primary exact 58 vs 52 pct (p 0.001); 2-hop 52.2 vs 51.2; signed 1-hop late 45 vs 49 (cancers), 42 vs 49 (subgroups) | test feedback at 24 h (phospho / time-course); secondary targets need more regulator-type targets |
| D Network sensitivity | Why do networks read a regulator differently? | **Resolved**: rewiring, explained by target overlap | overlap vs agreement 0.68-0.84 in 16/16 pair-windows; 83-95% of difference on unique targets; ESR1 overlap 1.4% | flag rewired regulators in reporting; keep per-cohort networks |
| E Program access | Are disease programs specially perturbable? | **Not resolved** beyond matched random sets; responses reproducible and drug-driven | reliability 0.61 / 0.76 vs pseudo 0.64 / 0.77; 3/474 above pseudo q95; line share 2.5% | programs stay interpretive; do not rank on them |
| F Synthesis | What evidence does each eligible pair carry? | 1,497 pairs assembled; no reranking | perturbation supported 27%; temporal supported 7-11% (not tested 43-59%); target supported 8-10% | use as the evidence card per candidate |

Not run, by amendment: PRISM/GDSC selective killing (A1.4). Not run, flagged for follow-up: B4 class enrichment on the
pseudo-replicate split and the high-confidence arm; genetic-knockdown controls through each network (D).
