# Exploration scripts (2026-11-02, not decisions)

One-off analyses behind the PROGRESS entry "exploration, can ARACNe edges be corrected?" and the CollecTRI edge test.
They are kept as a record, not as pipeline steps. Intermediate files were written to a session scratch directory
(the `SCR` path at the top of each script): set `SCR` to a writable folder before rerunning. Order:
`edge_explore.R` -> `edge_paired.R`, `edge_frac.R` -> `variant_eval.R` + `variant_summ.R` -> `hidden_factor.R` ->
`variant_eval_hf.R` -> `hf_regulator.R`; `collectri_edges.R` stands alone. `variant_eval*.R` use 28 cores (`NC`);
lower it (user: under half the machine's threads).
