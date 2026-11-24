# Final report — how much mechanism signal the production scorer resolves

Study `mechanistic_incremental_v1`, branch `study/mechanistic-incremental-v1`. Protocol `5cc2227` (spec Amendment 1 +
frozen config), stage settings frozen before each stage: B `2e55c43`, C `cbd2fb5`, D `d3763b7`, E `6583b74`, F
`6c29417`. Inputs: the 2026-11-23 full run of `scripts/run_pipeline.sh` only (Stage 0 verified run ids, checksums,
joins and exact decompositions). Scorer under study: production D211, target evidence T >= 50 then gene-level regulon
reversal `rev_gene`, early (<= 6 h) and late (>= 24 h) kept apart. No efficacy endpoint (no PRISM/GDSC). Production is
unchanged.

## The answer

**The perturbation data carry reproducible, drug-specific and mechanism-consistent signal at every step of the
production scorer, with a structured time course and network dependence that can be explained.** Given that 52% of eligible
pairs rest on a single condition and only 29% of (cancer, drug) pairs have both windows, that is a strong basis for
further work.

1. **Signal survives every processing step (A).** Two halves of one drug's conditions agree far more than two drugs in
   the same cancer (cosine gain 0.35 -> 0.43 early, 0.56 -> 0.62 late from gene profile to network activity); the
   production score is reliable (0.61 early, 0.74 late; ceilings 0.78 / 0.86).
2. **The signal is mechanistic (A, C).** Drugs move their on-target PROGENy pathway the expected way 76-77% of the time
   (null 64-65%); drugs sharing a target look alike (+0.28 to +0.44 cosine); a primary target that is itself a
   regulator ranks high in the drug's own response at 24 h (58th vs 52nd percentile), and its 2-hop neighbourhood moves
   at 6 h.
3. **Time matters, and by mechanism (B, C).** Under a calibrated noise model (95.3%), responses that appear only by 24 h
   exceed the replicate rate in every arm; PI3K inhibitors act late, CDK inhibitors persist, HDAC and HSP90 inhibitors
   act early (FDR 0.023). At 24 h, downstream regulators move against the signalling network's expected direction
   (45 vs 49, cancers; 42 vs 49, subgroups), the pattern adaptive feedback would produce and the one the anti-estrogens
   showed in MCF7 (clomifene in luminal A breast is a confident sign reversal).
4. **Network differences are rewiring, not noise (D).** Where two networks read a regulator differently, they give it
   different targets (overlap vs agreement 0.68-0.84 in all 16 comparisons, on held-out drugs); ESR1 between luminal A
   and basal breast shares 1.4% of its targets. This supports one network per cohort and a per-regulator caveat.

## What is not resolved at current data depth

- Disease programs are reproducible and drug-driven but not more perturbable than matched random regulator sets, and
  target links do not predict which programs a drug moves (E). Programs stay interpretive.
- Secondary targets with measured binding: only 15-16 are scored regulators; downstream recovery reaches the null
  (C). Annotation-only targets carry no recoverable signal.
- Confident sign reversals in the strict matched-dose arm; subgroup arms with few units (B, A).
- The D211 biomarker reading (`rev_gene` within the eligible set ranks EGFR inhibitors higher in LUAD than LUSC at
  24 h for 2 of 8 compounds; the biomarker signal sits in the filter) is unchanged by this study, which measured
  mechanism signal, not ranking quality.

## What it does not show

No result here is an efficacy claim, and none ranks drugs better or worse than another score. The time-course and
feedback readings are hypotheses for an experiment (time-resolved transcriptomics or phosphoproteomics), not
established signalling.

## Evidence per candidate (F)

`F_mechanism_evidence.tsv.gz` gives each of the 1,497 eligible (cancer, drug, window) pairs a status per layer
(disease, perturbation, temporal, target, network, program). Perturbation is supported for 27% of pairs and conflicting
for 24%; temporal evidence is supported for 7% (early) and 11% (late) and not tested for 43-59% (no second window);
target evidence is supported for 8-10%. `F_mechanism_cases.md` shows one rule-chosen case per category and window; the
frozen tie rule (alphabetical at equal percentile) means several categories show drugs at the 100th percentile with
names early in the alphabet.

## Reproduce

```
git checkout study/mechanistic-incremental-v1
python scripts/mechanistic_incremental_v1/00_audit.py
python scripts/mechanistic_incremental_v1/A_resolvable.py
python scripts/mechanistic_incremental_v1/B_temporal.py
python scripts/mechanistic_incremental_v1/C_target_recovery.py
python scripts/mechanistic_incremental_v1/D_network_sensitivity.py
python scripts/mechanistic_incremental_v1/E0_programs.py
python scripts/mechanistic_incremental_v1/E_program_access.py
python scripts/mechanistic_incremental_v1/F_synthesis.py
```

Interpreter `~/miniforge3/envs/celligner/bin/python`; each was run through `scripts/runlog.sh`; inputs are read from
`input_root` in `config/mechanistic_incremental_v1.yaml`. Implementation fixes made during the study (each recorded in
its stage report, none changing a rule or threshold): PROGENy resolved rule aligned to the frozen rule (A), odds-ratio
correction removed (B), two-sided check of C's summary, sorted 2-hop cap (C), empty-family guard (D).
