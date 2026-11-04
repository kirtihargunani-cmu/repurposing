# Disease side in plain language, with results

Written 2026-11-02 from the v2 results. This is the non-technical version of [`PIPELINE_V2.md`](PIPELINE_V2.md), which
has the commands, files and thresholds. The drug side continues in
[`DRUG_SIDE_PLAIN_LANGUAGE.md`](DRUG_SIDE_PLAIN_LANGUAGE.md). Decision numbers (D..) point to
[`DECISIONS.md`](DECISIONS.md).

## What we are trying to do

For each cancer, find which **control switches** (transcription factors: proteins that turn groups of genes on or
off) are stuck on or off compared with healthy tissue of the same organ. These switches, the cancer's **master
regulators**, are what the drug side later tries to push back toward healthy.

## How it works, step by step

1. **Gather patients.** We use public RNA-seq for tumors (TCGA, TARGET) and healthy organs from deceased donors (GTEx),
   all processed the same way. Mixed cancers are split into known subtypes, for example breast cancer into luminal A,
   luminal B, HER2 and basal. Groups with fewer than 30 patients are dropped. That leaves **49 groups and 10,211
   patients**.
2. **Keep only genes we can measure reliably.** Some genes look so much like other genes that sequencing reads cannot
   tell them apart. Others are "readthrough" transcripts that run across two genes. Both create false signals.
   **2,013 of 19,797 protein-coding genes are dropped; 17,784 are kept.**
3. **Prepare a fair healthy comparison.**
   - **Use the right part of the organ.** Colon samples that are mostly muscle wall and kidney samples from the wrong
     region are removed: colon goes from 165 to 109 samples, kidney from 28 to 21.
   - **Correct for how the tissue was collected.** Donor tissue decays depending on how the person died and how long
     the tissue waited. We adjust every sample to look like the best-preserved donors.
4. **Measure what changed in each tumor (the signature).** For every gene, how much higher or lower is it in the
   tumors than in healthy tissue? Age, sex and tumor purity are taken into account. Tumors vary more than healthy
   samples (1.2-1.8 times), so each side gets its own noise model (D69). **41 groups get a usable signature.**
5. **Learn which genes each switch controls (networks).** In each cancer, genes that consistently rise and fall
   together with a switch across patients are taken as its targets (ARACNe). We build **29 networks** and use only
   real transcription factors as switches, about 1,100-1,360 per network. Links between look-alike genes are removed,
   11-15% of links, because shared reads make them look connected.
6. **Score each switch (VIPER).** If a switch's targets are all shifted the same way in the tumor, the switch is
   probably on (or off). Each switch gets a score (NES) per cancer and per patient. To know what counts as chance, we
   rerun the scoring 1,000 times on random halves of the healthy samples.
7. **Check the results make biological sense.** Switches with known behavior should come out right. For example,
   cell-division switches (E2F1, FOXM1) should be up in tumors and estrogen-receptor switches up in luminal breast
   cancer.
8. **Look for subgroups of patients within each cancer.** Patients are clustered by their switch scores. A split is
   kept only if it reproduces in random halves of the patients and each subgroup has at least 30 patients. It also
   must not be explained by lab batch, tumor purity or RNA quality (D71).
9. **Give each kept subgroup its own signature and switch scores**, so drugs can be matched to the subgroup as well as
   to the whole cancer.

## Results

### Master regulators per cancer

- **Typical cancer: 151 switches clearly changed** (median; 116 at the stricter cutoff). The "confirmed" list has a
  median of 117 per cancer. A few cancers have almost none: cervical adenocarcinoma 0; bladder, stomach and uterine
  carcinosarcoma 1 each; head and neck 4.
- **Known-behavior checks:**
  - 77 come out right and clearly significant.
  - 55 come out in the right direction but below the cutoff.
  - 4 come out in the wrong direction, all with scores near zero, so none is a real contradiction.
- **Which switches pass a cutoff is fragile, but their order is not.** Changing the chance model moved colon cancer
  from 425 to 0 significant switches, while the ranking of switches stayed the same (correlation at least 0.98). So the
  drug side uses every switch, weighted by its score, instead of a cut list (D74).
- **Two ways of defining targets agree poorly.** The tumor-learned targets and a curated list of direct targets
  (CollecTRI) rank switches differently (median correlation 0.18). The drug side therefore averages both (D75).

### What could distort the signatures, and what we found

| Concern | Finding | What we did |
|---|---|---|
| Decay of donor tissue | After correction the post-mortem signal is closer to zero in 28 of 39 groups | Correction kept (D57) |
| Different cell mix in tumors vs healthy organ | Cell-type make-up explains a median 19% of each signature (11-38%) | Measured, not corrected (open) |
| A recurring block of ~50-80 switches that moves together | Present in 36 of 41 cancers; tracks lab batch, not biology | Left out when defining subgroups (D79; PI to confirm) |
| Network links that only reflect tumor type ("identity", not signaling) | 0-7% of links overall, but about half for lineage switches: breast ESR1 51%, esophagus HNF1A / FOXA2 up to 89% | Lists without them tested: small gain, see below (D82, D83) |
| A separate luminal-only breast network | Changes luminal A / B scores no more than random networks of the same size | Kept the pooled breast network (D81) |

Removing the identity links helped anti-estrogens only a little on the drug side: luminal A rose from the 9.7th to
the 12.3rd percentile and luminal B from the 23.4th to the 26.0th. Whether to switch to these lists is still open.

### Subgroups

20 of 41 cancers show significant structure. Six pass every keep rule, and colon cancer is kept as an exception: it
narrowly misses (0.748 vs 0.75) but matches the known MSI vs CIN split and replicated in two external cohorts in v1.
That gives **14 subgroups in 7 cancers**:

| Cancer | Subgroups (patients) |
|---|---|
| Bladder | basal-squamous-like (172) / luminal-like (235) |
| Head and neck | basal-mesenchymal-like (363) / atypical-classical-like (155) |
| Low-grade glioma, IDH-mutant non-codeleted | LGr2 (89) / LGr3 (152) |
| Pheochromocytoma | pseudohypoxia-Wnt (69) / kinase-signaling (77) |
| Thyroid | BRAF-like (296) / RAS-like (208) |
| Seminoma | C1 (32) / C2 (30) |
| Colon (exception) | CIN-like (194) / MSI-like (93) |

**Subgroups are more uniform than whole cancers.** 13 of 14 subgroups vary less from patient to patient than their
parent cancer (53-94% of its spread; the exception is glioma LGr2, 107%). 11 of 14 are as uniform as the most
uniform quarter of all 41 cancers. Patients resemble their own subgroup more than the sibling one in 85-100% of
cases. In the two seminoma subgroups and colon CIN-like, the main remaining difference between patients partly
follows lab batch.

The names come from the known TCGA subtypes of their patients. The glioma and seminoma pairs barely differ in drug
response, so the drug side reports one list for each of them.

### Separate networks for different cell types (2026-11-03)

Each cancer used to get one network, learned from all its tumors together, even when the cancer contains clearly
different kinds of cells (luminal vs basal breast cancer, squamous vs adenocarcinoma of the esophagus, the glioma
subtypes). In such a mix the network mostly learns which genes go with which cell type. Now every subtype with 85 or
more tumors gets its own network (12 subtypes in 8 cancers); smaller subtypes keep the shared one (D97). The
master-regulator lists of the big subtypes barely change (luminal A / B, glioma: 0.97-0.99 agreement with before);
those of mixed or small ones change more (basal breast 0.90). In basal breast cancer, low "luminal" regulators
(ESR1, GATA3) no longer dominate its list, which is what had made anti-estrogens look useful there. Smaller networks
are noisier: the overlap between two halves of a network's runs is 0.40 at 85 tumors vs about 0.6 for the largest
networks. The patient-level analyses and the subgroup clustering have not yet been rerun on the new networks.

## What is not done yet

- **Replication of the new subgroups** in other patient cohorts (cohorts to choose).
- **For the PI:** the recurring block (D26) and the kidney reference, which has only 21 donors (D60).
- **For the user:**
  - whether to adopt the identity-free target lists (D83);
  - whether to keep the "confirmed" rule (D18), whose curated part adds almost nothing;
  - the glioma purity measure (A37, tabled).
