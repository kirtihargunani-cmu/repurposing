# Drug side: caveats to revisit during ranking

Open issues on the drug side that are noted but not acted on. Revisit each one if it turns out to matter when
drugs are ranked. Newest first.

## LumA PC3 low end gets its drug data mostly from HER2-type cell lines (2026-10-12)

**Status:** noted, no change. Revisit if it distorts the breast stratum rankings (user, 2026-10-12).

**The axis.** LumA PC3 splits luminal A patients by proliferation:
- High end: fast-growing tumors (AURKA, MYBL2 high).
- Low end: slow-growing tumors with strong estrogen-receptor activity (ESR1, GATA3 high), the most classically
  ER-driven luminal A tumors.

**The issue.** Each stratum's drug data comes from the cell lines its patients most resemble (Celligner weights,
`scripts/74b_stratum_line_weights.py`). Only six breast lines have LINCS drug data: MCF7, BT-474, SK-BR-3, BT-20,
MDA-MB-231 and Hs 578T. Among those six:

| Stratum | Leans mainly on |
|---|---|
| `BRCA-LumA-PC3-high` (fast-growing) | MCF7 (46%), ER-positive |
| `BRCA-LumA-PC3-low` (ER-driven) | BT-474 and SK-BR-3 (80%), both HER2-amplified; SK-BR-3 is ER-negative |

So the most ER-driven patients get their drug data mostly from HER2-type cells, one of them without the estrogen
receptor. That is the opposite of what the biology suggests.

**Likely cause.** The low end's closest lines overall are mostly lines without drug data (MFM-223, MDA-MB-453,
CAMA-1). Only 8% of its in-tissue weight falls on the six drug-tested lines (15% for the high end), so the split
among them rests on a thin slice of similarity. MCF7 is a fast-growing line, which likely ties it to the high end.

**How much it matters.** A drug's score only uses the lines where that drug was tested. Many drugs were tested only
in MCF7, so both ends use MCF7 for them and nothing changes. It matters for drugs tested in several breast lines:
for those, the low end's score leans on HER2-type cells.

**Possible fix if needed.** Let ER-driven strata use only ER-positive drug-tested lines (MCF7, BT-474).

Numbers: `results/celligner/stratum_line_weights_long.tsv` (`in_key` lines), `stratum_line_summary.tsv`.
