# v1 Performance Profile

Baseline profiling of trnascan-py (thin Infernal binding) vs reference
tRNAscan-SE 2.0, showing where wall-clock goes and how the isotype optimization
made it faster than the reference.

Reproduce with:

```bash
python scripts/profile_v1.py /path/to/genome.fa --domain bact
```

## Benchmark

- **Genome:** *Escherichia coli* K-12 MG1655 (`NC_000913.3`), 4,641,652 bp, 1 contig
- **Domain:** bacterial (`-B` / `bact` models)
- **Hardware:** shared dev box, default Infernal threading (`--cpu` unset)
- **Date:** 2026-06-17

## Initial baseline (before optimization)

> Superseded — see [After](#after-e-coli-k-12-same-box) below for the current
> result. Kept to document where the time went.

| Metric | trnascan-py v1 | reference tRNAscan-SE 2.0 |
|---|---|---|
| tRNAs called | **89** | **89** |
| End-to-end wall-clock | 30.85 s | 25.31 s |
| Throughput | 150 kbp/s | 183 kbp/s |

The initial v1 was **0.82×** the reference's speed (~22% slower) — it ran the
first-pass and isotype scans as two independent whole-genome passes.
**Detection count is identical (89 = 89).**

### Stage breakdown (standalone timing)

| Stage | Time | Notes |
|---|---|---|
| `cmsearch` first-pass + confirm + align (1 general CM) | **1.76 s** | HMM-filtered, single model — cheap |
| `cmscan` isotype classification (~22 isotype CMs) | **34.50 s** | whole-genome × all isotype models — **the hot stage** |

> The first-pass general-model search is fast; the **isotype multi-model scan
> dominates wall-clock**. Essentially **100% of runtime is covariance-model DP**
> (both stages are CM search) — the isotype stage is just ~95% of that.

## Genome-scale extrapolation (v1, this hardware)

| Genome | Extrapolated `scan()` |
|---|---|
| 100 Mb | ~11 min |
| 1 Gb | ~1.85 h |
| Human 3.1 Gb | **~5.7 h** |

## Optimization applied — v1 now beats the reference

The isotype stage was the bottleneck. Two changes removed it:

1. **Anticodon-derived isotype.** A tRNA's amino acid is determined by its
   anticodon, which we already extract from the first-pass alignment. Translating
   it via the genetic code is essentially free and replaces the whole-genome
   isotype `cmscan`.
2. **CM refinement only for ambiguous loci.** Anticodons the genetic code can't
   disambiguate (`CAT` → Met/fMet/Ile2/iMet, `TCA` → SeC) are refined with a
   small CM scan of just those ~handful of loci.

This is the default `isotype_method="hybrid"`. (`"anticodon"` = fastest,
`"cm"` = scan every locus.)

### After (E. coli K-12, same box)

| Metric | v1 (hybrid) | reference |
|---|---|---|
| tRNAs called | 89 | 89 |
| Wall-clock | **15.4 s** | 26.1 s |
| **Speed ratio** | **1.70× faster** | 1.0× |
| Isotype concordance | 87/89* | — |

*The 2 differences are loci where the **reference** reports `Undet` (its own
anticodon detection returned `NNN`) while v1 extracts a valid anticodon and
assigns the correct isotype — i.e. v1 is at least as complete. The genuine
ambiguous cases (`CAT` initiator/Ile2) are resolved correctly by the CM
refinement.

`isotype_method="anticodon"` (skipping all CM) runs the isotype step in ~0s
(scan ≈ 1.2 s total, ~20× faster than reference) at 81/89 concordance.

## Implications

After this optimization the dominant remaining cost is the first-pass
general-model `cmsearch` — pure covariance-model CYK/Inside DP. A native GPU
re-implementation of that kernel was evaluated as a possible future direction and
**parked**: the CPU speedup here already covers typical genome-scale use, and the
remaining hot path is Infernal's already-tuned CM search (see the README
Roadmap). It would only be worth revisiting for routine scanning of many
multi-gigabase genomes where even the first pass dominates wall-clock.
