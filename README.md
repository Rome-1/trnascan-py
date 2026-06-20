# trnascan-py

**A fast Python interface for tRNA gene detection.**

Transfer-RNA genes are short (~70–90 nt) structured genes scattered across
genomes; finding them accurately means scoring candidate loci against a
*covariance model* (a profile that captures both sequence and the tRNA cloverleaf
secondary structure). The standard tool for this is
[tRNAscan-SE 2.0](http://lowelab.ucsc.edu/tRNAscan-SE/) (UCSC Lowe Lab), which
drives [Infernal](http://eddylab.org/infernal/) covariance-model search.

trnascan-py is a clean Python CLI and library that drives the same Infernal search
with the same covariance models, but with a leaner pipeline. Its calls
(coordinates, isotype, anticodon, introns, pseudogene flags) are validated
differentially against reference tRNAscan-SE 2.0 across bacteria, archaea,
eukaryotic-nuclear, vertebrate-mito, and plant-organellar genomes, with high
concordance (a few characterized edge-case differences, mostly a few-bp boundary
shifts — see [Performance](#performance) and [docs/BENCHMARK.md](docs/BENCHMARK.md)),
and runs **2.6–9.4× faster on genomes ≥ 150 kb**.

> **Status:** Alpha (v0.1). API may change before 1.0.

## What makes it faster

Same engine (Infernal), fewer passes. Two changes account for most of the
speedup:

1. **Anticodon-derived isotype.** A tRNA's amino acid is fixed by its anticodon,
   which we read directly from the first-pass `cmsearch` alignment and translate
   via the genetic code. That replaces reference tRNAscan-SE's separate
   *whole-genome* isotype covariance scan — the dominant cost — with a near-free
   table lookup. Only the few genetically-ambiguous anticodons (`CAT` →
   Met/fMet/Ile2/iMet, `TCA` → SeC) get a small covariance scan of just those loci.
2. **A single lean first pass.** One `cmsearch` over the genome with the general
   model, with candidate loci pulled via a streaming FASTA index (the genome is
   never held in memory). Reference bacterial/archaeal modes run in
   max-sensitivity (no HMM filter), which is why large prokaryotic genomes show
   the largest speedups.

## Performance

Differential benchmark vs reference tRNAscan-SE 2.0 (shared dev box, default
threading). Detection counts match within one call per genome and isotype/anticodon
concordance is high; the few differences are characterized (mostly a few-bp
boundary shift just outside the match window, plus cases where the reference itself
returns `Undet`). Full table, methodology, and the discrepancy analysis are in
[docs/BENCHMARK.md](docs/BENCHMARK.md).

| Genome | Size | tRNAs (ours / ref) | Speed |
|---|---|---|---|
| *S. cerevisiae* chr I (euk) | 230 kb | 4 / 4 | **9.4×** |
| *M. jannaschii* (archaea) | 1.66 Mb | 38 / 37 | **3.3×** |
| *E. coli* K-12 (bacteria) | 4.64 Mb | 89 / 89 | **2.6×** |
| *Arabidopsis* chloroplast | 154 kb | 30 / 30 | 1.2× |
| Human mtDNA (vert-mito) | 16.5 kb | 22 / 22 | 0.6× † |

† On a tiny (~16 kb) mitochondrial genome the mito two-pass has fixed overhead
that dominates; it is not a scaling effect. Median speedup across the corpus is
**2.6×**.

## How it works

1. **First pass + CM confirmation** — `cmsearch` with the general (isotype-merged)
   covariance model; hits at or above the bit-score cutoff are kept.
2. **Anticodon extraction** — the anticodon triplet is read from the `cmsearch`
   alignment (the anticodon loop is the second cloverleaf hairpin).
3. **Isotype assignment** — anticodon → amino acid via the genetic code, with a
   small covariance scan to disambiguate `CAT`/`TCA`. Selectable via
   `--isotype-method {hybrid,anticodon,cm}` (default `hybrid`).
4. **Introns, suppressors, pseudogenes** — introns are read from the alignment
   insert columns; stop-decoding anticodons are reported as `Sup`; low-scoring
   structurally-degenerate hits are flagged `pseudo` using the reference's exact
   primary/secondary score criterion.

## Requirements

trnascan-py is pure-Python glue; the search engine and covariance models are
**external runtime dependencies** (not bundled — see [License](#license)):

- Python ≥ 3.10
- [Infernal](http://eddylab.org/infernal/) 1.1.x (`cmsearch` on `PATH`)
- [tRNAscan-SE 2.0](http://lowelab.ucsc.edu/tRNAscan-SE/) — provides the bundled
  tRNA covariance models trnascan-py searches with, and serves as the differential
  oracle.

Install the engines (Debian/Ubuntu), or via conda:

```bash
sudo apt-get install -y infernal trnascan-se
# or: conda install -c bioconda infernal trnascan-se
```

If the models are not auto-discovered, set `TRNASCAN_MODELS_DIR` to the directory
containing `TRNAinf-*.cm` (e.g. `/usr/share/trnascan-se/models`).

## Install

```bash
pip install git+https://github.com/Rome-1/trnascan-py.git   # latest
# or, from a checkout:
pip install .
```

For development (editable install with test/lint tooling), see
[Development](#development).

## Usage

```bash
trnascan-py scan genome.fa --domain bact            # scan for tRNAs
trnascan-py scan genome.fa --domain bact --json     # JSON output
trnascan-py scan genome.fa -d bact --isotype-method anticodon   # fastest
trnascan-py version                                 # show engine versions
```

```python
from trnascan_py import scan

for hit in scan("genome.fa", domain="bact"):
    print(hit.seq_id, hit.start, hit.end, hit.strand.value,
          hit.isotype, hit.anticodon, hit.score, hit.note)
```

`--domain` selects the model set: `euk`, `bact`, `arch`, or a mitochondrial
lineage `mito-vert` / `mito-mammal`.

## Scope & known limitations (v1)

- **Engines required.** Without Infernal + the tRNAscan-SE covariance models,
  trnascan-py cannot run — it orchestrates these tools.
- **Isotype edge cases.** `hybrid`/`cm` resolve initiator vs elongator Met, Ile2,
  and SeC; `anticodon` mode does not (it reports `Met` for all `CAT`).
- **Mito anticodons** are reported only for Leu/Ser (encoded in the model name);
  general mito anticodon extraction is a follow-up. Detection/coordinates/isotype
  match the reference. Only vertebrate/mammalian mito lineages are wired up.
- Introns are detected and reported (coordinates match the reference); candidate
  loci are streamed from a FASTA index (no whole-genome memory load).

Open work and follow-ups are tracked in the GitHub issue tracker.

## Development

```bash
pip install -e ".[dev]"            # editable install + ruff/mypy/pytest
ruff check . && mypy
pytest -m "not differential"       # pure unit tests (no external tools)
pytest                             # full suite; runs differential tests too when
                                   # cmsearch + tRNAscan-SE are on PATH
python scripts/corpus_benchmark.py # multi-genome differential + speed benchmark
```

See [CONTRIBUTING.md](CONTRIBUTING.md).

## Roadmap

v1 — a correct, fast, well-tested interface and the differential oracle — is the
project. We evaluated a native GPU re-implementation of the covariance-model
dynamic-programming kernels and **decided not to pursue it**: v1's algorithmic CPU
speedups already cover typical genome-scale use, and the remaining hot path is
Infernal's already-tuned CM search. That exploration is parked unless a concrete
need arises for throughput well beyond what v1 + Infernal deliver on CPU (e.g.
routine scanning of many multi-gigabase genomes).

## License

trnascan-py is licensed **GPL-3.0-or-later** (see [LICENSE](LICENSE)). It derives
some logic from the reference tRNAscan-SE 2.0 (GPL-3) and depends at runtime on
Infernal (GPL-3) and tRNAscan-SE 2.0, which you install separately; their
covariance models are used at runtime and are not bundled or redistributed here.

## Acknowledgements

This project stands on [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/)
(Todd Lowe, Patricia Chan, and Sean Eddy) and
[Infernal](http://eddylab.org/infernal/) (Eddy–Rivas labs). Please cite their
papers if you use trnascan-py in research:

- Chan, Lin, Chen, Huang, Hann, Lowe. *tRNAscan-SE 2.0.* (2021).
- Nawrocki & Eddy. *Infernal 1.1: 100-fold faster RNA homology searches.*
  Bioinformatics (2013).
