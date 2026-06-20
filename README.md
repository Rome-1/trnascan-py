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

- Python ≥ 3.10
- The Infernal **`cmsearch`** binary on `PATH`
  ([Infernal](http://eddylab.org/infernal/) 1.1.x) — the search engine.
- The tRNA **covariance models**. These are **bundled with the package**, so a
  plain `pip install` is self-contained; you only need the `cmsearch` binary.
  (A system [tRNAscan-SE 2.0](http://lowelab.ucsc.edu/tRNAscan-SE/) install, if
  present, takes precedence and is also what the differential test suite uses as
  its oracle.) See [License](#license) for the bundled models' provenance.

To override model discovery, set `TRNASCAN_MODELS_DIR` to a directory containing
`TRNAinf-*.cm`.

## Install

**conda (one command — pulls the engine + models via tRNAscan-SE):**

```bash
conda install -c bioconda trnascan-py
```

**pip (models bundled; install the `cmsearch` binary separately):**

```bash
pip install git+https://github.com/Rome-1/trnascan-py.git   # latest; or: pip install .
# then ensure cmsearch is available, e.g.:
sudo apt-get install -y infernal        # or: conda install -c bioconda infernal
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

Copyright © 2026 Rome Thorstenson. trnascan-py is licensed **GPL-3.0-or-later**
(see [LICENSE](LICENSE)). It derives some logic from the reference tRNAscan-SE 2.0
(GPL-3) and drives the Infernal `cmsearch` binary (GPL-3), which you install
separately. The tRNA **covariance models bundled in this package are
redistributed unmodified from tRNAscan-SE 2.0 under GPL-3** — see
[`src/trnascan_py/data/models/NOTICE.md`](src/trnascan_py/data/models/NOTICE.md).
The Infernal / tRNAscan-SE compiled executables are **not** bundled.

## Citing

If you use trnascan-py in your research, please cite it **together with** the
tools it builds on (it drives them; it does not replace them). See
[`CITATION.cff`](CITATION.cff) or GitHub's "Cite this repository" button for the
full citation metadata.

- **trnascan-py** — Thorstenson, R. *trnascan-py: a fast Python interface for tRNA
  gene detection.* https://github.com/Rome-1/trnascan-py
- **tRNAscan-SE 2.0** — Chan, Lin, Chen, Huang, Hann, Lowe.
  *tRNAscan-SE 2.0: improved detection and functional classification of transfer
  RNA genes.* (2021).
- **Infernal 1.1** — Nawrocki & Eddy. *Infernal 1.1: 100-fold faster RNA homology
  searches.* Bioinformatics (2013).

## Acknowledgements

This project stands on [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/)
(Todd Lowe, Patricia Chan, and Sean Eddy) and
[Infernal](http://eddylab.org/infernal/) (Eddy–Rivas labs) — the engines and
covariance models that do the core science here.
