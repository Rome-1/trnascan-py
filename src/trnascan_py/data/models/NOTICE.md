# Bundled covariance models — provenance and license

The `TRNAinf-*` covariance-model files in this directory are **redistributed
unmodified from [tRNAscan-SE 2.0](http://lowelab.ucsc.edu/tRNAscan-SE/)** (UCSC
Lowe Lab — Todd Lowe, Patricia Chan, Sean Eddy and contributors).

They are licensed under the **GNU General Public License v3**, the same license
as this package, and are redistributed here under those terms. They are included
so that `trnascan-py` can run with only the Infernal `cmsearch` binary installed,
without a full tRNAscan-SE installation.

- Upstream: http://lowelab.ucsc.edu/tRNAscan-SE/
- License: GPL-3 (see the repository `LICENSE`).
- These files are **not** part of trnascan-py's own authored code; copyright in
  them remains with their original authors.

Bundled subset (the models trnascan-py searches with):

- `TRNAinf-{euk,bact,arch}.cm` — general (isotype-merged) models, first pass.
- `TRNAinf-{euk,bact,arch}-iso` — per-isotype model sets (used via `cmsearch`).
- `TRNAinf-{euk,bact,arch}-ns.cm` — no-secondary-structure models (pseudogene scoring).
- `TRNAinf-mito-vert`, `TRNAinf-mito-mammal` — mitochondrial model sets.

The `cmpress` index files (`.i1f/.i1i/.i1m/.i1p`) are intentionally **not**
bundled: trnascan-py drives `cmsearch` (which reads the model files directly), not
`cmscan`, so the press indices are unnecessary.
