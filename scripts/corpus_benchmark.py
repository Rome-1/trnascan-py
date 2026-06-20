#!/usr/bin/env python3
"""Differential + speed benchmark across a representative genome corpus.

For each (genome, domain) it runs trnascan-py ``scan`` and reference tRNAscan-SE,
then reports detection agreement, coordinate/isotype/anticodon concordance, and
the speed ratio. Used to validate that faithfulness and the ~1.7x speedup hold
across domains and genome sizes.

Usage:
    python scripts/corpus_benchmark.py            # default corpus under /tmp/corpus
    python scripts/corpus_benchmark.py manifest.tsv   # custom: path<TAB>domain<TAB>label
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

from trnascan_py.models import TRNAHit
from trnascan_py.oracle import run_trnascan
from trnascan_py.pipeline import scan

COORD_TOL = 5

# (path, our-domain, reference-domain, label)
DEFAULT_CORPUS = [
    ("/tmp/corpus/human_mito.fa", "mito-vert", "mito-vert", "vert-mito (16.5 kb)"),
    ("/tmp/corpus/arabidopsis_chloro.fa", "bact", "bact", "plant-chloroplast (154 kb)"),
    ("/tmp/corpus/yeast_chrI_euk.fa", "euk", "euk", "euk-nuclear yeast chrI (230 kb)"),
    ("/tmp/corpus/mjannaschii_arch.fa", "arch", "arch", "archaea M.jann (1.66 Mb)"),
    ("/tmp/corpus/ecoli_bact.fa", "bact", "bact", "bacteria E.coli (4.64 Mb)"),
]


def _match(o: TRNAHit, r: TRNAHit) -> bool:
    return (
        o.strand is r.strand
        and abs(min(o.start, o.end) - min(r.start, r.end)) <= COORD_TOL
        and abs(max(o.start, o.end) - max(r.start, r.end)) <= COORD_TOL
    )


def _bp(path: Path) -> int:
    return sum(len(line.strip()) for line in path.open() if not line.startswith(">"))


def bench(path: str, our_domain: str, ref_domain: str, label: str) -> dict:
    p = Path(path)
    size = _bp(p)
    t0 = time.perf_counter()
    ours = scan(p, domain=our_domain)
    t_ours = time.perf_counter() - t0
    t0 = time.perf_counter()
    ref = run_trnascan(p, domain=ref_domain)
    t_ref = time.perf_counter() - t0

    matched = iso_ok = ac_cmp = ac_ok = 0
    for r in ref:
        m = next((o for o in ours if _match(o, r)), None)
        if m is None:
            continue
        matched += 1
        iso_ok += int(m.isotype == r.isotype)
        if m.anticodon and r.anticodon:
            ac_cmp += 1
            ac_ok += int(m.anticodon == r.anticodon)
    return {
        "label": label,
        "size": size,
        "ours_n": len(ours),
        "ref_n": len(ref),
        "matched": matched,
        "iso": f"{iso_ok}/{matched}",
        "ac": f"{ac_ok}/{ac_cmp}" if ac_cmp else "n/a",
        "t_ours": t_ours,
        "t_ref": t_ref,
        "ratio": (t_ref / t_ours) if t_ours else float("nan"),
    }


def main() -> int:
    if len(sys.argv) > 1:
        corpus = []
        for line in Path(sys.argv[1]).read_text().splitlines():
            if line.strip() and not line.startswith("#"):
                pth, dom, label = line.split("\t")
                corpus.append((pth, dom, dom, label))
    else:
        corpus = DEFAULT_CORPUS

    rows = []
    for path, od, rd, label in corpus:
        if not Path(path).exists():
            print(f"SKIP {label}: {path} missing", flush=True)
            continue
        print(f"running {label} ...", flush=True)
        rows.append(bench(path, od, rd, label))

    print("\n" + "=" * 100)
    hdr = f"{'genome':<32}{'bp':>10} {'ours':>5}{'ref':>5}{'match':>6} {'isotype':>9}{'codon':>8}"
    hdr += f"{'ours_s':>9}{'ref_s':>8}{'ratio':>7}"
    print(hdr)
    print("-" * 100)
    for r in rows:
        print(
            f"{r['label']:<32}{r['size']:>10,} {r['ours_n']:>5}{r['ref_n']:>5}{r['matched']:>6} "
            f"{r['iso']:>9}{r['ac']:>8}{r['t_ours']:>9.2f}{r['t_ref']:>8.2f}{r['ratio']:>6.2f}x",
            flush=True,
        )
    ratios = [r["ratio"] for r in rows if r["ratio"] == r["ratio"]]
    if ratios:
        print("-" * 100)
        print(f"speed ratio: min {min(ratios):.2f}x  median {sorted(ratios)[len(ratios)//2]:.2f}x  "
              f"max {max(ratios):.2f}x  (>1 = trnascan-py faster)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
