#!/usr/bin/env python3
"""Profile trnascan-py v1 vs reference tRNAscan-SE on a genome.

Times the end-to-end ``scan()`` (and the first-pass ``cmsearch`` stage that
dominates the optimized pipeline), compares against reference tRNAscan-SE, and
prints a genome-scale extrapolation. Feeds the v2 GPU-acceleration decision.

Usage:
    python scripts/profile_v1.py GENOME.fa --domain bact [--cpu N]
"""

from __future__ import annotations

import argparse
import tempfile
import time
from pathlib import Path

from trnascan_py.infernal import run_cmsearch
from trnascan_py.models_registry import resolve_domain
from trnascan_py.oracle import run_trnascan
from trnascan_py.pipeline import scan


def _genome_size(fasta: Path) -> int:
    total = 0
    with fasta.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def _timed(label: str, fn):  # type: ignore[no-untyped-def]
    t0 = time.perf_counter()
    result = fn()
    dt = time.perf_counter() - t0
    print(f"  {label:<46} {dt:8.2f}s", flush=True)
    return result, dt


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta", type=Path)
    ap.add_argument("-d", "--domain", default="bact", choices=("euk", "bact", "arch"))
    ap.add_argument("--cpu", type=int, default=None)
    ap.add_argument("--method", default="hybrid", choices=("hybrid", "anticodon", "cm"))
    ap.add_argument("--skip-reference", action="store_true")
    args = ap.parse_args()

    size = _genome_size(args.fasta)
    print(f"Genome: {args.fasta}  ({size:,} bp)  domain={args.domain}  cpu={args.cpu}")
    general_cm, _ = resolve_domain(args.domain)

    print("\n[trnascan-py v1]")

    def _first_pass():  # type: ignore[no-untyped-def]
        with tempfile.TemporaryDirectory() as tmp:
            aln = Path(tmp) / "aln.sto"
            return run_cmsearch(general_cm, args.fasta, cpu=args.cpu, alignment_out=aln)

    _, t_search = _timed("first-pass cmsearch + align (general CM)", _first_pass)
    hits, t_scan = _timed(
        f"scan() total (isotype_method={args.method})",
        lambda: scan(args.fasta, domain=args.domain, isotype_method=args.method, cpu=args.cpu),
    )
    print(f"  -> {len(hits)} tRNAs", flush=True)

    t_ref = None
    if not args.skip_reference:
        print("\n[reference tRNAscan-SE 2.0]")
        ref, t_ref = _timed(
            "tRNAscan-SE total", lambda: run_trnascan(args.fasta, domain=args.domain)
        )
        print(f"  -> {len(ref)} tRNAs", flush=True)

    print("\n[summary]")
    print(f"  first-pass cmsearch share of scan(): {100 * t_search / t_scan:.0f}%  "
          "(the v2 GPU target)")
    print(f"  v1 throughput: {size / t_scan / 1000:,.1f} kbp/s  ({t_scan:.1f}s for {size:,} bp)")
    if t_ref:
        print(f"  reference throughput: {size / t_ref / 1000:,.1f} kbp/s")
        print(f"  v1 vs reference: {t_ref / t_scan:.2f}x  (>1 = v1 faster)")
    for gb, label in [(3.1, "human 3.1 Gb"), (1.0, "1 Gb"), (0.1, "100 Mb")]:
        scaled = t_scan * (gb * 1e9) / size
        mins, hrs = scaled / 60, scaled / 3600
        print(f"  extrapolated scan() for {label:<14}: {mins:8.1f} min  ({hrs:.2f} h)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
