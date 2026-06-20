"""Minimal FASTA reading and subsequence extraction.

Used to pull confirmed candidate loci out of a genome so the isotype ``cmscan``
runs on a handful of short sequences instead of the whole genome.

This loads sequences into memory; fine for the bacterial/small-eukaryote genomes
v1 targets. Very large genomes (multi-Gb) would want indexed random access
(e.g. ``esl-sfetch``); tracked for later.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path


def read_fasta(path: str | Path) -> dict[str, str]:
    """Read a FASTA file into a ``{seq_id: sequence}`` dict.

    The sequence id is the first whitespace-delimited token of the header (the
    same convention Infernal uses for its ``target name`` column).
    """
    seqs: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    with Path(path).open() as fh:
        for line in fh:
            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(chunks)
                current = line[1:].split(None, 1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if current is not None:
        seqs[current] = "".join(chunks)
    return seqs


_COMPLEMENT = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement (DNA), preserving case and passing through N."""
    return seq.translate(_COMPLEMENT)[::-1]


def extract(seq: str, start: int, end: int, pad: int = 0) -> str:
    """Extract the forward-strand subsequence spanning 1-based inclusive coords.

    ``start``/``end`` may be in either order (minus-strand hits report
    ``start > end``); the span ``[min, max]`` is returned in forward orientation
    plus ``pad`` bases on each side, clamped to the sequence ends. Strand is
    irrelevant here because downstream ``cmscan`` searches both strands.
    """
    lo = max(0, min(start, end) - 1 - pad)
    hi = min(len(seq), max(start, end) + pad)
    return seq[lo:hi]


def write_fasta(
    records: Iterator[tuple[str, str]] | list[tuple[str, str]], path: str | Path
) -> None:
    """Write ``(name, sequence)`` records to a FASTA file (60-col wrapped)."""
    with Path(path).open("w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


__all__ = ["read_fasta", "extract", "write_fasta"]
