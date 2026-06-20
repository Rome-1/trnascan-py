"""Intron detection from Infernal covariance-model alignments.

A tRNA intron is extra sequence (most often in the anticodon loop) that is not
part of the conserved tRNA body. In a ``cmsearch -A`` Stockholm alignment it
shows up as a run of **insert** columns — columns the model does not have a
consensus state for (``#=GC RF`` is ``.``), carrying lowercase residues in the
aligned hit.

We find the largest such insert run (≥ :data:`MIN_INTRON_LEN`) and map its
residue offsets back to genomic coordinates using the hit's ``seqid/from-to``
name, reporting the intron in the same orientation as the tRNA bounds (so a
minus-strand intron has ``start > end``, matching reference tRNAscan-SE).
"""

from __future__ import annotations

from trnascan_py.anticodon import StockholmHit

# Minimum insert-run length (nt) to call an intron, to avoid spurious 1-2 nt
# insertions. Real tRNA introns are typically ~10-60 nt.
MIN_INTRON_LEN = 4


def _parse_coords(name: str) -> tuple[int, int] | None:
    """Parse ``"seqid/from-to"`` -> ``(from, to)`` (1-based, genomic)."""
    if "/" not in name:
        return None
    _, span = name.rsplit("/", 1)
    if "-" not in span:
        return None
    a, _, b = span.rpartition("-")
    try:
        return int(a), int(b)
    except ValueError:
        return None


def intron_from_alignment(hit: StockholmHit) -> tuple[int, int] | None:
    """Return the intron's genomic ``(start, end)`` or ``None`` if there is none.

    Coordinates are 1-based inclusive and oriented like the tRNA (minus-strand
    introns have ``start > end``).
    """
    coords = _parse_coords(hit.name)
    if coords is None:
        return None
    g_from, g_to = coords
    minus = g_from > g_to

    # Walk the aligned hit, tracking the genomic offset (count of consumed
    # residues), and collect maximal runs of inserted residues.
    runs: list[tuple[int, int]] = []
    cur: list[int] | None = None
    residue_idx = 0
    for col, base in enumerate(hit.aligned_seq):
        if base in ".-":  # gap in this sequence — consumes no genomic base
            continue
        is_insert = col < len(hit.rf) and hit.rf[col] in ".-"
        if is_insert:
            if cur is None:
                cur = [residue_idx, residue_idx]
            else:
                cur[1] = residue_idx
        elif cur is not None:
            runs.append((cur[0], cur[1]))
            cur = None
        residue_idx += 1
    if cur is not None:
        runs.append((cur[0], cur[1]))

    if not runs:
        return None
    a, b = max(runs, key=lambda r: r[1] - r[0])
    if b - a + 1 < MIN_INTRON_LEN:
        return None

    if minus:
        return g_from - a, g_from - b
    return g_from + a, g_from + b


__all__ = ["intron_from_alignment", "MIN_INTRON_LEN"]
