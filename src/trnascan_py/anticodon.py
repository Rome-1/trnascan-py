"""Anticodon extraction from Infernal covariance-model alignments.

The tRNA cloverleaf has a fixed consensus secondary structure. In an Infernal
``cmsearch -A`` Stockholm alignment, the ``#=GC SS_cons`` line gives that
structure in WUSS notation. The anticodon loop is the **second hairpin loop**
(the D-arm loop is first, the anticodon-arm loop second), and the anticodon is
the middle three nucleotides of that 7-nt loop (standard positions 34-36).

We read those three residues from the aligned hit sequence, restricting to
*consensus* columns (``#=GC RF`` != ``.``) so that intron insertions — which
appear as non-consensus insert columns inside the anticodon loop — do not shift
the reading frame. Residues are returned in DNA convention (U→T, uppercase),
matching reference tRNAscan-SE output.
"""

from __future__ import annotations

from dataclasses import dataclass

# WUSS open/close bracket pairs used in Infernal SS_cons.
_OPEN = "<([{"
_CLOSE = ">)]}"
_PAIR = {"<": ">", "(": ")", "[": "]", "{": "}"}


@dataclass(frozen=True, slots=True)
class StockholmHit:
    """One aligned sequence from a Stockholm block, with shared structure lines."""

    name: str
    """Sequence name as written by Infernal, e.g. ``"contig/121-193"``."""

    aligned_seq: str
    ss_cons: str
    rf: str


def _pair_map(ss_cons: str) -> dict[int, int]:
    """Map each paired column index to its partner using a per-family stack."""
    stacks: dict[str, list[int]] = {c: [] for c in _OPEN}
    partner: dict[int, int] = {}
    # Map a closing char back to its opening char.
    close_to_open = {v: k for k, v in _PAIR.items()}
    for i, ch in enumerate(ss_cons):
        if ch in _OPEN:
            stacks[ch].append(i)
        elif ch in _CLOSE:
            opener = close_to_open[ch]
            if stacks[opener]:
                j = stacks[opener].pop()
                partner[i] = j
                partner[j] = i
    return partner


def _hairpin_loops(ss_cons: str) -> list[tuple[int, int]]:
    """Return ``(start, end)`` column spans of hairpin loops, ordered 5'→3'.

    A hairpin loop is a maximal run of unpaired columns whose flanking columns
    are base-paired *to each other* (the innermost pair of the enclosing stem).
    """
    partner = _pair_map(ss_cons)
    loops: list[tuple[int, int]] = []
    n = len(ss_cons)
    i = 0
    while i < n:
        if i in partner:
            i += 1
            continue
        j = i
        while j + 1 < n and (j + 1) not in partner:
            j += 1
        left, right = i - 1, j + 1
        if left >= 0 and right < n and partner.get(left) == right:
            loops.append((i, j))
        i = j + 1
    return loops


def anticodon_from_alignment(hit: StockholmHit) -> str | None:
    """Extract the anticodon (DNA letters) from a single aligned tRNA hit.

    Returns ``None`` if the structure does not have at least two hairpin loops
    or the anticodon loop yields fewer than three consensus residues.
    """
    loops = _hairpin_loops(hit.ss_cons)
    if len(loops) < 2:
        return None
    start, end = loops[1]  # second hairpin = anticodon loop

    # Keep only consensus columns (RF != '.'/'-'); these are the canonical
    # loop positions, excluding any intron insert columns.
    residues: list[str] = []
    for col in range(start, end + 1):
        if col >= len(hit.rf) or col >= len(hit.aligned_seq):
            break
        if hit.rf[col] in ".-":
            continue
        base = hit.aligned_seq[col]
        if base in ".-":
            continue
        residues.append(base.upper())

    if len(residues) < 3:
        return None
    mid = len(residues) // 2
    triplet = residues[mid - 1 : mid + 2]
    return "".join(triplet).replace("U", "T")


def parse_stockholm(text: str) -> list[StockholmHit]:
    """Parse a (possibly interleaved) Stockholm alignment from ``cmsearch -A``.

    Concatenates per-sequence rows and the ``#=GC SS_cons`` / ``#=GC RF`` lines
    across wrapped blocks. Returns one :class:`StockholmHit` per aligned sequence.
    """
    seqs: dict[str, list[str]] = {}
    order: list[str] = []
    ss_parts: list[str] = []
    rf_parts: list[str] = []

    for raw in text.splitlines():
        line = raw.rstrip("\n")
        if not line or line == "//":
            continue
        if line.startswith("#=GC SS_cons"):
            ss_parts.append(line.split(None, 2)[2])
            continue
        if line.startswith("#=GC RF"):
            rf_parts.append(line.split(None, 2)[2])
            continue
        if line.startswith("#"):
            continue  # other annotation (#=GF, #=GS, #=GR)
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        name, seq = parts
        if name not in seqs:
            seqs[name] = []
            order.append(name)
        seqs[name].append(seq)

    ss_cons = "".join(ss_parts)
    rf = "".join(rf_parts)
    return [
        StockholmHit(name=name, aligned_seq="".join(seqs[name]), ss_cons=ss_cons, rf=rf)
        for name in order
    ]


__all__ = ["StockholmHit", "anticodon_from_alignment", "parse_stockholm"]
