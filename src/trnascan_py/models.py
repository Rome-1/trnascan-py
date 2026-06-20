"""Typed result records for tRNA hits.

These mirror the columns of a tRNAscan-SE / Infernal output row closely enough
to support differential comparison against the reference tool, while staying a
clean Python API.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass


class Strand(enum.Enum):
    """Genomic strand a hit was found on."""

    PLUS = "+"
    MINUS = "-"

    @classmethod
    def from_str(cls, value: str) -> Strand:
        v = value.strip()
        if v in ("+", "plus", "PLUS"):
            return cls.PLUS
        if v in ("-", "minus", "MINUS"):
            return cls.MINUS
        raise ValueError(f"unrecognized strand: {value!r}")


@dataclass(frozen=True, slots=True)
class TRNAHit:
    """A single predicted tRNA gene.

    Coordinates are 1-based and inclusive, following tRNAscan-SE / Infernal
    conventions. ``start``/``end`` are reported in genomic orientation, so for a
    minus-strand hit ``start > end`` (as the reference tool reports it).
    """

    seq_id: str
    """Sequence (contig/chromosome) identifier the hit was found on."""

    start: int
    """1-based start coordinate (genomic orientation)."""

    end: int
    """1-based end coordinate (genomic orientation)."""

    strand: Strand
    """Strand the hit is on."""

    isotype: str
    """Predicted amino-acid isotype (e.g. ``"Lys"``, ``"Pseudo"``, ``"Undet"``)."""

    anticodon: str
    """Predicted anticodon (e.g. ``"CTT"``); empty string if undetermined."""

    score: float
    """Infernal covariance-model bit score."""

    intron_start: int | None = None
    """1-based intron start, or ``None`` if no intron was predicted."""

    intron_end: int | None = None
    """1-based intron end, or ``None`` if no intron was predicted."""

    note: str = ""
    """Free-text note (e.g. ``"pseudo"``), as the reference tool emits."""

    @property
    def length(self) -> int:
        """Span length in nucleotides (inclusive)."""
        return abs(self.end - self.start) + 1

    def overlaps(self, other: TRNAHit) -> bool:
        """Whether this hit overlaps ``other`` on the same sequence and strand."""
        if self.seq_id != other.seq_id or self.strand is not other.strand:
            return False
        a_lo, a_hi = sorted((self.start, self.end))
        b_lo, b_hi = sorted((other.start, other.end))
        return a_lo <= b_hi and b_lo <= a_hi
