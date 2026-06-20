"""Streaming FASTA index for random-access subsequence extraction.

For multi-Gb assemblies we must not hold the genome in memory. This builds a
``samtools faidx``-style index in a single streaming pass (O(1) memory beyond the
per-sequence metadata) and then extracts candidate loci by seeking, rather than
reading the whole genome into a dict.

The index requires the common FASTA convention of a uniform line width within
each record (the last line of a record may be shorter). Files that violate this
raise :class:`UnindexableFastaError`; callers may fall back to
:func:`trnascan_py.fasta.read_fasta`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


class UnindexableFastaError(ValueError):
    """Raised when a FASTA record does not have a uniform line width."""


@dataclass(frozen=True, slots=True)
class FaidxEntry:
    """Index record for one sequence (byte geometry, like a ``.fai`` line)."""

    length: int
    """Number of bases in the sequence."""

    offset: int
    """Byte offset of the first base."""

    linebases: int
    """Bases per full line."""

    linewidth: int
    """Bytes per full line, including the line terminator."""


class FastaIndex:
    """A built index over a FASTA file, supporting ``fetch`` by coordinate."""

    def __init__(self, path: str | Path, entries: dict[str, FaidxEntry]) -> None:
        self.path = Path(path)
        self.entries = entries

    def __contains__(self, seqid: str) -> bool:
        return seqid in self.entries

    @classmethod
    def build(cls, path: str | Path) -> FastaIndex:
        """Build an index by streaming the file once (low memory)."""
        path = Path(path)
        entries: dict[str, FaidxEntry] = {}
        name: str | None = None
        offset = 0
        length = 0
        linebases = 0
        linewidth = 0
        short_line_seen = False  # a non-final line shorter than the first is illegal
        pos = 0

        def finalize() -> None:
            if name is not None:
                entries[name] = FaidxEntry(length, offset, linebases, linewidth)

        with path.open("rb") as fh:
            for raw in fh:
                if raw.startswith(b">"):
                    finalize()
                    name = raw[1:].split()[0].decode() if raw[1:].strip() else ""
                    pos += len(raw)
                    offset = pos
                    length = 0
                    linebases = 0
                    linewidth = 0
                    short_line_seen = False
                    continue
                stripped = raw.rstrip(b"\r\n")
                nbases = len(stripped)
                if linebases == 0:
                    linebases = nbases
                    linewidth = len(raw)
                else:
                    if short_line_seen:
                        # A previous line was short but the record continued —
                        # line width is not uniform, so byte math would be wrong.
                        raise UnindexableFastaError(
                            f"non-uniform line width in record {name!r}"
                        )
                    if nbases < linebases:
                        short_line_seen = True
                    elif nbases > linebases:
                        raise UnindexableFastaError(
                            f"non-uniform line width in record {name!r}"
                        )
                length += nbases
                pos += len(raw)
        finalize()
        return cls(path, entries)

    def fetch(self, seqid: str, start: int, end: int) -> str:
        """Return the forward-strand subsequence for 1-based inclusive coords.

        ``start``/``end`` may be given in either order; the span ``[min, max]``
        is returned (clamped to the sequence ends). Raises ``KeyError`` if the
        sequence is unknown.
        """
        e = self.entries[seqid]
        lo = max(1, min(start, end))
        hi = min(e.length, max(start, end))
        if hi < lo or e.linebases == 0:
            return ""

        def byte_of(base1: int) -> int:
            z = base1 - 1
            return e.offset + (z // e.linebases) * e.linewidth + (z % e.linebases)

        read_start = byte_of(lo)
        read_end = byte_of(hi) + 1  # include the last base's byte
        with self.path.open("rb") as fh:
            fh.seek(read_start)
            raw = fh.read(read_end - read_start)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode()


__all__ = ["FastaIndex", "FaidxEntry", "UnindexableFastaError"]
