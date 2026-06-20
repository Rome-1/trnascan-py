"""Unit tests for the streaming FASTA index (no external tools)."""

from __future__ import annotations

from pathlib import Path

import pytest

from trnascan_py.faidx import FastaIndex, UnindexableFastaError
from trnascan_py.fasta import read_fasta


def _write(p: Path, text: str) -> Path:
    p.write_text(text)
    return p


def test_fetch_matches_read_fasta(tmp_path: Path) -> None:
    # 10 bases per line, two records, last lines short.
    fa = _write(
        tmp_path / "g.fa",
        ">chr1 desc\nACGTACGTAC\nGGGGCCCCAA\nTTT\n>chr2\nAAAACCCCGG\nTT\n",
    )
    idx = FastaIndex.build(fa)
    mem = read_fasta(fa)

    for seqid in ("chr1", "chr2"):
        full = mem[seqid]
        # fetch the whole sequence
        assert idx.fetch(seqid, 1, len(full)) == full
        # fetch every sub-span and compare to the in-memory slice
        for a in range(1, len(full) + 1):
            for b in range(a, len(full) + 1):
                assert idx.fetch(seqid, a, b) == full[a - 1 : b], (seqid, a, b)


def test_fetch_spans_line_boundary(tmp_path: Path) -> None:
    fa = _write(tmp_path / "g.fa", ">s\nAAAAACCCCC\nGGGGGTTTTT\n")
    idx = FastaIndex.build(fa)
    # bases 8..13 cross the line break (CCC | GGG)
    assert idx.fetch("s", 8, 13) == "CCCGGG"


def test_fetch_coords_either_order_and_clamp(tmp_path: Path) -> None:
    fa = _write(tmp_path / "g.fa", ">s\nACGTACGTAC\n")
    idx = FastaIndex.build(fa)
    assert idx.fetch("s", 6, 3) == idx.fetch("s", 3, 6) == "GTAC"
    assert idx.fetch("s", 1, 999) == "ACGTACGTAC"  # clamps to end


def test_crlf_line_endings(tmp_path: Path) -> None:
    fa = _write(tmp_path / "g.fa", ">s\r\nACGTAC\r\nGTACGT\r\n")
    idx = FastaIndex.build(fa)
    assert idx.fetch("s", 1, 12) == "ACGTACGTACGT"
    assert idx.fetch("s", 5, 8) == "ACGT"


def test_unindexable_nonuniform_width(tmp_path: Path) -> None:
    # A short non-final line followed by more sequence is not indexable.
    fa = _write(tmp_path / "g.fa", ">s\nACG\nACGTAC\n")
    with pytest.raises(UnindexableFastaError):
        FastaIndex.build(fa)


def test_contains(tmp_path: Path) -> None:
    fa = _write(tmp_path / "g.fa", ">a\nACGT\n>b\nTTTT\n")
    idx = FastaIndex.build(fa)
    assert "a" in idx and "b" in idx and "c" not in idx
