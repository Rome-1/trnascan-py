"""Unit tests for the minimal FASTA reader / extractor (no external tools)."""

from __future__ import annotations

from pathlib import Path

from trnascan_py.fasta import extract, read_fasta, reverse_complement, write_fasta


def test_reverse_complement() -> None:
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAGG") == "CCTT"
    assert reverse_complement("acgtN") == "Nacgt"


def test_read_fasta_multi(tmp_path: Path) -> None:
    p = tmp_path / "g.fa"
    p.write_text(">chr1 desc here\nACGT\nACGT\n>chr2\nTTTT\n")
    seqs = read_fasta(p)
    assert seqs == {"chr1": "ACGTACGT", "chr2": "TTTT"}


def test_extract_forward() -> None:
    seq = "AAAAACCCCCGGGGG"  # 1-based: 6..10 = CCCCC
    assert extract(seq, 6, 10) == "CCCCC"


def test_extract_minus_coords_same_span() -> None:
    seq = "AAAAACCCCCGGGGG"
    # start > end (minus strand) returns the same forward span.
    assert extract(seq, 10, 6) == "CCCCC"


def test_extract_pad_clamps_at_ends() -> None:
    seq = "AAAAACCCCCGGGGG"
    assert extract(seq, 1, 3, pad=100) == seq  # clamps to full sequence


def test_write_then_read_roundtrip(tmp_path: Path) -> None:
    p = tmp_path / "out.fa"
    write_fasta([("a", "ACGTACGTAC" * 8), ("b", "TT")], p)
    seqs = read_fasta(p)
    assert seqs["a"] == "ACGTACGTAC" * 8
    assert seqs["b"] == "TT"
