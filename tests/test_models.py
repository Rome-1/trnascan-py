"""Unit tests for the typed result records (no external tools)."""

from __future__ import annotations

import pytest

from trnascan_py.models import Strand, TRNAHit


def test_strand_from_str() -> None:
    assert Strand.from_str("+") is Strand.PLUS
    assert Strand.from_str("-") is Strand.MINUS
    assert Strand.from_str("plus") is Strand.PLUS
    with pytest.raises(ValueError):
        Strand.from_str("?")


def test_length_plus_strand() -> None:
    hit = TRNAHit("chr1", 121, 193, Strand.PLUS, "Phe", "GAA", 78.3)
    assert hit.length == 73


def test_length_minus_strand() -> None:
    hit = TRNAHit("chr1", 193, 121, Strand.MINUS, "Phe", "GAA", 78.3)
    assert hit.length == 73


def test_overlaps_same_strand() -> None:
    a = TRNAHit("chr1", 100, 200, Strand.PLUS, "Phe", "GAA", 50.0)
    b = TRNAHit("chr1", 150, 250, Strand.PLUS, "Lys", "CTT", 40.0)
    assert a.overlaps(b)


def test_no_overlap_different_strand() -> None:
    a = TRNAHit("chr1", 100, 200, Strand.PLUS, "Phe", "GAA", 50.0)
    b = TRNAHit("chr1", 150, 250, Strand.MINUS, "Lys", "CTT", 40.0)
    assert not a.overlaps(b)


def test_no_overlap_different_seq() -> None:
    a = TRNAHit("chr1", 100, 200, Strand.PLUS, "Phe", "GAA", 50.0)
    b = TRNAHit("chr2", 150, 250, Strand.PLUS, "Lys", "CTT", 40.0)
    assert not a.overlaps(b)
