"""Unit tests for intron extraction from alignments (no external tools)."""

from __future__ import annotations

from trnascan_py.anticodon import StockholmHit
from trnascan_py.intron import intron_from_alignment


def _hit(name: str, aligned: str, rf: str) -> StockholmHit:
    # ss_cons is unused by intron detection; pad it to the same length.
    return StockholmHit(name=name, aligned_seq=aligned, ss_cons="." * len(aligned), rf=rf)


# cols 0-3 consensus, 4-9 insert (6 nt), 10-13 consensus
ALIGNED = "ACGTaattttGCGT"
RF = "====......===="


def test_intron_plus_strand() -> None:
    # from=100 (+ strand): inserted residues at offsets 4..9 -> genomic 104..109
    assert intron_from_alignment(_hit("chr1/100-180", ALIGNED, RF)) == (104, 109)


def test_intron_minus_strand() -> None:
    # from=180 (- strand, from>to): offsets map downward -> 176..171
    assert intron_from_alignment(_hit("chr1/180-100", ALIGNED, RF)) == (176, 171)


def test_no_intron_when_all_consensus() -> None:
    assert intron_from_alignment(_hit("c/1-8", "ACGTACGT", "========")) is None


def test_short_insert_below_threshold_ignored() -> None:
    # only a 2-nt insert -> below MIN_INTRON_LEN
    assert intron_from_alignment(_hit("c/1-10", "ACGTaaGCGT", "====..====")) is None


def test_gaps_do_not_consume_genomic_position() -> None:
    # A deletion ('-') at a consensus column consumes no genomic base, so it must
    # not shift the genomic mapping of the downstream insert run.
    aligned = "AC-GTaaaaaaGT"  # the '-' at col 2 is a deletion (no residue)
    rf = "====.......=="  # cols 4-10 are insert columns (7 nt)
    res = intron_from_alignment(_hit("c/50-90", aligned, rf))
    assert res is not None
    start, end = res
    assert end - start + 1 == 7  # 7-nt insert run
    # residues before the run: A,C,G (3, the gap consumed none) -> offset 3
    assert start == 50 + 3
