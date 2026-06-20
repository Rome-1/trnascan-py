"""Unit tests for anticodon extraction from Stockholm alignments (no tools)."""

from __future__ import annotations

from trnascan_py.anticodon import (
    _hairpin_loops,
    anticodon_from_alignment,
    parse_stockholm,
)

# Captured `cmsearch -A` Stockholm output for the test tRNA-Phe (Infernal 1.1.5).
STOCKHOLM = """\
# STOCKHOLM 1.0
#=GF AU Infernal 1.1.5

#=GS test_contig/121-193 DE partial yeast tRNA-Phe gene with flanking sequence

test_contig/121-193         GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAG------------------GuCCUGUGUUCGAUCCACAGAAUUCGCA
#=GR test_contig/121-193 PP *********************************************..................****************************
#=GC SS_cons                (((((((,,<<<<________>>>>,<<<<<_______>>>>>,,<<<<<<<____>>>>>>>,.,<<<<<_______>>>>>))))))):
#=GC RF                     =======++====++++++++====+=====++***++=====++=======++++=======+.+=====+++++++============+
//
"""


def test_parse_stockholm_single_hit() -> None:
    hits = parse_stockholm(STOCKHOLM)
    assert len(hits) == 1
    h = hits[0]
    assert h.name == "test_contig/121-193"
    assert len(h.aligned_seq) == len(h.ss_cons) == len(h.rf)


def test_hairpin_loops_finds_three() -> None:
    hits = parse_stockholm(STOCKHOLM)
    loops = _hairpin_loops(hits[0].ss_cons)
    # tRNA cloverleaf: D-loop, anticodon-loop, variable-loop, T-loop.
    assert len(loops) >= 2


def test_anticodon_is_gaa() -> None:
    hits = parse_stockholm(STOCKHOLM)
    assert anticodon_from_alignment(hits[0]) == "GAA"


def test_anticodon_none_without_structure() -> None:
    from trnascan_py.anticodon import StockholmHit

    flat = StockholmHit(name="x", aligned_seq="ACGUACGU", ss_cons="........", rf="========")
    assert anticodon_from_alignment(flat) is None


def test_parse_stockholm_interleaved() -> None:
    """Wrapped (multi-block) Stockholm rows concatenate per sequence."""
    text = (
        "# STOCKHOLM 1.0\n"
        "seqA          GCGG\n"
        "#=GC SS_cons  ((((\n"
        "#=GC RF       ====\n"
        "\n"
        "seqA          CCGC\n"
        "#=GC SS_cons  ))))\n"
        "#=GC RF       ====\n"
        "//\n"
    )
    hits = parse_stockholm(text)
    assert len(hits) == 1
    assert hits[0].aligned_seq == "GCGGCCGC"
    assert hits[0].ss_cons == "(((())))"
