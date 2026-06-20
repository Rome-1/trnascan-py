"""Unit tests for the output parsers using captured fixtures (no external tools)."""

from __future__ import annotations

from trnascan_py.infernal import parse_tblout
from trnascan_py.models import Strand
from trnascan_py.oracle import parse_trnascan_out

# Captured `cmsearch --tblout` output (Infernal 1.1.5) for the test tRNA-Phe.
CMSEARCH_TBLOUT = """\
#target name         accession query name           accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
#------------------- --------- -------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- ---------------------
test_contig          -         euk-031616           -          cm        1       90      121      193      +    no    1 0.53   0.0   78.3   6.1e-22 !   partial yeast tRNA-Phe gene
"""

# Captured `cmscan --tblout` output (target = isotype model, query = sequence).
CMSCAN_TBLOUT = """\
#target name         accession query name           accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
euk-Phe              euk-Phe-r2 test_contig          -          cm        1       73      121      193      +    no    1 0.53   0.0  112.5   2.3e-35 !   -
"""

# Captured reference tRNAscan-SE tabular output.
TRNASCAN_OUT = """\
Sequence    \t\ttRNA\tBounds\ttRNA\tAnti\tIntron Bounds\tInf
Name        \ttRNA #\tBegin\tEnd\tType\tCodon\tBegin\tEnd\tScore\tNote
--------    \t------\t-----\t------\t----\t-----\t-----\t----\t------\t------
test_contig \t1\t121\t193\tPhe\tGAA\t0\t0\t78.4
"""


def test_parse_cmsearch_tblout() -> None:
    hits = parse_tblout(CMSEARCH_TBLOUT)
    assert len(hits) == 1
    h = hits[0]
    assert h.target_name == "test_contig"
    assert h.query_name == "euk-031616"
    assert h.seq_from == 121
    assert h.seq_to == 193
    assert h.strand is Strand.PLUS
    assert h.score == 78.3
    assert h.inc == "!"


def test_parse_cmscan_tblout_swaps_target_query() -> None:
    hits = parse_tblout(CMSCAN_TBLOUT)
    assert len(hits) == 1
    h = hits[0]
    # For cmscan, target is the model and query is the sequence.
    assert h.target_name == "euk-Phe"
    assert h.query_name == "test_contig"
    assert h.score == 112.5


def test_parse_cmsearch_skips_comments_and_blank() -> None:
    assert parse_tblout("#only a comment\n\n") == []


def test_parse_trnascan_out() -> None:
    hits = parse_trnascan_out(TRNASCAN_OUT)
    assert len(hits) == 1
    h = hits[0]
    assert h.seq_id == "test_contig"
    assert h.start == 121
    assert h.end == 193
    assert h.strand is Strand.PLUS
    assert h.isotype == "Phe"
    assert h.anticodon == "GAA"
    assert h.score == 78.4
    assert h.intron_start is None
    assert h.intron_end is None


def test_parse_trnascan_minus_strand_from_coords() -> None:
    row = "ctg\t1\t500\t430\tLys\tCTT\t0\t0\t60.2\n"
    hits = parse_trnascan_out(row)
    assert len(hits) == 1
    assert hits[0].strand is Strand.MINUS


def test_parse_trnascan_with_intron() -> None:
    row = "ctg\t1\t100\t200\tTyr\tGTA\t140\t160\t55.0\n"
    hits = parse_trnascan_out(row)
    assert hits[0].intron_start == 140
    assert hits[0].intron_end == 160
