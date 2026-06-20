"""Unit tests for anticodon -> isotype translation (no external tools)."""

from __future__ import annotations

from trnascan_py.genetic_code import AMBIGUOUS_ANTICODONS, isotype_from_anticodon


def test_standard_anticodons() -> None:
    assert isotype_from_anticodon("GAA") == "Phe"  # codon TTC
    assert isotype_from_anticodon("CTT") == "Lys"  # codon AAG
    assert isotype_from_anticodon("TGC") == "Ala"  # codon GCA
    assert isotype_from_anticodon("GCC") == "Gly"  # codon GGC


def test_rna_letters_accepted() -> None:
    assert isotype_from_anticodon("GAA") == isotype_from_anticodon("gaa")
    assert isotype_from_anticodon("UUC") == isotype_from_anticodon("TTC")


def test_special_cases() -> None:
    assert isotype_from_anticodon("TCA") == "SeC"  # recoded stop
    assert isotype_from_anticodon("CAT") == "Met"  # initiator/Ile2 need CM to refine


def test_suppressor_anticodons() -> None:
    assert isotype_from_anticodon("CTA") == "Sup"  # amber: decodes TAG
    assert isotype_from_anticodon("TTA") == "Sup"  # ochre: decodes TAA


def test_invalid() -> None:
    assert isotype_from_anticodon("NNN") == "Undet"
    assert isotype_from_anticodon("") == "Undet"
    assert isotype_from_anticodon("AC") == "Undet"


def test_ambiguous_set() -> None:
    assert "CAT" in AMBIGUOUS_ANTICODONS
    assert "TCA" in AMBIGUOUS_ANTICODONS
    assert "GAA" not in AMBIGUOUS_ANTICODONS
