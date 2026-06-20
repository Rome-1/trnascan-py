"""Mitochondrial scan: unit (model parsing) + end-to-end + differential."""

from __future__ import annotations

from pathlib import Path

import pytest

from trnascan_py.pipeline import _parse_mito_model, scan


def test_parse_mito_model_degenerate() -> None:
    assert _parse_mito_model("LeuTAA") == ("Leu", "TAA")
    assert _parse_mito_model("LeuTAG") == ("Leu", "TAG")
    assert _parse_mito_model("SerGCT") == ("Ser", "GCT")
    assert _parse_mito_model("SerTGA") == ("Ser", "TGA")


def test_parse_mito_model_plain() -> None:
    assert _parse_mito_model("Phe") == ("Phe", "")
    assert _parse_mito_model("Met") == ("Met", "")
    assert _parse_mito_model("Cys-no-darm") == ("Cys", "")


@pytest.mark.requires_infernal
def test_mito_scan_finds_22_human_trnas(mtdna_fasta: Path) -> None:
    hits = scan(mtdna_fasta, domain="mito-vert")
    # Human mtDNA encodes 22 tRNAs.
    assert len(hits) == 22
    isos = {h.isotype for h in hits}
    # Includes the structurally divergent Tyr (recovered by glocal re-scoring).
    assert "Tyr" in isos
    assert {"Leu", "Ser"} <= isos
    # Leu/Ser anticodons come from the model name.
    leu = [h for h in hits if h.isotype == "Leu"]
    assert all(h.anticodon in ("TAA", "TAG") for h in leu)
