"""End-to-end pipeline + CLI smoke tests (require Infernal + bundled models)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from trnascan_py.cli import main
from trnascan_py.models import Strand
from trnascan_py.pipeline import scan


@pytest.mark.requires_infernal
def test_scan_finds_phe(yeast_phe_fasta: Path) -> None:
    hits = scan(yeast_phe_fasta, domain="euk")
    assert len(hits) == 1
    h = hits[0]
    assert h.strand is Strand.PLUS
    assert h.isotype == "Phe"
    assert h.anticodon == "GAA"
    assert h.score > 20.0


@pytest.mark.requires_infernal
def test_scan_no_isotype_pass(yeast_phe_fasta: Path) -> None:
    hits = scan(yeast_phe_fasta, domain="euk", classify_isotype=False)
    assert len(hits) == 1
    assert hits[0].isotype == "Undet"


@pytest.mark.requires_infernal
def test_scan_anticodon_method(yeast_phe_fasta: Path) -> None:
    # GAA is unambiguous, so the genetic-code path alone yields Phe (no CM scan).
    hits = scan(yeast_phe_fasta, domain="euk", isotype_method="anticodon")
    assert len(hits) == 1
    assert hits[0].isotype == "Phe"
    assert hits[0].anticodon == "GAA"


@pytest.mark.requires_infernal
def test_real_trna_not_flagged_pseudo(yeast_phe_fasta: Path) -> None:
    # A real, high-scoring tRNA must never be flagged pseudo.
    hits = scan(yeast_phe_fasta, domain="euk")
    assert len(hits) == 1
    assert hits[0].note == ""


@pytest.mark.requires_infernal
def test_scan_detects_intron(intron_trna_fasta: Path) -> None:
    hits = scan(intron_trna_fasta, domain="euk")
    assert len(hits) == 1
    h = hits[0]
    assert h.isotype == "Phe"
    assert h.intron_start is not None and h.intron_end is not None
    # 16-nt synthetic intron
    assert h.intron_end - h.intron_start + 1 == 16


@pytest.mark.requires_infernal
def test_cli_scan_json(yeast_phe_fasta: Path, capsys: pytest.CaptureFixture[str]) -> None:
    rc = main(["scan", str(yeast_phe_fasta), "--domain", "euk", "--json"])
    assert rc == 0
    out = json.loads(capsys.readouterr().out)
    assert len(out) == 1
    assert out[0]["isotype"] == "Phe"
    assert out[0]["strand"] == "+"
    assert "note" in out[0]  # pseudogene flag field is exposed


@pytest.mark.requires_infernal
def test_cli_scan_missing_file_errors_cleanly(capsys: pytest.CaptureFixture[str]) -> None:
    # A missing input must exit nonzero with a concise stderr message, not a traceback.
    rc = main(["scan", "/no/such/file.fa", "--domain", "euk"])
    assert rc == 1
    err = capsys.readouterr().err
    assert "error" in err.lower()
    assert "Traceback" not in err


def test_cli_version(capsys: pytest.CaptureFixture[str]) -> None:
    rc = main(["version"])
    assert rc == 0
    assert "trnascan-py" in capsys.readouterr().out
