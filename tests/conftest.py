"""Shared pytest fixtures and tool-availability gating."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent / "data"

_HAVE_INFERNAL = shutil.which("cmsearch") is not None and shutil.which("cmscan") is not None
_HAVE_TRNASCAN = shutil.which("tRNAscan-SE") is not None


@pytest.fixture
def data_dir() -> Path:
    return DATA_DIR


@pytest.fixture
def yeast_phe_fasta() -> Path:
    return DATA_DIR / "yeast_trna_phe.fa"


@pytest.fixture
def intron_trna_fasta() -> Path:
    return DATA_DIR / "intron_trna.fa"


@pytest.fixture
def mtdna_fasta() -> Path:
    return DATA_DIR / "human_mtdna.fa"


@pytest.fixture
def amber_suppressor_fasta() -> Path:
    return DATA_DIR / "amber_suppressor.fa"


def pytest_runtest_setup(item: pytest.Item) -> None:
    """Skip tool-dependent tests when the underlying binaries are absent."""
    if item.get_closest_marker("requires_infernal") and not _HAVE_INFERNAL:
        pytest.skip("Infernal (cmsearch/cmscan) not on PATH")
    if item.get_closest_marker("requires_trnascan") and not _HAVE_TRNASCAN:
        pytest.skip("reference tRNAscan-SE not on PATH")
