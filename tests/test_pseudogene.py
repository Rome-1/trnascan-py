"""Unit tests for the faithful pseudogene decision (no external tools)."""

from __future__ import annotations

from trnascan_py.pseudogene import (
    MIN_HMM_SCORE,
    MIN_SS_SCORE,
    PSEUDO_FILTER_SCORE,
    is_pseudogene,
    secondary_structure_score,
)


def test_secondary_structure_score() -> None:
    assert secondary_structure_score(78.4, 49.6) == 78.4 - 49.6


def test_high_total_never_pseudo() -> None:
    # At/above the filter score a hit is never pseudo, even with a low ss/hmm split.
    assert is_pseudogene(PSEUDO_FILTER_SCORE, 0.0) is False
    assert is_pseudogene(80.0, 79.0) is False  # ss=1 < MIN_SS_SCORE but score high


def test_low_total_intact_structure_not_pseudo() -> None:
    # total 50 (< 55), hmm 30 (>= 10), ss = 20 (>= 5) -> real tRNA
    assert is_pseudogene(50.0, 30.0) is False


def test_low_secondary_structure_is_pseudo() -> None:
    # total 50, hmm 48 -> ss = 2 < MIN_SS_SCORE -> pseudo
    assert is_pseudogene(50.0, 48.0) is True


def test_low_primary_score_is_pseudo() -> None:
    # total 50, hmm 8 (< MIN_HMM_SCORE) -> pseudo (ss=42 is fine, but hmm too low)
    assert is_pseudogene(50.0, 8.0) is True


def test_thresholds_are_reference_defaults() -> None:
    assert (PSEUDO_FILTER_SCORE, MIN_SS_SCORE, MIN_HMM_SCORE) == (55.0, 5.0, 10.0)
