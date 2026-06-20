"""Differential tests: trnascan-py vs reference tRNAscan-SE 2.0.

These require both Infernal and the reference tRNAscan-SE on PATH and are the
correctness oracle for v1. They assert that our Infernal-backed pipeline agrees
with the reference tool on detection, coordinates, strand, isotype, and score
(within tolerance) for small test genomes.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from trnascan_py.models import TRNAHit
from trnascan_py.oracle import run_trnascan
from trnascan_py.pipeline import scan

# tRNAscan-SE rescores with the isotype model, so its reported "Inf Score"
# differs slightly from the general-model bit score. Allow a small tolerance.
SCORE_TOL = 2.0
# Coordinates from the first-pass general-model search can differ from the
# reference's refined bounds by a few bp at the ends.
COORD_TOL = 5


def _match(ours: TRNAHit, ref: TRNAHit) -> bool:
    return (
        ours.seq_id == ref.seq_id
        and ours.strand is ref.strand
        and abs(min(ours.start, ours.end) - min(ref.start, ref.end)) <= COORD_TOL
        and abs(max(ours.start, ours.end) - max(ref.start, ref.end)) <= COORD_TOL
    )


@pytest.mark.differential
@pytest.mark.requires_infernal
@pytest.mark.requires_trnascan
def test_yeast_phe_matches_reference(yeast_phe_fasta: Path) -> None:
    ours = scan(yeast_phe_fasta, domain="euk")
    ref = run_trnascan(yeast_phe_fasta, domain="euk")

    # Same number of tRNA calls.
    assert len(ours) == len(ref) == 1, f"ours={ours!r} ref={ref!r}"

    o, r = ours[0], ref[0]
    assert _match(o, r), f"coordinate/strand mismatch: ours={o!r} ref={r!r}"
    assert o.isotype == r.isotype, f"isotype mismatch: {o.isotype} != {r.isotype}"
    assert o.anticodon == r.anticodon, f"anticodon mismatch: {o.anticodon} != {r.anticodon}"
    assert abs(o.score - r.score) <= SCORE_TOL, f"score mismatch: {o.score} vs {r.score}"


@pytest.mark.differential
@pytest.mark.requires_infernal
@pytest.mark.requires_trnascan
def test_detection_count_agrees(yeast_phe_fasta: Path) -> None:
    """Every reference call is matched by one of ours (and vice versa)."""
    ours = scan(yeast_phe_fasta, domain="euk")
    ref = run_trnascan(yeast_phe_fasta, domain="euk")

    for r in ref:
        assert any(_match(o, r) for o in ours), f"reference call not found by us: {r!r}"
    for o in ours:
        assert any(_match(o, r) for r in ref), f"our call not in reference: {o!r}"


@pytest.mark.differential
@pytest.mark.requires_infernal
@pytest.mark.requires_trnascan
def test_pseudo_score_decomposition_matches_reference(yeast_phe_fasta: Path) -> None:
    """Our primary(HMM)/secondary score split matches reference --breakdown.

    This is the basis of faithful pseudogene calling: hmm_score from the
    no-secondary-structure model, ss_score = total - hmm.
    """
    from trnascan_py.models_registry import resolve_domain, resolve_ns
    from trnascan_py.oracle import run_trnascan_breakdown
    from trnascan_py.pipeline import (
        DEFAULT_CM_CUTOFF,
        _anticodons_from_cmsearch,
        _glocal_locus_scores,
    )

    general_cm, _ = resolve_domain("euk")
    raw, _ac, _intr = _anticodons_from_cmsearch(general_cm, yeast_phe_fasta, None)
    cand = [h for h in raw if h.score >= DEFAULT_CM_CUTOFF]
    hmm = _glocal_locus_scores(yeast_phe_fasta, cand, resolve_ns("euk"), None)
    ref = run_trnascan_breakdown(yeast_phe_fasta, domain="euk")
    assert len(cand) == len(ref) == 1

    our_total, our_hmm = cand[0].score, hmm[0]
    our_ss = our_total - our_hmm
    _rhit, ref_hmm, ref_ss = ref[0]
    assert abs(our_hmm - ref_hmm) <= 2.0, f"HMM: ours={our_hmm} ref={ref_hmm}"
    assert abs(our_ss - ref_ss) <= 2.0, f"2'Str: ours={our_ss} ref={ref_ss}"


@pytest.mark.differential
@pytest.mark.requires_infernal
@pytest.mark.requires_trnascan
def test_suppressor_matches_reference(amber_suppressor_fasta: Path) -> None:
    ours = scan(amber_suppressor_fasta, domain="euk")
    ref = run_trnascan(amber_suppressor_fasta, domain="euk")
    assert len(ours) == len(ref) == 1, f"ours={ours!r} ref={ref!r}"
    o, r = ours[0], ref[0]
    assert _match(o, r), f"coordinate/strand mismatch: ours={o!r} ref={r!r}"
    # Reference calls an amber suppressor "Sup" with anticodon CTA.
    assert r.isotype == "Sup", f"reference type was {r.isotype}"
    assert o.isotype == r.isotype, f"isotype: ours={o.isotype} ref={r.isotype}"
    assert o.anticodon == r.anticodon, f"anticodon: ours={o.anticodon} ref={r.anticodon}"


@pytest.mark.differential
@pytest.mark.requires_infernal
@pytest.mark.requires_trnascan
def test_mito_matches_reference(mtdna_fasta: Path) -> None:
    """Mitochondrial scan agrees with reference -M on detection, coords, isotype.

    (Anticodon for non-Leu/Ser mito tRNAs is a documented v1 limitation, so it is
    not asserted here.)
    """
    ours = scan(mtdna_fasta, domain="mito-vert")
    ref = run_trnascan(mtdna_fasta, domain="mito-vert")
    assert len(ours) == len(ref) == 22, f"ours={len(ours)} ref={len(ref)}"

    for r in ref:
        matches = [o for o in ours if _match(o, r)]
        assert matches, f"reference call not found by us: {r!r}"
        assert matches[0].isotype == r.isotype, (
            f"isotype @{r.start}: ours={matches[0].isotype} ref={r.isotype}"
        )


@pytest.mark.differential
@pytest.mark.requires_infernal
@pytest.mark.requires_trnascan
def test_intron_coords_match_reference(intron_trna_fasta: Path) -> None:
    ours = scan(intron_trna_fasta, domain="euk")
    ref = run_trnascan(intron_trna_fasta, domain="euk")
    assert len(ours) == len(ref) == 1, f"ours={ours!r} ref={ref!r}"
    o, r = ours[0], ref[0]
    assert _match(o, r), f"coordinate/strand mismatch: ours={o!r} ref={r!r}"
    assert r.intron_start is not None and r.intron_end is not None, "reference reported no intron"
    assert o.intron_start is not None and o.intron_end is not None, "we reported no intron"
    assert abs(o.intron_start - r.intron_start) <= COORD_TOL, (
        f"intron start: ours={o.intron_start} ref={r.intron_start}"
    )
    assert abs(o.intron_end - r.intron_end) <= COORD_TOL, (
        f"intron end: ours={o.intron_end} ref={r.intron_end}"
    )
