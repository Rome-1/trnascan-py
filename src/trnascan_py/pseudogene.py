"""Faithful pseudogene flagging, matching reference tRNAscan-SE 2.0.

Reference tRNAscan-SE (``tRNAscanSE/CM.pm::is_pseudo_gene``) decomposes a tRNA's
covariance-model bit score into a **primary** ("HMM") component and a
**secondary-structure** component, then flags low-scoring hits whose structure or
sequence score falls below a floor:

    hmm_score = score of the hit against the no-secondary-structure (NS) model
    ss_score  = total_cm_score - hmm_score
    pseudogene  iff  total_cm_score < PSEUDO_FILTER_SCORE
                     and (ss_score < MIN_SS_SCORE or hmm_score < MIN_HMM_SCORE)

The NS model is the general CM with base pairs removed; scoring against it (glocal
``cmsearch``) reproduces the reference's "HMM Score" column exactly. The pipeline
supplies ``total_score`` and ``hmm_score``; this module holds the thresholds and
the decision. Values are tRNAscan-SE 2.0 defaults (``min_cmsearch_pseudo_filter_score``
= 55, ``min_ss_score`` = 5, ``min_hmm_score`` = 10).
"""

from __future__ import annotations

PSEUDO_FILTER_SCORE = 55.0
MIN_SS_SCORE = 5.0
MIN_HMM_SCORE = 10.0


def secondary_structure_score(total_score: float, hmm_score: float) -> float:
    """Secondary-structure score = total CM score minus the primary (NS) score."""
    return total_score - hmm_score


def is_pseudogene(total_score: float, hmm_score: float) -> bool:
    """Whether a hit is a pseudogene, per reference tRNAscan-SE's criteria.

    Args:
        total_score: the full covariance-model bit score.
        hmm_score: the bit score against the no-secondary-structure model.
    """
    if total_score >= PSEUDO_FILTER_SCORE:
        return False
    ss_score = secondary_structure_score(total_score, hmm_score)
    return ss_score < MIN_SS_SCORE or hmm_score < MIN_HMM_SCORE


__all__ = [
    "is_pseudogene",
    "secondary_structure_score",
    "PSEUDO_FILTER_SCORE",
    "MIN_SS_SCORE",
    "MIN_HMM_SCORE",
]
