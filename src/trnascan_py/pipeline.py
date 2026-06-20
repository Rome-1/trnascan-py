"""The tRNAscan-SE-style detection pipeline, driven via Infernal.

v1 mirrors the reference pipeline's shape:

1. **First pass + CM confirmation** — ``cmsearch`` with the general (isotype-merged)
   covariance model for the domain; hits at or above ``cm_cutoff`` bits are kept.
2. **Anticodon extraction** — the ``cmsearch -A`` Stockholm alignment is parsed and
   the anticodon triplet is read from the anticodon loop (see :mod:`trnascan_py.anticodon`).
3. **Isotype classification** — by default the anticodon determines the isotype via
   the genetic code (:mod:`trnascan_py.genetic_code`), which is fast and matches the
   reference for standard cases. Anticodons the code cannot disambiguate (``CAT`` →
   Met/fMet/Ile2/iMet, ``TCA`` → SeC/Sup) are refined with a small CM scan of just
   those loci. See ``isotype_method``.

The reported ``TRNAHit.score`` is the general-model bit score, which matches the
reference tool's "Inf Score" column within tolerance.
"""

from __future__ import annotations

import tempfile
from collections.abc import Callable
from pathlib import Path
from typing import Literal

from trnascan_py.anticodon import anticodon_from_alignment, parse_stockholm
from trnascan_py.faidx import FastaIndex, UnindexableFastaError
from trnascan_py.fasta import extract, read_fasta, reverse_complement, write_fasta
from trnascan_py.genetic_code import AMBIGUOUS_ANTICODONS, isotype_from_anticodon
from trnascan_py.infernal import CmsearchHit, run_cmsearch
from trnascan_py.intron import intron_from_alignment
from trnascan_py.models import Strand, TRNAHit
from trnascan_py.models_registry import (
    MITO_MODELS,
    NS_MODELS,
    resolve_domain,
    resolve_mito,
    resolve_ns,
)
from trnascan_py.pseudogene import PSEUDO_FILTER_SCORE, is_pseudogene

IsotypeMethod = Literal["hybrid", "anticodon", "cm"]

# Default reporting cutoff (bit score), matching tRNAscan-SE's ``cm_cutoff``.
DEFAULT_CM_CUTOFF = 20.0
# Reporting cutoff for organellar/mito tRNAs (tRNAscan-SE ``organelle_cm_cutoff``).
ORGANELLE_CM_CUTOFF = 15.0
# Flanking bases kept around each candidate locus when extracting for the
# isotype scan, so the isotype model has full context to align against.
ISOTYPE_LOCUS_PAD = 20
# Fast glocal isotype-scan flags, mirroring reference tRNAscan-SE: no HMM filter
# (candidates are already real tRNAs, so the filter only adds overhead), top
# strand only (candidates are extracted in tRNA orientation), no truncation.
_ISOTYPE_CMSEARCH_ARGS = ["-g", "--nohmm", "--toponly", "--notrunc"]


def _isotype_from_model(model_name: str) -> str:
    """``"euk-Phe"`` -> ``"Phe"``; ``"bact-iMet"`` -> ``"iMet"``."""
    for prefix in ("euk-", "bact-", "arch-"):
        if model_name.startswith(prefix):
            return model_name[len(prefix) :]
    return model_name


def _locus_fetcher(fasta: str | Path) -> Callable[[str, int, int], str]:
    """Return a ``fetch(seqid, start, end) -> forward subsequence`` function.

    Prefers a streaming :class:`~trnascan_py.faidx.FastaIndex` (O(1) memory — does
    not hold the genome in RAM), which matters for multi-Gb assemblies. Falls
    back to an in-memory read for FASTA that cannot be index-seeked (non-uniform
    line widths).
    """
    try:
        idx = FastaIndex.build(fasta)
    except UnindexableFastaError:
        genome = read_fasta(fasta)

        def fetch_mem(seqid: str, start: int, end: int) -> str:
            seq = genome.get(seqid)
            return extract(seq, start, end) if seq else ""

        return fetch_mem

    def fetch_idx(seqid: str, start: int, end: int) -> str:
        return idx.fetch(seqid, start, end) if seqid in idx else ""

    return fetch_idx


def _classify_isotypes(
    fasta: str | Path,
    candidates: list[CmsearchHit],
    iso_db: Path,
    cpu: int | None,
) -> dict[int, str]:
    """Assign an isotype to each candidate by scanning only the candidate loci.

    Rather than ``cmscan`` over the whole genome (the dominant cost in early v1 —
    see docs/PROFILE.md), we extract each confirmed candidate **in tRNA
    orientation** (reverse-complemented for minus-strand hits) into a tiny FASTA,
    named ``cand_{i}`` by candidate index, and run a single ``cmsearch`` of all
    isotype models against it in fast glocal mode (no HMM filter, top strand
    only). The best-scoring model per ``cand_{i}`` gives that candidate's
    isotype — an exact index match, no coordinate-overlap heuristic.

    Loci are pulled via a streaming FASTA index so the genome is never held in
    memory (scales to multi-Gb assemblies).
    """
    fetch = _locus_fetcher(fasta)
    records: list[tuple[str, str]] = []
    for i, c in enumerate(candidates):
        lo = min(c.seq_from, c.seq_to) - ISOTYPE_LOCUS_PAD
        hi = max(c.seq_from, c.seq_to) + ISOTYPE_LOCUS_PAD
        sub = fetch(c.target_name, lo, hi)
        if not sub:
            continue
        if c.strand is Strand.MINUS:
            sub = reverse_complement(sub)
        records.append((f"cand_{i}", sub))
    if not records:
        return {}

    with tempfile.TemporaryDirectory() as tmp:
        loci_fa = Path(tmp) / "candidates.fa"
        write_fasta(records, loci_fa)
        iso_hits = run_cmsearch(
            iso_db, loci_fa, cpu=cpu, extra_args=_ISOTYPE_CMSEARCH_ARGS
        )

    best_score: dict[int, float] = {}
    assignment: dict[int, str] = {}
    for iso in iso_hits:
        # cmsearch: target_name = sequence ("cand_{i}"), query_name = isotype model.
        if not iso.target_name.startswith("cand_"):
            continue
        idx = int(iso.target_name[len("cand_") :])
        if idx not in best_score or iso.score > best_score[idx]:
            best_score[idx] = iso.score
            assignment[idx] = _isotype_from_model(iso.query_name)
    return assignment


def _anticodons_from_cmsearch(
    general_cm: Path, fasta: str | Path, cpu: int | None
) -> tuple[list[CmsearchHit], dict[str, str], dict[str, tuple[int, int]]]:
    """Run ``cmsearch`` with alignment capture; extract per-hit anticodon + intron.

    Returns the raw hits, a map ``"seqid/from-to" -> anticodon``, and a map
    ``"seqid/from-to" -> (intron_start, intron_end)`` for hits with an intron.
    """
    with tempfile.TemporaryDirectory() as tmp:
        aln_path = Path(tmp) / "hits.sto"
        raw = run_cmsearch(general_cm, fasta, cpu=cpu, alignment_out=aln_path)
        anticodons: dict[str, str] = {}
        introns: dict[str, tuple[int, int]] = {}
        if aln_path.exists():
            for sto in parse_stockholm(aln_path.read_text()):
                ac = anticodon_from_alignment(sto)
                if ac is not None:
                    anticodons[sto.name] = ac
                intron = intron_from_alignment(sto)
                if intron is not None:
                    introns[sto.name] = intron
    return raw, anticodons, introns


def _glocal_locus_scores(
    fasta: str | Path, candidates: list[CmsearchHit], cm: Path, cpu: int | None
) -> dict[int, float]:
    """Glocally score each candidate locus against ``cm``; return best per index.

    Mirrors reference tRNAscan-SE's ``cmsearch_scoring`` (glocal, no HMM filter) on
    the extracted candidate. Against the general model this reproduces the refined
    "Inf Score"; against the no-secondary-structure model it reproduces the
    "HMM Score".
    """
    fetch = _locus_fetcher(fasta)
    records: list[tuple[str, str]] = []
    for i, c in enumerate(candidates):
        lo = min(c.seq_from, c.seq_to) - ISOTYPE_LOCUS_PAD
        hi = max(c.seq_from, c.seq_to) + ISOTYPE_LOCUS_PAD
        sub = fetch(c.target_name, lo, hi)
        if not sub:
            continue
        if c.strand is Strand.MINUS:
            sub = reverse_complement(sub)
        records.append((f"cand_{i}", sub))
    if not records:
        return {}

    with tempfile.TemporaryDirectory() as tmp:
        loci_fa = Path(tmp) / "loci.fa"
        write_fasta(records, loci_fa)
        hits = run_cmsearch(cm, loci_fa, cpu=cpu, extra_args=_ISOTYPE_CMSEARCH_ARGS)

    scores: dict[int, float] = {}
    for h in hits:
        if not h.target_name.startswith("cand_"):
            continue
        idx = int(h.target_name[len("cand_") :])
        if idx not in scores or h.score > scores[idx]:
            scores[idx] = h.score
    return scores


def _hit_key(c: CmsearchHit) -> str:
    """Stockholm sequence name Infernal assigns a hit: ``"seqid/from-to"``."""
    return f"{c.target_name}/{c.seq_from}-{c.seq_to}"


def _assign_isotypes(
    fasta: str | Path,
    candidates: list[CmsearchHit],
    anticodons: dict[str, str],
    iso_db: Path,
    cpu: int | None,
    method: IsotypeMethod,
) -> dict[int, str]:
    """Assign an isotype to each candidate index using the requested method.

    * ``"anticodon"`` — genetic-code translation of the extracted anticodon (fast).
    * ``"cm"`` — covariance-model scan of every candidate locus (most faithful,
      slowest; ~all of v1 wall-clock — see docs/PROFILE.md).
    * ``"hybrid"`` (default) — anticodon for all, then a CM scan of *only* the
      ambiguous loci (``CAT``/``TCA``) and any the anticodon left ``Undet``.
    """
    if method == "cm":
        return _classify_isotypes(fasta, candidates, iso_db, cpu)

    base = {
        i: (isotype_from_anticodon(ac) if (ac := anticodons.get(_hit_key(c))) else "Undet")
        for i, c in enumerate(candidates)
    }
    if method == "anticodon":
        return base

    # hybrid: refine only the loci the anticodon can't resolve.
    refine = [
        i
        for i, c in enumerate(candidates)
        if anticodons.get(_hit_key(c), "").upper().replace("U", "T") in AMBIGUOUS_ANTICODONS
        or base[i] == "Undet"
    ]
    if refine:
        sub = [candidates[i] for i in refine]
        cm = _classify_isotypes(fasta, sub, iso_db, cpu)
        for local_i, global_i in enumerate(refine):
            if local_i in cm:
                base[global_i] = cm[local_i]
    return base


def _dedup_by_locus(hits: list[CmsearchHit]) -> list[CmsearchHit]:
    """Keep the highest-scoring hit per overlapping same-strand locus."""
    kept: list[CmsearchHit] = []
    for h in sorted(hits, key=lambda x: -x.score):
        lo, hi = sorted((h.seq_from, h.seq_to))
        if any(
            h.strand is k.strand
            and lo <= max(k.seq_from, k.seq_to)
            and min(k.seq_from, k.seq_to) <= hi
            for k in kept
        ):
            continue
        kept.append(h)
    return kept


def _parse_mito_model(name: str) -> tuple[str, str]:
    """``"LeuTAA"`` -> ``("Leu", "TAA")``; ``"Phe"`` -> ``("Phe", "")``.

    Mito isotype models encode the anticodon only for the degenerate amino acids
    (Leu, Ser); other anticodons are not recoverable from the model name.
    """
    base = name.split("-")[0]  # drop suffixes like "-no-darm"
    suffix = base[-3:]
    if len(base) > 3 and suffix.isupper() and all(b in "ACGT" for b in suffix):
        return base[:-3], suffix
    return base, ""


def _scan_mito(
    fasta: str | Path, mito_db: Path, cm_cutoff: float, cpu: int | None
) -> list[TRNAHit]:
    """Detect mitochondrial tRNAs using the lineage-specific mito CM database.

    Mito tRNAs are structurally divergent, so a single local search underscores
    them. We use a two-pass approach mirroring reference tRNAscan-SE's organellar
    mode: (1) a local ``cmsearch`` of the multi-model mito database finds candidate
    loci; (2) each locus is re-scored *glocally* against every mito model, and the
    best model gives the isotype and the authoritative score. Glocal scoring lets
    divergent tRNAs (e.g. vertebrate mito Tyr) align full-length and win their
    locus.
    """
    pass1 = [h for h in run_cmsearch(mito_db, fasta, cpu=cpu) if h.score >= cm_cutoff]
    loci = _dedup_by_locus(pass1)
    if not loci:
        return []

    fetch = _locus_fetcher(fasta)
    records: list[tuple[str, str]] = []
    index_of: list[int] = []
    for i, c in enumerate(loci):
        lo = min(c.seq_from, c.seq_to) - ISOTYPE_LOCUS_PAD
        hi = max(c.seq_from, c.seq_to) + ISOTYPE_LOCUS_PAD
        sub = fetch(c.target_name, lo, hi)
        if not sub:
            continue
        if c.strand is Strand.MINUS:
            sub = reverse_complement(sub)
        records.append((f"cand_{i}", sub))
        index_of.append(i)

    best: dict[int, tuple[str, float]] = {}
    if records:
        with tempfile.TemporaryDirectory() as tmp:
            loci_fa = Path(tmp) / "loci.fa"
            write_fasta(records, loci_fa)
            rescored = run_cmsearch(
                mito_db, loci_fa, cpu=cpu, extra_args=_ISOTYPE_CMSEARCH_ARGS
            )
        for h in rescored:
            if not h.target_name.startswith("cand_"):
                continue
            idx = int(h.target_name[len("cand_") :])
            if idx not in best or h.score > best[idx][1]:
                best[idx] = (h.query_name, h.score)

    results: list[TRNAHit] = []
    for i, c in enumerate(loci):
        model, score = best.get(i, (c.query_name, c.score))
        isotype, anticodon = _parse_mito_model(model)
        results.append(
            TRNAHit(
                seq_id=c.target_name,
                start=c.seq_from,
                end=c.seq_to,
                strand=c.strand,
                isotype=isotype,
                anticodon=anticodon,
                score=score,
            )
        )
    results.sort(key=lambda h: (h.seq_id, min(h.start, h.end)))
    return results


def scan(
    fasta: str | Path,
    *,
    domain: str = "euk",
    models_dir: str | Path | None = None,
    cm_cutoff: float = DEFAULT_CM_CUTOFF,
    classify_isotype: bool = True,
    isotype_method: IsotypeMethod = "hybrid",
    detect_pseudo: bool = True,
    cpu: int | None = None,
) -> list[TRNAHit]:
    """Scan a FASTA file for tRNA genes using the Infernal-backed pipeline.

    Args:
        fasta: Path to the input sequence FASTA file.
        domain: Model set — ``"euk"``, ``"bact"``, ``"arch"``, or a mitochondrial
            lineage ``"mito-vert"`` / ``"mito-mammal"``.
        models_dir: Override for the covariance-model directory (else auto-detected).
        cm_cutoff: Minimum general-model bit score to report a hit.
        classify_isotype: Whether to assign isotypes at all (else all ``"Undet"``).
        isotype_method: ``"hybrid"`` (default; anticodon + CM refinement of
            ambiguous loci), ``"anticodon"`` (genetic code only, fastest), or
            ``"cm"`` (covariance-model scan of every locus, slowest). Ignored for
            mito domains (isotype comes from the matched mito model).
        detect_pseudo: flag low-scoring, structurally-degenerate hits with
            ``note="pseudo"`` (conservative; never flags real tRNAs).
        cpu: ``--cpu`` worker count passed through to Infernal.

    Suppressor tRNAs (anticodons that decode a stop codon) are reported with
    isotype ``"Sup"``, matching reference tRNAscan-SE.

    Returns:
        Confirmed tRNA hits, sorted by sequence then start coordinate.

    Raises:
        InfernalNotFoundError: if Infernal is not installed.
        ModelsNotFoundError: if the covariance models cannot be located.
    """
    if domain in MITO_MODELS:
        mito_cutoff = cm_cutoff if cm_cutoff != DEFAULT_CM_CUTOFF else ORGANELLE_CM_CUTOFF
        return _scan_mito(fasta, resolve_mito(domain, models_dir), mito_cutoff, cpu)

    general_cm, iso_db = resolve_domain(domain, models_dir)

    raw, anticodons, introns = _anticodons_from_cmsearch(general_cm, fasta, cpu)
    candidates = [h for h in raw if h.score >= cm_cutoff]

    isotypes: dict[int, str] = {}
    if classify_isotype and candidates:
        isotypes = _assign_isotypes(fasta, candidates, anticodons, iso_db, cpu, isotype_method)

    # Pseudogene check (faithful to reference). Only hits below the pseudo-filter
    # score are checked. The local first-pass bounds can clip a few bp and under-
    # score a true tRNA, so we first glocally re-score the locus against the general
    # model — this recovers the reference's refined "Inf Score" — and use that as
    # the total for both the pseudo gate and the reported score. The primary
    # ("HMM") score comes from the no-secondary-structure model.
    pseudo_idx: set[int] = set()
    refined_score: dict[int, float] = {}
    if detect_pseudo and domain in NS_MODELS:
        low = [(i, c) for i, c in enumerate(candidates) if c.score < PSEUDO_FILTER_SCORE]
        if low:
            ns_cm = resolve_ns(domain, models_dir)
            sub = [c for _, c in low]
            total = _glocal_locus_scores(fasta, sub, general_cm, cpu)
            hmm = _glocal_locus_scores(fasta, sub, ns_cm, cpu)
            for local_i, (global_i, c) in enumerate(low):
                if local_i not in hmm:
                    continue
                t = max(c.score, total.get(local_i, c.score))
                refined_score[global_i] = t
                if is_pseudogene(t, hmm[local_i]):
                    pseudo_idx.add(global_i)

    results: list[TRNAHit] = []
    for i, c in enumerate(candidates):
        key = _hit_key(c)
        intron = introns.get(key)
        note = "pseudo" if i in pseudo_idx else ""
        results.append(
            TRNAHit(
                # For cmsearch output, the sequence name is the *target*.
                seq_id=c.target_name,
                start=c.seq_from,
                end=c.seq_to,
                strand=c.strand,
                isotype=isotypes.get(i, "Undet"),
                anticodon=anticodons.get(key, ""),
                score=refined_score.get(i, c.score),
                intron_start=intron[0] if intron else None,
                intron_end=intron[1] if intron else None,
                note=note,
            )
        )

    results.sort(key=lambda h: (h.seq_id, min(h.start, h.end)))
    return results


__all__ = ["scan", "DEFAULT_CM_CUTOFF", "IsotypeMethod", "Strand", "TRNAHit"]
