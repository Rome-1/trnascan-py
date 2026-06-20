"""Wrapper around the reference tRNAscan-SE 2.0 tool — the differential oracle.

v1 is validated by comparing :func:`trnascan_py.scan` output against the
reference implementation on test genomes. This module shells out to the
reference ``tRNAscan-SE`` and parses its tabular output into :class:`TRNAHit`
records so the two can be compared field by field.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from trnascan_py.models import Strand, TRNAHit

TRNASCAN_BIN = "tRNAscan-SE"

# Domain key -> reference tool search-mode flags.
_DOMAIN_FLAG = {
    "euk": ["-E"],
    "bact": ["-B"],
    "arch": ["-A"],
    "general": ["-G"],
    "mito-vert": ["-M", "vert"],
    "mito-mammal": ["-M", "mammal"],
}


class TrnascanNotFoundError(RuntimeError):
    """Raised when the reference ``tRNAscan-SE`` executable is not on ``PATH``."""


class TrnascanError(RuntimeError):
    """Raised when a reference ``tRNAscan-SE`` invocation exits non-zero."""


def trnascan_path() -> str:
    """Return the resolved path to ``tRNAscan-SE`` or raise."""
    path = shutil.which(TRNASCAN_BIN)
    if path is None:
        raise TrnascanNotFoundError(
            f"{TRNASCAN_BIN!r} not found on PATH. Install reference tRNAscan-SE 2.0 "
            "(e.g. `apt-get install trnascan-se`)."
        )
    return path


def trnascan_version() -> str:
    """Return the reference tool's version string."""
    out = subprocess.run(  # noqa: S603 - fixed argv, resolved binary
        [trnascan_path(), "-h"],
        capture_output=True,
        text=True,
        check=False,
    )
    for line in (out.stdout + out.stderr).splitlines():
        if line.lower().startswith("trnascan-se"):
            return line.strip()
    return "unknown"


def run_trnascan(
    fasta: str | Path,
    *,
    domain: str = "euk",
    extra_args: list[str] | None = None,
    timeout: float | None = None,
) -> list[TRNAHit]:
    """Run the reference tRNAscan-SE and parse its tabular output.

    Args:
        fasta: Input sequence FASTA path.
        domain: Search mode — ``"euk"``, ``"bact"``, ``"arch"``, or ``"general"``.
        extra_args: Extra raw arguments appended to the command.
        timeout: Optional subprocess timeout in seconds.

    Returns:
        The reference tool's tRNA calls as :class:`TRNAHit` records.
    """
    fasta = Path(fasta)
    if not fasta.exists():
        raise FileNotFoundError(f"sequence file not found: {fasta}")
    if domain not in _DOMAIN_FLAG:
        raise KeyError(f"unknown domain {domain!r}; expected one of {sorted(_DOMAIN_FLAG)}")

    # tRNAscan-SE blocks on an interactive "overwrite?" prompt if its -o output
    # file already exists, so we point it at a path that does NOT exist yet (and
    # close stdin) rather than pre-creating a temp file.
    with tempfile.TemporaryDirectory() as tmp:
        out_path = Path(tmp) / "trnascan.out"
        cmd: list[str] = [trnascan_path(), *_DOMAIN_FLAG[domain], "-q", "-o", str(out_path)]
        if extra_args:
            cmd += extra_args
        cmd += [str(fasta)]

        proc = subprocess.run(  # noqa: S603 - argv from resolved binary + paths
            cmd,
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout,
        )
        if proc.returncode != 0:
            raise TrnascanError(
                f"tRNAscan-SE exited {proc.returncode}: "
                f"{proc.stderr.strip() or proc.stdout.strip()}"
            )
        return parse_trnascan_out(out_path.read_text())


def parse_trnascan_out(text: str) -> list[TRNAHit]:
    """Parse tRNAscan-SE tabular output into :class:`TRNAHit` records.

    The tabular format has a three-line header followed by one row per tRNA::

        Name  tRNA#  Begin  End  Type  Codon  IntronBegin  IntronEnd  Score  [Note]

    Begin/End are in genomic orientation; ``Begin > End`` denotes the minus
    strand. Intron coordinates of ``0`` mean "no intron".
    """
    hits: list[TRNAHit] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        fields = line.split("\t") if "\t" in raw else line.split()
        fields = [f.strip() for f in fields if f.strip() != ""]
        # Skip header / separator rows: data rows have a numeric tRNA# in col 1
        # and numeric coordinates.
        if len(fields) < 9:
            continue
        try:
            begin = int(fields[2])
            end = int(fields[3])
            int_begin = int(fields[6])
            int_end = int(fields[7])
            score = float(fields[8])
        except (ValueError, IndexError):
            continue

        strand = Strand.PLUS if begin <= end else Strand.MINUS
        note = fields[9] if len(fields) > 9 else ""
        hits.append(
            TRNAHit(
                seq_id=fields[0],
                start=begin,
                end=end,
                strand=strand,
                isotype=fields[4],
                anticodon=fields[5],
                score=score,
                intron_start=int_begin or None,
                intron_end=int_end or None,
                note=note,
            )
        )
    hits.sort(key=lambda h: (h.seq_id, min(h.start, h.end)))
    return hits


def run_trnascan_breakdown(
    fasta: str | Path, *, domain: str = "euk", timeout: float | None = None
) -> list[tuple[TRNAHit, float, float]]:
    """Run reference tRNAscan-SE with ``--breakdown`` and parse score components.

    Returns ``(hit, hmm_score, ss_score)`` per tRNA, where ``hmm_score`` is the
    reference's "HMM Score" (primary structure) and ``ss_score`` is its "2'Str
    Score" (secondary structure). Used to validate the pseudogene decomposition.
    """
    fasta = Path(fasta)
    if not fasta.exists():
        raise FileNotFoundError(f"sequence file not found: {fasta}")
    if domain not in _DOMAIN_FLAG:
        raise KeyError(f"unknown domain {domain!r}; expected one of {sorted(_DOMAIN_FLAG)}")

    with tempfile.TemporaryDirectory() as tmp:
        out_path = Path(tmp) / "trnascan.out"
        cmd: list[str] = [
            trnascan_path(), *_DOMAIN_FLAG[domain], "-q", "--breakdown", "-o", str(out_path)
        ]
        cmd += [str(fasta)]
        proc = subprocess.run(  # noqa: S603 - argv from resolved binary + paths
            cmd,
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout,
        )
        if proc.returncode != 0:
            raise TrnascanError(
                f"tRNAscan-SE exited {proc.returncode}: "
                f"{proc.stderr.strip() or proc.stdout.strip()}"
            )
        return _parse_breakdown(out_path.read_text())


def _parse_breakdown(text: str) -> list[tuple[TRNAHit, float, float]]:
    """Parse ``--breakdown`` output: ... Score HMMScore 2'StrScore [Note]."""
    out: list[tuple[TRNAHit, float, float]] = []
    for raw in text.splitlines():
        fields = [f.strip() for f in raw.split("\t") if f.strip()] if "\t" in raw else raw.split()
        if len(fields) < 11:
            continue
        try:
            begin, end = int(fields[2]), int(fields[3])
            int_begin, int_end = int(fields[6]), int(fields[7])
            score, hmm, ss = float(fields[8]), float(fields[9]), float(fields[10])
        except (ValueError, IndexError):
            continue
        hit = TRNAHit(
            seq_id=fields[0],
            start=begin,
            end=end,
            strand=Strand.PLUS if begin <= end else Strand.MINUS,
            isotype=fields[4],
            anticodon=fields[5],
            score=score,
            intron_start=int_begin or None,
            intron_end=int_end or None,
            note=fields[11] if len(fields) > 11 else "",
        )
        out.append((hit, hmm, ss))
    return out


__all__ = [
    "run_trnascan",
    "run_trnascan_breakdown",
    "parse_trnascan_out",
    "trnascan_version",
    "TRNASCAN_BIN",
]
