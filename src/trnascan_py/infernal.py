"""Thin subprocess wrapper around Infernal's ``cmsearch``.

This is the v1 engine binding: we shell out to ``cmsearch``, run a covariance
model against a sequence database, and parse the ``--tblout`` table into typed
records. v2 will replace this with a native accelerated DP core; the parsed
record shape is intended to stay stable across that transition.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from trnascan_py.models import Strand

CMSEARCH_BIN = "cmsearch"
CMSCAN_BIN = "cmscan"


class InfernalNotFoundError(RuntimeError):
    """Raised when an Infernal executable cannot be located on ``PATH``."""


class CmsearchError(RuntimeError):
    """Raised when a ``cmsearch``/``cmscan`` invocation exits non-zero."""


@dataclass(frozen=True, slots=True)
class CmsearchHit:
    """One row of ``cmsearch --tblout`` output (the per-CM tabular format).

    Coordinates are 1-based inclusive in genomic orientation, matching Infernal:
    on the minus strand ``seq_from > seq_to``.
    """

    target_name: str
    query_name: str
    seq_from: int
    seq_to: int
    strand: Strand
    score: float
    evalue: float
    inc: str
    """Inclusion marker: ``"!"`` (passes threshold) or ``"?"`` (reported only)."""

    mdl_from: int = 0
    mdl_to: int = 0
    trunc: str = "-"
    description: str = ""


def _resolve(binary: str) -> str:
    path = shutil.which(binary)
    if path is None:
        raise InfernalNotFoundError(
            f"{binary!r} not found on PATH. Install Infernal "
            "(e.g. `apt-get install infernal` or `brew install infernal`)."
        )
    return path


def cmsearch_path() -> str:
    """Return the resolved path to ``cmsearch`` or raise :class:`InfernalNotFoundError`."""
    return _resolve(CMSEARCH_BIN)


def cmscan_path() -> str:
    """Return the resolved path to ``cmscan`` or raise :class:`InfernalNotFoundError`."""
    return _resolve(CMSCAN_BIN)


def cmsearch_version() -> str:
    """Return the Infernal version string reported by ``cmsearch -h``."""
    out = subprocess.run(  # noqa: S603 - fixed argv, resolved binary
        [cmsearch_path(), "-h"],
        capture_output=True,
        text=True,
        check=False,
    )
    for line in out.stdout.splitlines():
        if "INFERNAL" in line.upper():
            return line.strip("# ").strip()
    return "unknown"


def run_cmsearch(
    cm_path: str | Path,
    seq_path: str | Path,
    *,
    cpu: int | None = None,
    incE: float | None = None,
    extra_args: list[str] | None = None,
    timeout: float | None = None,
    alignment_out: str | Path | None = None,
) -> list[CmsearchHit]:
    """Run ``cmsearch`` for one covariance model against one sequence file.

    Args:
        cm_path: Path to the covariance model (``.cm``) file.
        seq_path: Path to the target sequence FASTA file.
        cpu: ``--cpu`` worker count; ``None`` leaves Infernal's default.
        incE: ``--incE`` inclusion E-value threshold; ``None`` leaves the default.
        extra_args: Additional raw arguments appended to the command.
        timeout: Optional subprocess timeout in seconds.
        alignment_out: If given, write the Stockholm alignment of hits (``-A``)
            to this path (used for anticodon extraction).

    Returns:
        Parsed hits from the ``--tblout`` table.

    Raises:
        InfernalNotFoundError: if ``cmsearch`` is not on ``PATH``.
        CmsearchError: if ``cmsearch`` exits non-zero.
        FileNotFoundError: if ``cm_path`` or ``seq_path`` does not exist.
    """
    cm_path = Path(cm_path)
    seq_path = Path(seq_path)
    if not cm_path.exists():
        raise FileNotFoundError(f"covariance model not found: {cm_path}")
    if not seq_path.exists():
        raise FileNotFoundError(f"sequence file not found: {seq_path}")

    import tempfile

    with tempfile.NamedTemporaryFile("r", suffix=".tblout", delete=True) as tbl:
        cmd: list[str] = [cmsearch_path(), "--tblout", tbl.name]
        if alignment_out is not None:
            cmd += ["-A", str(alignment_out)]
        if cpu is not None:
            cmd += ["--cpu", str(cpu)]
        if incE is not None:
            cmd += ["--incE", str(incE)]
        if extra_args:
            cmd += extra_args
        cmd += [str(cm_path), str(seq_path)]

        proc = subprocess.run(  # noqa: S603 - argv built from resolved binary + paths
            cmd,
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout,
        )
        if proc.returncode != 0:
            raise CmsearchError(
                f"cmsearch exited {proc.returncode}: {proc.stderr.strip() or proc.stdout.strip()}"
            )
        return parse_tblout(Path(tbl.name).read_text())


def run_cmscan(
    cm_db_path: str | Path,
    seq_path: str | Path,
    *,
    cpu: int | None = None,
    extra_args: list[str] | None = None,
    timeout: float | None = None,
) -> list[CmsearchHit]:
    """Run ``cmscan`` of a sequence file against a pressed CM database.

    Used for isotype-specific classification. The returned :class:`CmsearchHit`
    records reuse the same shape, but for ``cmscan`` output ``target_name`` is the
    **model** name (e.g. ``"euk-Phe"``) and ``query_name`` is the **sequence**
    name. ``seq_from``/``seq_to`` are still the coordinates within the sequence.

    The CM database must be pressed (``cmpress``); the tRNAscan-SE ``*-iso``
    databases ship pre-pressed.
    """
    cm_db_path = Path(cm_db_path)
    seq_path = Path(seq_path)
    if not cm_db_path.exists():
        raise FileNotFoundError(f"CM database not found: {cm_db_path}")
    if not seq_path.exists():
        raise FileNotFoundError(f"sequence file not found: {seq_path}")

    import tempfile

    with tempfile.NamedTemporaryFile("r", suffix=".tblout", delete=True) as tbl:
        cmd: list[str] = [cmscan_path(), "--tblout", tbl.name]
        if cpu is not None:
            cmd += ["--cpu", str(cpu)]
        if extra_args:
            cmd += extra_args
        cmd += [str(cm_db_path), str(seq_path)]

        proc = subprocess.run(  # noqa: S603 - argv built from resolved binary + paths
            cmd,
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout,
        )
        if proc.returncode != 0:
            raise CmsearchError(
                f"cmscan exited {proc.returncode}: {proc.stderr.strip() or proc.stdout.strip()}"
            )
        return parse_tblout(Path(tbl.name).read_text())


def parse_tblout(text: str) -> list[CmsearchHit]:
    """Parse the whitespace-delimited ``cmsearch --tblout`` table.

    The per-CM ``--tblout`` format has these leading columns (Infernal 1.1.x)::

        target_name accession query_name accession mdl mdl_from mdl_to
        seq_from seq_to strand trunc pass gc bias score E-value inc <desc...>

    Lines beginning with ``#`` are comments and are skipped.
    """
    hits: list[CmsearchHit] = []
    for raw in text.splitlines():
        line = raw.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 17:
            continue
        try:
            hit = CmsearchHit(
                target_name=fields[0],
                query_name=fields[2],
                mdl_from=int(fields[5]),
                mdl_to=int(fields[6]),
                seq_from=int(fields[7]),
                seq_to=int(fields[8]),
                strand=Strand.from_str(fields[9]),
                trunc=fields[10],
                score=float(fields[14]),
                evalue=float(fields[15]),
                inc=fields[16],
                description=" ".join(fields[17:]) if len(fields) > 17 else "",
            )
        except (ValueError, IndexError):
            # Malformed row — skip rather than abort the whole parse.
            continue
        hits.append(hit)
    return hits
