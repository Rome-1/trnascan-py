"""Locating the bundled tRNA covariance models that ship with tRNAscan-SE.

v1 reuses the covariance models distributed with tRNAscan-SE 2.0 rather than
shipping its own copies. The models directory is discovered from (in order):

1. the ``TRNASCAN_MODELS_DIR`` environment variable,
2. the standard package install location (``/usr/share/trnascan-se/models``),
3. a sibling ``models`` directory next to a ``tRNAscan-SE`` on ``PATH``.
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path

_DEFAULT_DIRS = (
    "/usr/share/trnascan-se/models",
    "/usr/local/share/trnascan-se/models",
    "/opt/homebrew/share/trnascan-se/models",
)


class ModelsNotFoundError(RuntimeError):
    """Raised when the tRNA covariance-model directory cannot be located."""


@dataclass(frozen=True, slots=True)
class DomainModels:
    """Covariance models for one taxonomic domain.

    Attributes:
        general: General (isotype-merged) CM used for first-pass + confirmation.
        isotype_db: Pressed CM database of isotype-specific models for ``cmscan``.
    """

    general: str
    isotype_db: str


# Map a domain key to its model filenames (relative to the models directory).
DOMAIN_MODELS: dict[str, DomainModels] = {
    "euk": DomainModels(general="TRNAinf-euk.cm", isotype_db="TRNAinf-euk-iso"),
    "bact": DomainModels(general="TRNAinf-bact.cm", isotype_db="TRNAinf-bact-iso"),
    "arch": DomainModels(general="TRNAinf-arch.cm", isotype_db="TRNAinf-arch-iso"),
}

# Mitochondrial domains use a single multi-model (per-isotype) CM database for
# both the candidate search and the glocal isotype re-score — there is no merged
# "general" mito CM. Models are lineage-specific.
MITO_MODELS: dict[str, str] = {
    "mito-vert": "TRNAinf-mito-vert",
    "mito-mammal": "TRNAinf-mito-mammal",
}

# No-secondary-structure ("NS") models: the general CM with base pairs removed.
# Scoring against these gives the primary-structure ("HMM") score; the secondary-
# structure score is (total CM score - NS score). Used for pseudogene detection.
NS_MODELS: dict[str, str] = {
    "euk": "TRNAinf-euk-ns.cm",
    "bact": "TRNAinf-bact-ns.cm",
    "arch": "TRNAinf-arch-ns.cm",
}


def find_models_dir(override: str | Path | None = None) -> Path:
    """Locate the directory holding the bundled tRNA covariance models."""
    if override is not None:
        p = Path(override)
        if p.is_dir():
            return p
        raise ModelsNotFoundError(f"TRNASCAN models dir not found: {p}")

    env = os.environ.get("TRNASCAN_MODELS_DIR")
    candidates: list[Path] = []
    if env:
        candidates.append(Path(env))
    candidates += [Path(d) for d in _DEFAULT_DIRS]

    trnascan = shutil.which("tRNAscan-SE")
    if trnascan:
        prefix = Path(trnascan).resolve().parent.parent
        candidates.append(prefix / "share" / "trnascan-se" / "models")

    for c in candidates:
        if c.is_dir():
            return c
    raise ModelsNotFoundError(
        "Could not locate tRNA covariance models. Set TRNASCAN_MODELS_DIR or "
        "install tRNAscan-SE 2.0 (e.g. `apt-get install trnascan-se`)."
    )


def resolve_domain(domain: str, models_dir: str | Path | None = None) -> tuple[Path, Path]:
    """Return ``(general_cm_path, isotype_db_path)`` for ``domain``.

    Raises:
        KeyError: if ``domain`` is not one of ``euk``/``bact``/``arch``.
        ModelsNotFoundError: if the models directory or a model file is missing.
    """
    if domain not in DOMAIN_MODELS:
        raise KeyError(f"unknown domain {domain!r}; expected one of {sorted(DOMAIN_MODELS)}")
    base = find_models_dir(models_dir)
    spec = DOMAIN_MODELS[domain]
    general = base / spec.general
    isotype = base / spec.isotype_db
    if not general.exists():
        raise ModelsNotFoundError(f"general model missing: {general}")
    if not isotype.exists():
        raise ModelsNotFoundError(f"isotype model DB missing: {isotype}")
    return general, isotype


def resolve_mito(domain: str, models_dir: str | Path | None = None) -> Path:
    """Return the multi-model mitochondrial CM database path for ``domain``.

    Raises:
        KeyError: if ``domain`` is not a known mito domain.
        ModelsNotFoundError: if the database file is missing.
    """
    if domain not in MITO_MODELS:
        raise KeyError(f"unknown mito domain {domain!r}; expected one of {sorted(MITO_MODELS)}")
    db = find_models_dir(models_dir) / MITO_MODELS[domain]
    if not db.exists():
        raise ModelsNotFoundError(f"mito model DB missing: {db}")
    return db


def resolve_ns(domain: str, models_dir: str | Path | None = None) -> Path:
    """Return the no-secondary-structure CM path for ``domain`` (for pseudo scoring).

    Raises:
        KeyError: if ``domain`` has no NS model.
        ModelsNotFoundError: if the NS model file is missing.
    """
    if domain not in NS_MODELS:
        raise KeyError(f"no NS model for domain {domain!r}; expected one of {sorted(NS_MODELS)}")
    ns = find_models_dir(models_dir) / NS_MODELS[domain]
    if not ns.exists():
        raise ModelsNotFoundError(f"NS model missing: {ns}")
    return ns
