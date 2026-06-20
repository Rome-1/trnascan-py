"""The covariance models bundled in the package are present and resolvable."""

from __future__ import annotations

from trnascan_py.models_registry import (
    _BUNDLED_MODELS_DIR,
    resolve_domain,
    resolve_mito,
    resolve_ns,
)


def test_bundled_models_dir_exists() -> None:
    assert _BUNDLED_MODELS_DIR.is_dir(), f"bundled models missing at {_BUNDLED_MODELS_DIR}"
    assert (_BUNDLED_MODELS_DIR / "NOTICE.md").exists()


def test_bundled_general_and_iso_resolve() -> None:
    for domain in ("euk", "bact", "arch"):
        general, iso = resolve_domain(domain, models_dir=_BUNDLED_MODELS_DIR)
        assert general.exists() and iso.exists(), domain


def test_bundled_ns_models_resolve() -> None:
    for domain in ("euk", "bact", "arch"):
        assert resolve_ns(domain, models_dir=_BUNDLED_MODELS_DIR).exists(), domain


def test_bundled_mito_models_resolve() -> None:
    for domain in ("mito-vert", "mito-mammal"):
        assert resolve_mito(domain, models_dir=_BUNDLED_MODELS_DIR).exists(), domain


def test_bundled_dir_is_a_discovery_fallback() -> None:
    # The bundled directory is registered as the final discovery fallback.
    from trnascan_py import models_registry

    assert models_registry._BUNDLED_MODELS_DIR.name == "models"
