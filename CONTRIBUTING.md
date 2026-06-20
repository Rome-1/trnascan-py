# Contributing to trnascan-py

Thanks for your interest! trnascan-py is a thin, well-tested Python interface to
Infernal-driven tRNA detection. Contributions that keep it correct, fast, and
honest about its scope are very welcome.

## Development setup

```bash
# 1. Install the external engines (provide cmsearch + the tRNA covariance models)
sudo apt-get install -y infernal trnascan-se     # or: conda install -c bioconda infernal trnascan-se

# 2. Install the package with dev extras
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

## Quality gates (run before every PR)

```bash
ruff check .            # lint + import order
mypy                    # static types (strict)
pytest                  # unit tests — no external tools needed
pytest -m differential  # differential tests vs reference tRNAscan-SE (needs the engines)
```

All four must pass. CI runs the same gates (see `.github/workflows/ci.yml`),
including a job that installs the engines for the differential tests.

## Conventions

- **Test-first & differential.** Correctness is defined by agreement with
  reference tRNAscan-SE 2.0. New behavior should come with a unit test and, where
  it affects calls, a differential assertion. Unit tests must not require the
  external engines — gate engine-dependent tests with the `requires_infernal` /
  `requires_trnascan` / `differential` markers.
- **Typed and documented.** Public functions are fully type-annotated (mypy
  strict) and have docstrings. Keep line length ≤ 100.
- **Thin v1.** v1 orchestrates Infernal; it does not re-implement covariance-model
  DP. Native/accelerated kernels belong to v2 — please discuss in an issue first.
- **Honest scope.** If a change has limitations, document them in the README
  "Scope & known limitations" section.

## Reporting issues

Include: input (a small FASTA if possible), the exact command, `trnascan-py
version` output, and Infernal / tRNAscan-SE versions.

## License of contributions

trnascan-py is licensed **GPL-3.0-or-later**. By contributing you agree your
contributions are licensed under the same terms. Logic derived from the reference
tRNAscan-SE 2.0 (itself GPL-3) is permitted; please attribute the source routine
in a comment when you port one.
