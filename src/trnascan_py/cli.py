"""Command-line interface for trnascan-py.

Mirrors the spirit of the reference ``tRNAscan-SE`` CLI for the v1 feature set:
a ``scan`` subcommand that runs the Infernal-backed pipeline, plus a ``version``
subcommand reporting the versions of the underlying engines.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence

from trnascan_py import __version__
from trnascan_py.infernal import CmsearchError, InfernalNotFoundError
from trnascan_py.models import TRNAHit
from trnascan_py.models_registry import ModelsNotFoundError
from trnascan_py.pipeline import DEFAULT_CM_CUTOFF, scan


def _hit_to_dict(hit: TRNAHit) -> dict[str, object]:
    return {
        "seq_id": hit.seq_id,
        "start": hit.start,
        "end": hit.end,
        "strand": hit.strand.value,
        "isotype": hit.isotype,
        "anticodon": hit.anticodon,
        "score": hit.score,
        "intron_start": hit.intron_start,
        "intron_end": hit.intron_end,
        "note": hit.note,
    }


def _print_table(hits: list[TRNAHit]) -> None:
    header = [
        "Name", "Begin", "End", "Strand", "Type", "Codon",
        "IntronBegin", "IntronEnd", "Score", "Note",
    ]
    print("\t".join(header))
    for h in hits:
        print(
            "\t".join(
                [
                    h.seq_id,
                    str(h.start),
                    str(h.end),
                    h.strand.value,
                    h.isotype,
                    h.anticodon or "-",
                    str(h.intron_start or 0),
                    str(h.intron_end or 0),
                    f"{h.score:.1f}",
                    h.note or "",
                ]
            )
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="trnascan-py",
        description="Python interface to tRNAscan-SE (v1: Infernal-backed).",
    )
    parser.add_argument("--version", action="version", version=f"trnascan-py {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    scan_p = sub.add_parser("scan", help="scan a FASTA file for tRNA genes")
    scan_p.add_argument("fasta", help="input sequence FASTA file")
    scan_p.add_argument(
        "-d",
        "--domain",
        choices=("euk", "bact", "arch", "mito-vert", "mito-mammal"),
        default="euk",
        help="model set: euk/bact/arch or mitochondrial mito-vert/mito-mammal (default: euk)",
    )
    scan_p.add_argument(
        "--cm-cutoff",
        type=float,
        default=DEFAULT_CM_CUTOFF,
        help=f"minimum bit score to report a hit (default: {DEFAULT_CM_CUTOFF})",
    )
    scan_p.add_argument(
        "--models-dir", default=None, help="override the covariance-model directory"
    )
    scan_p.add_argument(
        "--no-isotype", action="store_true", help="skip the isotype classification pass"
    )
    scan_p.add_argument(
        "--isotype-method",
        choices=("hybrid", "anticodon", "cm"),
        default="hybrid",
        help="isotype assignment: hybrid (default; anticodon + CM refine), "
        "anticodon (genetic code only, fastest), or cm (CM scan every locus)",
    )
    scan_p.add_argument("--cpu", type=int, default=None, help="number of Infernal worker threads")
    scan_p.add_argument("--json", action="store_true", help="emit JSON instead of a table")

    sub.add_parser("version", help="show versions of trnascan-py and underlying engines")
    return parser


def _cmd_version() -> int:
    from trnascan_py.infernal import InfernalNotFoundError, cmsearch_version
    from trnascan_py.oracle import TrnascanNotFoundError, trnascan_version

    print(f"trnascan-py {__version__}")
    try:
        print(f"infernal: {cmsearch_version()}")
    except InfernalNotFoundError as exc:
        print(f"infernal: not found ({exc})")
    try:
        print(f"reference tRNAscan-SE: {trnascan_version()}")
    except TrnascanNotFoundError as exc:
        print(f"reference tRNAscan-SE: not found ({exc})")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "version":
        return _cmd_version()

    if args.command == "scan":
        try:
            hits = scan(
                args.fasta,
                domain=args.domain,
                models_dir=args.models_dir,
                cm_cutoff=args.cm_cutoff,
                classify_isotype=not args.no_isotype,
                isotype_method=args.isotype_method,
                cpu=args.cpu,
            )
        except (
            FileNotFoundError,
            ModelsNotFoundError,
            InfernalNotFoundError,
            CmsearchError,
        ) as exc:
            print(f"trnascan-py: error: {exc}", file=sys.stderr)
            return 1
        if args.json:
            json.dump([_hit_to_dict(h) for h in hits], sys.stdout, indent=2)
            sys.stdout.write("\n")
        else:
            _print_table(hits)
        return 0

    parser.error(f"unknown command: {args.command}")
    return 2  # unreachable; parser.error exits


if __name__ == "__main__":
    raise SystemExit(main())
