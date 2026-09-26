from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .frozen_v1 import VERSION
from .frozen_v1.filtering import SelectionRequiredError, run_filter
from .frozen_v1.wgs import run_wgs
from .mapping import MappingConfig
from .workflow.mapping import run_mapping


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cinch", description="Cinch WGS v1 dependency discovery and Cinch Filter v1 recurrence filtering")
    parser.add_argument("--version", action="version", version=f"Cinch {VERSION}")
    sub = parser.add_subparsers(dest="command", required=True)
    mapping = sub.add_parser("map", help="indexed, parallel, resumable genome mapping and profile generation")
    mapping.add_argument("genomes", nargs="+", type=Path, help="input genome FASTA files")
    mapping.add_argument("-r", "--reference", required=True, type=Path, help="reference CDS FASTA")
    mapping.add_argument("-o", "--output", required=True, type=Path, help="mapping-stage output directory")
    mapping.add_argument("-t", "--threads", type=int, default=1, help="number of genome worker processes")
    mapping.add_argument("--min-identity", type=float, default=.90)
    mapping.add_argument("--detection-coverage", type=float, default=.60)
    mapping.add_argument("--callable-coverage", type=float, default=.95)
    mapping.add_argument("--min-mapq", type=int, default=0)
    wgs = sub.add_parser("wgs", help="build a complete unfiltered PP/PT/TP/TT dependency landscape from genomes")
    wgs.add_argument("genomes", nargs="+", type=Path, help="input genome FASTA files")
    wgs.add_argument("-r", "--reference", required=True, type=Path, help="reference CDS FASTA with unique locus IDs")
    wgs.add_argument("-p", "--prefix", required=True, help="run prefix; default output is PREFIX.cinch")
    wgs.add_argument("-t", "--threads", type=int, default=1, help="requested worker threads (recorded; deterministic output)")
    wgs.add_argument("-o", "--output", type=Path, help="explicit WGS result directory")
    wgs.add_argument("--phenotype", type=Path, help="optional phenotype TSV/CSV; validated and frozen only")
    wgs.add_argument("--min-identity", type=float, default=.90, help="minimum full-length nucleotide identity for a locus call [0.90]")
    wgs.add_argument("--min-informative", type=int, default=20, help="minimum pairwise-complete samples [20]")
    wgs.add_argument("--min-state-count", type=int, default=3, help="minimum marginal state support after rare pooling [3]")
    wgs.add_argument("--min-order-observations", type=int, default=5, help="minimum same-contig observations for order distance [5]")
    filt = sub.add_parser("filter", help="filter one complete Cinch WGS result directory into network-ready candidate dependencies")
    filt.add_argument("--wgs_results", required=True, type=Path, help="path to a complete Cinch WGS v1 result directory")
    filt.add_argument("--hc", type=int, help="dataset-specific HC allele-distance level; HCx is not x clusters")
    filt.add_argument("--order-threshold", type=float, help="dataset-specific long-range gene-order threshold")
    filt.add_argument("--minimum-populations", type=int, default=3, help="minimum HC blocks carrying the same driver [3, frozen v1]")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "map":
            config = MappingConfig(
                minimum_identity=args.min_identity,
                detection_coverage=args.detection_coverage,
                callable_coverage=args.callable_coverage,
                minimum_mapq=args.min_mapq,
            )
            output = run_mapping(args.reference, args.genomes, args.output, args.threads, config)
        elif args.command == "wgs":
            output = run_wgs(args.reference, args.genomes, args.prefix, args.threads, args.output,
                             args.phenotype, args.min_identity, args.min_informative,
                             args.min_state_count, args.min_order_observations)
        else:
            output = run_filter(args.wgs_results, args.hc, args.order_threshold, args.minimum_populations)
        print(f"RESULT_DIR={output}")
        return 0
    except SelectionRequiredError as error:
        print(f"SELECTION_REQUIRED: {error}", file=sys.stderr)
        return 3
    except Exception as error:
        print(f"ERROR: {type(error).__name__}: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
