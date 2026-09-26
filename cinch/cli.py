from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .frozen_v1 import VERSION
from .frozen_v1.filtering import SelectionRequiredError, run_filter
from .frozen_v1.wgs import run_wgs
from .mapping import MappingConfig
from .workflow.mapping import run_mapping
from .workflow.profile import run_profile_conversion
from .workflow.association import run_association
from .workflow.advanced_filter import run_advanced_filter
from .workflow.reporting import run_report


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
    mapping.add_argument("--no-hit-policy", choices=["unresolved", "absence"], default="unresolved",
                         help="interpret zero alignment hits explicitly; absence enables PP but assumes adequate assembly/search sensitivity")
    profile = sub.add_parser("profile", help="convert a legacy allele table using an explicit missingness policy")
    profile.add_argument("source", type=Path, help="sample-by-locus TSV/CSV; first column is sample ID")
    profile.add_argument("-o", "--output", required=True, type=Path)
    profile.add_argument("--missing-policy", required=True, choices=["absence", "unresolved"])
    profile.add_argument("--separator", help="explicit field separator; default auto-detect")
    associate = sub.add_parser("associate", help="restartable blockwise PP/PT/TP/TT association stage")
    associate.add_argument("--profile", required=True, type=Path, help="validated PROFILE_V2/UNIFIED_STATES npz")
    associate.add_argument("--coordinates", required=True, type=Path, help="mapping coordinates TSV or Parquet")
    associate.add_argument("-o", "--output", required=True, type=Path)
    associate.add_argument("--min-informative", type=int, default=20)
    associate.add_argument("--min-state-count", type=int, default=3)
    associate.add_argument("--min-order-observations", type=int, default=5)
    associate.add_argument("--block-pairs", type=int, default=100_000)
    advanced = sub.add_parser("advanced-filter", help="SHC permutation, BH, and ARACNE filtering for PP edges")
    advanced.add_argument("--association", required=True, type=Path, help="association stage directory containing blocks/")
    advanced.add_argument("--profile", required=True, type=Path)
    advanced.add_argument("--shc", required=True, type=Path, help="table with sample_id and shc")
    advanced.add_argument("--weights", type=Path, help="optional table with sample_id and weight")
    advanced.add_argument("-o", "--output", required=True, type=Path)
    advanced.add_argument("--permutations", type=int, default=999)
    advanced.add_argument("--min-cluster-size", type=int, default=4)
    advanced.add_argument("--min-informative-clusters", type=int, default=2)
    advanced.add_argument("--alpha", type=float, default=.05)
    advanced.add_argument("--seed", type=int, default=20260926)
    report = sub.add_parser("report", help="render report-only outputs from a staged advanced-filter result")
    report.add_argument("--filter-results", required=True, type=Path)
    report.add_argument("-o", "--output", required=True, type=Path)
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
                no_hit_policy=args.no_hit_policy,
            )
            output = run_mapping(args.reference, args.genomes, args.output, args.threads, config)
        elif args.command == "profile":
            output = run_profile_conversion(
                args.source, args.output, missing_policy=args.missing_policy, separator=args.separator,
            )
        elif args.command == "associate":
            output = run_association(
                args.profile, args.coordinates, args.output,
                minimum_informative=args.min_informative,
                minimum_state_count=args.min_state_count,
                minimum_order_observations=args.min_order_observations,
                block_pairs=args.block_pairs,
            )
        elif args.command == "advanced-filter":
            output = run_advanced_filter(
                args.association, args.profile, args.shc, args.output,
                weights_path=args.weights, permutations=args.permutations,
                minimum_cluster_size=args.min_cluster_size,
                minimum_informative_clusters=args.min_informative_clusters,
                alpha=args.alpha, seed=args.seed,
            )
        elif args.command == "report":
            output = run_report(args.filter_results, args.output)
        elif args.command == "wgs":
            output = run_wgs(args.reference, args.genomes, args.prefix, args.threads, args.output,
                             args.phenotype, args.min_identity, args.min_informative,
                             args.min_state_count, args.min_order_observations)
        elif args.command == "filter":
            output = run_filter(args.wgs_results, args.hc, args.order_threshold, args.minimum_populations)
        else:  # pragma: no cover - argparse enforces known commands
            raise ValueError(f"unknown command {args.command}")
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
