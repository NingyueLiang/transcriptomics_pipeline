#!/usr/bin/env python3
"""
Entry point for the transcriptomics analysis pipeline.

Run this script from the repository root inside the `ma4bio` conda environment.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Ensure repo root is in path so 'scripts' package can be imported
REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.mas import TranscriptomicsPipelineMAS
from src.console_utils import Colors, print_warning, print_error


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Python agent for orchestrating the R transcriptomics pipeline."
    )
    parser.add_argument(
        "--experiment-name",
        default="xen_tran_2024_12",
        help="Experiment identifier used for results folders.",
    )
    parser.add_argument(
        "--assay-file",
        default="data/PLX073248_assay_data.csv",
        help="Path to the assay counts CSV (relative to repo root).",
    )
    parser.add_argument(
        "--sample-file",
        default="data/PLX073248_sample_data.csv",
        help="Path to the sample metadata CSV (relative to repo root).",
    )
    parser.add_argument(
        "--xenbase-file",
        default="data/Xenbase_Xenopus_Orthology_Predictions.txt",
        help="Path to the Xenbase orthology file.",
    )
    parser.add_argument(
        "--hcop-file",
        default="data/HCOP_Xenopus_Orthology_Predictions.txt",
        help="Path to the HCOP orthology file.",
    )
    parser.add_argument(
        "--suffix-mode",
        default="all",
        help="Suffix handling mode forwarded to the R mapper.",
    )
    parser.add_argument(
        "--collections",
        default="H,C2,C5",
        help="Comma-separated MSigDB collections to use for GSEA.",
    )
    parser.add_argument(
        "--reference-condition",
        default="Vehicle",
        help="Reference condition supplied to DESeq2.",
    )
    parser.add_argument(
        "--reference-timepoint",
        default="T",
        help="Reference timepoint supplied to DESeq2.",
    )
    parser.add_argument(
        "--group-cols",
        default="condition,timepoint,dose,drug,replicate",
        help="Comma-separated metadata columns for grouping.",
    )
    parser.add_argument(
        "--comparison-cols",
        default="condition,timepoint",
        help="Comma-separated columns for DESeq2 comparison.",
    )
    parser.add_argument(
        "--skip-mapping",
        action="store_true",
        help="Skip the gene mapping step.",
    )
    parser.add_argument(
        "--skip-gsea",
        action="store_true",
        help="Skip the GSEA step.",
    )
    parser.add_argument(
        "--summary-path",
        default=None,
        help="Override the default summary JSON path.",
    )
    parser.add_argument(
        "--report-path",
        default=None,
        help="Where to write the LLM report (defaults to results/<experiment>/llm_report.md).",
    )
    parser.add_argument(
        "--skip-r-setup",
        action="store_true",
        help="Skip the R environment bootstrap step.",
    )
    parser.add_argument(
        "--skip-llm",
        action="store_true",
        help="Skip the LLM report generation phase.",
    )
    parser.add_argument(
        "--user-query",
        default="Provide a comprehensive analysis of the transcriptomics data, highlighting key findings and biological insights.",
        help="User query that guides report generation and verification.",
    )
    parser.add_argument(
        "--num-reports",
        type=int,
        default=1,
        help="Number of LLM reports to generate.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    
    # Configure logging to be less verbose (we have our own pretty output)
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format=f'{Colors.DIM}%(levelname)s: %(message)s{Colors.RESET}'
    )
    
    repo_root = REPO_ROOT

    if args.summary_path:
        summary_path = Path(args.summary_path)
        if not summary_path.is_absolute():
            summary_path = repo_root / summary_path
    else:
        summary_path = repo_root / "results" / args.experiment_name / "pipeline_summary.json"

    if args.report_path:
        report_path = Path(args.report_path)
        if not report_path.is_absolute():
            report_path = repo_root / report_path
    else:
        report_path = repo_root / "results" / args.experiment_name / "llm_report.md"

    try:
        mas = TranscriptomicsPipelineMAS(
            repo_root=repo_root,
            experiment_name=args.experiment_name,
            assay_file=repo_root / args.assay_file,
            sample_file=repo_root / args.sample_file,
            xenbase_file=repo_root / args.xenbase_file,
            hcop_file=repo_root / args.hcop_file,
            suffix_mode=args.suffix_mode,
            collections=args.collections,
            reference_condition=args.reference_condition,
            reference_timepoint=args.reference_timepoint,
            summary_path=summary_path,
            report_path=report_path,
            skip_r_setup=args.skip_r_setup,
            skip_llm=args.skip_llm,
            group_cols=args.group_cols,
            comparison_cols=args.comparison_cols,
            skip_mapping=args.skip_mapping,
            skip_gsea=args.skip_gsea,
        )
        mas.run(user_query=args.user_query, num_reports=args.num_reports)
    except KeyboardInterrupt:
        print_warning("Pipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        print_error(f"Pipeline failed: {e}")
        logging.exception("Pipeline error details:")
        sys.exit(1)


if __name__ == "__main__":
    main()

