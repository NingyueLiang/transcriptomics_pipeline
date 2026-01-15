#!/usr/bin/env python3
"""
Entry point for the transcriptomics analysis pipeline.

Supports two agent modes via --agent-type:
- qa: Question-Answering agent (input: question + database → output: report)
- se: Self-Exploration agent (input: database → output: research questions + results)

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

from src.question_answering_agent.agent import QuestionAnsweringAgent
from src.self_exploration_agent.agent import SelfExplorationAgent
from src.utils.console_utils import Colors, print_warning, print_error


def parse_args() -> argparse.Namespace:
    """Create the argument parser."""
    parser = argparse.ArgumentParser(
        description="Transcriptomics Multi-Agent System",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run Question-Answering agent (default)
  python main.py --agent-type qa --experiment-name scenario_11a \\
    --assay-file scenarios/scenario_11a/scenario_11a_assay_data.csv \\
    --sample-file scenarios/scenario_11a/scenario_11a_sample_data.csv \\
    --user-query "What genes are differentially expressed?"

  # Run Self-Exploration agent
  python main.py --agent-type se --experiment-name scenario_11a \\
    --assay-file scenarios/scenario_11a/scenario_11a_assay_data.csv \\
    --sample-file scenarios/scenario_11a/scenario_11a_sample_data.csv \\
    --num-rqs 10
        """,
    )

    # Agent selection
    parser.add_argument(
        "--agent-type",
        choices=["qa", "se"],
        default="qa",
        help="Agent type: 'qa' (Question-Answering) or 'se' (Self-Exploration). Default: qa",
    )

    # Common arguments
    parser.add_argument(
        "--experiment-name",
        default="experiment",
        help="Experiment identifier used for results folders.",
    )
    parser.add_argument(
        "--assay-file",
        required=True,
        help="Path to the assay counts CSV (relative to repo root).",
    )
    parser.add_argument(
        "--sample-file",
        required=True,
        help="Path to the sample metadata CSV (relative to repo root).",
    )
    parser.add_argument(
        "--reference-condition",
        default="Control",
        help="Reference condition for DESeq2.",
    )
    parser.add_argument(
        "--group-cols",
        default="condition,replicate,batch",
        help="Comma-separated metadata columns for grouping.",
    )
    parser.add_argument(
        "--comparison-cols",
        default="condition",
        help="Comma-separated columns for DESeq2 comparison.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity.",
    )

    # QA agent specific arguments
    qa_group = parser.add_argument_group("QA Agent Options")
    qa_group.add_argument(
        "--xenbase-file",
        default="data/Xenbase_Xenopus_Orthology_Predictions.txt",
        help="Path to the Xenbase orthology file.",
    )
    qa_group.add_argument(
        "--hcop-file",
        default="data/HCOP_Xenopus_Orthology_Predictions.txt",
        help="Path to the HCOP orthology file.",
    )
    qa_group.add_argument(
        "--suffix-mode",
        default="all",
        help="Suffix handling mode forwarded to the R mapper.",
    )
    qa_group.add_argument(
        "--collections",
        default="H,C2,C5",
        help="Comma-separated MSigDB collections to use for GSEA.",
    )
    qa_group.add_argument(
        "--reference-timepoint",
        default="T",
        help="Reference timepoint supplied to DESeq2.",
    )
    qa_group.add_argument(
        "--skip-mapping",
        action="store_true",
        help="Skip the gene mapping step.",
    )
    qa_group.add_argument(
        "--skip-gsea",
        action="store_true",
        help="Skip the GSEA step.",
    )
    qa_group.add_argument(
        "--summary-path",
        default=None,
        help="Override the default summary JSON path.",
    )
    qa_group.add_argument(
        "--report-path",
        default=None,
        help="Where to write the LLM report.",
    )
    qa_group.add_argument(
        "--skip-r-setup",
        action="store_true",
        help="Skip the R environment bootstrap step.",
    )
    qa_group.add_argument(
        "--skip-llm",
        action="store_true",
        help="Skip the LLM report generation phase.",
    )
    qa_group.add_argument(
        "--user-query",
        default="Provide a comprehensive analysis of the transcriptomics data, highlighting key findings and biological insights.",
        help="User query that guides report generation.",
    )
    qa_group.add_argument(
        "--num-reports",
        type=int,
        default=1,
        help="Number of LLM reports to generate.",
    )

    # SE agent specific arguments
    se_group = parser.add_argument_group("SE Agent Options")
    se_group.add_argument(
        "--num-rqs",
        type=int,
        default=5,
        help="Number of research questions to generate and analyze.",
    )
    se_group.add_argument(
        "--output-dir",
        default=None,
        help="Directory to save exploration results (defaults to results/<experiment>/exploration/).",
    )

    return parser.parse_args()


def run_qa_agent(args: argparse.Namespace, repo_root: Path) -> None:
    """Run the Question-Answering agent."""
    # QA agent outputs go to results/<experiment>/qa_agent/
    qa_output_dir = repo_root / "results" / args.experiment_name / "qa_agent"

    if args.summary_path:
        summary_path = Path(args.summary_path)
        if not summary_path.is_absolute():
            summary_path = repo_root / summary_path
    else:
        summary_path = qa_output_dir / "pipeline_summary.json"

    if args.report_path:
        report_path = Path(args.report_path)
        if not report_path.is_absolute():
            report_path = repo_root / report_path
    else:
        report_path = qa_output_dir / "report.md"

    mas = QuestionAnsweringAgent(
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


def run_se_agent(args: argparse.Namespace, repo_root: Path) -> None:
    """Run the Self-Exploration agent."""
    # SE agent outputs go to results/<experiment>/se_agent/
    if args.output_dir:
        output_dir = Path(args.output_dir)
        if not output_dir.is_absolute():
            output_dir = repo_root / output_dir
    else:
        output_dir = repo_root / "results" / args.experiment_name / "se_agent"

    agent = SelfExplorationAgent(
        repo_root=repo_root,
        experiment_name=args.experiment_name,
        assay_file=repo_root / args.assay_file,
        sample_file=repo_root / args.sample_file,
        output_dir=output_dir,
        reference_condition=args.reference_condition,
        group_cols=args.group_cols,
        comparison_cols=args.comparison_cols,
    )
    agent.run(num_rqs=args.num_rqs)


def main() -> None:
    args = parse_args()

    # Configure logging
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format=f'{Colors.DIM}%(levelname)s: %(message)s{Colors.RESET}'
    )

    repo_root = REPO_ROOT

    try:
        if args.agent_type == "qa":
            run_qa_agent(args, repo_root)
        elif args.agent_type == "se":
            run_se_agent(args, repo_root)
        else:
            print_error(f"Unknown agent type: {args.agent_type}")
            sys.exit(1)
    except KeyboardInterrupt:
        print_warning("Agent interrupted by user")
        sys.exit(1)
    except Exception as e:
        print_error(f"Agent failed: {e}")
        logging.exception("Agent error details:")
        sys.exit(1)


if __name__ == "__main__":
    main()
