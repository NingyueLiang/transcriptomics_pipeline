#!/usr/bin/env python3
"""
Pipeline agent that orchestrates the R-based transcriptomics analysis and
generates an LLM-backed report. Run this script from the repository root
inside the `ma4bio` conda environment so that Python dependencies (litellm)
and CLI tools are available.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import textwrap
from pathlib import Path
from typing import Any, Dict

from litellm import completion
from dotenv import load_dotenv
load_dotenv()

class TranscriptomicsPipelineAgent:
    """Thin orchestrator for the R pipeline and LLM-based report generation."""

    def __init__(
        self,
        repo_root: Path,
        experiment_name: str,
        assay_file: Path,
        sample_file: Path,
        xenbase_file: Path,
        hcop_file: Path,
        suffix_mode: str,
        collections: str,
        reference_condition: str,
        reference_timepoint: str,
        summary_path: Path,
        report_path: Path,
        skip_r_setup: bool = False,
        skip_llm: bool = False,
    ) -> None:
        self.repo_root = repo_root
        self.experiment_name = experiment_name
        self.assay_file = assay_file
        self.sample_file = sample_file
        self.xenbase_file = xenbase_file
        self.hcop_file = hcop_file
        self.suffix_mode = suffix_mode
        self.collections = collections
        self.reference_condition = reference_condition
        self.reference_timepoint = reference_timepoint
        self.summary_path = summary_path
        self.report_path = report_path
        self.skip_r_setup = skip_r_setup
        self.skip_llm = skip_llm

    def run(self) -> None:
        if not self.skip_r_setup:
            self.setup_r_environment()
        else:
            logging.info("Skipping R environment setup as requested.")

        self.run_r_pipeline()
        summary = self.load_summary()

        if self.skip_llm:
            logging.info("Skipping LLM report generation as requested.")
            return

        report = self.generate_llm_report(summary)
        self.write_report(report)

    def setup_r_environment(self) -> None:
        """Bootstrap R dependencies and directory scaffolding."""
        cmd = [
            "Rscript",
            "scripts/setup_r_env.R",
            "--base-dir",
            str(self.repo_root),
        ]
        logging.info("Setting up R environment via %s", " ".join(cmd))
        self._run_command(cmd)

    def run_r_pipeline(self) -> None:
        """Execute the R pipeline that performs the transcriptomics analysis."""
        cmd = [
            "Rscript",
            "scripts/run_r_pipeline.R",
            "--experiment-name",
            self.experiment_name,
            "--assay-file",
            str(self.assay_file),
            "--sample-file",
            str(self.sample_file),
            "--xenbase-file",
            str(self.xenbase_file),
            "--hcop-file",
            str(self.hcop_file),
            "--suffix-mode",
            self.suffix_mode,
            "--collections",
            self.collections,
            "--reference-condition",
            self.reference_condition,
            "--reference-timepoint",
            self.reference_timepoint,
            "--summary-file",
            str(self.summary_path),
        ]
        logging.info("Running R pipeline via %s", " ".join(cmd))
        self._run_command(cmd)

    def load_summary(self) -> Dict[str, Any]:
        """Load the JSON summary emitted by the R pipeline."""
        if not self.summary_path.exists():
            raise FileNotFoundError(
                f"Pipeline summary not found at {self.summary_path}. "
                "Ensure the R pipeline completed successfully."
            )
        with self.summary_path.open("r", encoding="utf-8") as handle:
            summary = json.load(handle)
        logging.info(
            "Loaded pipeline summary for experiment %s",
            summary.get("experiment", {}).get("name", self.experiment_name),
        )
        return summary

    def generate_llm_report(self, summary: Dict[str, Any]) -> str:
        """Call the Together-hosted Qwen model via litellm to draft a report."""
        prompt = textwrap.dedent(
            f"""
            You are an experienced computational biologist. Review the following
            transcriptomics analysis summary (JSON) and produce a concise
            narrative report covering:
            1. Experimental context and inputs.
            2. Key differential expression findings (mention contrasts, counts,
               and notable genes).
            3. Ortholog mapping coverage/quality.
            4. Major GSEA themes and their biological implications.
            5. Recommended follow-up analyses or validations.

            Respond in markdown with short sections and bullet lists where
            helpful. Keep the tone factual and avoid speculation beyond the
            data provided.

            Summary JSON:
            {json.dumps(summary, indent=2)}
            """
        ).strip()

        messages = [{"role": "user", "content": prompt}]
        logging.info("Requesting LLM report from together_ai/Qwen/Qwen3-Next-80B-A3B-Instruct")
        response = completion(
            model="together_ai/Qwen/Qwen3-Next-80B-A3B-Instruct",
            messages=messages,
        )
        content = response["choices"][0]["message"]["content"]
        logging.info("Received LLM response (%d characters).", len(content))
        return content

    def write_report(self, report: str) -> None:
        """Persist the LLM report to disk."""
        self.report_path.parent.mkdir(parents=True, exist_ok=True)
        with self.report_path.open("w", encoding="utf-8") as handle:
            handle.write(report)
        logging.info("LLM report written to %s", self.report_path)

    def _run_command(self, cmd: list[str]) -> None:
        """Execute a command in the repository root with error handling."""
        env = os.environ.copy()
        subprocess.run(
            cmd,
            cwd=self.repo_root,
            env=env,
            check=True,
        )


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
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))
    repo_root = Path(__file__).resolve().parents[1]

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

    agent = TranscriptomicsPipelineAgent(
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
    )
    agent.run()


if __name__ == "__main__":
    main()
