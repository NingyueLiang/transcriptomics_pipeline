#!/usr/bin/env python3
"""
Transcriptomics Pipeline Multi-Agent System (MAS).

Orchestrates the R-based transcriptomics analysis and generates an LLM-backed report.
"""

from __future__ import annotations

import json
import logging
import os
import re
import subprocess
import textwrap
import time
from pathlib import Path
from typing import Any, Dict, List

from dotenv import load_dotenv

from .hvd_api import call_hvd_api
from .console_utils import (
    print_banner,
    print_section,
    print_step,
    print_substep,
    print_info,
    print_error,
    print_config,
    print_timing,
    print_completion_banner,
    print_results_summary,
)

load_dotenv()


class TranscriptomicsPipelineMAS:
    """Thin orchestrator for the R pipeline and LLM-based report generation."""

    TOTAL_STEPS = 5  # Total number of main pipeline steps

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
        group_cols: str = "condition,timepoint,dose,drug,replicate",
        comparison_cols: str = "condition,timepoint",
        skip_mapping: bool = False,
        skip_gsea: bool = False,
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
        self.group_cols = group_cols
        self.comparison_cols = comparison_cols
        self.skip_mapping = skip_mapping
        self.skip_gsea = skip_gsea
        self.pipeline_start_time = None

    def run(self, user_query: str, num_reports: int = 1) -> None:
        """
        Run the transcriptomics pipeline.
        
        Args:
            user_query: The user's query/request that guides report generation and verification.
            num_reports: Number of LLM reports to generate.
        """
        self.pipeline_start_time = time.time()
        self.user_query = user_query
        self.num_reports = num_reports
        
        # Print startup banner
        print_banner()
        
        # Print configuration overview
        print_section("Pipeline Configuration", "⚙️")
        print_config({
            "Experiment": self.experiment_name,
            "Assay File": self.assay_file.name,
            "Sample File": self.sample_file.name,
            "Reference Condition": self.reference_condition,
            "Reference Timepoint": self.reference_timepoint,
            "GSEA Collections": self.collections,
            "Summary Output": str(self.summary_path),
            "Report Output": str(self.report_path),
            "User Query": user_query[:80] + "..." if len(user_query) > 80 else user_query,
            "Number of Reports": str(num_reports),
        })
        
        print_section("Pipeline Execution", "🚀")
        
        #! Step 1: R Environment Setup
        if not self.skip_r_setup:
            print_step(1, self.TOTAL_STEPS, "Setting up R environment")
            self.setup_r_environment()
        else:
            print_step(1, self.TOTAL_STEPS, "R environment setup (skipped)")
            print_substep("User requested --skip-r-setup", "skip")

        #! Step 2: R Pipeline
        print_step(2, self.TOTAL_STEPS, "Running R transcriptomics pipeline")
        self.run_r_pipeline()
        
        #! Step 3: Load Summary
        print_step(3, self.TOTAL_STEPS, "Loading pipeline summary")
        summary = self.load_summary()

        #! Step 4: LLM Report Generation
        reports = []
        if self.skip_llm:
            print_step(4, self.TOTAL_STEPS, "LLM report generation (skipped)")
            print_substep("User requested --skip-llm", "skip")
        else:
            print_step(4, self.TOTAL_STEPS, f"Generating {num_reports} LLM-powered report(s)")
            reports = self.generate_llm_reports(summary, user_query, num_reports)
            self.write_reports(reports)
        
        #! Step 5: Result Verification
        if self.skip_llm or not reports:
            print_step(5, self.TOTAL_STEPS, "Report verification (skipped)")
            print_substep("No reports to verify", "skip")
        else:
            print_step(5, self.TOTAL_STEPS, "Verifying generated reports")
            self.verify_and_save_reports(reports, user_query)
        
        # Final summary
        self._print_completion_summary()

    def setup_r_environment(self) -> None:
        """Bootstrap R dependencies and directory scaffolding."""
        start_time = time.time()
        
        print_substep("Installing R package dependencies...", "running")
        print_substep("Creating directory structure...", "running")
        
        cmd = [
            "Rscript",
            "src/setup_r_env.R",
            "--base-dir",
            str(self.repo_root),
        ]
        logging.info("Setting up R environment via %s", " ".join(cmd))
        self._run_command(cmd)
        
        print_substep("R environment ready", "done")
        print_timing(start_time, "R setup")

    def run_r_pipeline(self) -> None:
        """Execute the R pipeline that performs the transcriptomics analysis."""
        start_time = time.time()
        
        print_substep("Importing and validating input data...", "running")
        print_substep("Cleaning data for differential expression...", "running")
        print_substep("Running DESeq2 analysis...", "running")
        print_substep("Mapping Xenopus genes → Human orthologs...", "running")
        print_substep("Performing GSEA pathway analysis...", "running")
        
        cmd = [
            "Rscript",
            "src/run_r_pipeline.R",
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
            "--group-cols",
            self.group_cols,
            "--comparison-cols",
            self.comparison_cols,
            "--summary-file",
            str(self.summary_path),
        ]
        if self.skip_mapping:
            cmd.append("--skip-mapping")
        if self.skip_gsea:
            cmd.append("--skip-gsea")
        logging.info("Running R pipeline via %s", " ".join(cmd))
        self._run_command(cmd)
        
        print_substep("R pipeline completed successfully", "done")
        print_timing(start_time, "R pipeline")

    def load_summary(self) -> Dict[str, Any]:
        """Load the JSON summary emitted by the R pipeline."""
        print_substep(f"Reading summary from {self.summary_path.name}...", "running")
        
        if not self.summary_path.exists():
            print_error(f"Summary file not found: {self.summary_path}")
            raise FileNotFoundError(
                f"Pipeline summary not found at {self.summary_path}. "
                "Ensure the R pipeline completed successfully."
            )
        with self.summary_path.open("r", encoding="utf-8") as handle:
            summary = json.load(handle)
        
        # Print summary highlights
        de_results = summary.get("differential_expression", [])
        mapping = summary.get("mapping", {})
        gsea = summary.get("gsea", {})

        print_substep("Summary loaded", "done")
        print_info("DE Contrasts", str(len(de_results)))
        print_info("Mapped Genes", str(mapping.get("mapped_genes", "N/A") if isinstance(mapping, dict) else "N/A"))
        gsea_pathways = gsea.get("significant_pathways", "N/A") if isinstance(gsea, dict) else "Skipped"
        print_info("GSEA Pathways", str(gsea_pathways))
        
        logging.info(
            "Loaded pipeline summary for experiment %s",
            summary.get("experiment", {}).get("name", self.experiment_name),
        )
        return summary

    def generate_llm_reports(
        self, summary: Dict[str, Any], user_query: str, num_reports: int
    ) -> list[str]:
        """
        Call the Harvard-hosted GPT-5-mini model via HVD API to draft multiple reports.
        
        Args:
            summary: The pipeline summary data.
            user_query: The user's query/request that guides report generation.
            num_reports: Number of reports to generate.
            
        Returns:
            A list of generated report strings.
        """
        start_time = time.time()
        reports = []
        
        print_substep("Preparing prompt for LLM...", "running")
        
        prompt_template = textwrap.dedent(
            f"""
            You are an experienced computational biologist. A user has requested
            the following analysis:
            
            === USER QUERY ===
            {user_query}
            === END USER QUERY ===
            
            Review the following transcriptomics analysis summary (JSON) and produce 
            a concise narrative report that addresses the user's query. Your report 
            should cover:
            1. Experimental context and inputs.
            2. Key differential expression findings (mention contrasts, counts,
               and notable genes).
            3. Ortholog mapping coverage/quality.
            4. Major GSEA themes and their biological implications.
            5. Recommended follow-up analyses or validations.
            6. Specific answers or insights addressing the user's query.

            Respond in markdown with short sections and bullet lists where
            helpful. Keep the tone factual and avoid speculation beyond the
            data provided. Make sure your report directly addresses the user's
            query.

            Summary JSON:
            {json.dumps(summary, indent=2)}
            """
        ).strip()

        print_info("Model", "gpt-5-mini")
        print_info("User Query", user_query[:60] + "..." if len(user_query) > 60 else user_query)
        print_info("Prompt Length", f"{len(prompt_template):,} chars")
        print_info("Reports to Generate", str(num_reports))
        
        for i in range(num_reports):
            print_substep(f"Generating report {i + 1}/{num_reports}...", "running")
            logging.info(f"Requesting LLM report {i + 1}/{num_reports} from Harvard API (gpt-5-mini)")
            
            content = call_hvd_api(prompt_template, model="gpt-5-mini")
            reports.append(content)
            
            print_substep(f"Report {i + 1} received ({len(content):,} chars)", "done")
            logging.info(f"Received LLM response {i + 1} ({len(content)} characters).")
        
        print_timing(start_time, "LLM generation")
        
        return reports

    def write_reports(self, reports: list[str]) -> None:
        """Persist multiple LLM reports to disk."""
        self.report_path.parent.mkdir(parents=True, exist_ok=True)
        
        for i, report in enumerate(reports):
            if len(reports) == 1:
                report_file = self.report_path
            else:
                # For multiple reports, add index to filename
                stem = self.report_path.stem
                suffix = self.report_path.suffix
                report_file = self.report_path.parent / f"{stem}_{i + 1}{suffix}"
            
            print_substep(f"Writing report {i + 1} to {report_file.name}...", "running")
            
            with report_file.open("w", encoding="utf-8") as handle:
                handle.write(report)
            
            print_substep(f"Report {i + 1} saved", "done")
            print_info("Location", str(report_file))
            
            logging.info("LLM report %d written to %s", i + 1, report_file)

    def verify_llm_report(self, report: str, user_query: str) -> Dict[str, Any]:
        """
        Verify if a generated report satisfies the user query.
        
        Args:
            report: The generated report content.
            user_query: The original user query to verify against.
            
        Returns:
            A dictionary containing verification results.
        """
        verification_prompt = textwrap.dedent(
            f"""
            You are a quality assurance specialist for scientific reports.
            
            A user requested the following analysis:
            
            === USER QUERY ===
            {user_query}
            === END USER QUERY ===
            
            The following report was generated in response:
            
            === GENERATED REPORT ===
            {report}
            === END GENERATED REPORT ===
            
            Please evaluate whether the report adequately addresses the user's query.
            
            Respond in the following JSON format ONLY (no additional text):
            {{
                "satisfies_query": true/false,
                "confidence_score": 0.0-1.0,
                "addressed_aspects": ["list of aspects from the query that were addressed"],
                "missing_aspects": ["list of aspects from the query that were NOT addressed"],
                "quality_score": 0.0-1.0,
                "strengths": ["list of report strengths"],
                "weaknesses": ["list of report weaknesses"],
                "summary": "Brief summary of the verification result"
            }}
            """
        ).strip()
        
        logging.info("Verifying report against user query via LLM")
        response = call_hvd_api(verification_prompt, model="gpt-5-mini")
        
        # Parse the JSON response
        try:
            # Try to extract JSON from the response (in case there's extra text)
            json_match = re.search(r'\{[\s\S]*\}', response)
            if json_match:
                verification_result = json.loads(json_match.group())
            else:
                verification_result = json.loads(response)
        except json.JSONDecodeError as e:
            logging.warning(f"Failed to parse verification response as JSON: {e}")
            verification_result = {
                "satisfies_query": None,
                "confidence_score": 0.0,
                "addressed_aspects": [],
                "missing_aspects": [],
                "quality_score": 0.0,
                "strengths": [],
                "weaknesses": [],
                "summary": f"Verification parsing failed: {str(e)}",
                "raw_response": response,
            }
        
        return verification_result

    def verify_and_save_reports(self, reports: list[str], user_query: str) -> None:
        """
        Verify all generated reports and save results alongside reports.
        
        Args:
            reports: List of generated report strings.
            user_query: The original user query.
        """
        start_time = time.time()
        self.report_path.parent.mkdir(parents=True, exist_ok=True)
        
        for i, report in enumerate(reports):
            print_substep(f"Verifying report {i + 1}/{len(reports)}...", "running")
            
            # Verify the report
            verification_result = self.verify_llm_report(report, user_query)
            
            # Determine file paths
            if len(reports) == 1:
                json_file = self.report_path.with_suffix(".json")
            else:
                stem = self.report_path.stem
                json_file = self.report_path.parent / f"{stem}_{i + 1}_verified.json"
            
            # Create combined output with report and verification
            combined_output = {
                "user_query": user_query,
                "report_index": i + 1,
                "report_content": report,
                "verification": verification_result,
            }
            
            # Save to JSON file
            with json_file.open("w", encoding="utf-8") as handle:
                json.dump(combined_output, handle, indent=2)
            
            # Print verification summary
            satisfies = verification_result.get("satisfies_query", "unknown")
            confidence = verification_result.get("confidence_score", 0.0)
            quality = verification_result.get("quality_score", 0.0)
            
            status_icon = "✓" if satisfies else "✗" if satisfies is False else "?"
            print_substep(f"Report {i + 1} verification: {status_icon} (conf: {confidence:.2f}, qual: {quality:.2f})", "done")
            print_info(f"Verification saved", str(json_file))
            
            logging.info(
                "Report %d verification: satisfies=%s, confidence=%.2f, quality=%.2f",
                i + 1, satisfies, confidence, quality
            )
        
        print_timing(start_time, "Verification")
    
    def _print_completion_summary(self) -> None:
        """Print a beautiful completion summary."""
        total_time = time.time() - self.pipeline_start_time
        print_completion_banner()
        print_results_summary(
            experiment_name=self.experiment_name,
            summary_path=str(self.summary_path),
            report_path=str(self.report_path),
            total_time=total_time,
        )

    def _run_command(self, cmd: list[str]) -> None:
        """Execute a command in the repository root with error handling."""
        env = os.environ.copy()
        subprocess.run(
            cmd,
            cwd=self.repo_root,
            env=env,
            check=True,
        )

