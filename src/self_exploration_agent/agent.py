#!/usr/bin/env python3
"""
Self-Exploration Agent for Transcriptomics Data.

Autonomously explores transcriptomics databases, proposes research questions,
and validates hypotheses using R-based analysis tools.
"""

from __future__ import annotations

import json
import logging
import re
import subprocess
import textwrap
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

from dotenv import load_dotenv

from ..utils.hvd_api import call_hvd_api
from ..utils.console_utils import (
    print_banner,
    print_section,
    print_step,
    print_substep,
    print_info,
    print_error,
    print_config,
    print_timing,
    print_completion_banner,
)

load_dotenv()


@dataclass
class ResearchQuestion:
    """A research question with its analysis results."""
    question_id: int
    question: str
    hypothesis: str
    analysis_type: str  # 'differential_expression', 'pca', 'gsea', 'correlation', etc.
    analysis_params: Dict[str, Any] = field(default_factory=dict)
    result: Optional[Dict[str, Any]] = None
    interpretation: Optional[str] = None
    validated: bool = False
    validation_score: float = 0.0
    timestamp: str = ""


class SelfExplorationAgent:
    """
    Agent that autonomously explores transcriptomics data.

    Workflow:
    1. Load and explore the database
    2. Generate research questions based on data characteristics
    3. Execute appropriate R analyses
    4. Validate and interpret results
    5. Loop to generate more questions
    """

    def __init__(
        self,
        repo_root: Path,
        experiment_name: str,
        assay_file: Path,
        sample_file: Path,
        output_dir: Path,
        reference_condition: str = "Control",
        group_cols: str = "condition,replicate,batch",
        comparison_cols: str = "condition",
    ) -> None:
        self.repo_root = repo_root
        self.experiment_name = experiment_name
        self.assay_file = assay_file
        self.sample_file = sample_file
        self.output_dir = output_dir
        self.reference_condition = reference_condition
        self.group_cols = group_cols
        self.comparison_cols = comparison_cols

        # State tracking
        self.data_summary: Dict[str, Any] = {}
        self.research_questions: List[ResearchQuestion] = []
        self.exploration_history: List[Dict[str, Any]] = []
        self.question_counter = 0
        self.start_time: Optional[float] = None

    def run(self, num_rqs: int = 5) -> List[ResearchQuestion]:
        """
        Run the self-exploration loop.

        Args:
            num_rqs: Number of research questions to generate and analyze.

        Returns:
            List of research questions with results.
        """
        self.start_time = time.time()

        print_banner()
        print_section("Self-Exploration Agent", "🔬")
        print_config({
            "Experiment": self.experiment_name,
            "Assay File": self.assay_file.name,
            "Sample File": self.sample_file.name,
            "Target Questions": str(num_rqs),
            "Output Directory": str(self.output_dir),
        })

        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Phase 1: Data Exploration
        print_section("Phase 1: Data Exploration", "📊")
        print_step(1, 3, "Exploring database structure")
        self.explore_database()

        # Phase 2: Generate all research questions at once
        print_section("Phase 2: Hypothesis Generation & Testing", "🧪")
        print_step(2, 3, f"Generating {num_rqs} research questions")
        questions = self.generate_research_questions(num_questions=num_rqs)

        # Phase 3: Execute and validate each question
        for i, rq in enumerate(questions, 1):
            print_substep(f"[{i}/{len(questions)}] {rq.question[:50]}...", "running")

            # Execute analysis
            self.execute_analysis(rq)

            # Validate results
            self.validate_result(rq)

            # Store result
            self.research_questions.append(rq)

            status = "✓" if rq.validated else "○"
            print_substep(f"{status} Q{i} complete (score: {rq.validation_score:.2f})", "done")

        # Phase 4: Save results
        print_section("Phase 3: Saving Results", "💾")
        print_step(3, 3, "Saving exploration results")
        self.save_results()

        # Print summary
        self._print_summary()

        return self.research_questions

    def explore_database(self) -> None:
        """
        Explore the database to understand its structure and characteristics.
        """
        start_time = time.time()

        print_substep("Reading assay data...", "running")
        assay_info = self._analyze_assay_file()

        print_substep("Reading sample metadata...", "running")
        sample_info = self._analyze_sample_file()

        # Combine into data summary
        self.data_summary = {
            "experiment_name": self.experiment_name,
            "assay": assay_info,
            "samples": sample_info,
            "exploration_timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        }

        # Use LLM to generate initial data insights
        print_substep("Generating data insights via LLM...", "running")
        self.data_summary["llm_insights"] = self._generate_data_insights()

        print_substep("Database exploration complete", "done")
        print_info("Genes", str(assay_info.get("num_genes", "N/A")))
        print_info("Samples", str(sample_info.get("num_samples", "N/A")))
        print_info("Conditions", str(sample_info.get("conditions", "N/A")))
        print_timing(start_time, "Data exploration")

    def _analyze_assay_file(self) -> Dict[str, Any]:
        """Analyze the assay data file."""
        import csv

        with open(self.assay_file, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)

            # Count genes
            gene_count = sum(1 for _ in reader)

        # Get sample columns (all except first which is gene_symbol)
        sample_cols = header[1:]

        return {
            "num_genes": gene_count,
            "sample_columns": sample_cols,
            "file_path": str(self.assay_file),
        }

    def _analyze_sample_file(self) -> Dict[str, Any]:
        """Analyze the sample metadata file."""
        import csv

        with open(self.sample_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        if not rows:
            return {"num_samples": 0, "columns": [], "conditions": []}

        columns = list(rows[0].keys())

        # Extract unique values for key columns
        conditions = list(set(r.get("condition", "") for r in rows if r.get("condition")))
        batches = list(set(r.get("batch", "") for r in rows if r.get("batch")))

        return {
            "num_samples": len(rows),
            "columns": columns,
            "conditions": conditions,
            "batches": batches,
            "sample_data": rows,
            "file_path": str(self.sample_file),
        }

    def _generate_data_insights(self) -> str:
        """Use LLM to generate insights about the data."""
        prompt = textwrap.dedent(f"""
            You are a computational biologist exploring a transcriptomics dataset.

            Dataset Summary:
            - Experiment: {self.experiment_name}
            - Number of genes: {self.data_summary.get('assay', {}).get('num_genes', 'Unknown')}
            - Number of samples: {self.data_summary.get('samples', {}).get('num_samples', 'Unknown')}
            - Conditions: {self.data_summary.get('samples', {}).get('conditions', [])}
            - Batches: {self.data_summary.get('samples', {}).get('batches', [])}
            - Sample columns: {self.data_summary.get('samples', {}).get('columns', [])}

            Provide a brief (2-3 paragraph) summary of:
            1. What this dataset appears to contain
            2. What experimental design seems to be used
            3. What types of biological questions could be addressed

            Be concise and factual.
        """).strip()

        return call_hvd_api(prompt, model="gpt-5-mini")

    def generate_research_questions(self, num_questions: int = 3) -> List[ResearchQuestion]:
        """
        Generate research questions based on the data and exploration history.
        """
        start_time = time.time()

        # Build context from exploration history
        history_context = ""
        if self.research_questions:
            history_context = "\n\nPreviously explored questions:\n"
            for rq in self.research_questions[-5:]:  # Last 5 questions
                history_context += f"- {rq.question} (analysis: {rq.analysis_type})\n"

        prompt = textwrap.dedent(f"""
            You are a computational biologist generating research questions for a transcriptomics dataset.

            Dataset Summary:
            - Experiment: {self.experiment_name}
            - Number of genes: {self.data_summary.get('assay', {}).get('num_genes', 'Unknown')}
            - Conditions: {self.data_summary.get('samples', {}).get('conditions', [])}
            - Batches: {self.data_summary.get('samples', {}).get('batches', [])}
            - Available metadata: {self.data_summary.get('samples', {}).get('columns', [])}

            Data Insights:
            {self.data_summary.get('llm_insights', 'No insights available')}
            {history_context}

            Available analysis types:
            1. differential_expression - Compare gene expression between conditions using DESeq2
            2. pca - Principal component analysis to explore sample clustering
            3. gene_statistics - Basic statistics on gene expression levels
            4. batch_effect - Analyze potential batch effects in the data

            Generate exactly {num_questions} NEW research questions that can be answered with this data.
            DO NOT repeat questions from the history.

            Respond in JSON format ONLY:
            {{
                "questions": [
                    {{
                        "question": "The research question",
                        "hypothesis": "The expected outcome or hypothesis",
                        "analysis_type": "One of: differential_expression, pca, gene_statistics, batch_effect",
                        "params": {{"key": "value"}}  // Optional parameters for the analysis
                    }}
                ]
            }}
        """).strip()

        print_substep(f"Generating {num_questions} research questions...", "running")
        response = call_hvd_api(prompt, model="gpt-5-mini")

        # Parse response
        questions = []
        try:
            json_match = re.search(r'\{[\s\S]*\}', response)
            if json_match:
                data = json.loads(json_match.group())
                for q_data in data.get("questions", []):
                    self.question_counter += 1
                    rq = ResearchQuestion(
                        question_id=self.question_counter,
                        question=q_data.get("question", ""),
                        hypothesis=q_data.get("hypothesis", ""),
                        analysis_type=q_data.get("analysis_type", "differential_expression"),
                        analysis_params=q_data.get("params", {}),
                        timestamp=time.strftime("%Y-%m-%dT%H:%M:%S"),
                    )
                    questions.append(rq)
                    print_substep(f"Q{rq.question_id}: {rq.question[:50]}...", "done")
        except json.JSONDecodeError as e:
            logging.warning(f"Failed to parse research questions: {e}")
            print_error(f"Failed to parse LLM response: {e}")

        print_timing(start_time, "Question generation")
        return questions

    def execute_analysis(self, rq: ResearchQuestion) -> None:
        """
        Execute the appropriate R analysis for a research question.
        """
        start_time = time.time()

        print_substep(f"Running {rq.analysis_type} analysis...", "running")

        if rq.analysis_type == "differential_expression":
            rq.result = self._run_differential_expression()
        elif rq.analysis_type == "pca":
            rq.result = self._run_pca_analysis()
        elif rq.analysis_type == "gene_statistics":
            rq.result = self._run_gene_statistics()
        elif rq.analysis_type == "batch_effect":
            rq.result = self._run_batch_effect_analysis()
        else:
            rq.result = {"error": f"Unknown analysis type: {rq.analysis_type}"}
            print_error(f"Unknown analysis type: {rq.analysis_type}")
            return

        print_substep("Analysis complete", "done")
        print_timing(start_time, rq.analysis_type)

    def _run_differential_expression(self) -> Dict[str, Any]:
        """Run differential expression analysis using the R pipeline."""
        output_file = self.output_dir / f"de_results_{int(time.time())}.json"

        cmd = [
            "Rscript",
            "src/scripts/run_R_pipeline.R",
            "--experiment-name", self.experiment_name,
            "--assay-file", str(self.assay_file),
            "--sample-file", str(self.sample_file),
            "--reference-condition", self.reference_condition,
            "--group-cols", self.group_cols,
            "--comparison-cols", self.comparison_cols,
            "--skip-mapping",
            "--skip-gsea",
            "--summary-file", str(output_file),
        ]

        try:
            subprocess.run(
                cmd,
                cwd=self.repo_root,
                check=True,
                capture_output=True,
                text=True,
            )

            if output_file.exists():
                with open(output_file, 'r') as f:
                    return json.load(f)
            else:
                return {"error": "Output file not created"}
        except subprocess.CalledProcessError as e:
            return {"error": str(e), "stderr": e.stderr}

    def _run_pca_analysis(self) -> Dict[str, Any]:
        """Run PCA analysis."""
        # Create a simple R script for PCA
        pca_script = textwrap.dedent(f"""
            library(jsonlite)

            # Read data
            assay_data <- read.csv("{self.assay_file}", row.names=1)
            sample_data <- read.csv("{self.sample_file}")

            # Log transform and run PCA
            log_counts <- log2(assay_data + 1)
            pca_result <- prcomp(t(log_counts), scale=TRUE)

            # Extract results
            variance_explained <- summary(pca_result)$importance[2, 1:min(5, ncol(pca_result$x))]
            pc_coords <- as.data.frame(pca_result$x[, 1:min(5, ncol(pca_result$x))])
            pc_coords$sample_id <- rownames(pc_coords)

            # Merge with sample info
            pc_coords <- merge(pc_coords, sample_data, by="sample_id")

            result <- list(
                variance_explained = as.list(variance_explained),
                pc_coordinates = pc_coords,
                num_samples = nrow(pc_coords)
            )

            cat(toJSON(result, auto_unbox=TRUE, pretty=TRUE))
        """)

        try:
            result = subprocess.run(
                ["Rscript", "-e", pca_script],
                cwd=self.repo_root,
                capture_output=True,
                text=True,
                check=True,
            )
            return json.loads(result.stdout)
        except (subprocess.CalledProcessError, json.JSONDecodeError) as e:
            return {"error": str(e)}

    def _run_gene_statistics(self) -> Dict[str, Any]:
        """Calculate basic gene expression statistics."""
        stats_script = textwrap.dedent(f"""
            library(jsonlite)

            assay_data <- read.csv("{self.assay_file}", row.names=1)

            # Calculate statistics
            gene_means <- rowMeans(assay_data)
            gene_vars <- apply(assay_data, 1, var)
            gene_max <- apply(assay_data, 1, max)

            # Top expressed genes
            top_genes <- names(sort(gene_means, decreasing=TRUE)[1:20])

            # Most variable genes
            most_variable <- names(sort(gene_vars, decreasing=TRUE)[1:20])

            result <- list(
                total_genes = nrow(assay_data),
                mean_expression = mean(gene_means),
                median_expression = median(gene_means),
                top_expressed_genes = top_genes,
                most_variable_genes = most_variable,
                zero_count_genes = sum(gene_max == 0)
            )

            cat(toJSON(result, auto_unbox=TRUE, pretty=TRUE))
        """)

        try:
            result = subprocess.run(
                ["Rscript", "-e", stats_script],
                cwd=self.repo_root,
                capture_output=True,
                text=True,
                check=True,
            )
            return json.loads(result.stdout)
        except (subprocess.CalledProcessError, json.JSONDecodeError) as e:
            return {"error": str(e)}

    def _run_batch_effect_analysis(self) -> Dict[str, Any]:
        """Analyze potential batch effects."""
        batch_script = textwrap.dedent(f"""
            library(jsonlite)

            assay_data <- read.csv("{self.assay_file}", row.names=1)
            sample_data <- read.csv("{self.sample_file}")

            # Check if batch column exists
            if (!"batch" %in% colnames(sample_data)) {{
                cat(toJSON(list(error = "No batch column found"), auto_unbox=TRUE))
                quit()
            }}

            # Log transform
            log_counts <- log2(assay_data + 1)

            # Calculate sample means by batch
            sample_means <- colMeans(log_counts)
            batch_means <- tapply(sample_means, sample_data$batch, mean)

            # Simple batch effect test using ANOVA
            batch_anova <- summary(aov(sample_means ~ sample_data$batch))
            p_value <- batch_anova[[1]][["Pr(>F)"]][1]

            result <- list(
                batches = unique(sample_data$batch),
                batch_means = as.list(batch_means),
                batch_effect_pvalue = p_value,
                significant_batch_effect = p_value < 0.05
            )

            cat(toJSON(result, auto_unbox=TRUE, pretty=TRUE))
        """)

        try:
            result = subprocess.run(
                ["Rscript", "-e", batch_script],
                cwd=self.repo_root,
                capture_output=True,
                text=True,
                check=True,
            )
            return json.loads(result.stdout)
        except (subprocess.CalledProcessError, json.JSONDecodeError) as e:
            return {"error": str(e)}

    def validate_result(self, rq: ResearchQuestion) -> None:
        """
        Validate the analysis result using LLM.
        """
        if not rq.result or "error" in rq.result:
            rq.validated = False
            rq.validation_score = 0.0
            rq.interpretation = f"Analysis failed: {rq.result.get('error', 'Unknown error')}"
            return

        print_substep("Validating results via LLM...", "running")

        prompt = textwrap.dedent(f"""
            You are a computational biologist validating analysis results.

            Research Question: {rq.question}
            Hypothesis: {rq.hypothesis}
            Analysis Type: {rq.analysis_type}

            Results:
            {json.dumps(rq.result, indent=2, default=str)[:3000]}

            Evaluate the results and provide:
            1. Whether the analysis successfully addresses the question
            2. An interpretation of the findings
            3. A quality/confidence score (0.0-1.0)

            Respond in JSON format ONLY:
            {{
                "valid": true/false,
                "score": 0.0-1.0,
                "interpretation": "Your interpretation of the results",
                "key_findings": ["Finding 1", "Finding 2"],
                "limitations": ["Limitation 1"]
            }}
        """).strip()

        try:
            response = call_hvd_api(prompt, model="gpt-5-mini")
            json_match = re.search(r'\{[\s\S]*\}', response)
            if json_match:
                validation = json.loads(json_match.group())
                rq.validated = validation.get("valid", False)
                rq.validation_score = validation.get("score", 0.0)
                rq.interpretation = validation.get("interpretation", "")

                # Store additional validation info in result
                rq.result["validation"] = validation

                print_substep(f"Validation score: {rq.validation_score:.2f}", "done")
            else:
                rq.validated = False
                rq.validation_score = 0.0
        except (json.JSONDecodeError, Exception) as e:
            logging.warning(f"Validation failed: {e}")
            rq.validated = False
            rq.validation_score = 0.0

    def save_results(self) -> None:
        """Save all results to the output directory with per-question folders."""
        print_substep("Saving exploration results...", "running")

        # Create questions directory
        questions_dir = self.output_dir / "questions"
        questions_dir.mkdir(parents=True, exist_ok=True)

        # Save data exploration results
        data_file = self.output_dir / "data_exploration.json"
        with open(data_file, 'w') as f:
            json.dump(self.data_summary, f, indent=2, default=str)
        print_substep("Saved data exploration", "done")

        # Save each research question in its own folder
        validated_count = 0
        for rq in self.research_questions:
            rq_dir = questions_dir / f"rq_{rq.question_id:03d}"
            rq_dir.mkdir(parents=True, exist_ok=True)

            # Save question details
            question_data = {
                "question_id": rq.question_id,
                "question": rq.question,
                "hypothesis": rq.hypothesis,
                "analysis_type": rq.analysis_type,
                "analysis_params": rq.analysis_params,
                "timestamp": rq.timestamp,
            }
            with open(rq_dir / "question.json", 'w') as f:
                json.dump(question_data, f, indent=2)

            # Save analysis results
            with open(rq_dir / "result.json", 'w') as f:
                json.dump(rq.result, f, indent=2, default=str)

            # Save validation results
            validation_data = {
                "validated": rq.validated,
                "validation_score": rq.validation_score,
                "interpretation": rq.interpretation,
            }
            with open(rq_dir / "validation.json", 'w') as f:
                json.dump(validation_data, f, indent=2)

            if rq.validated:
                validated_count += 1

            print_substep(f"Saved rq_{rq.question_id:03d}/", "done")

        # Save overall summary
        summary_data = {
            "experiment_name": self.experiment_name,
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "total_questions": len(self.research_questions),
            "validated_questions": validated_count,
            "questions": [
                {
                    "id": rq.question_id,
                    "question": rq.question,
                    "analysis_type": rq.analysis_type,
                    "validated": rq.validated,
                    "score": rq.validation_score,
                }
                for rq in self.research_questions
            ],
        }
        with open(self.output_dir / "summary.json", 'w') as f:
            json.dump(summary_data, f, indent=2)

        # Save markdown summary
        self._save_markdown_summary()

        print_info("Output directory", str(self.output_dir))

    def _save_markdown_summary(self) -> None:
        """Save a human-readable markdown summary."""
        summary_file = self.output_dir / "summary.md"
        validated_count = sum(1 for rq in self.research_questions if rq.validated)

        lines = [
            f"# Self-Exploration Results: {self.experiment_name}",
            f"\n**Generated:** {time.strftime('%Y-%m-%d %H:%M:%S')}",
            f"\n## Data Summary",
            f"\n- **Genes:** {self.data_summary.get('assay', {}).get('num_genes', 'N/A')}",
            f"- **Samples:** {self.data_summary.get('samples', {}).get('num_samples', 'N/A')}",
            f"- **Conditions:** {', '.join(self.data_summary.get('samples', {}).get('conditions', []))}",
            f"\n### Data Insights\n",
            self.data_summary.get('llm_insights', 'No insights available'),
            f"\n## Research Questions & Findings",
            f"\n**Total questions:** {len(self.research_questions)}",
            f"**Validated:** {validated_count}",
            "\n| # | Question | Type | Score | Status |",
            "|---|----------|------|-------|--------|",
        ]

        for rq in self.research_questions:
            status = "✓" if rq.validated else "○"
            q_short = rq.question[:50] + "..." if len(rq.question) > 50 else rq.question
            lines.append(f"| {rq.question_id} | {q_short} | {rq.analysis_type} | {rq.validation_score:.2f} | {status} |")

        lines.append("\n## Detailed Findings\n")

        for rq in self.research_questions:
            lines.extend([
                f"### RQ {rq.question_id}: {rq.question}",
                f"\n**Hypothesis:** {rq.hypothesis}",
                f"\n**Analysis:** {rq.analysis_type}",
                f"\n**Score:** {rq.validation_score:.2f}",
                f"\n**Interpretation:**\n{rq.interpretation or 'N/A'}",
                f"\n*Results saved to:* `questions/rq_{rq.question_id:03d}/`",
                "\n---\n",
            ])

        with open(summary_file, 'w') as f:
            f.write('\n'.join(lines))

        print_info("Summary file", str(summary_file))

    def _print_summary(self) -> None:
        """Print final exploration summary."""
        total_time = time.time() - self.start_time if self.start_time else 0

        print_completion_banner()
        print_section("Exploration Summary", "📋")
        print_info("Total questions generated", str(self.question_counter))
        print_info("Questions validated", str(len(self.research_questions)))
        print_info("Total time", f"{total_time:.1f}s")
        print_info("Results saved to", str(self.output_dir))
