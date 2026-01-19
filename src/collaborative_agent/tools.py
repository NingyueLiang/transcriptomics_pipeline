"""
Tools available to the collaborative agent.

This module provides a flexible tool system that can be extended for different domains.
Default tools are provided for transcriptomics analysis, but custom tools can be registered.
"""

from __future__ import annotations

import json
import subprocess
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional


@dataclass
class Tool:
    """Represents a tool that the agent can use."""
    name: str
    description: str
    parameters: Dict[str, Any]
    function: Callable
    requires_confirmation: bool = False  # If True, agent should confirm with user before running
    category: str = "general"  # Tool category for organization


class ToolRegistry:
    """
    Registry of tools available to the agent.

    Tools can be registered dynamically for different domains.
    """

    def __init__(self, repo_root: Path, domain: str = "general"):
        """
        Initialize the tool registry.

        Args:
            repo_root: Path to the repository root.
            domain: Domain for default tools. Options: "general", "transcriptomics"
        """
        self.repo_root = repo_root
        self.domain = domain
        self.tools: Dict[str, Tool] = {}

        # Register default tools based on domain
        self._register_general_tools()
        if domain == "transcriptomics":
            self._register_transcriptomics_tools()

    def _register_general_tools(self) -> None:
        """Register general-purpose tools."""

        # Tool: Read file
        self.register(Tool(
            name="read_file",
            description="Read the contents of a file. Useful for examining data, code, or documentation.",
            parameters={
                "file_path": {"type": "string", "description": "Path to the file to read"},
                "max_lines": {"type": "integer", "description": "Maximum lines to read (default: 100)", "default": 100},
            },
            function=self._read_file,
            requires_confirmation=False,
            category="general",
        ))

        # Tool: List directory
        self.register(Tool(
            name="list_directory",
            description="List files and directories in a path. Useful for exploring data structure.",
            parameters={
                "path": {"type": "string", "description": "Directory path to list"},
                "pattern": {"type": "string", "description": "Optional glob pattern to filter", "default": "*"},
            },
            function=self._list_directory,
            requires_confirmation=False,
            category="general",
        ))

        # Tool: Run Python code
        self.register(Tool(
            name="run_python",
            description="Execute Python code and return the result. Useful for data analysis and computations.",
            parameters={
                "code": {"type": "string", "description": "Python code to execute"},
            },
            function=self._run_python,
            requires_confirmation=True,
            category="general",
        ))

        # Tool: Search in files
        self.register(Tool(
            name="search_files",
            description="Search for a pattern in files. Useful for finding specific content.",
            parameters={
                "pattern": {"type": "string", "description": "Text or regex pattern to search"},
                "path": {"type": "string", "description": "Directory to search in", "default": "."},
                "file_pattern": {"type": "string", "description": "File glob pattern", "default": "*"},
            },
            function=self._search_files,
            requires_confirmation=False,
            category="general",
        ))

    def _register_transcriptomics_tools(self) -> None:
        """Register transcriptomics-specific tools."""

        # Tool: Describe data
        self.register(Tool(
            name="describe_data",
            description="Load and describe the structure of a transcriptomics dataset. Returns information about genes, samples, conditions, and basic statistics.",
            parameters={
                "assay_file": {"type": "string", "description": "Path to the assay data CSV"},
                "sample_file": {"type": "string", "description": "Path to the sample metadata CSV"},
            },
            function=self._describe_data,
            requires_confirmation=False,
            category="transcriptomics",
        ))

        # Tool: Run differential expression analysis
        self.register(Tool(
            name="run_differential_expression",
            description="Run DESeq2 differential expression analysis comparing conditions. Returns significantly differentially expressed genes.",
            parameters={
                "assay_file": {"type": "string", "description": "Path to the assay data CSV"},
                "sample_file": {"type": "string", "description": "Path to the sample metadata CSV"},
                "reference_condition": {"type": "string", "description": "Reference condition for comparison"},
                "comparison_cols": {"type": "string", "description": "Columns to use for comparison (comma-separated)"},
            },
            function=self._run_differential_expression,
            requires_confirmation=True,
            category="transcriptomics",
        ))

        # Tool: Run PCA analysis
        self.register(Tool(
            name="run_pca",
            description="Run Principal Component Analysis to visualize sample clustering and identify batch effects.",
            parameters={
                "assay_file": {"type": "string", "description": "Path to the assay data CSV"},
                "sample_file": {"type": "string", "description": "Path to the sample metadata CSV"},
            },
            function=self._run_pca,
            requires_confirmation=False,
            category="transcriptomics",
        ))

        # Tool: Get gene statistics
        self.register(Tool(
            name="get_gene_statistics",
            description="Calculate expression statistics for genes, including top expressed and most variable genes.",
            parameters={
                "assay_file": {"type": "string", "description": "Path to the assay data CSV"},
                "top_n": {"type": "integer", "description": "Number of top genes to return", "default": 20},
            },
            function=self._get_gene_statistics,
            requires_confirmation=False,
            category="transcriptomics",
        ))

        # Tool: Analyze specific gene
        self.register(Tool(
            name="analyze_gene",
            description="Get detailed information about a specific gene's expression across conditions.",
            parameters={
                "assay_file": {"type": "string", "description": "Path to the assay data CSV"},
                "sample_file": {"type": "string", "description": "Path to the sample metadata CSV"},
                "gene_symbol": {"type": "string", "description": "Gene symbol to analyze"},
            },
            function=self._analyze_gene,
            requires_confirmation=False,
            category="transcriptomics",
        ))

        # Tool: Check batch effects
        self.register(Tool(
            name="check_batch_effects",
            description="Analyze potential batch effects in the data.",
            parameters={
                "assay_file": {"type": "string", "description": "Path to the assay data CSV"},
                "sample_file": {"type": "string", "description": "Path to the sample metadata CSV"},
            },
            function=self._check_batch_effects,
            requires_confirmation=False,
            category="transcriptomics",
        ))

    def register(self, tool: Tool) -> None:
        """Register a new tool."""
        self.tools[tool.name] = tool

    def unregister(self, name: str) -> bool:
        """Unregister a tool by name. Returns True if tool was removed."""
        if name in self.tools:
            del self.tools[name]
            return True
        return False

    def get_tool(self, name: str) -> Optional[Tool]:
        """Get a tool by name."""
        return self.tools.get(name)

    def list_tools(self, category: Optional[str] = None) -> List[str]:
        """List all tool names, optionally filtered by category."""
        if category:
            return [name for name, tool in self.tools.items() if tool.category == category]
        return list(self.tools.keys())

    def execute(self, tool_name: str, params: Dict[str, Any]) -> Dict[str, Any]:
        """Execute a tool with given parameters."""
        tool = self.get_tool(tool_name)
        if tool is None:
            return {"error": f"Unknown tool: {tool_name}", "available_tools": self.list_tools()}

        try:
            return tool.function(**params)
        except Exception as e:
            return {"error": str(e), "tool": tool_name, "params": params}

    def get_tools_description(self, category: Optional[str] = None) -> str:
        """Get a formatted description of tools, optionally filtered by category."""
        lines = []
        for name, tool in self.tools.items():
            if category and tool.category != category:
                continue

            params_str = ", ".join(
                f"{k}: {v.get('type', 'any')}"
                for k, v in tool.parameters.items()
            )
            confirm_note = " [REQUIRES USER CONFIRMATION]" if tool.requires_confirmation else ""
            lines.append(f"- {name}({params_str}): {tool.description}{confirm_note}")
        return "\n".join(lines)

    # ==========================================================================
    # General Tool Implementations
    # ==========================================================================

    def _read_file(self, file_path: str, max_lines: int = 100) -> Dict[str, Any]:
        """Read a file's contents."""
        try:
            path = Path(file_path)
            if not path.is_absolute():
                path = self.repo_root / path

            if not path.exists():
                return {"error": f"File not found: {file_path}"}

            with open(path, 'r', encoding='utf-8', errors='replace') as f:
                lines = f.readlines()

            total_lines = len(lines)
            content = ''.join(lines[:max_lines])

            return {
                "success": True,
                "file_path": str(path),
                "total_lines": total_lines,
                "lines_returned": min(max_lines, total_lines),
                "content": content,
                "truncated": total_lines > max_lines,
            }
        except Exception as e:
            return {"error": str(e)}

    def _list_directory(self, path: str, pattern: str = "*") -> Dict[str, Any]:
        """List directory contents."""
        try:
            dir_path = Path(path)
            if not dir_path.is_absolute():
                dir_path = self.repo_root / dir_path

            if not dir_path.exists():
                return {"error": f"Directory not found: {path}"}

            if not dir_path.is_dir():
                return {"error": f"Not a directory: {path}"}

            items = list(dir_path.glob(pattern))
            files = [str(p.relative_to(dir_path)) for p in items if p.is_file()]
            dirs = [str(p.relative_to(dir_path)) + "/" for p in items if p.is_dir()]

            return {
                "success": True,
                "path": str(dir_path),
                "directories": sorted(dirs),
                "files": sorted(files),
                "total_items": len(items),
            }
        except Exception as e:
            return {"error": str(e)}

    def _run_python(self, code: str) -> Dict[str, Any]:
        """Execute Python code."""
        try:
            # Create a restricted execution environment
            local_vars = {}
            exec(code, {"__builtins__": __builtins__}, local_vars)

            # Try to capture a result
            result = local_vars.get('result', local_vars.get('output', None))

            return {
                "success": True,
                "result": result,
                "local_vars": {k: str(v)[:500] for k, v in local_vars.items() if not k.startswith('_')},
            }
        except Exception as e:
            return {"error": str(e), "code": code[:200]}

    def _search_files(self, pattern: str, path: str = ".", file_pattern: str = "*") -> Dict[str, Any]:
        """Search for pattern in files."""
        try:
            search_path = Path(path)
            if not search_path.is_absolute():
                search_path = self.repo_root / search_path

            matches = []
            for file_path in search_path.rglob(file_pattern):
                if file_path.is_file():
                    try:
                        with open(file_path, 'r', encoding='utf-8', errors='replace') as f:
                            for i, line in enumerate(f, 1):
                                if pattern.lower() in line.lower():
                                    matches.append({
                                        "file": str(file_path.relative_to(search_path)),
                                        "line_number": i,
                                        "line": line.strip()[:200],
                                    })
                                    if len(matches) >= 50:  # Limit results
                                        break
                    except:
                        pass
                if len(matches) >= 50:
                    break

            return {
                "success": True,
                "pattern": pattern,
                "search_path": str(search_path),
                "matches": matches,
                "total_matches": len(matches),
                "truncated": len(matches) >= 50,
            }
        except Exception as e:
            return {"error": str(e)}

    # ==========================================================================
    # Transcriptomics Tool Implementations
    # ==========================================================================

    def _describe_data(self, assay_file: str, sample_file: str) -> Dict[str, Any]:
        """Describe the dataset structure."""
        import csv

        try:
            # Read assay data
            with open(assay_file, 'r') as f:
                reader = csv.reader(f)
                header = next(reader)
                gene_count = sum(1 for _ in reader)

            sample_cols = header[1:]

            # Read sample data
            with open(sample_file, 'r') as f:
                reader = csv.DictReader(f)
                rows = list(reader)

            if not rows:
                return {"error": "No sample data found"}

            columns = list(rows[0].keys())
            conditions = list(set(r.get("condition", "") for r in rows if r.get("condition")))
            batches = list(set(r.get("batch", "") for r in rows if r.get("batch")))

            return {
                "success": True,
                "assay_info": {
                    "num_genes": gene_count,
                    "num_samples": len(sample_cols),
                    "sample_ids": sample_cols,
                },
                "sample_info": {
                    "num_samples": len(rows),
                    "columns": columns,
                    "conditions": conditions,
                    "num_conditions": len(conditions),
                    "batches": batches,
                    "num_batches": len(batches),
                },
            }
        except Exception as e:
            return {"error": str(e)}

    def _run_differential_expression(
        self,
        assay_file: str,
        sample_file: str,
        reference_condition: str,
        comparison_cols: str = "condition",
    ) -> Dict[str, Any]:
        """Run differential expression analysis."""
        import tempfile

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            output_file = f.name

        cmd = [
            "Rscript",
            "src/scripts/run_R_pipeline.R",
            "--experiment-name", "collab_de_temp",
            "--assay-file", assay_file,
            "--sample-file", sample_file,
            "--reference-condition", reference_condition,
            "--comparison-cols", comparison_cols,
            "--skip-mapping",
            "--skip-gsea",
            "--summary-file", output_file,
        ]

        try:
            subprocess.run(
                cmd,
                cwd=self.repo_root,
                check=True,
                capture_output=True,
                text=True,
            )

            with open(output_file, 'r') as f:
                result = json.load(f)

            Path(output_file).unlink(missing_ok=True)
            return {"success": True, "results": result}
        except subprocess.CalledProcessError as e:
            return {"error": str(e), "stderr": e.stderr}

    def _run_pca(self, assay_file: str, sample_file: str) -> Dict[str, Any]:
        """Run PCA analysis."""
        pca_script = textwrap.dedent(f"""
            library(jsonlite)
            assay_data <- read.csv("{assay_file}", row.names=1)
            sample_data <- read.csv("{sample_file}")
            log_counts <- log2(assay_data + 1)
            pca_result <- prcomp(t(log_counts), scale=TRUE)
            variance_explained <- summary(pca_result)$importance[2, 1:min(5, ncol(pca_result$x))]
            pc_coords <- as.data.frame(pca_result$x[, 1:min(5, ncol(pca_result$x))])
            pc_coords$sample_id <- rownames(pc_coords)
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
            return {"success": True, "results": json.loads(result.stdout)}
        except (subprocess.CalledProcessError, json.JSONDecodeError) as e:
            return {"error": str(e)}

    def _get_gene_statistics(self, assay_file: str, top_n: int = 20) -> Dict[str, Any]:
        """Get gene expression statistics."""
        stats_script = textwrap.dedent(f"""
            library(jsonlite)
            assay_data <- read.csv("{assay_file}", row.names=1)
            gene_means <- rowMeans(assay_data)
            gene_vars <- apply(assay_data, 1, var)
            gene_max <- apply(assay_data, 1, max)
            top_genes <- names(sort(gene_means, decreasing=TRUE)[1:{top_n}])
            most_variable <- names(sort(gene_vars, decreasing=TRUE)[1:{top_n}])
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
            return {"success": True, "results": json.loads(result.stdout)}
        except (subprocess.CalledProcessError, json.JSONDecodeError) as e:
            return {"error": str(e)}

    def _analyze_gene(self, assay_file: str, sample_file: str, gene_symbol: str) -> Dict[str, Any]:
        """Analyze a specific gene."""
        gene_script = textwrap.dedent(f"""
            library(jsonlite)
            assay_data <- read.csv("{assay_file}", row.names=1)
            sample_data <- read.csv("{sample_file}")
            gene_symbol <- "{gene_symbol}"
            if (!gene_symbol %in% rownames(assay_data)) {{
                cat(toJSON(list(error = paste("Gene not found:", gene_symbol)), auto_unbox=TRUE))
                quit()
            }}
            gene_expr <- as.numeric(assay_data[gene_symbol, ])
            names(gene_expr) <- colnames(assay_data)
            expr_df <- data.frame(sample_id = names(gene_expr), expression = gene_expr)
            expr_df <- merge(expr_df, sample_data, by="sample_id")
            if ("condition" %in% colnames(expr_df)) {{
                condition_stats <- aggregate(expression ~ condition, data=expr_df,
                    FUN=function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
            }} else {{
                condition_stats <- NULL
            }}
            result <- list(
                gene_symbol = gene_symbol,
                mean_expression = mean(gene_expr),
                sd_expression = sd(gene_expr),
                min_expression = min(gene_expr),
                max_expression = max(gene_expr),
                expression_by_sample = as.list(gene_expr),
                condition_stats = condition_stats
            )
            cat(toJSON(result, auto_unbox=TRUE, pretty=TRUE))
        """)

        try:
            result = subprocess.run(
                ["Rscript", "-e", gene_script],
                cwd=self.repo_root,
                capture_output=True,
                text=True,
                check=True,
            )
            return {"success": True, "results": json.loads(result.stdout)}
        except (subprocess.CalledProcessError, json.JSONDecodeError) as e:
            return {"error": str(e)}

    def _check_batch_effects(self, assay_file: str, sample_file: str) -> Dict[str, Any]:
        """Check for batch effects."""
        batch_script = textwrap.dedent(f"""
            library(jsonlite)
            assay_data <- read.csv("{assay_file}", row.names=1)
            sample_data <- read.csv("{sample_file}")
            if (!"batch" %in% colnames(sample_data)) {{
                cat(toJSON(list(error = "No batch column found in sample data"), auto_unbox=TRUE))
                quit()
            }}
            log_counts <- log2(assay_data + 1)
            sample_means <- colMeans(log_counts)
            batch_means <- tapply(sample_means, sample_data$batch, mean)
            batch_anova <- summary(aov(sample_means ~ sample_data$batch))
            p_value <- batch_anova[[1]][["Pr(>F)"]][1]
            result <- list(
                batches = unique(sample_data$batch),
                num_batches = length(unique(sample_data$batch)),
                batch_means = as.list(batch_means),
                batch_effect_pvalue = p_value,
                significant_batch_effect = p_value < 0.05,
                recommendation = ifelse(p_value < 0.05,
                    "Significant batch effect detected. Consider batch correction.",
                    "No significant batch effect detected.")
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
            return {"success": True, "results": json.loads(result.stdout)}
        except (subprocess.CalledProcessError, json.JSONDecodeError) as e:
            return {"error": str(e)}
