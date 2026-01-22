# Transcriptomics Multi-Agent System

A Multi-Agent System (MAS) for Xenopus laevis transcriptomics analysis, combining R-based bioinformatics tools with LLM-powered analysis and report generation.

## Overview

This system provides automated workflows for differential expression analysis, gene mapping to human orthologs, and pathway enrichment. It was developed to analyze transcriptomics data from Xenopus laevis experiments, particularly focusing on anesthetic drug responses.

## Key Features

- **Differential Expression**: DESeq2-based analysis with automated contrast generation
- **Gene Mapping**: Xenopus to human ortholog mapping using Xenbase and HCOP databases
- **Pathway Analysis**: GSEA using MSigDB gene sets with human gene symbols
- **LLM-Powered Reports**: Automated report generation with verification
- **Multi-Agent Architecture**: Three specialized agents for different analysis workflows

## Agent Types

| Agent | Purpose | Input | Output |
|-------|---------|-------|--------|
| **QA Agent** | Question-Answering | Question + Database | Report |
| **SE Agent** | Self-Exploration | Database | Research questions + Results |
| **Collaborative Agent** | Human-AI Research | Interactive conversation | Analysis results + Memory |

## Installation

```bash
# Clone the repository
git clone https://github.com/agrossberg/transcriptomics_pipeline.git
cd transcriptomics_pipeline

# Set up Python environment
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Set API key
export HVD_API_KEY=<your-key>
```

## Quick Start

### Question-Answering Agent

Runs the full R pipeline and generates an LLM report based on your query:

```bash
python main.py --agent-type qa \
  --experiment-name my_experiment \
  --assay-file data/assay_data.csv \
  --sample-file data/sample_data.csv \
  --user-query "What genes are differentially expressed between conditions?"
```

### Self-Exploration Agent

Autonomously generates and investigates research questions:

```bash
python main.py --agent-type se \
  --experiment-name my_experiment \
  --assay-file data/assay_data.csv \
  --sample-file data/sample_data.csv \
  --num-rqs 10
```

### Collaborative Agent

Interactive human-AI research collaboration with persistent memory:

```bash
bash jobs/launch_collab.sh

# Or with custom parameters:
MODE=simulated SCENARIO=scenario_2 MAX_TURNS=15 bash jobs/launch_collab.sh
```

### Using Launch Scripts

Launch scripts set up the environment and API keys automatically:

```bash
bash jobs/launch_qa.sh   # QA agent
bash jobs/launch_se.sh   # SE agent
```

## Data Requirements

### Input Files

**Assay data (CSV):**
- First column: gene symbols (used as row names)
- Remaining columns: sample expression counts

**Sample metadata (CSV):**
- `sample_id`: Unique identifiers matching assay column names
- `condition`: Treatment conditions (e.g., Vehicle, Ketamine_100)
- `timepoint`: Time points (e.g., T, R)
- Optional: `batch`, `replicate`, `dose`, `drug`

**Orthology files** (for gene mapping):
- Xenbase orthology predictions
- HCOP orthology predictions

## Output Structure

```
results/
└── experiment_name/
    ├── qa_agent/              # QA agent outputs
    │   ├── pipeline_summary.json
    │   └── report.md
    ├── se_agent/              # SE agent outputs
    │   ├── summary.json
    │   ├── summary.md
    │   └── questions/
    │       └── rq_XXX/
    ├── collab_agent/          # Collaborative agent outputs
    ├── de/                    # Differential expression
    │   ├── csv/
    │   ├── tables/
    │   └── plots/
    ├── mapping/               # Gene mapping results
    ├── pca/                   # PCA analysis
    └── gsea/                  # Pathway enrichment
```

## Command Line Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--agent-type` | Agent type: `qa` or `se` | `qa` |
| `--experiment-name` | Experiment identifier | `experiment` |
| `--assay-file` | Path to assay counts CSV | Required |
| `--sample-file` | Path to sample metadata CSV | Required |
| `--user-query` | Question for QA agent | General analysis |
| `--num-rqs` | Number of research questions (SE) | `5` |
| `--reference-condition` | Reference for DESeq2 | `Control` |
| `--skip-r-setup` | Skip R environment bootstrap | `false` |
| `--skip-llm` | Skip LLM report generation | `false` |
| `--skip-mapping` | Skip gene mapping step | `false` |
| `--skip-gsea` | Skip GSEA step | `false` |

## Architecture

```
main.py                          # Entry point
src/
├── question_answering_agent/    # QA agent
├── self_exploration_agent/      # SE agent
├── collaborative_agent/         # Collaborative agent
│   ├── agent.py                 # Core agent logic
│   ├── tools.py                 # Analysis tools
│   ├── memory_system.py         # Persistent memory
│   └── run_experiment.py        # Experiment runner
├── utils/                       # Shared utilities
├── scripts/                     # R pipeline scripts
└── R_tools/                     # R analysis functions
```

## R Dependencies

Core packages: `DESeq2`, `fgsea`, `msigdbr`, `jsonlite`, `optparse`, `here`, `dplyr`, `tidyr`, `ggplot2`, `ComplexHeatmap`, `viridis`

The R environment is set up automatically on first run, or manually via:
```bash
Rscript src/scripts/setup_R_env.R
```

## Contact

**Allison Grossberg**
Wyss Institute for Biologically Inspired Engineering
Harvard University
Email: allison.grossberg@wyss.harvard.edu

## License

MIT License - see LICENSE file for details.
