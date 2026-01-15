#!/bin/bash
set -e

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

source .venv/bin/activate

export HVD_API_KEY=jzmwbIL1QJG1vT2AC8yTZzJJFA80xeBGVtrTxfkTHmGu5GIE

# Run Self-Exploration agent
python main.py \
  --agent-type se \
  --experiment-name scenario_2 \
  --assay-file scenarios/scenario_2/scenario_2_assay_data.csv \
  --sample-file scenarios/scenario_2/scenario_2_sample_data.csv \
  --reference-condition Control \
  --group-cols "condition,replicate,batch" \
  --comparison-cols "condition" \
  --num-rqs 5 \
  "$@"
