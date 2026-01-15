#!/bin/bash
set -e

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

source .venv/bin/activate

export HVD_API_KEY=jzmwbIL1QJG1vT2AC8yTZzJJFA80xeBGVtrTxfkTHmGu5GIE

python main.py \
  --experiment-name scenario_11a \
  --assay-file scenarios/scenario_11a/scenario_11a_assay_data.csv \
  --sample-file scenarios/scenario_11a/scenario_11a_sample_data.csv \
  --reference-condition Control \
  --group-cols "condition,replicate,batch" \
  --comparison-cols "condition" \
  --skip-mapping \
  --skip-gsea \
  "$@"
