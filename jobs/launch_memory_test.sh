#!/bin/bash
# Simple test script for memory system development
# Uses scripted human responses (no LLM needed for user)

set -e

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

source .venv/bin/activate

# Set API key for agent
export HVD_API_KEY="${HVD_API_KEY:-jzmwbIL1QJG1vT2AC8yTZzJJFA80xeBGVtrTxfkTHmGu5GIE}"

echo "=============================================="
echo "Memory System Test"
echo "=============================================="
echo "Mode: scripted (predefined user responses)"
echo "Scenario: memory_test"
echo "=============================================="

python -m src.collaborative_agent.run_experiment \
  --mode scripted \
  --scenario memory_test \
  --max-turns 6 \
  --agent-llm gpt-5-mini \
  --agent-api hvd \
  --output-dir results/collab_agent \
  "$@"
