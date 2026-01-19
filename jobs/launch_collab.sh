#!/bin/bash
set -e

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

source .venv/bin/activate

# Set API key
export HVD_API_KEY="${HVD_API_KEY:-jzmwbIL1QJG1vT2AC8yTZzJJFA80xeBGVtrTxfkTHmGu5GIE}"

# Default parameters
MODE="${MODE:-simulated}"
SCENARIO="${SCENARIO:-scenario_2}"
HUMAN_PROFILE="${HUMAN_PROFILE:-domain_expert}"
HUMAN_BELIEFS="${HUMAN_BELIEFS:-exploratory}"
MAX_TURNS="${MAX_TURNS:-15}"
AGENT_LLM="${AGENT_LLM:-gpt-5-mini}"
AGENT_API="${AGENT_API:-hvd}"
USER_LLM="${USER_LLM:-gpt-5-mini}"
USER_API="${USER_API:-hvd}"

echo "=============================================="
echo "Collaborative Agent Experiment"
echo "=============================================="
echo "Mode: $MODE"
echo "Scenario: $SCENARIO"
echo "Human Profile: $HUMAN_PROFILE"
echo "Human Beliefs: $HUMAN_BELIEFS"
echo "Max Turns: $MAX_TURNS"
echo "Agent: $AGENT_LLM ($AGENT_API)"
echo "User: $USER_LLM ($USER_API)"
echo "=============================================="

python -m src.collaborative_agent.run_experiment \
  --mode "$MODE" \
  --scenario "$SCENARIO" \
  --human-profile "$HUMAN_PROFILE" \
  --human-beliefs "$HUMAN_BELIEFS" \
  --max-turns "$MAX_TURNS" \
  --agent-llm "$AGENT_LLM" \
  --agent-api "$AGENT_API" \
  --user-llm "$USER_LLM" \
  --user-api "$USER_API" \
  --output-dir "results/collab_agent" \
  "$@"
