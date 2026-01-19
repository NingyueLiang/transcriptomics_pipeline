# Collaborative Agent

A research platform for studying human-AI collaboration in scientific research tasks.

## Overview

This system enables studying how AI agents can effectively collaborate with human researchers on any domain. The agent operates as a teammate, not just a tool—asking clarifying questions, proposing plans, challenging assumptions when appropriate, and learning from failures.

### Research Questions

1. **Collaboration Formation**: How does the agent form the notion of "collaboration"?
   - Asking clarification questions when uncertain
   - Proposing plans for complex tasks
   - Challenging assumptions when appropriate
   - Reflecting on failures to learn

2. **User Modeling**: How does the agent mentally model the user?
   - Stores observations in persistent memory
   - LLM reasons about user preferences from memory content
   - No numeric belief scores—qualitative understanding

3. **Memory Accumulation**: How does the agent accumulate memory?
   - File-based markdown memory (human-readable)
   - Categories: users, episodes, skills, outcomes
   - Persists across all experiments (shared learning)
   - On-demand retrieval via `search_memory` tool

## Quick Start

```bash
# Set up environment
source .venv/bin/activate
export HVD_API_KEY=<your-key>

# Run a simulated experiment
python -m src.collaborative_agent.run_experiment \
    --mode simulated \
    --scenario scenario_2 \
    --human-profile domain_expert \
    --human-beliefs exploratory \
    --max-turns 10

# Run with different LLMs for agent and simulated user
python -m src.collaborative_agent.run_experiment \
    --mode simulated \
    --agent-llm gpt-5-mini \
    --agent-api hvd \
    --user-llm gpt-5-mini \
    --user-api hvd \
    --scenario scenario_2

# Run an interactive session (real human)
python -m src.collaborative_agent.run_experiment \
    --mode interactive \
    --scenario general \
    --task-description "Help me analyze this dataset"
```

## Architecture

```
collaborative_agent/
├── agent.py           # CollaborativeAgent - the AI agent with collaboration actions
├── simulated_human.py # SimulatedHuman - LLM-powered human simulation
├── environment.py     # CollaborationEnvironment - orchestrates interaction
├── tools.py           # Extensible tool system
├── prompts.py         # System prompts + user profiles + hidden beliefs
├── memory_system.py   # File-based persistent memory (markdown files)
├── scenarios.py       # Predefined experimental scenarios
└── run_experiment.py  # Entry point for experiments

memory/                # Persistent memory storage (created at runtime)
├── users/             # User profiles and preferences
├── episodes/          # Session histories and experiences
├── skills/            # Knowledge about tools and techniques
└── outcomes/          # Successes and failures with lessons
    ├── successes/
    └── failures/
```

## Command Line Arguments

### LLM Configuration

| Argument | Description |
|----------|-------------|
| `--agent-llm` | Model name for the agent (default: `gpt-5-mini`) |
| `--agent-api` | API for agent: `hvd`, `openai`, `anthropic` |
| `--user-llm` | Model name for simulated user (default: `gpt-5-mini`) |
| `--user-api` | API for simulated user: `hvd`, `openai`, `anthropic` |

### Experiment Settings

| Argument | Description |
|----------|-------------|
| `--mode` | `simulated` or `interactive` |
| `--scenario` | Scenario name or `general` |
| `--task-description` | Custom task description |
| `--human-profile` | User profile for simulation |
| `--human-beliefs` | Hidden beliefs for simulation |
| `--max-turns` | Maximum conversation turns |
| `--output-dir` | Output directory (default: `results/collab_agent`) |

## Components

### CollaborativeAgent

The main agent that:
- Responds to user input with structured actions
- Can ask clarifications, propose plans, execute tools
- Tracks its internal state for research analysis

**Collaboration Actions:**
- `ASK_CLARIFICATION` - Ask when information is missing or ambiguous
- `PROPOSE_PLAN` - Propose a plan when task has multiple steps
- `CHALLENGE` - Challenge assumptions when they may be incorrect
- `REFLECT` - Analyze failures to learn from them
- `SUMMARIZE` - Checkpoint mutual understanding
- `DEFER` - Admit when uncertain how to proceed
- `EXECUTE_TOOL` - Take action with tools

### SimulatedHuman

LLM-powered simulation of a human researcher:
- Has configurable profile (expertise level, communication style)
- Has hidden beliefs/goals that influence behavior
- Tracks internal state (satisfaction, confusion, goal progress)

**Predefined Profiles:**
- `domain_expert` - Senior researcher, technical, direct
- `domain_novice` - Different field, needs explanations
- `student` - Learning, curious, needs guidance
- `busy_pi` - Wants quick results, impatient
- `skeptical_reviewer` - Questions everything

**Hidden Beliefs:**
- `has_hypothesis` - Has specific hypothesis to test
- `skeptical_of_ai` - Doesn't fully trust AI
- `time_pressured` - Under deadline pressure
- `exploratory` - Wants open-ended exploration
- `perfectionist` - High standards for quality
- `collaborative` - Enjoys working together

### Tools

**General Tools (always available):**
- `read_file` - Read file contents
- `list_directory` - List directory contents
- `run_python` - Execute Python code (requires confirmation)
- `search_files` - Search for patterns in files

**Transcriptomics Tools (domain-specific):**
- `describe_data` - Describe dataset structure
- `run_differential_expression` - DESeq2 analysis
- `run_pca` - Principal component analysis
- `get_gene_statistics` - Expression statistics
- `analyze_gene` - Analyze specific gene
- `check_batch_effects` - Batch effect analysis

### Memory System

File-based persistent memory stored as human-readable markdown files.

**Memory Categories:**
- **users/** - Information about users (preferences, expertise, communication style)
- **episodes/** - What happened in sessions (decisions, outcomes, context)
- **skills/** - Knowledge about tools, techniques, what works
- **outcomes/** - Successes and failures with lessons learned

**Memory Tools (available to agent):**
- `add_to_memory(category, content, topic, subcategory, tags)` - Store information
- `search_memory(query, category)` - Grep-like search across memory files

## Experiment Output

Results save to `results/collab_agent/` with:

```json
{
    "experiment_name": "...",
    "llm_config": {
        "agent_llm": "gpt-5-mini",
        "agent_api": "hvd",
        "user_llm": "gpt-5-mini",
        "user_api": "hvd"
    },
    "config": {...},
    "metrics": {
        "total_turns": 10,
        "agent_questions_asked": 3,
        "agent_plans_proposed": 1,
        "final_human_satisfaction": 0.85,
        "final_goal_progress": 0.9
    },
    "conversation_log": [...],
    "human_internal_thoughts": [...],
    "agent_state": {...}
}
```

## Extending the System

### Adding New Tools

```python
from src.collaborative_agent.tools import ToolRegistry, Tool

registry = ToolRegistry(repo_root, domain="general")

registry.register(Tool(
    name="my_tool",
    description="What the tool does",
    parameters={"param1": {"type": "string", "description": "..."}},
    function=my_function,
    requires_confirmation=False,
    category="custom",
))
```

### Adding New User Profiles

```python
# In prompts.py
USER_PROFILES["my_profile"] = {
    "profile": "Description of who this user is",
    "knowledge_level": "What they know",
    "communication_style": "How they communicate",
}
```

### Using Different APIs

```python
from src.collaborative_agent.run_experiment import create_llm_call

# Harvard API
agent_llm = create_llm_call("hvd", "gpt-5-mini")

# OpenAI
agent_llm = create_llm_call("openai", "gpt-4")

# Anthropic
agent_llm = create_llm_call("anthropic", "claude-sonnet-4-20250514")
```

## Environment Variables

| Variable | Required For |
|----------|--------------|
| `HVD_API_KEY` | Harvard OpenAI Direct API |
| `OPENAI_API_KEY` | OpenAI API |
| `ANTHROPIC_API_KEY` | Anthropic API |
