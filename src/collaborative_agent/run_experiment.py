#!/usr/bin/env python3
"""
Run collaborative agent experiments.

Examples:
    # Run a simulated experiment with same LLM for agent and user
    python -m src.collaborative_agent.run_experiment \
        --mode simulated \
        --scenario scenario_2 \
        --human-profile domain_expert \
        --human-beliefs has_hypothesis \
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
        --scenario scenario_2
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Callable, Optional

# Add repo root to path
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from dotenv import load_dotenv

from src.collaborative_agent.environment import (
    CollaborationEnvironment,
    InteractiveEnvironment,
    ExperimentConfig,
)
from src.collaborative_agent.prompts import USER_PROFILES, HIDDEN_BELIEFS_TEMPLATES

load_dotenv()


# =============================================================================
# LLM API Wrappers
# =============================================================================

def create_hvd_llm_call(model: str = "gpt-5-mini") -> Callable[[str, str], str]:
    """Create a function that calls the Harvard OpenAI Direct API."""
    from src.utils.hvd_api import call_hvd_api

    def llm_call(prompt: str, system_prompt: str) -> str:
        full_prompt = f"""<system>
{system_prompt}
</system>

<user>
{prompt}
</user>"""
        return call_hvd_api(full_prompt, model=model)

    return llm_call


def create_openai_llm_call(model: str = "gpt-4") -> Callable[[str, str], str]:
    """Create a function that calls the OpenAI API directly."""
    try:
        import openai
    except ImportError:
        raise ImportError("openai package not installed. Run: pip install openai")

    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise ValueError("OPENAI_API_KEY environment variable not set")

    client = openai.OpenAI(api_key=api_key)

    def llm_call(prompt: str, system_prompt: str) -> str:
        response = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": prompt},
            ],
        )
        return response.choices[0].message.content

    return llm_call


def create_anthropic_llm_call(model: str = "claude-sonnet-4-20250514") -> Callable[[str, str], str]:
    """Create a function that calls the Anthropic API."""
    try:
        import anthropic
    except ImportError:
        raise ImportError("anthropic package not installed. Run: pip install anthropic")

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        raise ValueError("ANTHROPIC_API_KEY environment variable not set")

    client = anthropic.Anthropic(api_key=api_key)

    def llm_call(prompt: str, system_prompt: str) -> str:
        response = client.messages.create(
            model=model,
            max_tokens=4096,
            system=system_prompt,
            messages=[{"role": "user", "content": prompt}],
        )
        return response.content[0].text

    return llm_call


def create_llm_call(api: str, model: str) -> Callable[[str, str], str]:
    """
    Create an LLM call function based on API type and model.

    Args:
        api: API type - "hvd", "openai", or "anthropic"
        model: Model name to use

    Returns:
        Function that takes (prompt, system_prompt) and returns response string.
    """
    if api == "hvd":
        return create_hvd_llm_call(model)
    elif api == "openai":
        return create_openai_llm_call(model)
    elif api == "anthropic":
        return create_anthropic_llm_call(model)
    else:
        raise ValueError(f"Unknown API type: {api}. Supported: hvd, openai, anthropic")


# =============================================================================
# Scenario Configuration
# =============================================================================

def get_scenario_config(scenario_name: str) -> dict:
    """Get configuration for a predefined scenario."""
    scenarios = {
        # Simple test scenarios (for memory development)
        "simple_hello": {
            "task_description": "Simple conversation test.",
            "data_files": {},
            "script_name": "simple_hello",
        },
        "memory_test": {
            "task_description": "Test memory system with simple questions.",
            "data_files": {},
            "script_name": "memory_test",
        },
        # Research scenarios
        "scenario_2": {
            "task_description": "Analyze transcriptomics data to identify differentially expressed genes between treatment and control conditions.",
            "data_files": {
                "assay_file": f"scenarios/{scenario_name}/{scenario_name}_assay_data.csv",
                "sample_file": f"scenarios/{scenario_name}/{scenario_name}_sample_data.csv",
            },
        },
        "scenario_11a": {
            "task_description": "Analyze Xenopus transcriptomics data from anesthetic experiments to understand drug response mechanisms.",
            "data_files": {
                "assay_file": f"scenarios/{scenario_name}/{scenario_name}_assay_data.csv",
                "sample_file": f"scenarios/{scenario_name}/{scenario_name}_sample_data.csv",
            },
        },
        "general": {
            "task_description": "General research collaboration task. The specific task will be defined by the user.",
            "data_files": {},
        },
    }

    if scenario_name in scenarios:
        return scenarios[scenario_name]
    else:
        # Try to construct a generic config
        return {
            "task_description": f"Research task: {scenario_name}",
            "data_files": {},
        }


# =============================================================================
# Argument Parsing
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run collaborative agent experiments",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Mode
    parser.add_argument(
        "--mode",
        choices=["simulated", "interactive", "scripted"],
        default="simulated",
        help="Experiment mode: 'simulated' uses LLM-simulated human, 'interactive' uses real human, 'scripted' uses predefined responses",
    )

    # Scenario
    parser.add_argument(
        "--scenario",
        default="general",
        help="Scenario name (use 'general' for open-ended tasks)",
    )

    parser.add_argument(
        "--task-description",
        default=None,
        help="Custom task description (overrides scenario default)",
    )

    # Agent LLM configuration
    agent_group = parser.add_argument_group("Agent LLM Configuration")
    agent_group.add_argument(
        "--agent-llm",
        default="gpt-5-mini",
        help="Model name for the agent LLM",
    )
    agent_group.add_argument(
        "--agent-api",
        choices=["hvd", "openai", "anthropic"],
        default="hvd",
        help="API to use for agent LLM",
    )

    # Simulated User LLM configuration
    user_group = parser.add_argument_group("Simulated User LLM Configuration")
    user_group.add_argument(
        "--user-llm",
        default="gpt-5-mini",
        help="Model name for the simulated user LLM",
    )
    user_group.add_argument(
        "--user-api",
        choices=["hvd", "openai", "anthropic"],
        default="hvd",
        help="API to use for simulated user LLM",
    )

    # Human simulation settings
    human_group = parser.add_argument_group("Simulated Human Settings")
    human_group.add_argument(
        "--human-profile",
        choices=list(USER_PROFILES.keys()),
        default="domain_expert",
        help="Profile for simulated human",
    )
    human_group.add_argument(
        "--human-beliefs",
        choices=list(HIDDEN_BELIEFS_TEMPLATES.keys()),
        default="exploratory",
        help="Hidden beliefs for simulated human",
    )

    # Experiment settings
    exp_group = parser.add_argument_group("Experiment Settings")
    exp_group.add_argument(
        "--max-turns",
        type=int,
        default=15,
        help="Maximum number of conversation turns",
    )
    exp_group.add_argument(
        "--output-dir",
        default="results/collab_agent",
        help="Directory to save experiment results",
    )
    exp_group.add_argument(
        "--experiment-name",
        default=None,
        help="Name for this experiment (auto-generated if not provided)",
    )

    return parser.parse_args()


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()

    # Validate API keys based on selected APIs
    if args.agent_api == "hvd" or (args.mode == "simulated" and args.user_api == "hvd"):
        if not os.environ.get("HVD_API_KEY"):
            print("Error: HVD_API_KEY environment variable not set")
            print("Please set it: export HVD_API_KEY=<your-key>")
            sys.exit(1)

    if args.agent_api == "openai" or (args.mode == "simulated" and args.user_api == "openai"):
        if not os.environ.get("OPENAI_API_KEY"):
            print("Error: OPENAI_API_KEY environment variable not set")
            sys.exit(1)

    if args.agent_api == "anthropic" or (args.mode == "simulated" and args.user_api == "anthropic"):
        if not os.environ.get("ANTHROPIC_API_KEY"):
            print("Error: ANTHROPIC_API_KEY environment variable not set")
            sys.exit(1)

    # Get scenario config
    scenario_config = get_scenario_config(args.scenario)

    # Override task description if provided
    if args.task_description:
        scenario_config["task_description"] = args.task_description

    # Create experiment name
    if args.experiment_name:
        experiment_name = args.experiment_name
    else:
        import time
        experiment_name = f"exp_{args.scenario}_{args.human_profile}_{time.strftime('%Y%m%d_%H%M%S')}"

    # Create experiment config
    # Memory persists across experiments in a shared directory
    memory_dir = REPO_ROOT / "memory"

    # Determine user_id based on mode
    if args.mode == "scripted":
        user_id = "scripted_user"
    elif args.mode == "simulated":
        user_id = "simulated_user"
    else:
        user_id = "interactive_user"

    config = ExperimentConfig(
        experiment_name=experiment_name,
        max_turns=args.max_turns,
        task_description=scenario_config["task_description"],
        data_files=scenario_config.get("data_files", {}),
        human_profile=args.human_profile,
        human_beliefs=args.human_beliefs,
        memory_dir=memory_dir,
        session_id=experiment_name,
        user_id=user_id,
        # Scripted mode settings
        use_scripted_human=(args.mode == "scripted"),
        script_name=scenario_config.get("script_name", "default"),
    )

    # Create LLM call functions
    agent_llm_call = create_llm_call(args.agent_api, args.agent_llm)

    if args.mode == "simulated":
        user_llm_call = create_llm_call(args.user_api, args.user_llm)
    else:
        user_llm_call = None  # Not needed for interactive or scripted mode

    print("\n" + "=" * 70)
    print("   COLLABORATIVE AGENT EXPERIMENT")
    print("=" * 70)
    print(f"[CONFIG] Mode: {args.mode}")
    print(f"[CONFIG] Scenario: {args.scenario}")
    print(f"[CONFIG] Task: {scenario_config['task_description'][:60]}...")
    print(f"[CONFIG] Experiment name: {experiment_name}")
    print("-" * 70)
    print(f"[LLM] Agent: {args.agent_llm} (API: {args.agent_api})")
    if args.mode == "simulated":
        print(f"[LLM] Simulated User: {args.user_llm} (API: {args.user_api})")
        print(f"[HUMAN SIM] Profile: {args.human_profile}")
        print(f"[HUMAN SIM] Hidden beliefs: {args.human_beliefs}")
    elif args.mode == "scripted":
        print(f"[SCRIPTED] Script: {scenario_config.get('script_name', 'default')}")
        print(f"[SCRIPTED] No LLM needed for user responses")
    print("-" * 70)
    print(f"[SETTINGS] Max turns: {args.max_turns}")
    print(f"[SETTINGS] Memory dir: {memory_dir}")
    print(f"[SETTINGS] Output dir: {args.output_dir}")
    print("=" * 70)
    print("\n[STARTUP] Creating LLM call functions...")
    print(f"[STARTUP] ✓ Agent LLM ready")

    # Run experiment
    if args.mode == "interactive":
        print(f"[STARTUP] ✓ Interactive mode - no simulated user needed")
        print("\n[STARTUP] Creating InteractiveEnvironment...")
        env = InteractiveEnvironment(
            llm_call_fn=agent_llm_call,
            repo_root=REPO_ROOT,
            config=config,
        )
        print("[STARTUP] ✓ Environment ready, starting interactive session...")
        results = env.run_interactive()
    elif args.mode == "scripted":
        print(f"[STARTUP] ✓ Scripted mode - using predefined responses")
        print("\n[STARTUP] Creating CollaborationEnvironment with ScriptedHuman...")
        env = CollaborationEnvironment(
            llm_call_fn=agent_llm_call,
            repo_root=REPO_ROOT,
            config=config,
            user_llm_call_fn=None,  # Not needed for scripted
        )
        print("[STARTUP] ✓ Environment ready, starting experiment...")
        results = env.run_experiment()
    else:  # simulated
        print(f"[STARTUP] ✓ Simulated User LLM ready")
        print("\n[STARTUP] Creating CollaborationEnvironment...")
        env = CollaborationEnvironment(
            llm_call_fn=agent_llm_call,
            repo_root=REPO_ROOT,
            config=config,
            user_llm_call_fn=user_llm_call,
        )
        print("[STARTUP] ✓ Environment ready, starting experiment...")
        results = env.run_experiment()

    # Add LLM config to results
    results["llm_config"] = {
        "agent_llm": args.agent_llm,
        "agent_api": args.agent_api,
        "user_llm": args.user_llm if args.mode == "simulated" else None,
        "user_api": args.user_api if args.mode == "simulated" else None,
    }

    # Save session summary to memory
    try:
        memory_path = env.agent.save_session_summary()
        if memory_path:
            print(f"\nSession summary saved to memory: {memory_path}")
    except Exception as e:
        print(f"\nWarning: Failed to save session summary: {e}")

    # Save results
    output_dir = REPO_ROOT / args.output_dir
    output_path = env.save_results(results, output_dir)

    print(f"\n{'='*60}")
    print("EXPERIMENT COMPLETE")
    print(f"{'='*60}")
    print(f"Results saved to: {output_path}")
    print(f"\nKey metrics:")
    print(f"  Total turns: {results['metrics']['total_turns']}")
    print(f"  Questions asked: {results['metrics']['agent_questions_asked']}")
    print(f"  Plans proposed: {results['metrics']['agent_plans_proposed']}")
    print(f"  Tools used: {len(results['metrics']['agent_tools_used'])}")

    if args.mode == "simulated":
        print(f"  Final satisfaction: {results['metrics'].get('final_human_satisfaction', 'N/A'):.2f}")
        print(f"  Final goal progress: {results['metrics'].get('final_goal_progress', 'N/A'):.2f}")


if __name__ == "__main__":
    main()
