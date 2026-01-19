"""
Collaboration Environment for running human-AI collaboration experiments.

This module orchestrates the interaction between:
- CollaborativeAgent: The AI agent that responds to human input
- SimulatedHuman (or real human): Provides input to the agent

Execution Flow (simulated experiment):
1. Environment initializes agent and simulated human
2. Human sends initial message
3. Loop for max_turns:
   a. Agent responds to human message
   b. Check end conditions (goal achieved, satisfaction too low)
   c. Human responds to agent message
   d. Record conversation and metrics
4. Save results to file
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from .agent import CollaborativeAgent
from .simulated_human import SimulatedHuman, ScriptedHuman


@dataclass
class ExperimentConfig:
    """Configuration for a collaboration experiment."""
    experiment_name: str
    max_turns: int = 20
    task_description: str = ""
    data_files: Dict[str, str] = field(default_factory=dict)
    human_profile: str = "domain_expert"
    human_beliefs: str = "exploratory"
    custom_human_profile: Optional[Dict[str, str]] = None
    custom_human_beliefs: Optional[str] = None
    # Memory settings
    memory_dir: Optional[Path] = None
    session_id: Optional[str] = None
    user_id: Optional[str] = None
    # Testing mode
    use_scripted_human: bool = False  # If True, use ScriptedHuman instead of LLM
    script_name: str = "default"      # Script to use for ScriptedHuman


@dataclass
class ExperimentMetrics:
    """Metrics collected during an experiment."""
    total_turns: int = 0
    agent_questions_asked: int = 0
    agent_plans_proposed: int = 0
    agent_tools_used: List[str] = field(default_factory=list)
    human_satisfaction_history: List[float] = field(default_factory=list)
    human_confusion_history: List[float] = field(default_factory=list)
    human_goal_progress_history: List[float] = field(default_factory=list)
    turn_timestamps: List[float] = field(default_factory=list)
    tool_execution_results: List[Dict[str, Any]] = field(default_factory=list)


class CollaborationEnvironment:
    """
    Environment for running human-AI collaboration experiments.

    Orchestrates the interaction between a CollaborativeAgent and a
    SimulatedHuman (or real human via callbacks).
    """

    def __init__(
        self,
        llm_call_fn: Callable,
        repo_root: Path,
        config: ExperimentConfig,
        user_llm_call_fn: Optional[Callable] = None,
    ):
        """
        Initialize the collaboration environment.

        Args:
            llm_call_fn: Function to call the LLM for the agent.
            repo_root: Path to the repository root.
            config: Experiment configuration.
            user_llm_call_fn: Function to call the LLM for simulated user.
                              If None, uses llm_call_fn for both.
        """
        print("\n" + "=" * 60)
        print("[ENV INIT] Initializing CollaborationEnvironment...")
        print("=" * 60)
        print(f"[ENV INIT] Experiment: {config.experiment_name}")
        print(f"[ENV INIT] Task: {config.task_description[:80]}{'...' if len(config.task_description) > 80 else ''}")
        print(f"[ENV INIT] Max turns: {config.max_turns}")

        self.agent_llm_call = llm_call_fn
        self.user_llm_call = user_llm_call_fn or llm_call_fn
        self.repo_root = repo_root
        self.config = config

        # Build task context from config
        task_context = self._build_task_context()

        # Initialize agent (uses agent LLM)
        print("\n[ENV INIT] Creating CollaborativeAgent...")
        self.agent = CollaborativeAgent(
            llm_call_fn=self.agent_llm_call,
            repo_root=repo_root,
            task_context=task_context,
            memory_dir=config.memory_dir,
            session_id=config.session_id or config.experiment_name,
            user_id=config.user_id,
        )

        # Initialize human (scripted or LLM-simulated)
        if config.use_scripted_human:
            print(f"[ENV INIT] Creating ScriptedHuman (script={config.script_name})...")
            self.human = ScriptedHuman(script_name=config.script_name)
            print("[ENV INIT] ✓ ScriptedHuman created (no LLM needed for user)")
        else:
            print(f"[ENV INIT] Creating SimulatedHuman (profile={config.human_profile}, beliefs={config.human_beliefs})...")
            self.human = SimulatedHuman(
                llm_call_fn=self.user_llm_call,
                profile_name=config.human_profile,
                hidden_beliefs_name=config.human_beliefs,
                task_description=config.task_description,
                custom_profile=config.custom_human_profile,
                custom_beliefs=config.custom_human_beliefs,
            )
            print("[ENV INIT] ✓ SimulatedHuman created")

        # Metrics for tracking experiment progress
        self.metrics = ExperimentMetrics()

        # Conversation log for analysis
        self.conversation_log: List[Dict[str, Any]] = []

        # Callbacks for real human interaction
        self.human_input_callback: Optional[Callable[[str], str]] = None
        self.display_callback: Optional[Callable[[str, str], None]] = None

        print("[ENV INIT] ✓ Environment ready")
        print("=" * 60 + "\n")

    def _build_task_context(self) -> str:
        """Build task context from config."""
        context_parts = [f"Task: {self.config.task_description}"]

        if self.config.data_files:
            context_parts.append("\nAvailable data files:")
            for name, path in self.config.data_files.items():
                context_parts.append(f"- {name}: {path}")

        return "\n".join(context_parts)

    def set_human_input_callback(self, callback: Callable[[str], str]) -> None:
        """
        Set a callback for getting real human input.

        When set, the environment will use this instead of the simulated human.

        Args:
            callback: Function that takes agent's message and returns human's response.
        """
        self.human_input_callback = callback

    def set_display_callback(self, callback: Callable[[str, str], None]) -> None:
        """
        Set a callback for displaying messages.

        Args:
            callback: Function that takes (role, message) and displays it.
        """
        self.display_callback = callback

    def _display(self, role: str, message: str) -> None:
        """Display a message."""
        if self.display_callback:
            self.display_callback(role, message)
        else:
            print(f"\n[{role.upper()}]: {message}")

    def _get_human_response(self, agent_message: str, agent_action: str) -> Dict[str, Any]:
        """Get response from human (real or simulated)."""
        if self.human_input_callback:
            # Real human
            response_text = self.human_input_callback(agent_message)
            return {
                "internal_state": {
                    "satisfaction": None,
                    "confusion": None,
                    "goal_progress": None,
                    "thoughts": "Real human - no simulated state",
                },
                "response_to_agent": response_text,
            }
        else:
            # Simulated human
            return self.human.respond(agent_message, agent_action)

    def run_turn(self, human_message: str) -> Dict[str, Any]:
        """
        Run a single turn of interaction.

        Args:
            human_message: The human's message to the agent.

        Returns:
            Dict with agent_response and human_response details.
        """
        turn_start = time.time()
        self.metrics.total_turns += 1
        self.metrics.turn_timestamps.append(turn_start)

        # Agent responds to human
        agent_response = self.agent.respond(human_message)
        agent_message = agent_response.get("message_to_user", "")
        agent_action = agent_response.get("action", "")

        # Log agent turn
        self.conversation_log.append({
            "turn": self.metrics.total_turns,
            "role": "agent",
            "action": agent_action,
            "message": agent_message,
            "thinking": agent_response.get("thinking", ""),
            "action_input": agent_response.get("action_input", {}),
            "timestamp": turn_start,
        })

        # Display agent response
        self._display("Agent", agent_message)

        # Handle tool execution if needed
        tool_result = None
        if agent_action == "EXECUTE_TOOL":
            tool_info = agent_response.get("action_input", {})
            tool_name = tool_info.get("tool", "")
            tool_params = tool_info.get("params", {})

            # Check if confirmation is needed
            if agent_response.get("awaiting_confirmation"):
                self._display("System", f"Tool '{tool_name}' requires confirmation.")
            else:
                # Execute tool directly
                tool_result = self.agent.tools.execute(tool_name, tool_params)
                self.metrics.tool_execution_results.append({
                    "tool": tool_name,
                    "params": tool_params,
                    "result": tool_result,
                })
                self._display("System", f"Tool '{tool_name}' executed. Result: {json.dumps(tool_result, indent=2)[:500]}...")

        # Update agent metrics
        self.metrics.agent_questions_asked = self.agent.state.clarifications_asked
        self.metrics.agent_plans_proposed = self.agent.state.plans_proposed
        self.metrics.agent_tools_used = self.agent.state.tools_executed.copy()

        return {
            "agent_response": agent_response,
            "agent_message": agent_message,
            "agent_action": agent_action,
            "tool_result": tool_result,
        }

    def run_experiment(self) -> Dict[str, Any]:
        """
        Run a full experiment with simulated human.

        This is the main experiment loop:
        1. Get initial message from human
        2. For each turn: agent responds, check end conditions, human responds
        3. Compile and return results

        Returns:
            Dict containing experiment results and metrics.
        """
        print("\n" + "#" * 60)
        print("[EXPERIMENT] Starting simulated experiment...")
        print("#" * 60)

        experiment_start = time.time()

        # Get initial message from human
        print("\n[EXPERIMENT] Getting initial message from simulated human...")
        human_message = self.human.get_initial_message()
        self._display("Human", human_message)

        self.conversation_log.append({
            "turn": 0,
            "role": "human",
            "message": human_message,
            "internal_state": {
                "satisfaction": self.human.state.satisfaction,
                "confusion": self.human.state.confusion,
                "goal_progress": self.human.state.goal_progress,
            },
            "timestamp": experiment_start,
        })

        # Main interaction loop
        print(f"\n[EXPERIMENT] Starting main loop (max {self.config.max_turns} turns)...")
        for turn in range(self.config.max_turns):
            print(f"\n{'='*60}")
            print(f"[EXPERIMENT] === TURN {turn + 1}/{self.config.max_turns} ===")
            print(f"{'='*60}")

            # Agent turn
            turn_result = self.run_turn(human_message)

            # Record human state
            self.metrics.human_satisfaction_history.append(self.human.state.satisfaction)
            self.metrics.human_confusion_history.append(self.human.state.confusion)
            self.metrics.human_goal_progress_history.append(self.human.state.goal_progress)

            # Check for experiment end conditions
            if self.human.state.goal_progress >= 0.95:
                print(f"\n[EXPERIMENT] ★ Goal achieved (progress={self.human.state.goal_progress:.2f})! Ending experiment.")
                self._display("System", "Goal achieved! Ending experiment.")
                break

            if self.human.state.satisfaction < 0.1:
                print(f"\n[EXPERIMENT] ✗ Human satisfaction too low ({self.human.state.satisfaction:.2f}). Ending experiment.")
                self._display("System", "Human satisfaction too low. Ending experiment.")
                break

            # Human responds
            print(f"\n[EXPERIMENT] Getting response from simulated human...")
            human_response = self._get_human_response(
                turn_result["agent_message"],
                turn_result["agent_action"],
            )

            human_message = human_response.get("response_to_agent", "")
            internal = human_response.get("internal_state", {})
            print(f"[EXPERIMENT] Human internal state: satisfaction={internal.get('satisfaction', 'N/A'):.2f}, "
                  f"confusion={internal.get('confusion', 'N/A'):.2f}, "
                  f"goal_progress={internal.get('goal_progress', 'N/A'):.2f}")
            self._display("Human", human_message)

            # Record human turn
            self.conversation_log.append({
                "turn": turn + 1,
                "role": "human",
                "message": human_message,
                "internal_state": human_response.get("internal_state", {}),
                "timestamp": time.time(),
            })

            # Record in agent's history
            self.agent.record_human_response(
                human_message,
                human_response.get("internal_state"),
            )

        experiment_end = time.time()
        print(f"\n[EXPERIMENT] Experiment completed in {experiment_end - experiment_start:.2f} seconds")

        # Compile results
        results = {
            "experiment_name": self.config.experiment_name,
            "config": {
                "max_turns": self.config.max_turns,
                "task_description": self.config.task_description,
                "human_profile": self.config.human_profile,
                "human_beliefs": self.config.human_beliefs,
            },
            "metrics": {
                "total_turns": self.metrics.total_turns,
                "agent_questions_asked": self.metrics.agent_questions_asked,
                "agent_plans_proposed": self.metrics.agent_plans_proposed,
                "agent_tools_used": self.metrics.agent_tools_used,
                "final_human_satisfaction": self.human.state.satisfaction,
                "final_human_confusion": self.human.state.confusion,
                "final_goal_progress": self.human.state.goal_progress,
                "human_satisfaction_history": self.metrics.human_satisfaction_history,
                "human_confusion_history": self.metrics.human_confusion_history,
                "human_goal_progress_history": self.metrics.human_goal_progress_history,
                "duration_seconds": experiment_end - experiment_start,
            },
            "conversation_log": self.conversation_log,
            "human_internal_thoughts": self.human.get_internal_log(),
            "agent_state": self.agent.get_state_summary(),
        }

        return results

    def save_results(self, results: Dict[str, Any], output_dir: Path) -> Path:
        """
        Save experiment results to a file.

        Args:
            results: The experiment results dict.
            output_dir: Directory to save results.

        Returns:
            Path to the saved file.
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        timestamp = time.strftime("%Y%m%d_%H%M%S")
        filename = f"{self.config.experiment_name}_{timestamp}.json"
        output_path = output_dir / filename

        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)

        return output_path


class InteractiveEnvironment(CollaborationEnvironment):
    """
    Environment for interactive (real human) collaboration.

    Provides a simple CLI interface for human-agent interaction.
    """

    def __init__(
        self,
        llm_call_fn: Callable,
        repo_root: Path,
        config: ExperimentConfig,
    ):
        super().__init__(llm_call_fn, repo_root, config)

        # Set up real human input
        self.set_human_input_callback(self._cli_input)

    def _cli_input(self, agent_message: str) -> str:
        """Get input from CLI."""
        return input("\nYour response: ").strip()

    def run_interactive(self) -> Dict[str, Any]:
        """Run an interactive session with a real human."""
        print("\n" + "=" * 60)
        print("COLLABORATIVE AGENT - Interactive Session")
        print("=" * 60)
        print(f"\nTask: {self.config.task_description}")
        print("\nType 'quit' or 'exit' to end the session.")
        print("Type 'status' to see current metrics.")
        print("-" * 60)

        # Get initial message from user
        human_message = input("\nYou: ").strip()

        while human_message.lower() not in ['quit', 'exit']:
            if human_message.lower() == 'status':
                print(f"\n[STATUS] Turns: {self.metrics.total_turns}, "
                      f"Questions asked: {self.metrics.agent_questions_asked}, "
                      f"Tools used: {len(self.metrics.agent_tools_used)}")
                human_message = input("\nYou: ").strip()
                continue

            if not human_message:
                human_message = input("\nYou: ").strip()
                continue

            # Run turn
            self.run_turn(human_message)

            # Get next human input
            human_message = input("\nYou: ").strip()

        print("\n" + "=" * 60)
        print("Session ended.")

        return {
            "experiment_name": self.config.experiment_name,
            "metrics": {
                "total_turns": self.metrics.total_turns,
                "agent_questions_asked": self.metrics.agent_questions_asked,
                "agent_plans_proposed": self.metrics.agent_plans_proposed,
                "agent_tools_used": self.metrics.agent_tools_used,
            },
            "conversation_log": self.conversation_log,
            "agent_state": self.agent.get_state_summary(),
        }
