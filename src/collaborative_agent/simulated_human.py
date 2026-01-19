"""
Simulated Human for testing human-AI collaboration.

Two types of simulated humans:
1. ScriptedHuman - Uses predefined responses (for testing/development)
2. SimulatedHuman - Uses LLM to generate responses (for research)
"""

from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from .prompts import (
    SIMULATED_HUMAN_SYSTEM_PROMPT,
    USER_PROFILES,
    HIDDEN_BELIEFS_TEMPLATES,
)


# =============================================================================
# SHARED STATE CLASS
# =============================================================================

@dataclass
class HumanState:
    """Tracks the simulated human's internal state."""
    satisfaction: float = 0.5
    confusion: float = 0.0
    goal_progress: float = 0.0
    turn_count: int = 0
    thoughts_history: List[str] = field(default_factory=list)


# =============================================================================
# SCRIPTED HUMAN (for simple testing)
# =============================================================================

# Predefined conversation scripts for testing
CONVERSATION_SCRIPTS = {
    "simple_hello": [
        "Hello! I'm here to test the system.",
        "That sounds good. Can you tell me more?",
        "Interesting. What else can you do?",
        "Thanks for the explanation. That's all I needed.",
    ],
    "memory_test": [
        "Hi! My name is Alice and I'm a biologist.",
        "Can you remember my name?",
        "What's my profession?",
        "Great! I prefer detailed explanations over brief ones.",
        "Do you remember my preferences?",
        "Thanks! That's all for now.",
    ],
    "default": [
        "Hello, I need help with something.",
        "That makes sense. Please continue.",
        "I understand. What's next?",
        "Good. Is there anything else?",
        "Thank you. That's helpful.",
    ],
}


class ScriptedHuman:
    """
    A scripted human for simple testing.

    Uses predefined responses instead of LLM calls.
    Good for testing memory system without LLM costs/latency.
    """

    def __init__(self, script_name: str = "default"):
        """
        Initialize with a conversation script.

        Args:
            script_name: Name of predefined script to use.
        """
        self.script = CONVERSATION_SCRIPTS.get(script_name, CONVERSATION_SCRIPTS["default"])
        self.turn_index = 0
        self.state = HumanState()
        print(f"[SCRIPTED HUMAN] Using script: {script_name} ({len(self.script)} messages)")

    def respond(self, agent_message: str, agent_action: str = "") -> Dict[str, Any]:
        """Return the next scripted response."""
        self.state.turn_count += 1
        self.turn_index += 1

        # Get next response from script (loop if exceeded)
        if self.turn_index < len(self.script):
            response_text = self.script[self.turn_index]
        else:
            response_text = "Thanks, I think we're done here."
            self.state.goal_progress = 1.0

        # Simulate progress
        self.state.goal_progress = min(1.0, self.turn_index / len(self.script))
        self.state.satisfaction = 0.7

        print(f"[SCRIPTED HUMAN] Turn {self.state.turn_count}: {response_text[:50]}...")

        return {
            "internal_state": {
                "satisfaction": self.state.satisfaction,
                "confusion": 0.0,
                "goal_progress": self.state.goal_progress,
                "thoughts": "Scripted response",
            },
            "response_to_agent": response_text,
        }

    def get_initial_message(self) -> str:
        """Return the first message from the script."""
        initial = self.script[0] if self.script else "Hello!"
        print(f"[SCRIPTED HUMAN] Initial message: {initial}")
        return initial

    def get_state_summary(self) -> Dict[str, Any]:
        """Get state summary."""
        return {
            "turn_count": self.state.turn_count,
            "satisfaction": self.state.satisfaction,
            "confusion": self.state.confusion,
            "goal_progress": self.state.goal_progress,
            "profile": "scripted",
        }

    def get_internal_log(self) -> List[str]:
        """Get internal thoughts log."""
        return ["Scripted human - no internal thoughts"]


# =============================================================================
# LLM-BASED SIMULATED HUMAN (for research)
# =============================================================================


class SimulatedHuman:
    """
    A simulated human for testing the collaborative agent.

    The simulated human has:
    - A profile defining their expertise level
    - Hidden beliefs/goals that influence their behavior
    - Internal state tracking (satisfaction, confusion, goal progress)
    """

    def __init__(
        self,
        llm_call_fn,
        profile_name: str = "expert_biologist",
        hidden_beliefs_name: str = "exploratory",
        task_description: str = "",
        custom_profile: Optional[Dict[str, str]] = None,
        custom_beliefs: Optional[str] = None,
    ):
        """
        Initialize the simulated human.

        Args:
            llm_call_fn: Function to call the LLM. Signature: (prompt, system_prompt) -> str
            profile_name: Name of predefined profile to use.
            hidden_beliefs_name: Name of predefined hidden beliefs to use.
            task_description: Description of the task the human wants to accomplish.
            custom_profile: Custom profile dict (overrides profile_name).
            custom_beliefs: Custom hidden beliefs string (overrides hidden_beliefs_name).
        """
        self.llm_call = llm_call_fn
        self.task_description = task_description

        # Set up profile
        if custom_profile:
            self.profile = custom_profile
        elif profile_name in USER_PROFILES:
            self.profile = USER_PROFILES[profile_name]
        else:
            self.profile = USER_PROFILES["expert_biologist"]

        # Set up hidden beliefs
        if custom_beliefs:
            self.hidden_beliefs = custom_beliefs
        elif hidden_beliefs_name in HIDDEN_BELIEFS_TEMPLATES:
            self.hidden_beliefs = HIDDEN_BELIEFS_TEMPLATES[hidden_beliefs_name]["beliefs"]
        else:
            self.hidden_beliefs = HIDDEN_BELIEFS_TEMPLATES["exploratory"]["beliefs"]

        # State tracking
        self.state = HumanState()
        self.conversation_history: List[Dict[str, Any]] = []

    def _build_system_prompt(self) -> str:
        """Build the system prompt for the simulated human."""
        history_str = self._format_conversation_history()

        return SIMULATED_HUMAN_SYSTEM_PROMPT.format(
            user_profile=self.profile.get("profile", ""),
            hidden_beliefs=self.hidden_beliefs,
            knowledge_level=self.profile.get("knowledge_level", ""),
            communication_style=self.profile.get("communication_style", ""),
            task_description=self.task_description or "General transcriptomics analysis",
            conversation_history=history_str or "Conversation just started.",
        )

    def _format_conversation_history(self, max_turns: int = 10) -> str:
        """Format recent conversation history."""
        if not self.conversation_history:
            return ""

        recent = self.conversation_history[-max_turns:]
        lines = []

        for turn in recent:
            role = turn.get("role", "")
            content = turn.get("content", "")
            lines.append(f"[{role.capitalize()}]: {content}")

        return "\n".join(lines)

    def _parse_response(self, response: str) -> Dict[str, Any]:
        """Parse the simulated human's JSON response."""
        try:
            return json.loads(response)
        except json.JSONDecodeError:
            pass

        json_match = re.search(r'\{[\s\S]*\}', response)
        if json_match:
            try:
                return json.loads(json_match.group())
            except json.JSONDecodeError:
                pass

        # Fallback
        return {
            "internal_state": {
                "satisfaction": 0.5,
                "confusion": 0.5,
                "goal_progress": self.state.goal_progress,
                "thoughts": "Response parsing failed",
            },
            "response_to_agent": response,
        }

    def respond(self, agent_message: str, agent_action: str = "") -> Dict[str, Any]:
        """
        Generate a response to the agent's message.

        Args:
            agent_message: The agent's message to respond to.
            agent_action: The type of action the agent took.

        Returns:
            Dict with internal_state and response_to_agent.
        """
        self.state.turn_count += 1

        # Record agent's message
        self.conversation_history.append({
            "role": "agent",
            "action": agent_action,
            "content": agent_message,
        })

        # Build prompt
        system_prompt = self._build_system_prompt()

        user_prompt = f"""The agent said: "{agent_message}"

The agent's action type was: {agent_action}

Generate your response as the simulated human. Remember to include your internal state (satisfaction, confusion, goal_progress, thoughts) and your response to the agent."""

        # Call LLM
        raw_response = self.llm_call(user_prompt, system_prompt)

        # Parse response
        response = self._parse_response(raw_response)

        # Update state
        internal = response.get("internal_state", {})
        self.state.satisfaction = internal.get("satisfaction", self.state.satisfaction)
        self.state.confusion = internal.get("confusion", self.state.confusion)
        self.state.goal_progress = internal.get("goal_progress", self.state.goal_progress)

        thoughts = internal.get("thoughts", "")
        if thoughts:
            self.state.thoughts_history.append(thoughts)

        # Record human's response
        self.conversation_history.append({
            "role": "human",
            "content": response.get("response_to_agent", ""),
        })

        return response

    def get_initial_message(self) -> str:
        """Generate an initial message to start the conversation."""
        system_prompt = self._build_system_prompt()

        user_prompt = """You are starting a new conversation with an AI research assistant.
Generate your opening message to initiate the collaboration. This should introduce your task/question.

Respond with JSON including your internal_state and response_to_agent (which will be your opening message)."""

        raw_response = self.llm_call(user_prompt, system_prompt)
        response = self._parse_response(raw_response)

        # Record the initial message
        initial_msg = response.get("response_to_agent", "")
        self.conversation_history.append({
            "role": "human",
            "content": initial_msg,
        })

        return initial_msg

    def get_state_summary(self) -> Dict[str, Any]:
        """Get a summary of the human's current state."""
        return {
            "turn_count": self.state.turn_count,
            "satisfaction": self.state.satisfaction,
            "confusion": self.state.confusion,
            "goal_progress": self.state.goal_progress,
            "profile": self.profile.get("profile", ""),
        }

    def get_internal_log(self) -> List[str]:
        """Get the history of internal thoughts (for research analysis)."""
        return self.state.thoughts_history.copy()
