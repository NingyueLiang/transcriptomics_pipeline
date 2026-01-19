"""
Collaborative Agent - Human-AI Research Collaboration System.

This module implements a general-purpose system where an AI agent can actively
collaborate with humans (or simulated humans) on research tasks across any domain.

Research Questions:
- How does the agent form the notion of "collaboration"?
- How does the agent mentally model the user?
- How does the agent accumulate memory?
"""

from .agent import CollaborativeAgent, AgentState, ConversationTurn, COLLABORATION_ACTIONS
from .simulated_human import SimulatedHuman, ScriptedHuman, HumanState
from .environment import CollaborationEnvironment, InteractiveEnvironment, ExperimentConfig
from .tools import ToolRegistry, Tool
from .memory_system import MemorySystem, MemoryEntry, create_memory_tools
from .scenarios import Scenario, get_scenario, list_scenarios

__all__ = [
    # Core components
    "CollaborativeAgent",
    "AgentState",
    "ConversationTurn",
    "COLLABORATION_ACTIONS",
    "SimulatedHuman",
    "ScriptedHuman",
    "HumanState",
    "CollaborationEnvironment",
    "InteractiveEnvironment",
    "ExperimentConfig",
    # Tools
    "ToolRegistry",
    "Tool",
    # Memory
    "MemorySystem",
    "MemoryEntry",
    "create_memory_tools",
    # Scenarios
    "Scenario",
    "get_scenario",
    "list_scenarios",
]
