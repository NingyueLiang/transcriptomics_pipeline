"""
Collaborative Agent for human-AI research collaboration.

The agent operates as a teammate, not just a tool:
- Asks clarifying questions when uncertain
- Proposes plans for complex tasks
- Challenges assumptions when appropriate
- Reflects on failures to learn
- Uses memory to persist learning across sessions

Execution Flow:
1. Agent receives user input via respond()
2. Agent builds system prompt with task context and conversation history
3. Agent calls LLM to decide on action (one of COLLABORATION_ACTIONS)
4. Agent parses LLM response and validates action
5. If action is EXECUTE_TOOL, agent executes the tool
6. Agent updates internal state and returns response
"""

from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Callable

from .prompts import COLLABORATIVE_AGENT_SYSTEM_PROMPT
from .tools import ToolRegistry
from .memory_system import MemorySystem, create_memory_tools


# =============================================================================
# COLLABORATION ACTIONS
# These are the actions the agent can take during collaboration.
# The LLM decides which action to take based on the conversation context.
# =============================================================================
COLLABORATION_ACTIONS = [
    "ASK_CLARIFICATION",  # Ask user for more information
    "PROPOSE_PLAN",       # Propose a plan for user approval
    "CHALLENGE",          # Challenge user's assumption
    "REFLECT",            # Reflect on what happened (usually after failure)
    "SUMMARIZE",          # Summarize current understanding
    "DEFER",              # Admit uncertainty, ask for help
    "EXECUTE_TOOL",       # Execute a tool to take action
]


@dataclass
class AgentState:
    """Tracks the agent's internal state for research analysis."""
    turn_count: int = 0
    actions_taken: List[Dict[str, Any]] = field(default_factory=list)
    clarifications_asked: int = 0
    plans_proposed: int = 0
    challenges_made: int = 0
    reflections_made: int = 0
    tools_executed: List[str] = field(default_factory=list)
    memory_writes: int = 0
    memory_reads: int = 0


@dataclass
class ConversationTurn:
    """A single turn in the conversation."""
    turn_id: int
    role: str  # "agent" or "human"
    action: str
    content: Dict[str, Any]
    timestamp: float
    reasoning: Optional[str] = None  # Agent's thinking (for analysis)


class CollaborativeAgent:
    """
    An agent that collaborates with humans on research tasks.

    Key features:
    - Multiple collaboration actions (ask, plan, challenge, reflect, summarize, defer, execute)
    - Persistent memory system
    - Full logging for research analysis
    """

    def __init__(
        self,
        llm_call_fn: Callable[[str, str], str],
        repo_root: Path,
        task_context: str = "",
        memory_dir: Optional[Path] = None,
        session_id: Optional[str] = None,
        user_id: Optional[str] = None,
    ):
        """
        Initialize the collaborative agent.

        Args:
            llm_call_fn: Function to call the LLM. Signature: (prompt, system_prompt) -> str
            repo_root: Path to the repository root.
            task_context: Description of the current task.
            memory_dir: Directory for persistent memory (default: repo_root/memory)
            session_id: Identifier for this session
            user_id: Identifier for the user
        """
        print("\n" + "=" * 60)
        print("[AGENT INIT] Initializing CollaborativeAgent...")
        print("=" * 60)

        self.llm_call = llm_call_fn
        self.repo_root = repo_root
        self.task_context = task_context

        # Initialize memory system
        # Memory persists across sessions in markdown files
        self.memory_dir = memory_dir or (repo_root / "memory")
        self.session_id = session_id or f"session_{int(time.time())}"
        self.user_id = user_id or "default_user"

        print(f"[AGENT INIT] Memory directory: {self.memory_dir}")
        print(f"[AGENT INIT] Session ID: {self.session_id}")
        print(f"[AGENT INIT] User ID: {self.user_id}")

        self.memory = MemorySystem(
            memory_dir=self.memory_dir,
            session_id=self.session_id,
            user_id=self.user_id,
        )
        print("[AGENT INIT] ✓ Memory system initialized")

        # Initialize tools (includes memory tools)
        # Tools available: read_file, list_directory, run_python, search_files,
        #                  add_to_memory, search_memory
        self.tools = ToolRegistry(repo_root)
        self._register_memory_tools()
        print(f"[AGENT INIT] ✓ Tools registered: {list(self.tools.tools.keys())}")

        # State tracking for research analysis
        self.state = AgentState()
        self.conversation_history: List[ConversationTurn] = []
        print("[AGENT INIT] ✓ Agent ready for collaboration")
        print("=" * 60 + "\n")

    def _register_memory_tools(self) -> None:
        """Register memory tools with the tool registry."""
        from .tools import Tool

        memory_tools = create_memory_tools(self.memory)

        # Register add_to_memory
        self.tools.register(Tool(
            name="add_to_memory",
            description="Store information in persistent memory. Categories: 'users' (user info), 'episodes' (what happened), 'skills' (how to do things), 'outcomes' (successes/failures)",
            parameters={
                "category": {"type": "string", "description": "One of: users, episodes, skills, outcomes"},
                "content": {"type": "string", "description": "Information to remember (markdown)"},
                "topic": {"type": "string", "description": "Brief topic for the entry"},
                "subcategory": {"type": "string", "description": "Optional: for skills (skill name) or outcomes (successes/failures)"},
                "tags": {"type": "string", "description": "Optional: comma-separated tags"},
            },
            function=memory_tools["add_to_memory"],
            requires_confirmation=False,
            category="memory",
        ))

        # Register search_memory
        self.tools.register(Tool(
            name="search_memory",
            description="Search memory for relevant information. Use this to recall past interactions, user preferences, lessons learned, etc.",
            parameters={
                "query": {"type": "string", "description": "Search term (case-insensitive)"},
                "category": {"type": "string", "description": "Optional: limit search to category (users, episodes, skills, outcomes)"},
            },
            function=memory_tools["search_memory"],
            requires_confirmation=False,
            category="memory",
        ))

    def set_task_context(self, context: str) -> None:
        """Set or update the task context."""
        self.task_context = context

    def _build_system_prompt(self) -> str:
        """Build the system prompt with current context."""
        history_str = self._format_conversation_history()

        return COLLABORATIVE_AGENT_SYSTEM_PROMPT.format(
            tools_description=self.tools.get_tools_description(),
            task_context=self.task_context or "No specific task context provided.",
            conversation_history=history_str or "No conversation history yet.",
        )

    def _format_conversation_history(self, max_turns: int = 10) -> str:
        """Format recent conversation history for the prompt."""
        if not self.conversation_history:
            return ""

        recent = self.conversation_history[-max_turns:]
        lines = []

        for turn in recent:
            if turn.role == "agent":
                msg = turn.content.get("message_to_user", "")
                action = turn.action
                lines.append(f"[Agent - {action}]: {msg}")
            else:
                msg = turn.content.get("response_to_agent", turn.content.get("message", ""))
                lines.append(f"[Human]: {msg}")

        return "\n".join(lines)

    def _parse_agent_response(self, response: str) -> Dict[str, Any]:
        """Parse the agent's JSON response."""
        # Try to extract JSON from the response
        try:
            return json.loads(response)
        except json.JSONDecodeError:
            pass

        # Try to find JSON in the response
        json_match = re.search(r'\{[\s\S]*\}', response)
        if json_match:
            try:
                return json.loads(json_match.group())
            except json.JSONDecodeError:
                pass

        # Fallback: return as a message
        return {
            "thinking": "Failed to parse response as JSON",
            "action": "DEFER",
            "action_input": {"reason": "Response parsing failed", "what_would_help": "Clearer formatting"},
            "message_to_user": response,
        }

    def respond(self, user_input: str) -> Dict[str, Any]:
        """
        Generate a response to user input.

        This is the main entry point for agent responses. The agent:
        1. Builds a system prompt with context
        2. Calls the LLM to decide on an action
        3. Parses and validates the response
        4. Executes tools if needed
        5. Updates internal state

        Args:
            user_input: The user's message.

        Returns:
            Dict containing the agent's response with action, content, and message.
        """
        self.state.turn_count += 1

        print("\n" + "-" * 60)
        print(f"[AGENT RESPOND] Turn {self.state.turn_count}")
        print("-" * 60)
        print(f"[AGENT RESPOND] Received input: {user_input[:100]}{'...' if len(user_input) > 100 else ''}")

        # Build the prompt with task context and conversation history
        system_prompt = self._build_system_prompt()
        print("[AGENT RESPOND] ✓ Built system prompt with context")

        user_prompt = f"""User's message: {user_input}

Decide what action to take and respond in JSON format with thinking, action, action_input, and message_to_user."""

        # Call the LLM to decide on action
        print("[AGENT RESPOND] Calling LLM for action decision...")
        raw_response = self.llm_call(user_prompt, system_prompt)
        print("[AGENT RESPOND] ✓ LLM responded")

        # Parse the response (extract JSON from LLM output)
        response = self._parse_agent_response(raw_response)

        # Extract and validate action
        action = response.get("action", "DEFER")
        if action not in COLLABORATION_ACTIONS:
            print(f"[AGENT RESPOND] ⚠ Invalid action '{action}', defaulting to DEFER")
            action = "DEFER"
            response["action"] = action

        print(f"[AGENT RESPOND] ★ Action chosen: {action}")
        if response.get("thinking"):
            print(f"[AGENT RESPOND] Reasoning: {response['thinking'][:150]}{'...' if len(response.get('thinking', '')) > 150 else ''}")

        # Record the turn in conversation history
        turn = ConversationTurn(
            turn_id=self.state.turn_count,
            role="agent",
            action=action,
            content=response,
            timestamp=time.time(),
            reasoning=response.get("thinking", ""),
        )
        self.conversation_history.append(turn)

        # Update state based on action (counters for research analysis)
        self._update_state_from_action(response)

        # Handle tool execution if the action is EXECUTE_TOOL
        if action == "EXECUTE_TOOL":
            tool_result = self._handle_tool_execution(response.get("action_input", {}))
            response["tool_result"] = tool_result

        # Record action for analysis
        self.state.actions_taken.append({
            "turn": self.state.turn_count,
            "action": action,
            "reasoning": response.get("thinking", ""),
            "timestamp": time.time(),
        })

        print(f"[AGENT RESPOND] Message to user: {response.get('message_to_user', '')[:100]}{'...' if len(response.get('message_to_user', '')) > 100 else ''}")
        print("-" * 60)

        return response

    def _update_state_from_action(self, response: Dict[str, Any]) -> None:
        """
        Update internal state based on the action taken.
        This tracks metrics for research analysis.
        """
        action = response.get("action", "")

        if action == "ASK_CLARIFICATION":
            self.state.clarifications_asked += 1
            print(f"[AGENT STATE] Clarifications asked: {self.state.clarifications_asked}")
        elif action == "PROPOSE_PLAN":
            self.state.plans_proposed += 1
            print(f"[AGENT STATE] Plans proposed: {self.state.plans_proposed}")
        elif action == "CHALLENGE":
            self.state.challenges_made += 1
            print(f"[AGENT STATE] Challenges made: {self.state.challenges_made}")
        elif action == "REFLECT":
            self.state.reflections_made += 1
            print(f"[AGENT STATE] Reflections made: {self.state.reflections_made}")
        elif action == "EXECUTE_TOOL":
            tool_name = response.get("action_input", {}).get("tool", "")
            self.state.tools_executed.append(tool_name)
            print(f"[AGENT STATE] Tool executed: {tool_name}")

            # Track memory operations specifically
            if tool_name == "add_to_memory":
                self.state.memory_writes += 1
                print(f"[AGENT STATE] Memory writes: {self.state.memory_writes}")
            elif tool_name == "search_memory":
                self.state.memory_reads += 1
                print(f"[AGENT STATE] Memory reads: {self.state.memory_reads}")

    def _handle_tool_execution(self, action_input: Dict[str, Any]) -> Dict[str, Any]:
        """
        Execute a tool and return the result.
        Called when the agent chooses EXECUTE_TOOL action.
        """
        tool_name = action_input.get("tool", "")
        params = action_input.get("params", {})

        if not tool_name:
            print("[AGENT TOOL] ✗ Error: No tool specified")
            return {"error": "No tool specified"}

        print(f"[AGENT TOOL] Executing tool: {tool_name}")
        print(f"[AGENT TOOL] Parameters: {json.dumps(params, indent=2)[:200]}")

        result = self.tools.execute(tool_name, params)

        # Show abbreviated result
        result_str = str(result)
        print(f"[AGENT TOOL] ✓ Result: {result_str[:200]}{'...' if len(result_str) > 200 else ''}")

        return result

    def record_human_response(self, response: str, internal_state: Optional[Dict] = None) -> None:
        """Record a human's response in the conversation history."""
        turn = ConversationTurn(
            turn_id=self.state.turn_count,
            role="human",
            action="response",
            content={
                "response_to_agent": response,
                "internal_state": internal_state or {},
            },
            timestamp=time.time(),
        )
        self.conversation_history.append(turn)

    def get_state_summary(self) -> Dict[str, Any]:
        """Get a summary of the agent's current state."""
        return {
            "turn_count": self.state.turn_count,
            "clarifications_asked": self.state.clarifications_asked,
            "plans_proposed": self.state.plans_proposed,
            "challenges_made": self.state.challenges_made,
            "reflections_made": self.state.reflections_made,
            "tools_executed": self.state.tools_executed,
            "memory_writes": self.state.memory_writes,
            "memory_reads": self.state.memory_reads,
            "conversation_length": len(self.conversation_history),
            "action_distribution": self._get_action_distribution(),
        }

    def _get_action_distribution(self) -> Dict[str, int]:
        """Get distribution of actions taken."""
        distribution = {action: 0 for action in COLLABORATION_ACTIONS}
        for action_record in self.state.actions_taken:
            action = action_record.get("action", "")
            if action in distribution:
                distribution[action] += 1
        return distribution

    def get_conversation_log(self) -> List[Dict[str, Any]]:
        """Get the full conversation log for analysis."""
        return [
            {
                "turn_id": turn.turn_id,
                "role": turn.role,
                "action": turn.action,
                "content": turn.content,
                "reasoning": turn.reasoning,
                "timestamp": turn.timestamp,
            }
            for turn in self.conversation_history
        ]

    def save_session_summary(self) -> str:
        """Save a summary of this session to episodic memory."""
        print("\n[AGENT MEMORY] Saving session summary to episodic memory...")
        summary = f"""## Session Summary

**Session ID**: {self.session_id}
**User ID**: {self.user_id}
**Total Turns**: {self.state.turn_count}

### Actions Taken
- Clarifications asked: {self.state.clarifications_asked}
- Plans proposed: {self.state.plans_proposed}
- Challenges made: {self.state.challenges_made}
- Reflections made: {self.state.reflections_made}
- Tools executed: {len(self.state.tools_executed)}
- Memory writes: {self.state.memory_writes}
- Memory reads: {self.state.memory_reads}

### Action Distribution
{json.dumps(self._get_action_distribution(), indent=2)}

### Task Context
{self.task_context}
"""
        result = self.memory.add_to_memory(
            category="episodes",
            content=summary,
            topic="session_summary",
            tags=["summary", "session_end"],
        )

        file_path = result.get("file_path", "")
        print(f"[AGENT MEMORY] ✓ Session summary saved to: {file_path}")
        return file_path
