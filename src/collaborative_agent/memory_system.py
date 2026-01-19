"""
File-based Memory System for Collaborative Agent.

Memory is stored as human-readable markdown files organized by category.
The agent interacts with memory through two tools:
- add_to_memory(): Store new information
- search_memory(): Retrieve relevant information (grep-like)

Memory Structure:
    memory/
    ├── users/                    # User profiles and observations
    │   └── {user_id}/
    │       └── {timestamp}_{topic}.md
    ├── episodes/                 # Session histories and experiences
    │   └── {session_id}/
    │       └── {timestamp}_{event}.md
    ├── skills/                   # Knowledge about tools and techniques
    │   └── {skill_name}/
    │       └── {timestamp}_{detail}.md
    └── outcomes/                 # Successes and failures
        ├── successes/
        │   └── {timestamp}_{summary}.md
        └── failures/
            └── {timestamp}_{summary}.md
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional


@dataclass
class MemoryEntry:
    """A single memory entry with metadata."""
    content: str
    category: str
    subcategory: Optional[str] = None
    tags: List[str] = field(default_factory=list)
    timestamp: str = ""
    session_id: Optional[str] = None
    user_id: Optional[str] = None
    file_path: Optional[str] = None

    def __post_init__(self):
        if not self.timestamp:
            self.timestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

    def to_markdown(self) -> str:
        """Convert entry to markdown format."""
        lines = [
            f"# Memory Entry",
            f"",
            f"**Timestamp**: {self.timestamp}",
            f"**Category**: {self.category}",
        ]

        if self.subcategory:
            lines.append(f"**Subcategory**: {self.subcategory}")
        if self.session_id:
            lines.append(f"**Session**: {self.session_id}")
        if self.user_id:
            lines.append(f"**User**: {self.user_id}")
        if self.tags:
            lines.append(f"**Tags**: {', '.join(self.tags)}")

        lines.extend([
            f"",
            f"---",
            f"",
            self.content,
        ])

        return "\n".join(lines)

    @classmethod
    def from_markdown(cls, content: str, file_path: str) -> "MemoryEntry":
        """Parse a markdown file into a MemoryEntry."""
        # Extract metadata
        timestamp_match = re.search(r'\*\*Timestamp\*\*:\s*(.+)', content)
        category_match = re.search(r'\*\*Category\*\*:\s*(.+)', content)
        subcategory_match = re.search(r'\*\*Subcategory\*\*:\s*(.+)', content)
        session_match = re.search(r'\*\*Session\*\*:\s*(.+)', content)
        user_match = re.search(r'\*\*User\*\*:\s*(.+)', content)
        tags_match = re.search(r'\*\*Tags\*\*:\s*(.+)', content)

        # Extract content (everything after ---)
        parts = content.split("---", 1)
        main_content = parts[1].strip() if len(parts) > 1 else content

        return cls(
            content=main_content,
            category=category_match.group(1).strip() if category_match else "unknown",
            subcategory=subcategory_match.group(1).strip() if subcategory_match else None,
            tags=[t.strip() for t in tags_match.group(1).split(",")] if tags_match else [],
            timestamp=timestamp_match.group(1).strip() if timestamp_match else "",
            session_id=session_match.group(1).strip() if session_match else None,
            user_id=user_match.group(1).strip() if user_match else None,
            file_path=file_path,
        )


class MemorySystem:
    """
    File-based memory system for the collaborative agent.

    Provides two main operations:
    - add_to_memory(): Store new information as markdown files
    - search_memory(): Grep-like search across memory files
    """

    CATEGORIES = ["users", "episodes", "skills", "outcomes"]
    OUTCOME_SUBCATEGORIES = ["successes", "failures"]

    def __init__(self, memory_dir: Path, session_id: Optional[str] = None, user_id: Optional[str] = None):
        """
        Initialize the memory system.

        Args:
            memory_dir: Root directory for memory storage
            session_id: Current session identifier
            user_id: Current user identifier
        """
        print(f"[MEMORY INIT] Initializing MemorySystem...")
        self.memory_dir = Path(memory_dir)
        self.session_id = session_id or f"session_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        self.user_id = user_id or "default_user"

        print(f"[MEMORY INIT] Memory directory: {self.memory_dir}")
        print(f"[MEMORY INIT] Categories: {self.CATEGORIES}")

        # Create directory structure
        self._initialize_directories()
        print(f"[MEMORY INIT] ✓ Directory structure created")

    def _initialize_directories(self) -> None:
        """Create the memory directory structure."""
        self.memory_dir.mkdir(parents=True, exist_ok=True)

        for category in self.CATEGORIES:
            (self.memory_dir / category).mkdir(exist_ok=True)

        # Outcomes has subcategories
        for subcat in self.OUTCOME_SUBCATEGORIES:
            (self.memory_dir / "outcomes" / subcat).mkdir(parents=True, exist_ok=True)

    def _generate_filename(self, topic: str) -> str:
        """Generate a filename from timestamp and topic."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        # Sanitize topic for filename
        safe_topic = re.sub(r'[^\w\s-]', '', topic)[:50].strip().replace(' ', '_').lower()
        return f"{timestamp}_{safe_topic}.md"

    def _get_category_path(self, category: str, subcategory: Optional[str] = None) -> Path:
        """Get the path for a category, creating subdirectories as needed."""
        if category == "users":
            path = self.memory_dir / "users" / self.user_id
        elif category == "episodes":
            path = self.memory_dir / "episodes" / self.session_id
        elif category == "skills":
            subcat = subcategory or "general"
            path = self.memory_dir / "skills" / subcat
        elif category == "outcomes":
            subcat = subcategory or "general"
            if subcat not in ["successes", "failures", "general"]:
                subcat = "general"
            path = self.memory_dir / "outcomes" / subcat
        else:
            path = self.memory_dir / category

        path.mkdir(parents=True, exist_ok=True)
        return path

    def add_to_memory(
        self,
        category: str,
        content: str,
        topic: str = "note",
        subcategory: Optional[str] = None,
        tags: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """
        Add a new memory entry.

        Args:
            category: One of "users", "episodes", "skills", "outcomes"
            content: The information to store (markdown formatted)
            topic: Brief topic for the filename
            subcategory: For skills (skill name) or outcomes (successes/failures)
            tags: Optional tags for the entry

        Returns:
            Dict with status and file path
        """
        print(f"\n[MEMORY WRITE] Adding memory entry...")
        print(f"[MEMORY WRITE] Category: {category}, Topic: {topic}")
        if subcategory:
            print(f"[MEMORY WRITE] Subcategory: {subcategory}")
        if tags:
            print(f"[MEMORY WRITE] Tags: {tags}")

        if category not in self.CATEGORIES:
            print(f"[MEMORY WRITE] ✗ Error: Invalid category '{category}'")
            return {
                "success": False,
                "error": f"Invalid category '{category}'. Must be one of: {self.CATEGORIES}",
            }

        # Create the entry
        entry = MemoryEntry(
            content=content,
            category=category,
            subcategory=subcategory,
            tags=tags or [],
            session_id=self.session_id,
            user_id=self.user_id,
        )

        # Determine file path
        category_path = self._get_category_path(category, subcategory)
        filename = self._generate_filename(topic)
        file_path = category_path / filename

        # Write to file
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(entry.to_markdown())

            print(f"[MEMORY WRITE] ✓ Saved to: {file_path.relative_to(self.memory_dir)}")
            print(f"[MEMORY WRITE] Content preview: {content[:100]}{'...' if len(content) > 100 else ''}")

            return {
                "success": True,
                "file_path": str(file_path.relative_to(self.memory_dir)),
                "category": category,
                "timestamp": entry.timestamp,
            }
        except Exception as e:
            print(f"[MEMORY WRITE] ✗ Error: {e}")
            return {
                "success": False,
                "error": str(e),
            }

    def search_memory(
        self,
        query: str,
        category: Optional[str] = None,
        max_results: int = 10,
    ) -> Dict[str, Any]:
        """
        Search memory for relevant entries (grep-like).

        Args:
            query: Search term (case-insensitive substring match)
            category: Optional category to limit search
            max_results: Maximum number of results to return

        Returns:
            Dict with matching entries
        """
        print(f"\n[MEMORY SEARCH] Searching for: '{query}'")
        if category:
            print(f"[MEMORY SEARCH] Category filter: {category}")

        results = []
        search_dirs = []

        # Determine which directories to search
        if category:
            if category not in self.CATEGORIES:
                print(f"[MEMORY SEARCH] ✗ Error: Invalid category '{category}'")
                return {
                    "success": False,
                    "error": f"Invalid category '{category}'",
                    "results": [],
                }
            search_dirs = [self.memory_dir / category]
        else:
            search_dirs = [self.memory_dir / cat for cat in self.CATEGORIES]

        print(f"[MEMORY SEARCH] Searching in: {[str(d.name) for d in search_dirs]}")

        # Search through files
        query_lower = query.lower()
        files_searched = 0

        for search_dir in search_dirs:
            if not search_dir.exists():
                continue

            for file_path in search_dir.rglob("*.md"):
                files_searched += 1
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        content = f.read()

                    if query_lower in content.lower():
                        entry = MemoryEntry.from_markdown(content, str(file_path))

                        # Find matching lines for context
                        matching_lines = []
                        for i, line in enumerate(content.split('\n')):
                            if query_lower in line.lower():
                                matching_lines.append({
                                    "line_number": i + 1,
                                    "text": line.strip()[:200],
                                })

                        results.append({
                            "file": str(file_path.relative_to(self.memory_dir)),
                            "category": entry.category,
                            "timestamp": entry.timestamp,
                            "tags": entry.tags,
                            "matching_lines": matching_lines[:5],  # Limit context
                            "content_preview": entry.content[:300] + "..." if len(entry.content) > 300 else entry.content,
                        })

                        if len(results) >= max_results:
                            break
                except Exception:
                    continue

            if len(results) >= max_results:
                break

        # Sort by timestamp (most recent first)
        results.sort(key=lambda x: x.get("timestamp", ""), reverse=True)

        print(f"[MEMORY SEARCH] ✓ Searched {files_searched} files, found {len(results)} matches")
        if results:
            print(f"[MEMORY SEARCH] Top result: {results[0]['file']}")

        return {
            "success": True,
            "query": query,
            "category_filter": category,
            "total_results": len(results),
            "results": results[:max_results],
        }

    def get_user_profile(self, user_id: Optional[str] = None) -> Dict[str, Any]:
        """
        Get all memory entries for a specific user.

        Args:
            user_id: User to get profile for (defaults to current user)

        Returns:
            Dict with user's memory entries
        """
        target_user = user_id or self.user_id
        user_dir = self.memory_dir / "users" / target_user

        if not user_dir.exists():
            return {
                "success": True,
                "user_id": target_user,
                "entries": [],
                "message": "No memory entries found for this user",
            }

        entries = []
        for file_path in sorted(user_dir.glob("*.md"), reverse=True):
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.read()
                entry = MemoryEntry.from_markdown(content, str(file_path))
                entries.append({
                    "file": file_path.name,
                    "timestamp": entry.timestamp,
                    "tags": entry.tags,
                    "content": entry.content,
                })
            except Exception:
                continue

        return {
            "success": True,
            "user_id": target_user,
            "total_entries": len(entries),
            "entries": entries,
        }

    def get_session_history(self, session_id: Optional[str] = None) -> Dict[str, Any]:
        """
        Get all memory entries for a specific session.

        Args:
            session_id: Session to get history for (defaults to current session)

        Returns:
            Dict with session's memory entries
        """
        target_session = session_id or self.session_id
        session_dir = self.memory_dir / "episodes" / target_session

        if not session_dir.exists():
            return {
                "success": True,
                "session_id": target_session,
                "entries": [],
                "message": "No memory entries found for this session",
            }

        entries = []
        for file_path in sorted(session_dir.glob("*.md")):
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.read()
                entry = MemoryEntry.from_markdown(content, str(file_path))
                entries.append({
                    "file": file_path.name,
                    "timestamp": entry.timestamp,
                    "tags": entry.tags,
                    "content": entry.content,
                })
            except Exception:
                continue

        return {
            "success": True,
            "session_id": target_session,
            "total_entries": len(entries),
            "entries": entries,
        }

    def list_categories(self) -> Dict[str, Any]:
        """List all categories and their contents summary."""
        summary = {}

        for category in self.CATEGORIES:
            category_path = self.memory_dir / category
            if category_path.exists():
                file_count = len(list(category_path.rglob("*.md")))
                subdirs = [d.name for d in category_path.iterdir() if d.is_dir()]
                summary[category] = {
                    "file_count": file_count,
                    "subdirectories": subdirs,
                }
            else:
                summary[category] = {
                    "file_count": 0,
                    "subdirectories": [],
                }

        return {
            "success": True,
            "memory_dir": str(self.memory_dir),
            "categories": summary,
        }


# =============================================================================
# Tool Functions (for agent integration)
# =============================================================================

def create_memory_tools(memory_system: MemorySystem) -> Dict[str, callable]:
    """
    Create tool functions that wrap the memory system.

    These are the functions that get exposed to the agent.
    """

    def add_to_memory(
        category: str,
        content: str,
        topic: str = "note",
        subcategory: str = "",
        tags: str = "",
    ) -> str:
        """
        Store information in memory.

        Args:
            category: One of "users", "episodes", "skills", "outcomes"
            content: Information to remember (markdown formatted)
            topic: Brief topic (used in filename)
            subcategory: For skills (skill name) or outcomes ("successes"/"failures")
            tags: Comma-separated tags

        Returns:
            Confirmation of what was stored
        """
        tag_list = [t.strip() for t in tags.split(",") if t.strip()] if tags else []

        result = memory_system.add_to_memory(
            category=category,
            content=content,
            topic=topic,
            subcategory=subcategory if subcategory else None,
            tags=tag_list,
        )

        if result["success"]:
            return f"Memory stored: {result['file_path']} (category: {category}, timestamp: {result['timestamp']})"
        else:
            return f"Failed to store memory: {result['error']}"

    def search_memory(query: str, category: str = "") -> str:
        """
        Search memory for relevant information.

        Args:
            query: Search term (case-insensitive)
            category: Optional category to search ("users", "episodes", "skills", "outcomes")

        Returns:
            Matching memory entries
        """
        result = memory_system.search_memory(
            query=query,
            category=category if category else None,
        )

        if not result["success"]:
            return f"Search failed: {result.get('error', 'Unknown error')}"

        if not result["results"]:
            return f"No results found for query: '{query}'"

        # Format results
        lines = [f"Found {result['total_results']} result(s) for '{query}':", ""]

        for i, r in enumerate(result["results"], 1):
            lines.append(f"### Result {i}: {r['file']}")
            lines.append(f"**Category**: {r['category']} | **Time**: {r['timestamp']}")
            if r["tags"]:
                lines.append(f"**Tags**: {', '.join(r['tags'])}")
            lines.append(f"\n{r['content_preview']}\n")

        return "\n".join(lines)

    return {
        "add_to_memory": add_to_memory,
        "search_memory": search_memory,
    }
