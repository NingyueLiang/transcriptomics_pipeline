#!/usr/bin/env python3
"""
Beautiful console output utilities for the transcriptomics pipeline.

Provides colorful, structured terminal output with ANSI escape codes.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Dict


# ═══════════════════════════════════════════════════════════════════════════════
# ANSI COLOR CODES
# ═══════════════════════════════════════════════════════════════════════════════

class Colors:
    """ANSI color codes for beautiful terminal output."""
    HEADER = '\033[95m'      # Magenta
    BLUE = '\033[94m'        # Blue
    CYAN = '\033[96m'        # Cyan
    GREEN = '\033[92m'       # Green
    YELLOW = '\033[93m'      # Yellow
    RED = '\033[91m'         # Red
    BOLD = '\033[1m'         # Bold
    DIM = '\033[2m'          # Dim
    ITALIC = '\033[3m'       # Italic
    UNDERLINE = '\033[4m'    # Underline
    RESET = '\033[0m'        # Reset


# ═══════════════════════════════════════════════════════════════════════════════
# BANNER & SECTION HEADERS
# ═══════════════════════════════════════════════════════════════════════════════

def print_banner():
    """Print a beautiful startup banner."""
    banner = f"""
{Colors.CYAN}{Colors.BOLD}
  ╔══════════════════════════════════════════════════════════════════════════╗
  ║                                                                          ║
  ║   ████████╗██████╗  █████╗ ███╗   ██╗███████╗ ██████╗██████╗ ██╗██████╗  ║
  ║   ╚══██╔══╝██╔══██╗██╔══██╗████╗  ██║██╔════╝██╔════╝██╔══██╗██║██╔══██╗ ║
  ║      ██║   ██████╔╝███████║██╔██╗ ██║███████╗██║     ██████╔╝██║██████╔╝ ║
  ║      ██║   ██╔══██╗██╔══██║██║╚██╗██║╚════██║██║     ██╔══██╗██║██╔═══╝  ║
  ║      ██║   ██║  ██║██║  ██║██║ ╚████║███████║╚██████╗██║  ██║██║██║      ║
  ║      ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝╚══════╝ ╚═════╝╚═╝  ╚═╝╚═╝╚═╝      ║
  ║                                                                          ║
  ╚══════════════════════════════════════════════════════════════════════════╝
{Colors.RESET}"""
    print(banner)


def print_section(title: str, icon: str = "▶"):
    """Print a section header."""
    width = 70
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'─' * width}{Colors.RESET}")
    print(f"{Colors.BOLD}{Colors.CYAN}  {icon}  {title}{Colors.RESET}")
    print(f"{Colors.BOLD}{Colors.BLUE}{'─' * width}{Colors.RESET}\n")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP & SUBSTEP INDICATORS
# ═══════════════════════════════════════════════════════════════════════════════

def print_step(step_num: int, total: int, description: str):
    """Print a step indicator."""
    progress = f"[{step_num}/{total}]"
    print(f"{Colors.BOLD}{Colors.GREEN}  ● {progress}{Colors.RESET} {description}")


def print_substep(description: str, status: str = "running"):
    """Print a substep with status.
    
    Args:
        description: Text to display
        status: One of "running", "done", "skip", or other (default arrow)
    """
    if status == "running":
        symbol = f"{Colors.YELLOW}⟳{Colors.RESET}"
    elif status == "done":
        symbol = f"{Colors.GREEN}✓{Colors.RESET}"
    elif status == "skip":
        symbol = f"{Colors.DIM}○{Colors.RESET}"
    else:
        symbol = f"{Colors.BLUE}→{Colors.RESET}"
    print(f"      {symbol} {Colors.DIM}{description}{Colors.RESET}")


# ═══════════════════════════════════════════════════════════════════════════════
# INFO & STATUS MESSAGES
# ═══════════════════════════════════════════════════════════════════════════════

def print_info(label: str, value: str):
    """Print an info line with label and value."""
    print(f"      {Colors.DIM}├─{Colors.RESET} {Colors.BOLD}{label}:{Colors.RESET} {Colors.CYAN}{value}{Colors.RESET}")


def print_success(message: str):
    """Print a success message."""
    print(f"\n  {Colors.GREEN}✅ {Colors.BOLD}{message}{Colors.RESET}")


def print_warning(message: str):
    """Print a warning message."""
    print(f"\n  {Colors.YELLOW}⚠️  {Colors.BOLD}{message}{Colors.RESET}")


def print_error(message: str):
    """Print an error message."""
    print(f"\n  {Colors.RED}❌ {Colors.BOLD}{message}{Colors.RESET}")


# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION & TIMING
# ═══════════════════════════════════════════════════════════════════════════════

def print_config(config: Dict[str, Any]):
    """Print configuration in a nice format."""
    print(f"\n  {Colors.BOLD}{Colors.HEADER}📋 Configuration:{Colors.RESET}")
    for key, value in config.items():
        if isinstance(value, Path):
            value = str(value)
        print(f"      {Colors.DIM}├─{Colors.RESET} {key}: {Colors.CYAN}{value}{Colors.RESET}")


def print_timing(start_time: float, step_name: str):
    """Print elapsed time for a step."""
    elapsed = time.time() - start_time
    if elapsed < 60:
        time_str = f"{elapsed:.1f}s"
    else:
        mins = int(elapsed // 60)
        secs = elapsed % 60
        time_str = f"{mins}m {secs:.1f}s"
    print(f"      {Colors.DIM}└─ Completed in {Colors.YELLOW}{time_str}{Colors.RESET}")


def format_duration(seconds: float) -> str:
    """Format a duration in seconds to a human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f} seconds"
    else:
        mins = int(seconds // 60)
        secs = seconds % 60
        return f"{mins} min {secs:.1f} sec"


# ═══════════════════════════════════════════════════════════════════════════════
# COMPLETION BANNERS
# ═══════════════════════════════════════════════════════════════════════════════

def print_completion_banner():
    """Print a success completion banner."""
    print(f"""
{Colors.GREEN}{Colors.BOLD}
  ╔══════════════════════════════════════════════════════════════════════════╗
  ║                                                                          ║
  ║   ✅  PIPELINE COMPLETED SUCCESSFULLY                                    ║
  ║                                                                          ║
  ╚══════════════════════════════════════════════════════════════════════════╝
{Colors.RESET}""")


def print_results_summary(
    experiment_name: str,
    summary_path: str,
    report_path: str,
    total_time: float
):
    """Print the final results summary."""
    time_str = format_duration(total_time)
    
    print(f"  {Colors.BOLD}📊 Results Summary:{Colors.RESET}")
    print(f"      {Colors.DIM}├─{Colors.RESET} Experiment: {Colors.CYAN}{experiment_name}{Colors.RESET}")
    print(f"      {Colors.DIM}├─{Colors.RESET} Summary: {Colors.CYAN}{summary_path}{Colors.RESET}")
    print(f"      {Colors.DIM}├─{Colors.RESET} Report: {Colors.CYAN}{report_path}{Colors.RESET}")
    print(f"      {Colors.DIM}└─{Colors.RESET} Total Time: {Colors.YELLOW}{time_str}{Colors.RESET}")
    print()

