"""
Predefined scenarios for collaborative agent experiments.

Each scenario defines:
- Task description
- Data files
- Expected user goals (for simulated experiments)
- Evaluation criteria
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class Scenario:
    """A predefined experimental scenario."""
    name: str
    description: str
    task_description: str
    data_files: Dict[str, str]
    expected_analyses: List[str]  # What analyses the user likely wants
    success_criteria: Dict[str, Any]  # How to evaluate success
    difficulty: str = "medium"  # easy, medium, hard
    tags: List[str] = field(default_factory=list)


# Predefined scenarios
SCENARIOS = {
    # ==========================================================================
    # SIMPLE TEST SCENARIOS (for memory system development)
    # ==========================================================================
    "simple_hello": Scenario(
        name="simple_hello",
        description="Simple greeting test",
        task_description="A simple conversation to test basic agent functionality.",
        data_files={},
        expected_analyses=[],
        success_criteria={"responded": True},
        difficulty="easy",
        tags=["simple", "test"],
    ),

    "memory_test": Scenario(
        name="memory_test",
        description="Test memory system with simple questions",
        task_description="Test the agent's ability to remember and recall information.",
        data_files={},
        expected_analyses=[],
        success_criteria={"used_memory": True},
        difficulty="easy",
        tags=["memory", "test"],
    ),

    # ==========================================================================
    # RESEARCH SCENARIOS
    # ==========================================================================
    "scenario_2": Scenario(
        name="scenario_2",
        description="Basic differential expression analysis with treatment vs control",
        task_description="Analyze transcriptomics data to identify differentially expressed genes between treatment and control conditions. The goal is to understand what genes change in response to treatment.",
        data_files={
            "assay_file": "scenarios/scenario_2/scenario_2_assay_data.csv",
            "sample_file": "scenarios/scenario_2/scenario_2_sample_data.csv",
        },
        expected_analyses=[
            "differential_expression",
            "gene_statistics",
            "pca",
        ],
        success_criteria={
            "identified_de_genes": True,
            "provided_interpretation": True,
            "answered_user_questions": True,
        },
        difficulty="easy",
        tags=["differential_expression", "basic"],
    ),

    "scenario_11a": Scenario(
        name="scenario_11a",
        description="Xenopus anesthetic response analysis",
        task_description="Analyze transcriptomics data from Xenopus laevis experiments studying anesthetic drug responses. The goal is to understand molecular mechanisms underlying anesthetic effects.",
        data_files={
            "assay_file": "scenarios/scenario_11a/scenario_11a_assay_data.csv",
            "sample_file": "scenarios/scenario_11a/scenario_11a_sample_data.csv",
        },
        expected_analyses=[
            "differential_expression",
            "pca",
            "batch_effect",
            "gene_statistics",
        ],
        success_criteria={
            "identified_de_genes": True,
            "checked_batch_effects": True,
            "provided_biological_interpretation": True,
        },
        difficulty="medium",
        tags=["differential_expression", "xenopus", "pharmacology"],
    ),

    "exploratory_analysis": Scenario(
        name="exploratory_analysis",
        description="Open-ended data exploration",
        task_description="Explore this transcriptomics dataset to find interesting patterns and generate hypotheses. No specific question - just see what's there.",
        data_files={
            "assay_file": "scenarios/scenario_2/scenario_2_assay_data.csv",
            "sample_file": "scenarios/scenario_2/scenario_2_sample_data.csv",
        },
        expected_analyses=[
            "describe_data",
            "gene_statistics",
            "pca",
            "batch_effect",
        ],
        success_criteria={
            "explored_data_structure": True,
            "identified_patterns": True,
            "generated_hypotheses": True,
        },
        difficulty="hard",
        tags=["exploratory", "open-ended"],
    ),
}


def get_scenario(name: str) -> Optional[Scenario]:
    """Get a scenario by name."""
    return SCENARIOS.get(name)


def list_scenarios() -> List[str]:
    """List all available scenario names."""
    return list(SCENARIOS.keys())


def get_scenario_by_tags(tags: List[str]) -> List[Scenario]:
    """Get scenarios that match any of the given tags."""
    results = []
    for scenario in SCENARIOS.values():
        if any(tag in scenario.tags for tag in tags):
            results.append(scenario)
    return results


# Experimental conditions for systematic studies
EXPERIMENTAL_CONDITIONS = {
    "user_expertise": [
        ("expert_biologist", "Expert user who knows what they want"),
        ("wet_lab_scientist", "Domain expert but limited computational skills"),
        ("graduate_student", "Learning, needs guidance"),
        ("pi_busy", "Wants quick results, impatient"),
    ],

    "user_beliefs": [
        ("has_hypothesis", "User has specific hypothesis to test"),
        ("skeptical_of_ai", "User doesn't fully trust AI"),
        ("time_pressured", "User is under deadline pressure"),
        ("exploratory", "User wants open-ended exploration"),
    ],

    "task_complexity": [
        ("simple", "Single analysis, clear goal"),
        ("moderate", "Multiple analyses, some ambiguity"),
        ("complex", "Open-ended, requires iteration"),
    ],
}


def generate_experiment_matrix() -> List[Dict[str, Any]]:
    """
    Generate a matrix of experimental conditions.

    Returns a list of experiment configurations for systematic study.
    """
    experiments = []

    for scenario_name, scenario in SCENARIOS.items():
        for profile_name, profile_desc in EXPERIMENTAL_CONDITIONS["user_expertise"]:
            for beliefs_name, beliefs_desc in EXPERIMENTAL_CONDITIONS["user_beliefs"]:
                experiments.append({
                    "scenario": scenario_name,
                    "human_profile": profile_name,
                    "human_beliefs": beliefs_name,
                    "description": f"{scenario.description} with {profile_desc} who {beliefs_desc.lower()}",
                    "expected_difficulty": scenario.difficulty,
                })

    return experiments
