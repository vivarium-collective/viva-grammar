"""The ABM-agnostic grammar engine: signals -> behaviors via Hill math.

Clean-room reimplementation of the cell-behavior hypothesis grammar
(Johnson et al., *Cell* 2025).  See ``references/PROVENANCE.md``.
"""
from __future__ import annotations

from .behaviors import Behavior, BehaviorRegistry, ResolvedBehavior, default_registry, parse_behavior
from .cell_definitions import CellDefinitions
from .evaluator import BehaviorValue, GrammarEvaluator
from .response import combine, hill_response, hill_term, linear_response, saturating_response, step_response
from .rules import DECREASES, INCREASES, Rule, RuleSet
from .signals import SignalContext, resolve_signal

__all__ = [
    "Behavior", "BehaviorRegistry", "ResolvedBehavior", "default_registry", "parse_behavior",
    "CellDefinitions",
    "BehaviorValue", "GrammarEvaluator",
    "combine", "hill_response", "hill_term", "linear_response", "saturating_response", "step_response",
    "Rule", "RuleSet", "INCREASES", "DECREASES",
    "SignalContext", "resolve_signal",
]
