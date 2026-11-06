"""viva-grammar: a generic, ABM-agnostic implementation of the cell-behavior
hypothesis grammar (Johnson et al., *Cell* 2025) as a process-bigraph module.

Clean-room reimplementation of the grammar's signals->behaviors Hill math;
not PhysiCell.  See ``references/PROVENANCE.md``.
"""
from __future__ import annotations

from .core import build_core
from .engine import (
    CellDefinitions, GrammarEvaluator, Rule, RuleSet, SignalContext, default_registry,
)
from .steps import BlueprintGrammarProcess, GrammarRulesProcess
from .reference import ReferenceSubstrateProcess

__all__ = [
    "build_core",
    "GrammarRulesProcess", "BlueprintGrammarProcess", "ReferenceSubstrateProcess",
    "CellDefinitions", "GrammarEvaluator", "Rule", "RuleSet", "SignalContext", "default_registry",
]
