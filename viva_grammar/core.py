"""build_core(): register viva-grammar's processes on a process-bigraph core.

Called by every ``sims/run.py`` and by the dashboard so the workspace's own
processes resolve regardless of how the editable install was made (see the
viva-expert Auto-Discovery convention).
"""
from __future__ import annotations

from process_bigraph import allocate_core

from .reference.substrate import ReferenceSubstrateProcess
from .steps.adapters.blueprint import BlueprintGrammarProcess
from .steps.grammar_process import GrammarRulesProcess

__all__ = ["build_core"]


def build_core(core=None):
    if core is None:
        core = allocate_core()
    core.register_link("GrammarRulesProcess", GrammarRulesProcess)
    core.register_link("BlueprintGrammarProcess", BlueprintGrammarProcess)
    core.register_link("ReferenceSubstrateProcess", ReferenceSubstrateProcess)
    return core
