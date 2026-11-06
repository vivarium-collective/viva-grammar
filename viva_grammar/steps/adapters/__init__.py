"""Host adapters isolating ABM coupling (blueprint port surface, reference helpers)."""
from __future__ import annotations

from .blueprint import BlueprintGrammarProcess
from .reference import context_from_dict

__all__ = ["BlueprintGrammarProcess", "context_from_dict"]
