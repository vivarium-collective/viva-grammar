"""CellDefinitions: the base/min behavior values the rules CSV does not carry.

The 8-column PhysiCell rules CSV stores only each rule's *saturated* value
(``bM`` for an "increases" rule, ``bm`` for a "decreases" rule).  The base
value ``b0`` (behavior in the absence of signal) and, where a behavior is only
up-regulated, its min ``bm`` live in the cell **definition** (PhysiCell's cell
XML), not the rules table — so this module supplies them, mirroring that split.

A :class:`CellDefinitions` maps ``cell_type -> {canonical_behavior -> {base, min}}``.
Lookups fall back to the behavior registry's ``default_base`` / ``default_min``
when a cell_def omits a behavior.  For the blueprint adapter (which has no cell
XML) the same config is provided inline.

See ``references/PROVENANCE.md`` for the hypoxia baseline table.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

from .behaviors import BehaviorRegistry, default_registry

__all__ = ["CellDefinitions"]


@dataclass
class _BehaviorBase:
    base: float = 0.0
    min: float = 0.0


class CellDefinitions:
    """Per-(cell_type, behavior) base/min values with registry fallbacks."""

    def __init__(self, definitions: Optional[dict] = None, registry: Optional[BehaviorRegistry] = None):
        self.registry = registry or default_registry()
        # cell_type -> canonical_behavior -> _BehaviorBase
        self._defs: Dict[str, Dict[str, _BehaviorBase]] = {}
        for ct, behaviors in (definitions or {}).items():
            self._defs[ct] = {}
            for behavior, vals in (behaviors or {}).items():
                canonical = self.registry.resolve(behavior).canonical
                if isinstance(vals, (int, float)):
                    vals = {"base": float(vals)}
                self._defs[ct][canonical] = _BehaviorBase(
                    base=float(vals.get("base", 0.0)),
                    min=float(vals.get("min", 0.0)),
                )

    def _fallback(self, canonical_behavior: str) -> _BehaviorBase:
        # Resolve the behavior to pull its registry defaults.
        try:
            b = self.registry.resolve(canonical_behavior).behavior
            return _BehaviorBase(base=b.default_base, min=b.default_min)
        except KeyError:
            return _BehaviorBase()

    def base(self, cell_type: str, canonical_behavior: str) -> float:
        entry = self._defs.get(cell_type, {}).get(canonical_behavior)
        if entry is not None:
            return entry.base
        return self._fallback(canonical_behavior).base

    def min(self, cell_type: str, canonical_behavior: str) -> float:
        entry = self._defs.get(cell_type, {}).get(canonical_behavior)
        if entry is not None:
            return entry.min
        return self._fallback(canonical_behavior).min

    def to_dict(self) -> dict:
        return {
            ct: {b: {"base": v.base, "min": v.min} for b, v in behaviors.items()}
            for ct, behaviors in self._defs.items()
        }
