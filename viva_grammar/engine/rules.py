"""Rule and RuleSet: the parsed form of a cell-behavior hypothesis.

A *rule* is one grammar statement:

    In [cell_type], [signal] increases|decreases [behavior]
    (with a saturated value, half-max, Hill power, applies-to-dead flag)

i.e. one row of the PhysiCell 8-column rules CSV.  A :class:`RuleSet` groups
rules by ``(cell_type, behavior)`` so the evaluator can assemble the
multivariate up/down responses of eq 4.

Clean-room reimplementation; see ``references/PROVENANCE.md``.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from .behaviors import BehaviorRegistry, ResolvedBehavior, default_registry

__all__ = ["Rule", "RuleSet", "INCREASES", "DECREASES", "HILL", "LINEAR", "STEP"]

INCREASES = "increases"
DECREASES = "decreases"

HILL = "hill"
LINEAR = "linear"
STEP = "step"


@dataclass(frozen=True)
class Rule:
    """One parsed grammar rule (one CSV row)."""

    cell_type: str
    signal: str
    direction: str          # INCREASES or DECREASES
    behavior: str           # raw behavior string (may be a family like "transform to X")
    saturated_value: float  # bM for increases, bm for decreases
    half_max: float
    hill_power: float
    applies_to_dead: bool = False
    response_type: str = HILL

    def __post_init__(self):
        if self.direction not in (INCREASES, DECREASES):
            raise ValueError(
                f"rule direction must be {INCREASES!r} or {DECREASES!r}, got {self.direction!r}"
            )
        if self.response_type not in (HILL, LINEAR, STEP):
            raise ValueError(f"unknown response_type {self.response_type!r}")

    @property
    def increases(self) -> bool:
        return self.direction == INCREASES


class RuleSet:
    """A validated collection of rules, grouped for evaluation.

    Validation resolves every behavior name against the supplied
    :class:`BehaviorRegistry` (unknown behaviors raise early with a helpful
    message), and records each behavior's resolved metadata + target.
    """

    def __init__(self, rules: List[Rule], registry: Optional[BehaviorRegistry] = None):
        self.registry = registry or default_registry()
        self.rules: List[Rule] = list(rules)
        # (cell_type, canonical_behavior) -> {"resolved": ResolvedBehavior, "rules": [Rule...]}
        self._groups: Dict[Tuple[str, str], dict] = {}
        self._resolved_cache: Dict[str, ResolvedBehavior] = {}
        for rule in self.rules:
            resolved = self._resolve(rule.behavior)
            key = (rule.cell_type, resolved.canonical)
            grp = self._groups.setdefault(key, {"resolved": resolved, "rules": []})
            grp["rules"].append(rule)

    def _resolve(self, behavior: str) -> ResolvedBehavior:
        if behavior not in self._resolved_cache:
            self._resolved_cache[behavior] = self.registry.resolve(behavior)
        return self._resolved_cache[behavior]

    @property
    def cell_types(self) -> List[str]:
        seen = []
        for (ct, _), _ in self._groups.items():
            if ct not in seen:
                seen.append(ct)
        return seen

    def groups_for(self, cell_type: str) -> List[dict]:
        """All ``{resolved, rules}`` behavior-groups that apply to ``cell_type``."""
        return [
            {"canonical": canon, **grp}
            for (ct, canon), grp in self._groups.items()
            if ct == cell_type
        ]

    def __len__(self) -> int:
        return len(self.rules)

    def __iter__(self):
        return iter(self.rules)
