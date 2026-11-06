"""The grammar's dictionary of **behaviors** (the verbs: what a cell does).

A *behavior* is a cell-scale phenotype parameter the grammar can modulate:
a rate (cycling, death, transformation, secretion, attack) or a continuous
parameter (migration speed/bias, adhesion, chemotactic sensitivity).  Each
behavior carries metadata the evaluator and the host adapters need:

- ``kind``   — ``"rate"`` (per-time, becomes a stochastic event via
  ``Pr = 1 - exp(-rate * interval)``) or ``"parameter"`` (a continuous value
  the host reads directly, e.g. migration speed).
- ``default_base`` / ``default_min`` — fallbacks for ``b0`` / ``bm`` when a
  :class:`~viva_grammar.engine.cell_definitions.CellDefinitions` entry omits
  them (the rules CSV carries only the *saturated* value ``bM``).
- ``parametric`` — ``True`` for family behaviors like ``"transform to <type>"``
  whose target is parsed out of the behavior string.

This is a clean-room reimplementation of a subset of PhysiCell's behavior
dictionary (Methods S1), covering what the hypoxia example and the blueprint
ports need, with a registry anyone can extend.  See ``references/PROVENANCE.md``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

__all__ = ["Behavior", "BehaviorRegistry", "default_registry", "parse_behavior"]

RATE = "rate"
PARAMETER = "parameter"


@dataclass(frozen=True)
class Behavior:
    """Metadata for one behavior name."""

    name: str
    kind: str  # RATE or PARAMETER
    default_base: float = 0.0
    default_min: float = 0.0
    units: str = ""
    parametric: bool = False  # True for "transform to <type>" style families
    description: str = ""

    @property
    def is_rate(self) -> bool:
        return self.kind == RATE


# Family behaviors: the behavior string is "<prefix> <target>" and the target
# (a cell type / substrate name) is parsed off.  Keyed by prefix.
_FAMILIES = {
    "transform to": Behavior(
        "transform to", RATE, 0.0, 0.0, "1/min", parametric=True,
        description="Rate of transition into the named target cell type.",
    ),
    "transition to": Behavior(
        "transition to", RATE, 0.0, 0.0, "1/min", parametric=True,
        description="Alias of 'transform to'.",
    ),
}


def _base_registry() -> dict:
    items = [
        Behavior("cycle entry", RATE, 0.0, 0.0, "1/min",
                 description="Rate of entering the cell cycle (division propensity)."),
        Behavior("apoptosis", RATE, 0.0, 0.0, "1/min",
                 description="Programmed (apoptotic) death rate."),
        Behavior("necrosis", RATE, 0.0, 0.0, "1/min",
                 description="Necrotic death rate."),
        Behavior("death", RATE, 0.0, 0.0, "1/min",
                 description="Generic death rate (maps to a death event)."),
        Behavior("migration speed", PARAMETER, 0.0, 0.0, "micron/min",
                 description="Random-motility speed."),
        Behavior("migration bias", PARAMETER, 0.0, 0.0, "dimensionless",
                 description="Directional bias of migration in [0,1]."),
        Behavior("migration persistence time", PARAMETER, 0.0, 0.0, "min",
                 description="Persistence time of the biased random walk."),
        Behavior("chemotactic response", PARAMETER, 0.0, 0.0, "dimensionless",
                 description="Chemotactic sensitivity toward a substrate gradient."),
        Behavior("cell-cell adhesion", PARAMETER, 0.0, 0.0, "dimensionless",
                 description="Homotypic adhesion strength."),
        Behavior("cell-cell repulsion", PARAMETER, 0.0, 0.0, "dimensionless",
                 description="Repulsion strength."),
        Behavior("secretion", RATE, 0.0, 0.0, "1/min", parametric=True,
                 description="Secretion rate of the named substrate."),
        Behavior("uptake", RATE, 0.0, 0.0, "1/min", parametric=True,
                 description="Uptake rate of the named substrate."),
    ]
    return {b.name: b for b in items}


class BehaviorRegistry:
    """Resolvable set of behaviors; extend with :meth:`register`."""

    def __init__(self, behaviors: Optional[dict] = None, families: Optional[dict] = None):
        self._behaviors = dict(behaviors if behaviors is not None else _base_registry())
        self._families = dict(families if families is not None else _FAMILIES)

    def register(self, behavior: Behavior) -> None:
        if behavior.parametric:
            self._families[behavior.name] = behavior
        else:
            self._behaviors[behavior.name] = behavior

    def resolve(self, behavior_name: str) -> "ResolvedBehavior":
        """Resolve a raw behavior string to its :class:`Behavior` + target.

        Raises ``KeyError`` with a helpful message for an unknown behavior.
        """
        name = behavior_name.strip()
        if name in self._behaviors:
            return ResolvedBehavior(self._behaviors[name], target=None)
        # Try family prefixes (longest prefix first for determinism).
        for prefix in sorted(self._families, key=len, reverse=True):
            if name == prefix:
                return ResolvedBehavior(self._families[prefix], target=None)
            if name.startswith(prefix + " "):
                target = name[len(prefix) + 1:].strip()
                return ResolvedBehavior(self._families[prefix], target=target)
        raise KeyError(
            f"unknown behavior {behavior_name!r}; known behaviors: "
            f"{sorted(self._behaviors)} and families {sorted(self._families)}"
        )

    def __contains__(self, name: str) -> bool:
        try:
            self.resolve(name)
            return True
        except KeyError:
            return False


@dataclass(frozen=True)
class ResolvedBehavior:
    """A behavior resolved from a rule string, plus any parsed target."""

    behavior: Behavior
    target: Optional[str]

    @property
    def is_rate(self) -> bool:
        return self.behavior.is_rate

    @property
    def canonical(self) -> str:
        """Canonical key used to group rules and index outputs."""
        if self.target is not None:
            return f"{self.behavior.name} {self.target}"
        return self.behavior.name


def default_registry() -> BehaviorRegistry:
    return BehaviorRegistry()


def parse_behavior(behavior_name: str, registry: Optional[BehaviorRegistry] = None) -> ResolvedBehavior:
    return (registry or default_registry()).resolve(behavior_name)
