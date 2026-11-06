"""GrammarEvaluator: signals -> behavior parameters, per cell (eq 4).

Given a :class:`~viva_grammar.engine.rules.RuleSet` and
:class:`~viva_grammar.engine.cell_definitions.CellDefinitions`, evaluate the
behavior values for one cell's :class:`~viva_grammar.engine.signals.SignalContext`.

For each ``(cell_type, behavior)`` group the evaluator:
  1. splits the group's rules into up-regulating ("increases") and
     down-regulating ("decreases");
  2. forms the multivariate up/down responses ``U`` and ``D`` from each rule's
     signal value (eqs 2-3);
  3. reads ``b0`` from the cell definition, ``bM`` from the up-rules' saturated
     value, ``bm`` from the down-rules' saturated value;
  4. returns ``b = combine(b0, bM, bm, U, D)`` (eq 4).

Mixed response types (hill / capped-linear / step) compose in the multivariate
sum by converting each rule's single-signal response ``r`` to an equivalent
Hill term ``r / (1 - r)``; a pure-Hill group reduces exactly to eqs 2-4.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

from .behaviors import ResolvedBehavior
from .cell_definitions import CellDefinitions
from .response import (
    combine,
    hill_term,
    linear_response,
    saturating_response,
    step_response,
)
from .rules import HILL, LINEAR, STEP, Rule, RuleSet
from .signals import SignalContext, resolve_signal

__all__ = ["BehaviorValue", "GrammarEvaluator"]

# A near-saturating term used to represent a fired step / fully-on response.
_BIG_TERM = 1.0e9


@dataclass(frozen=True)
class BehaviorValue:
    """One evaluated behavior for a cell."""

    canonical: str          # e.g. "cycle entry" or "transform to motile tumor"
    value: float
    is_rate: bool
    behavior_name: str      # the base behavior, e.g. "transform to"
    target: Optional[str]   # parsed target, e.g. "motile tumor" (None if not a family)


def _rule_term(rule: Rule, signal_value: float) -> float:
    """The additive multivariate term contributed by one rule's signal."""
    if rule.response_type == HILL:
        return hill_term(signal_value, rule.half_max, rule.hill_power)
    if rule.response_type == LINEAR:
        r = linear_response(signal_value, rule.half_max)
    elif rule.response_type == STEP:
        r = step_response(signal_value, rule.half_max)
    else:  # pragma: no cover - guarded at Rule construction
        raise ValueError(f"unknown response_type {rule.response_type!r}")
    r = min(max(r, 0.0), 1.0 - 1.0e-12)
    return r / (1.0 - r)


class GrammarEvaluator:
    """Evaluate behavior values from signals for a single cell type's rules."""

    def __init__(self, ruleset: RuleSet, cell_definitions: Optional[CellDefinitions] = None):
        self.ruleset = ruleset
        self.cell_definitions = cell_definitions or CellDefinitions(registry=ruleset.registry)

    def evaluate(self, context: SignalContext) -> Dict[str, BehaviorValue]:
        """Return ``{canonical_behavior -> BehaviorValue}`` for ``context``."""
        out: Dict[str, BehaviorValue] = {}
        for group in self.ruleset.groups_for(context.cell_type):
            resolved: ResolvedBehavior = group["resolved"]
            canonical: str = group["canonical"]
            rules: List[Rule] = group["rules"]

            up_terms: List[float] = []
            down_terms: List[float] = []
            up_sat: List[float] = []
            down_sat: List[float] = []

            for rule in rules:
                if context.dead and not rule.applies_to_dead:
                    continue
                sig = resolve_signal(rule.signal, context)
                term = _rule_term(rule, sig)
                if rule.increases:
                    up_terms.append(term)
                    up_sat.append(rule.saturated_value)
                else:
                    down_terms.append(term)
                    down_sat.append(rule.saturated_value)

            b0 = self.cell_definitions.base(context.cell_type, canonical)
            # bM: up-regulated target (saturated max). None present -> no up effect.
            b_max = max(up_sat) if up_sat else b0
            # bm: down-regulated target (saturated min). None present -> fall back.
            b_min = min(down_sat) if down_sat else self.cell_definitions.min(context.cell_type, canonical)

            up = saturating_response(up_terms)
            down = saturating_response(down_terms)
            value = combine(b0, b_max, b_min, up, down)

            out[canonical] = BehaviorValue(
                canonical=canonical,
                value=value,
                is_rate=resolved.is_rate,
                behavior_name=resolved.behavior.name,
                target=resolved.target,
            )
        return out
