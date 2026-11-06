"""Blueprint adapter: the grammar engine as a drop-in for RulesProcess.

``viva-multicell-blueprint`` separates its shared mechanism (fields, rules,
immune, ODE) from swappable mechanical substrates (CPM / vertex / center-based)
behind one ``MechanicsBackend`` port.  Its built-in rule backend is
``RulesProcess``, which reads the backend's cohort readouts
(``types``/``neighbors``/``shared_surface_by_type``/``positions``/``sizes``/
``alive``/``field_grid``) and writes ``fates``/``scalar``/``divisions``/
``deaths``.

:class:`BlueprintGrammarProcess` is :class:`GrammarRulesProcess` with the same
port surface, so it wires into the blueprint composite exactly where
``RulesProcess`` does — but expresses behavior with the *published Macklin
grammar* (Hill-function signals -> behaviors) instead of the blueprint's ad-hoc
``rate_expr`` grammar.

Parity note
-----------
The blueprint's ``RulesProcess`` applies ``Pr = 1 - exp(-rate)`` and ignores the
step interval (its ``rate_expr`` values are already per-step).  So this adapter
defaults to ``rate_mode="per_step"`` to reproduce that convention exactly.
Confirm the blueprint composite's effective step duration before a numeric
parity run; switch to ``rate_mode="per_time"`` to use physical per-time rates
with the real interval.

It does **not** emit ``activations`` — that is a MechanicsBackend *input* fed by
``ImmuneRecruitmentProcess``, not a rule-engine output.

Expressiveness boundary
-----------------------
The grammar combines rules independently per ``(cell_type, behavior)`` via the
up/down bilinear math (eq 4); it has **no rule ordering** and **no
scalar-written-then-read chaining**, unlike ``RulesProcess`` (ordered rules +
``resistance`` scalars + ``const``).  Use this adapter for models expressible in
the grammar (signal -> transform/divide/die, continuous parameters); it is not a
general replacement for ordered/resistance rule chains.
"""
from __future__ import annotations

from ..grammar_process import GrammarRulesProcess

__all__ = ["BlueprintGrammarProcess"]


class BlueprintGrammarProcess(GrammarRulesProcess):
    """GrammarRulesProcess configured as a blueprint RulesProcess drop-in."""

    config_schema = {
        **GrammarRulesProcess.config_schema,
        # Match RulesProcess's per-step `1 - exp(-rate)` convention by default.
        "rate_mode": {"_type": "string", "_default": "per_step"},
    }
