"""Reference-substrate helpers.

The reference substrate (``viva_grammar.reference``) emits the same cohort
readout ports :class:`GrammarRulesProcess` consumes, so it needs no bespoke
adapter Process — it wires directly.  This module provides a small helper to
build a :class:`SignalContext` from a plain per-cell dict, handy for tests,
notebooks, and driving the engine outside a Composite.
"""
from __future__ import annotations

from typing import Optional

from ...engine import SignalContext

__all__ = ["context_from_dict"]


def context_from_dict(cell: dict) -> SignalContext:
    """Build a :class:`SignalContext` from a flat cell dict.

    Keys: ``type`` (required), ``fields``, ``gradients``, ``neighbor_count``,
    ``contact_area``, ``volume``, ``pressure``, ``damage``, ``dead``, ``time``.
    """
    return SignalContext(
        cell_type=cell["type"],
        fields=dict(cell.get("fields", {})),
        gradients=dict(cell.get("gradients", {})),
        neighbor_count=dict(cell.get("neighbor_count", {})),
        contact_area=dict(cell.get("contact_area", {})),
        volume=float(cell.get("volume", 1.0)),
        pressure=float(cell.get("pressure", 0.0)),
        damage=float(cell.get("damage", 0.0)),
        dead=bool(cell.get("dead", False)),
        time=float(cell.get("time", 0.0)),
    )
