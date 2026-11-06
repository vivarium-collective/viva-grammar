"""The grammar's dictionary of **signals** (the nouns a cell senses).

A *signal* is a stimulus a cell can read to drive its behaviors: a substrate
concentration or gradient, contact with cells of a given type, mechanical
pressure, volume, damage, dead status, or simulation time.

The engine is ABM-agnostic, so a host supplies a :class:`SignalContext` of raw
per-cell readouts; a signal *resolver* maps a grammar signal name (the string
in a rule) to a scalar value from that context.  Host adapters
(reference substrate, blueprint) are responsible only for *populating* the
context from their own state — the signal vocabulary lives here.

Clean-room reimplementation of a subset of PhysiCell's signal dictionary
(Methods S1).  See ``references/PROVENANCE.md``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional

__all__ = ["SignalContext", "resolve_signal"]


@dataclass
class SignalContext:
    """Raw per-cell readouts a host ABM provides for one cell, one step.

    All fields are optional; a signal referencing a missing readout resolves to
    0.0 (treated as "signal absent").  Cell ids are the host's own; this context
    describes a single cell already identified by the caller.
    """

    cell_type: str
    fields: Dict[str, float] = field(default_factory=dict)        # substrate conc at cell
    gradients: Dict[str, float] = field(default_factory=dict)     # |grad substrate| at cell
    neighbor_count: Dict[str, int] = field(default_factory=dict)  # contacts, by neighbor type
    contact_area: Dict[str, float] = field(default_factory=dict)  # shared surface, by type
    volume: float = 1.0
    pressure: float = 0.0
    damage: float = 0.0
    dead: bool = False
    time: float = 0.0
    custom: Dict[str, float] = field(default_factory=dict)        # escape hatch


# Signal names with fixed meaning (not substrate/contact families).
_SCALAR_SIGNALS = {
    "pressure": lambda c: float(c.pressure),
    "volume": lambda c: float(c.volume),
    "damage": lambda c: float(c.damage),
    "dead": lambda c: 1.0 if c.dead else 0.0,
    "time": lambda c: float(c.time),
}


def resolve_signal(signal_name: str, context: SignalContext) -> float:
    """Resolve a grammar signal string to its scalar value for ``context``.

    Resolution order:
      1. fixed scalar signals (pressure, volume, damage, dead, time);
      2. ``"contact with <type>"``        -> neighbor count of that type;
      3. ``"<substrate> gradient"``       -> gradient magnitude of that field;
      4. ``"intracellular <substrate>"``  -> custom[...] (host-provided);
      5. a bare name                      -> substrate concentration ``fields[name]``;
      6. finally ``custom[name]``, else 0.0.
    """
    name = signal_name.strip()

    if name in _SCALAR_SIGNALS:
        return _SCALAR_SIGNALS[name](context)

    if name.startswith("contact with "):
        ctype = name[len("contact with "):].strip()
        return float(context.neighbor_count.get(ctype, 0))

    if name.endswith(" gradient"):
        sub = name[: -len(" gradient")].strip()
        return float(context.gradients.get(sub, 0.0))

    if name.startswith("intracellular "):
        sub = name[len("intracellular "):].strip()
        return float(context.custom.get(f"intracellular {sub}", context.custom.get(sub, 0.0)))

    if name in context.fields:
        return float(context.fields[name])

    if name in context.custom:
        return float(context.custom[name])

    # Unknown / absent signal -> treated as 0 (no effect).
    return 0.0
