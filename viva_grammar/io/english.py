"""Convenience parser for the grammar's English behavioral statements.

Parses the human-readable form the paper leads with:

    In [cell type T], [signal S] increases|decreases [behavior B].

e.g. ``"In malignant epithelial cells, cisplatin increases apoptosis."`` or
``"In naive T cells, IL-10 decreases transition to CD8+ T cell."``

The CSV remains the canonical interchange format (it carries the numeric
parameters); this parser produces :class:`Rule` objects with the numeric slots
left to defaults or supplied via ``defaults=``, for quickly sketching rules in
prose.  Optional-argument grammar (the paper's "[with optional arguments]") is
out of scope for v1.
"""
from __future__ import annotations

import re
from typing import List, Optional

from ..engine.rules import DECREASES, INCREASES, Rule

__all__ = ["parse_statement", "parse_statements"]

# "In <T>, <S> increases|decreases <B>."  (leading "In" optional; trailing period optional)
_PATTERN = re.compile(
    r"^\s*(?:in\s+)?(?P<cell>.+?)\s*,\s*(?P<signal>.+?)\s+(?P<dir>increases|decreases)\s+(?P<behavior>.+?)\s*\.?\s*$",
    re.IGNORECASE,
)

# Light singular-ization of the cell-type phrase ("cells" -> "cell"); the
# grammar's type names are singular ("malignant epithelial cell").
_PLURAL_TAIL = re.compile(r"\bcells\b\s*$", re.IGNORECASE)


def _normalize_cell(phrase: str) -> str:
    phrase = phrase.strip()
    phrase = _PLURAL_TAIL.sub("cell", phrase)
    return phrase


def parse_statement(
    text: str,
    *,
    saturated_value: float = 0.0,
    half_max: float = 0.5,
    hill_power: float = 4.0,
    applies_to_dead: bool = False,
    normalize_cell_type: bool = True,
) -> Rule:
    """Parse one English statement into a :class:`Rule`.

    Numeric slots (saturated value, half-max, Hill power, dead flag) default to
    the keyword arguments — prose carries no numbers.
    """
    m = _PATTERN.match(text)
    if not m:
        raise ValueError(
            f"could not parse grammar statement: {text!r}\n"
            "expected: 'In <cell type>, <signal> increases|decreases <behavior>'"
        )
    cell = m.group("cell")
    if normalize_cell_type:
        cell = _normalize_cell(cell)
    return Rule(
        cell_type=cell,
        signal=m.group("signal").strip(),
        direction=m.group("dir").lower(),
        behavior=m.group("behavior").strip(),
        saturated_value=saturated_value,
        half_max=half_max,
        hill_power=hill_power,
        applies_to_dead=applies_to_dead,
    )


def parse_statements(text: str, **defaults) -> List[Rule]:
    """Parse multiple statements, one per non-empty / non-comment line."""
    rules: List[Rule] = []
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("//"):
            continue
        rules.append(parse_statement(s, **defaults))
    return rules
