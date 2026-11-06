"""Read/write the PhysiCell 8-column cell-behavior rules CSV.

Column layout (the canonical grammar interchange format):

    cell_type, signal, response, behavior, saturated_value, half_max, hill_power, applies_to_dead

``response`` is ``increases`` or ``decreases``.  ``//`` lines and blank lines
are comments (as in the reference ``grammar_samples`` CSVs).  See
``references/PROVENANCE.md``.
"""
from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import List, Optional, Union

from ..engine.rules import DECREASES, INCREASES, Rule

__all__ = ["parse_rules_csv", "load_rules_csv", "rules_to_csv"]

_N_COLS = 8


def _parse_bool(token: str) -> bool:
    t = token.strip().lower()
    return t in ("1", "true", "yes", "y")


def parse_rules_csv(text: str) -> List[Rule]:
    """Parse rules from CSV text; skip ``//`` comments and blank lines."""
    rules: List[Rule] = []
    reader = csv.reader(io.StringIO(text))
    for lineno, row in enumerate(reader, start=1):
        if not row:
            continue
        first = row[0].strip()
        if not first or first.startswith("//") or first.startswith("#"):
            continue
        cells = [c.strip() for c in row]
        if len(cells) < _N_COLS:
            raise ValueError(
                f"rules CSV line {lineno}: expected {_N_COLS} columns, got {len(cells)}: {row!r}"
            )
        cell_type, signal, response, behavior = cells[0:4]
        response = response.lower()
        if response not in (INCREASES, DECREASES):
            raise ValueError(
                f"rules CSV line {lineno}: response must be 'increases'/'decreases', got {response!r}"
            )
        try:
            saturated_value = float(cells[4])
            half_max = float(cells[5])
            hill_power = float(cells[6])
        except ValueError as exc:
            raise ValueError(f"rules CSV line {lineno}: non-numeric parameter: {exc}") from exc
        rules.append(
            Rule(
                cell_type=cell_type,
                signal=signal,
                direction=response,
                behavior=behavior,
                saturated_value=saturated_value,
                half_max=half_max,
                hill_power=hill_power,
                applies_to_dead=_parse_bool(cells[7]),
            )
        )
    return rules


def load_rules_csv(path: Union[str, Path]) -> List[Rule]:
    """Load rules from a CSV file path."""
    return parse_rules_csv(Path(path).read_text())


def rules_to_csv(rules: List[Rule]) -> str:
    """Serialize rules back to the 8-column CSV format."""
    buf = io.StringIO()
    writer = csv.writer(buf)
    for r in rules:
        writer.writerow([
            r.cell_type, r.signal, r.direction, r.behavior,
            r.saturated_value, r.half_max,
            int(r.hill_power) if float(r.hill_power).is_integer() else r.hill_power,
            int(r.applies_to_dead),
        ])
    return buf.getvalue()
