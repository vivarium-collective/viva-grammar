"""GrammarRulesProcess: the ABM-agnostic grammar rule engine as a pbg Process.

Each step it reads a per-step **snapshot of the whole cell cohort** (types,
neighbors, contacts, positions, sizes, alive, and the shared substrate fields),
builds a :class:`SignalContext` per cell, evaluates the grammar
(:class:`GrammarEvaluator`), and emits the behavioral outcome on four ports:

- ``fates``     — ``cell_id -> new_type``   (from "transform to <type>" behaviors)
- ``divisions`` — ``[{"cell_id": id}, ...]`` (from the "cycle entry" rate)
- ``deaths``    — ``[{"cell_id": id}, ...]`` (from apoptosis / necrosis / death rates)
- ``scalar``    — ``cell_id -> {behavior -> value}`` (continuous parameters:
  migration speed/bias, secretion/uptake rates, …)

These are **exactly** the output ports of ``viva_multicell_blueprint``'s
``RulesProcess`` (``fates`` / ``scalar`` / ``divisions`` / ``deaths``), and the
input ports are its readout surface, so this is a drop-in alternative rule
backend for any MechanicsBackend — and for any ABM that produces the same
readouts (e.g. this package's reference substrate).

Rate behaviors become stochastic events via ``Pr = 1 - exp(-rate * interval)``.
``rate_mode="per_step"`` treats the interval as 1 (matching ``RulesProcess``'s
``1 - exp(-rate)`` convention exactly); ``"per_time"`` (default) uses the real
``interval`` with physical (per-time) rates.

Clean-room reimplementation of the grammar; see ``references/PROVENANCE.md``.
"""
from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
from process_bigraph import Process

from ..engine import CellDefinitions, GrammarEvaluator, Rule, RuleSet, SignalContext
from ..io import parse_rules_csv

__all__ = ["GrammarRulesProcess"]

# Behavior names (canonical) that map to lineage events rather than scalars.
_DEATH_BEHAVIORS = {"apoptosis", "necrosis", "death"}
_DIVISION_BEHAVIORS = {"cycle entry"}
_TRANSFORM_PREFIXES = ("transform to", "transition to")


def _rules_from_config(config: dict) -> List[Rule]:
    rules_csv = (config.get("rules_csv") or "").strip()
    if rules_csv:
        # A path has no newline/comma and names an existing file; otherwise it's CSV text.
        looks_like_path = ("\n" not in rules_csv and "," not in rules_csv and len(rules_csv) < 1024)
        if looks_like_path and Path(rules_csv).exists():
            return parse_rules_csv(Path(rules_csv).read_text())
        return parse_rules_csv(rules_csv)
    inline = config.get("rules") or []
    out: List[Rule] = []
    for r in inline:
        out.append(Rule(
            cell_type=r["cell_type"], signal=r["signal"], direction=r["direction"],
            behavior=r["behavior"], saturated_value=float(r["saturated_value"]),
            half_max=float(r["half_max"]), hill_power=float(r["hill_power"]),
            applies_to_dead=bool(r.get("applies_to_dead", False)),
            response_type=r.get("response_type", "hill"),
        ))
    return out


class GrammarRulesProcess(Process):
    config_schema = {
        "rules_csv": {"_type": "string", "_default": ""},        # CSV text or a file path
        "rules": {"_type": "list[node]", "_default": []},         # inline rule dicts (alt.)
        "cell_definitions": {"_type": "node", "_default": {}},    # cell_type -> {behavior -> {base,min}}
        "seed": {"_type": "integer", "_default": 0},
        "rate_mode": {"_type": "string", "_default": "per_time"}, # "per_time" | "per_step"
        "default_interval": {"_type": "float", "_default": 1.0},
        # world->grid index mapping for field sampling: {"origin":[x0,y0], "spacing":dx}
        "field_sampling": {"_type": "node", "_default": {}},
    }

    def initialize(self, config):
        self._rng = np.random.default_rng(int(config.get("seed", 0)))
        ruleset = RuleSet(_rules_from_config(config))
        cell_defs = CellDefinitions(config.get("cell_definitions") or {}, registry=ruleset.registry)
        self._evaluator = GrammarEvaluator(ruleset, cell_defs)
        self._rate_mode = config.get("rate_mode", "per_time")
        self._default_interval = float(config.get("default_interval", 1.0))
        fs = config.get("field_sampling") or {}
        self._origin = list(fs.get("origin", [])) or None
        self._spacing = float(fs["spacing"]) if fs.get("spacing") else None

    def inputs(self):
        return {
            "positions": "map[list[float]]",
            "sizes": "map[float]",
            "types": "map[string]",
            "shared_surface_by_type": "map[map[float]]",
            "neighbors": "map[list[string]]",
            "alive": "map[boolean]",
            "field_grid": "map[node]",
            "pressure": "map[float]",   # optional; absent -> neighbor-density proxy
            "time": "float",
        }

    def outputs(self):
        # Same surface as viva_multicell_blueprint RulesProcess; overwrite so the
        # per-step result maps replace the shared stores whole (no stale carry-over).
        return {
            "fates": "overwrite[map[string]]",
            "scalar": "overwrite[map[map[float]]]",
            "divisions": "list[node]",
            "deaths": "list[node]",
        }

    # -- field sampling --------------------------------------------------------
    def _site_index(self, grid: np.ndarray, pos) -> tuple:
        if self._origin is not None and self._spacing:
            coords = [(float(p) - o) / self._spacing for p, o in zip(pos, self._origin)]
        else:
            coords = [float(p) for p in pos]
        return tuple(int(min(max(round(c), 0), n - 1)) for c, n in zip(coords, grid.shape))

    def _field_value(self, grid: np.ndarray, pos) -> float:
        return float(grid[self._site_index(grid, pos)])

    def _gradient_mag(self, grid: np.ndarray, pos) -> float:
        idx = self._site_index(grid, pos)
        g2 = 0.0
        for ax, n in enumerate(grid.shape):
            lo = list(idx); hi = list(idx)
            lo[ax] = max(idx[ax] - 1, 0); hi[ax] = min(idx[ax] + 1, n - 1)
            span = hi[ax] - lo[ax]
            if span > 0:
                d = (float(grid[tuple(hi)]) - float(grid[tuple(lo)])) / span
                g2 += d * d
        return math.sqrt(g2)

    # -- main ------------------------------------------------------------------
    def update(self, inputs, interval=None):
        types = inputs.get("types") or {}
        neighbors = inputs.get("neighbors") or {}
        shared = inputs.get("shared_surface_by_type") or {}
        positions = inputs.get("positions") or {}
        sizes = inputs.get("sizes") or {}
        alive = inputs.get("alive") or {}
        pressure_in = inputs.get("pressure") or {}
        sim_time = float(inputs.get("time", 0.0) or 0.0)
        grids = {k: np.asarray(v, dtype=float) for k, v in (inputs.get("field_grid") or {}).items()}

        if interval is None:
            interval = self._default_interval
        eff_interval = 1.0 if self._rate_mode == "per_step" else float(interval)

        cids = sorted(c for c in types if alive.get(c, True))

        fates: Dict[str, str] = {}
        scalar: Dict[str, Dict[str, float]] = {}
        divisions: List[dict] = []
        deaths: List[dict] = []
        decided = set()  # cells that already took a lineage/fate outcome this step

        for cid in cids:
            pos = positions.get(cid, [0.0, 0.0])
            ncount: Dict[str, int] = defaultdict(int)
            for j in neighbors.get(cid) or []:
                if alive.get(j, True) and j in types:
                    ncount[types[j]] += 1
            # pressure: host-provided, else a simple neighbor-crowding proxy.
            if cid in pressure_in:
                pressure = float(pressure_in[cid])
            else:
                pressure = float(sum(ncount.values()))
            ctx = SignalContext(
                cell_type=types[cid],
                fields={n: self._field_value(g, pos) for n, g in grids.items()},
                gradients={n: self._gradient_mag(g, pos) for n, g in grids.items()},
                neighbor_count=dict(ncount),
                contact_area=dict(shared.get(cid) or {}),
                volume=float(sizes.get(cid, 1.0)),
                pressure=pressure,
                dead=not alive.get(cid, True),
                time=sim_time,
            )
            evaluated = self._evaluator.evaluate(ctx)

            fired_death = fired_transform = fired_division = None
            for canonical, bv in evaluated.items():
                if not bv.is_rate:
                    scalar.setdefault(cid, {})[canonical] = bv.value
                    continue
                rate = max(bv.value, 0.0)
                prob = 1.0 - math.exp(-rate * eff_interval)
                fired = self._rng.random() < prob
                if not fired:
                    continue
                if bv.behavior_name in _DEATH_BEHAVIORS:
                    fired_death = canonical
                elif bv.behavior_name in _TRANSFORM_PREFIXES and bv.target:
                    fired_transform = bv.target
                elif canonical in _DIVISION_BEHAVIORS:
                    fired_division = canonical

            # One lineage/fate outcome per cell: death > transform > division.
            if fired_death is not None:
                # `cause` (necrosis/apoptosis/death) is extra metadata; the
                # blueprint's deaths port ignores it, hosts that care can use it.
                deaths.append({"cell_id": cid, "cause": fired_death}); decided.add(cid)
            elif fired_transform is not None:
                fates[cid] = fired_transform; decided.add(cid)
            elif fired_division is not None:
                divisions.append({"cell_id": cid}); decided.add(cid)

        return {"fates": fates, "scalar": scalar, "divisions": divisions, "deaths": deaths}
