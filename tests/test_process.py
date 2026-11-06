"""Unit tests for GrammarRulesProcess: signal sampling, rate->event routing,
and the {fates, scalar, divisions, deaths} output surface.
"""
from __future__ import annotations

import numpy as np
import pytest
from process_bigraph import allocate_core

from viva_grammar.steps import BlueprintGrammarProcess, GrammarRulesProcess


def _proc(rule_rows, cell_defs=None, **cfg):
    core = allocate_core()
    config = {"rules": rule_rows, "cell_definitions": cell_defs or {}, "seed": 0, **cfg}
    return GrammarRulesProcess(config=config, core=core)


def _rule(cell_type, signal, direction, behavior, sat, half, power=1, dead=0):
    return {
        "cell_type": cell_type, "signal": signal, "direction": direction,
        "behavior": behavior, "saturated_value": sat, "half_max": half,
        "hill_power": power, "applies_to_dead": dead,
    }


# One cell "0" of type A at a grid site reading oxygen = 1000 (far above half-max).
def _cohort(cell_type="A", oxygen=1000.0):
    grid = np.full((3, 3), oxygen, dtype=float)
    return {
        "types": {"0": cell_type},
        "positions": {"0": [1, 1]},
        "sizes": {"0": 1.0},
        "neighbors": {"0": []},
        "shared_surface_by_type": {"0": {}},
        "alive": {"0": True},
        "field_grid": {"oxygen": grid},
        "pressure": {"0": 0.0},
        "time": 0.0,
    }


def test_cycle_entry_fires_division():
    p = _proc([_rule("A", "oxygen", "increases", "cycle entry", 1.0, 1.0)])
    out = p.update(_cohort(), interval=10000.0)  # huge interval -> Pr ~ 1
    assert out["divisions"] == [{"cell_id": "0"}]
    assert out["deaths"] == [] and out["fates"] == {}


def test_necrosis_fires_death():
    p = _proc([_rule("A", "oxygen", "increases", "necrosis", 1.0, 1.0)])
    out = p.update(_cohort(), interval=10000.0)
    assert out["deaths"] == [{"cell_id": "0", "cause": "necrosis"}]
    assert out["divisions"] == []


def test_transform_fires_fate_with_target():
    p = _proc([_rule("A", "oxygen", "increases", "transform to B", 1.0, 1.0)])
    out = p.update(_cohort(), interval=10000.0)
    assert out["fates"] == {"0": "B"}


def test_parameter_behavior_emits_scalar_not_event():
    p = _proc([_rule("A", "oxygen", "increases", "migration speed", 2.5, 1.0)])
    out = p.update(_cohort(oxygen=1000.0), interval=1.0)
    # migration speed is a continuous parameter -> scalar, no lineage event.
    # Hill: U = 1000/(1+1000); value = saturated * U.
    assert out["divisions"] == [] and out["deaths"] == [] and out["fates"] == {}
    assert out["scalar"]["0"]["migration speed"] == pytest.approx(2.5 * 1000.0 / 1001.0, rel=1e-9)


def test_death_preempts_division_same_step():
    # Both a certain death and a certain division; death wins (one outcome/cell).
    p = _proc([
        _rule("A", "oxygen", "increases", "necrosis", 1.0, 1.0),
        _rule("A", "oxygen", "increases", "cycle entry", 1.0, 1.0),
    ])
    out = p.update(_cohort(), interval=10000.0)
    assert out["deaths"] == [{"cell_id": "0", "cause": "necrosis"}]
    assert out["divisions"] == []


def test_low_signal_no_events():
    p = _proc([_rule("A", "oxygen", "increases", "cycle entry", 1.0, 1.0)])
    out = p.update(_cohort(oxygen=0.0), interval=10000.0)  # oxygen 0 -> rate 0
    assert out["divisions"] == [] and out["deaths"] == []


def test_dead_cells_skipped():
    p = _proc([_rule("A", "oxygen", "increases", "cycle entry", 1.0, 1.0)])
    c = _cohort(); c["alive"]["0"] = False
    out = p.update(c, interval=10000.0)
    assert out["divisions"] == []


def test_rate_mode_per_step_ignores_interval():
    # per_step: Pr = 1 - exp(-rate); a tiny rate stays tiny regardless of interval.
    p = _proc([_rule("A", "oxygen", "increases", "cycle entry", 1e-6, 1.0)],
              rate_mode="per_step")
    fired = 0
    for seed in range(200):
        pp = _proc([_rule("A", "oxygen", "increases", "cycle entry", 1e-6, 1.0)],
                   rate_mode="per_step", seed=seed)
        if pp.update(_cohort(), interval=10000.0)["divisions"]:
            fired += 1
    assert fired == 0  # per-step prob ~1e-6, never fires across 200 seeds


def test_blueprint_adapter_defaults_per_step():
    core = allocate_core()
    bp = BlueprintGrammarProcess(config={"rules": [_rule("A", "oxygen", "increases", "cycle entry", 1.0, 1.0)]}, core=core)
    assert bp.config["rate_mode"] == "per_step"
    # same output surface as RulesProcess
    assert set(bp.outputs()) == {"fates", "scalar", "divisions", "deaths"}


def test_field_sampling_world_to_grid():
    # positions in world coords; origin/spacing map them onto grid indices.
    grid = np.zeros((5, 5)); grid[2, 2] = 100.0
    p = _proc([_rule("A", "oxygen", "increases", "migration speed", 1.0, 1.0)],
              field_sampling={"origin": [-50.0, -50.0], "spacing": 20.0})
    cohort = {
        "types": {"0": "A"}, "positions": {"0": [-10.0, -10.0]},  # -> index (2,2)
        "sizes": {"0": 1.0}, "neighbors": {"0": []},
        "shared_surface_by_type": {"0": {}}, "alive": {"0": True},
        "field_grid": {"oxygen": grid}, "pressure": {"0": 0.0}, "time": 0.0,
    }
    out = p.update(cohort, interval=1.0)
    # oxygen sampled as 100 at index (2,2); Hill U = 100/101, value = 1.0 * U.
    assert out["scalar"]["0"]["migration speed"] == pytest.approx(100.0 / 101.0, rel=1e-9)
