"""Unit tests for the grammar engine — the math is checked against hand values
and against the paper's hypoxia rules (``references/physicell_hypoxia_cell_rules.csv``).
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from viva_grammar.engine import (
    CellDefinitions,
    GrammarEvaluator,
    RuleSet,
    SignalContext,
    combine,
    hill_response,
    hill_term,
    linear_response,
    saturating_response,
    step_response,
)
from viva_grammar.io import load_rules_csv, parse_rules_csv, parse_statement, rules_to_csv

REFCSV = Path(__file__).resolve().parents[1] / "references" / "physicell_hypoxia_cell_rules.csv"


# ---------------------------------------------------------------- response math
def test_hill_half_max_is_one_half():
    assert hill_response(5.0, 5.0, 4) == pytest.approx(0.5)
    assert hill_response(0.0, 5.0, 4) == pytest.approx(0.0)


def test_hill_monotonic_and_saturating():
    vals = [hill_response(v, 5.0, 4) for v in (0, 2.5, 5, 10, 1000)]
    assert vals == sorted(vals)
    assert vals[-1] == pytest.approx(1.0, abs=1e-6)


def test_hill_power_sharpens():
    # Higher Hill power -> sharper switch: smaller response below half-max.
    assert hill_response(3.0, 5.0, 8) < hill_response(3.0, 5.0, 2)


def test_linear_and_step():
    assert linear_response(5.0, 5.0) == pytest.approx(0.5)
    assert linear_response(10.0, 5.0) == pytest.approx(1.0)
    assert linear_response(100.0, 5.0) == pytest.approx(1.0)  # capped
    assert step_response(4.9, 5.0) == 0.0
    assert step_response(5.0, 5.0) == 1.0


def test_saturating_response_reduces_to_single_hill():
    t = hill_term(7.0, 5.0, 4)
    assert saturating_response([t]) == pytest.approx(hill_response(7.0, 5.0, 4))


def test_combine_single_up():
    # b0=0, bM=1, up=0.5, no down -> 0.5
    assert combine(0.0, 1.0, 0.0, up=0.5, down=0.0) == pytest.approx(0.5)


def test_combine_down_pulls_to_min():
    # full down -> min, regardless of up
    assert combine(2.0, 5.0, 0.3, up=1.0, down=1.0) == pytest.approx(0.3)


# ---------------------------------------------------------------- CSV parsing
def test_parse_reference_hypoxia_csv():
    rules = load_rules_csv(REFCSV)
    assert len(rules) == 8
    cts = {r.cell_type for r in rules}
    assert cts == {"tumor", "motile tumor"}
    # round-trips back to CSV and re-parses to the same rules
    reparsed = parse_rules_csv(rules_to_csv(rules))
    assert [r.behavior for r in reparsed] == [r.behavior for r in rules]


def test_csv_skips_comments():
    text = "// a comment\n\ntumor,oxygen,increases,cycle entry,7.2e-3,21.5,4,0\n"
    rules = parse_rules_csv(text)
    assert len(rules) == 1 and rules[0].behavior == "cycle entry"


def test_english_parser():
    r = parse_statement("In malignant epithelial cells, cisplatin increases apoptosis.")
    assert r.cell_type == "malignant epithelial cell"
    assert r.signal == "cisplatin"
    assert r.direction == "increases"
    assert r.behavior == "apoptosis"


def test_unknown_behavior_raises():
    with pytest.raises(KeyError):
        RuleSet(parse_rules_csv("tumor,oxygen,increases,teleport,0,1,2,0"))


# ---------------------------------------------------------------- evaluator / hypoxia
def _hypoxia_cell_defs():
    # Baselines from the hypoxia PhysiCell_settings.xml (references/PROVENANCE.md).
    base = {
        "cycle entry": {"base": 7.2e-4},
        "apoptosis": {"base": 5.31667e-5},
        "necrosis": {"base": 2.8e-3},
    }
    return CellDefinitions({
        "tumor": {**base,
                  "migration speed": {"base": 0.0},
                  "transform to motile tumor": {"base": 0.001}},
        "motile tumor": {**base,
                         "migration speed": {"base": 0.47},
                         "transform to tumor": {"base": 0.0}},
    })


def _hypoxia_evaluator():
    rs = RuleSet(load_rules_csv(REFCSV))
    return GrammarEvaluator(rs, _hypoxia_cell_defs())


def test_cycle_entry_at_half_max_oxygen():
    ev = _hypoxia_evaluator()
    ctx = SignalContext(cell_type="tumor", fields={"oxygen": 21.5}, pressure=0.0)
    out = ev.evaluate(ctx)
    # U=0.5 (oxygen at half-max, power 4), D=0 (pressure 0). base 7.2e-4, max 7.2e-3.
    expected = 0.5 * 7.2e-4 + 0.5 * 7.2e-3
    assert out["cycle entry"].value == pytest.approx(expected, rel=1e-9)
    assert out["cycle entry"].is_rate


def test_necrosis_high_at_low_oxygen():
    ev = _hypoxia_evaluator()
    anoxic = ev.evaluate(SignalContext(cell_type="tumor", fields={"oxygen": 0.0}))
    normoxic = ev.evaluate(SignalContext(cell_type="tumor", fields={"oxygen": 100.0}))
    # necrosis = base at zero oxygen, -> ~0 at high oxygen (oxygen decreases necrosis)
    assert anoxic["necrosis"].value == pytest.approx(2.8e-3, rel=1e-6)
    assert normoxic["necrosis"].value < 1e-6
    # half-max check
    half = ev.evaluate(SignalContext(cell_type="tumor", fields={"oxygen": 3.75}))
    assert half["necrosis"].value == pytest.approx(0.5 * 2.8e-3, rel=1e-6)


def test_transform_to_motile_high_at_low_oxygen():
    ev = _hypoxia_evaluator()
    out = ev.evaluate(SignalContext(cell_type="tumor", fields={"oxygen": 0.0}))
    bv = out["transform to motile tumor"]
    assert bv.value == pytest.approx(0.001, rel=1e-6)   # base, full at anoxia
    assert bv.behavior_name == "transform to"
    assert bv.target == "motile tumor"
    assert bv.is_rate
    # vanishes at high oxygen
    hi = ev.evaluate(SignalContext(cell_type="tumor", fields={"oxygen": 100.0}))
    assert hi["transform to motile tumor"].value < 1e-6


def test_motile_speed_is_baseline_property():
    # The grammar has no migration-speed rule for the hypoxia model; the motile
    # type's speed is a cell-definition baseline, surfaced via evaluation only if
    # a rule references it. Here we assert the baseline is sourced correctly.
    cd = _hypoxia_cell_defs()
    assert cd.base("motile tumor", "migration speed") == pytest.approx(0.47)
    assert cd.base("tumor", "migration speed") == pytest.approx(0.0)


def test_pressure_suppresses_cycle_entry():
    ev = _hypoxia_evaluator()
    lowp = ev.evaluate(SignalContext(cell_type="tumor", fields={"oxygen": 100.0}, pressure=0.0))
    highp = ev.evaluate(SignalContext(cell_type="tumor", fields={"oxygen": 100.0}, pressure=1.0))
    # pressure decreases cycle entry -> high pressure suppresses it
    assert highp["cycle entry"].value < lowp["cycle entry"].value
