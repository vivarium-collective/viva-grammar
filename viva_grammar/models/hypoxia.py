"""The paper's hypoxia -> motility model, as reusable grammar objects.

Bundles the verbatim PhysiCell hypoxia rules (``references/physicell_hypoxia_cell_rules.csv``)
with the cell-definition baselines extracted from that model's
``PhysiCell_settings.xml`` (see ``references/PROVENANCE.md``).  Both the
standalone reference study and the blueprint adapter import these.

Encoded faithfully as the paper describes it: *transition-based* — "oxygen
decreases transform to motile tumor" drives a `tumor -> motile tumor` transition
at low oxygen; the elevated motility of the motile phenotype is a baseline
property of the ``motile tumor`` cell type (base migration speed 0.47 vs 0.0),
not a migration-speed rule.
"""
from __future__ import annotations

from pathlib import Path
from typing import List

from ..engine import CellDefinitions, Rule, RuleSet
from ..io import load_rules_csv

__all__ = ["hypoxia_rules", "hypoxia_cell_definitions", "hypoxia_ruleset", "HYPOXIA_RULES_CSV"]

HYPOXIA_RULES_CSV = (
    Path(__file__).resolve().parents[2] / "references" / "physicell_hypoxia_cell_rules.csv"
)


def hypoxia_rules() -> List[Rule]:
    """The 8 verbatim hypoxia rules from the paper's reference implementation."""
    return load_rules_csv(HYPOXIA_RULES_CSV)


def hypoxia_cell_definitions() -> CellDefinitions:
    """Baseline phenotype values (b0 / bm) from the hypoxia cell definitions."""
    shared = {
        "cycle entry": {"base": 7.2e-4},
        "apoptosis": {"base": 5.31667e-5},
        "necrosis": {"base": 2.8e-3},
    }
    return CellDefinitions({
        "tumor": {
            **shared,
            "migration speed": {"base": 0.0},
            "migration persistence time": {"base": 0.0},
            "migration bias": {"base": 0.0},
            "transform to motile tumor": {"base": 0.001},
        },
        "motile tumor": {
            **shared,
            "migration speed": {"base": 0.47},
            "migration persistence time": {"base": 15.0},
            "migration bias": {"base": 0.18},
            "transform to tumor": {"base": 0.0},
        },
    })


def hypoxia_ruleset() -> RuleSet:
    return RuleSet(hypoxia_rules())
