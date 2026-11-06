"""Composite generator: the paper's hypoxia -> motility model on the reference substrate.

Wires the minimal reference substrate (oxygen field + point agents) to the
generic :class:`GrammarRulesProcess`, driven by the verbatim PhysiCell hypoxia
rules (``references/physicell_hypoxia_cell_rules.csv``) and the hypoxia cell
definitions.  Demonstrates the transient post-hypoxic motile phenotype end to
end — with the grammar engine as a *separate, pluggable* process, exactly as it
will attach to the MISM blueprint's ABMs.

Oxygen schedule: a dense tumor makes a hypoxic core via uptake (necrotic core +
transition to motile at the rim); late, cellular uptake is relieved so the core
reoxygenates and motile cells revert — the transience the paper reports.
"""
from __future__ import annotations

from process_bigraph.composite_generator import composite_generator

from ..io import rules_to_csv
from ..models.hypoxia import hypoxia_cell_definitions, hypoxia_rules

_EMIT = ["stores/populations", "stores/radial", "stores/time"]


def _cell_defs_dict():
    return hypoxia_cell_definitions().to_dict()


@composite_generator(
    name="hypoxia_motility",
    description="Paper's oxygen->motility tumor model: grammar rules on a minimal "
                "oxygen-field + point-agent reference substrate (transient post-hypoxic migration).",
    parameters={
        "n_cells": {"type": "integer", "default": 350,
                    "description": "Initial tumor cells seeded in a disk."},
        "seed_radius": {"type": "float", "default": 130.0,
                        "description": "Radius (µm) of the initial tumor disk."},
        "interval": {"type": "float", "default": 30.0,
                     "description": "Composite tick (min); rules use per-time rates."},
        "uptake": {"type": "float", "default": 1.0,
                   "description": "Per-cell oxygen uptake (controls core hypoxia / gradient)."},
        "reoxygenation_time": {"type": "float", "default": 3000.0,
                               "description": "Minute at which uptake is relieved (reoxygenation)."},
        "seed": {"type": "integer", "default": 0, "description": "RNG seed."},
    },
    emitters=[{"address": "local:ParquetEmitter", "paths": _EMIT}],
)
def hypoxia_motility(core=None, *, n_cells=350, seed_radius=130.0, interval=30.0,
                     uptake=50.0, reoxygenation_time=3000.0, seed=0):
    cell_defs = _cell_defs_dict()
    rules_csv = rules_to_csv(hypoxia_rules())

    # Reoxygenation: full uptake until reoxygenation_time, then ramp uptake -> 0.
    uptake_schedule = [
        [0.0, 1.0],
        [float(reoxygenation_time), 1.0],
        [float(reoxygenation_time) + 600.0, 0.05],
        [float(reoxygenation_time) + 3000.0, 0.05],
    ]

    substrate = {
        "_type": "process",
        "address": "local:ReferenceSubstrateProcess",
        "config": {
            "domain": {"x_min": -400.0, "x_max": 400.0, "y_min": -400.0, "y_max": 400.0},
            "dx": 20.0,
            "cell_radius": 8.4,
            "max_cells": 1600,
            # D is reduced from PhysiCell's 1e5 µm²/min to give this small toy
            # substrate a ~130 µm diffusion length, so a hypoxic core forms at
            # toy scale. The grammar rules/parameters are verbatim; only the
            # substrate's PDE constants are toy-scaled (see README).
            "oxygen": {"D": 5.0e3, "decay": 0.1, "boundary": 38.0,
                       "uptake": float(uptake), "relax_iters": 40,
                       "uptake_schedule": uptake_schedule},
            "cell_definitions": cell_defs,
            "initial_cells": {"n": int(n_cells), "radius": float(seed_radius), "type": "tumor"},
            "seed": int(seed),
        },
        "interval": float(interval),
        "inputs": {
            "fates": ["stores", "fates"],
            "divisions": ["stores", "divisions"],
            "deaths": ["stores", "deaths"],
            "scalar": ["stores", "scalar"],
        },
        "outputs": {
            "positions": ["stores", "positions"],
            "sizes": ["stores", "sizes"],
            "types": ["stores", "types"],
            "shared_surface_by_type": ["stores", "shared_surface_by_type"],
            "neighbors": ["stores", "neighbors"],
            "alive": ["stores", "alive"],
            "field_grid": ["stores", "field_grid"],
            "pressure": ["stores", "pressure"],
            "populations": ["stores", "populations"],
            "radial": ["stores", "radial"],
            "time": ["stores", "time"],
        },
    }

    grammar = {
        "_type": "process",
        "address": "local:GrammarRulesProcess",
        "config": {
            "rules_csv": rules_csv,
            "cell_definitions": cell_defs,
            "seed": int(seed) + 1,
            "rate_mode": "per_time",
            # map world coords -> oxygen grid indices (same as the substrate's grid)
            "field_sampling": {"origin": [-400.0, -400.0], "spacing": 20.0},
        },
        "interval": float(interval),
        "inputs": {
            "positions": ["stores", "positions"],
            "sizes": ["stores", "sizes"],
            "types": ["stores", "types"],
            "shared_surface_by_type": ["stores", "shared_surface_by_type"],
            "neighbors": ["stores", "neighbors"],
            "alive": ["stores", "alive"],
            "field_grid": ["stores", "field_grid"],
            "pressure": ["stores", "pressure"],
            "time": ["stores", "time"],
        },
        "outputs": {
            "fates": ["stores", "fates"],
            "scalar": ["stores", "scalar"],
            "divisions": ["stores", "divisions"],
            "deaths": ["stores", "deaths"],
        },
    }

    return {
        "substrate": substrate,
        "grammar": grammar,
        "stores": {
            "fates": {}, "scalar": {}, "divisions": [], "deaths": [],
            "positions": {}, "sizes": {}, "types": {}, "shared_surface_by_type": {},
            "neighbors": {}, "alive": {}, "field_grid": {}, "pressure": {},
            "populations": {}, "radial": {}, "time": 0.0,
        },
        "emitter": {
            "_type": "step",
            "address": "local:RAMEmitter",
            "config": {"emit": {"populations": "map[float]", "radial": "map[float]", "time": "float"}},
            "inputs": {
                "populations": ["stores", "populations"],
                "radial": ["stores", "radial"],
                "time": ["stores", "time"],
            },
        },
    }
