"""Integration tests: the hypoxia composite runs, the generator is registered,
and the paper's qualitative phenotype (induction, plume, reversion, core) emerges.
"""
from __future__ import annotations

import pytest
from process_bigraph import Composite, gather_emitter_results

from viva_grammar.composites.hypoxia import hypoxia_motility
from viva_grammar.core import build_core


def test_generator_is_registered():
    from process_bigraph.composite_generator import _REGISTRY
    matches = [eid for eid in _REGISTRY if eid.endswith(".hypoxia_motility")]
    assert matches, f"hypoxia_motility missing; have {list(_REGISTRY)[:5]}"


def test_composite_builds_and_runs_one_step():
    core = build_core()
    doc = hypoxia_motility(core, n_cells=80, seed_radius=90.0, interval=30.0)
    sim = Composite({"state": doc}, core=core)
    sim.run(30.0)
    rows = gather_emitter_results(sim)[("emitter",)]
    assert rows and "populations" in rows[-1]


@pytest.mark.slow
def test_hypoxia_phenotype_hits_acceptance_bands():
    """Full qualitative reproduction of the paper's hypoxia -> motility dynamics."""
    core = build_core()
    interval = 30.0
    doc = hypoxia_motility(core, n_cells=350, seed_radius=130.0, interval=interval,
                           uptake=1.0, reoxygenation_time=3000.0, seed=0)
    sim = Composite({"state": doc}, core=core)
    sim.run(150 * interval)
    rows = gather_emitter_results(sim)[("emitter",)]

    def frac(r):
        p = r["populations"]
        live = p.get("live_tumor", 0) + p.get("live_motile_tumor", 0)
        return (p.get("live_motile_tumor", 0) / live) if live else 0.0

    peak_i = max(range(len(rows)), key=lambda i: frac(rows[i]))
    peak_frac = frac(rows[peak_i])
    peak_row = rows[peak_i]

    # Band 1: induction — peak motile fraction >= 0.15
    assert peak_frac >= 0.15, f"peak motile fraction {peak_frac:.3f} < 0.15"

    # Band 2: invasive plume — motile cells further out than non-motile at the peak
    rad = peak_row["radial"]
    assert rad["mean_motile_radius"] >= 1.2 * rad["mean_nonmotile_radius"], \
        f"plume: motile r {rad['mean_motile_radius']:.1f} vs nonmotile {rad['mean_nonmotile_radius']:.1f}"

    # Band 3: transience — reverts to <= 0.02 after reoxygenation (final frame)
    assert frac(rows[-1]) <= 0.02, f"final motile fraction {frac(rows[-1]):.3f} did not revert"

    # Band 4: necrotic core central — necrotic radius < live radius once a core exists
    core_rows = [r for r in rows if r["populations"].get("necrotic", 0) > 10]
    assert core_rows, "no necrotic core formed"
    for r in core_rows:
        rd = r["radial"]
        if rd["mean_nonmotile_radius"] > 0:
            assert rd["mean_necrotic_radius"] < rd["mean_nonmotile_radius"]
