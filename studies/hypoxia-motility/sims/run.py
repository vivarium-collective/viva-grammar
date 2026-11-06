#!/usr/bin/env python3
"""Canonical run for the ``hypoxia-motility`` study.

Drives the ``hypoxia_motility`` composite (reference substrate + grammar engine),
captures the spatial + population timecourse, scores the four acceptance bands,
renders the interactive visualization, and records the run in ``.pbg/runs.jsonl``.
"""
from __future__ import annotations

import json
import time
import uuid
from pathlib import Path

STUDY_DIR = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = STUDY_DIR.parents[1]

import sys
sys.path.insert(0, str(WORKSPACE_ROOT))

from process_bigraph import Composite, gather_emitter_results  # noqa: E402

from viva_grammar.core import build_core  # noqa: E402
from viva_grammar.composites.hypoxia import hypoxia_motility  # noqa: E402
from viva_grammar.viz import render_hypoxia_html  # noqa: E402

SPEC_ID = "viva_grammar.composites.hypoxia.hypoxia_motility"
STUDY_SLUG = "hypoxia-motility"
INVESTIGATION_SLUG = "grammar-showcase"

INTERVAL = 30.0
N_TICKS = 150
REOXY_TIME = 3000.0
CAPTURE_EVERY = 4
DOMAIN = {"x_min": -400.0, "x_max": 400.0, "y_min": -400.0, "y_max": 400.0}


def _append_run_event(event: dict) -> None:
    """Record a run event, preferring vivarium_workbench, else a local fallback."""
    try:
        from vivarium_workbench.lib.run_log import append_run_event
        append_run_event(WORKSPACE_ROOT, event)
        return
    except Exception:
        pbg = WORKSPACE_ROOT / ".pbg"
        pbg.mkdir(exist_ok=True)
        with (pbg / "runs.jsonl").open("a") as fh:
            fh.write(json.dumps(event) + "\n")


def _frac(pops):
    live = pops.get("live_tumor", 0.0) + pops.get("live_motile_tumor", 0.0)
    return (pops.get("live_motile_tumor", 0.0) / live) if live else 0.0


def _build_frames(rows):
    frames = []
    for i, r in enumerate(rows):
        if i % CAPTURE_EVERY and i != len(rows) - 1:
            continue
        positions = r.get("positions", {})
        types = r.get("types", {})
        alive = r.get("alive", {})
        xy, state = [], []
        for cid, p in positions.items():
            xy.append([float(p[0]), float(p[1])])
            state.append(types.get(cid, "tumor") if alive.get(cid, True) else "necrotic")
        frames.append({
            "time": r.get("time", i * INTERVAL),
            "oxygen": r["field_grid"]["oxygen"],
            "cell_xy": xy, "cell_state": state,
            "populations": r["populations"], "radial": r["radial"],
        })
    return frames


def _score_bands(rows):
    fr = [_frac(r["populations"]) for r in rows]
    peak_i = max(range(len(rows)), key=lambda i: fr[i])
    peak = rows[peak_i]
    rad = peak["radial"]
    plume_ratio = (rad["mean_motile_radius"] / rad["mean_nonmotile_radius"]
                   if rad["mean_nonmotile_radius"] else 0.0)
    final = fr[-1]
    core_rows = [r for r in rows if r["populations"].get("necrotic", 0) > 10]
    core_ok = all(r["radial"]["mean_necrotic_radius"] < r["radial"]["mean_nonmotile_radius"]
                  for r in core_rows if r["radial"]["mean_nonmotile_radius"] > 0) and bool(core_rows)
    cards = [
        {"label": "Phenotype induction", "value": f"{fr[peak_i]:.2f}",
         "ok": fr[peak_i] >= 0.15, "detail": "peak motile fraction ≥ 0.15"},
        {"label": "Invasive plume", "value": f"{plume_ratio:.2f}×",
         "ok": plume_ratio >= 1.2, "detail": "motile radius ≥ 1.2× non-motile (at peak)"},
        {"label": "Transience / reversion", "value": f"{final:.3f}",
         "ok": final <= 0.02, "detail": "motile fraction ≤ 0.02 after reoxygenation"},
        {"label": "Necrotic core", "value": "central" if core_ok else "—",
         "ok": core_ok, "detail": "necrotic radius < live radius throughout"},
    ]
    return {"cards": cards, "all_pass": all(c["ok"] for c in cards)}


def main() -> int:
    run_id = uuid.uuid4().hex
    _append_run_event({
        "run_id": run_id, "event": "started", "spec_id": SPEC_ID, "label": STUDY_SLUG,
        "started_at": time.time(), "status": "running", "n_steps": N_TICKS, "emitter": "ram",
        "origin": "canonical_run", "study_slug": STUDY_SLUG,
        "investigation_slug": INVESTIGATION_SLUG, "params": {"interval": INTERVAL, "n_ticks": N_TICKS},
    })
    try:
        core = build_core()
        doc = hypoxia_motility(core, n_cells=350, seed_radius=130.0, interval=INTERVAL,
                               uptake=1.0, reoxygenation_time=REOXY_TIME, seed=0)
        # augment the emitter to also capture spatial state for the animation
        doc["emitter"]["config"]["emit"].update({
            "positions": "map[list[float]]", "types": "map[string]",
            "alive": "map[boolean]", "field_grid": "map[node]",
        })
        doc["emitter"]["inputs"].update({
            "positions": ["stores", "positions"], "types": ["stores", "types"],
            "alive": ["stores", "alive"], "field_grid": ["stores", "field_grid"],
        })
        sim = Composite({"state": doc}, core=core)
        sim.run(N_TICKS * INTERVAL)
        rows = gather_emitter_results(sim)[("emitter",)]

        bands = _score_bands(rows)
        frames = _build_frames(rows)
        html = render_hypoxia_html(frames, DOMAIN, REOXY_TIME, bands)
        viz_dir = STUDY_DIR / "viz"
        viz_dir.mkdir(parents=True, exist_ok=True)
        (viz_dir / "hypoxia_motility.html").write_text(html)

        peak = max(_frac(r["populations"]) for r in rows)
        print(f"captured {len(rows)} ticks, {len(frames)} viz frames; peak motile fraction {peak:.3f}")
        for c in bands["cards"]:
            print(f"  [{'PASS' if c['ok'] else 'FAIL'}] {c['label']}: {c['value']} ({c['detail']})")
        print(f"viz -> {viz_dir / 'hypoxia_motility.html'}")
    except Exception:
        _append_run_event({"run_id": run_id, "event": "completed", "completed_at": time.time(),
                           "n_steps": 0, "status": "failed"})
        raise

    _append_run_event({"run_id": run_id, "event": "completed", "completed_at": time.time(),
                       "n_steps": N_TICKS, "status": "completed",
                       "outcomes": {"all_bands_pass": bands["all_pass"]}})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
