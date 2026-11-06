"""A deliberately-minimal 2D agent-based substrate for standalone grammar studies.

This is **not** a serious ABM — it exists only so the grammar engine has
*something* to drive inside ``viva-grammar`` with no external simulator, and so
the engine's ABM-agnosticism is demonstrated (this reference substrate is ABM
#1; the MISM blueprint's CPM/vertex/munk backends are ABM #2+).

It holds the authoritative cell cohort + oxygen field on the instance, and each
step:
  1. applies the grammar's events read on its input ports — ``divisions`` (spawn
     a daughter), ``deaths`` (apoptosis removes; necrosis leaves a persistent
     necrotic corpse), ``fates`` (change cell type), ``scalar`` (per-cell
     parameters, e.g. an overriding migration speed);
  2. advances a quasi-steady oxygen field (Jacobi relaxation of
     ``D∇²O - decay·O - uptake·ρ_live·O = 0`` with a time-scheduled Dirichlet
     boundary) — cheap and stable despite oxygen's fast diffusion;
  3. moves motile cells by a persistent walk biased up the oxygen gradient, at
     the cell type's baseline migration speed (from the cell definitions) unless
     the grammar supplied one;
  4. emits a fresh cohort snapshot on the readout ports the grammar consumes,
     plus populations / radial metrics for the study.

Migration speed is a *baseline property of the cell type* here (tumor 0.0,
motile tumor 0.47 µm/min), matching the paper's transition-based hypoxia model.
"""
from __future__ import annotations

import math
from collections import defaultdict
from typing import Dict, List, Optional

import numpy as np
from process_bigraph import Process

from ..engine import CellDefinitions

__all__ = ["ReferenceSubstrateProcess"]

LIVE, NECROTIC, APOPTOTIC = "live", "necrotic", "apoptotic"


def _schedule_value(schedule, t: float, default: float) -> float:
    """Piecewise-linear interpolation of a [(time, value), ...] schedule."""
    if not schedule:
        return default
    pts = sorted((float(a), float(b)) for a, b in schedule)
    if t <= pts[0][0]:
        return pts[0][1]
    if t >= pts[-1][0]:
        return pts[-1][1]
    for (t0, v0), (t1, v1) in zip(pts, pts[1:]):
        if t0 <= t <= t1:
            f = 0.0 if t1 == t0 else (t - t0) / (t1 - t0)
            return v0 + f * (v1 - v0)
    return default


class ReferenceSubstrateProcess(Process):
    config_schema = {
        "domain": {"_type": "node", "_default": {"x_min": -400.0, "x_max": 400.0,
                                                 "y_min": -400.0, "y_max": 400.0}},
        "dx": {"_type": "float", "_default": 20.0},
        "cell_radius": {"_type": "float", "_default": 8.4},
        "max_cells": {"_type": "integer", "_default": 1500},
        "oxygen": {"_type": "node", "_default": {}},  # D, decay, boundary, uptake, relax_iters, schedule
        "cell_definitions": {"_type": "node", "_default": {}},
        "initial_cells": {"_type": "node", "_default": {}},  # {n, radius, type} OR {"cells":[...]}
        "seed": {"_type": "integer", "_default": 0},
    }

    # -- lifecycle -------------------------------------------------------------
    def initialize(self, config):
        self._rng = np.random.default_rng(int(config.get("seed", 0)))
        d = config["domain"]
        self._x0, self._x1 = float(d["x_min"]), float(d["x_max"])
        self._y0, self._y1 = float(d["y_min"]), float(d["y_max"])
        self._dx = float(config["dx"])
        self._nx = max(int(round((self._x1 - self._x0) / self._dx)), 3)
        self._ny = max(int(round((self._y1 - self._y0) / self._dx)), 3)
        self._radius = float(config["cell_radius"])
        self._max_cells = int(config["max_cells"])
        self._cdefs = CellDefinitions(config.get("cell_definitions") or {})

        ox = config.get("oxygen") or {}
        self._oD = float(ox.get("D", 1.0e5))
        self._odecay = float(ox.get("decay", 0.1))
        self._oboundary = float(ox.get("boundary", 38.0))
        self._ouptake = float(ox.get("uptake", 10.0))
        self._orelax = int(ox.get("relax_iters", 25))
        self._oschedule = ox.get("schedule") or []  # [(time, boundary_value)]
        # multiplier on cellular uptake over time; dropping it late reoxygenates
        # the core (the controllable "reoxygenation" lever for the study).
        self._ouptake_schedule = ox.get("uptake_schedule") or []
        self._uptake_mult = 1.0

        self._field = np.full((self._ny, self._nx), self._oboundary, dtype=float)
        self._time = 0.0
        self._next_id = 0
        self._cells: Dict[str, dict] = {}
        self._heading: Dict[str, np.ndarray] = {}
        self._seed_cells(config.get("initial_cells") or {})
        self._relax_oxygen(self._oboundary)

    def _new_id(self) -> str:
        cid = str(self._next_id)
        self._next_id += 1
        return cid

    def _seed_cells(self, spec: dict):
        if spec.get("cells"):
            for c in spec["cells"]:
                self._add_cell(c["type"], float(c["x"]), float(c["y"]))
            return
        n = int(spec.get("n", 300))
        radius = float(spec.get("radius", 150.0))
        ctype = spec.get("type", "tumor")
        for _ in range(n):
            # uniform in a disk
            r = radius * math.sqrt(self._rng.random())
            th = 2.0 * math.pi * self._rng.random()
            self._add_cell(ctype, r * math.cos(th), r * math.sin(th))

    def _add_cell(self, ctype: str, x: float, y: float) -> str:
        cid = self._new_id()
        self._cells[cid] = {"type": ctype, "pos": np.array([x, y], float), "state": LIVE}
        self._heading[cid] = self._rng.normal(size=2)
        self._heading[cid] /= (np.linalg.norm(self._heading[cid]) + 1e-12)
        return cid

    # -- ports -----------------------------------------------------------------
    def inputs(self):
        return {
            "fates": "map[string]",
            "divisions": "list[node]",
            "deaths": "list[node]",
            "scalar": "map[map[float]]",
        }

    def outputs(self):
        return {
            "positions": "overwrite[map[list[float]]]",
            "sizes": "overwrite[map[float]]",
            "types": "overwrite[map[string]]",
            "shared_surface_by_type": "overwrite[map[map[float]]]",
            "neighbors": "overwrite[map[list[string]]]",
            "alive": "overwrite[map[boolean]]",
            "field_grid": "overwrite[map[node]]",
            "pressure": "overwrite[map[float]]",
            "populations": "overwrite[map[float]]",
            "radial": "overwrite[map[float]]",
            "time": "overwrite[float]",
        }

    def initial_state(self):
        return self._snapshot()

    # -- grid helpers ----------------------------------------------------------
    def _site(self, pos) -> tuple:
        ix = int(min(max(round((pos[0] - self._x0) / self._dx), 0), self._nx - 1))
        iy = int(min(max(round((pos[1] - self._y0) / self._dx), 0), self._ny - 1))
        return iy, ix  # row, col

    def _live_density(self) -> np.ndarray:
        rho = np.zeros((self._ny, self._nx), dtype=float)
        for c in self._cells.values():
            if c["state"] == LIVE:
                iy, ix = self._site(c["pos"])
                rho[iy, ix] += 1.0
        return rho

    def _relax_oxygen(self, boundary: float):
        """A few Jacobi sweeps toward the steady diffusion-reaction state."""
        rho = self._live_density()
        O = self._field
        O[0, :] = O[-1, :] = O[:, 0] = O[:, -1] = boundary  # Dirichlet
        coef_diff = self._oD / (self._dx * self._dx)
        # steady: coef*(sum_neighbors - 4*O) - (decay + uptake*rho)*O = 0
        denom = 4.0 * coef_diff + self._odecay + self._ouptake * self._uptake_mult * rho
        for _ in range(self._orelax):
            nbr = np.zeros_like(O)
            nbr[1:-1, 1:-1] = O[:-2, 1:-1] + O[2:, 1:-1] + O[1:-1, :-2] + O[1:-1, 2:]
            newO = coef_diff * nbr / denom
            newO[0, :] = newO[-1, :] = newO[:, 0] = newO[:, -1] = boundary
            O[1:-1, 1:-1] = newO[1:-1, 1:-1]
        self._field = O

    def _oxygen_grad(self, pos) -> np.ndarray:
        iy, ix = self._site(pos)
        gx = (self._field[iy, min(ix + 1, self._nx - 1)] - self._field[iy, max(ix - 1, 0)])
        gy = (self._field[min(iy + 1, self._ny - 1), ix] - self._field[max(iy - 1, 0), ix])
        return np.array([gx, gy], float)

    # -- neighbor search (spatial binning) ------------------------------------
    def _neighbors(self):
        contact = 2.0 * self._radius
        bin_size = contact
        bins: Dict[tuple, List[str]] = defaultdict(list)
        for cid, c in self._cells.items():
            b = (int(c["pos"][0] // bin_size), int(c["pos"][1] // bin_size))
            bins[b].append(cid)
        neighbors: Dict[str, List[str]] = {cid: [] for cid in self._cells}
        shared: Dict[str, Dict[str, float]] = {cid: defaultdict(float) for cid in self._cells}
        for cid, c in self._cells.items():
            bx, by = int(c["pos"][0] // bin_size), int(c["pos"][1] // bin_size)
            for dxb in (-1, 0, 1):
                for dyb in (-1, 0, 1):
                    for other in bins.get((bx + dxb, by + dyb), ()):
                        if other == cid:
                            continue
                        dist = float(np.linalg.norm(c["pos"] - self._cells[other]["pos"]))
                        if dist < contact:
                            neighbors[cid].append(other)
                            overlap = max(contact - dist, 0.0)
                            shared[cid][self._cells[other]["type"]] += overlap
        return neighbors, {k: dict(v) for k, v in shared.items()}

    # -- snapshot --------------------------------------------------------------
    def _snapshot(self) -> dict:
        neighbors, shared = self._neighbors()
        positions = {cid: [float(c["pos"][0]), float(c["pos"][1])] for cid, c in self._cells.items()}
        types = {cid: c["type"] for cid, c in self._cells.items()}
        sizes = {cid: math.pi * self._radius ** 2 for cid in self._cells}
        alive = {cid: (c["state"] == LIVE) for cid, c in self._cells.items()}
        pressure = {cid: float(len(neighbors[cid])) for cid in self._cells}

        pops: Dict[str, float] = defaultdict(float)
        for c in self._cells.values():
            if c["state"] == LIVE:
                pops[f"live_{c['type'].replace(' ', '_')}"] += 1.0
            elif c["state"] == NECROTIC:
                pops["necrotic"] += 1.0
        pops["total_live"] = sum(v for k, v in pops.items() if k.startswith("live_"))
        pops["min_oxygen"] = float(self._field.min())
        pops["mean_oxygen"] = float(self._field.mean())

        radial = self._radial_metrics()
        return {
            "positions": positions, "sizes": sizes, "types": types,
            "shared_surface_by_type": shared, "neighbors": neighbors,
            "alive": alive, "field_grid": {"oxygen": self._field.tolist()},
            "pressure": pressure, "populations": dict(pops), "radial": radial,
            "time": self._time,
        }

    def _radial_metrics(self) -> Dict[str, float]:
        live_r, motile_r, necro_r = [], [], []
        for c in self._cells.values():
            r = float(np.linalg.norm(c["pos"]))
            if c["state"] == NECROTIC:
                necro_r.append(r)
            elif c["state"] == LIVE:
                live_r.append(r)
                if c["type"] == "motile tumor":
                    motile_r.append(r)
        mean = lambda xs: float(np.mean(xs)) if xs else 0.0
        nonmotile_r = [r for c, r in
                       ((cc, float(np.linalg.norm(cc["pos"]))) for cc in self._cells.values())
                       if c["state"] == LIVE and c["type"] != "motile tumor"]
        return {
            "mean_live_radius": mean(live_r),
            "mean_motile_radius": mean(motile_r),
            "mean_nonmotile_radius": mean(nonmotile_r),
            "mean_necrotic_radius": mean(necro_r),
            "n_motile": float(len(motile_r)),
        }

    # -- main step -------------------------------------------------------------
    def update(self, inputs, interval=1.0):
        interval = float(interval)
        # 1) apply grammar events ------------------------------------------------
        for cid, new_type in (inputs.get("fates") or {}).items():
            if cid in self._cells and self._cells[cid]["state"] == LIVE:
                self._cells[cid]["type"] = new_type
        for ev in (inputs.get("deaths") or []):
            cid = ev["cell_id"] if isinstance(ev, dict) else ev
            cause = ev.get("cause", "apoptosis") if isinstance(ev, dict) else "apoptosis"
            if cid in self._cells:
                if cause == "necrosis":
                    self._cells[cid]["state"] = NECROTIC  # persists -> forms the core
                else:
                    self._cells.pop(cid, None); self._heading.pop(cid, None)
        for ev in (inputs.get("divisions") or []):
            cid = ev["cell_id"] if isinstance(ev, dict) else ev
            if cid in self._cells and self._cells[cid]["state"] == LIVE and len(self._cells) < self._max_cells:
                parent = self._cells[cid]
                off = self._rng.normal(scale=self._radius, size=2)
                self._add_cell(parent["type"], parent["pos"][0] + off[0], parent["pos"][1] + off[1])
        scalar = inputs.get("scalar") or {}

        # 2) advance oxygen (scheduled Dirichlet boundary + uptake multiplier) ---
        self._time += interval
        boundary = _schedule_value(self._oschedule, self._time, self._oboundary)
        self._uptake_mult = _schedule_value(self._ouptake_schedule, self._time, 1.0)
        self._relax_oxygen(boundary)

        # 3) move motile cells ---------------------------------------------------
        for cid, c in self._cells.items():
            if c["state"] != LIVE:
                continue
            speed = float((scalar.get(cid) or {}).get("migration speed",
                                                       self._cdefs.base(c["type"], "migration speed")))
            if speed <= 0.0:
                continue
            bias = self._cdefs.base(c["type"], "migration bias") or 0.18
            g = self._oxygen_grad(c["pos"])
            gn = np.linalg.norm(g)
            grad_dir = g / gn if gn > 1e-9 else np.zeros(2)
            rand = self._rng.normal(size=2)
            rand /= (np.linalg.norm(rand) + 1e-12)
            h = bias * grad_dir + (1.0 - bias) * rand + 0.4 * self._heading[cid]
            hn = np.linalg.norm(h)
            if hn > 1e-9:
                h = h / hn
            self._heading[cid] = h
            new = c["pos"] + speed * interval * h
            new[0] = min(max(new[0], self._x0), self._x1)
            new[1] = min(max(new[1], self._y0), self._y1)
            c["pos"] = new

        return self._snapshot()
