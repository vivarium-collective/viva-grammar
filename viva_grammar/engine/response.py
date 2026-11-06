"""Response-function mathematics of the cell-behavior hypothesis grammar.

Clean-room reimplementation of the response math described in:

    Johnson, Bergman, …, Fertig, Macklin (2025). "Human interpretable grammar
    encodes multicellular systems biology models to democratize virtual cell
    laboratories." *Cell* 188, 4711-4733.  https://doi.org/10.1016/j.cell.2025.06.048

This module is **not** PhysiCell; it re-derives the published equations
(eqs 1-4, Methods S1) in portable NumPy so the grammar can drive any
agent-based model.  See ``references/PROVENANCE.md``.

Equations
---------
Single rule (eq 1):      ``b(s) = b0 + (bM - b0) * R(s)``
Multivariate up (eq 2):  ``U = S_u / (1 + S_u)``,  ``S_u = Σ_i (u_i / u*_i)**p_i``
Multivariate down (eq 3):``D = S_d / (1 + S_d)``,  ``S_d = Σ_i (d_i / d*_i)**q_i``
Combined (eq 4):         ``b = (1 - D) * [(1 - U) * b0 + U * bM] + D * bm``

``U`` and ``D`` each vary in ``[0, 1)``.  With a single up signal and no down
signals, eq 4 collapses to the familiar single Hill curve ``b(u)``.
"""
from __future__ import annotations

from typing import Iterable, Sequence

__all__ = [
    "hill_term",
    "hill_response",
    "linear_response",
    "step_response",
    "saturating_response",
    "combine",
]


def hill_term(value: float, half_max: float, hill_power: float) -> float:
    """A single normalized Hill term ``(value / half_max)**hill_power``.

    These terms are summed to form the multivariate up/down responses
    (``S_u`` / ``S_d`` in eqs 2-3).  Negative ``value`` is clamped to 0;
    a non-positive ``half_max`` is treated as "no response" (term 0).
    """
    v = max(float(value), 0.0)
    hm = float(half_max)
    if hm <= 0.0:
        # A zero half-max would divide by zero; interpret as an inactive signal.
        return 0.0
    return (v / hm) ** float(hill_power)


def hill_response(value: float, half_max: float, hill_power: float) -> float:
    """Single-signal Hill response ``t / (1 + t)`` with ``t = (value/half_max)**n``.

    Varies monotonically in ``[0, 1)``; equals 0.5 at ``value == half_max``.
    """
    t = hill_term(value, half_max, hill_power)
    return t / (1.0 + t)


def linear_response(value: float, half_max: float) -> float:
    """Capped-linear response in ``[0, 1]`` (Fig S2 alternative to Hill).

    Reaches 0.5 at ``value == half_max`` and saturates at 1 for
    ``value >= 2 * half_max``.
    """
    hm = float(half_max)
    if hm <= 0.0:
        return 0.0
    return min(max(float(value), 0.0) / (2.0 * hm), 1.0)


def step_response(value: float, half_max: float) -> float:
    """Heaviside step response (Fig S2 alternative): 1 if ``value >= half_max`` else 0."""
    return 1.0 if float(value) >= float(half_max) else 0.0


def saturating_response(terms: Iterable[float]) -> float:
    """Multivariate response ``S / (1 + S)`` for ``S = Σ terms`` (eqs 2-3).

    ``terms`` are the per-signal Hill terms from :func:`hill_term`.  An empty
    iterable yields 0 (no signals of this sign -> no response).
    """
    s = 0.0
    for t in terms:
        s += max(float(t), 0.0)
    return s / (1.0 + s)


def combine(b0: float, b_max: float, b_min: float, up: float, down: float) -> float:
    """Bilinear combine of up/down responses into a behavior value (eq 4).

    Parameters
    ----------
    b0 : base value (no signal).
    b_max : saturated value under up-regulating signals.
    b_min : saturated value under down-regulating signals.
    up, down : the ``U`` and ``D`` multivariate responses, each in ``[0, 1]``.
    """
    b0 = float(b0)
    target = (1.0 - up) * b0 + up * float(b_max)
    return (1.0 - down) * target + down * float(b_min)
