"""Composite generators (dashboard-visible entry points)."""
from __future__ import annotations

from . import hypoxia  # noqa: F401  (side-effect: registers the generator)
from .hypoxia import hypoxia_motility

__all__ = ["hypoxia_motility"]
