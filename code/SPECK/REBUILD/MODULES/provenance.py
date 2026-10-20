"""Provenance tiers — every number in this engine carries one.

The tier vocabulary is spectre's (ORESTES_CLASS_REBUILD.txt §0.2), extended
by the KINETIC tier the plan's §5 architecture adds:

  KINETIC  first-principles export from PIC_MCC_3D, read as DATA (never an
           import); ranks above LOKI and is SCORED against it, not swapped in
           silently (S4).
  LOKI     electron-impact k(<eps>) from the curated LoKI-MC mixture tables
           (copied under MODULES/loki_tables/ with provenance).
  CHEN     Chen 2025 Table III verbatim (SI as printed; typography flags
           preserved in chen2025.py).
  REPO     value + citation lifted from phantom11/spectre config — copied,
           never re-tuned; the citation string names the source line.
  LIT      a value fetched first-hand from the literature by this repo
           (cite names paper, table, and fetch date).
  EST      no citation yet: stated analogy, flagged, and sensitivity-
           bracketed ×1/3..×3 by tests. S5 drives this count to zero.

Ordering is trust, highest first. ``TIER_ORDER.index(tier)`` sorts a
provenance report; an unknown tier is an error, not a default.
"""

from __future__ import annotations

from dataclasses import dataclass

TIER_ORDER = ("KINETIC", "LOKI", "CHEN", "REPO", "LIT", "EST")

# EST sensitivity bracket (multiplicative), the spectre convention.
EST_BRACKET = (1.0 / 3.0, 3.0)


def check_tier(tier: str) -> str:
    if tier not in TIER_ORDER:
        raise ValueError(f"unknown provenance tier {tier!r}; use {TIER_ORDER}")
    return tier


@dataclass(frozen=True)
class Sourced:
    """A value that knows where it came from."""

    value: float
    tier: str
    cite: str

    def __post_init__(self):
        check_tier(self.tier)
        if not self.cite:
            raise ValueError("a Sourced value with no citation is the defect "
                             "this class exists to prevent")


def tier_census(items) -> dict:
    """Count items per tier. ``items`` yields objects with a ``.tier``."""
    out = {t: 0 for t in TIER_ORDER}
    for it in items:
        out[check_tier(it.tier)] += 1
    return out


def census_line(census: dict) -> str:
    total = sum(census.values())
    parts = ", ".join(f"{t} {n}" for t, n in census.items() if n)
    return f"{total} rows: {parts}"
