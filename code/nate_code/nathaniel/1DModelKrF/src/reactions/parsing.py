"""Parse reactant/product operands into (multiplicity, species-name) terms.

Accepts either a tuple of tokens (``('Kr+', '2e')``) or a single string
(``'Kr+ + 2e'`` / ``'Kr+ + 2e-'``). A token is an optional integer coefficient
followed by a species name; species names keep their charge/excitation markers
(``Kr+``, ``Kr*``, ``F-``). The reaction-level ``+`` separator must be surrounded by
whitespace (`` + ``) so the charge ``+`` in ``Kr+`` is never mistaken for a separator.
Common source-notation aliases are normalized to registry names.
"""

from __future__ import annotations

import re
from collections.abc import Iterable

# Source notation -> registry species name. Bare ``hv`` defaults to the 248 nm KrF
# line; prefer the explicit ``hv248``/``hv193``/``hv95`` in data files.
ALIASES: dict[str, str] = {
    "e": "e-",
    "e-": "e-",
    "hv": "hv248",
    "hnu": "hv248",
}

_SEP = re.compile(r"\s+\+\s+")  # the reaction '+' separator (whitespace on both sides)
_COEFF = re.compile(r"^\s*(\d+)\s*(.+?)\s*$")  # leading integer coefficient + name


def normalize_name(name: str) -> str:
    """Map a source token to its registry species name (identity if unknown)."""
    return ALIASES.get(name, name)


def parse_side(operand: str | Iterable[str]) -> tuple[tuple[int, str], ...]:
    """Return ``((count, species_name), ...)`` for one reaction side.

    ``'Kr+ + 2e'`` -> ``((1, 'Kr+'), (2, 'e-'))``; ``('Kr+', '2e')`` -> same.
    Repeated species are kept as separate terms (the caller/compiler sums orders).
    """
    if isinstance(operand, str):
        tokens = [t for t in _SEP.split(operand.strip()) if t]
    else:
        tokens = [str(t).strip() for t in operand if str(t).strip()]

    terms: list[tuple[int, str]] = []
    for token in tokens:
        m = _COEFF.match(token)
        if m:
            count, name = int(m.group(1)), m.group(2).strip()
        else:
            count, name = 1, token
        if count < 1:
            raise ValueError(f"non-positive multiplicity in reaction token {token!r}")
        terms.append((count, normalize_name(name)))
    return tuple(terms)
