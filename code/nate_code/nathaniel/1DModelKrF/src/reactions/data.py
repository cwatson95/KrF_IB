"""Load a reaction network from a TOML data file (stdlib ``tomllib``).

Each ``[[reaction]]`` row is dispatched to :meth:`ReactionNetwork.add` by its
``rate_form``. Tagged rate constants may be published in cgs or SI and are converted
to internal SI at load. See ``data/reactions/krf.toml`` for the schema.
"""

from __future__ import annotations

import tomllib
from copy import deepcopy
from pathlib import Path

from particles.species import SpeciesRegistry

from .network import ReactionNetwork
from .types import RateForm


def _add_row(net: ReactionNetwork, row: dict) -> None:
    form = RateForm(row["rate_form"])
    common = dict(
        id=str(row.get("id", "")),
        reference=str(row.get("reference", "")),
        enabled=bool(row.get("enabled", True)),
        photon_mode=row.get("photon_mode"),
    )
    reactants, products = row["reactants"], row["products"]

    if form is RateForm.CONSTANT:
        net.add(reactants, products, r=row["k"], r_units=row["k_units"], **common)
    elif form is RateForm.RADIATIVE:
        net.add(reactants, products, tau_s=float(row["tau_s"]), **common)
    elif form is RateForm.TE_DEPENDENT:
        # E0_eV controls the published analytic rate expression; energy_loss_eV is the
        # independently declared electron-energy transfer.  Keep both.  Treating E0 as the
        # loss (or dropping the explicit loss here) gives the chemistry and transport closure
        # different reaction-energy tables.
        params = {
            k: row[k]
            for k in ("E0_eV", "te_power", "energy_loss_eV", "energy_loss_ev")
            if k in row
        }
        net.add(
            reactants, products, r=row["K0"], r_units=row.get("k_units", "cm3/s"),
            rate_form=RateForm.TE_DEPENDENT, params=params, **common,
        )
    elif form is RateForm.ELECTRON_IMPACT:
        params = {
            k: row[k]
            for k in ("threshold_eV", "threshold_ev", "energy_loss_eV", "energy_loss_ev")
            if k in row
        }
        net.add(
            reactants, products, sigma=row.get("cross_section"),
            rate_form=RateForm.ELECTRON_IMPACT, params=params, **common,
        )
    elif form is RateForm.PHOTO:
        if "sigma_m2" not in row:
            raise ValueError(
                f"photo reaction {common['id'] or reactants!r} requires sigma_m2"
            )
        sig = float(row["sigma_m2"])
        params = {"sigma_m2": sig}
        # Photoelectron birth energy is hν minus the binding or detachment energy.
        if "electron_energy_eV" in row:
            params["electron_energy_eV"] = float(row["electron_energy_eV"])
        net.add(
            reactants, products, r=sig, r_units="m2", rate_form=RateForm.PHOTO,
            params=params, **common,
        )
    elif form is RateForm.STIMULATED:
        # sigma_eff (cm^2) -> m^2 is handled where gain is applied; store both raw + SI.
        params = {"sigma_eff_cm2": row["sigma_eff_cm2"], "sigma_eff_m2": row["sigma_eff_cm2"] * 1e-4}
        net.add(
            reactants, products, rate_form=RateForm.STIMULATED, params=params, **common,
        )
    else:  # pragma: no cover - RateForm is exhaustive
        raise ValueError(f"unhandled rate_form {form!r}")


def load_reaction_network(
    path: str | Path,
    registry: SpeciesRegistry | None = None,
) -> ReactionNetwork:
    """Build a :class:`ReactionNetwork` from a TOML reaction file."""
    source = Path(path)
    doc = _load_network_document(source, seen=set())
    net = ReactionNetwork(registry)
    for row in doc.get("reaction", []):
        _add_row(net, row)
    return net


def _load_network_document(path: Path, *, seen: set[Path]) -> dict:
    """Load a full TOML network or a provenance-marked reproduction overlay."""
    resolved = path.resolve()
    if resolved in seen:
        raise ValueError(f"cyclic base_network reference at {path}")
    seen.add(resolved)
    with open(path, "rb") as fh:
        doc = tomllib.load(fh)
    if "base_network" not in doc:
        return doc

    base_path = path.parent / str(doc["base_network"])
    base_doc = _load_network_document(base_path, seen=seen)
    rows = deepcopy(base_doc.get("reaction", []))
    by_id = {str(row.get("id", "")): row for row in rows}
    for reaction_id in doc.get("disable_reactions", []):
        if reaction_id not in by_id:
            raise KeyError(f"reproduction overlay disables unknown reaction {reaction_id!r}")
        by_id[reaction_id]["enabled"] = False
    for override in doc.get("override", []):
        reaction_id = str(override.get("id", ""))
        if reaction_id not in by_id:
            raise KeyError(f"reproduction overlay overrides unknown reaction {reaction_id!r}")
        patch = {key: value for key, value in override.items() if key != "id"}
        by_id[reaction_id].update(patch)
    rows.extend(deepcopy(doc.get("reaction", [])))
    return {**doc, "reaction": rows}
