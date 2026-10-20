"""Checkpoint / restart -- save a run, resume it BITWISE-EXACTLY.

ADAPTED FROM ``projects/MHD`` (``flare/io.py``), which recreates FLASH's
checkpoint discipline (guide ch. 9): a checkpoint holds the full
double-precision state plus every runtime parameter and its provenance, and a
restart reproduces the run exactly -- flare gates that with
``test_restart_is_bitwise_exact``, and so does this module.

WHY THIS PACKAGE NEEDED IT.  There was no persistence at all: a run existed
only inside one Python process. Three consequences, all now closed:

  * a long 3-D run could not survive an interruption, so the affordable
    problem size was bounded by "how long dare I leave it";
  * a result could not be re-examined after the fact -- only whatever the
    history dict happened to record survived;
  * and there was no way to test reproducibility, which is the check that
    actually catches hidden state. Bitwise-exact restart is a strong gate
    precisely because ANY unsaved mutable state breaks it.

THE RNG IS PART OF THE STATE.  A PIC-MCC run's identity includes its random
stream, so the bit generator's full internal state is saved and restored --
not the seed, which would rewind the stream and silently change the answer.
This is the piece that makes exactness testable rather than approximate.

WHAT IS *NOT* SERIALISED: the MCC operator and its cross-section tables.
Those are large, immutable, and rebuildable from the gas description, which
IS saved. :func:`restore` takes the rebuilt operator as an argument. Pickling
a live object graph would have been less code and much worse: it would freeze
the tables into every checkpoint and break the moment a module was edited.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import es_sim as ES
from . import species3d as SP

FORMAT_VERSION = 1
_PARTICLE_FIELDS = ("x", "y", "z", "vx", "vy", "vz", "w")


def _rng_state_json(rng: np.random.Generator) -> str:
    """Serialise a Generator's bit-generator state (ints may exceed 64 bits)."""
    return json.dumps(rng.bit_generator.state, default=str)


def _rng_from_json(text: str) -> np.random.Generator:
    st = json.loads(text)
    name = st.get("bit_generator", "PCG64")
    bg = getattr(np.random, name)()
    # ints round-tripped through str() by the encoder above
    def _ints(o):
        if isinstance(o, dict):
            return {k: _ints(v) for k, v in o.items()}
        if isinstance(o, str) and (o.isdigit() or
                                   (o.startswith("-") and o[1:].isdigit())):
            return int(o)
        return o
    bg.state = _ints(st)
    return np.random.Generator(bg)


def save(path, sim, gas_description: dict = None, params=None,
         extra: dict = None) -> Path:
    """Write the full state of ``sim`` to ``path`` (a ``.npz``).

    ``gas_description`` should carry whatever is needed to rebuild the MCC
    (e.g. ``{"gas": "Ar", "pressure_atm": 1.0, "T_gas": 300.0}``); it is
    stored verbatim and handed back by :func:`load`.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out = {}

    meta = {
        "format_version": FORMAT_VERSION,
        "shape": list(sim.shape), "h": sim.h, "dt": sim.dt, "bc": sim.bc,
        "time": sim.time, "step_count": sim.step_count,
        "rho_background": sim.rho_background,
        "wall_charge": sim.wall_charge, "wall_energy_J": sim.wall_energy_J,
        "gas_reservoir_eV": getattr(sim, "gas_reservoir_eV", 0.0),
        "mode_axis": sim._mode_axis, "mode_number": sim._mode_number,
        "n_species": len(sim.species),
        "gas": gas_description or {},
        "params": params.as_dict() if params is not None else {},
        "provenance": params.provenance() if params is not None else {},
        "extra": extra or {},
    }
    # phi_lo/phi_hi may be arrays (an (nx,ny) profile) or scalars
    for tag in ("phi_lo", "phi_hi"):
        v = getattr(sim, tag)
        out[f"__{tag}"] = np.asarray(v, dtype=float)
    out["__meta"] = np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    out["__rng"] = np.frombuffer(_rng_state_json(sim.rng).encode(),
                                 dtype=np.uint8)
    out["rho_heavy"] = sim.rho_heavy
    # The FIELDS are saved rather than recomputed on restore.  They are a
    # cache -- step() re-solves before using them -- but they are the cache
    # from BEFORE the last push/move, so re-solving after a restart would
    # produce a different (post-move) field and the state would not match.
    # Saving them also means a checkpoint can be loaded and plotted without
    # re-running the solver.
    out["phi"] = sim.phi
    out["Ex"], out["Ey"], out["Ez"] = sim.E

    if sim.collision_tally is not None:
        out["collision_tally"] = sim.collision_tally
        out["collision_energy_eV"] = sim.collision_energy_eV

    for k, v in sim.history.items():
        out[f"hist__{k}"] = np.asarray(v)

    for i, pop in enumerate(sim.species):
        n = pop.n
        sp = pop.species
        out[f"sp{i}__meta"] = np.frombuffer(json.dumps({
            "name": sp.name, "charge": sp.charge, "mass": sp.mass,
            "n": n, "subcycle": int(getattr(pop, "subcycle", 1)),
            "is_mcc_target": pop is sim.mcc_species,
        }).encode(), dtype=np.uint8)
        for f in _PARTICLE_FIELDS:
            out[f"sp{i}__{f}"] = getattr(pop, f)[:n].copy()

    np.savez(path, **out)
    return path


def load(path) -> dict:
    """Read a checkpoint into a plain dict (no simulation is built)."""
    z = np.load(Path(path), allow_pickle=False)
    meta = json.loads(bytes(z["__meta"]).decode())
    if meta.get("format_version") != FORMAT_VERSION:
        raise ValueError(
            f"checkpoint format {meta.get('format_version')} != "
            f"{FORMAT_VERSION}: written by a different version of this module")
    species = []
    for i in range(meta["n_species"]):
        sm = json.loads(bytes(z[f"sp{i}__meta"]).decode())
        sm["arrays"] = {f: z[f"sp{i}__{f}"] for f in _PARTICLE_FIELDS}
        species.append(sm)
    hist = {k[6:]: z[k] for k in z.files if k.startswith("hist__")}
    return {
        "meta": meta, "species": species, "history": hist,
        "rho_heavy": z["rho_heavy"],
        "phi": z["phi"], "E": (z["Ex"], z["Ey"], z["Ez"]),
        "phi_lo": z["__phi_lo"], "phi_hi": z["__phi_hi"],
        "rng_state": bytes(z["__rng"]).decode(),
        "collision_tally": z["collision_tally"] if "collision_tally" in z.files
        else None,
        "collision_energy_eV": z["collision_energy_eV"]
        if "collision_energy_eV" in z.files else None,
    }


def restore(path, mcc=None) -> ES.ESSim:
    """Rebuild an :class:`ESSim` from a checkpoint.

    ``mcc`` is the collision operator rebuilt from ``meta["gas"]`` -- see the
    module docstring for why it is not serialised. Pass ``None`` for a
    collisionless restart.
    """
    d = load(path)
    m = d["meta"]

    pops, mcc_target = [], None
    for sm in d["species"]:
        sp = SP.Species(sm["name"], sm["charge"], sm["mass"])
        pop = SP.ParticleSpecies(sp, capacity=max(int(sm["n"]) * 2, 16))
        a = sm["arrays"]
        if sm["n"]:
            pop.spawn(a["x"], a["y"], a["z"], a["vx"], a["vy"], a["vz"], a["w"])
        pop.subcycle = int(sm["subcycle"])
        pops.append(pop)
        if sm.get("is_mcc_target"):
            mcc_target = pop

    def _scalar_or_array(arr):
        arr = np.asarray(arr)
        return float(arr) if arr.ndim == 0 else arr

    nx, ny, nz = m["shape"]
    sim = ES.ESSim(nx, ny, nz, m["h"], SP.SpeciesSet(pops), m["dt"],
                   bc=m["bc"],
                   phi_lo=_scalar_or_array(d["phi_lo"]),
                   phi_hi=_scalar_or_array(d["phi_hi"]),
                   background_charge_density=m["rho_background"],
                   mcc=mcc, mcc_species=mcc_target,
                   rng=_rng_from_json(d["rng_state"]))

    sim.time = m["time"]
    sim.step_count = m["step_count"]
    sim.wall_charge = m["wall_charge"]
    sim.wall_energy_J = m["wall_energy_J"]
    sim.gas_reservoir_eV = m["gas_reservoir_eV"]
    sim.set_mode(m["mode_axis"], m["mode_number"])
    sim.rho_heavy = np.array(d["rho_heavy"])
    sim.phi = np.array(d["phi"])
    sim.E = tuple(np.array(a) for a in d["E"])
    if d["collision_tally"] is not None:
        sim.collision_tally = np.array(d["collision_tally"])
        sim.collision_energy_eV = np.array(d["collision_energy_eV"])
    sim.history = {k: list(v) for k, v in d["history"].items()}
    return sim


def state_digest(sim) -> dict:
    """A compact fingerprint of a run's state -- the restart-equality check.

    Deliberately built from RAW ARRAYS, not from summary scalars: a digest of
    means and totals would call two runs identical while individual particles
    differed. Compare with ``==`` on the returned dict.
    """
    parts = {
        "time": repr(sim.time), "step": sim.step_count,
        "wall_charge": repr(sim.wall_charge),
        "wall_energy": repr(sim.wall_energy_J),
        "rho_heavy": _hash(sim.rho_heavy),
        "phi": _hash(sim.phi),
    }
    for i, pop in enumerate(sim.species):
        parts[f"sp{i}_n"] = pop.n
        for f in _PARTICLE_FIELDS:
            parts[f"sp{i}_{f}"] = _hash(getattr(pop, f)[:pop.n])
    if sim.collision_tally is not None:
        parts["tally"] = _hash(sim.collision_tally.astype(float))
    parts["rng"] = _rng_state_json(sim.rng)
    return parts


def _hash(arr) -> str:
    import hashlib
    a = np.ascontiguousarray(np.asarray(arr, dtype=float))
    return hashlib.sha256(a.tobytes()).hexdigest()[:16]
