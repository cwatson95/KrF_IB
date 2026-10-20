"""S1 — the 76-entry history converted to gates, each citing its entry.

ORESTES_CLASS_REBUILD.txt §6 S1: "convert every recorded defect into a
test... every Known-bad item either has a test or a written statement of why
it cannot apply." The statements-of-inapplicability live at the bottom.
"""

from pathlib import Path

import pytest

MODULES = Path(__file__).resolve().parents[1] / "MODULES"


def _source_of(name):
    return (MODULES / name).read_text()


def _code_only(path):
    """Strings/comments stripped: prose describing a defect is not the
    defect."""
    import tokenize
    out = []
    with open(path, "rb") as f:
        for tok in tokenize.tokenize(f.readline):
            if tok.type in (tokenize.STRING, tokenize.COMMENT):
                continue
            out.append(tok.string)
    return " ".join(out)


def test_plan72_kannari_relaxation_unconditional():
    """PLAN (72) / root Known-bad: phantom11's K_PS_RELAX defaults OFF as a
    parity artifact (tau ~0.45 ns at 1 atm vs 25 ns radiative — 3.4x in
    baseline light, W-bound breached). Here: rows present, no flag, and no
    LIT_QUENCH-style gate exists anywhere in the engine."""
    from REBUILD.MODULES.network import build_network
    rows = build_network("krf", {"Kr*", "Kr**", "Kr+", "F-", "F2(a3)",
                                 "F2(A1)", "F2+", "Ne*", "Ne**", "Ne+"})
    by_rid = {r.rid: r for r in rows}
    assert "KPS_Ar" in by_rid and "APS_Ar" in by_rid
    for f in MODULES.glob("*.py"):
        assert "LIT_QUENCH" not in _code_only(f), f.name


def test_plan70_class_no_unledgered_conversion():
    """PLAN (70): phantom11's electron-heating -> light chain is internally
    inconsistent by ~18x WITH NO INVARIANT TO CATCH IT. The structural fix
    is the ledger; its closure and teeth gates live in test_ledger.py. This
    test asserts the invariant EXISTS in the RHS path (the channel rows are
    integrated state, not post-hoc bookkeeping)."""
    src = _source_of("engine.py")
    assert "i_ledger0" in src and "drive_rate" in src
    from REBUILD.MODULES.ledger import CHANNELS
    assert CHANNELS == ("drive_in", "laser_out", "cavity_loss",
                        "radiated_escape")


def test_r115_class_killed_by_mechanical_balance():
    """Root Known-bad: the parity harness was structurally blind to
    Ar_g-weighted terms — a sign error on r115 (40% of dArF) survived.
    Here NO row can carry an element/charge imbalance (test_network teeth)
    and every RHS evaluation conserves each element to machine precision
    (test_chemistry.test_elemental_mass_conserved_exactly). This entry
    records the mapping; the teeth live in those files."""
    from REBUILD.MODULES.network import NetworkError, Reaction, check_reaction
    with pytest.raises(NetworkError):
        check_reaction(Reaction("R115LIKE", ("Ar2F*", "Ar"),
                                ("ArF*", "Ar"), ("const", 1e-16), "EST",
                                "drops an Ar — the class of typo r115 was"))


def test_plan74_class_no_module_globals_mutated():
    """PLAN (74): krf0d mutating phantom11 globals cost 10 of 17 suite
    failures. Engine modules must not hold mutable module-level state that
    a run rewrites; test_architecture diffs it mechanically. Here: the
    known offender pattern (a module-level CHOW_EN_TD_DEFAULT-style
    read-write global) must not exist."""
    for f in MODULES.glob("*.py"):
        code = _code_only(f)
        assert "CHOW_EN_TD_DEFAULT" not in code, f.name
        assert "DIS_SCALE" not in code, f.name


def test_scipy_tolerance_known_bad_enforced():
    """Root Known-bad: never integrate at loose tolerance and compare
    tight. The rtol declaration refuses anything looser than 1e-3."""
    from REBUILD.MODULES.params import ParamError, engine_parameters
    P = engine_parameters()
    with pytest.raises(ParamError):
        P.set("rtol", 5e-3)


def test_ne_is_not_y0():
    """Root Known-bad: 'ne is not y[0]' — electron density comes from
    quasineutrality, never from a state slot."""
    from REBUILD.MODULES.species import KRF_MODE_SPECIES, StateLayout
    import numpy as np
    lay = StateLayout(KRF_MODE_SPECIES, 1)
    y = np.zeros(lay.size)
    y[0] = 7.777e20                              # first species slot (Ne)
    dens = lay.density_block(y)
    ne = lay.electron_density(dens)
    assert ne[0] == 1.0                          # floored: no ions present


def test_no_engine_import_of_sibling_packages():
    """PIC_MCC_3D Known-bad (owner warning 2026-08-03): export and compare,
    never inject. The engine must not import PIC_MCC_3D, phantom11, or
    spectre — the KINETIC coupling is a data file. Checked on the AST, so
    a docstring saying 'adapted from spectre' is documentation, not a
    dependency."""
    import ast
    for f in MODULES.glob("*.py"):
        tree = ast.parse(f.read_text())
        for node in ast.walk(tree):
            mods = []
            if isinstance(node, ast.Import):
                mods = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                mods = [node.module]
            for mod in mods:
                assert mod.split(".")[0] not in (
                    "PIC_MCC_3D", "phantom11", "spectre"), (f.name, mod)


def test_beam_w_never_defaulted():
    """PIC_MCC_3D beam3d rule: a W-value input must be explicit, so no run
    can quietly turn W into an input and be quoted for it. Here the
    fallback route is recorded on the yields object and printed."""
    from REBUILD.MODULES.beam import w_table_yields
    y = w_table_yields({"Ar": 1.0})
    assert y.mode == "w_table" and "Lee" in y.cite


# --- Known-bad items that CANNOT apply here, with the reason ---------------
INAPPLICABLE = {
    "og-pyfiles parity gates": "no frozen reference exists for this "
        "package; D1 is closed by construction (no parity harness at all)",
    "phantom11 benchmark orchestration (KrF_phantom11.py)": "never copied; "
        "the calibrated path is banned and no calibration knob can be "
        "declared (params._BANNED_KNOB)",
    "LoKI-MC BSR-Ne beyond 150 Td / Phelps Ne": "this package runs no "
        "LoKI-MC sweeps; it consumes curated outputs whose provenance "
        "README records the Biagi choice",
    "np.add.at deposition slowness": "no particle deposition exists here; "
        "the chemistry scatter is per-row indexed adds",
    "CPML / DC-bearing sources / overdense-mirror": "no FDTD in this "
        "package (ARPAE's instrument owns the EM path)",
}


def test_inapplicability_statements_exist():
    assert len(INAPPLICABLE) >= 5
    for k, why in INAPPLICABLE.items():
        assert len(why) > 20, k
