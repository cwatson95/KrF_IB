"""chen2025.py — the verbatim copy stays verbatim (pins, not re-typing).

0.1.25: TABLE_III remains the AS-PRINTED record — the adjudicated
corrections live ONLY in CHEN_ERRATA, and both are pinned here."""

from REBUILD.MODULES.chen2025 import (CHEN_ERRATA, CHEN_PARAMS, TABLE_III,
                                      te_fit)


def test_row_count_is_41():
    assert len(TABLE_III) == 41


def _row(rid):
    for row in TABLE_III:
        if row[0] == rid:
            return row
    raise KeyError(rid)


def test_spot_pins_as_printed():
    assert _row("R9")[4] == (1.56e-13, 0.71, -2.63)
    assert _row("R22")[4] == 4e-14
    assert _row("R27")[4] == 2.5e-43
    assert _row("R37")[4] == 2.381e8
    assert _row("R38")[4] == 4e-20
    assert _row("R40")[4] == 1e-22
    assert _row("R41")[4] == 1e-21


def test_typography_flags_preserved():
    """The R19 exponent and R26 product adjudications (PLAN (46)) must stay
    visible on the rows — deleting the flags would delete the audit."""
    assert "Te^+0.5" in _row("R19")[6] or "PRINTED" in _row("R19")[6]
    assert "endothermic" in _row("R26")[6]


def test_errata_pins_and_table_stays_as_printed():
    """CHEN_ERRATA (0.1.25): the five numeric transcription errors of the
    printed table, adjudicated against primary sources — R13/R14 cm3
    misprint /1e6, R23 = their own R27 pasted onto the excimer row, R24
    x1e6 slip (0.4 ps Ne* lifetime as printed), R36 exponent typo x100.
    TABLE_III itself must keep the PRINTED values (the record of the
    paper); the corrections ride only in the errata block."""
    assert set(CHEN_ERRATA) == {"R13", "R14", "R23", "R24", "R36"}
    for rid, (fixed, printed, basis) in CHEN_ERRATA.items():
        row = _row(rid)
        as_printed = row[4][0] if isinstance(row[4], tuple) else row[4]
        assert as_printed == printed, rid       # the table stays as printed
        assert fixed != printed and basis
    assert CHEN_ERRATA["R13"][0] == 2.08e-14        # /1e6 of the print
    assert CHEN_ERRATA["R14"][0] == 5.61e-15
    assert CHEN_ERRATA["R23"][0] == 1.1e-44         # Ma [13]/SRI/Kannari
    assert CHEN_ERRATA["R24"][0] == 4.0e-46         # R&K/Phelps/Lee triple
    assert CHEN_ERRATA["R36"][0] == 3.5e-43         # Lee/Kannari/Ma x100
    # R23's provenance: the printed value IS their R27 byte-for-byte
    assert CHEN_ERRATA["R23"][1] == _row("R27")[4]


def test_chen_parameters():
    assert CHEN_PARAMS["mixture"] == {"Ne": 0.964, "Ar": 0.035,
                                      "F2": 0.001, "Xe": 0.0001}
    assert CHEN_PARAMS["pressure_atm"] == 3.4
    assert CHEN_PARAMS["R_OC"] == 0.08


def test_te_fit_form():
    # A Te^B exp(C/Te) at Te = 1 -> A * exp(C)
    import math
    assert abs(te_fit(2.0, 0.5, -1.0, 1.0) - 2.0 * math.exp(-1.0)) < 1e-15
    assert te_fit(2.0, 0.5, -1.0, 0.0) == 0.0
