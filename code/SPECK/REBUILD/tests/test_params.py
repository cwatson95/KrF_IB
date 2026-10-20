"""params.py gates: precedence, validation, and the D4 knob ban."""

import pytest

from REBUILD.MODULES.params import (Param, ParamError, ParameterSet,
                                    engine_parameters)


def test_unknown_name_is_an_error_not_a_noop():
    P = engine_parameters()
    with pytest.raises(ParamError):
        P.set("pressure_atmm", 2.0)          # the misspelled-override defect


def test_out_of_range_fails_at_assignment():
    P = engine_parameters()
    with pytest.raises(ParamError):
        P.set("pressure_atm", 99.0)
    with pytest.raises(ParamError):
        P.set("rtol", 1e-2)                  # the loose-tolerance Known-bad


def test_precedence_lower_cannot_clobber_higher():
    P = engine_parameters()
    P.set("t_end_ns", 100.0, source="api")
    P.set("t_end_ns", 55.0, source="problem")   # weaker: ignored
    assert P["t_end_ns"] == 100.0
    assert P.source("t_end_ns") == "api"


def test_provenance_answers_why():
    P = engine_parameters()
    assert P.source("beam_kV") == "default"
    P.set("beam_kV", 500.0, source="cli")
    assert P.source("beam_kV") == "cli"
    assert "[CHANGED: cli]" in P.format_table(only_changed=True)


def test_banned_calibration_knob_is_refused():
    """D4 (ORESTES_CLASS_REBUILD.txt §2): phantom11's DIS_SCALE — a factor
    auto-tuned to hit a target density — must be UNDECLARABLE here."""
    P = ParameterSet()
    with pytest.raises(ParamError):
        P.declare("dis_scale", Param(1.0, "output scale factor"))
    with pytest.raises(ParamError):
        P.declare("boost", Param(1.0, "calibrated against Chow"))
    with pytest.raises(ParamError):
        P.declare("k_fudge_factor", Param(1.0, "adjust rates"))


def test_no_declared_engine_param_is_a_knob():
    # engine_parameters() itself would raise at declaration if one crept in;
    # this asserts the set builds and stays banned-pattern-free
    P = engine_parameters()
    assert len(P.names()) > 15


def test_bool_coercion_and_choices():
    P = ParameterSet()
    P.declare("flag", Param(True, "a boolean"))
    assert P.set("flag", "false") is False
    with pytest.raises(ParamError):
        P.set("flag", "maybe")
    P2 = engine_parameters()
    with pytest.raises(ParamError):
        P2.set("mode", "xef")
