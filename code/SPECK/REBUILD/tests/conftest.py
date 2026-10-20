"""Shared fixtures for the REBUILD gates."""

import numpy as np
import pytest

NE_FILL = {"Ne": 0.949, "Kr": 0.05, "F2": 0.001}
AR_FILL = {"Ar": 0.949, "Kr": 0.05, "F2": 0.001}


@pytest.fixture(scope="session")
def loki_ar_fill():
    from REBUILD.MODULES.loki_tier import rates_for_mode
    return rates_for_mode("krf", AR_FILL)


@pytest.fixture(scope="session")
def loki_arf():
    from REBUILD.MODULES.loki_tier import rates_for_mode
    return rates_for_mode("arf_chen", {})


@pytest.fixture(scope="session")
def krf_device():
    """A small 0-D KrF beam device (session-cached; treat as read-only)."""
    from REBUILD.MODULES.engine import Device
    from REBUILD.MODULES.params import engine_parameters
    p = engine_parameters()
    p.set_many({"mode": "krf", "drive": "beam", "pressure_atm": 3.0,
                "t_end_ns": 120.0, "nx": 1, "ny": 1, "nz": 1,
                "Lx_m": 0.0141, "Ly_m": 0.0141, "Lz_m": 0.75,
                "beam_axis": 0}, source="problem")
    return Device(p)


@pytest.fixture(scope="session")
def krf_result(krf_device):
    return krf_device.run(n_save=31)


def rng():
    return np.random.default_rng(42)
