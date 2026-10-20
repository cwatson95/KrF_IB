"""The dataset registry: complete entries, and the predictive-only guard."""

import pytest

from ml import ROOT, data, fit


def test_every_dataset_is_complete():
    names = set()
    for d in data.all_datasets():
        assert d.name not in names
        names.add(d.name)
        assert d.role in data.ROLES and d.source and d.sigma > 0
        assert d.points and all(v > 0 for v in d.values)
        if d.file:
            assert (ROOT / d.file).is_file(), d.file


def test_device_data_never_train():
    """Predictive-only policy (KrF_and_LoKI/agents/CLAUDE.md): laser
    performance scores the model and never trains it. Changing this is the
    owner's decision, made in the policy text first."""
    for d in data.all_datasets():
        if d.kind in data.DEVICE_KINDS:
            assert d.role == "device", d.name
    with pytest.raises(ValueError):
        fit.fit(data.by_role("device"))


def test_train_set_is_the_consistent_kinetic_set():
    assert [d.name for d in data.by_role("train")] == [
        "eden_kr", "eden_f2", "ksu_fig2", "ksu_fig1", "ys_ir300"]
