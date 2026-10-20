"""The neural-network estimator trains cleanly and learns an identifiable
channel (a small, fast version of run.py nn)."""

from ml import data, nn, params


def test_network_learns_k_eq():
    model = nn.train(data.by_role("train"), n_train=2500, n_val=500,
                     epochs=25, width=64)
    i = params.NAMES.index("K_eq")
    assert model["r2"][i] > 0.6
    m, s = nn.apply(model)
    assert (s > 0).all()
