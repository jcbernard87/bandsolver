"""The worked examples double as validation tests (observed convergence orders)."""
import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "examples"))

import fd_jacobian  # noqa: E402
import nonlinear_bvp  # noqa: E402
import reaction_diffusion  # noqa: E402
import transient_diffusion  # noqa: E402
import xy_boundary  # noqa: E402


def orders(rows):
    return np.array([r["order"] for r in rows[1:]])


def test_nonlinear_bvp_second_order(backend):
    rows = nonlinear_bvp.run(backend, grids=(41, 81, 161))
    assert np.all(np.abs(orders(rows) - 2) < 0.05)
    assert all(r["iterations"] <= 8 for r in rows)


def test_reaction_diffusion_second_order_and_quadratic(backend):
    rows = reaction_diffusion.run(backend, grids=(41, 81, 161))
    assert np.all(np.abs(orders(rows) - 2) < 0.05)
    s = rows[1]["steps"]
    assert s[-1] < 1e-12 and s[-2] < 10 * s[-3] ** 2


def test_transient_implicit_euler_first_order(backend):
    rows = transient_diffusion.run(backend, nj=201, dts=(0.01, 0.005, 0.0025))
    assert np.all(np.abs(orders(rows) - 1) < 0.1)


def test_xy_closure_restores_second_order(backend):
    out = xy_boundary.run(backend, grids=(321, 641, 1281))
    assert np.all(np.abs(orders(out["second_order_XY"]) - 2) < 0.05)
    assert np.all(np.abs(orders(out["first_order"]) - 1) < 0.05)


@pytest.mark.parametrize("module", [nonlinear_bvp, reaction_diffusion])
def test_backends_give_identical_errors(module):
    a = module.run("cpp", grids=(41,))
    b = module.run("fortran", grids=(41,))
    assert a[0]["iterations"] == b[0]["iterations"]
    assert abs(a[0]["error"] - b[0]["error"]) <= 1e-14


def test_fd_example_matches_analytic_and_flags_bug(backend):
    rows, good, bad = fd_jacobian.run(backend, grids=(41, 81, 161))
    for r in rows:
        assert r["it_fd"] == r["it_analytic"] and r["diff"] < 1e-12 and r["evals"] == 10 * r["it_fd"]
    assert np.all(np.abs(orders(rows) - 2) < 0.05)
    assert good.max_error < 1e-5
    name, m = bad.worst()
    assert (name, m.row, m.col) == ("B", 1, 0) and m.error > 1e-2 and np.sign(m.user) != np.sign(m.fd)
