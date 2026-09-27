"""P6: bandsolver.integrate (adaptive BDF 1-2 DAE integrator) from Python."""
import pathlib
import sys

import numpy as np
import pytest

import bandsolver as bs

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "examples"))
import adaptive_integration as ai  # noqa: E402

LAM = np.array([[1.0], [3.0], [10.0]])


def res(t, c, cd):
    return cd + LAM * c


def jac(t, c, cd, alpha):
    return np.zeros((3, 1, 1)), (LAM + alpha)[:, :, None], np.zeros((3, 1, 1))


def err(r):
    return np.abs(r.y[-1, :, 0] - np.exp(-LAM[:, 0])).max()


@pytest.mark.parametrize("order,ratio", [(1, 2.0), (2, 4.0)])
def test_fixed_step_order(order, ratio):
    e = [err(bs.integrate(res, np.ones((3, 1)), [1.0], jacobian=jac, adaptive=False, dt=dt, max_order=order))
         for dt in (0.02, 0.01)]
    assert abs(e[0] / e[1] - ratio) < 0.35


def test_adaptive_outputs_and_accuracy():
    r = bs.integrate(res, np.ones((3, 1)), [0.25, 0.5, 1.0], jacobian=jac, rtol=1e-6, atol=1e-8)
    np.testing.assert_allclose(r.t, [0.25, 0.5, 1.0], rtol=0, atol=1e-15)
    assert r.y.shape == (3, 3, 1)
    assert err(r) < 1e-4 and r.stats["steps"] > 0


def test_fd_jacobian_path_matches_analytic():
    a = bs.integrate(res, np.ones((3, 1)), [1.0], jacobian=jac, rtol=1e-6)
    f = bs.integrate(res, np.ones((3, 1)), [1.0], rtol=1e-6)
    np.testing.assert_allclose(f.y, a.y, rtol=0, atol=1e-10)


def test_binary_electrolyte_dae_reuse_on_off():
    r_on, _, e_on = ai.run(rtol=1e-6, atol=1e-6)
    r_off, _, e_off = ai.run(rtol=1e-6, atol=1e-6, jacobian_reuse=False)
    assert e_on < 2e-3 and e_off < 2e-3                  # vs the analytic series (incl. spatial error)
    assert r_on.stats["factorizations"] < r_off.stats["factorizations"] / 3
    np.testing.assert_allclose(r_on.y, r_off.y, rtol=0, atol=5e-3)
    r_fix, _, e_fix = ai.run(adaptive=False, dt=0.01)
    assert e_fix < 1e-4 and r_fix.stats["factorizations"] <= 5    # constant alpha: one Jacobian serves all steps


def test_fixed_step_large_step_inconsistent_algebraic_start():
    # Regression: fixed step 0.1 s from phi = 0 (inconsistent) needs more than max_newton_iter
    # Newton iterations on the first step; fixed mode must retry with more iterations.
    r, _, e = ai.run(adaptive=False, dt=0.1)
    assert r.stats["steps"] == 50 and e < 5e-2


class Boom(Exception):
    pass


def test_errors():
    def bad(t, c, cd):
        raise Boom("residual")

    with pytest.raises(Boom):
        bs.integrate(bad, np.ones((3, 1)), [1.0])
    with pytest.raises(ValueError, match="dt"):
        bs.integrate(res, np.ones((3, 1)), [1.0], adaptive=False)
    with pytest.raises(ValueError, match="max_order"):
        bs.integrate(res, np.ones((3, 1)), [1.0], max_order=3)
    with pytest.raises(ValueError):
        bs.integrate(res, np.ones((3, 1)), [1.0, 0.5])            # t_out not increasing
    with pytest.raises(bs.IntegrationError) as e:
        bs.integrate(res, np.ones((3, 1)), [1.0], max_steps=3)
    assert e.value.result.stats["steps"] == 3


@pytest.mark.parametrize("dt,t_end,steps", [(5e-3, 5.0, 1000), (0.1, 3.0, 30), (1e-3, 0.7, 700)])
def test_fixed_step_count_has_no_roundoff_sliver(dt, t_end, steps):
    # Regression: accumulated round-off in t used to leave a tiny extra step before t_end.
    r = bs.integrate(res, np.ones((3, 1)), [t_end], jacobian=jac, adaptive=False, dt=dt)
    assert r.stats["steps"] == steps and r.t[-1] == t_end
