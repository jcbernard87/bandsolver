import numpy as np
import pytest

import bandsolver as bs
from conftest import assemble_dense, random_system


@pytest.mark.parametrize("pivot", ["partial", "legacy"])
@pytest.mark.parametrize("xy", [False, True])
@pytest.mark.parametrize("n,nj", [(1, 3), (2, 3), (3, 4), (6, 30), (4, 200)])
def test_matches_dense_solve(rng, backend, pivot, xy, n, nj):
    A, B, D, G, X, Y = random_system(rng, n, nj, xy)
    dc = bs.solve(A, B, D, G, X, Y, pivot=pivot, backend=backend)
    assert dc.shape == (nj, n)
    K = assemble_dense(A, B, D, X, Y)
    ref = np.linalg.solve(K, G.ravel()).reshape(nj, n)
    np.testing.assert_allclose(dc, ref, rtol=0, atol=1e-12 * np.abs(ref).max())
    berr = np.abs(K @ dc.ravel() - G.ravel()).max() / (np.abs(K).sum(1).max() * np.abs(dc).max() + np.abs(G).max())
    assert berr < 1e-14


def test_backends_agree(rng):
    A, B, D, G, X, Y = random_system(rng, 5, 60, True)
    for pivot in ("partial", "legacy"):
        a = bs.solve(A, B, D, G, X, Y, pivot=pivot, backend="cpp")
        b = bs.solve(A, B, D, G, X, Y, pivot=pivot, backend="fortran")
        np.testing.assert_allclose(a, b, rtol=0, atol=64 * np.finfo(float).eps * np.abs(a).max())


def test_inputs_not_modified_and_layouts_accepted(rng, backend):
    A, B, D, G, X, Y = random_system(rng, 3, 10, True)
    copies = [a.copy() for a in (A, B, D, G, X, Y)]
    ref = bs.solve(A, B, D, G, X, Y, backend=backend)
    for a, c in zip((A, B, D, G, X, Y), copies):
        np.testing.assert_array_equal(a, c)
    # Fortran-ordered and strided views give the same answer.
    Bf = np.asfortranarray(B)
    Gs = np.repeat(G, 2, axis=0)[::2]
    np.testing.assert_array_equal(bs.solve(A, Bf, D, Gs, X, Y, backend=backend), ref)


@pytest.mark.parametrize("pivot", ["partial", "legacy"])
@pytest.mark.parametrize("node", [0, 4, 8])
def test_singular_block_reports_node(rng, backend, pivot, node):
    A, B, D, G, _, _ = random_system(rng, 2, 9, False)
    A[node] = 0
    B[node, 1] = 2 * B[node, 0]
    with pytest.raises(bs.SingularBlockError) as e:
        bs.solve(A, B, D, G, pivot=pivot, backend=backend)
    assert e.value.node == node


def test_non_finite(rng, backend):
    A, B, D, G, _, _ = random_system(rng, 2, 5, False)
    G[2, 1] = np.nan
    with pytest.raises(bs.NonFiniteError):
        bs.solve(A, B, D, G, backend=backend)


def test_argument_validation(rng):
    A, B, D, G, X, Y = random_system(rng, 2, 5, True)
    with pytest.raises(ValueError, match="shape"):
        bs.solve(A[:4], B, D, G)
    with pytest.raises(ValueError, match="shape"):
        bs.solve(A, B, D, G[:, :1])
    with pytest.raises(ValueError, match="shape"):
        bs.solve(A, B, D, G, X=np.eye(3))
    with pytest.raises(ValueError, match="nodes"):
        bs.solve(A[:2], B[:2], D[:2], G[:2])
    with pytest.raises(ValueError, match="backend"):
        bs.solve(A, B, D, G, backend="julia")
    with pytest.raises(ValueError, match="pivot"):
        bs.solve(A, B, D, G, pivot="full")


@pytest.mark.parametrize("xy", [False, True])
@pytest.mark.parametrize("n,nj", [(1, 3), (3, 4), (7, 60)])
def test_fortran_kernels_bit_identical(rng, xy, n, nj):
    A, B, D, G, X, Y = random_system(rng, n, nj, xy)
    fast = bs.solve(A, B, D, G, X, Y, backend="fortran", kernel="fast")
    ref = bs.solve(A, B, D, G, X, Y, backend="fortran", kernel="reference")
    np.testing.assert_array_equal(fast, ref)
    with pytest.raises(ValueError, match="kernel"):
        bs.solve(A, B, D, G, backend="fortran", kernel="turbo")
