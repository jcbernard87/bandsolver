import numpy as np
import pytest


def assemble_dense(A, B, D, X=None, Y=None):
    """Dense K for the Appendix C system, built independently of the solvers."""
    nj, n, _ = B.shape
    K = np.zeros((nj * n, nj * n))
    for j in range(nj):
        r = slice(j * n, (j + 1) * n)
        K[r, j * n:(j + 1) * n] += B[j]
        if j > 0:
            K[r, (j - 1) * n:j * n] += A[j]
        if j < nj - 1:
            K[r, (j + 1) * n:(j + 2) * n] += D[j]
    if X is not None:
        K[0:n, 2 * n:3 * n] += X
    if Y is not None:
        K[(nj - 1) * n:nj * n, (nj - 3) * n:(nj - 2) * n] += Y
    return K


def random_system(rng, n, nj, xy):
    A = rng.uniform(-1, 1, (nj, n, n))
    B = rng.uniform(-1, 1, (nj, n, n)) + (3 * n + 3) * np.eye(n)
    D = rng.uniform(-1, 1, (nj, n, n))
    G = rng.uniform(-1, 1, (nj, n))
    X = rng.uniform(-1, 1, (n, n)) if xy else None
    Y = rng.uniform(-1, 1, (n, n)) if xy else None
    return A, B, D, G, X, Y


@pytest.fixture
def rng():
    return np.random.default_rng(20260926)


@pytest.fixture(params=["cpp", "fortran"])
def backend(request):
    return request.param
