# bandsolver

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22997640.svg)](https://doi.org/10.5281/zenodo.22997640)
[![CI](https://github.com/jcbernard87/bandsolver/actions/workflows/ci.yml/badge.svg)](https://github.com/jcbernard87/bandsolver/actions/workflows/ci.yml)

A standalone implementation of John Newman's **BAND** algorithm for coupled, nonlinear one-dimensional boundary-value and transient problems. The algorithm is described in Appendix C of *Electrochemical Systems* (Newman & Thomas-Alyea). The library provides two interchangeable native cores and one Python interface:

| Layer | Language | Entry points |
|---|---|---|
| Fortran core | Fortran 2008 | `bandsolver_kernel` (`band_solve`), `bandsolver_newton` (`band_newton`, `band_problem`) |
| C interface to the Fortran core | C | `fortran/include/bandsolver_f.h` (`bandsolver_f_solve`, `bandsolver_f_newton`) |
| C++ core | C++17 | `<bandsolver/band.hpp>`, `<bandsolver/newton.hpp>` |
| Python | numpy | `bandsolver.solve`, `bandsolver.newton`, `backend="cpp"` or `"fortran"` |

The library contains no battery or electrochemistry model. You write the residual and Jacobian for your own equations; the library assembles the block system, solves it, and iterates with Newton's method.

## What it solves

With `nj ≥ 3` nodes and `n` unknowns per node, each Newton step solves a block-tridiagonal system. The first and last rows may reach one extra node, through the blocks `X` and `Y`:

```
node 0:        B₀ Δc₀ + D₀ Δc₁ + X Δc₂                         = G₀
interior j:    Aⱼ Δcⱼ₋₁ + Bⱼ Δcⱼ + Dⱼ Δcⱼ₊₁                    = Gⱼ
node nj−1:     Y Δc_{nj−3} + A Δc_{nj−2} + B Δc_{nj−1}          = G_{nj−1}
```

Here `G = −F(c)` is the negative residual and `A, B, D` are its Jacobian blocks with respect to the neighbouring, local and next unknowns. The cost is O(nj·n³) time and O(nj·n²) memory. See [docs/math.md](docs/math.md).

## Build and install

Requirements: CMake ≥ 3.20, a Fortran 2008 compiler (tested with gfortran 14), and a C++17 compiler (tested with Apple clang 16). The Python package also needs Python ≥ 3.9 and numpy.

```sh
# Fortran + C++ libraries, tests, and the quickstart programs
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build

# Python package from source (builds both backends into one extension)
python -m pip install .
python -m pytest
```

**Prebuilt wheels** for Linux (x86_64, aarch64), macOS (arm64, x86_64) and Windows (x86_64) are attached to each [GitHub release](https://github.com/jcbernard87/bandsolver/releases). They cover CPython 3.9 and later and bundle the Fortran runtime, so no compiler is needed:

```sh
pip install bandsolver-<version>-<python>-<platform>.whl
```

## Quickstart

**Python**:
```python
import numpy as np
import bandsolver as bs

nj = 11; h = 1 / (nj - 1)

def fill(c):                     # -c'' = 1, c(0) = c(1) = 0
    A = np.zeros((nj, 1, 1)); B = np.zeros((nj, 1, 1)); D = np.zeros((nj, 1, 1)); G = np.zeros((nj, 1))
    A[1:-1] = D[1:-1] = -1 / h**2
    B[1:-1] = 2 / h**2
    G[1:-1, 0] = (c[2:, 0] - 2 * c[1:-1, 0] + c[:-2, 0]) / h**2 + 1   # G = -F(c)
    B[0] = B[-1] = 1
    G[0, 0], G[-1, 0] = -c[0, 0], -c[-1, 0]
    return A, B, D, G

r = bs.newton(fill, np.zeros((nj, 1)), backend="fortran")   # or backend="cpp"
print(r.converged, r.iterations, r.c[:, 0])
```

If you only want to write the residual `F(c)`, let the library build the Jacobian by finite differences. This takes 3n + 1 residual calls per iteration, independent of the grid size:

```python
r = bs.newton_fd(lambda c: -fill(c)[3], np.zeros((nj, 1)))   # residual only
print(bs.check_jacobian(fill, r.c).max_error)                # verify a hand-written Jacobian
```

**Fortran**: see [examples/fortran_quickstart.f90](examples/fortran_quickstart.f90). Extend `band_problem`, implement `fill`, and call `band_newton`; or call `band_solve` directly.

**C++**: see [examples/cpp_quickstart.cpp](examples/cpp_quickstart.cpp). Fill a `BlockSystem` and call `bandsolver::solve`, or pass a lambda to `bandsolver::newton`.

More complete examples are in [examples/](examples): a nonlinear BVP, a coupled 3-unknown differential-algebraic system, transient diffusion with implicit Euler, second-order Neumann boundaries through `X`/`Y`, and residual-only solving with Jacobian checking (`fd_jacobian.py`). Each reports its observed convergence order.

## Citing

See [CITATION.cff](CITATION.cff); GitHub's "Cite this repository" button uses it. Archived on Zenodo: [doi:10.5281/zenodo.22997640](https://doi.org/10.5281/zenodo.22997640), which resolves to the latest version. Version 0.1.1 is [doi:10.5281/zenodo.22997641](https://doi.org/10.5281/zenodo.22997641). Please also cite Newman's method: J. Newman, *Ind. Eng. Chem. Fundam.* 7, 514 (1968), and *Electrochemical Systems*, 3rd ed., Appendix C.

## Tutorial notebooks

Worked, fully explained examples live in [`notebooks/`](notebooks). Each one sets out the equations in LaTeX, derives the discretization, implements it, and validates it against an analytic solution:

| Notebook | What it shows |
|---|---|
| [1 · Getting started](notebooks/01_getting_started.ipynb) | The BAND block structure; `solve`; a nonlinear BVP with a hand-written Jacobian and quadratic Newton convergence; `newton_fd`; how `check_jacobian` catches Jacobian bugs; second-order accuracy |
| [2 · Binary electrolyte](notebooks/02_binary_electrolyte.ipynb) | **Coupled transient PDEs**: Nernst–Planck transport of concentration and potential under constant current, implicit Euler, validated against the analytic series solution and steady-state relations |
| [3 · Porous electrode](notebooks/03_porous_electrode.ipynb) | The Newman–Tobias current distribution: solid and electrolyte potentials coupled by Butler–Volmer kinetics, solved from the residual alone, validated against the linear-kinetics analytic solution |
| [4 · Finite differences vs finite volumes](notebooks/04_fd_vs_fv.ipynb) | The same electrode/separator diffusion problem discretized both ways, with the same solver. Nodal FD converges to the wrong answer (+8 % salt); FD with interface matching is first order; FV is second order and conserves exactly. Also covers a graded mesh, when FD is fine, and a decision guide |

To run them locally: `pip install bandsolver[notebooks]` (or `pip install jupyterlab matplotlib`), then `jupyter lab notebooks/`. CI re-executes every notebook.

## Documentation

- [docs/math.md](docs/math.md): formulation, elimination algorithm, pivoting, and Newton convergence.
- [docs/api.md](docs/api.md): Fortran, C, C++ and Python API reference, layouts, and error codes.
- [docs/validation.md](docs/validation.md): measured accuracy and convergence, cross-backend and legacy agreement.
- [docs/provenance.md](docs/provenance.md): origin of the algorithm, the frozen legacy oracle, and the open licensing questions.

## Relationship to the historical code

The kernel reimplements the `BAND`/`MATINV` routines used in the author's PhD battery-model Fortran sources, which follow Appendix C. With `pivot="legacy"` it reproduces a frozen, byte-exact copy of the archival kernel **bit for bit** on the test platform. That copy is kept privately as a test oracle and is not distributed, because the archival routines closely follow the copyrighted textbook listing. The default `pivot="partial"` uses standard row partial pivoting. Two deliberate fixes are documented in [docs/math.md](docs/math.md#differences-from-the-archival-kernel): a singular-block threshold, and a correct 3-node X+Y case.

## License

BSD 3-Clause; see [LICENSE](LICENSE). If you use bandsolver in published work, please cite Newman's Appendix C (see [docs/provenance.md](docs/provenance.md)).

## Status and limitations

- The current release is 0.1.1 (see [CHANGELOG.md](CHANGELOG.md)); the source and prebuilt wheels are on GitHub. It is not yet on PyPI.
- CI builds and tests on every push: Linux x86_64/aarch64 (gfortran 14), macOS arm64/x86_64 (gfortran 14 + Apple clang), and Windows x86_64 (Intel ifx 2025 + MSVC).
- Time integration is left to user code (see `examples/transient_diffusion.py`). No adaptive step control or DAE integrator is included.

