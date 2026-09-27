// Quickstart: solve -c'' = 1 on [0,1], c(0) = c(1) = 0 (exact: x(1-x)/2), first as a
// linear solve, then via Newton with a lambda fill function.
#include <bandsolver/band.hpp>
#include <bandsolver/newton.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

int main() {
    const int n = 1, nj = 11;
    const double h = 1.0 / (nj - 1);

    // 1) Linear solve with an owning BlockSystem.
    bandsolver::BlockSystem sys(n, nj);
    sys.B(0, 0, 0) = 1;
    sys.B(nj - 1, 0, 0) = 1;
    for (int j = 1; j < nj - 1; ++j) {
        sys.A(j, 0, 0) = -1 / (h * h);
        sys.B(j, 0, 0) = 2 / (h * h);
        sys.D(j, 0, 0) = -1 / (h * h);
        sys.G(j, 0) = 1;
    }
    std::vector<double> dc = bandsolver::solve(sys);  // throws bandsolver::Error on failure
    double err = 0;
    for (int j = 0; j < nj; ++j) err = std::fmax(err, std::abs(dc[j] - j * h * (1 - j * h) / 2));
    std::printf("solve: max error=%.2e\n", err);

    // 2) Newton: fill receives c and a zeroed BlockSystem; G = -F(c).
    auto fill = [h](const double* c, bandsolver::BlockSystem& s) {
        const int m = s.nj();
        s.B(0, 0, 0) = 1;     s.G(0, 0) = -c[0];
        s.B(m - 1, 0, 0) = 1; s.G(m - 1, 0) = -c[m - 1];
        for (int j = 1; j < m - 1; ++j) {
            s.A(j, 0, 0) = -1 / (h * h); s.B(j, 0, 0) = 2 / (h * h); s.D(j, 0, 0) = -1 / (h * h);
            s.G(j, 0) = (c[j + 1] - 2 * c[j] + c[j - 1]) / (h * h) + 1;
        }
    };
    std::vector<double> c(nj, 0.0);
    bandsolver::NewtonResult r = bandsolver::newton(n, nj, fill, c.data());
    err = 0;
    for (int j = 0; j < nj; ++j) err = std::fmax(err, std::abs(c[j] - j * h * (1 - j * h) / 2));
    std::printf("newton: status=%s converged=%d iterations=%d max error=%.2e\n",
                bandsolver::to_string(r.status), r.converged, r.iterations, err);
    return r.status == bandsolver::Status::ok ? 0 : 1;
}
