# Provenance, attribution, and open release decisions

## Origin

- **Algorithm:** J. Newman's BAND method for coupled ordinary differential equations. It is published as Appendix C of J. Newman and K. E. Thomas-Alyea, *Electrochemical Systems*, 3rd ed. (Wiley-Interscience, 2004; ISBN 0-471-47756-7), Appendix C, printed pp. 611–622. The edition was confirmed by the author and matches the title and copyright pages of the supplied PDF. The book is © 2004 John Wiley & Sons, all rights reserved. Its §C.4 offers the Fortran 77 subroutines `BAND(J)` and `MATINV` for readers to call from their own programs, but grants no redistribution license. That appendix includes a Fortran listing of `BAND(J)` and `MATINV`. The method goes back to Newman's earlier work (Newman, *Ind. Eng. Chem. Fundam.* 7, 514 (1968)). **Verify the 1968 citation against the original before publication.**
- **Archival implementation:** the author's PhD-era battery-model Fortran sources in `~/Desktop/Research` and related folders. Their shared `BAND`/`MATINV` text was inventoried and compared in `../../BAND_SOURCE_INVENTORY.md` and `../../NUMERICAL_SOLVER_AUDIT.md` (§8.18 compares it with Appendix C).
- **This library:** a new implementation (Fortran 2008, C++17 and Python). It has explicit interfaces, no global state, and status reporting. The elimination and the legacy pivot rule deliberately follow the operation order of the archival and Appendix C routines, so that historical results can be reproduced.

## Legacy test oracle (private, not distributed)

The archival kernel is textually close to the copyrighted Appendix C listing, so it is **not part of this repository or any release**. The author decided this on 2026-09-26. A byte-exact extract is kept privately, next to the repository, in `../bandsolver-legacy-oracle/`. That folder also holds its own `PROVENANCE.md` and the wrapper that drives it. CMake builds the two comparison tests (`tests/fortran/legacy/`) only when that folder is found; the location can be overridden with `-DBANDSOLVER_LEGACY_DIR=...`. Without it, all other tests still build and run.

Facts recorded from the oracle, which contain no code:

| Item | Value |
|---|---|
| Archive source | `Research/Modelling/MnO2/Old Work/ZnMn02_v8.f95`, SHA-256 `8abba38490cc8feefd282841b441a9c77d4f49677a6c90fda76a13a702cc191a` |
| Extract | lines 1983–2100 (`MATINV` 1983–2037, `BAND` 2041–2100), SHA-256 `99eb0fb89443aba84debb2329d1478c48036239df85860268288a1f564b4ef8f` |
| Normalized routine hashes | `MATINV` `3194efd4…48f8a3`, `BAND` `e8e1ce3a…958b8c`; both match the archive inventory |
| Library agreement | `pivot="legacy"` matches bit for bit (0 ulp); `pivot="partial"` differs by ≤ 9.6e-16 relative (see `validation.md`) |
| Archival defect | with nj = 3 and both X and Y nonzero, the result is wrong (backward error ≈ 3e-3); see `math.md` |

## Decisions for the author (not made by the implementer)

1. **License — decided 2026-09-26: BSD-3-Clause** (`LICENSE`). It covers this repository's own code: the reimplementation, tests, examples and docs. The BAND algorithm itself is a published method and is cited, not licensed.
2. **Shipping the legacy extract — decided 2026-09-26: not shipped.** See above. To publish it later would need permission from Wiley (permreq@wiley.com).
3. **Institutional rights — resolved 2026-09-26 by the author.** Columbia is not involved in this library's development. The underlying algorithm is Newman's published method, and this implementation is the author's own work, so no institutional claim applies. The author approved public release.
4. **Attribution and citation text**, and whether to add `CITATION.cff` or a Zenodo DOI at release.
5. **Package name** on PyPI/conda: `bandsolver` has not been checked for availability.
6. **Publication target**, for example JOSS or an electrochemistry-methods venue, and the scope of the accompanying paper. Stage 7 benchmarks and Track 2 model reproductions would support that paper.

No push, upload, or package release has been made. Each requires explicit authorization (see `../../AGENT_INSTRUCTIONS.md`, Stage 8).
