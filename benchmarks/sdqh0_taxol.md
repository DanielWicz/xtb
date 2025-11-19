# SDQH0 Benchmark – Taxol (GFN2-xTB)

- **Date**: 2025-11-19
- **Host**: NEBULA (same node as development environment)
- **Command**:
  ```bash
  XTB_SDQH0_TIMING=1 OMP_NUM_THREADS=4 build/xtb assets/inputs/xyz/taxol.xyz --gfn2
  ```
- **Relevant xtb settings**: default integral cutoff from accuracy 1.0, SCF converged normally.
- **Key timings (from `build_SDQH0` diagnostics)**:

| Revision | Translation (s) | Transform (s) | Screening (s) | Diagonal (s) |
| --- | --- | --- | --- | --- |
| Baseline (pre-change) | 0.029564 | 0.003730 | 0.000384 | 0.001667 |
| After skip+rowstart optimizations | 0.025725 | 0.003465 | 0.000403 | 0.001383 |
- **Notes**: run produces standard scratch files in the build directory. Delete or re-use them before rerunning to avoid accidental file churn.

This measurement is the performance reference for the upcoming optimizations in `src/xtb/hamiltonian.F90`.
