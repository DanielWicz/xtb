# SDQH0 Benchmarks for 19 Nov 2025 Commits

All runs executed on NEBULA with `meson setup build --buildtype=release`, `ninja -C build`, `ninja -C build test`, and the benchmark command:

```bash
XTB_SDQH0_TIMING=1 OMP_NUM_THREADS=4 build/xtb assets/inputs/xyz/taxol.xyz --gfn2
```

| Commit | Description | Build/Test | Translation (s) | Transform (s) | Screening (s) | Diagonal (s) | SCC total (s) | Benchmark Log |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 14ebcad | Stream Coulomb gradients without large buffers | ✅ | n/a¹ | n/a¹ | n/a¹ | n/a¹ | 0.364 | [log](runs/benchmark_14ebcad.log) |
| c3fe6aa | Free SCF integral tensors earlier | ✅ | n/a¹ | n/a¹ | n/a¹ | n/a¹ | 0.360 | [log](runs/benchmark_c3fe6aa.log) |
| 0d1805f | Fix self-energy shell loops | ✅ | n/a¹ | n/a¹ | n/a¹ | n/a¹ | 0.366 | [log](runs/benchmark_0d1805f.log) |
| 0ff2702 | Parallelize integral counters | ✅ | n/a¹ | n/a¹ | n/a¹ | n/a¹ | 0.360 | [log](runs/benchmark_0ff2702.log) |
| 1b8857c | Tile shell pair accumulation | ✅ | n/a¹ | n/a¹ | n/a¹ | n/a¹ | 0.365 | [log](runs/benchmark_1b8857c.log) |
| f09ee64 | Tile diagonal shell accumulation | ✅ | n/a¹ | n/a¹ | n/a¹ | n/a¹ | 0.360 | [log](runs/benchmark_f09ee64.log) |
| 8fbd559 | Optimize SDQH0 loops and record benchmark | ✅ | 0.030528 | 0.003767 | 0.000387 | 0.001660 | 0.368 | [log](runs/benchmark_8fbd559.log) |
| 93bfc74 | Use packed H0 in gradient build | ✅ | 0.030480 | 0.003803 | 0.000388 | 0.001742 | 0.365 | [log](runs/benchmark_93bfc74.log) |

See `2025-11-19_walltime.md` for per-system wall-time data (including `--bhess` runs).

¹These older commits predate the `XTB_SDQH0_TIMING` diagnostics, so the benchmark log does not contain per-phase timings even though the full xtb run is recorded.

`✅` indicates `ninja -C build test` completed with only the expected failures listed in the log.
