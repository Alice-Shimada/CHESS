# Simone Zoia F3: Mathematica versus FORM/FLINT batch evaluation

Date: 2026-08-17

## Problem and timing boundary

This reproduces the large non-canonical example preserved under
`ttW_physical_variables_2608_03746`: the 131-dimensional F3 master-integral
system supplied with arXiv:2504.13011, transported from X1 to X2 along the
physical-coordinate path.  Its connection is

```text
B0(t) + epsilon B1(t) + epsilon^2 B2(t)
```

and the propagated boundary contains epsilon orders 0 through 4.  Both routes
use 96 Lobatto intervals, 100 decimal digits, and eight workers/threads.  Both
are compared with the independently loaded X2 reference values.

The expensive part is the evaluation of 639 large irreducible polynomials at
all ordinary nodes.  One-time symbolic pullback, FORM optimisation, and SLP
construction are excluded from both routes.  The fresh Mathematica worker
wall time includes worker startup and loading the prepared DAG and is bound to
a separately hashed `/usr/bin/time` log; its summed per-node timers do not.
The FLINT node-backend timer includes path-coordinate
evaluation, point serialization, the evaluator process, binary decoding, and
sparse matrix assembly, but excludes package/example setup.

## Fresh result

The all-Mathematica matrices were regenerated in a new `/tmp` directory; no
preserved matrix file was reused or overwritten.

| quantity | all Mathematica | FORM/FLINT run 1 | run 2 | run 3 |
|---|---:|---:|---:|---:|
| 96-node worker wall / complete node backend | 435.870 s | 39.492236 s | 41.106568 s | 36.154701 s |
| slowest worker's summed node timers | 427.247814 s | -- | -- | -- |
| FLINT evaluator process only | -- | 24.533456 s | 26.251546 s | 21.332926 s |
| FLINT interface including I/O and decode | -- | 25.685369 s | 27.413958 s | 22.508912 s |
| sparse matrix assembly | included above | 13.678793 s | 13.566877 s | 13.518641 s |
| fresh matrix import | 1.277295 s | -- | -- | -- |
| CHESSv2 propagation including node backend | -- | 49.205035 s | 50.767542 s | 45.860153 s |
| propagation after importing precomputed Mathematica nodes | 8.780763 s | -- | -- | -- |
| maximum absolute X2 error | 5.19528e-32 | 5.19528e-32 | 5.19528e-32 | 5.19528e-32 |

Using the actual eight-process wall time, the complete FLINT node backend is
10.60--12.06 times faster, with a median speedup of 11.04.  Using only the
slowest worker's summed node timers gives 10.39--11.82 times.  Comparing the
fresh Mathematica worker wall plus matrix import and propagation with the three
FLINT propagation calls gives an observed end-to-end improvement of
8.78--9.72 times, with median 9.06.

The historical result is therefore reproduced.  The earlier records reported
462.1 s for Mathematica node evaluation, 38.769--39.010 s for the experimental
FLINT node backend, and 48.536--48.918 s end to end.  The frozen CHESSv2 runs
give a median 39.49 s for the complete node backend and 49.21 s through the
unified `SpectralPropagate` entry.  The 36.15--41.11 s node range also exposes
the real process-load variability instead of presenting one favourable run.

## What the result means

The roughly ninefold end-to-end gain does not come from replacing the
non-canonical spectral solver.  After node matrices are available, the
Mathematica solve is only 8.79 s.  The gain comes from evaluating the 790 MiB
FORM-optimised straight-line program once for all nodes with FLINT/OpenMP;
Mathematica then handles only compact roots/one-forms and sparse assembly.

This benchmark also checks the new interface boundary: the FLINT module is
named `FLINTBatchEvaluation` and exposes `CHESSFLINT*`; `NativeBackend` is
reserved for the separate C++ spectral-propagation backend.  A polynomial
`CHESSNodeEvaluator` is sent through the unified entry by providing the
explicit positive `"EpsilonDegree" -> 2` option.  Every returned coefficient is
checked to be a numerical 131 by 131 matrix before it enters the solver.

Host: Intel Core i9-13950HX (24 cores, 32 logical CPUs), Mathematica 14.3.0,
g++ 11.4.0, FLINT 3.4.0, eight FLINT/OpenMP threads or eight Mathematica worker
processes.  The evaluator links to the hashed FLINT library listed below.

## Reproduction

FLINT route:

```bash
chess_f3_run_dir=$(mktemp -d /tmp/chessv2-f3-run.XXXXXX)

CHESSV2_F3_BENCHMARK_OUTPUT="${chess_f3_run_dir}/flint.wl" \
  wolframscript -file CHESSv2/Benchmarks/f3_zoia_flint_benchmark.wls
```

Fresh Mathematica node generation, using a new directory:

```bash
chess_f3_mma_dir="${chess_f3_run_dir}/mma-matrices"
chess_f3_mma_time="${chess_f3_run_dir}/mma-time.log"
mkdir -p "${chess_f3_mma_dir}"

/usr/bin/time -f 'CHESSV2_F3_MMA_WORKER_WALL_SECONDS=%e' \
  -o "${chess_f3_mma_time}" \
  sh -c "seq 8 | xargs -P 8 -I '{}' wolframscript -file \
  CHESSv2/Benchmarks/f3_zoia_mathematica_worker.wls \
  /home/liuyuanche/Program/Chebyshev_nonUT/ttW_physical_variables_2608_03746 \
  96 100 8 '{}' ${chess_f3_mma_dir}"

wolframscript -file \
  CHESSv2/Benchmarks/f3_zoia_mathematica_assemble.wls \
  /home/liuyuanche/Program/Chebyshev_nonUT/ttW_physical_variables_2608_03746 \
  96 100 "${chess_f3_mma_dir}" \
  "${chess_f3_run_dir}/mathematica.wl" 435.87 \
  "${chess_f3_mma_time}"
```

The wall-time argument should be replaced by the value in the current log.  The
assembler verifies that value against the log and records the log hash.  The
fresh 35 MiB matrix directory is intentionally retained until verification;
it may be removed after the raw record and time log have been archived.

## Provenance

External evaluator artifacts:

```text
43c702b11e8aac02254931e66a22795e41740302dcd1d888c8984f06ec5ff8cd  chess_slp_evaluator
6000b52692d0903ac4d833300b3d47371d96071c9f58a9bbffd66212ca79251d  chess_slp_evaluator.cpp
4915222391a69a8c708a505832b660ece833ef8af819128352f460aa178c9ca0  FLINT 3.4.0 libflint.so.22
d8b6c2f596a64544f696e364ed5f5eb7741435b8a420940b9ca02982a02f1511  F3_polynomials.chslp
08101639929f31874c576cb2147d0180864da1dda176de06e356d9c13550a51b  prepared_F3_backend.mx
81bee65753609d03479ce13455a4b7c0cd3792f9898a809f14f122e5a9d39a14  F3PhysicalBackend.wl
29b86f7b57c2bef35c1e914814258c76f2b11ef57ba0d5516c23b376f6c528af  F3_X1.m
44923eee9850055dcbee7ca658ea5a8e793f9abaa1f6e65b33d74bf32cfe285e  F3_X2.m
f422580f7a74d62d241c77c7a3158c3395a284679098f0c2a30bcf4c90255fde  prepared_F3_physical_DAG_X1_X2.mx
9931da8456ec886ed73e6cb0a0488a71e83bded6c45a4b1e34fba63b80db5d9e  f3_physical_dag_evaluator.wl
48294b347777fca190b062b53e2fe2235de9fe05660d9ab0957a008150e54035  legacy ChessNonCanonical.wl
```

Measured CHESSv2 sources:

```text
90016b6e47ec6258e6f198b5a9d9c4f1e6d271a23781c0208644c537fb14b96f  CHESSv2.wl
036a670f5983992dbe70c6b335ff6b498a0e1a868b8de9824b9918092c97d583  Core/Canonical.wl
03242f6b33c5cdc6cba0043393b0f806c151e4d4dc37601068ba21002929e3d9  Core/NonCanonical.wl
720fd96dcb777b61dd824a3c0d6f643d7a368eada00fb0c154d3efd864a81a4d  Core/Dispatch.wl
871d99f7abc3671fe4a8526138a863aa83b5123db5fda2d49bb3ea0d30bf117f  Core/FLINTBatchEvaluation.wl
e42be07dc117af951fb8d791212c65c3d6281cf415e853b647f894a8a1599629  Benchmarks/f3_zoia_flint_benchmark.wls
dfb1ce4a7e325775bca7481dce61b2fe2c105fe24c28aea43811a7285325a573  Benchmarks/f3_zoia_mathematica_worker.wls
f44ad3855c969ad5c299300383cf6dd7d72f2fae4488794b557f8f804230dbdd  Benchmarks/f3_zoia_mathematica_assemble.wls
```

Raw records were written outside the repository:

```text
9a4ddd044c199300fc74dd7fe22820f1e923a5843ebe842b6518e2311dc7091b  /tmp/chessv2-f3-zoia-flint-final-1.wl
477ae46bbc2ad1a17ec43709033a65ec9028edffb952b8390741b3e67ab6a885  /tmp/chessv2-f3-zoia-flint-final-2.wl
d0cac3c46507937506cf32f95529a671729628bfa83320ea4a8e52179a92f65a  /tmp/chessv2-f3-zoia-flint-final-3.wl
2a98d15f884abe4d000a02ed2e855868665f18a164e54847cd921b6febf37fca  /tmp/chessv2-f3-zoia-mathematica-final.wl
83bf975566405ba527b7f7164045e1f0c27776745b49e46d37b986d526416cc7  /tmp/chessv2-f3-mma-final-fbx3YM.time
```
