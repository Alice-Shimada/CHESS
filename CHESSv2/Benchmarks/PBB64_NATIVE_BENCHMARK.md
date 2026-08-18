# PBB64 native-versus-Mathematica propagation benchmark

Date: 2026-08-17

## Timing boundary

This is a numerical-transport benchmark, not a full Feynman-integral or IBP
benchmark. Both routes propagate the same public 316-dimensional PBB equation
on the regular interval `{-17/10,-19/10}` at the physical parameter
`2275901/2275800`.

Data loading, rational path construction, and exact `Atilde` decomposition are
outside both timed calls. For the reused native handle only, process startup,
81-node matrix sampling, setup serialization, and native LU factorization are
also outside and reported separately. The complete Mathematica public call
still constructs its node operators; the native public call still includes
handle-definition fingerprint validation, boundary formatting, WSTP transfer,
kernel work, and endpoint parsing. Thus the
comparison is the repeated-call API cost available to a user, not two kernels
with identical preprocessing removed. Both routes use 64 decimal digits,
matrix precision 90, 82 Lobatto points (`"Nodes" -> 81`), and fixed auxiliary-
delta order 56. The native route uses eight OpenMP threads and
`"ResultData" -> "Endpoint"`.

## Result

| quantity | value |
|---|---:|
| Mathematica propagation times | 19.741337, 21.170378, 21.880414, 21.871567, 22.959160 s |
| Native propagation times | 2.591288, 2.279936, 2.271004, 2.307695, 2.274918 s |
| Mathematica median | 21.871567 s |
| Native median | 2.279936 s |
| Last native kernel time | 1.398000036 s |
| Native preparation | 17.283158 s |
| of which native load/factorization | 0.698352 s |
| Cached propagation speedup | 9.593062x |
| Mathematica warm median / (prepare + native warm median) | 1.118001x |
| Observed Mathematica cold / (prepare + native cold) | 1.005748x |
| Scaled endpoint difference | 3.99402e-57 |

The cached C++ route is therefore about 9.59 times faster for this repeated
propagation task. Setup is substantial, so the supported conclusion is about a
reused handle. The warm-median setup ratio is an amortization estimate, not an
observed cold one-shot pair. The separately reported cold ratio uses the actual
two cold calls. Both happen to be slightly above one on this large task;
neither crossover should be generalized to smaller systems.

The difference is measured as

```text
max_i |native_i - Mathematica_i| / (1 + |Mathematica_i|).
```

## Reproduction

```bash
make -C CHESSv2/Native/Propagation

OMP_NUM_THREADS=8 CHESSV2_BENCHMARK_REPEATS=5 \
  wolframscript -file CHESSv2/Benchmarks/pbb64_native_vs_mathematica.wls
```

The harness writes the complete endpoints to
`/tmp/chessv2-pbb64-native-vs-mma.wl` by default; the large generated result is
not committed.

Host: Intel Core i9-13950HX, Mathematica 14.3, g++ 11.4.0, Boost 1.83,
MPFR/MPC/GMP, eight OpenMP threads. Public PBB data commit:
`3e5862a2c28bd1bd080d141894db85a010bce2c1`.

Source hashes for the measured run:

```text
8060a1062e857c9fa2f3bbea79d375ffbb819ee32d5d342b9ad66c3fadbd9e0f  Native/Propagation/src/backend.cpp
445438d9947c299584533739cb2e80eb84b894cbc497526de3631846e6f75a21  Core/NativeBackend.wl
00724534ef906da1d70b2f16e68077156c14361e2e3e9a08dec7e474804f6627  Benchmarks/pbb64_native_vs_mathematica.wls
```
