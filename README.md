# CHESS

CHESS is a high-precision Chebyshev--Lobatto spectral transport package for
systems of differential equations arising in Feynman-integral calculations.
The `experiment` branch contains CHESSv2, which unifies canonical,
non-canonical, and auxiliary-parameter propagation behind one Wolfram Language
entry point and provides optional C++ and FORM/FLINT numerical backends.

The original release package remains available as `Chess.wl`. New development,
tests, and backend integration live under `CHESSv2/`.

## Highlights

- one public `SpectralPropagate` function for canonical and polynomial-
  epsilon systems;
- automatic routing for pure `B0`, pure `B1`, and mixed
  `B0 + epsilon B1 + ...` equations;
- endpoint regularization through the verified Mathematica cores;
- automatic auxiliary-delta order selection for regular `B0`-only systems;
- arbitrary-precision C++ transport using MPFR/MPC and OpenMP;
- direct FLINT `nfloat_complex` factorization for regular mixed active blocks;
- FORM/FLINT multipoint evaluation for large generated expressions;
- reusable Native handles for repeated boundary transports;
- fail-closed dimension, precision, protocol, and convergence checks.

## Requirements

The current package and test suite are validated with Mathematica 14.3 through
WolframScript 1.13.0. The optional Native backend additionally requires a C++17
compiler, WSTP, Boost.Multiprecision, MPFR, MPC, GMP, FLINT, and OpenMP. Build
paths are configurable in `CHESSv2/Native/Propagation/Makefile`.

## Quick start

Load CHESSv2 directly:

```wl
Get["/path/to/CHESS/CHESSv2/CHESSv2.wl"];
```

No context prefix is required:

```wl
result = SpectralPropagate[
  {B0, B1, B2},
  boundary,
  {0, 1},
  "Nodes" -> 48,
  "Precision" -> 100,
  "WorkingPrecisionA" -> 140
];
```

The equation convention is

```text
dY(t,epsilon)/dt =
  (B0(t) + epsilon B1(t) + epsilon^2 B2(t) + ...) Y(t,epsilon).
```

| Matrix specification | Route |
|---|---|
| `B1` or `{0,B1}` | canonical sequential-epsilon solver |
| `{B0}` | auxiliary-delta recursion |
| `{B0,B1,...}` | active-support non-canonical solver |

A bare evaluator is interpreted as `B1` for compatibility with canonical CHESS
notebooks. Routing is structural: CHESSv2 never samples a function and guesses
that it is zero.

## Matrix evaluators

Coefficient functions may use any established CHESS calling convention:

```wl
Bp[t_, "Precision" -> p_] := ...
Bp[t_, p_] := ...
Bp[t_] := ...
```

For expensive multipoint evaluation, provide both scalar and batch functions:

```wl
CHESSNodeEvaluator[scalarEvaluator, batchEvaluator]
```

The batch signature is
`batchEvaluator[nodes, precision, threads]`. If each batch item is a complete
coefficient list `{B0,B1,...}`, also set a positive `"EpsilonDegree"`.

## Optional Native backend

Build the propagation executable:

```bash
make -C CHESSv2/Native/Propagation
```

Enable it on a supported regular path:

```wl
nativeResult = SpectralPropagate[
  {B0}, boundary, {0, 1},
  "Nodes" -> 96,
  "Precision" -> 100,
  "WorkingPrecisionA" -> 140,
  "NumericalBackend" -> "Native",
  "ResultData" -> "Endpoint"
];
```

Native arithmetic uses 20 guard digits by default. Thus `"Precision" -> p`
keeps the public result at `p` digits while the base C++ transport uses
`p + 20`. Set `"NativeGuardDigits" -> 0` to reproduce the historical
no-guard behavior, or increase it for a difficult system. Low-precision
nonzero matrices, boundaries, and interval endpoints are rejected rather than
padded to a misleading arbitrary precision.

For repeated transports, prepare and reuse a `"NativeHandle"`; setup includes
matrix sampling, serialization, and LU factorization and can dominate a
one-shot calculation.

## FORM/FLINT expression evaluation

`CHESSFlintRunSLP` evaluates a prepared FORM straight-line program at all path
nodes in one FLINT process. `CHESSNativeFlintMatrixAdapter` combines that batch
evaluator with problem-specific coordinate and matrix-assembly functions.
Generated process-specific SLPs are intentionally not embedded in the package.

## Supported Native routes

- one regular canonical interval;
- regular `B0`-only propagation, including regular segmentation;
- one regular mixed polynomial batch with nonempty `B0` active support;
- multiple physical boundary columns;
- full node data or endpoint-only output.

Singular endpoint regularization and unsupported mixed/segmented combinations
remain on the Mathematica route or fail closed.

## Repository layout

```text
Chess.wl                       original release package
CHESSv2/CHESSv2.wl             modular loader
CHESSv2/Core/                  Wolfram routing and numerical modules
CHESSv2/Native/Propagation/    C++/MPFR/MPC/FLINT backend
CHESSv2/FLINT/                 FORM/FLINT protocol documentation
CHESSv2/Tests/                 regression and backend tests
CHESSv2/Benchmarks/            reproducible performance studies
```

## Tests

```bash
wolframscript -file CHESSv2/Tests/run_all.wls
wolframscript -file CHESSv2/Tests/test_native_backend.wls
wolframscript -file CHESSv2/Tests/test_flint_evaluation.wls
```

The Native suite requires a compiled backend. The real F3 FLINT smoke test is
enabled with `CHESSV2_FLINT_EVALUATOR` and `CHESSV2_FLINT_SLP`.

## Documentation

- [CHESSv2 reference and usage](CHESSv2/README.md)
- [FORM/FLINT batch interface](CHESSv2/FLINT/README.md)
- [C++ propagation backend](CHESSv2/Native/Propagation/README.md)
- [PBB64 Native benchmark](CHESSv2/Benchmarks/PBB64_NATIVE_BENCHMARK.md)
- [Simone Zoia F3 benchmark](CHESSv2/Benchmarks/F3_ZOIA_FLINT_BENCHMARK.md)

Benchmark reports state their timing boundaries explicitly. They measure
matrix evaluation or differential-equation transport, not complete IBP or
boundary-generation workflows.
