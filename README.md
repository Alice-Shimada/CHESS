# CHESS

CHESS is a Wolfram Language package for Chebyshev pseudo-spectral transport of
canonical Feynman-integral differential equations along one-dimensional
kinematic paths.

## Contents

- `Chess.wl` - the package.
- `examples/run_example.wls` - runnable examples matching the paper benchmarks.
- `examples/run_tests.wls` - lightweight pass/fail smoke tests.
- `examples/BHB_family_example.nb` - a detailed notebook for the 3L5P BHB
  workflow.
- `examples/expected/` - sample outputs for the smoke-test settings.
- `data/` - matrices, prepared one-dimensional letters, boundary values, and
  references needed by the examples.

No benchmark logs, paper source files, development scripts, or temporary files
are included in this release archive.

## Requirements

- Wolfram Mathematica / WolframScript 13.0 or newer.

The bundled examples have been smoke-tested with Mathematica 14.3.  Earlier
13.x versions should work because the package uses standard Wolfram Language
functions, and the previous release smoke tests were checked on 13.0.1.

## Examples

From this directory:

```bash
wolframscript -file examples/run_example.wls DP
```

Arguments are:

```bash
wolframscript -file examples/run_example.wls CASE NODES KERNELS PRECISION PRECISION_A
```

Available cases:

- `DP` - two-loop six-point DP physical-region transport, no endpoint
  regularization, checked against the bundled AMFlow value.
- `PBB`, `BPB`, `BHB`, `PBP` - three-loop five-point planar families, using
  left-endpoint regularization at the maximally symmetric Euclidean point.
- `BHABHA` - two-loop massive Bhabha-scattering planar family using a direct
  pulled-back matrix evaluator, checked against bundled AMFlow endpoint values.

Small smoke-style runs:

```bash
wolframscript -file examples/run_example.wls DP 8 1 60 80
wolframscript -file examples/run_example.wls PBB 8 1 60 80
wolframscript -file examples/run_example.wls BPB 8 1 60 80
wolframscript -file examples/run_example.wls BHB 8 1 60 80
wolframscript -file examples/run_example.wls PBP 8 1 60 80
wolframscript -file examples/run_example.wls BHABHA 24 1 80 100
```

Paper-scale runs use larger node counts and more kernels.  For example, the
paper uses 72, 96, and 120 nodes for the DP physical-region convergence table,
96 and 120 nodes for the 3L5P node-convergence comparison, and 96 nodes for the
Bhabha direct-matrix smoke comparison.  The Bhabha high-precision table in the
paper also reports larger node counts up to 512.

The runner prints matrix dimensions, endpoint mode, runtime, endpoint
diagnostics, a final-state checksum, and a reference error when a reference is
bundled.

## Smoke Tests

Run all lightweight release checks with:

```bash
wolframscript -file examples/run_tests.wls
```

A successful run terminates with exit code 0 and prints
`CHESS_RELEASE_TESTS_PASSED`.

The smoke tests use one kernel and low node counts.  They verify numerical
agreement for the DP and BHABHA examples and successful endpoint-regularized
completion for the four 3L5P families.

## Input Convention

For a standard dlog problem, provide:

- `Atilde`: a matrix whose entries are linear combinations of `Log[W[i]]` or
  `logW[i]`;
- `dLettersLine`: a list where `dLettersLine[[i]]` is
  `d log(W_i(x(t))) / dt`;
- `boundaryValues`: an `nIntegrals x nEpsilonOrders` matrix of epsilon
  coefficients at the left endpoint;
- a path interval, usually `{0, 1}`.

The package core supports the unified `W`/`logW` convention.  Case-specific
letter heads, such as legacy `Wtilde` or `What`, are converted in the examples
before calling `CHESSAtildeLinearData[]`.

The Bhabha example shows the more general interface: users may supply any
matrix-valued function `nAfun[t]`, not necessarily one assembled from prepared
dlog letters.  In that example the runner directly evaluates
`(3/10) A_s(s(x),t(x),1) - (7/10) A_t(s(x),t(x),1)` along the path
`s(x)=3 x/10`, `t(x)=-7 x/10`.
