# CHESSv2

CHESSv2 unifies the original canonical CHESS implementation and the
`CHESS_non_canonical` solver behind one public function:

```wl
SpectralPropagate[matrixSpec, boundary, {t0, t1}, options]
```

The default numerical backend is Mathematica. Optional FORM/FLINT and C++
modules accelerate large expression evaluation and regular spectral transport
without changing the public entry point.

## Loading

```wl
Get["/path/to/CHESS/CHESSv2/CHESSv2.wl"];
```

Use `SpectralPropagate` directly. The experimental version string is stored in
`$CHESSv2Version`.

Function names use PascalCase while preserving `CHESS` as an uppercase project
acronym and spelling the library name as `Flint`, for example
`CHESSFlintRunSLP`. The former `CHESSFLINT*` spellings are intentionally not
kept as compatibility aliases on the experimental branch.

## Equation convention and routing

CHESSv2 solves

```text
dY(t,eps)/dt = (B0(t) + eps B1(t) + eps^2 B2(t) + ...) Y(t,eps).
```

The item at position `p+1` in an explicit coefficient list is `Bp`:

```wl
SpectralPropagate[{B0, B1, B2}, boundary, {t0, t1}, options]
```

| Non-zero coefficients | Algorithm |
|---|---|
| `B0` only | auxiliary-delta recursion |
| `B1` only | original canonical sequential-epsilon solver |
| mixed `B0,B1,...` or higher powers | original active-support non-canonical solver |

Only literal zero and exact zero constant matrices are classified as zero.
CHESSv2 does not sample a function and guess that it vanishes.

For compatibility with canonical notebooks, a bare evaluator is interpreted as
`B1`:

```wl
SpectralPropagate[B1, boundary, {t0, t1}, options]
```

## Matrix evaluators

Each coefficient may be a constant numerical matrix or a function using one of
these conventions:

```wl
Bp[t_, "Precision" -> p_] := ...
Bp[t_, p_] := ...
Bp[t_] := ...
```

Every result must be a numerical square matrix matching the number of rows in
`boundary`.

For expensive generated expressions, provide a scalar and a batch evaluator:

```wl
CHESSNodeEvaluator[scalarEvaluator, batchEvaluator]
```

The batch function has the form `batchEvaluator[nodes, precision, threads]`.
If one batch item is the complete list `{B0,B1,...}`, give an explicit positive
epsilon degree:

```wl
SpectralPropagate[
  CHESSNodeEvaluator[scalarFallback, polynomialBatch],
  boundary,
  {0, 1},
  "EpsilonDegree" -> 2
]
```

## Native matrix adapters

`Core/NativeMatrixAdapter.wl` provides two independent helpers. Both return a
`CHESSNodeEvaluator` accepted by the same `SpectralPropagate` call.

For an existing Wolfram coefficient list, use:

```wl
nativeEvaluator = CHESSNativeMatrixAdapter[{B0, B1, B2}, dimension];

SpectralPropagate[
  nativeEvaluator, boundary, {0, 1},
  "EpsilonDegree" -> 2,
  "NumericalBackend" -> "Native"
]
```

The list entries may be constant matrices or any evaluator convention described
above. This form accelerates transport, but the node matrices are still
evaluated in Mathematica. If the matrices are explicit expressions in one path
variable, first `Clear[t]`, then use
`CHESSNativeMatrixAdapter[{B0expr,B1expr,...}, t, dimension]`.
The built-in Wolfram batch wrapper maps nodes serially; a user evaluator may
implement its own parallelism when useful.

To move large node expression evaluation to FORM/FLINT, use:

```wl
nativeEvaluator = CHESSNativeFlintMatrixAdapter[
  evaluatorExecutable,
  preparedSLP,
  coordinateFunction,
  assembleFunction,
  dimension
];
```

`coordinateFunction[point,precision]` returns one real SLP input row.
`assembleFunction[point,values,precision]` converts one FLINT output row into
`{B0,B1,...}`. The helper batches all requested nodes in one FLINT process.
If the executable or SLP file is replaced at the same path, prepare a new
Native handle; file contents are external state and are not part of the Wolfram
definition fingerprint.

## FORM/FLINT expression evaluation

`Core/FLINTBatchEvaluation.wl` exposes:

```wl
CHESSFlintRunSLP[executable, slpFile, pointRows, precision, threads]
```

FORM may optimize large exact expressions into a straight-line program. FLINT
then evaluates every requested point with `nfloat` arithmetic and returns the
lossless `NFLOAT01` binary encoding of those finite-precision values.

The default 20 guard digits are a heuristic, not a rigorous error bound. Near a
pole or under severe cancellation, increase `"GuardDigits"` or set
`"WorkingBits"` explicitly and repeat the calculation at higher precision.
Non-finite values and invalid headers, signs, exponents, or zero mantissas fail
closed.

See [FLINT/README.md](FLINT/README.md) for the protocol and adapter example.

## Native spectral propagation

Build the optional backend with:

```bash
make -C CHESSv2/Native/Propagation
```

The development Makefile uses WSTP, Boost.Multiprecision, MPFR/MPC/GMP, FLINT,
and OpenMP. Override these paths when necessary:

```bash
make -C CHESSv2/Native/Propagation \
  BOOST_INCLUDE=/path/to/boost \
  WSTP_DIR=/path/to/WSTP/CompilerAdditions \
  FLINT_INCLUDE=/path/to/flint/include \
  FLINT_DIR=/path/to/flint/lib
```

Enable it through the unified entry:

```wl
result = SpectralPropagate[
  {B0}, boundary, {0, 1},
  "Nodes" -> 48,
  "Precision" -> 100,
  "WorkingPrecisionA" -> 120,
  "NumericalBackend" -> "Native",
  "ResultData" -> "Endpoint"
];
```

Supported native routes are:

- regular canonical propagation;
- regular `B0`-only auxiliary-delta propagation;
- one regular mixed polynomial `CHESSNodeEvaluator` batch with an explicit
  positive `"EpsilonDegree"` and non-empty `B0` active support.

The mixed route directly constructs
`D tensor I_active - diag(B0)` and factorizes it with FLINT
`nfloat_complex`; it does not use an auxiliary-delta series. Remaining inactive
components reuse the scalar Lobatto LU.

This direct mixed native capability is batch-only. An explicit Native request
for an ordinary `{B0,B1,...}` coefficient list fails rather than silently
changing algorithms.

Prepared handles may be reused through `"NativeHandle"`. A handle is bound to
the evaluator definitions, interval, dimension, node count, requested
precision, and Native guard precision. Ordinary Wolfram definition changes
invalidate it automatically.

Native propagation uses 20 guard digits by default: `"Precision" -> p`
returns a result at the requested precision while the base Native transport,
Lobatto matrix, boundary serialization, and result parsing use `p + 20`
digits. Matrix sampling uses the larger of `"WorkingPrecisionA"` and this
Native working precision plus ten sampling digits; the direct mixed FLINT LU
retains an additional internal margin above that sampling precision.
Set `"NativeGuardDigits" -> 0` to recover the historical no-guard behavior, or
increase it for an ill-conditioned problem. This is a working-precision
heuristic rather than a rigorous error bound, so precision-sensitive results
should still be repeated at a higher setting.

Native preparation fails closed if a nonzero matrix entry, interval endpoint,
or boundary value carries insufficient actual precision. In particular, a
machine-precision evaluator is never padded with decimal zeros and presented as
an arbitrary-precision result.

| Option | Default | Meaning |
|---|---:|---|
| `"NumericalBackend"` | `"Mathematica"` | select `"Native"` explicitly |
| `"NativeHandle"` | `Automatic` | reuse a prepared native system |
| `"NativeGuardDigits"` | `20` | extra C++/FLINT working digits |
| `"ResultData"` | `"Full"` | use `"Endpoint"` to avoid node-state transfer |

## Auxiliary delta for `B0`-only systems

For

```text
Y'(t) = B0(t) Y(t),
```

CHESSv2 introduces `delta` and expands

```text
Y(t,delta) = C0(t) + delta C1(t) + delta^2 C2(t) + ...,
C0'(t) = 0,
Ck'(t) = B0(t) C(k-1)(t).
```

After propagation, the coefficients are summed at `delta=1`. Automatic order
selection uses the last three consecutive coefficient ratios and one guard
layer. This is an empirical stopping diagnostic, not a strict bound on the
unseen tail. Important results should be repeated with a tighter tolerance or a
higher fixed order.

| Option | Default | Meaning |
|---|---:|---|
| `"DeltaOrder"` | `Automatic` | automatic selection or fixed non-negative order |
| `"DeltaTolerance"` | `Automatic` | target empirical tail size |
| `"MaxDeltaOrder"` | `128` | maximum automatic search order |
| `"Segments"` | `1` | regular-path segmentation |

Failure to satisfy the criterion by `"MaxDeltaOrder"` returns `$Failed`.

## Result structure

All routes return the historical six-part form:

```wl
{nodes, nodeValues, finalValue, endpointBlocks, endpointInfo, methodMetadata}
```

`finalValue` is the transported endpoint coefficient matrix. On a native route,
`"ResultData" -> "Endpoint"` makes `nodeValues`
`Missing["NotRequested"]`. Mathematica routes retain their historical full
node data.

## Current limitations

- `B0`-only auxiliary-delta propagation currently requires regular endpoints.
- Direct mixed native propagation requires one regular interval, an explicit
  polynomial batch evaluator, and non-empty `B0` active support.
- Singular endpoint regularization remains in the Mathematica solvers.
- Segmented mixed propagation with repeated endpoint constraints fails closed.
- Invalid dimensions, malformed backend data, unsupported option combinations,
  and unconverged automatic delta searches return `$Failed`.

## Source layout

```text
CHESSv2.wl                    loader
Core/Canonical.wl             canonical numerical core
Core/NonCanonical.wl          polynomial-epsilon numerical core
Core/MatrixAdapters.wl        evaluator normalization and dimension checks
Core/FakeDelta.wl             B0-only auxiliary-delta logic
Core/FLINTBatchEvaluation.wl  FORM/FLINT batch protocol
Core/NativeMatrixAdapter.wl   Wolfram and FORM/FLINT matrix input helpers
Core/NativeBackend.wl         Wolfram/native interface and handle management
Core/Dispatch.wl              unified SpectralPropagate routing
Native/Propagation/           C++/FLINT propagation backend
  src/backend.cpp             WSTP protocol and transport orchestration
  src/flint_coupled_solver.*  FLINT contexts, conversion, direct LU and solve
  src/backend_types.hpp       shared numerical POD types
FLINT/                        expression-evaluator documentation
Tests/                        current regression tests
Benchmarks/                   reproducible performance studies
```

## Tests

```bash
wolframscript -file CHESSv2/Tests/run_all.wls
wolframscript -file CHESSv2/Tests/test_native_backend.wls
wolframscript -file CHESSv2/Tests/test_flint_evaluation.wls
```

`test_native_backend.wls` requires a compiled backend. The optional real F3 SLP
smoke test is enabled by setting `CHESSV2_FLINT_EVALUATOR` and
`CHESSV2_FLINT_SLP`.

## Benchmarks

- `Benchmarks/PBB64_NATIVE_BENCHMARK.md`: 316-dimensional, 64-digit regular
  propagation with a reused native handle.
- `Benchmarks/F3_ZOIA_FLINT_BENCHMARK.md`: 131-dimensional F3 expression
  evaluation and non-canonical propagation at 100 digits.

Each report states its timing boundary. Do not interpret these transport and
evaluation measurements as full IBP-workflow timings.
