# C++ spectral-propagation backend

This optional backend moves the regular-point floating-point hot loop out of
Mathematica. It caches one collocation problem, then accepts one or more
boundary columns through WSTP.

The kernel uses:

- Boost.Multiprecision `mpfr_float` / `mpc_complex` for arbitrary-precision
  real and complex arithmetic;
- MPFR, MPC, and GMP as the numerical libraries;
- FLINT `nfloat_complex` for the directly factorized mixed-system active block;
- OpenMP across independent components and node source products;
- one cached LU factorization of the scalar Lobatto block;
- sparse node operators;
- a real-only fast path for the sequential canonical/fake-delta kernel.

The separate `FLINTBatchEvaluation` module remains the public interface for
large FORM-generated expressions. NativeBackend uses FLINT internally only for
the direct dense LU of a small active-support collocation block; these are
distinct capabilities and APIs.

## Build

The defaults match the development machine and can be overridden:

```bash
make \
  BOOST_INCLUDE=/path/to/boost_1_83_0 \
  WSTP_DIR=/path/to/Mathematica/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/CompilerAdditions \
  FLINT_INCLUDE=/path/to/flint/include \
  FLINT_DIR=/path/to/flint/lib
```

The generated `chess_native_link`, object directory, and `wsprep` output are
ignored by Git. A package checkout therefore never ships an opaque executable.

Source ownership is intentionally narrow: `backend_types.hpp` contains shared
numeric and sparse data types, `flint_coupled_solver.*` owns every FLINT context/conversion/LU,
and `backend.cpp` owns only the WSTP protocol and transport orchestration.

## Supported routes

- one regular canonical interval;
- regular B0-only propagation through the auxiliary-delta recurrence;
- one regular mixed polynomial `CHESSNodeEvaluator` batch interval through a
  directly factorized active-support operator, without an auxiliary-delta
  series;
- multiple boundary columns in one call;
- full node data or endpoint-only transfer;
- fixed or automatically estimated delta order (selection remains in the
  Mathematica interface).

Segmented B0 propagation prepares one handle per segment. Direct mixed
transport currently requires one regular interval, an explicit positive
`EpsilonDegree`, and a `CHESSNodeEvaluator` batch function. Ordinary
coefficient lists remain on the Mathematica route; an explicit Native request
for such a list fails. Singular endpoint regularization also remains in
Mathematica.

## Performance contract

`CHESSNativePrepare` samples matrices, serializes them, and factors the Lobatto
block once. Reuse its handle through the public `"NativeHandle"` option when
transporting several boundaries. `"ResultData" -> "Endpoint"` avoids caching
and transferring all node states and is the intended high-throughput mode.
The preparation signatures are

```wl
CHESSNativePrepare[matrixSpec, dimension, interval, optionRules]
CHESSNativePolynomialPrepare[evaluator, boundary, interval, optionRules]
```

Pass the returned Association back through `"NativeHandle"`; the mixed form
also requires a positive `"EpsilonDegree"` in `optionRules`.

The public `"Precision"` remains the requested output precision; Native setup
and propagation use `"NativeGuardDigits" -> 20` extra decimal digits by
default. The sampled matrix precision is raised when necessary so it is never
requested below the guarded Native precision plus ten sampling digits; returned
nonzero entries must retain at least the guarded Native precision. Handles
record both values and cannot be reused with a different guard setting.
Nonzero machine-precision matrix, interval, or boundary input fails closed
instead of being padded to the guarded precision. The direct mixed FLINT LU
keeps its existing additional internal margin above the sampled matrix
precision.
`CHESSNativePolynomialPrepare` additionally constructs the active component
list. C++ then assembles and factorizes the Kronecker product of `D` with the
active-space identity minus the block-diagonal matrix whose blocks are
`B0(t_j)[[active,active]]` at ordinary, non-boundary Lobatto nodes; the boundary
block is unshifted. Positive epsilon powers remain sparse RHS operators. This
preparation is substantial, so the direct mixed route is most useful when the
handle is reused.
The handle is bound to the evaluator expression, the transitive definitions of
its Wolfram helper functions, and the numerical setup. Ordinary helper
redefinitions invalidate it automatically. Mutable external state which is not
visible in Wolfram definitions still requires preparing a new handle.
