# C++ spectral-propagation backend

This optional backend moves the regular-point floating-point hot loop out of
Mathematica. It caches one collocation problem, then accepts one or more
boundary columns through WSTP.

The kernel uses:

- Boost.Multiprecision `mpfr_float` / `mpc_complex` for arbitrary-precision
  real and complex arithmetic;
- MPFR, MPC, and GMP as the numerical libraries;
- OpenMP across independent components and node source products;
- one cached LU factorization of the scalar Lobatto block;
- sparse node operators;
- a real-only fast path when both operators and boundary data are real.

FLINT is not used for this linear-algebra kernel. The separate
`Native/Evaluation` module uses FLINT `nfloat` where it is advantageous: large
FORM-generated straight-line expressions evaluated at many nodes.

## Build

The defaults match the development machine and can be overridden:

```bash
make \
  BOOST_INCLUDE=/path/to/boost_1_83_0 \
  WSTP_DIR=/path/to/Mathematica/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/CompilerAdditions
```

The generated `chess_native_link`, object directory, and `wsprep` output are
ignored by Git. A package checkout therefore never ships an opaque executable.

## Supported routes

- one regular canonical interval;
- regular B0-only propagation through the auxiliary-delta recurrence;
- multiple boundary columns in one call;
- full node data or endpoint-only transfer;
- fixed or automatically estimated delta order (selection remains in the
  Mathematica interface).

Segmented B0 propagation prepares one handle per segment. Singular endpoint
regularization and genuinely mixed polynomial-epsilon systems remain on the
existing Mathematica algorithms. Unsupported native requests fail explicitly.

## Performance contract

`CHESSNativePrepare` samples matrices, serializes them, and factors the Lobatto
block once. Reuse its handle through the public `"NativeHandle"` option when
transporting several boundaries. `"ResultData" -> "Endpoint"` avoids caching
and transferring all node states and is the intended high-throughput mode.
The handle is bound to the evaluator expression, the transitive definitions of
its Wolfram helper functions, and the numerical setup. Ordinary helper
redefinitions invalidate it automatically. Mutable external state which is not
visible in Wolfram definitions still requires preparing a new handle.
