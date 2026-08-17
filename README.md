# CHESS experiment branch

This branch contains CHESSv2, the modular development version of CHESS for
high-precision Chebyshev spectral transport of canonical and non-canonical
Feynman-integral differential equations.

The old release examples and smoke-test runner have been removed from this
branch. Current regression tests and reproducible performance studies live
inside `CHESSv2/Tests` and `CHESSv2/Benchmarks`.

## Quick start

```wl
Get["/path/to/CHESS/CHESSv2/CHESSv2.wl"];

result = SpectralPropagate[
  {B0, B1, B2},
  boundary,
  {0, 1},
  "Nodes" -> 48,
  "Precision" -> 100,
  "WorkingPrecisionA" -> 120
];
```

`SpectralPropagate` is loaded directly into the user context; callers do not
need a `CHESS`` or `CHESSv2`` prefix.

## Optional numerical backends

CHESSv2 keeps symbolic routing in Wolfram Language and moves expensive
floating-point work to optional native modules:

- FORM and FLINT batch evaluation for large generated expressions;
- C++/MPFR/MPC sequential spectral propagation;
- FLINT `nfloat_complex` factorization of regular mixed-system active blocks.

Build the propagation backend with:

```bash
make -C CHESSv2/Native/Propagation
```

The default backend remains Mathematica. Native execution is enabled with
`"NumericalBackend" -> "Native"` on supported regular paths. Direct mixed
native transport is intentionally batch-only: it requires
`CHESSNodeEvaluator[scalar,batch]` and an explicit positive `"EpsilonDegree"`.
Ordinary mixed coefficient lists remain on the Mathematica route.

## Documentation and validation

- [CHESSv2 documentation](CHESSv2/README.md)
- [FORM/FLINT batch interface](CHESSv2/FLINT/README.md)
- [native propagation backend](CHESSv2/Native/Propagation/README.md)

Run the current tests with:

```bash
wolframscript -file CHESSv2/Tests/run_all.wls
wolframscript -file CHESSv2/Tests/test_native_backend.wls
wolframscript -file CHESSv2/Tests/test_flint_evaluation.wls
```

The native tests require a compiled backend. The real F3 FLINT smoke test also
requires the external evaluator and prepared SLP paths described in the test
file.
