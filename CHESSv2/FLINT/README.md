# FORM/FLINT batch-expression interface

This directory documents the optional expression-evaluation side of CHESSv2.
It is deliberately separate from `Native/Propagation`: a problem-specific
FORM straight-line program may be evaluated by FLINT while spectral transport
still uses Mathematica, or the same evaluator may feed the C++ transport
backend.

The package function

```wl
CHESSFLINTRunSLP[executable, slpFile, pointRows, precision, threads]
```

implements the FLINT batch protocol already exercised by the experimental
non-canonical phenomenology backend:

```text
executable evaluate SLP POINTS BITS THREADS OUTPUT
```

`pointRows` is one real coordinate row per path node. The FLINT evaluator writes
the little-endian `NFLOAT01` binary format; the result Association contains its
arbitrary-precision values under `"Values"`.

The default `"GuardDigits" -> 20` is only a working-precision heuristic. It
cannot bound cancellation, proximity to a pole, or intermediate expression
growth. Callers can raise `"GuardDigits"`, set `"WorkingBits"` explicitly, and
should repeat sensitive evaluations at higher precision. Non-finite FLINT
records (`+Infinity`, `-Infinity`, or `NaN`) and exponents outside FLINT's
documented finite `nfloat` interval are rejected rather than converted into
matrix entries. `"Timeout" -> seconds` places a finite bound on the
external call; its default is `Infinity` because production SLP sizes vary by
orders of magnitude.

The FORM program, its input-variable names, and matrix assembly remain local to
the physics example. A typical adapter is:

```wl
batch[nodes_, precision_, threads_] := Module[{coordinates, flintResult},
  coordinates = pathCoordinates /@ nodes;
  flintResult = CHESSFLINTRunSLP[
    evaluatorExecutable, preparedSLP, coordinates, precision, threads
  ];
  If[flintResult === $Failed,
    $Failed,
    MapThread[
      assembleCoefficientMatrices,
      {nodes, flintResult["Values"]}
    ]
  ]
];

flintEvaluator = CHESSNodeEvaluator[scalarFallback, batch];

result = SpectralPropagate[
  flintEvaluator, boundary, {0, 1},
  "EpsilonDegree" -> 2
];
```

The explicit positive epsilon degree tells the unified entry that each batch
item is a polynomial matrix list rather than one canonical matrix.

No compiled process-specific expression is shipped in the package. This keeps
large generated FORM code and its variable convention out of CHESS, while the
stable multipoint transport interface remains reusable.
