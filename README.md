# CHESS

CHESS is a Wolfram Language package for Chebyshev pseudo-spectral transport of
canonical Feynman-integral differential equations along a one-dimensional path.

## Contents

- `Chess.wl` - the package.
- `examples/run_example.wls` - runnable examples.
- `data/` - matrices, prepared letters, boundary values, and reference values
  needed by the examples.

No benchmark logs, paper-generation scripts, or development files are included
in this release archive.

## Requirements

- Wolfram Mathematica / WolframScript 13.0 or newer.

The bundled examples have been smoke-tested with Mathematica 13.0.1 and 14.3.

## Run

From this directory:

```bash
wolframscript -file examples/run_example.wls DP
```

The command above runs the DP example with 24 Chebyshev nodes, one kernel,
`"Precision" -> 80`, and `"WorkingPrecisionA" -> 100`.

Arguments are:

```bash
wolframscript -file examples/run_example.wls CASE NODES KERNELS PRECISION PRECISION_A
```

Available cases:

- `DP` - no endpoint regularization.
- `HB` - right-endpoint regularization.
- `BHB` - left-endpoint regularization.

For example:

```bash
wolframscript -file examples/run_example.wls HB 24 2 80 100
wolframscript -file examples/run_example.wls BHB 24 2 80 100
```

The script prints the matrix dimensions, endpoint mode, runtime, endpoint
diagnostics, a final-state checksum, and a reference error when a reference
value is included.
