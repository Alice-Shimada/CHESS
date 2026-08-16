#!/usr/bin/env wolframscript

(*
  Fresh-kernel regression suite for the first minimal CHESSv2 merge.

  The tests are intentionally small enough to run on every edit.  They check
  routing and mathematical equivalence, not production-scale performance.
*)

ClearAll["Global`*"];
$HistoryLength = 0;
$MaxExtraPrecision = 3000;

testRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{testRoot, "CHESSv2.wl"}]];

failures = {};
assertionCount = 0;
expectedAssertionCount = 23;
assert[name_String, condition_] := If[
  assertionCount++;
  TrueQ[condition],
  Print["PASS: ", name],
  AppendTo[failures, name];
  Print["FAIL: ", name]
];

maxError[a_, b_] := Max[Abs[Flatten[N[a - b, 50]]]];
metadataValue[metadata_List, key_String] := Module[{position},
  position = FirstPosition[metadata, key, Missing["NotFound"]];
  If[MissingQ[position] || position[[1]] >= Length[metadata],
    Missing["NotFound"],
    metadata[[position[[1]] + 1]]
  ]
];

assert[
  "single-loader-and-unified-entry",
  $CHESSv2Version === "2.0.0-experiment" &&
  Context[SpectralPropagate] === "Global`" &&
  DownValues[SpectralPropagate] =!= {} &&
  MemberQ[First /@ Options[SpectralPropagate], "DeltaOrder"] &&
  MemberQ[First /@ Options[SpectralPropagate], "EpsilonDegree"]
];

(* The compatibility adapter must stop after the first valid convention.  This
   guards against silently evaluating an expensive or stateful matrix function
   two additional times at every collocation node. *)
ClearAll[lazyMatrix, lazyCallCount];
lazyCallCount = 0;
lazyMatrix[x_, "Precision" -> p_] := (
  lazyCallCount++;
  N[{{1}}, p]
);
lazyMatrixValue = CHESSUnifiedEvaluateMatrix[lazyMatrix, 1/3, 60];
assert[
  "matrix-adapter-stops-after-first-valid-convention",
  CHESSNumericMatrixQ[lazyMatrixValue] && lazyCallCount === 1
];
assert[
  "option-filter-preserves-rule-delayed",
  CHESSFilterOptionRules[
    {"Nodes" :> 12, "NotAnOption" :> 99},
    $CHESSCanonicalOptions
  ] === {"Nodes" :> 12}
];

(* Verify the thin canonical failure boundary independently of a particular
   LinearSolve failure trigger.  The copied core can emit solvefail and then
   have Return swallowed by Do; the wrapper must convert that message to a
   top-level $Failed result. *)
Quiet[checkedCanonicalFailure = Block[{SpectralPropagateSequential},
  SpectralPropagateSequential[___] := (
    Message[SpectralPropagate::solvefail, 1, $Failed, {}];
    "partial-success-shape"
  );
  CHESSCheckedCanonicalPropagate[Null, {{1}}, {0, 1}, {}]
]];
assert[
  "canonical-solvefail-message-is-promoted-to-failure",
  checkedCanonicalFailure === $Failed
];

(* ---------------------------------------------------------------------- *)
(* Canonical B1 routes                                                    *)
(* ---------------------------------------------------------------------- *)

ClearAll[canonicalB1, twoArgumentB1];
Options[canonicalB1] = {"Precision" -> 80};
canonicalB1[x_, OptionsPattern[]] :=
  N[{{1/2}}, OptionValue["Precision"]];
twoArgumentB1[x_, precision_] := N[{{1/2}}, precision];

canonicalBoundary = N[{{1, 0, 0, 0}}, 70];
canonicalExpected = N[{{1, 1/2, 1/8, 1/48}}, 70];

bareCanonical = SpectralPropagate[
  canonicalB1,
  canonicalBoundary,
  {0, 1},
  "Nodes" -> 22,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False
];
explicitCanonical = SpectralPropagate[
  {0, twoArgumentB1},
  canonicalBoundary,
  {0, 1},
  "Nodes" -> 22,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False
];

assert[
  "bare-b1-uses-canonical",
  bareCanonical =!= $Failed &&
  bareCanonical[[6]] === "SequentialEpsilon" &&
  maxError[bareCanonical[[3]], canonicalExpected] < 10^-55
];
assert[
  "explicit-only-b1-uses-canonical",
  explicitCanonical =!= $Failed &&
  explicitCanonical[[6]] === "SequentialEpsilon" &&
  maxError[explicitCanonical[[3]], bareCanonical[[3]]] < 10^-60
];

(* SparseArray is not a List expression.  The generic legacy route must still
   recognize it as a constant B1 matrix instead of trying to call it as a
   function. *)
sparseConstantCanonical = SpectralPropagate[
  SparseArray[{{1, 1} -> 1/2}, {1, 1}],
  canonicalBoundary,
  {0, 1},
  "Nodes" -> 22,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False
];
assert[
  "sparse-constant-b1-uses-canonical-adapter",
  sparseConstantCanonical =!= $Failed &&
  maxError[sparseConstantCanonical[[3]], canonicalExpected] < 10^-55
];

(* Every unified route must reject a square operator of the wrong physical
   dimension before the historical coefficient loops can form partial output. *)
ClearAll[wrongDimensionMatrix];
wrongDimensionMatrix[x_, p_] := N[IdentityMatrix[2], p];
Quiet[wrongBareConstant = SpectralPropagate[
  IdentityMatrix[2], N[{{1, 0}}, 50], {0, 1},
  "Nodes" -> 10, "Precision" -> 50, "ParallelEvaluation" -> False
]];
Quiet[wrongOnlyB1 = SpectralPropagate[
  {0, wrongDimensionMatrix}, N[{{1, 0}}, 50], {0, 1},
  "Nodes" -> 10, "Precision" -> 50, "ParallelEvaluation" -> False
]];
Quiet[wrongOnlyB0 = SpectralPropagate[
  {wrongDimensionMatrix}, N[{{1}}, 50], {0, 1},
  "Nodes" -> 10, "Precision" -> 50, "ParallelEvaluation" -> False,
  "DeltaOrder" -> 4
]];
Quiet[wrongMixed = SpectralPropagate[
  {wrongDimensionMatrix, {{1/3}}},
  N[{{1, 0}}, 50], {0, 1},
  "Nodes" -> 10, "Precision" -> 50, "ParallelEvaluation" -> False
]];
assert[
  "operator-dimensions-checked-on-every-route",
  wrongBareConstant === $Failed && wrongOnlyB1 === $Failed &&
  wrongOnlyB0 === $Failed && wrongMixed === $Failed
];

(* Regular canonical segmentation is implemented only in the thin wrapper;
   each segment still uses the original canonical solver unchanged. *)
segmentedCanonical = SpectralPropagate[
  {0, twoArgumentB1},
  canonicalBoundary,
  {0, 1},
  "Segments" -> 2,
  "Nodes" -> 18,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False
];
assert[
  "canonical-regular-segmentation",
  segmentedCanonical =!= $Failed &&
  segmentedCanonical[[6]] === {"CanonicalSegments", 2} &&
  maxError[segmentedCanonical[[3]], canonicalExpected] < 10^-50
];

(* ---------------------------------------------------------------------- *)
(* B0-only fake-delta routes                                              *)
(* ---------------------------------------------------------------------- *)

ClearAll[scalarB0, diagonalB0];
scalarB0[x_, precision_] := N[{{1/3}}, precision];
diagonalB0[x_, precision_] := N[
  DiagonalMatrix[{1/3, -1/5}], precision
];

autoDelta = SpectralPropagate[
  {scalarB0},
  N[{{1}}, 70],
  {0, 1},
  "Nodes" -> 22,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False,
  "DeltaTolerance" -> 10^-45,
  "MaxDeltaOrder" -> 64
];
fixedDelta = SpectralPropagate[
  {scalarB0},
  N[{{1}}, 70],
  {0, 1},
  "Nodes" -> 22,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False,
  "DeltaOrder" -> 48,
  "DeltaTolerance" -> 10^-45
];

assert[
  "b0-only-auto-delta-analytic",
  autoDelta =!= $Failed &&
  autoDelta[[6, 1]] === "FakeDelta" &&
  autoDelta[[6, 2]] === "Automatic" &&
  IntegerQ[metadataValue[autoDelta[[6]], "SelectedOrder"]] &&
  metadataValue[autoDelta[[6]], "TailEstimate"] < 10^-45 &&
  Abs[autoDelta[[3, 1, 1]] - Exp[N[1/3, 70]]] < 10^-45
];
assert[
  "b0-only-auto-versus-higher-fixed-order",
  fixedDelta =!= $Failed &&
  maxError[autoDelta[[3]], fixedDelta[[3]]] < 10^-45
];

(* A time-dependent matrix whose values at different points do not commute is
   a more representative spectral-transport test than a scalar exponential.
   There is no elementary closed form to lean on here, so the automatic result
   is compared with a deliberately higher fixed delta order while keeping the
   same collocation grid.  The explicit commutator check prevents this fixture
   from accidentally degenerating into another commuting example. *)
ClearAll[noncommutingB0];
noncommutingB0[x_, precision_] := N[
  {{0, 1 + x}, {-x/3, 1/5}}, precision
];
noncommutingBoundary = N[{{1}, {-2/3}}, 70];
noncommutingAuto = SpectralPropagate[
  {noncommutingB0},
  noncommutingBoundary,
  {0, 1},
  "Nodes" -> 26,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False,
  "DeltaTolerance" -> 10^-35,
  "MaxDeltaOrder" -> 64
];
noncommutingFixed = SpectralPropagate[
  {noncommutingB0},
  noncommutingBoundary,
  {0, 1},
  "Nodes" -> 26,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False,
  "DeltaOrder" -> 48,
  "DeltaTolerance" -> 10^-35
];
assert[
  "b0-only-noncommuting-auto-versus-fixed",
  noncommutingAuto =!= $Failed && noncommutingFixed =!= $Failed &&
  noncommutingB0[0, 70] . noncommutingB0[1, 70] =!=
    noncommutingB0[1, 70] . noncommutingB0[0, 70] &&
  maxError[noncommutingAuto[[3]], noncommutingFixed[[3]]] < 10^-35
];

(* Several physical boundary columns must remain independent under B0.  This
   test also verifies that the public result has the original n x q shape and
   not the internal n x (K+1) dummy-coefficient shape. *)
multiBoundary = N[{{1, 2}, {3, 4}}, 70];
multiExpected = N[
  DiagonalMatrix[{Exp[1/3], Exp[-1/5]}] . multiBoundary,
  70
];
multiDelta = SpectralPropagate[
  {diagonalB0},
  multiBoundary,
  {0, 1},
  "Segments" -> 2,
  "Nodes" -> 20,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False,
  "DeltaTolerance" -> 10^-40,
  "MaxDeltaOrder" -> 64
];
assert[
  "b0-only-multiple-columns-and-segments",
  multiDelta =!= $Failed &&
  Dimensions[multiDelta[[3]]] === Dimensions[multiBoundary] &&
  multiDelta[[6, 1]] === "FakeDeltaSegments" &&
  maxError[multiDelta[[3]], multiExpected] < 10^-38
];

(* B0=0 is a useful limiting case: every higher dummy coefficient vanishes and
   the physical state must be exactly unchanged. *)
zeroDelta = SpectralPropagate[
  {0},
  multiBoundary,
  {0, 1},
  "Nodes" -> 12,
  "Precision" -> 60,
  "WorkingPrecisionA" -> 80,
  "ParallelEvaluation" -> False,
  "DeltaTolerance" -> 10^-35
];
assert[
  "b0-zero-limiting-case",
  zeroDelta =!= $Failed && zeroDelta[[3]] === SetPrecision[multiBoundary, 60]
];

(* The automatic selector must fail closed when the maximum order is plainly
   insufficient; returning a partial sum here would defeat automatic control. *)
ClearAll[largeB0];
largeB0[x_, precision_] := N[{{20}}, precision];
Quiet[insufficientDelta = SpectralPropagate[
  {largeB0},
  N[{{1}}, 60],
  {0, 1},
  "Nodes" -> 16,
  "Precision" -> 60,
  "WorkingPrecisionA" -> 80,
  "ParallelEvaluation" -> False,
  "DeltaTolerance" -> 10^-40,
  "MaxDeltaOrder" -> 4
]];
assert["b0-auto-order-fails-closed", insufficientDelta === $Failed];

(* ---------------------------------------------------------------------- *)
(* General polynomial-epsilon route                                      *)
(* ---------------------------------------------------------------------- *)

ClearAll[mixedB0, mixedB1];
mixedB0[x_, precision_] := N[{{1/3}}, precision];
mixedB1[x_, precision_] := N[{{2/5}}, precision];
mixedBoundary = N[{{1, 0, 0, 0}}, 70];
mixedExpected = N[{{
  Exp[1/3],
  Exp[1/3] (2/5),
  Exp[1/3] (2/5)^2/2,
  Exp[1/3] (2/5)^3/6
}}, 70];

mixedUnified = SpectralPropagate[
  {mixedB0, mixedB1},
  mixedBoundary,
  {0, 1},
  "Nodes" -> 22,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False
];
mixedDirect = SpectralPropagatePolynomialEpsilon[
  {mixedB0, mixedB1},
  mixedBoundary,
  {0, 1},
  "Nodes" -> 22,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False
];
assert[
  "mixed-route-preserves-noncanonical-algorithm",
  mixedUnified =!= $Failed && mixedDirect =!= $Failed &&
  ListQ[mixedUnified[[6]]] &&
  mixedUnified[[6, 1]] === "PolynomialEpsilonActiveSupport" &&
  maxError[mixedUnified[[3]], mixedDirect[[3]]] < 10^-60 &&
  (* At 22 Lobatto nodes the copied non-canonical core has about 42 correct
     decimal digits on this fixture.  The tighter equality above is the actual
     regression check for the new dispatch layer. *)
  maxError[mixedUnified[[3]], mixedExpected] < 10^-40
];

(* Literal constant matrices exercise the normalization adapter. *)
constantMixed = SpectralPropagate[
  {{{1/3}}, {{2/5}}},
  mixedBoundary,
  {0, 1},
  "Nodes" -> 22,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False
];
assert[
  "constant-polynomial-matrix-adapter",
  constantMixed =!= $Failed &&
  maxError[constantMixed[[3]], mixedExpected] < 10^-40
];

(* Mixed routes must accept the canonical keyword-precision convention too;
   this is handled by MatrixAdapters rather than by changing the old solver. *)
ClearAll[keywordB0, keywordB1];
keywordB0[x_, "Precision" -> p_] := N[{{1/3}}, p];
keywordB1[x_, "Precision" -> p_] := N[{{2/5}}, p];
keywordMixed = SpectralPropagate[
  {keywordB0, keywordB1},
  mixedBoundary,
  {0, 1},
  "Nodes" -> 22,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False
];
assert[
  "mixed-route-accepts-keyword-precision-evaluators",
  keywordMixed =!= $Failed &&
  maxError[keywordMixed[[3]], mixedExpected] < 10^-40
];

(* A coefficient that fails only in the second segment checks the historical
   non-canonical loop boundary.  The whole propagation must fail, rather than
   returning the first segment with success-shaped metadata. *)
ClearAll[secondSegmentB0, secondSegmentB1];
secondSegmentB0[x_, p_] := If[
  TrueQ[x > 11/20], $Failed, N[{{1/3}}, p]
];
secondSegmentB1[x_, p_] := N[{{2/5}}, p];
Quiet[secondSegmentFailure = SpectralPropagate[
  {secondSegmentB0, secondSegmentB1},
  mixedBoundary,
  {0, 1},
  "Segments" -> 2,
  "Nodes" -> 16,
  "Precision" -> 60,
  "WorkingPrecisionA" -> 80,
  "ParallelEvaluation" -> False
]];
assert[
  "noncanonical-later-segment-failure-propagates",
  secondSegmentFailure === $Failed
];

(* RuleDelayed is valid option syntax.  The outer mixed dispatcher preserves
   it, while the recursive segment call must remove the old value before adding
   Segments->1; otherwise the solver recurses indefinitely. *)
mixedDelayedSegments = SpectralPropagate[
  {mixedB0, mixedB1},
  mixedBoundary,
  {0, 1},
  "Segments" :> 2,
  "Nodes" -> 18,
  "Precision" -> 60,
  "WorkingPrecisionA" -> 80,
  "ParallelEvaluation" -> False
];
assert[
  "noncanonical-rule-delayed-segments-terminate",
  mixedDelayedSegments =!= $Failed &&
  mixedDelayedSegments[[6]] === "PolynomialEpsilonSegments"
];

(* ---------------------------------------------------------------------- *)
(* Fail-closed API boundaries and provenance                             *)
(* ---------------------------------------------------------------------- *)

Quiet[badSegments = SpectralPropagate[
  {scalarB0}, {{1}}, {0, 1}, "Segments" -> 0
]];
Quiet[badDeltaOrder = SpectralPropagate[
  {scalarB0}, {{1}}, {0, 1}, "DeltaOrder" -> -1
]];
Quiet[b0Endpoint = SpectralPropagate[
  {scalarB0}, {{1}}, {0, 1}, "RegularizedEndpoints" -> {"Right"}
]];
Quiet[b0EndpointData = SpectralPropagate[
  {scalarB0}, {{1}}, {0, 1}, "EndpointData" -> {1 -> {}}
]];
Quiet[canonicalEndpointData = SpectralPropagate[
  {0, twoArgumentB1}, canonicalBoundary, {0, 1},
  "EndpointData" -> {1 -> {}}
]];
Quiet[mixedOrphanEndpointData = SpectralPropagate[
  {mixedB0, mixedB1}, mixedBoundary, {0, 1},
  "EndpointData" -> {1 -> {}}
]];
Quiet[mixedSegmentedEndpoint = SpectralPropagate[
  {mixedB0, mixedB1}, mixedBoundary, {0, 1},
  "Segments" -> 2,
  "RegularizedEndpoints" -> {"Right"},
  "EndpointData" -> {1 -> {}}
]];
Quiet[mixedOutsideEndpoint = SpectralPropagate[
  {mixedB0, mixedB1}, mixedBoundary, {0, 1},
  "RegularizedEndpoints" -> {2},
  "EndpointData" -> {2 -> {}}
]];
assert[
  "invalid-and-unsupported-routes-fail-closed",
  badSegments === $Failed && badDeltaOrder === $Failed &&
  b0Endpoint === $Failed && b0EndpointData === $Failed &&
  canonicalEndpointData === $Failed &&
  mixedOrphanEndpointData === $Failed &&
  mixedSegmentedEndpoint === $Failed &&
  mixedOutsideEndpoint === $Failed
];

hashHex[file_] := IntegerString[FileHash[file, "SHA256"], 16, 64];
assert[
  "original-sources-remain-untouched",
  hashHex["/home/liuyuanche/Program/Chebyshev/Chess_release/Chess.wl"] ===
    "036a670f5983992dbe70c6b335ff6b498a0e1a868b8de9824b9918092c97d583" &&
  hashHex["/home/liuyuanche/Program/Chebyshev_nonUT/ChessNonCanonical.wl"] ===
    "48294b347777fca190b062b53e2fe2235de9fe05660d9ab0957a008150e54035"
];
assert[
  "canonical-core-is-verbatim-copy",
  hashHex[FileNameJoin[{testRoot, "Core", "Canonical.wl"}]] ===
    "036a670f5983992dbe70c6b335ff6b498a0e1a868b8de9824b9918092c97d583"
];

(* A malformed assert call otherwise remains unevaluated and can make a key
   check disappear without affecting the process exit code. *)
If[assertionCount =!= expectedAssertionCount,
  AppendTo[failures, "assertion-count"];
  Print[
    "FAIL: assertion-count (expected ", expectedAssertionCount,
    ", observed ", assertionCount, ")"
  ]
];

If[failures === {},
  Print["CHESSV2_ALL_TESTS_PASSED"];
  Exit[0],
  Print[InputForm[{"CHESSV2_TEST_FAILURES", failures}]];
  Exit[1]
];
