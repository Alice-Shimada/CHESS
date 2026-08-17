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
expectedAssertionCount = 25;
CheckAssertion[name_String, condition_] := If[
  assertionCount++;
  TrueQ[condition],
  Print["PASS: ", name],
  AppendTo[failures, name];
  Print["FAIL: ", name]
];

MaxError[a_, b_] := Max[Abs[Flatten[N[a - b, 50]]]];
MetadataValue[metadata_List, key_String] := Module[{position},
  position = FirstPosition[metadata, key, Missing["NotFound"]];
  If[MissingQ[position] || position[[1]] >= Length[metadata],
    Missing["NotFound"],
    metadata[[position[[1]] + 1]]
  ]
];

CheckAssertion[
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
ClearAll[LazyMatrix, lazyCallCount];
lazyCallCount = 0;
LazyMatrix[x_, "Precision" -> p_] := (
  lazyCallCount++;
  N[{{1}}, p]
);
lazyMatrixValue = CHESSUnifiedEvaluateMatrix[LazyMatrix, 1/3, 60];
CheckAssertion[
  "matrix-adapter-stops-after-first-valid-convention",
  CHESSNumericMatrixQ[lazyMatrixValue] && lazyCallCount === 1
];
CheckAssertion[
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
CheckAssertion[
  "canonical-solvefail-message-is-promoted-to-failure",
  checkedCanonicalFailure === $Failed
];

(* ---------------------------------------------------------------------- *)
(* Canonical B1 routes                                                    *)
(* ---------------------------------------------------------------------- *)

ClearAll[CanonicalB1, TwoArgumentB1];
Options[CanonicalB1] = {"Precision" -> 80};
CanonicalB1[x_, OptionsPattern[]] :=
  N[{{1/2}}, OptionValue["Precision"]];
TwoArgumentB1[x_, precision_] := N[{{1/2}}, precision];

canonicalBoundary = N[{{1, 0, 0, 0}}, 70];
canonicalExpected = N[{{1, 1/2, 1/8, 1/48}}, 70];

bareCanonical = SpectralPropagate[
  CanonicalB1,
  canonicalBoundary,
  {0, 1},
  "Nodes" -> 22,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False
];
explicitCanonical = SpectralPropagate[
  {0, TwoArgumentB1},
  canonicalBoundary,
  {0, 1},
  "Nodes" -> 22,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False
];

CheckAssertion[
  "bare-b1-uses-canonical",
  bareCanonical =!= $Failed &&
  bareCanonical[[6]] === "SequentialEpsilon" &&
  MaxError[bareCanonical[[3]], canonicalExpected] < 10^-55
];
CheckAssertion[
  "explicit-only-b1-uses-canonical",
  explicitCanonical =!= $Failed &&
  explicitCanonical[[6]] === "SequentialEpsilon" &&
  MaxError[explicitCanonical[[3]], bareCanonical[[3]]] < 10^-60
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
CheckAssertion[
  "sparse-constant-b1-uses-canonical-adapter",
  sparseConstantCanonical =!= $Failed &&
  MaxError[sparseConstantCanonical[[3]], canonicalExpected] < 10^-55
];

(* Every unified route must reject a square operator of the wrong physical
   dimension before the historical coefficient loops can form partial output. *)
ClearAll[WrongDimensionMatrix];
WrongDimensionMatrix[x_, p_] := N[IdentityMatrix[2], p];
Quiet[wrongBareConstant = SpectralPropagate[
  IdentityMatrix[2], N[{{1, 0}}, 50], {0, 1},
  "Nodes" -> 10, "Precision" -> 50, "ParallelEvaluation" -> False
]];
Quiet[wrongOnlyB1 = SpectralPropagate[
  {0, WrongDimensionMatrix}, N[{{1, 0}}, 50], {0, 1},
  "Nodes" -> 10, "Precision" -> 50, "ParallelEvaluation" -> False
]];
Quiet[wrongOnlyB0 = SpectralPropagate[
  {WrongDimensionMatrix}, N[{{1}}, 50], {0, 1},
  "Nodes" -> 10, "Precision" -> 50, "ParallelEvaluation" -> False,
  "DeltaOrder" -> 4
]];
Quiet[wrongMixed = SpectralPropagate[
  {WrongDimensionMatrix, {{1/3}}},
  N[{{1, 0}}, 50], {0, 1},
  "Nodes" -> 10, "Precision" -> 50, "ParallelEvaluation" -> False
]];
CheckAssertion[
  "operator-dimensions-checked-on-every-route",
  wrongBareConstant === $Failed && wrongOnlyB1 === $Failed &&
  wrongOnlyB0 === $Failed && wrongMixed === $Failed
];

(* Regular canonical segmentation is implemented only in the thin wrapper;
   each segment still uses the original canonical solver unchanged. *)
segmentedCanonical = SpectralPropagate[
  {0, TwoArgumentB1},
  canonicalBoundary,
  {0, 1},
  "Segments" -> 2,
  "Nodes" -> 18,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False
];
CheckAssertion[
  "canonical-regular-segmentation",
  segmentedCanonical =!= $Failed &&
  segmentedCanonical[[6]] === {"CanonicalSegments", 2} &&
  MaxError[segmentedCanonical[[3]], canonicalExpected] < 10^-50
];

(* ---------------------------------------------------------------------- *)
(* B0-only fake-delta routes                                              *)
(* ---------------------------------------------------------------------- *)

ClearAll[ScalarB0, DiagonalB0];
ScalarB0[x_, precision_] := N[{{1/3}}, precision];
DiagonalB0[x_, precision_] := N[
  DiagonalMatrix[{1/3, -1/5}], precision
];

autoDelta = SpectralPropagate[
  {ScalarB0},
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
  {ScalarB0},
  N[{{1}}, 70],
  {0, 1},
  "Nodes" -> 22,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False,
  "DeltaOrder" -> 48,
  "DeltaTolerance" -> 10^-45
];

CheckAssertion[
  "b0-only-auto-delta-analytic",
  autoDelta =!= $Failed &&
  autoDelta[[6, 1]] === "FakeDelta" &&
  autoDelta[[6, 2]] === "Automatic" &&
  IntegerQ[MetadataValue[autoDelta[[6]], "SelectedOrder"]] &&
  MetadataValue[autoDelta[[6]], "TailEstimate"] < 10^-45 &&
  Abs[autoDelta[[3, 1, 1]] - Exp[N[1/3, 70]]] < 10^-45
];
CheckAssertion[
  "b0-only-auto-versus-higher-fixed-order",
  fixedDelta =!= $Failed &&
  MaxError[autoDelta[[3]], fixedDelta[[3]]] < 10^-45
];

(* A time-dependent matrix whose values at different points do not commute is
   a more representative spectral-transport test than a scalar exponential.
   There is no elementary closed form to lean on here, so the automatic result
   is compared with a deliberately higher fixed delta order while keeping the
   same collocation grid.  The explicit commutator check prevents this fixture
   from accidentally degenerating into another commuting example. *)
ClearAll[NoncommutingB0];
NoncommutingB0[x_, precision_] := N[
  {{0, 1 + x}, {-x/3, 1/5}}, precision
];
noncommutingBoundary = N[{{1}, {-2/3}}, 70];
noncommutingAuto = SpectralPropagate[
  {NoncommutingB0},
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
  {NoncommutingB0},
  noncommutingBoundary,
  {0, 1},
  "Nodes" -> 26,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False,
  "DeltaOrder" -> 48,
  "DeltaTolerance" -> 10^-35
];
CheckAssertion[
  "b0-only-noncommuting-auto-versus-fixed",
  noncommutingAuto =!= $Failed && noncommutingFixed =!= $Failed &&
  NoncommutingB0[0, 70] . NoncommutingB0[1, 70] =!=
    NoncommutingB0[1, 70] . NoncommutingB0[0, 70] &&
  MaxError[noncommutingAuto[[3]], noncommutingFixed[[3]]] < 10^-35
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
  {DiagonalB0},
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
CheckAssertion[
  "b0-only-multiple-columns-and-segments",
  multiDelta =!= $Failed &&
  Dimensions[multiDelta[[3]]] === Dimensions[multiBoundary] &&
  multiDelta[[6, 1]] === "FakeDeltaSegments" &&
  MaxError[multiDelta[[3]], multiExpected] < 10^-38
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
CheckAssertion[
  "b0-zero-limiting-case",
  zeroDelta =!= $Failed && zeroDelta[[3]] === SetPrecision[multiBoundary, 60]
];

(* The automatic selector must fail closed when the maximum order is plainly
   insufficient; returning a partial sum here would defeat automatic control. *)
ClearAll[LargeB0];
LargeB0[x_, precision_] := N[{{20}}, precision];
Quiet[insufficientDelta = SpectralPropagate[
  {LargeB0},
  N[{{1}}, 60],
  {0, 1},
  "Nodes" -> 16,
  "Precision" -> 60,
  "WorkingPrecisionA" -> 80,
  "ParallelEvaluation" -> False,
  "DeltaTolerance" -> 10^-40,
  "MaxDeltaOrder" -> 4
]];
CheckAssertion["b0-auto-order-fails-closed", insufficientDelta === $Failed];

(* ---------------------------------------------------------------------- *)
(* General polynomial-epsilon route                                      *)
(* ---------------------------------------------------------------------- *)

ClearAll[MixedB0, MixedB1];
MixedB0[x_, precision_] := N[{{1/3}}, precision];
MixedB1[x_, precision_] := N[{{2/5}}, precision];
mixedBoundary = N[{{1, 0, 0, 0}}, 70];
mixedExpected = N[{{
  Exp[1/3],
  Exp[1/3] (2/5),
  Exp[1/3] (2/5)^2/2,
  Exp[1/3] (2/5)^3/6
}}, 70];

mixedUnified = SpectralPropagate[
  {MixedB0, MixedB1},
  mixedBoundary,
  {0, 1},
  "Nodes" -> 22,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False
];
mixedDirect = SpectralPropagatePolynomialEpsilon[
  {MixedB0, MixedB1},
  mixedBoundary,
  {0, 1},
  "Nodes" -> 22,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False
];
CheckAssertion[
  "mixed-route-preserves-noncanonical-algorithm",
  mixedUnified =!= $Failed && mixedDirect =!= $Failed &&
  ListQ[mixedUnified[[6]]] &&
  mixedUnified[[6, 1]] === "PolynomialEpsilonActiveSupport" &&
  MaxError[mixedUnified[[3]], mixedDirect[[3]]] < 10^-60 &&
  (* At 22 Lobatto nodes the copied non-canonical core has about 42 correct
     decimal digits on this fixture.  The tighter equality above is the actual
     regression check for the new dispatch layer. *)
  MaxError[mixedUnified[[3]], mixedExpected] < 10^-40
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
CheckAssertion[
  "constant-polynomial-matrix-adapter",
  constantMixed =!= $Failed &&
  MaxError[constantMixed[[3]], mixedExpected] < 10^-40
];

(* Mixed routes must accept the canonical keyword-precision convention too;
   this is handled by MatrixAdapters rather than by changing the old solver. *)
ClearAll[KeywordB0, KeywordB1];
KeywordB0[x_, "Precision" -> p_] := N[{{1/3}}, p];
KeywordB1[x_, "Precision" -> p_] := N[{{2/5}}, p];
keywordMixed = SpectralPropagate[
  {KeywordB0, KeywordB1},
  mixedBoundary,
  {0, 1},
  "Nodes" -> 22,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False
];
CheckAssertion[
  "mixed-route-accepts-keyword-precision-evaluators",
  keywordMixed =!= $Failed &&
  MaxError[keywordMixed[[3]], mixedExpected] < 10^-40
];

(* A FORM/FLINT-style batch evaluator returns every epsilon coefficient matrix
   in one call.  Explicit EpsilonDegree must route it through the unified
   non-canonical entry without splitting it into repeated evaluator calls. *)
ClearAll[PolynomialScalar, PolynomialBatch];
polynomialBatchCallCount = 0;
PolynomialScalar[x_, requestedPrecision_] := N[
  {{{1/5}}, {{2/7}}}, requestedPrecision
];
PolynomialBatch[points_, requestedPrecision_, threads_] := (
  polynomialBatchCallCount++;
  PolynomialScalar[#, requestedPrecision] & /@ points
);
polynomialNodeEvaluator = CHESSNodeEvaluator[
  PolynomialScalar, PolynomialBatch
];
unifiedPolynomialBatch = SpectralPropagate[
  polynomialNodeEvaluator, {{1}}, {0, 1},
  "Nodes" -> 12,
  "Precision" -> 60,
  "WorkingPrecisionA" -> 70,
  "EpsilonDegree" -> 1,
  "ParallelEvaluation" -> True,
  "ParallelKernels" -> 2
];
unifiedPolynomialBatchCallCount = polynomialBatchCallCount;
directPolynomialBatch = SpectralPropagatePolynomialEpsilon[
  polynomialNodeEvaluator, {{1}}, {0, 1},
  "Nodes" -> 12,
  "Precision" -> 60,
  "WorkingPrecisionA" -> 70,
  "EpsilonDegree" -> 1,
  "ParallelEvaluation" -> True,
  "ParallelKernels" -> 2
];
ClearAll[WrongDimensionPolynomialBatch];
WrongDimensionPolynomialBatch[points_, requestedPrecision_, threads_] :=
  ConstantArray[
    N[{{{1/5}}, (2/7) IdentityMatrix[2]}, requestedPrecision],
    Length[points]
  ];
wrongDimensionPolynomialResult = Quiet @ SpectralPropagate[
  CHESSNodeEvaluator[PolynomialScalar, WrongDimensionPolynomialBatch],
  IdentityMatrix[2],
  {0, 1},
  "Nodes" -> 8,
  "Precision" -> 50,
  "WorkingPrecisionA" -> 60,
  "EpsilonDegree" -> 1
];
CheckAssertion[
  "unified-polynomial-batch-evaluator-route",
  unifiedPolynomialBatch =!= $Failed && directPolynomialBatch =!= $Failed &&
  unifiedPolynomialBatchCallCount === 1 &&
  polynomialBatchCallCount === 2 &&
  wrongDimensionPolynomialResult === $Failed &&
  MaxError[unifiedPolynomialBatch[[3]], directPolynomialBatch[[3]]] < 10^-50
];

(* A coefficient that fails only in the second segment checks the historical
   non-canonical loop boundary.  The whole propagation must fail, rather than
   returning the first segment with success-shaped metadata. *)
ClearAll[SecondSegmentB0, SecondSegmentB1];
SecondSegmentB0[x_, p_] := If[
  TrueQ[x > 11/20], $Failed, N[{{1/3}}, p]
];
SecondSegmentB1[x_, p_] := N[{{2/5}}, p];
Quiet[secondSegmentFailure = SpectralPropagate[
  {SecondSegmentB0, SecondSegmentB1},
  mixedBoundary,
  {0, 1},
  "Segments" -> 2,
  "Nodes" -> 16,
  "Precision" -> 60,
  "WorkingPrecisionA" -> 80,
  "ParallelEvaluation" -> False
]];
CheckAssertion[
  "noncanonical-later-segment-failure-propagates",
  secondSegmentFailure === $Failed
];

(* RuleDelayed is valid option syntax.  The outer mixed dispatcher preserves
   it, while the recursive segment call must remove the old value before adding
   Segments->1; otherwise the solver recurses indefinitely. *)
mixedDelayedSegments = SpectralPropagate[
  {MixedB0, MixedB1},
  mixedBoundary,
  {0, 1},
  "Segments" :> 2,
  "Nodes" -> 18,
  "Precision" -> 60,
  "WorkingPrecisionA" -> 80,
  "ParallelEvaluation" -> False
];
CheckAssertion[
  "noncanonical-rule-delayed-segments-terminate",
  mixedDelayedSegments =!= $Failed &&
  mixedDelayedSegments[[6]] === "PolynomialEpsilonSegments"
];

(* ---------------------------------------------------------------------- *)
(* Fail-closed API boundaries and provenance                             *)
(* ---------------------------------------------------------------------- *)

Quiet[badSegments = SpectralPropagate[
  {ScalarB0}, {{1}}, {0, 1}, "Segments" -> 0
]];
Quiet[badDeltaOrder = SpectralPropagate[
  {ScalarB0}, {{1}}, {0, 1}, "DeltaOrder" -> -1
]];
Quiet[b0Endpoint = SpectralPropagate[
  {ScalarB0}, {{1}}, {0, 1}, "RegularizedEndpoints" -> {"Right"}
]];
Quiet[b0EndpointData = SpectralPropagate[
  {ScalarB0}, {{1}}, {0, 1}, "EndpointData" -> {1 -> {}}
]];
Quiet[canonicalEndpointData = SpectralPropagate[
  {0, TwoArgumentB1}, canonicalBoundary, {0, 1},
  "EndpointData" -> {1 -> {}}
]];
Quiet[mixedOrphanEndpointData = SpectralPropagate[
  {MixedB0, MixedB1}, mixedBoundary, {0, 1},
  "EndpointData" -> {1 -> {}}
]];
Quiet[mixedSegmentedEndpoint = SpectralPropagate[
  {MixedB0, MixedB1}, mixedBoundary, {0, 1},
  "Segments" -> 2,
  "RegularizedEndpoints" -> {"Right"},
  "EndpointData" -> {1 -> {}}
]];
Quiet[mixedOutsideEndpoint = SpectralPropagate[
  {MixedB0, MixedB1}, mixedBoundary, {0, 1},
  "RegularizedEndpoints" -> {2},
  "EndpointData" -> {2 -> {}}
]];
CheckAssertion[
  "invalid-and-unsupported-routes-fail-closed",
  badSegments === $Failed && badDeltaOrder === $Failed &&
  b0Endpoint === $Failed && b0EndpointData === $Failed &&
  canonicalEndpointData === $Failed &&
  mixedOrphanEndpointData === $Failed &&
  mixedSegmentedEndpoint === $Failed &&
  mixedOutsideEndpoint === $Failed
];

(* User-facing conversion helpers must produce the same polynomial batch
   contract from function lists and from explicit path-variable expressions. *)
ClearAll[AdapterB0, AdapterB1, AdapterFailure, adapterT];
AdapterB0[x_, p_] := N[{{1/3 + x/10}}, p];
AdapterB1[x_, p_] := N[{{2/5 - x/20}}, p];
AdapterFailure[x_, p_] := If[x > 2/5, IdentityMatrix[2], N[{{1}}, p]];
functionAdapter = CHESSNativeMatrixAdapter[{AdapterB0, AdapterB1}, 1];
expressionAdapter = CHESSNativeMatrixAdapter[
  {{{1/3 + adapterT/10}}, {{2/5 - adapterT/20}}}, adapterT, 1
];
functionAdapterValues = CHESSNodeEvaluatorBatch[functionAdapter][
  {1/4, 1/2}, 60, 1
];
expressionAdapterValues = CHESSNodeEvaluatorBatch[expressionAdapter][
  {1/4, 1/2}, 60, 1
];
badExpressionAdapter = Quiet @ CHESSNativeMatrixAdapter[
  {{{1}}, IdentityMatrix[2]}, adapterT, 1
];
failingAdapter = CHESSNativeMatrixAdapter[{AdapterFailure}, 1];
failingAdapterValues = Quiet @ CHESSNodeEvaluatorBatch[failingAdapter][
  {1/4, 1/2}, 60, 1
];
badFlintAdapter = Quiet @ CHESSNativeFlintMatrixAdapter[
  "missing", "missing", Identity, Identity, 0
];
CheckAssertion[
  "native-matrix-adapter-function-and-expression-formats",
  CHESSNodeEvaluatorQ[functionAdapter] &&
  CHESSNodeEvaluatorQ[expressionAdapter] &&
  functionAdapterValues =!= $Failed &&
  expressionAdapterValues =!= $Failed &&
  MaxError[functionAdapterValues, expressionAdapterValues] < 10^-55 &&
  badExpressionAdapter === $Failed && failingAdapterValues === $Failed &&
  badFlintAdapter === $Failed
];

HashHex[file_] := IntegerString[FileHash[file, "SHA256"], 16, 64];
CheckAssertion[
  "original-sources-remain-untouched",
  HashHex["/home/liuyuanche/Program/Chebyshev/Chess_release/Chess.wl"] ===
    "036a670f5983992dbe70c6b335ff6b498a0e1a868b8de9824b9918092c97d583" &&
  HashHex["/home/liuyuanche/Program/Chebyshev_nonUT/ChessNonCanonical.wl"] ===
    "48294b347777fca190b062b53e2fe2235de9fe05660d9ab0957a008150e54035"
];
CheckAssertion[
  "canonical-core-is-verbatim-copy",
  HashHex[FileNameJoin[{testRoot, "Core", "Canonical.wl"}]] ===
    "036a670f5983992dbe70c6b335ff6b498a0e1a868b8de9824b9918092c97d583"
];

(* A malformed CheckAssertion call otherwise remains unevaluated and can make a key
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
