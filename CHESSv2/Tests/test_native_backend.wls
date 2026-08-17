#!/usr/bin/env wolframscript

(* Focused test for the optional C++ regular-point propagation module.  This
   file is separate from run_all.wls because it requires a compiled executable:

       make -C CHESSv2/Native/Propagation
*)

ClearAll["Global`*"];
$HistoryLength = 0;
$MaxExtraPrecision = 3000;

testRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{testRoot, "CHESSv2.wl"}]];

failures = {};
assert[name_String, condition_] := If[
  TrueQ[condition],
  Print["PASS: ", name],
  AppendTo[failures, name];
  Print["FAIL: ", name]
];
maxError[a_, b_] := Max[Abs[Flatten[N[a - b, 50]]]];

If[!CHESSNativeAvailableQ[],
  Print["FAIL: native-executable-missing"];
  Exit[2]
];

ClearAll[scalarOperator];
scalarOperator[x_, precision_] := N[{{3/7}}, precision];
precision = 70;
degree = 20;
rules = {
  "Nodes" -> 32,
  "Precision" -> precision,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False,
  "ParallelKernels" -> 1
};

handle = CHESSNativePrepare[scalarOperator, 1, {0, 1}, rules];
assert["native-prepare", AssociationQ[handle]];

canonicalBoundary = N[{{1} ~Join~ ConstantArray[0, degree]}, precision];
canonicalRun = CHESSNativeRun[handle, {Transpose[canonicalBoundary]}];
canonicalResult = If[
  AssociationQ[canonicalRun],
  CHESSNativeCanonicalResult[handle, canonicalRun, "Full"],
  $Failed
];
canonicalExpected = N[
  {Table[(3/7)^order/order!, {order, 0, degree}]},
  precision
];
assert[
  "native-canonical-full-result",
  canonicalResult =!= $Failed &&
  Dimensions[canonicalResult[[2]]] === {33, degree + 1} &&
  maxError[canonicalResult[[3]], canonicalExpected] < 10^-60
];

publicNative = SpectralPropagate[
  {0, scalarOperator},
  canonicalBoundary,
  {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Native",
  "NativeHandle" -> handle,
  "ResultData" -> "Full"
];
publicEndpointOnly = SpectralPropagate[
  {0, scalarOperator},
  canonicalBoundary,
  {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Native",
  "NativeHandle" -> handle,
  "ResultData" -> "Endpoint"
];
assert[
  "unified-spectral-propagate-native-canonical",
  publicNative =!= $Failed && publicEndpointOnly =!= $Failed &&
  publicNative[[6, 2]] === "Native" &&
  MissingQ[publicEndpointOnly[[2]]] &&
  maxError[publicNative[[3]], canonicalExpected] < 10^-60 &&
  maxError[publicEndpointOnly[[3]], canonicalExpected] < 10^-60
];

(* Two independent physical columns exercise the batched B0/fake-delta data
   layout.  Only C0 is nonzero at the left boundary. *)
physicalBoundary = N[{{1, -2}}, precision];
fakeTensor = ConstantArray[0, {2, degree + 1, 1}];
fakeTensor[[All, 1, 1]] = First[physicalBoundary];
fakeRun = CHESSNativeRun[handle, fakeTensor];
fakeFetch = If[
  AssociationQ[fakeRun], CHESSNativeFetch[handle, 1, degree], $Failed
];
fakeEndpoint = If[
  AssociationQ[fakeFetch],
  Lookup[fakeFetch, "Data"][[All, -1, 1]],
  $Failed
];
fakeExpected = N[
  First[physicalBoundary] Sum[(3/7)^order/order!, {order, 0, degree}],
  precision
];
assert[
  "native-batched-physical-columns-delta-sum",
  fakeEndpoint =!= $Failed &&
  maxError[fakeEndpoint, fakeExpected] < 10^-60
];

publicFakeDelta = SpectralPropagate[
  {scalarOperator},
  physicalBoundary,
  {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Native",
  "NativeHandle" -> handle,
  "DeltaTolerance" -> 10^-50,
  "MaxDeltaOrder" -> 64,
  "ResultData" -> "Full"
];
publicFakeDeltaEndpoint = SpectralPropagate[
  {scalarOperator},
  physicalBoundary,
  {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Native",
  "NativeHandle" -> handle,
  "DeltaTolerance" -> 10^-50,
  "MaxDeltaOrder" -> 64,
  "ResultData" -> "Endpoint"
];
mathematicaFakeDelta = SpectralPropagate[
  {scalarOperator},
  physicalBoundary,
  {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Mathematica",
  "DeltaTolerance" -> 10^-50,
  "MaxDeltaOrder" -> 64
];
assert[
  "unified-native-b0-auto-order",
  publicFakeDelta =!= $Failed && publicFakeDeltaEndpoint =!= $Failed &&
  mathematicaFakeDelta =!= $Failed &&
  publicFakeDelta[[6, 3]] === "Native" &&
  MissingQ[publicFakeDeltaEndpoint[[2]]] &&
  maxError[publicFakeDelta[[3]], mathematicaFakeDelta[[3]]] < 10^-55 &&
  maxError[publicFakeDeltaEndpoint[[3]], mathematicaFakeDelta[[3]]] < 10^-55
];

(* A prepared handle is an opaque cache for one sampled equation.  Reusing it
   with another evaluator must fail before any native propagation occurs. *)
ClearAll[wrongOperator];
wrongOperator[x_, requestedPrecision_] := N[{{1/2}}, requestedPrecision];
wrongHandleResult = Quiet @ SpectralPropagate[
  {wrongOperator}, physicalBoundary, {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Native",
  "NativeHandle" -> handle,
  "DeltaOrder" -> degree,
  "ResultData" -> "Endpoint"
];
invalidMaximumResult = Quiet @ SpectralPropagate[
  {scalarOperator}, physicalBoundary, {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Native",
  "NativeHandle" -> handle,
  "DeltaOrder" -> 4,
  "MaxDeltaOrder" -> 3.5,
  "ResultData" -> "Endpoint"
];
assert[
  "native-handle-system-binding-and-option-failure",
  wrongHandleResult === $Failed && invalidMaximumResult === $Failed
];

(* Reloading the loader must retain ownership of the installed process and its
   generation counter; otherwise a later generation-1 handle could resurrect a
   stale generation-1 handle while leaking the original WSTP process. *)
linksBeforeReload = Length[Links[]];
generationBeforeReload = $CHESSNativeGeneration;
Get[FileNameJoin[{testRoot, "CHESSv2.wl"}]];
assert[
  "native-package-reload-preserves-link-ownership",
  Length[Links[]] === linksBeforeReload &&
  $CHESSNativeGeneration === generationBeforeReload &&
  CHESSNativeHandleValidQ[handle]
];

(* Reusing the same symbol name after changing its DownValues must also be
   rejected; hashing only HoldComplete[scalarOperator] would miss this case. *)
Clear[scalarOperator];
scalarOperator[x_, requestedPrecision_] := N[{{1/2}}, requestedPrecision];
redefinedSymbolResult = Quiet @ SpectralPropagate[
  {scalarOperator}, physicalBoundary, {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Native",
  "NativeHandle" -> handle,
  "DeltaOrder" -> degree,
  "ResultData" -> "Endpoint"
];
assert[
  "native-handle-detects-symbol-redefinition",
  redefinedSymbolResult === $Failed
];
Clear[scalarOperator];
scalarOperator[x_, requestedPrecision_] := N[{{3/7}}, requestedPrecision];

(* The dependency fingerprint must be transitive.  Redefining a helper called
   by an otherwise unchanged outer evaluator changes its sampled matrices and
   therefore must invalidate the prepared handle. *)
ClearAll[helperOperator, outerOperator];
helperOperator[x_, requestedPrecision_] := N[1/3, requestedPrecision];
outerOperator[x_, requestedPrecision_] :=
  {{helperOperator[x, requestedPrecision]}};
helperHandle = CHESSNativePrepare[outerOperator, 1, {0, 1}, rules];
Clear[helperOperator];
helperOperator[x_, requestedPrecision_] := N[1/2, requestedPrecision];
helperRedefinitionResult = If[
  AssociationQ[helperHandle],
  Quiet @ SpectralPropagate[
    {outerOperator}, physicalBoundary, {0, 1},
    Sequence @@ rules,
    "NumericalBackend" -> "Native",
    "NativeHandle" -> helperHandle,
    "DeltaOrder" -> degree,
    "ResultData" -> "Endpoint"
  ],
  Missing["PrepareFailed"]
];
assert[
  "native-handle-detects-helper-redefinition",
  AssociationQ[helperHandle] && helperRedefinitionResult === $Failed
];

(* Options, DefaultValues, and NValues are definition state too.  In
   particular, OptionsPattern can change a sampled matrix while leaving all
   ordinary DownValues text unchanged. *)
ClearAll[optionOperator];
Options[optionOperator] = {"Coefficient" -> 1/3};
optionOperator[x_, requestedPrecision_, OptionsPattern[]] :=
  N[{{OptionValue["Coefficient"]}}, requestedPrecision];
optionHandle = CHESSNativePrepare[optionOperator, 1, {0, 1}, rules];
Options[optionOperator] = {"Coefficient" -> 1/2};
optionRedefinitionResult = If[
  AssociationQ[optionHandle],
  Quiet @ SpectralPropagate[
    {optionOperator}, physicalBoundary, {0, 1},
    Sequence @@ rules,
    "NumericalBackend" -> "Native",
    "NativeHandle" -> optionHandle,
    "DeltaOrder" -> degree,
    "ResultData" -> "Endpoint"
  ],
  Missing["PrepareFailed"]
];
assert[
  "native-handle-detects-option-redefinition",
  AssociationQ[optionHandle] && optionRedefinitionResult === $Failed
];

(* A separate complex operator forces the MPC branch.  Preparing it also
   invalidates the earlier handle, exercising the single-cache generation
   contract instead of accidentally reusing stale native matrices. *)
ClearAll[complexOperator];
complexOperator[x_, requestedPrecision_] :=
  N[{{2/5 + I/7}}, requestedPrecision];
complexHandle = CHESSNativePrepare[complexOperator, 1, {0, 1}, rules];
complexBoundary = N[{{1} ~Join~ ConstantArray[0, degree]}, precision];
complexResult = If[
  AssociationQ[complexHandle],
  SpectralPropagate[
    {0, complexOperator}, complexBoundary, {0, 1},
    Sequence @@ rules,
    "NumericalBackend" -> "Native",
    "NativeHandle" -> complexHandle,
    "ResultData" -> "Endpoint"
  ],
  $Failed
];
complexExpected = N[
  {Table[(2/5 + I/7)^order/order!, {order, 0, degree}]},
  precision
];
assert[
  "native-complex-mpc-and-handle-generation",
  complexResult =!= $Failed &&
  !CHESSNativeHandleValidQ[handle] &&
  maxError[complexResult[[3]], complexExpected] < 10^-60
];

(* A genuine B0+B1+B2 system exercises the direct active-support LU.  The
   complex B0 entry forces the FLINT nfloat-complex factorization; no fake-delta
   expansion is involved. *)
ClearAll[mixedB0, mixedB1, mixedB2, mixedScalar, mixedBatch];
mixedB0[x_, requestedPrecision_] := N[
  {
    {1/5 + x/20, 0, I (1 + x)/11},
    {0, 0, 0},
    {(1 - x/3)/13, 0, -1/7 + x/30}
  },
  requestedPrecision
];
mixedB1[x_, requestedPrecision_] := N[
  {{2/9, x/29, 0}, {1/17, -1/8, x/31}, {0, 1/37, 1/10}},
  requestedPrecision
];
mixedB2[x_, requestedPrecision_] := N[
  {{0, 1/19, 0}, {-1/23, 0, 1/41}, {1/43, 0, 0}},
  requestedPrecision
];
mixedScalar[x_, requestedPrecision_] := {
  mixedB0[x, requestedPrecision],
  mixedB1[x, requestedPrecision],
  mixedB2[x, requestedPrecision]
};
mixedBatch[points_, requestedPrecision_, threads_] :=
  mixedScalar[#, requestedPrecision] & /@ points;
mixedEvaluator = CHESSNodeEvaluator[mixedScalar, mixedBatch];
mixedBoundary = N[
  {{1, 1/3, 0, 0}, {-2, 0, 1/5, 0}, {1/7, -1/9, 0, 1/11}},
  precision
];
mixedRules = {
  "Nodes" -> 20,
  "Precision" -> precision,
  "WorkingPrecisionA" -> 90,
  "EpsilonDegree" -> 2,
  "ParallelEvaluation" -> False
};
mixedMathematica = SpectralPropagate[
  mixedEvaluator, mixedBoundary, {0, 1},
  Sequence @@ mixedRules,
  "NumericalBackend" -> "Mathematica"
];
mixedHandle = CHESSNativePolynomialPrepare[
  mixedEvaluator, mixedBoundary, {0, 1}, mixedRules
];
mixedNative = SpectralPropagate[
  mixedEvaluator, mixedBoundary, {0, 1},
  Sequence @@ mixedRules,
  "NumericalBackend" -> "Native",
  "NativeHandle" -> mixedHandle,
  "ResultData" -> "Endpoint"
];
mixedNativeFull = SpectralPropagate[
  mixedEvaluator, mixedBoundary, {0, 1},
  Sequence @@ mixedRules,
  "NumericalBackend" -> "Native",
  "NativeHandle" -> mixedHandle,
  "ResultData" -> "Full"
];
mixedNativeEndpointRequest = Quiet @ SpectralPropagate[
  mixedEvaluator, mixedBoundary, {0, 1},
  Sequence @@ mixedRules,
  "NumericalBackend" -> "Native",
  "RegularizedEndpoints" -> {0}
];
mixedListNative = Quiet @ SpectralPropagate[
  {mixedB0, mixedB1, mixedB2}, mixedBoundary, {0, 1},
  Sequence @@ mixedRules,
  "NumericalBackend" -> "Native"
];
assert[
  "native-direct-mixed-polynomial-complex",
  AssociationQ[mixedHandle] && mixedMathematica =!= $Failed &&
  mixedNative =!= $Failed && mixedNativeFull =!= $Failed &&
  mixedNative[[6, 4]] === "DirectCoupledLU" &&
  Dimensions[mixedNativeFull[[2]]] === {21, 12} &&
  mixedNativeEndpointRequest === $Failed && mixedListNative === $Failed &&
  maxError[mixedNativeFull[[2]], mixedMathematica[[2]]] < 10^-55 &&
  maxError[mixedNativeFull[[3]], mixedMathematica[[3]]] < 10^-55 &&
  maxError[mixedNative[[3]], mixedMathematica[[3]]] < 10^-55
];

CHESSNativeUninstall[];

If[failures === {},
  Print["CHESSV2_NATIVE_TESTS_PASSED"];
  Exit[0],
  Print[InputForm[{"CHESSV2_NATIVE_TEST_FAILURES", failures}]];
  Exit[1]
];
