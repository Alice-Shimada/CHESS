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
CheckAssertion[name_String, condition_] := If[
  TrueQ[condition],
  Print["PASS: ", name],
  AppendTo[failures, name];
  Print["FAIL: ", name]
];
MaxError[a_, b_] := Max[Abs[Flatten[N[a - b, 50]]]];

If[!CHESSNativeAvailableQ[],
  Print["FAIL: native-executable-missing"];
  Exit[2]
];

ClearAll[ScalarOperator];
ScalarOperator[x_, precision_] := N[{{3/7}}, precision];
precision = 70;
degree = 20;
rules = {
  "Nodes" -> 32,
  "Precision" -> precision,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False,
  "ParallelKernels" -> 1
};

handle = CHESSNativePrepare[ScalarOperator, 1, {0, 1}, rules];
CheckAssertion[
  "native-prepare",
  AssociationQ[handle] &&
  Lookup[handle, "Precision"] === precision &&
  Lookup[handle, "NativeGuardDigits"] === 20 &&
  Lookup[handle, "NativePrecision"] === precision + 20 &&
  Lookup[handle, "NativeWorkingPrecisionA"] === 100
];

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
CheckAssertion[
  "native-canonical-full-result",
  canonicalResult =!= $Failed &&
  Dimensions[canonicalResult[[2]]] === {33, degree + 1} &&
  MaxError[canonicalResult[[3]], canonicalExpected] < 10^-60
];

publicNative = SpectralPropagate[
  {0, ScalarOperator},
  canonicalBoundary,
  {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Native",
  "NativeHandle" -> handle,
  "ResultData" -> "Full"
];
publicEndpointOnly = SpectralPropagate[
  {0, ScalarOperator},
  canonicalBoundary,
  {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Native",
  "NativeHandle" -> handle,
  "ResultData" -> "Endpoint"
];
CheckAssertion[
  "unified-spectral-propagate-native-canonical",
  publicNative =!= $Failed && publicEndpointOnly =!= $Failed &&
  publicNative[[6, 2]] === "Native" &&
  MissingQ[publicEndpointOnly[[2]]] &&
  MaxError[publicNative[[3]], canonicalExpected] < 10^-60 &&
  MaxError[publicEndpointOnly[[3]], canonicalExpected] < 10^-60
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
CheckAssertion[
  "native-batched-physical-columns-delta-sum",
  fakeEndpoint =!= $Failed &&
  MaxError[fakeEndpoint, fakeExpected] < 10^-60
];

publicFakeDelta = SpectralPropagate[
  {ScalarOperator},
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
  {ScalarOperator},
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
  {ScalarOperator},
  physicalBoundary,
  {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Mathematica",
  "DeltaTolerance" -> 10^-50,
  "MaxDeltaOrder" -> 64
];
CheckAssertion[
  "unified-native-b0-auto-order",
  publicFakeDelta =!= $Failed && publicFakeDeltaEndpoint =!= $Failed &&
  mathematicaFakeDelta =!= $Failed &&
  publicFakeDelta[[6, 3]] === "Native" &&
  MissingQ[publicFakeDeltaEndpoint[[2]]] &&
  MaxError[publicFakeDelta[[3]], mathematicaFakeDelta[[3]]] < 10^-55 &&
  MaxError[publicFakeDeltaEndpoint[[3]], mathematicaFakeDelta[[3]]] < 10^-55
];

(* A prepared handle is an opaque cache for one sampled equation.  Reusing it
   with another evaluator must fail before any native propagation occurs. *)
ClearAll[WrongOperator];
WrongOperator[x_, requestedPrecision_] := N[{{1/2}}, requestedPrecision];
wrongHandleResult = Quiet @ SpectralPropagate[
  {WrongOperator}, physicalBoundary, {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Native",
  "NativeHandle" -> handle,
  "DeltaOrder" -> degree,
  "ResultData" -> "Endpoint"
];
invalidMaximumResult = Quiet @ SpectralPropagate[
  {ScalarOperator}, physicalBoundary, {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Native",
  "NativeHandle" -> handle,
  "DeltaOrder" -> 4,
  "MaxDeltaOrder" -> 3.5,
  "ResultData" -> "Endpoint"
];
CheckAssertion[
  "native-handle-system-binding-and-option-failure",
  wrongHandleResult === $Failed && invalidMaximumResult === $Failed
];

(* Reloading the loader must retain ownership of the installed process and its
   generation counter; otherwise a later generation-1 handle could resurrect a
   stale generation-1 handle while leaking the original WSTP process. *)
linksBeforeReload = Length[Links[]];
generationBeforeReload = $CHESSNativeGeneration;
Get[FileNameJoin[{testRoot, "CHESSv2.wl"}]];
CheckAssertion[
  "native-package-reload-preserves-link-ownership",
  Length[Links[]] === linksBeforeReload &&
  $CHESSNativeGeneration === generationBeforeReload &&
  CHESSNativeHandleValidQ[handle]
];

(* Reusing the same symbol name after changing its DownValues must also be
   rejected; hashing only HoldComplete[ScalarOperator] would miss this case. *)
Clear[ScalarOperator];
ScalarOperator[x_, requestedPrecision_] := N[{{1/2}}, requestedPrecision];
redefinedSymbolResult = Quiet @ SpectralPropagate[
  {ScalarOperator}, physicalBoundary, {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Native",
  "NativeHandle" -> handle,
  "DeltaOrder" -> degree,
  "ResultData" -> "Endpoint"
];
CheckAssertion[
  "native-handle-detects-symbol-redefinition",
  redefinedSymbolResult === $Failed
];
Clear[ScalarOperator];
ScalarOperator[x_, requestedPrecision_] := N[{{3/7}}, requestedPrecision];

(* The dependency fingerprint must be transitive.  Redefining a helper called
   by an otherwise unchanged outer evaluator changes its sampled matrices and
   therefore must invalidate the prepared handle. *)
ClearAll[HelperOperator, OuterOperator];
HelperOperator[x_, requestedPrecision_] := N[1/3, requestedPrecision];
OuterOperator[x_, requestedPrecision_] :=
  {{HelperOperator[x, requestedPrecision]}};
helperHandle = CHESSNativePrepare[OuterOperator, 1, {0, 1}, rules];
Clear[HelperOperator];
HelperOperator[x_, requestedPrecision_] := N[1/2, requestedPrecision];
helperRedefinitionResult = If[
  AssociationQ[helperHandle],
  Quiet @ SpectralPropagate[
    {OuterOperator}, physicalBoundary, {0, 1},
    Sequence @@ rules,
    "NumericalBackend" -> "Native",
    "NativeHandle" -> helperHandle,
    "DeltaOrder" -> degree,
    "ResultData" -> "Endpoint"
  ],
  Missing["PrepareFailed"]
];
CheckAssertion[
  "native-handle-detects-helper-redefinition",
  AssociationQ[helperHandle] && helperRedefinitionResult === $Failed
];

(* Options, DefaultValues, and NValues are definition state too.  In
   particular, OptionsPattern can change a sampled matrix while leaving all
   ordinary DownValues text unchanged. *)
ClearAll[OptionOperator];
Options[OptionOperator] = {"Coefficient" -> 1/3};
OptionOperator[x_, requestedPrecision_, OptionsPattern[]] :=
  N[{{OptionValue["Coefficient"]}}, requestedPrecision];
optionHandle = CHESSNativePrepare[OptionOperator, 1, {0, 1}, rules];
Options[OptionOperator] = {"Coefficient" -> 1/2};
optionRedefinitionResult = If[
  AssociationQ[optionHandle],
  Quiet @ SpectralPropagate[
    {OptionOperator}, physicalBoundary, {0, 1},
    Sequence @@ rules,
    "NumericalBackend" -> "Native",
    "NativeHandle" -> optionHandle,
    "DeltaOrder" -> degree,
    "ResultData" -> "Endpoint"
  ],
  Missing["PrepareFailed"]
];
CheckAssertion[
  "native-handle-detects-option-redefinition",
  AssociationQ[optionHandle] && optionRedefinitionResult === $Failed
];

(* A separate complex operator forces the MPC branch.  Preparing it also
   invalidates the earlier handle, exercising the single-cache generation
   contract instead of accidentally reusing stale native matrices. *)
ClearAll[ComplexOperator];
ComplexOperator[x_, requestedPrecision_] :=
  N[{{2/5 + I/7}}, requestedPrecision];
complexHandle = CHESSNativePrepare[ComplexOperator, 1, {0, 1}, rules];
complexBoundary = N[{{1} ~Join~ ConstantArray[0, degree]}, precision];
complexResult = If[
  AssociationQ[complexHandle],
  SpectralPropagate[
    {0, ComplexOperator}, complexBoundary, {0, 1},
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
CheckAssertion[
  "native-complex-mpc-and-handle-generation",
  complexResult =!= $Failed &&
  !CHESSNativeHandleValidQ[handle] &&
  MaxError[complexResult[[3]], complexExpected] < 10^-60
];

(* A genuine B0+B1+B2 system exercises the direct active-support LU.  The
   complex B0 entry forces the FLINT nfloat-complex factorization; no fake-delta
   expansion is involved. *)
ClearAll[MixedB0, MixedB1, MixedB2];
MixedB0[x_, requestedPrecision_] := N[
  {
    {1/5 + x/20, 0, I (1 + x)/11},
    {0, 0, 0},
    {(1 - x/3)/13, 0, -1/7 + x/30}
  },
  requestedPrecision
];
MixedB1[x_, requestedPrecision_] := N[
  {{2/9, x/29, 0}, {1/17, -1/8, x/31}, {0, 1/37, 1/10}},
  requestedPrecision
];
MixedB2[x_, requestedPrecision_] := N[
  {{0, 1/19, 0}, {-1/23, 0, 1/41}, {1/43, 0, 0}},
  requestedPrecision
];
mixedEvaluator = CHESSNativeMatrixAdapter[
  {MixedB0, MixedB1, MixedB2}, 3
];
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
  {MixedB0, MixedB1, MixedB2}, mixedBoundary, {0, 1},
  Sequence @@ mixedRules,
  "NumericalBackend" -> "Native"
];
CheckAssertion[
  "native-direct-mixed-polynomial-complex",
  AssociationQ[mixedHandle] && mixedMathematica =!= $Failed &&
  mixedNative =!= $Failed && mixedNativeFull =!= $Failed &&
  mixedNative[[6, 4]] === "DirectCoupledLU" &&
  Dimensions[mixedNativeFull[[2]]] === {21, 12} &&
  mixedNativeEndpointRequest === $Failed && mixedListNative === $Failed &&
  Lookup[mixedHandle, "NativeGuardDigits"] === 20 &&
  Lookup[mixedHandle, "NativePrecision"] === precision + 20 &&
  MemberQ[mixedNative[[6]], "NativeWorkingPrecision"] &&
  MaxError[mixedNativeFull[[2]], mixedMathematica[[2]]] < 10^-55 &&
  MaxError[mixedNativeFull[[3]], mixedMathematica[[3]]] < 10^-55 &&
  MaxError[mixedNative[[3]], mixedMathematica[[3]]] < 10^-55
];

(* NativeGuardDigits raises only the backend working precision.  The requested
   public Precision remains unchanged, participates in handle compatibility,
   and may be restored to the historical no-guard behavior explicitly. *)
zeroGuardRules = Join[rules, {"NativeGuardDigits" -> 0}];
zeroGuardHandle = CHESSNativePrepare[
  ScalarOperator, 1, {0, 1}, zeroGuardRules
];
zeroGuardResult = If[
  AssociationQ[zeroGuardHandle],
  SpectralPropagate[
    {0, ScalarOperator}, canonicalBoundary, {0, 1},
    Sequence @@ zeroGuardRules,
    "NumericalBackend" -> "Native",
    "NativeHandle" -> zeroGuardHandle,
    "ResultData" -> "Endpoint"
  ],
  $Failed
];
invalidGuardResult = Quiet @ SpectralPropagate[
  {0, ScalarOperator}, canonicalBoundary, {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Native",
  "NativeGuardDigits" -> -1,
  "ResultData" -> "Endpoint"
];
CheckAssertion[
  "native-guard-precision-contract",
  AssociationQ[zeroGuardHandle] && zeroGuardResult =!= $Failed &&
  Lookup[zeroGuardHandle, "Precision"] === precision &&
  Lookup[zeroGuardHandle, "NativePrecision"] === precision &&
  Lookup[zeroGuardHandle, "NativeGuardDigits"] === 0 &&
  !CHESSNativeHandleCompatibleQ[
    zeroGuardHandle, ScalarOperator, 1, {0, 1}, rules
  ] &&
  MemberQ[zeroGuardResult[[6]], "NativeWorkingPrecision"] &&
  invalidGuardResult === $Failed
];

(* Machine-precision nonzero input must fail closed.  Printing it with more
   digits would only manufacture a high-precision-looking C++ result.  Exact
   zeros remain valid because they carry no numerical uncertainty. *)
ClearAll[MachineOperator];
MachineOperator[x_, requestedPrecision_] := {{1.0/3.0}};
machineHandle = Quiet @ CHESSNativePrepare[
  MachineOperator, 1, {0, 1}, rules
];
machineBoundary = canonicalBoundary;
machineBoundary[[1, 1]] = 1.0;
machineBoundaryRun = Quiet @ CHESSNativeRun[
  zeroGuardHandle, {Transpose[machineBoundary]}, False
];
machineIntervalHandle = Quiet @ CHESSNativePrepare[
  ScalarOperator, 1, {0., 1.}, rules
];
machineMixedEvaluator = CHESSNativeMatrixAdapter[
  {MachineOperator, ScalarOperator}, 1
];
machineMixedRules = Join[
  rules, {"EpsilonDegree" -> 1}
];
machineMixedHandle = Quiet @ CHESSNativePolynomialPrepare[
  machineMixedEvaluator, canonicalBoundary, {0, 1}, machineMixedRules
];
CheckAssertion[
  "native-low-precision-input-fails-closed",
  machineHandle === $Failed && machineBoundaryRun === $Failed &&
  machineIntervalHandle === $Failed && machineMixedHandle === $Failed &&
  !CHESSNativeMatrixPrecisionQ[{{1.0/3.0}}, precision] &&
  CHESSNativeMatrixPrecisionQ[{{0., 1/3}}, precision]
];

(* Segmented B0 propagation must construct every subinterval at guarded
   precision and carry the unrounded endpoint internally to the next segment.
   Only the final public state is rounded back to the requested Precision. *)
ClearAll[SegmentOperator];
segmentSamplePrecisions = {};
SegmentOperator[x_, requestedPrecision_?NumericQ] := (
  AppendTo[segmentSamplePrecisions, {Precision[x], requestedPrecision}];
  N[{{1/3}}, requestedPrecision]
);
segmentedNative = SpectralPropagate[
  {SegmentOperator}, N[{{1}}, precision], {0, 1},
  "Nodes" -> 18,
  "Precision" -> precision,
  "WorkingPrecisionA" -> 90,
  "NativeGuardDigits" -> 20,
  "NumericalBackend" -> "Native",
  "Segments" -> 3,
  "DeltaOrder" -> 32,
  "MaxDeltaOrder" -> 32,
  "ResultData" -> "Endpoint",
  "ParallelEvaluation" -> False
];
CheckAssertion[
  "native-segmented-guard-precision",
  segmentedNative =!= $Failed && segmentSamplePrecisions =!= {} &&
  Min[segmentSamplePrecisions[[All, 1]]] > precision + 10 &&
  Union[segmentSamplePrecisions[[All, 2]]] === {
    Max[90, precision + 20 + 10]
  } &&
  Count[Flatten[segmentedNative[[6]]], "NativeWorkingPrecision"] === 3 &&
  Abs[segmentedNative[[3, 1, 1]] - Exp[N[1/3, precision]]] < 10^-40
];

CHESSNativeUninstall[];

If[failures === {},
  Print["CHESSV2_NATIVE_TESTS_PASSED"];
  Exit[0],
  Print[InputForm[{"CHESSV2_NATIVE_TEST_FAILURES", failures}]];
  Exit[1]
];
